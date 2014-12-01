##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

## May 2010 Dream3/4 pipeline (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
##  		     "Alex Greenfield" <ag1868@nyu.edu> 
## NYU - Center for Genomics and Systems Biology

##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

clr <- function(m) {
	repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}
	ur = apply(m, 1, mean)
	uc = apply(m, 2, mean)
	sr = apply(m, 1, sd)
	sc = apply(m, 2, sd)
	
	uc = t(repmat(uc, 1, nrow(m)))
	sc = t(repmat(sc, 1, nrow(m)))
	
	z1 = (m - ur) / sr
	
	z2 = (m - uc) / sc
	
	z1[which(z1 < 0)] = 0
	z2[which(z2 < 0)] = 0
	
	z = sqrt(z1 ^ 2 + z2 ^ 2) 
	
	return(z)
}

# input:
#	- m.stat matrix of static MI
#	- m.dyn matrix of dynamic MI
#
#
mixed_clr <- function(m.stat,m.dyn) {
	Z = matrix(0,nrow=dim(m.stat)[1],ncol=dim(m.stat)[2])
	stat.sd <- apply(m.stat, 2, sd)
	stat.mean <- apply(m.stat, 2, mean)
	dyn.sd <- apply(m.dyn, 1, sd)
	dyn.mean <- apply(m.dyn, 1, mean)	 
	ind = 1:dim(m.stat)[1]
	x = lapply(ind, mixed_clr_level1, m = m.dyn, dyn.mean = dyn.mean, dyn.sd = dyn.sd, stat.mean = stat.mean, stat.sd = stat.sd)
	for (i in 1:length(x)) {
		Z[i,] = x[[i]]
	}	
	cat(i)
	# coerce mixed clr results to matrix (1 row for each target)
	return(Z)
}

mixed_clr_level1 <- function(ix, m, dyn.mean, dyn.sd, stat.mean, stat.sd) {
	z_ix = numeric(length=dim(m)[2])
	ind = 1:dim(m)[2]
	z_ix = as.numeric(lapply(ind, mixed_clr_level2, dyn.mi.val = m[ix,], dyn.mean.val = dyn.mean[ix], dyn.sd.val = dyn.sd[ix], stat.mean = stat.mean, stat.sd = stat.sd))
	cat(".")
	return(z_ix)
}

mixed_clr_level2 <- function(ix, dyn.mi.val, dyn.mean.val, dyn.sd.val, stat.mean, stat.sd) {
	z.stat = max(0, (dyn.mi.val[ix] - stat.mean[ix])/stat.sd[ix] )
	z.dyn = max(0, (dyn.mi.val[ix] - dyn.mean.val)/dyn.sd.val )
	return(sqrt(z.stat ^ 2 + z.dyn ^ 2))
}

mixed_clr_Stouffer <- function(m.stat,m.dyn) {
	Z = matrix(0,nrow=dim(m.stat)[1],ncol=dim(m.stat)[2])
	stat.sd <- apply(m.stat, 2, sd)
	stat.mean <- apply(m.stat, 2, mean)
	dyn.sd <- apply(m.dyn, 1, sd)
	dyn.mean <- apply(m.dyn, 1, mean)	 
	ind = 1:dim(m.stat)[1]
	x = lapply(ind, mixed_clr_level1_Stouffer, m = m.dyn, dyn.mean = dyn.mean, dyn.sd = dyn.sd, stat.mean = stat.mean, stat.sd = stat.sd)
	for (i in 1:length(x)) {
		Z[i,] = x[[i]]
	}	
	cat(i)
	# coerce mixed clr results to matrix (1 row for each target)
	return(Z)
}

mixed_clr_level1_Stouffer <- function(ix, m, dyn.mean, dyn.sd, stat.mean, stat.sd) {
	z_ix = numeric(length=dim(m)[2])
	ind = 1:dim(m)[1]
	z_ix = as.numeric(lapply(ind, mixed_clr_level2_Stouffer, dyn.mi.val = m[ix,], dyn.mean.val = dyn.mean[ix], dyn.sd.val = dyn.sd[ix], stat.mean = stat.mean, stat.sd = stat.sd))
	cat(".")
	return(z_ix)
}

mixed_clr_level2_Stouffer <- function(ix, dyn.mi.val, dyn.mean.val, dyn.sd.val, stat.mean, stat.sd) {
	z.stat = (dyn.mi.val[ix] - stat.mean[ix])/stat.sd[ix]
	z.dyn = (dyn.mi.val[ix] - dyn.mean.val)/dyn.sd.val
	return(sum(z.stat,z.dyn)/sqrt(2))
}

mixed_clr_parallel <- function(m.stat,m.dyn, processorsNumber = 1) {
	if(processorsNumber == 1){
		Z = mixed_clr(m.stat,m.dyn)
		return(Z)
	}
	
	Z = matrix(0,nrow=dim(m.stat)[1],ncol=dim(m.stat)[2])
	stat.sd <- apply(m.stat, 2, sd)
	stat.mean <- apply(m.stat, 2, mean)
	dyn.sd <- apply(m.dyn, 1, sd)
	dyn.mean <- apply(m.dyn, 1, mean)	 
	ind = 1:dim(m.dyn)[1]
	x = mclapply(ind, mixed_clr_level1, m = m.dyn, dyn.mean = dyn.mean, dyn.sd = dyn.sd, stat.mean = stat.mean, stat.sd = stat.sd, mc.cores=processorsNumber)
	for (i in 1:length(x)) {
		Z[i,] = x[[i]]
	}	
	cat(i)
	# coerce mixed clr results to matrix (1 row for each target)
	return(Z)
}

mixed_clr_parallel_Stouffer <- function(m.stat,m.dyn, processorsNumber = 1) {
	Z = matrix(0,nrow=dim(m.stat)[1],ncol=dim(m.stat)[2])
	stat.sd <- apply(m.stat, 2, sd)
	stat.mean <- apply(m.stat, 2, mean)
	dyn.sd <- apply(m.dyn, 1, sd)
	dyn.mean <- apply(m.dyn, 1, mean)	 
	ind = 1:dim(m.stat)[1]
	x = mclapply(ind, mixed_clr_level1_Stouffer, m = m.dyn, dyn.mean = dyn.mean, dyn.sd = dyn.sd, stat.mean = stat.mean, stat.sd = stat.sd, mc.cores=processorsNumber)
	for (i in 1:length(x)) {
		Z[i,] = x[[i]]
	}	
	cat(i)
	# coerce mixed clr results to matrix (1 row for each target)
	return(Z)
}

mixed_clr_one_by_one <- function(dyn.mi.val, stat.mean.val, dyn.mean.val, stat.sd.val, dyn.sd.val ) {
	z.stat = max(0, (dyn.mi.val - stat.mean.val)/stat.sd.val )
	z.dyn = max(0, (dyn.mi.val - dyn.mean.val)/dyn.sd.val )
	return(sqrt(z.stat ^ 2 + z.dyn ^ 2))
}



# calc MI matrix in batches
calc_MI_inBatces <- function(Y,X,btch_size, n.bins=10) {
	M = matrix(0,nrow=nrow(Y),ncol=nrow(X))
	x= nrow(Y)
	y=1
	cat("mi calculation for: ")
	while(x > btch_size){ 
		cat(y,":",(y+(btch_size-1)),", ", sep="") ; 		flush.console()
		M[y:(y+(btch_size-1)),] = mi(Y[y:(y+(btch_size-1)),], X,n=n.bins)
		x = x - btch_size
		y=y+btch_size
	}
	cat(y,":",(y+x-1),"\n", sep="") ; 	flush.console()
	M[y:(y+x-1),] = mi(Y[y:(y+x-1),], X,n=n.bins)
	return(M)
}


# like calc MI one by one but parallalized
calc_MI_one_by_one_parallel <- function(Y, X, Pi, processorsNumber = 1, n.bins=10, verbose = F){
        mem.start <- memory.profile()
	M = matrix(0,nrow=nrow(Y),ncol=nrow(X))
	ind = as.vector(c(1:nrow(Y)))
	repeat {
		x = mclapply(ind, mi_parllel, Y, X, Pi, mc.cores = processorsNumber, n.bins=n.bins)
		if(length(which(lapply(x, is.null) == T)) == 0)
			break;
		cat("\nREPEATING\n")
		print(which(lapply(x, is.null) == T))
	}
	if(verbose){
    cat("\ncalc_MI: before copy (diff)\n")
    print(memory.profile() - mem.start)
	}
	for (i in 1:length(x)) {
		M[i,] = x[[i]]
	}
	cat(i)
	if(verbose){
    cat("\ncalc_MI: after copy (diff)\n")
    print(memory.profile() - mem.start)
	}
	return(M)
}


mi_parllel <- function(ix, Y_ind, X, Pi_ind, n.bins=10){
	#for debugging only:
	#	cat(dim(Y_ind),"\t")
	#	cat(dim(Pi_ind),"\t")
	#	cat(dim(X),"\n")
	
	# old call: 
	#vec = mi(t(as.matrix(Y_ind[ix,Pi_ind[ix,]])), X[,Pi_ind[ix,]], n=n.bins )
	
	# if we have few conditions, add some fuzz
	fuzz.n <- ceiling(165 / length(Pi_ind[ix,]))

	# odd numbers of fuzz work better
	if (fuzz.n %% 2 == 0) {
	  fuzz.n <- fuzz.n - 1
	}
	
	vec = miNew(Y_ind[ix, Pi_ind[ix,]], X[, Pi_ind[ix, ]], n.bins, fuzz.n = fuzz.n)
	
	cat(".")
	return(vec)
}

mi_parllel_old <- function(ix, Y_ind, X, sigma, n.bins=10){
	vec = mi(t(as.matrix(Y_ind[ix, ])), X, n=n.bins )
	cat(".")
	return(vec)
}


calc_MI_one_by_one_parallel_old <- function(Y, X, processorsNumber = 1, n.bins=10, verbose = F){
	mem.start <- memory.profile()
	M = matrix(0,nrow=nrow(Y),ncol=nrow(X))
	ind = as.vector(c(1:nrow(Y)))
	repeat {
		x = mclapply(ind, mi_parllel,Y[ind, ], X, mc.cores=processorsNumber,n.bins=n.bins)
		if(length(which(lapply(x, is.null) == T)) == 0)
			break;
		cat("\nREPEATING\n")
		print(which(lapply(x, is.null) == T))
	}
	if(verbose){
		cat("\ncalc_MI: before copy (diff)\n")
		print(memory.profile() - mem.start)
	}
	for (i in 1:length(x)) {
		M[i,] = x[[i]]
	}
	cat(i)
	if(verbose){
		cat("\ncalc_MI: after copy (diff)\n")
		print(memory.profile() - mem.start)
	}
	return(M)
}

## below this line: helper functions for MI by CH 03/19/2012
################################################################################

miNew <- function(x, Y, nbins, fuzz.sdf = 0.5, fuzz.n = 1) {
  # calculate MI of x to each column of Y
  # input:
  # x - vector
  # Y - matrix
  # nbins - number if bins to use for discretization
  # fuzz.sdf - sd of fuzz added
  # fuzz.n - number of fuzz points added
  # output:
  # vector of MI values
  
  # fuzzyfy
  if(length(x) <= 25){
		x <- fuzzyfy(x, fuzz.sdf, fuzz.n)
  	Y <- t(apply(Y, 1, fuzzyfy, fuzz.sdf, fuzz.n))
	}
  # discretize
  x <- discretize(x, nbins)
  Y <- t(apply(Y, 1, discretize, nbins))
  
  n <- nrow(Y)
  ret <- rep(0, n)
  for (i in 1:n) {
    ret[i] <- discMI(x, Y[i, ], nbins)
  }
  return(ret)
}

fuzzyfy <- function(x, sdf, n) {
  if (n <= 1) {
    return(x)
  }
  N <- length(x)
  sd <- sd(x) * sdf
  ret <- rep(0, N * n)
  fuzz <- nrnorm(n, 0, sd)
  
  for (i in 1:N) {
    ret[((i - 1) * n + 1):(i * n)] <- x[i] + fuzz
  }
  
  return(ret)
}

nrnorm <- function(K, mean = 0, sd = 1) {
  # non-random points from a normal distribution
  return(qnorm((1 / (2 * K)) + ((0:(K - 1)) / K)) * sd + mean)
}

discretize <- function(X, nbins, eps = 0.0001) {
  N <- length(X)
  X.min <- min(X)
  X.max <- max(X)
  X.disc <- ceiling((X - X.min + eps) / (X.max - X.min + eps) * nbins)
  return(as.integer(X.disc))
}

discMI <- function(X.disc, Y.disc, nlevels) {
  # return the MI of two integer vectors with 0 < values <= nlevels
  tmp <- discMIhelper(length(X.disc), nlevels, X.disc, Y.disc, 
                      rep(0, nlevels * nlevels), rep(0, nlevels * nlevels))
  return(sum(tmp$pxy * log(tmp$pxy / tmp$pxpy), na.rm = T))
}

discMIhelper.sig <- signature(n = 'integer', k = 'integer', x = 'integer', 
                              y = 'integer', pxy = 'numeric', 
                              pxpy = 'numeric')

# C code for computing the frequency and joint frequency of integer variables
discMIhelper.code <- '
  int i, j;
  double *px, *py;

  px = (double *)calloc(*k, sizeof(double));
  py = (double *)calloc(*k, sizeof(double));

  for (i = 0; i < *n; ++i) {
    ++pxy[(y[i]-1) * (*k) + (x[i]-1)];
    ++px[x[i]-1];
    ++py[y[i]-1];
  }
  for (i = 0; i < *k; ++i) {
    for (j = 0; j < *k; ++j) {
      pxy[j * (*k) + i] /= *n;
      pxpy[j * (*k) + i] = (px[i] / *n) * (py[j] / *n);
    }
  }

  free(px);
  free(py);'

discMIhelper <- cfunction(discMIhelper.sig, discMIhelper.code, convention=".C")


