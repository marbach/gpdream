load.net.avg <- function(netNum) {
  return(read.table(paste("Network", netNum, "_Data/net", netNum, "_expression_data_avg_t.tsv", sep="")))
}

load.net.full <- function(netNum) {
  return(read.table(paste("Network", netNum, "_Data/net", netNum, "_expression_data_t.tsv", sep="")))
}

load.tfs <- function(netNum) {
  return(as.character(read.table(paste("Network", netNum, "_Data/Network", netNum, "_transcription_factors.tsv", sep=""))[[1]]))
}

fix.colname <- function(m) {
  m2 <- m
  v2 <- array(length(m2))
  for(a in 1:length(m2))
    v2[a] <- sub("_.*$", paste("_", a, sep=""), colnames(m2[a]))
  colnames(m2) <- v2
  return(m2)
}

# shows a heat map for a symmetric matrix with reordering but without the actual dendrogram display
heatmap.crispy <- function(c, as.image=FALSE, do.reorder=TRUE, ...) {
  if(do.reorder) {
    dd <- as.dendrogram(hclust(dist(c)))
    dd2 <- reorder(dd, rowMeans(c))
    idx <- order.dendrogram(dd2)
  } else {
    idx <- 1:nrow(c)
  }
  if(as.image)
    image(c[idx,rev(idx)], axes=FALSE, ...)
  else
    heatmap(c[idx,idx], Rowv=NA, Colv=NA, revC=TRUE, scale="none", ...)
}

heat.net <- function(m, fix=TRUE, ...) {
  if(fix)
    m2 <- fix.colname(m)
  else
    m2 <- m
  m2c <- cor(m2)
  heatmap.crispy(m2c, ...)
  return(m2c)
}

get.expnum <- function(m) {
  num <- 0
  v <- array()
  for (a in 1:dim(m)[2]) {
    n <- sub("^.*([0-9]+)", "\\1", colnames(m)[a])
    if(n == 1)
      num <- num + 1
    v[a] <- num
  }
  return(factor(v))
}

s.to.n <- function(m, v=get.expnum(m), avg=NULL) {
  r <- t(apply(m, 1, function(x, v) tapply(x, v, sd), v))             # conditions with one rep will have sd of NA
  med <- apply(r, 1, median, na.rm=TRUE)  # so for each gene, take median of non-NA values
  if(is.null(avg))
    avg <- t(apply(m, 1, tmean <- function(x, v) tapply(x, v, mean), v))
  return(sd(t(avg)) / med)
}

# compares correlation coefficients between and among two sets of conditions
compare.cond <- function(m, g1, g2, exp.g1=NULL, exp.g2=NULL, ...) {
  m <- fix.colname(m)
  if(is.null(exp.g1))
    g1.idx <- g1
  else
    g1.idx <- grep(paste("^X", exp.g1, "_", sep=""), colnames(m))
  if(is.null(exp.g2))
    g2.idx <- g2
  else
    g2.idx <- grep(paste("^X", exp.g2, "_", sep=""), colnames(m))
  idx.sub <- c(g1.idx, g2.idx)
  m.sub <- m[,idx.sub]
  c <- heat.net(m.sub, fix=FALSE, ...)
  corr.mean(c)
}

compare.cond.2 <- function(m, exp1, exp2, ...) {
  g1.idx <- grep(paste("^X", exp1, "_", sep=""), colnames(m))
  g2.idx <- grep(paste("^X", exp2, "_", sep=""), colnames(m))
  idx.sub <- c(g1.idx, g2.idx)
  m.sub <- m[,idx.sub]
  colnames(m.sub) <- gsub("^([^_]*)_[^_]*_[^_]*_[^_]*_([^_]*)_.*", "\\1_\\2", colnames(m.sub))
  m.ord <- m.sub[,sort(colnames(m.sub))]
  c <- cor(m.ord)
  heatmap.crispy(c, do.reorder=FALSE, ...)
  N <- nrow(c) / 2
  c(corr.mean(c), mean(diag(c[1:N,(N+1):(2*N)])))
}

# find means between and within 2 same-sized groups of condition in a correlation matrix
corr.mean <- function(m) {
  N <- nrow(m) / 2
  # easy part: between groups
  bt <- mean(m[1:N,(N+1):(2*N)])
  # harder: within groups, exluding diagonal
  g1 <- mean((m[1:N,1:N])[lower.tri(m[1:N,1:N])])
  g2 <- mean((m[(N+1):(2*N),(N+1):(2*N)])[lower.tri(m[(N+1):(2*N),(N+1):(2*N)])])
  wi <- mean(c(g1, g2))
  return(c(wi, bt))
}

make.plots2 <- function() {
  net2.avg <- load.net.avg(2)
  net2.full <- load.net.full(2)
  # network 2: overall heatmap
  X11()
  c2 <- heat.net(net2.avg, main="Net 2 experiment correlation, all genes (2800)")
  # network 2: detailed experiment comparison, X23 & X25
  X11()
  net2.fix <- fix.colname(net2.avg)
  net2.small <- net2.fix[,c("X23_42","X23_43","X23_44","X23_45","X25_50","X25_51","X25_52","X25_53")]
  rn <- c("X23.WT","X23.KO1","X23.KO2","X23.DKO","X25.WT","X25.KO1","X25.KO2","X25.DKO")
  colnames(net2.small) <- rn
  c2.small <- heat.net(net2.small, fix=FALSE, main="Net 2 experiment vs. condition correlation")
  print(paste("inter-experiment mean:", mean(c2.small[5:8,1:4])))
  i <- lower.tri(c2.small)
  i[5:8,1:4] <- FALSE
  print(paste("intra-experiment mean:", mean(c2.small[i])))
  # network 2: s-to-n histogram
  X11()
  sn <- s.to.n(net2.full)
  hist(sn, nclass=500, main="Net 2 signal-to-noise")
  # network 2: above-mean s-to-n heatmap
  X11()
  q.thr <- .5
  sn2.keep <- names(sn2[sn2 > quantile(sn2, .5)])
  net2.high <- net2.avg[sn2.keep,]
  c2.high <- heat.net(net2.high, main="Net 2 experiment correlation, high s-to-n genes (1400)")
}
  
make.plots1 <- function() {
  net1.avg <- load.net.avg(1)
  net1.full <- load.net.full(1)
  sn <- s.to.n(net1.full)
  hist(sn, nclass=500, main="Net 1 signal-to-noise")
}

plot.comp1 <- function(as.image=FALSE, ...) {
  if(!as.image) X11()
  c <- compare.cond(net1.avg, exp.g1=70, exp.g2=87, as.image=as.image,
                    main="net1 70 vs 87", ...)
  print(c)
  if(!as.image) X11()
  c <- compare.cond(net1.avg, exp.g1=68, exp.g2=77, as.image=as.image,
                    main="net1 68 vs 77", ...)
  print(c)
  if(!as.image) X11()
  c <- compare.cond(net1.avg, exp.g1=72, exp.g2=95, as.image=as.image,
                    main="net1 72 vs 95", ...)
  print(c)
}

plot.comp13.2 <- function() {
  zl <- c(0.5, 1.0)
  par(mfrow=c(2,3), mar=c(1,1,1,1))
  print(compare.cond.2(net1.avg, 70, 87, as.image=TRUE, main="net1 70 vs 87", zlim=zl))
  print(compare.cond.2(net1.avg, 68, 77, as.image=TRUE, main="net1 68 vs 77", zlim=zl))
  print(compare.cond.2(net1.avg, 72, 95, as.image=TRUE, main="net1 72 vs 95", zlim=zl))
  print(compare.cond.2(net3.avg, 64, 86, as.image=TRUE, main="net3 64 vs 86", zlim=zl))
  print(compare.cond.2(net3.avg, 71, 90, as.image=TRUE, main="net3 71 vs 90", zlim=zl))
  print(compare.cond.2(net3.avg, 82, 92, as.image=TRUE, main="net3 82 vs 92", zlim=zl))
}

plot.comp13.2.high <- function() {
  zl <- c(0.5, 1.0)
  par(mfrow=c(2,3), mar=c(1,1,1,1))
  print(compare.cond.2(net1.high, 70, 87, as.image=TRUE, main="net1 70 vs 87", zlim=zl))
  print(compare.cond.2(net1.high, 68, 77, as.image=TRUE, main="net1 68 vs 77", zlim=zl))
  print(compare.cond.2(net1.high, 72, 95, as.image=TRUE, main="net1 72 vs 95", zlim=zl))
  print(compare.cond.2(net3.high, 64, 86, as.image=TRUE, main="net3 64 vs 86", zlim=zl))
  print(compare.cond.2(net3.high, 71, 90, as.image=TRUE, main="net3 71 vs 90", zlim=zl))
  print(compare.cond.2(net3.high, 82, 92, as.image=TRUE, main="net3 82 vs 92", zlim=zl))
}

plot.comp3 <- function(as.image=FALSE, ...) {
  if(!as.image) X11()
  c <- compare.cond(net3.avg, exp.g1=64, exp.g2=86, as.image=as.image,
                    main="net3 64 vs 86", ...)
  print(c)
  if(!as.image) X11()
  c <- compare.cond(net3.avg, exp.g1=71, exp.g2=90, as.image=as.image,
                    main="net3 71 vs 90", ...)
  print(c)
  if(!as.image) X11()
  c <- compare.cond(net3.avg, exp.g1=82, exp.g2=92, as.image=as.image,
                    main="net3 82 vs 92", ...)
  print(c)
}

plot.comp13 <- function() {
  zl <- c(0.5, 1.0)
  par(mfrow=c(2,3), mar=c(1,1,1,1))
  plot.comp1(as.image=TRUE, zlim=zl)
  plot.comp3(as.image=TRUE, zlim=zl)
}

assign.type <- function(f) {
  if(f["time"] != "NA")
    return(factor("TS", levels=c("WT","SS","TS")))
  else if(f["pert"] == "NA" && f["deleted"] == "NA" && f["overexpressed"] == "NA")
    return(factor("WT", levels=c("WT","SS","TS")))
  else
    return(factor("SS", levels=c("WT","SS","TS")))
}

# annotates DREAM5 data using known column naming convention
# adds "type" and "group" columns for condition type and group within type
get.features <- function(cnames) {
  if(is.data.frame(cnames) || is.matrix(cnames))
    cnames <- colnames(cnames)
  mat <- t(sapply(strsplit(cnames, "_"), rbind, simplify=T))
  feat <- as.data.frame(mat)
  colnames(feat) <- c("exp","pert","pert.label","treatment","deleted","overexpressed","time","rep")
  type <- apply(feat, 1, assign.type)
  suppressWarnings(feat[,7] <- as.numeric(mat[,7]))
  # assign sequential group numbers to wt and ss types
  group <- array(0, length(type))
  wt.idx <- which(type=="WT")
  group[wt.idx] <- 1:length(wt.idx)
  ss.idx <- which(type=="SS")
  group[ss.idx] <- 1:length(ss.idx)
  lastExp <- ""
  gid <- 0
  tsInds <- which(type=="TS")
  if(length(tsInds) > 0){	
		for (r in 1:length(tsInds)) {
    	if( (lastExp != feat[tsInds[r],"exp"]) || (feat[tsInds[r-1],"time"] > feat[tsInds[r],"time"]))
   	   gid <- gid + 1
   	 group[tsInds[r]] <- gid
   	 lastExp <- feat[tsInds[r],"exp"]
  	}
	}
  return(cbind(feat, type, group))
}

# splits columns of n into matrices by type. types are not mutually exclusive
# types:
#   ts   = all time points
#   flts = first and last only of the time series
#   wto  = wt only
#   wtt  = wt plus all 0 time points
#   pto  = perturbation, all steady-state perturbations (no time series)
#   ptt  = perturbation, all perturbations including time series
#   gno  = genetic, all KO/KD/overexpression minus time series and perturbations
#   gnt  = genetic, all KO/KD/overexpression including time series and perturbations
split.by.type <- function(m) {
  f <- get.features(colnames(m))
  ts <- m[,which(f$type == "TS")]
  wto <- m[,which(f$type == "WT")]
  wtt <- m[,which(f$type == "WT" | (!is.na(f$time) & f$time == 0))]
  pto <- m[,which(f$pert != "NA" & f$type != "TS")]
  ptt <- m[,which(f$pert != "NA")]
  gno <- m[,which((f$deleted != "NA" | f$overexpressed != "NA") & f$pert == "NA" & f$type != "TS")]
  gnt <- m[,which((f$deleted != "NA" | f$overexpressed != "NA") & f$pert == "NA")]
  
  #getting first and last of time-series
  group <- f[which(f[,"type"]=="TS"),"group"]
  isFirstLast <- c(1)
  for(i in 2:length(group)){
    if( i == length(group))
      isFirstLast <- c(isFirstLast,i)
    else if( (group[i-1] != group[i]) || (group[i+1] != group[i]))
      isFirstLast <- c(isFirstLast,i) 
  }
  #f[which(f[,"type"]=="TS")[isFirstLast],] #just for error checking
  flts <- m[,which(f[,"type"]=="TS")[isFirstLast]]
  
  return(list(ts=ts, flts=flts, wto=wto, wtt=wtt, pto=pto, ptt=ptt, gno=gno, gnt=gnt))
}

# finds putative KO genes
find.pko <- function(sd.set, co.set, co.q=.01, top=FALSE, abs=FALSE, graphics=FALSE) {
  # sd and median of sd.set
  sd1 <- sd(t(sd.set))
  med1 <- apply(sd.set, 1, median)
  # calculate z-score on co.set using sd.set sd and med
  z <- apply(co.set, 2, function(x) (x - med1) / sd1)
  if(abs)
    z <- abs(z)
  # find a z-score cutoff
  co <- quantile(z, co.q)
  # get the counts of genes meeting cutoff by condition
  if(top)
    y <- t(apply(z, 2, function(x) length(x[x>=co])))
  else
    y <- t(apply(z, 2, function(x) length(x[x<=co])))
  # optionally show the results
  if(graphics) {
    hist(y, nclass=max(y), main=paste("putative genes, co=", co, " (", co.q, ")", str=""))
    lines(c(5,5), c(0,120), col="red")
  }
  # find conditions with under-/over-expressed genes
  a <- colnames(y)[which(y>0 & y<5, arr.ind=TRUE)[,2]]
  z3 <- z[,a]
  conditions <- colnames(z3)
  # find genes over/under cutoff
  if(top)
    n <- which(z3>=co, arr.ind=T)
  else
    n <- which(z3<=co, arr.ind=T)
  if(length(n) > 0) {
    u <- unique(sort(n[,1]))
    genes <- rownames(z3)[u]
  } else {
    genes <- character()
  }
  return(list(n=n,genes=genes,conditions=conditions))
}

show.pko <- function(ko.in) {
  n <- ko.in$n
  u <- unique(sort(n[,1]))
  cond <- ko.in$conditions
  m <- matrix(0, length(u), length(cond))
  rownames(m) <- as.character(u)
  colnames(m) <- cond
  n2 <- n
  n2[,1] <- as.character(n[,1])
  n2[,2] <- cond[n[,2]]
  for (i in 1:nrow(n2))
    m[n2[i,1],n2[i,2]] <- T
  par(mar=c(3,1,1,12))
  image(m, zlim=c(0,1), axes=FALSE)
  axis(RIGHT<-4, at=seq(0,1,1/(ncol(m)-1)), labels=colnames(m), las=HORIZONTAL<-1, cex.axis=0.75)
  axis(BOTTOM<-1, at=seq(0,1,1/(nrow(m)-1)), labels=rownames(m), las=VERTICAL<-3, cex.axis=0.8)
}

get.samr.raw <- function(g1, g2) {
  m <- cbind(g1, g2)
	if(!is.data.frame(m))
		m <- as.data.frame(m)
  g1.idx <- 1:ncol(g1)
  g2.idx <- ncol(g1) + 1:ncol(g2)
  return(get.samr(m, g1.idx, g2.idx))
}

# Performs basic samr run
# m.in is overall matrix (e.g. net3.full)
# g1 and g2 are indices of groups (length > 1) of replicates to compare
get.samr <- function(m.in, g1.idx, g2.idx) {
  x <- m.in[c(g1.idx,g2.idx)]
  y <- c(rep(1,length(g1.idx)), rep(2,length(g2.idx)))
  data <- list(x=x, y=y, geneid=as.character(1:nrow(x)), genenames=rownames(x), logged2=TRUE)
	capture.output(samr.obj <- samr(data, resp.type="Two class unpaired", nperms=100))
  return(samr.obj)
}

samr.rank <- function(m.in, g1, g2, gene.id, use.abs=TRUE) {
  samr.obj <- get.samr(m.in, g1, g2)
  if(use.abs)
    tt <- sort(abs(samr.obj$tt), decreasing=TRUE)
  else
    tt <- sort(samr.obj$tt, decreasing=FALSE)
  return(list(tt=samr.obj$tt[gene.id], rank=which(names(tt)==gene.id)))
}

print.pairs.ko <- function(l.in) {
  for (l in l.in)
    cat(l$expr, ": ", l$wt.idx, " / ", l$cond.idx, " : ", l$gene, "\n")
}

print.pairs.pk <- function(l.in) {
  for (l in l.in)
    cat(l$expr, ": ", l$wt.idx, " / ", l$cond.idx, " : ", l$pert, " , ", l$level, " , ", l$treatment, "\n")
}

# finds samr ranks from find.ko.pairs.full output
samr.pairs <- function(m.in, l.in, ...) {
  r <- list()
  for (l in l.in) {
    if(length(grep("^G[0-9]*$", l$gene)) != 1)
      r <- append(r, "D")
    else if(length(l$wt.idx) == 1)
      r <- append(r, "A")
    else if(length(l$cond.idx) == 1)
      r <- append(r, "B")
    else {
      v <- samr.rank(m.in, l$wt.idx, l$cond.idx, l$gene, ...)
      r <- append(r, paste(v$tt, " : ", v$rank))
    }
  }
  for (i in 1:length(l.in)) {
    l <- l.in[[i]]
    cat(l$expr, ": ", l$wt.idx, " / ", l$cond.idx, " : ", l$gene, " : ", r[[i]], "\n")
  }
}

# computes z-scores, using samr or diff.exp method, from output of find.annot.pairs.full
#   m: expression matrix
#   l: list in style of find.annot.pairs.full output
#   tf: list of TFs
# returns z-scores as a matrix, z-score as a list of vectors, gene list, and diff and SAM scores
get.ko.z.scores <- function(m, l, tf, global.wt = NULL) {
# analysis for score threshold:
#  for (i in 1:17) y[i] <- length(which(z$scores.samr[[i]] > 2))
#  hist(y, nclass=length(y))
  if(is.null(global.wt))
    global.wt <- m
  sp.start <- "[-+]"
  sp.split <- "[-.+]"  # separator characters for gene name list
  scores.samr <- list()
  scores.diff <- list()
  gene.count <- 0
  genes <- list()
  mat.count <- 0
  mat <- list()
  for (a in l) {
    # parse gene list
    gn <- strsplit(sub(sp.start, "", a$gene), sp.split)
    gene.count <- gene.count + 1
    genes[[gene.count]] <- gn[[1]]
    # run SAM or z.diff
    #   accumulate z-scores matching cutoff criteria
    if(a$wt.idx[1] == 0){
      wt <- global.wt
		}	else{
      wt <- m[,a$wt.idx]
		}
		ko <- m[,a$cond.idx]
    #print(a$wt.idx)
    #print(a$cond.idx)
    if(length(a$wt.idx) >= 3 && length(a$cond.idx) >= 3) {
      # we never get here if we're using global.wt
      z <- get.samr.raw(wt, ko)$tt
      scores.samr[[length(scores.samr)+1]] <- z
    } else {
      z <- diff.exp.no.rep(m, wt=wt, ko=ko)
      scores.diff[[length(scores.diff)+1]] <- z
    }
    # store z-scores for TFs
    for (g in gn[[1]]) {
      if(length(which(tf == g)) > 0) {
        mat.count <- mat.count + 1
        mat[[mat.count]] <- list(gene=g, z=z)
      }
    }
  }
  # now we have mat, a list of lists of TF name and z-scores
  # we need to take the median over all z-scores for each gene for TFs with multiple KOs
  # zf will be this list
  zf <- list()
  for (r in mat) {
    if(is.null(zf[[r$gene]]))
      zf[[r$gene]] <- matrix(0, 0, length(r$z))
    zf[[r$gene]] <- rbind(zf[[r$gene]], r$z)
  }
  # now take the medians
  for (i in 1:length(zf)) {
    zf[[i]] <- apply(zf[[i]], 2, median)
    zf[[i]] <- abs(zf[[i]])
		#colnames(zf[[i]]) <- rownames(m)
    zf[[i]][names(zf)[i]] <- 0     # no self-regulation
  }
  # make a TF/gene matrix of all z-scores
  fmat <- matrix(0, nrow(m), length(tf))
  colnames(fmat) <- tf
  rownames(fmat) <- rownames(m)
  for (n in names(zf))
    fmat[,n] <- zf[[n]]
  # apply thresholds and cutoff
  return(list(z.mat=fmat, z.list.full=zf, genes=genes, scores.samr=scores.samr, scores.diff=scores.diff))
}

# computes z-scores, using samr or diff.exp method, from output of find.putative.pairs.full
#   m: expression matrix
#   l: list in style of find.annot.pairs.full output
#   min.z.samr: minimum z-score (absolute value) for inclusion using SAM
#   min.z.diff: minimum z-score (positive only) for inclusion using diff.exp
#   max.count: maximum number of genes to consider (either method)
# returns sparse matrix as list
# OBSOLETE 9/15
#get.putative.z.scores <- function(m, l, min.z.samr=2, min.z.diff=2, max.count=500) {
# analysis for score threshold:
#  for (i in 1:17) y[i] <- length(which(z$scores.samr[[i]] > 2))
#  hist(y, nclass=length(y))
##   sp <- "[-.+]"  # separator characters for gene name list
##   scores.samr <- list()
##   scores.diff <- list()
##   pert.count <- 0
##   perts <- list()
##   mat.count <- 0
##   mat <- list()
##   for (a in l) {
##     # parse gene list
##     gn <- strsplit(sub(sp, "", a$pert), sp)
##     pert.count <- pert.count + 1
##     perts[[pert.count]] <- gn[[1]]
##     # run SAM or z.diff
##     #   accumulate z-scores matching cutoff criteria
##     if(length(a$wt.idx) >= 3 && length(a$cond.idx) >= 3) {
##       z <- abs(get.samr(m, a$wt.idx, a$cond.idx)$tt)
##       scores.samr[[length(scores.samr)+1]] <- z
##       z <- z[z >= min.z.samr]
##     } else {
##       z <- diff.exp.no.rep(m, a$wt.idx, a$cond.idx)
##       scores.diff[[length(scores.diff)+1]] <- z
##       z <- z[z >= min.z.diff]
##     }
##     # store z-scores
##     ## for(g in gn[[1]]) {
##     ##   if(length(which(tf == g)) > 0) {
##     ##     mat.count <- mat.count + 1
##     ##     mat[[mat.count]] <- list(gene=g, z=z[1:min(max.count, length(z))])
##     ##   }
##     ## }
##   }
##   return(list(scores.samr=scores.samr, scores.diff=scores.diff))
## }

# converts the z-score lists from the $mat element of get.ko.z.scores output
# to a list of z-score vectors by TF - the format required for knockOutConfList
#   koz: full output of get.ko.z.scores, contains $mat element
#   m: expression matrix (used to get number of genes and names)
#   tf: transcription factors (used to get number of TFs and names)
# OBSOLETE 9/15
z.list.to.ko.conf <- function(koz, m, tf.in) {
  # first build a hashtable, more or less
  # list of arrays of z-scores for each TF/gene pair that has a value, index by key "TF_gene"
  tf.gene.list <- list()
  for (a in koz$mat) {
    t <- a$gene
    z <- a$z
    for (gn in names(z)) {
      ix <- paste(t, "_", gn, sep="")   # hash key
      if(is.null(tf.gene.list[[ix]]))
        tf.gene.list[[ix]] <- numeric()
      tf.gene.list[[ix]][length(tf.gene.list[[ix]])+1] <- z[gn]
    }
  }
  # now make the matrix from the value lists above
  z <- matrix(0, nrow(m), length(tf.in))
  rownames(z) <- rownames(m)
  colnames(z) <- tf.in
  for (n in names(tf.gene.list)) {
    # get TF & gene from name
    f <- strsplit(n, "_")
    t <- f[[1]][1]
    g <- f[[1]][2]
    z[g,t] <- sum(tf.gene.list[[n]]) / sqrt(length(tf.gene.list[[n]]))
  }
  # finally make a list for each non-zero TF containing its non-zero gene z-scores
  kl <- list()
  for (t in names(which(apply(z, 2, sum) > 0)))
    kl[[t]] <- z[which(z[,t] > 0), t]
  return(list(ko.list=kl, z=z, tf.gene.list=tf.gene.list))
}

# find pairs of KO/WT experiments
# OBSOLETE
## find.ko.pairs.avg <- function(m) {
##   feat <- get.features(m)
##   wt <- NULL
##   for (i in 1:nrow(feat)) {
##     f <- feat[i,]
##     if(f["type"] == "WT") {
##       wt <- f
##     } else if(!is.null(wt) && f["deleted"] != "NA" && f["pert"] == "NA" && f["overexpressed"] == "NA") {
##       if(f["exp"] == wt["exp"])
##         cat(wt[["exp"]], " : ", rownames(wt)[[1]], " ", rownames(f)[[1]], "\n")
##     }
##   }
## }

# finds all "deleted" and genes along with their local WT if any
# those without local WT are not reported
#   m: DREAM5 expression matrix
# TODO check wt-finding code (use get avg.wt or so?)
find.annot.pairs.full <- function(m.avg, m.full, inc.nowt=T) {
  f.avg <- get.features(m.avg)
  # first find list of first/last time series points
  group <- f.avg[which(f.avg[,"type"]=="TS"),"group"]
  isFirstLast <- c(1)
  for(i in 2:length(group)){
    if( i == length(group))
      isFirstLast <- c(isFirstLast,i)
    else if( (group[i-1] != group[i]) || (group[i+1] != group[i]))
      isFirstLast <- c(isFirstLast,i) 
  }
  ifl.avg <- as.numeric(rownames(f.avg[which(f.avg[,"type"]=="TS"),][isFirstLast,]))
  # now convert from avg idx to full idx
  f.full <- get.features(m.full)
  ifl.full <- unlist(sapply(ifl.avg, function(i) avg.to.full(f.avg, f.full, i)))
  
  # now get on with it
  l <- list()
  lc <- 0
  wt.start <- NULL
  for (i in 1:nrow(f.full)) {
    f <- f.full[i,]
    if(f["type"] == "WT") {
      if(f["rep"] == 1) {  # first wt in set of replicates
        wt.start <- f
        wt.idx <- i
      } else {  # following wt in set of replicates
        wt.idx <- c(wt.idx, i)
      }
    } else if((f["type"] == "SS" || i %in% ifl.full) && !is.null(wt.start) && f["deleted"] != "NA") {
      if(f["rep"] == 1)
        ko.idx <- i
      else
        ko.idx <- c(ko.idx, i)
      if(i == nrow(f.full) || f.full[i+1,"rep"] == 1) {   # last in set of KOs
        if(wt.start["exp"] == f["exp"] || inc.nowt) {     # have local wt, or use global wt flag is true
          if(wt.start["exp"] != f["exp"])
            wt.idx <- 0
          l.ko <- paste("-", as.character(f[["deleted"]]), sep="")
          if(l.ko == "-NA") l.ko <- ""
          lc <- lc + 1
          l[[lc]] <- list(expr=as.character(f[["exp"]]), wt.idx=wt.idx, cond.idx=ko.idx, gene=paste(l.ko, sep=""))
        }
      }
    }
  }
  return(l)
}

# finds all perturbed genes along with their local WT if any
# those without local WT are not reported
#   m: DREAM5 expression matrix
# OBSOLETE 9/15
## find.putative.pairs.full <- function(m) {
##   l <- list()
##   lc <- 0
##   feat <- get.features(m)
##   wt.start <- NULL
##   for (i in 1:nrow(feat)) {
##     f <- feat[i,]
##     if(f["type"] == "WT") {
##       if(f["rep"] == 1) {  # first wt in set of replicates
##         wt.start <- f
##         wt.idx <- i
##       } else {  # following wt in set of replicates
##         wt.idx <- c(wt.idx, i)
##       }
##     } else if(f["type"] == "SS" && !is.null(wt.start) && f["pert"] != "NA" &&
##         f["deleted"] == "NA" && f["overexpressed"] == "NA") {
##       if(f["rep"] == 1)
##         pk.idx <- i
##       else
##         pk.idx <- c(pk.idx, i)
##       if(wt.start["exp"] == f["exp"] && (i == nrow(feat) || feat[i+1,"rep"] == 1)) {   # last in set of KOs
##         lc <- lc + 1
##         l[[lc]] <- list(expr=as.character(wt.start[["exp"]]), wt.idx=wt.idx, cond.idx=pk.idx, pert=as.character(f[["pert"]]),
##             level=as.character(f[["pert.label"]]), treatment=as.character(f[["treatment"]]))
##       }
##     }
##   }
##   return(l)
## }

# computes z-score of differential expression between wt and KO using log difference and log fold-change
#   net: expression matrix
#   wt.idx: indices of wild-type columns in net
#   ko.idx: knockout column indices
#   logged2: (logical=TRUE) has data been log2 normalized?
#   do.plot: (logical=FALSE) should the final z-scores be plotted?
# TODO change mean to median
diff.exp.no.rep <- function(net, wt.idx=NULL, ko.idx=NULL, wt=NULL, ko=NULL, logged2=TRUE, do.plot=FALSE) {
  if((is.null(wt) && is.null(wt.idx)) || (is.null(ko) && is.null(ko.idx)))
    stop("must provide both wt/wt.idx and ko/ko.idx")
  if(is.null(wt))
    wt <- net[,wt.idx]
  if(is.null(ko))
    ko <- net[,ko.idx]
  # take means if wt or ko has more than one condition
  if(!is.null(ncol(wt)) && ncol(wt) > 1) {
    if(logged2)
      wt <- 2^wt
    wt <- apply(wt, 1, mean)
    if(logged2)
      wt <- log2(wt)
  }
  if(!is.null(ncol(ko)) && ncol(ko) > 1) {
    if(logged2)
      ko <- 2^ko
    ko <- apply(ko, 1, mean)
    if(logged2)
      ko <- log2(ko)
  }
  # log2 transform if it isn't done already
  if(!logged2) {
    wt <- log2(wt)
    ko <- log2(ko)
  }
  # z-score of log absolute difference in expression
  abs.diff <- log2(abs(2^wt - 2^ko) + .01)
  z.abs <- (abs.diff - median(abs.diff)) / sd(abs.diff)
  # z-score of absolute fold change
  fc.diff <- ko - wt
  z.fc <- (fc.diff - median(fc.diff)) / sd(fc.diff)
  # Stouffer combination of z-scores
  z.comb <- (z.abs + abs(z.fc)) / sqrt(2)
  names(z.comb) <- rownames(net)
  z.comb[which(z.comb < 0)] <- 0
  ix.net <- which(z.fc < 0)
  z.comb[ix.net] = -z.comb[ix.net]
  if(do.plot) {
    ix <- sort(z.comb, index.return=TRUE,decreasing=T)$ix
    plot(wt[ix],ko[ix])
    points(wt[ix[1:100]],ko[ix[1:100]],col="red",pch=3)
  }
  return(z.comb)
}

# computes fraction of common genes in the top <count> between samr and diff.exp z-score methods
#   net: expression data
#   gwt: indices of wild-type conditions
#   gko: indices of knockout (or other) conditions
compare.samr.diff <- function(net, gwt, gko, count=100) {
  so <- get.samr(net, gwt, gko)
  z.samr <- so$tt
  z.diff <- diff.exp.no.rep(net, gwt, ,gko)
  z.samr.ix <- sort(abs(z.samr), decreasing=T, index.return=T)$ix
  z.diff.ix <- sort(z.diff, decreasing=T, index.return=T)$ix
  i <- intersect(z.samr.ix[1:count], z.diff.ix[1:count])
  return(length(i)/count)
}

get.ranges <- function(l, offset=0) {
  x <- numeric()
  for (a in l)
    x <- c(x, range(a+offset))
  lo <- floor(range(x)[1])
  hi <- ceiling(range(x)[2])
  return(c(lo, hi)-offset)
}

get.hist.max <- function(l, b) {
  x <- numeric()
  for (a in l)
    x <- c(x, max(hist(a, breaks=b, plot=F)$counts))
  return(max(x))
}

plot.samr.diff.scores <- function(m, type, output.prefix=NULL, ...) {
  x <- get.ranges(m, offset=.5)
  b <- x[1]:x[2]
  y <- get.hist.max(m, b)
  for (i in 1:length(m)) {
    if(!is.null(output.prefix))
      png(file=paste(output.prefix, "_", i, ".png", sep=""), ...)
    hist(m[[i]], breaks=b, main=paste(type, i), ylim=c(0, y))
    if(!is.null(output.prefix))
      dev.off()
  }
  return(x)
}

# returns the rows of full corresponding to the given row of avg
avg.to.full <- function(f.avg, f.full, idx.avg) {
  # sum of # reps from 0:idx.avg-1
  l <- 1+sum(as.integer(as.character(f.avg$rep[0:(idx.avg-1)])))
  i <- l+1
  while(i <= nrow(f.full) && f.full[i,"rep"] != 1) {
    l <- c(l, i)
    i <- i + 1
  }
  return(l)
  ## return(which(f.full[,"exp"] == f.full[l,"exp"] & f.full[,"pert"] == f.full[l,"pert"] & f.full[,"pert.label"] == f.full[l,"pert.label"]
  ##                    & f.full[,"treatment"] == f.full[l,"treatment"] & f.full[,"deleted"] == f.full[l,"deleted"]
  ##                    & f.full[,"overexpressed"] == f.full[l,"overexpressed"]
  ##                    & (f.full[,"time"] == f.full[l,"time"] | is.na(f.full[,"time"]) & is.na(f.full[l,"time"]))))
}

# returns the WT(s) corresponding to the given row of avg; returned value may include idx.avg
get.avg.wt <- function(f.avg, idx.avg) {
  return(which(f.avg[,"exp"] == f.avg[idx.avg, "exp"] & f.avg[,"type"] == "WT"))
}
