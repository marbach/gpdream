
makeDataGP <- function(d.path, f.path, tfs.path){
	#loading data
	# tfs <- as.character(read.table(tfs.path)[[1]])
	
	data.full  <- t(read.table(d.path))
	rnames <- as.vector(data.full[, 1]); data.full <- data.full[, -1]
	#removing the quotes from the values, ie. making it numeric
	data.full <- matrix(as.numeric(data.full), nrow(data.full), ncol(data.full))
	rownames(data.full) <- rnames
	
	if(!is.na(tfs.path)){
		tfs <- as.character(read.table(tfs.path)[[1]])
	} else {
		tfs <- rnames
	}	

	## AM switch ##
	#getting the features table, which we need to make colnames
	feat <- c()
	if(is.na(f.path)){
		cat("\nNOTE: No meta information was given as input. 
				 Defaulting to using no time series and no knockout data.")
		#make default features table
		feat <- matrix(NA, ncol(data.full), 8)
		feat[, 1] <- 1:ncol(data.full)
		feat[, 2] <- paste("P", 1:ncol(data.full) , sep = "")
		feat[, 8] <- rep(1, ncol(data.full))
	}else{
		feat <- read.table(f.path)
	}
	##################
	
	
	
	#adding colnames
	n <- apply(feat, 1, function(i) paste(i, collapse = "_"))
	n <- unlist(lapply(n, function(i) gsub(" ", "", i)))
	n <- unlist(lapply(n, function(i) gsub(",", ".", i)))
	colnames(data.full) = n
	
	#averaging replicates
	data.avg <- combine.replicates(data.full, feat)
	
	## AM fix this 60 unit conversion
	l <- make.ratios(data.avg, 1)#mult.t)
	ratios <- l$ratios
	if(is.na(f.path)){
		# make one element of colMap
		 o <- list(isTs = F,
		        is1stLast = factor("e", levels=c("e","f","m","l")),
		        prevCol= as.logical("NA"),
		        del.t = as.logical("NA"),
		        condName = character())
		colMap <- list()
		for(j in 1:ncol(ratios)){
			colMap[[j]] <- o
			colMap[[j]][["condName"]] <- colnames(ratios)[j]
		}
	} else {
		colMap <-  make.colMap(ratios, l$features)
	}
	clusterStack <- make.clusterStack(ratios)
	
	#making the ko filter list
	ko <- tryCatch({make.knockOutConfList(data.avg, data.full, tfs )}, 
			           error = function(err) {
									             #stop("package samr not installed, install samr and its dependencies")
									             return(NULL)
														}
						 )
			         
	knockOutFilterList <- ko$lst
	knockOutCombine <- ko$mat
	
	if(is.null(knockOutFilterList) || is.null(knockOutCombine))
		cat("\nNOTE: Did not detect any knockout expirements.
				 Defaulting to steady state mode.
         To use knockout data make sure a meta file is included and  package SAMR is installed on the server!")
	 
	
	obj <- list()
	obj[["ratios"]] <- ratios
	obj[["cM"]]     <- colMap
	obj[["cS"]]     <- clusterStack
	obj[["tfs"]]    <- tfs
	obj[["koFilt"]] <- knockOutFilterList
	obj[["koComb"]] <- knockOutCombine
	
	return(obj)
}

# makes the inferelator "ratios" matrix from input, e.g. net2.avg
# mult.t is the factor by which to multiply time points such that they are all integers
# experiments are reordered as follows: wt, ss, ts
# experiments are named as follows: {WT,SS,TS}_Ng_Nx[_tNt] where Ng is the sequential number with the wt/ss/ts type,
#   Nx is the original sequential experiment (technically condition) number, and Nt is the integer time point
# returns the ratios matrix and a permuted features matrix to match the row numbers of ratios
make.ratios <- function(m, mult.t) {
  f <- get.features(colnames(m))
	
	mats.list <- list()

  wt.idx <- which(f$type == "WT")
  wt <- m[,wt.idx]
  cn <- character()
  for(i in 1:length(wt.idx))
    cn[i] <- paste("WT_", f[wt.idx[i], "group"], "_", wt.idx[i], sep="")
  if(ncol(wt) != 0)
		colnames(wt) <- cn
 	mats.list[["WT"]] = list()
	mats.list[["WT"]][["idx"]] <- wt.idx
	mats.list[["WT"]][["mat"]] <- wt
	
	ss.idx <- which(f$type == "SS")
  ss <- m[,ss.idx]
  cn <- character()
  for(i in 1:length(ss.idx))
    cn[i] <- paste("SS_", f[ss.idx[i], "group"], "_", ss.idx[i], sep="")
  if(ncol(ss) != 0)
		colnames(ss) <- cn
	mats.list[["SS"]]  <- list()
	mats.list[["SS"]][["idx"]] <- ss.idx
	mats.list[["SS"]][["mat"]] <- ss
	
	ts.idx <- which(f$type == "TS")
  ts <- m[,ts.idx]
  cn <- character()
  for(i in 1:length(ts.idx))
    cn[i] <- paste("TS_", f[ts.idx[i], "group"], "_", ts.idx[i], "_t",
        round(f[ts.idx[i], "time"]*mult.t), sep="")
  if(ncol(ts) != 0)
		colnames(ts) <- cn
	mats.list[["TS"]]  <- list()
	mats.list[["TS"]][["idx"]] <- ts.idx
	mats.list[["TS"]][["mat"]] <- ts
	
	num.idx <-	unlist(lapply(mats.list, function(i) length(i$idx)))
	non.zero.mats <- names(which(num.idx != 0))
	
	ratios   <- NULL
	for(i in non.zero.mats){
		ratios <- cbind(ratios, mats.list[[paste(i, sep = "")]]$mat)
	}
	
	features <- f[c(wt.idx,ss.idx,ts.idx),]
  features["time"] <- round(features["time"]*mult.t)
  return(list(ratios=ratios, features=features))
}

combine.replicates <- function(d, f){
	nd <- matrix(0, nrow(d), 0)
	idx = 1
	cnames = c()
	while (idx  < nrow(f)) {
		#creating a name from features
		#and using the name to see if there are duplicates
		n1 <- paste(f[idx, 1:7], collapse = "_")
		n2 <- paste(f[idx + 1, 1:7], collapse = "_")
				
		if (n1 != n2){
			nd  <- cbind(nd, d[, idx])
			cnames <- c(cnames, colnames(d)[idx])
			
			idx <- idx +  1
			next;
		}
		
		#if we get here then we have duplicates
		dup.idx <- idx
		idx = idx + 1
		while((idx <= nrow(f)) && (f[idx, 8] != 1)){
			dup.idx <- c(dup.idx, idx)
			idx = idx + 1
		}
		cnames <- c(cnames, colnames(d)[max(dup.idx)])
		
		cr <- apply(d[, dup.idx], 1, mean)
		nd <- cbind(nd, cr)
	}
	
	n1 <- paste(f[idx - 1, 1:7], collapse = "_")
	n2 <- paste(f[idx, 1:7], collapse = "_")
	#adding the last row
	if ((n1 != n2) && (idx == ncol(d))){
		nd  <- cbind(nd, d[, idx])
		cnames <- c(cnames, colnames(d)[idx])
		idx <- idx +  1
	}
	
	colnames(nd) <- cnames
	
	return(nd)
}

# makes the inferelator "colMap" list from ratios matrix with DREAM5 column names
#   and permuted features matrix from make.ratios
make.colMap <- function(ratios, features) {
  colMap <- list()
  cnames <- colnames(ratios)
  numTS <- 0
  for (i in 1:ncol(ratios)) {
    o <- list(isTs = F,
        is1stLast = factor("e", levels=c("e","f","m","l")),
        prevCol= as.logical("NA"),
        del.t = 0,
#              del.t = as.logical("NA"),
        condName = character())
    if(features[i,"type"] == "TS") {
      # condition is part of a time series
      o$isTs <- T
      numTS <- numTS + 1
      # which time point and group?
      tp.curr <- features[i, "time"]
      tp.grp <- features[i, "group"]
      # time point first/middle/last parsing logic (at this point we know the current condition is a time point)
      # start new time series when: condition is overall first condition; previous condition was not TS,
      #   was different group, or had later time point than current
      # end current time series when: condition is overall last condition; next condition is not TS,
      #   is different group, or has earlier time point than current
      if(i == 1) tp.prev <- NA else tp.prev <- features[i-1, "time"]
      if(is.na(tp.prev) || tp.curr < tp.prev || tp.grp != features[i-1, "group"]) {
        o$is1stLast <- factor("f", levels=c("e","f","m","l"))
      } else {
        # time point of next condition, or NA to force "end" if we're at the last condition
        if(i == ncol(ratios))
          tp.next <- NA
        else
          tp.next <- features[i+1, "time"]
        # find out if the next time point starts a new series: if so, we're last, otherwise middle
        if(is.na(tp.next) || tp.next < tp.curr || tp.grp != features[i+1, "group"])
          o$is1stLast <- factor("l", levels=c("e","f","m","l"))
        else
          o$is1stLast <- factor("m", levels=c("e","f","m","l"))
        # assign the name of the previous time point (can't get here if i==0)
        o$prevCol <- cnames[i-1]
        # calculate time step from previous point
        o$del.t <- tp.curr - tp.prev
      }
    }
    # in all cases, assign condition name
    o$condName <- cnames[i]
    
    # final assignments
    colMap[[ cnames[i] ]] <- o
  }
  #colMap[["numTS"]] <- numTS
  return(colMap)
}

# essentially a direct copy of DREAM[234] makeClusterStack from readTxtDataUtil.R
make.clusterStack <- function(ratios) {
  object <- list()
  object$cols <- colnames(ratios)
  object$ncols <- length(colnames(ratios))
  object$rows <- character()
  object$nrows <- 1
  object$resid <- as.logical("NA")
  object$k <- integer()
  object$redExp <- numeric()
  
  clusterStack <- list()
  for (i in 1:dim(ratios)[1]) {
    clusterStack[[i]] <- object
    clusterStack[[i]]$rows <- rownames(ratios)[i]
    clusterStack[[i]]$k <- i
    clusterStack[[i]]$redExp <- ratios[i,]
  }
  clusterStack[[(i+1)]] <- i
  return(clusterStack)
}

make.knockOutConfList <- function(net.avg, net.full, tf) {
  l <- find.annot.pairs.full(net.avg, net.full)
  koz <- get.ko.z.scores(net.full, l, tf)
  return(list(lst=koz$z.list.full, mat=koz$z.mat))
}

meta.make.knockOutConfList <- function(net) {
  net.avg <- load.net.avg(net)
  net.full <- load.net.full(net)
  tf <- load.tfs(net)
  ko <- make.knockOutConfList(net.avg, net.full, tf)
  knockOutFilterList <- ko$lst
  knockOutCombine <- ko$mat
  save(knockOutFilterList, file=paste("input/DREAM5/net", net, "/knockOutFilterList.RData", sep=""))
  save(knockOutCombine, file=paste("input/DREAM5/net", net, "/knockOutCombine.RData", sep=""))
}


# create putative knockout list from get.putative.kos and get.putative.annot.kos
# returns a list of z=z-score and pko=putative KO list (for making plots with plot.pkos
make.pkoConfList <- function(net.avg.orig, net.full.orig, tf, percentage, maxLow) {
  # this is awkward: to remove bogus conditions, we have to remove them from both average and full!
  f.avg.orig <- get.features(net.avg.orig)
  f.full.orig <- get.features(net.full.orig)
  bad.idx <- bogus.idx(net.avg.orig)
  bad.idx.full <- unlist(sapply(bad.idx, function(i) avg.to.full(f.avg.orig, f.full.orig, i)))
  net.avg <- net.avg.orig[,-bad.idx]
  net.full <- net.full.orig[,-bad.idx.full]
  f.avg <- get.features(net.avg)
  f.full <- get.features(net.full)
  # now we're ready to proceed with making the pko list
  netSplit <- split.by.type(net.avg)
  m <- as.matrix(cbind(netSplit$flts, netSplit$gno, netSplit$pto, netSplit$wto))
  pko.raw <- get.putative.kos(m, qtlz=1-percentage, qtlabs=percentage, tfNames=tf, maxLow=maxLow, matchWt=T)
  pko <- condense.pko.list(pko.raw[[1]])
  p <- convert.putative(pko, net.avg, f.avg, f.full)
  conf.pko <- get.ko.z.scores(net.full, p, tf)
  return(list(z.list=conf.pko$z.list.full,z.mat=conf.pko$z.mat,pko.raw=pko.raw,pko=pko))
}

meta.make.pkoConfList <- function(net) {
  # cutoff values for each network
  pc.filter  <- c(0.941, 0.855, 0.94, 0.941)
  lo.filter  <- c(14, 16, 24, 34)
  pc.combine <- c(0.947, .906, .928, .932)
  lo.combine <- c(10, 11, 17, 17)
  # used for both runs
  net.avg <- load.net.avg(net)
  net.full <- load.net.full(net)
  tf <- load.tfs(net)
  # filter list
  poz.filter <- make.pkoConfList(net.avg, net.full, tf, pc.filter[net], lo.filter[net])
  pkoFilterList <- poz.filter$z.list
  save(pkoFilterList, file=paste("input/DREAM5/net", net, "/pkoFilterList.RData", sep=""))
  # combination matrix
  poz.combine <- make.pkoConfList(net.avg, net.full, tf, pc.combine[net], lo.combine[net])
  pkoCombine <- poz.combine$z.mat
  save(pkoCombine, file=paste("input/DREAM5/net", net, "/pkoCombine.RData", sep=""))
}

make.input <- function(net, mult.t) {
  data <- load.net.avg(net)
  l <- make.ratios(data, mult.t)
  ratios <- l$ratios
  colMap <-  make.colMap(ratios, l$features)
  clusterStack <- make.clusterStack(ratios)
  tfNames <- load.tfs(net)
  # TODO update for pKO lists
  stop("this must be updated for knockout and pkos")
  ko <- make.knockOutConfList(data, tfNames)
  knockOutConfList <- ko$lst
  knockOutConfMatrix <- ko$mat
  save(ratios, file=paste("input/DREAM5/net", net, "/ratios.RData", sep=""))
  save(colMap, file=paste("input/DREAM5/net", net, "/colMap.RData", sep=""))
  save(clusterStack, file=paste("input/DREAM5/net", net, "/clusterStack.RData", sep=""))
  save(tfNames, file=paste("input/DREAM5/net", net, "/tfNames.RData", sep=""))
  save(knockOutConfList, file=paste("input/DREAM5/net", net, "/knockOutConfList.RData", sep=""))
  save(knockOutConfMatrix, file=paste("input/DREAM5/net", net, "/knockOutConfMatrix.RData", sep=""))
}

make.input.all <- function() {
  cat("network 1\n")
  make.input(1, 1)
  cat("network 2\n")
  make.input(2, 60)
  cat("network 3\n")
  make.input(3, 60)
  cat("network 4\n")
  make.input(4, 60)
}

# convert Alex G list of putative KOs into the same format used by find.annot.pairs.full
#   pko.list is named list: name is condition name; value is pKO gene
convert.putative <- function(pko.list, net.avg, f.avg, f.full) {
  r <- list()
  for (i in 1:length(pko.list)) {
    l <- NULL
    gene <- pko.list[i]
    cond.avg <- which(colnames(net.avg) == names(pko.list[i]))
    cond.full <- avg.to.full(f.avg, f.full, cond.avg)
    wt.avg <- get.avg.wt(f.avg, cond.avg)
    if(!any(wt.avg) || any(intersect(cond.avg, wt.avg)))   # if there's no wt or cond is itself a wt
      wt.full <- 0
    else      # we have a wt that is not the condition
      wt.full <- avg.to.full(f.avg, f.full, wt.avg)
    l <- list(expr=as.character(f.avg[cond.avg,"exp"]), wt.idx=wt.full, cond.idx=cond.full, gene=gene)
    if(!is.null(l))
      r[[length(r)+1]] <- l
  }
  return(r)
}

# convert wt.cond.all to list suitable for get.ko.z.scores input
# broken, do not use
all.zsc <- function(net.avg, net.full) {
  r <- list()
  for (i in 1:length(n)) {
    l <- NULL
    if(n[i] != 0) {
      gene <- "Gxxx"
      wt.full <- avg.to.full(f.avg, f.full, n[i])
      cond.full <- avg.to.full(f.avg, f.full, i)
      l <- list(expr=as.character(f.avg[i,"exp"]), wt.idx=wt.full, cond.idx=cond.full, gene=gene)
    }
    if(!is.null(l))
      r[[length(r)+1]] <- l
  }
  return(r)
}
