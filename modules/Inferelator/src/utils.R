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

# plots  a histogram with values from bg.scores and overlays conf.vals for specific set on the hitograms in red (many par() params)
plot.hits.in.hist <- function(bg.scores,conf.vals,color="red",linewidth=1,main="",yLable="Counts",whr.lgnd="topright",
                              main.title.name="clr scores",xlab.title.name = "clr pseudo z score",ylab.title.name = "Counts",qunt = .75,cx=.7) {
  h <- hist(bg.scores,nclass=1000,main=main.title.name,xlab = xlab.title.name,ylab = ylab.title.name)
  hits.height = quantile(h$counts,qunt)
  zero.cnt <- 0
  for(j in 1:length(conf.vals)){
    lines(c(conf.vals[j],conf.vals[j]),c(0,-hits.height),lwd=linewidth,col="red")
    if(conf.vals[j] == 0) {
      zero.cnt <- zero.cnt + 1
    }
  }
  ix <- sort(conf.vals,decreasing=TRUE,index.return=TRUE)$ix
  legend(whr.lgnd,paste(names(conf.vals[ix]),"=",round(conf.vals[ix],2)),lty=1,lwd=linewidth,col="red",cex=cx)
}

# tw samples unequal varaince t test
# input: xm1,xm2 - mean of each group
#        sd1,sd2 - standard deviations of each group
#        number of observations (equal size for both groups)
# output: zscores
zScoreTwoSampleUnequlVar <- function(xm1,xm2,sd1,sd2,n){
  (xm1-xm2)/sqrt((sd1^2+sd2^2)/n)
}



# calculate the sample standard deviation from the s_j sums
#input: s_j = \sigma_{k=1}^{N} x^j_k
# output sd.
stdBySums <- function(s0,s1,s2) {
  sqrt(((s0*s2)-s1^2)/(s0*(s0-1)))
}

getMedianNetworkFromBootstrapsExcludeZeros <- function(allRes, pipeline = "MCZ.MixCLR.Inf"){
	# 1) get all MCZ.MixCLR.Inf confidence scores matrices
	x <- lapply(allRes, function(i) i[[pipeline]])
	# 2) collapse x list into a 3 dim array
	all.conf.scores <- array(0,c(dim(x[[1]]),length(x)))
	for (i in 1:length(x)){
		all.conf.scores[,,i] <- x[[i]]
	}
	# 3) get median (and stndrd dev) matrix (median confidence score for each putative regulatory interaction)
	median.conf.scores.exclude.zeros <- matrix(0,nrow=dim(all.conf.scores)[1],ncol=dim(all.conf.scores)[2])
	for (i in 1: dim(all.conf.scores)[1]){
		for(j in 1:dim(all.conf.scores)[2]) {
                  non.zeros.ix <- which(all.conf.scores[i,j,]!=0)
                  if(length(non.zeros.ix)>0){
                    median.conf.scores.exclude.zeros[i,j] <- median(all.conf.scores[i,j,][non.zeros.ix])
                  }
		}
	}	
	return(median.conf.scores.exclude.zeros)
}

getFrequencyNetworkFromBootstraps <- function(allRes, pipeline = "MCZ.MixCLR.Inf"){
	# 1) get all MCZ.MixCLR.Inf confidence scores matrices
	x <- lapply(allRes, function(i) i[[pipeline]])
	# 2) collapse x list into a 3 dim array
        frequency.matrix <- matrix(0,nrow=dim(x[[1]])[1],ncol=dim(x[[1]])[2])
	for (i in 1:length(x)){
		frequency.matrix[which(x[[i]]!=0)] <- frequency.matrix[which(x[[i]]!=0)] + 1
	}
        frequency.matrix <- frequency.matrix/length(x)
	return(frequency.matrix)
}

getMedianNetworkFromBootstraps <- function(allRes, pipeline = "MCZ.MixCLR.Inf"){
	# 1) get all MCZ.MixCLR.Inf confidence scores matrices
	x <- lapply(allRes, function(i) i[[pipeline]])
	# 2) collapse x list into a 3 dim array
	all.conf.scores <- array(0,c(dim(x[[1]]),length(x)))
	for (i in 1:length(x)){
		all.conf.scores[,,i] <- x[[i]]
	}
	# 3) get median (and stndrd dev) matrix (median confidence score for each putative regulatory interaction)
	median.conf.scores <- matrix(0,nrow=dim(all.conf.scores)[1],ncol=dim(all.conf.scores)[2])
	for (i in 1: dim(all.conf.scores)[1]){
		for(j in 1:dim(all.conf.scores)[2]) {
			median.conf.scores[i,j] <- median(all.conf.scores[i,j,])
		}
	}	
	return(median.conf.scores)
}

getMeanNetworkFromBootstraps <- function(allRes, pipeline = "MCZ.MixCLR.Inf"){
	# 1) get all MCZ.MixCLR.Inf confidence scores matrices
	x <- lapply(allRes, function(i) i[[pipeline]])
	# 2) collapse x list into a 3 dim array
	all.conf.scores <- array(0,c(dim(x[[1]]),length(x)))
	for (i in 1:length(x)){
		all.conf.scores[,,i] <- x[[i]]
	}
	# 3) get median (and stndrd dev) matrix (median confidence score for each putative regulatory interaction)
	mean.conf.scores <- matrix(0,nrow=dim(all.conf.scores)[1],ncol=dim(all.conf.scores)[2])
	for (i in 1: dim(all.conf.scores)[1]){
		for(j in 1:dim(all.conf.scores)[2]) {
			mean.conf.scores[i,j] <- mean(all.conf.scores[i,j,])
		}
	}	
	return(mean.conf.scores)
}

makePredictions <- function(Z,pred.mat.bias,beta.mat,dream_ko,dream_wt, pred.mat.lnet, inCut=75,inDblKo, initChoice="koMean"){
	biasVec <- matrix(0,nrow(Z),1)
	unqModels <- unique(pred.mat.bias[,1])
	for( i in 1:length(unqModels)){
		biasVec[unqModels[i],1] <- pred.mat.bias[ which(pred.mat.bias[,1] == unqModels[i])[1], 3]
	}
	modelWeights <- cbind(biasVec,beta.mat)
	maxVec <- apply(dream_ko,2,max)
	predWeights <- apply(pred.mat.lnet,1,sum)
	wtWeights <- 1 - predWeights
	
	dblKo <- inDblKo #INPUT$general$dbl_ko
	dblKoPreds <- matrix(, nrow(dblKo),ncol(dream_ko))
	bestCutOffVal <- quantile(Z, inCut/100)
	
	rmsdKo1 <- c()
	rmsdKo2 <- c()
	rmsdWt <- c()
	rmsdInit <- c()
	
	for( dblInd in 1:nrow(dblKo)){
		curKos <- dblKo[dblInd,]
		curWt <- c()
		curZ <- apply( Z[,curKos],1,max )
		if( initChoice == "origWt"){
			curWt <- dream_wt
		}else if(initChoice == "koMean"){
			curWt <- apply( dream_ko[,inDblKo[dblInd,]], 1, median)
		}else if(initChoice == "combine"){
			curWt <- apply( dream_ko[,curKos], 1, median)
			if( any(curZ > bestCutOffVal) ){
				toChange <- which( curZ > bestCutOffVal )
				zOne <- Z[toChange,curKos[1]]
				zTwo <- Z[toChange,curKos[2]]
				curWt[toChange] <- (zOne*dream_ko[ toChange,curKos[1] ] + zTwo*dream_ko[ toChange,curKos[2] ])/(zOne + zTwo)
				#MAKE SURE THIS WORKS...combine based on Z-score, and divide by sum of z-scores
				#curWt[toChange] <- apply( dream_ko[toChange,curKos],1,max) 
				#				curWt[toChange] <- dream_ko[toChange,curKos[1]]*
			}	
		}
		
		curInitCond <- curWt#dream_wt    #setting initial condition to wild type
		curInitCond[ curKos ] <- 0 #knocking out two genes
		curInitCond <- c(1, curInitCond) #here initCond is just a column vector, we add a one to it to 
																		 #incorporate bias term into our model
		
		#idea below(commented out lines) was to make one prediction using column i of dream_ko
		#and column j (i,j corresponding to the knockouts), and then combine the values
		#curInitCond1 <- c(1, dream_ko[,curKos[1]])
		#curInitCond2 <- c(1, dream_ko[,curKos[2]])
		#curPrediction1 <- modelWeights%*%curInitCond1
		#curPrediction2 <- modelWeights%*%curInitCond2
		#curPrediction <- apply(cbind(curPrediction1,curPrediction2),1,median)#modelWeights%*%curInitCond
		
		curPrediction <- modelWeights%*%curInitCond
		
		#now we squash between zero and max we see
		curPrediction[ curPrediction > maxVec ] <- maxVec[ curPrediction > maxVec ]
		curPrediction[ curPrediction < 0] <- 0
		
		#now we do our weighting
		curPrediction <- curPrediction*predWeights + curWt*wtWeights
		
		#now filter
		curIndsToReset <- which( curZ < bestCutOffVal )
		curPrediction[ curIndsToReset ] <- dream_wt[ curIndsToReset ]
		
		#now set the predicted values for the genes we just knocked out to zero
		curPrediction[ curKos ] <- 0
		#calc rmsds
		rmsdKo1 <- c( rmsdKo1,( sum((dream_ko[ ,curKos[1]] - curPrediction)^2)/length(curPrediction)))
		rmsdKo2 <- c( rmsdKo2,( sum((dream_ko[ ,curKos[2]] - curPrediction)^2)/length(curPrediction)))
		rmsdWt <- c(rmsdWt, ( sum((curWt - curPrediction)^2)/length(curPrediction)))
		rmsdInit <- c(rmsdInit, ( sum((curInitCond[-1] - curPrediction)^2)/length(curPrediction)))
		#add to matrix
		dblKoPreds[ dblInd, ] <- curPrediction
	}
	return(dblKoPreds)
	#dblPredFName <-paste(chalName,"_Bonneau_InSilico_Size",numGenesInNet,"_",numNet,"dualknockouts.txt",sep="")
	#write.table(dblKoPreds,sep="\t",file=dblPredFName,row.names=F,col.names=F,quote=F)	
}

#function to test how good our double knockout predictions are
#--INPUT--#
#1) predictions  - matrix with the rows being the particular doublkockounts, and columns correspond to values of the other genes
#2) goldStandard - the true values of network states for all genes given particular knockout - rows correspond to rows of predictions matrix
#3) ko_inds      - matrix of indices of which genes are doubly knockedout...columns indiciate which gene is knocked out
#4) dream_ko     - matrix of single knockouts
#5) wt_meas      - what we use for wt values, could be the given wt values, or the ones we calculate (ie. median over all conditions)
#6) whichNull    - determines what we use as our null prediction...choices:
#                   "origWt" - use wt values(whichever way we calculate them) as null predictions
#                   "meanKo" - take columns i,j of dream_ko (where i and j are the genes being knocked out), and use the mean of them as null prediction  
#--OUTPUT--#
#one list containing
#1) rmsdPred - error(in terms of RMSD) of our prediction to the true answer
#2) rmsdNull - error(in terms of RMSD) of the null prediction(in whichever way we calculated it) to the true answer
calcDblKoRmsd <- function( predictions, goldStandard, ko_inds, dream_ko, wt_meas, whichNull = "meanKo"){
	#calculates RMSD between true dblk_ko and predicted dbl_ko
	#also RMSD between true dlb_ko, and using WT as dbl_ko pred
	rmsdPred <- matrix(0,nrow(predictions),1)
	rmsdNull <- matrix(0, nrow(predictions),1)
	for(i in 1:nrow(predictions)){
		if( (predictions[i,ko_inds[i,]] != 0) || (goldStandard[i,ko_inds[i,]] != 0) ){
			cat( "KOINDS NOT EQUAL ZERO AT PRED IND",i,"\n")
		}
		curNull <- c()
		#curWt
		if(whichNull == "origWt"){
			curNull <- wt_meas
		}else if(whichNull == "meanKo"){
			curNull <- apply( dream_ko[,ko_inds[i,]], 1, mean)#wt_meas
		}
		curNull[ ko_inds[i,] ] <- 0
		rmsdPred[i,1] <- sqrt( sum((predictions[i,] - goldStandard[i,])^2)/ncol(predictions) )
		rmsdNull[i,1] <- sqrt( sum((curNull - goldStandard[i,])^2)/ncol(predictions) )
	}
	names(rmsdPred) <- paste(ko_inds[,1],"_",ko_inds[,2],sep="")
	names(rmsdNull) <- paste(ko_inds[,1],"_",ko_inds[,2],sep="")
	return( list(as.vector(rmsdPred), as.vector(rmsdNull)) )
}


# description:
calcFoldChange <- function( perturbedData, wildType, epsZero ){
	epsZero <- epsZero #.0050
	foldChange <- matrix( 0, nrow(perturbedData), nrow(perturbedData))
	rownames( foldChange ) <- rownames( perturbedData )
	colnames( foldChange ) <- colnames( perturbedData )
	for( i in 1:nrow(perturbedData)){
	  foldChange[,i] <- t((perturbedData[,i] - wildType)/(wildType + epsZero))	
		#rownames(foldChange) <- rownames( perturbedData )
	}
	return( foldChange )
}
# description:
#otherData is the other dataset which is to be used as an estimate of sd
calcZscores <- function( inData, otherData, wtData, sigmaZero, numRem, calcSDOverAll = FALSE ){ #sigma-zero is a correction term
 	#eventually want a psuedocount
	sigmaZero <- sigmaZero#.00375 #can also try .005 or .0025 #is most useful for dream3, not so much for dream4
	
	zScores <- matrix(,nrow(inData), ncol(inData))
	rownames( zScores ) <- rownames( inData )
	colnames( zScores ) <- colnames( inData )
	
	allSds <- c()
	
	numRem <- numRem#c() #number of outliers to remove in each direction when calculating sd
#	if( nrow(inData) == 10){
#		numRem <- 0 #prev value was 1
#	}else{
#		numRem <- 0 #prev val was 5
#	}
	
	for( i in 1:nrow(inData)){
		highest <- sort(inData[i,], decreasing=T, index.return=T)$ix
		highest <- highest[ -which( highest == i ) ]
		lowest <- sort(inData[i,], decreasing=F, index.return=T)$ix
		lowest <- lowest[ -which( lowest == i ) ]
		self <- i
		toRemove <- unique( c(lowest[1:numRem], highest[1:numRem], self))
		
		if(calcSDOverAll){
			otherHighest <- sort(otherData[i,], decreasing=T, index.return=T)$ix
			otherHighest <- otherHighest[ -which( otherHighest == i ) ]
			otherLowest <- sort(otherData[i,], decreasing=F, index.return=T)$ix
			otherLowest <- otherLowest[ -which( otherLowest == i ) ]
			otherSelf <- i
			otherToRemove <- unique( c(otherLowest[1:numRem], otherHighest[1:numRem], otherSelf))
		}
		
		if(calcSDOverAll){
			allSds[i] <- sd( c(inData[i,-toRemove], otherData[i,-otherToRemove]) )
		}else{
			allSds[i] <- sd( inData[i,-toRemove] )
		}
	}
	
	for(i in 1:ncol(inData)){
		#cat("using sigma zero! \n")
		zScores[,i] <- (inData[,i] - wtData)/(allSds + sigmaZero)
	}
	rtrn.obj = list()
	rtrn.obj[[1]] = zScores
	rtrn.obj[[2]] = allSds
	return( rtrn.obj )
}

splitDreamDataByType <- function(inData){
	ratios <- c()
	if(is.null(inData)){
		source("scripts/init_util.R")
		dataSet <- let_usr_choose_dataset()
		allData <- get_usr_chosen_dataset(dataSet)
		ratios <- allData[[1]]	
		gs <- allData[[5]]
		dreamPred <- allData[[6]]
	}else{
		ratios <- inData
	}

	#return( list(dream4_ts,dream4_ko,dream4_kd,dream4_wt) )
	wtCols <- grep( "wt", colnames(ratios) )
	dream3_wt <-  ratios[,wtCols[1]]
	
	tsCols <- grep( "TS", colnames(ratios)) 
	dream3_ts <- ratios[ ,tsCols]
	
	koCols <- grep( "-/-", colnames(ratios) )
	dream3_ko <- ratios[,koCols]
	
	ssCols <- grep("SS",colnames(ratios))
	
	#nonKd <- c(wtCols,koCols,tsCols)
	if(any(ssCols)){
		kdCols <- c(1:ncol(ratios))[ -c(wtCols,koCols,tsCols,ssCols) ]
	}else{
		kdCols <- c(1:ncol(ratios))[ -c(wtCols,koCols,tsCols) ]
	}
	dream3_kd <- ratios[ ,kdCols]
	
	return ( list( dream3_ts, dream3_ko, dream3_kd, dream3_wt))	
}


# nicely print all parameters of a inf run
print_params <- function(PARAMS,file="") {
	cat("",file=file)
	for (i in 1:length(PARAMS)) {
		if(length(PARAMS[[i]]) > 0){
			cat(names(PARAMS[i]), ":\n",file=file,append = TRUE)
			for (j in 1:length(PARAMS[[i]]))
				cat("\t",names(PARAMS[[i]][j])," -> ", PARAMS[[i]][[j]],"\n",file=file,append = TRUE)
		} else {
			cat(names(PARAMS[i]),": empty list","\n",file=file,append = TRUE)
		}
	}
}

# Add a confidance measure for each non-zero wighted predictor
# Measure is how much does each predictor explain of the total variance (as part of the whole model)

add_weight_beta <- function(bL,model_errors,n,pS,pD,col=4,col_name = "prd_xpln_var"){
#	wBetaList = betaList
	for (j in 1:length(bL)) {
		if(length(bL[[j]]) > 0) {
			w = matrix(0,n,ncol=(pS+pD))
			x = unsparse(bL[[j]],w)
			sum_beta = apply(abs(x),1,sum)
			for (k in 1:n) {
				if(sum_beta[k] > 0){
					w[k,] = x[k,]/sum_beta[k]*(1-model_errors[k,j])
				}
			}
			bL[[j]] = cbind(bL[[j]],0)
			x= make_sparse2(w)
#		y = sort(x[,1],index.return=TRUE)
			for (i in 1:nrow(x)){
				r=which(bL[[j]][,1]==x[i,1] & bL[[j]][,2]==x[i,2])
				bL[[j]][r,col] = abs(x[i,3])
			}
			colnames(bL[[j]])[col] = col_name
		}
	}
	return(bL)
}

# Add a bias term for each non-zero wighted predictor	
# bT - bias term
# bL - beta list

add_bias_term <- function(bL,bT,col=7,col_name = "bias"){
  if(length(bL[[1]]) > 0) {
	for (j in 1:length(bL)) {
		bL[[j]] = cbind(bL[[j]],rep(0,nrow(bL[[j]])))
		colnames(bL[[j]])[length(colnames(bL[[j]]))] = col_name
		for(i in 1:nrow(bL[[j]])){
			bL[[j]][i,col_name] = bT[bL[[j]][i,"trgt"],j]
		}		
	}
  }
	return(bL)
}

# Combine the predictions from several L2 params into one list based on error

combine_l2_net_res <- function(bL,mE,col = "prd_xpln_var"){
	lnet.mat = matrix(0,ncol=3,nrow=0)
	colnames(lnet.mat) = c(colnames(bL[[1]])[1:2],col)
	for (i in 1:nrow(mE)) {
		l = which.min(mE[i,])
		rg_intrs = which(bL[[l]][,"trgt"]==i)
		if(length(rg_intrs)>0){
			lnet.mat = rbind(lnet.mat,bL[[l]][rg_intrs,colnames(lnet.mat)])
		}
	}
	if( col == "bias" ){
		unqModels <- unique(lnet.mat[,1])
                if(length(unqModels)>1){
                  new.mat = matrix(0,ncol=3,nrow=max(unqModels))
                  new.mat[,1] <- c(1:max(unqModels))
                  for( i in 1:length(unqModels)){
                    new.mat[unqModels[i],1] <- lnet.mat[ which(lnet.mat[,1] == unqModels[i])[1], 1]
                    new.mat[unqModels[i],2] <- lnet.mat[ which(lnet.mat[,1] == unqModels[i])[1], 2]
                    new.mat[unqModels[i],3] <- lnet.mat[ which(lnet.mat[,1] == unqModels[i])[1], 3]
                  }
                  lnet.mat <- new.mat
                }
	}
	return(lnet.mat)
}

# Combine two matrices M1 M2
# pre-process: assign M=M1.
# 1) function takes the highest rank (based on values) pixel in M2.
# 2) function finds what is the index of that pixel, denote it by j.
# 3) function adds the highest ranked value of M1 to M[j].

combine_mtrcs_Stouffer <- function(M1,M2){
	ix_m1 = sort(M1,decreasing = TRUE,index.return = TRUE)$ix
	ix_m2 = sort(M2,decreasing = TRUE,index.return = TRUE)$ix
	M = M1
        rm(M1)
	i=1
	while(M2[ix_m2[i]] > 0){
          if(M[ix_m2[i]]>0){
            M[ix_m2[i]] = sqrt(M[ix_m2[i]]^2+M[ix_m1[i]]^2)
          } else {
            M[ix_m2[i]] = sqrt(0+M[ix_m1[i]]^2)
          }
          	i = i+1
	}
	return(M)
}



# Combine two matrices with entries either >0 or 0
# map higest values in tmplt.mat (M1) with sml.mat (M2) that you want to combine to tmplt_mat

combine_mtrcs <- function(M1,M2,qtile=1){
	ix_m1 = sort(M1,decreasing = TRUE,index.return = TRUE)$ix
	ix_m2 = sort(M2,decreasing = TRUE,index.return = TRUE)$ix
        ix_m1 <- ix_m1[round(1+(1-qtile)*(length(ix_m1)-1)):length(ix_m1)]
        if(length(which(M1[ix_m1] > 0)) < length(which(M2 > 0)))
          warning("M2 has more non-zero entries than (possibly reduced) M1")
	M = M1
	i=1
	while(M2[ix_m2[i]] > 0){
		M[ix_m2[i]] = sqrt(M[ix_m2[i]]^2+M1[ix_m1[i]]^2)
		i = i+1
	}
	return(M)
}

combine_mtrcs_new <- function(M1,M2,base.vec=NULL,qtile=1){
  if(is.null(base.vec))
    stop("base.vec required")
	ix_m2 = sort(M2,decreasing = TRUE,index.return = TRUE)$ix
        base.vec <- base.vec[round(1+(1-qtile)*(length(base.vec)-1)):length(base.vec)]
        if(length(base.vec) < length(which(M2 > 0)))
          warning("M2 has more non-zero entries than (possibly reduced) M1")
	M = M1
	i=1
	while(M2[ix_m2[i]] > 0){
		M[ix_m2[i]] = sqrt(M[ix_m2[i]]^2+base.vec[i]^2)
		i = i+1
	}
	return(M)
}

#add_zscore <- function(bL,M1,M2=NULL,col=5,col_name = "clr_zs"){
#	zBL = bL
#	if(length(bL) > 0) {
#		for(i in 1:length(bL)) {
#			zBL[[i]] = cbind(zBL[[i]],0)
#		}
#	for (i in 1:length(bL)) {
#		if(length(bL) > 0) {
#		for(i in 1:length(bL)) {
#			r = bL[[i]][,1]
#			c = bL[[i]][,2]
#			inters = which(c > nrow(M1))
#			singles= which(c <= nrow(M1))
#			if (length(singles) > 0)
#			zBL[[i]][singles,col] = M1[cbind(r[singles],c[singles])]
#			if (length(inters) > 0)
#			zBL[[i]][inters,col] = M2[cbind(r[inters],c[inters]-nrow(M1))]
#			colnames(zBL[[i]])[col] = col_name
#		}
#	}
#	}
#	return(zBL)
#}

add_zscore <- function(bL,M1,M2=NULL,col=5,col_name = "clr_zs"){
	zBL = bL
	if(length(bL[[1]]) > 0) {
		for(i in 1:length(bL)) {
			zBL[[i]] = cbind(zBL[[i]],0)
		}
#	for (i in 1:length(bL)) {
#		if(length(bL) > 0) {
		for(i in 1:length(bL)) {
			r = bL[[i]][,1]
			c = bL[[i]][,2]
			inters = which(c > nrow(M1))
			singles= which(c <= nrow(M1))
			if (length(singles) > 0)
			zBL[[i]][singles,col] = M1[cbind(r[singles],c[singles])]
			if (length(inters) > 0)
			zBL[[i]][inters,col] = M2[cbind(r[inters],c[inters]-nrow(M1))]
			colnames(zBL[[i]])[col] = col_name
		}
	}
#	}
	return(zBL)
}

######################################
make_sparse2 <- function(M)
{	
	non_zero_idx = which(M != 0,arr.ind = TRUE)
	w = matrix(0,nrow(non_zero_idx),3)
	w[,1] = non_zero_idx[,1]
	w[,2] = non_zero_idx[,2]
	w[,3] = M[non_zero_idx]
	return(w)
}
make_sparse_pretty <- function(M) {
  w <- make_sparse2(M)
  w <- w[sort(w[,3], index.return=T, decreasing=T)$ix,]
  w <- w[,c(2,1,3)]
  colnames(w) <- c("tf", "gene", "score")
  return(w)
}
save.predictions.dream5 <- function(M, in.file, cut.off) {
	w <- make_sparse_pretty(M)
	w[,3] <- w[,3] / max(w[,3])
	w[,1] <- paste("G", w[,1], sep="")
	w[,2] <- paste("G", w[,2], sep="")
	w <- w[1:min(cut.off, nrow(w)),]
	write(t(w), file = in.file, sep="\t", append=F, ncolumns=3)
	return(w)
}
######################################
#INPUT
#1) sM - a matrix of dim N*3, where N is the number of interactions
#2) M -  the matrix into which the un-sparsed network go
#OUTPUT
#1) M - binary matrix of unsparsed network. Rows are targets and columns are regulators
unsparse <- function(sM,M) {
  if(nrow(sM)>0) {
    for (i in 1:nrow(sM)){
      M[sM[i,1],sM[i,2]] = as.numeric(sM[i,3])
    }
  }
  return(M)
}

formatOutFNameGP <- function(f.name){
	is.backslash = length(grep("/", f.name)) != 0
	if(!is.backslash)
		return(strsplit(f.name, "\\.")[[1]][1])
	
	#if we reach here then we have backslashes 
	#in the filename
	
	#find the element of the split string
	#that contains a period followed by .txt or .tsv (or.csv)
	x <- strsplit(f.name, "/")[[1]]
	# y <- grep("\\.t*", x)
	return(strsplit(x[length(x)], "\\.")[[1]][1])
}