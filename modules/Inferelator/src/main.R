##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## May th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 1- reads params, design and response matrices, found in PARAMS and INPUT list respectively
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# init PARAMS and INPUT

#params for main (object PARAMS defined in init.R)
b = 1 # this follow current iteration/bootstrap number
N_b = PARAMS$"general"$"numBoots" # number of bootstraps
btch_size = 10 # calculate this mumber of genes to all predictors MI scores in batches (to avoid running out of memory when calculating MI)
percentCoverage <- PARAMS[["general"]][["percentCoverage"]] # (usually 100) percent of matrix that we want to resample
lambda = PARAMS[["lars"]][["lambda"]] # set of l2 norm regularization weights to try in elastic net
cleanUp <- FALSE # clear functions and other intermediate variables at end of run (leaves important variables more visible for end users)
# response and design matrices for clr
Y_clr = INPUT[["clr"]][["response_matrix"]]

###########
X_clr = INPUT[["clr"]][["design_matrix"]] # single predictors
# response and design matrices for lars
Y_lars = INPUT[["lars"]][["response_matrix"]]
###########
X_lars = INPUT[["lars"]][["design_matrix"]] # single predictors
# store results (ODEs,Z scores, and error for each model for each bootstrap run, respectively)
betaList = vector("list", N_b)
modelErrorList = vector("list", N_b)
#startTime <- date() #times how long a run takes
allResults <- list() #list for storing all models

# check that size of dataset is OK
if(ncol(X_lars)<10){
	stop("Too few conditions.  Min number of conditions (experiments) required is 10! (30 or more for robust results). Bailing out...")
} 
if(nrow(X_lars)<5){
	stop("Too few regulators.  Min number of regulators (TFs) required is 5! Bailing out...")
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 2- calculate Median corrected Zscores based on KO data
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# skip for gene pattern version

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 3- setup for bootstrap: create Pi-perm_vector/matrix, Y^pi,X^pi
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

while (b <= N_b) {
  #create permutation matrix
  cat("bootstrap #: ",b,"\n")
  if(b == 1){
    #here we want the original permutation, ie. getOrigPerm = TRUE (i.e. first bootstrap is exact dataset, no resampling)
    Pi_s_clr=createPermMatrix(cS=INPUT[["general"]][["clusterStack"]], allConds = colnames(Y_clr), getOrigPerm = TRUE, percentCoverage = percentCoverage)
    Pi_s_lars=Pi_s_clr
  } else {
    Pi_s_clr=createPermMatrix(cS=INPUT[["general"]][["clusterStack"]], allConds = colnames(Y_clr), getOrigPerm = FALSE, percentCoverage = percentCoverage)
    Pi_s_lars=Pi_s_clr
  }
  #create bicluster specific permutation matrix (ie. read from Pi_g, algorithm described in method comments)
  #this should be changed to be general for both cases where we have only single genes and cases where we havee biclusters
  Y_clr_p = permuteCols(Y_clr,Pi_s_clr)
  X_clr_p = permuteCols(X_clr,Pi_s_clr)
  
  #----added by Alex----9/19-----------#
  #the code below speeds up the calculation of Ms and Ms_bg by
  #making a smaller design matrix
  if( PARAMS$clr$speedUp){
	if(PARAMS$clr$numGenes < nrow(Y_clr)){
	    tfIx <- which(rownames(X_clr_p) %in% INPUT$general$tf_names)
	    otherIx <- c(1:nrow(X_clr_p))[-tfIx]
	    rIx1 <- c(tfIx,sample( otherIx, PARAMS$clr$numGenes -length(tfIx),replace=F))
	} else {
		rIx1 <- 1:nrow(Y_clr)
	}
  }
  
  ##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
  # 4- pass one: fill M - mutual information matrix or correlation matrix
  ##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
  
  # dynamic MI scores stored here
  cat("calculating dynamic MI ")
  #if(PARAMS[["general"]][["processorsNumber"]] > 1){
  if( PARAMS$clr$speedUp){
     Ms <- calc_MI_one_by_one_parallel( Y_clr, X_clr[rIx1,], Pi_s_clr, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
   }
   else{
     Ms <- calc_MI_one_by_one_parallel( Y_clr, X_clr, Pi_s_clr, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
   }

  if(PARAMS$clr$speedUp & ( PARAMS$clr$numGenes < nrow(Y_clr) )){
		x <- cbind(rIx1,1:PARAMS$clr$numGenes)
    for(i in 1:nrow(x)){
    	Ms[x[i,1],x[i,2]] <- 0
    }
  } else {
    diag(Ms) = 0
  }
  
  cat("\n")
  # static MI scores stored here
  cat("calculating background MI ")
  #if(PARAMS[["general"]][["processorsNumber"]] > 1){
  if(PARAMS$clr$speedUp & ( PARAMS$clr$numGenes < nrow(Y_clr) ) ){
     Ms_bg <- calc_MI_one_by_one_parallel( X_clr, X_clr[rIx1,], Pi_s_clr, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
   }else{
     Ms_bg <- calc_MI_one_by_one_parallel( X_clr, X_clr, Pi_s_clr, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
   }

  if(PARAMS$clr$speedUp & ( PARAMS$clr$numGenes < nrow(Y_clr) )){
    x <- cbind(rIx1,1:PARAMS$clr$numGenes)
    for(i in 1:nrow(x)){
      Ms_bg[x[i,1],x[i,2]] <- 0
    }
  }else{
    diag(Ms_bg) = 0
  }
  
  cat("\n")
  
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 5- calculate mixed-CLR (or clr) matrix
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

  if(PARAMS[["general"]][["use_mixCLR"]]){
    cat("running mix-CLR ")
#		Z_nt_fltrd = mixed_clr(Ms_bg,Ms)
    Z_nt_fltrd = mixed_clr_parallel(Ms_bg,Ms,processorsNumber=PARAMS[["general"]][["processorsNumber"]])
  } else {
    cat("running CLR ")
    Z_nt_fltrd = clr(Ms)
  }
  
  cat("\n")
  if(PARAMS$clr$speedUp){
    colnames(Z_nt_fltrd) <- rownames(X_clr)[rIx1]
  }else{
    colnames(Z_nt_fltrd) <- rownames(X_clr)
  }
  rownames(Z_nt_fltrd) <- rownames(X_clr)
  Z <- Z_nt_fltrd[,INPUT[["general"]][["tf_names"]]]
  
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 6- apply MCZ filter -- i.e. remove unlikely reg inters from further consideration by mixedCLR (and thus from Inf)
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
  
  # filter cutoff
  # KOs first
  ct = PARAMS[["general"]][["MCZ_fltr_prcntile"]]
	
  if(!is.null(INPUT$general$knockOutFilterList)){
    x <- INPUT$general$knockOutFilterList
    tfs <- names(x)
    if(!is.null(tfs))
      for(i in 1:length(x)) {
        bad.trgts <- names(x[[ tfs[i] ]])[which(x[[ tfs[i] ]] < quantile(x[[ tfs[i] ]],ct))]
        Z[bad.trgts,tfs[i]] <- 0
      }
  }
  # make sure Z has at least two non-zero predictors for each target (without allowing self regulation)
	n.pred.per.trgt <- apply(Z,1,function(i) length(which(i!=0)))
	ix <- which(n.pred.per.trgt<2)
	if(length(ix)>0){
		min.z <- min(Z[which(Z!=0)])
		ix <- which(n.pred.per.trgt==1)
		if(length(ix)>0){
			for(k in 1:length(ix)){
				ix.replacable <- which(Z[ix[k],]==0)
				ix.bad <- which(names(ix.replacable) %in% names(ix[k]))
				if(length(ix.bad)>0){
					ix.replacable <- ix.replacable[-ix.bad]					
				}
				Z[ix[k],which(Z[ix[k],]==0)[ix.replacable[1]]] <- min.z
			}
		}
		ix <- which(n.pred.per.trgt==0)		
		if(length(ix)>0){
			for(k in 1:length(ix)){
				ix.replacable <- which(Z[ix[k],]==0)
				ix.bad <- which(names(ix.replacable) %in% names(ix[k]))
				if(length(ix.bad)>0){
					ix.replacable <- ix.replacable[-ix.bad]					
				}
				Z[ix[k],which(Z[ix[k],]==0)[ix.replacable[1:2]]] <- min.z
			}
		}
	}

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 7- run Inferelator (elastic net with ODE based modifications to response and design matrices) 
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
  cat("running elastic net ")
  nCv <- min(10,floor(ncol(Y_lars)/2))
  x = calc_ode_model_weights_parallel(Xs = X_lars,Y = Y_lars, Pi = Pi_s_lars, M1 = Z, nS = PARAMS[["lars"]][["max_single_preds"]], nCv = nCv,
	        lambda=lambda, processorsNumber = PARAMS[["general"]][["processorsNumber"]], plot.it = FALSE,
	         plot.file.name = "",verbose = FALSE)


  cat("\n")
  betaList[[b]] = x[[1]]
  modelErrorList[[b]] = t(x[[2]])
  betaList[[b]]=add_weight_beta(bL=betaList[[b]],model_errors=modelErrorList[[1]],n=nrow(Y_lars),pS=nrow(X_lars),pD=0,col=4,col_name = "prd_xpln_var" )
  betaList[[b]]=add_zscore(bL=betaList[[b]],M1=Z,M2=NULL,col=5,col_name = "clr_zs")
  betaList[[b]]=add_bias_term(bL=betaList[[b]],bT=t(x[[3]]),col=6,col_name = "bias")
  rm(x)
  
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 8- run heuristic to combine results from different methods (MCZ, mixCLR, and Inf)--- i.e. results from different pipelines
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
  
  # elastic net produces a model for each l2 norm weight (choose a model for each target from the l2 norm weight with minimum CV error)
  beta.mat = combine_l2_net_res(betaList[[b]],modelErrorList[[b]],col="beta")
  # beta list is a sparse matrix representation.  Turn it into a matrix
  beta.mat = unsparse(beta.mat ,matrix(0,dim(Z)[1],dim(Z)[2]) )
  # same as beta.mat only instead of having beta weight as values it has predictive value for each reg inter
  pred.mat.lnet = combine_l2_net_res(betaList[[b]],modelErrorList[[b]],col="prd_xpln_var")
  pred.mat.lnet = unsparse(pred.mat.lnet,matrix(0,dim(Z)[1],dim(Z)[2]) )
  # for each trgt get the bias term (needed to predict system's response to new perturbations)
  pred.mat.bias = combine_l2_net_res(betaList[[b]],modelErrorList[[b]],col="bias")
  # this is the heuristic described in DREAM3 and DREAM4 papers z = sqrt(z1^2+z2^2)^2
  #  first for DREAM3 pipeline (not additive with MCZ)
  base.vec <- sort(Z,decreasing=T)
  base.vec <- base.vec[which(base.vec>0)]
  pred.mat.lnet.mixCLR = combine_mtrcs_new(Z,pred.mat.lnet,base.vec=base.vec)
  #  second for DREAM5 pipeline
  #   apply threshold to combine z-scores from annotated KOs, then use to push up scores
  if(!is.null(INPUT$general$knockOutCombine)){
 		z.annot.ko <- INPUT$general$knockOutCombine
 	 	z.annot.ko[which(z.annot.ko < PARAMS[["general"]][["z_score_co"]])] <- 0
  	pred.mat.mixCLR.ko <- combine_mtrcs_new(pred.mat.lnet.mixCLR, z.annot.ko, base.vec=base.vec)
	}else {
		pred.mat.mixCLR.ko = NULL
	}

  
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 10- store current re-sampling results
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
  
  allResults[[b]] <- list()
#  allResults[[b]][["betaMat"]] <- beta.mat
#  allResults[[b]][["bias"]] <- pred.mat.bias
  allResults[[b]][["MixCLR.Inf"]] <- pred.mat.lnet.mixCLR
  allResults[[b]][["MixCLR.Inf.ko"]] <- pred.mat.mixCLR.ko
  
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 11- while b<N_b increament b by 1 adn repeat steps 2-6
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
  
  b = b + 1
}

#save(allResults,file=paste(PARAMS$general$saveToDir,"/all_results_",b-1,".RData",sep=""))

if(is.null(allResults[[1]][["MixCLR.Inf.ko"]])){
	median.conf.scores <- getMedianNetworkFromBootstraps(allResults, "MixCLR.Inf")
	dimnames(median.conf.scores) <- dimnames(allResults[[1]][["MixCLR.Inf"]])
}else{
	median.conf.scores <- getMedianNetworkFromBootstraps(allResults, "MixCLR.Inf.ko")
	dimnames(median.conf.scores) <- dimnames(allResults[[1]][["MixCLR.Inf.ko"]])
}

#write.table(median.conf.scores, file=paste("mix_clr_inf_median_scores_",b-1,".xls",sep=""),sep="\t")

#formatting the filename

file.name = paste(formatOutFNameGP(PARAMS$general$d.path), "_InferelatorPipeline_predictions.txt", sep = "")
x <- save.predictions.dream5(median.conf.scores, file.name, PARAMS[["general"]][["num.inters.out"]])

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 13- cleanup tmp variables functions
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

if (cleanUp) {
  rm(numGenesInNet,make_final_design_and_response_matrix,add_bias_term,add_weight_beta,add_zscore,calc_MI_inBatces,calc_ode_model_weights,
      calcDblKoRmsd,calcFoldChange,calcZscores,create_Pi_g,create_Pi_s,create_Xpi,createPermMatrix,fName,get_all_perms,get_best_preds_idx,
      get_usr_chosen_dataset,get_usr_chosen_design_matrix,get_usr_chosen_response,let_usr_choose_dataset,let_usr_choose_design_matrix,
      let_usr_choose_response,load_gold_standard,load_predictions,make_sparse2,makePredictions,modelErrorList,percentCoverage,permuteCols,
      Pi_s_clr,Pi_s_lars,saveInt,splitDreamDataByType,btch_size)
}
