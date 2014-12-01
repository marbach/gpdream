##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

## January 2009 Inferelator
## Bonneau lab - Aviv Madar
## NYU - Center for Genomics and Systems Biology

##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
let_usr_choose_dataset <- function() {
  data_is = NULL
  while(is.null(data_is)){
    data_is = switch( menu(c("DREAM5 Network 1","DREAM5 Network 2","DREAM5 Network 3","DREAM5 Network 4","Small DREAM5 Network 2"),
            graphics = FALSE, title = "Run for:"), "DREAM5.1", "DREAM5.2", "DREAM5.3", "DREAM5.4", "DREAM5.2s")
    if(!is.null(data_is)) { break }
    cat("You must choose a data-set. Click ctrl c to exit.")
  }
  return (data_is)
}
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

let_usr_choose_response <- function() {
  response_clr = NULL
  while(is.null(response_clr)){
    response_clr = switch( menu(c("y_k [trivial]","y_k - y_{k-1} [time difference]","(y_k - y_{k-1})/(t_k-t_{k-1}) [rate]",
		"(y_k - y_{k-1})/(t_k-t_{k-1}) + 1/tau*y_k [inf 1]", "(y_k - y_{k-k'})/(t_k-t_{k-k'}) + 1/tau*y_k' [inf 1 all ts intervals]"), 
	    graphics = FALSE, title = "Choose response for CLR:"), 
	"trivial","time_difference","rate","inf_1","inf_1_all_intervals")
    if(!is.null(response_clr)) { break }
    cat("You must choose a response. Click ctrl c to exit.")
  }
  response_lars = NULL
  while(is.null(response_lars)){
    response_lars = switch( menu(c("y_k [trivial]","y_k - y_{k-1} [time difference]","(y_k - y_{k-1})/(t_k-t_{k-1}) [rate]",
		"(y_k - y_{k-1})/(t_k-t_{k-1}) + 1/tau*y_k [inf 1]", "(y_k - y_{k-k'})/(t_k-t_{k-k'}) + 1/tau*y_k' [inf 1 all ts intervals]"), 
	    graphics = FALSE, title = "Choose response for LARS:"), 
	"trivial","time_difference","rate","inf_1","inf_1_all_intervals")
    if(!is.null(response_lars)) { break }
    cat("You must choose a response. Click ctrl c to exit.")
  }
  return( list(response_clr,response_lars) )
}

##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '


let_usr_choose_design_matrix <- function() {
  design_clr = NULL
  while(is.null(design_clr)){
    design_clr = switch( menu(c("X_k [trivial]","X_{k-1} [time delayed]"), 
	    graphics = FALSE, title = "Choose design matrix for CLR:"), 
	"trivial","time_delayed")
    if(!is.null(design_clr)) { break }
    cat("You must choose a design matrix. Click ctrl c to exit.")
  }
  design_lars = NULL
  while(is.null(design_lars)){
    design_lars = switch( menu(c("X_k [trivial]","X_{k-1} [time delayed]"),
	    graphics = FALSE, title = "Choose design matrix for LARS:"), 
	"trivial","time_delayed")
    if(!is.null(design_lars)) { break }
    cat("You must choose a design matrix. Click ctrl c to exit.")
  }
  return( list(design_clr,design_lars) )
}


##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input: cM - colMap
#			r - expression matrix (assumed already normalized)
#			param - a vector of 2 params:
#						1-  c1 - cutoff1: represent the maximal time interval allowed for time series (ts) data
#						2- 'trivial' or 'time_delayed' (param to choose the type of design matrix)
#						3- 'consecutive' or 'all_intervals' (determine if consecutive time measurements [no longer than c1], 
#							 or all permutations of time measurements [up to c1] respectively)
#						4- 'TRUE' or 'FALSE' use_t0_as_steady_state
#						5- 'TRUE' or 'FALSE' use_delt_bigger_than_cutoff_as_steady_state
#
# output: 
#			 steadyStateDesignMat, timeSeriesDesignMat

get_usr_chosen_design_matrix <- function(cM, r, params) {
  delT_min = as.numeric(params[1])
  delT_max = as.numeric(params[2])
  use_t0_as_steady_state = as.logical(params[5])
  use_delt_bigger_than_cutoff_as_steady_state = as.logical(params[6])
  
  delT_vec = sapply(cM, function(i) i$del.t)
  isTs_vec = sapply(cM, function(i) i$isTs)
  eq_idx = which(!isTs_vec) #
  ts_idx = which(isTs_vec)
  
  # data for steady state
  rSS = r[,eq_idx]
  # if we have time series experiments
  if(length(ts_idx)>0) { 
	  delT_vec = delT_vec[ts_idx]
	  # set delT_vec: 
	  #		0 - last time measurement in ts
	  #		>0 - first and middle time measurements in ts
	  # following line make 0s indicate first time measurement in ts
	  delT_vec[which(is.na(delT_vec))] = 0
	  delT_vec_trivial = delT_vec
	  # following 2 lines make 0s indicate last time measurement in ts
	  delT_vec[-length(delT_vec)] = delT_vec[-1]
	  delT_vec[length(delT_vec)] = 0
    # data for time series
    rTS = r[,-eq_idx]
    eq_idx_pseudo = numeric()
    
    # get ts starting conditions, we treat these as equilibrium
    if (use_t0_as_steady_state)
      eq_idx_pseudo = which(delT_vec_trivial == 0)
    
    # get ts conditions with larger than c1 delt, we treat these as equilibrium
    if (use_delt_bigger_than_cutoff_as_steady_state)
      eq_idx_pseudo = c(eq_idx_pseudo, which(delT_vec_trivial > delT_max))
    
    # create design matrix for steady state
    DesignMatSS = cbind(rSS, rTS[,eq_idx_pseudo])
    
    if (params[4] == 'all_intervals') { # all permutations time series
      x = get_all_perms(delT_vec, delT_min, delT_max)
      init_ind = x[[1]]
      boundary_ind = x[[2]]
      DesignMatTS = rTS[,init_ind]
    } else # consecutive time series 
    {
      if(params[3] == 'time_delayed') {
	DesignMatTS = rTS[,which(delT_vec != 0 & delT_vec <= delT_max)]
      } else {
	DesignMatTS = rTS[,which(delT_vec_trivial != 0 & delT_vec_trivial <= delT_max)]
      }
    }
  } else {
    DesignMatSS = rSS
    DesignMatTS = matrix()
  }
  return (list(DesignMatSS,DesignMatTS))
}

##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input: cM - colMap
#			r - ratios matrix (assumed already normalized)
#			param - a vector of 5 params:
#						1-  c1 - cutoff1: represent the maximal time interval allowed for time series (ts) data
#						2- 'trivial','time_difference','rate','inf_1', or 'inf_1_all_intervals'
#						3- tau
#						4- TRUE or FALSE use_t0_as_steady_state
#						5- TRUE or FALSE use_delt_bigger_than_cutoff_as_steady_state
#
# output:
#			 steadyStateResponseMat timeSeriesResponseMat

get_usr_chosen_response  <- function(cM, r, params) {
  delT_min = as.numeric(params[1])
  delT_max = as.numeric(params[2])
  tau = as.numeric(params[4])
  use_t0_as_steady_state = as.logical(params[5])
  use_delt_bigger_than_cutoff_as_steady_state = as.logical(params[6])
  
  delT_vec = sapply(cM, function(i) i$del.t)
  isTs_vec = sapply(cM, function(i) i$isTs)
  eq_idx = which(!isTs_vec)
  ts_idx = which(isTs_vec)
  
  rSS = r[,eq_idx]
  # if we have time series experiments
  if(length(ts_idx)>0) { 
		delT_vec = delT_vec[ts_idx]
		# set delT_vec: 
		#		0 - last time measurement in ts
		#		>0 - first and middle time measurements in ts
		# following line make 0's indicate first time measurement in ts
		delT_vec[which(is.na(delT_vec))] = 0
		delT_vec_trivial = delT_vec
		# following 2 lines make 0's indicate last time measurement in ts
		delT_vec[-length(delT_vec)] = delT_vec[-1]
		delT_vec[length(delT_vec)] = 0
    rTS = r[,-eq_idx]
    
    eq_idx_pseudo = numeric()
    # get ts starting conditions, we treat these as equilibrium
    if (use_t0_as_steady_state)
      eq_idx_pseudo = which(delT_vec_trivial == 0)
    
    # get ts conditions with larger than c1 delt, we treat these as equilibrium
    if (use_delt_bigger_than_cutoff_as_steady_state)
      eq_idx_pseudo = c(eq_idx_pseudo, which(delT_vec_trivial > delT_max))
    
    response_matrixSS = cbind(rSS, rTS[,eq_idx_pseudo])
    
    init_ind = which(delT_vec != 0 & delT_vec <= delT_max)
    boundary_ind = init_ind+1
    
    # finished response matrices now go on to response
    if (params[3] == 'trivial') {
      response_matrixTS = rTS[,boundary_ind]
    } 	else if (params[3] == 'time_difference') {
      response_matrixTS = (rTS[,boundary_ind] - rTS[,init_ind])
    }	else if (params[3] == 'rate') {
      response_matrixTS = t(1/delT_vec[init_ind] * t(rTS[,boundary_ind] - rTS[,init_ind]))
    }	else if (params[3] == 'inf_1') {
      response_matrixTS = t(tau/delT_vec[init_ind] * t(rTS[,boundary_ind] - rTS[,init_ind])) + (rTS[,init_ind])
    } else if (params[3] == 'inf_1_all_intervals') {
      x = get_all_perms(delT_vec, delT_min, delT_max)
      init_ind = x[[1]]
      boundary_ind = x[[2]]
      response_matrixTS = t(tau/delT_vec[init_ind] * t(rTS[,boundary_ind] - rTS[,init_ind])) + (rTS[,init_ind])
    } else {
      stop("unknown response read")
    }
  } else {
    response_matrixSS = rSS
    response_matrixTS = matrix()		
  }
  return(list(response_matrixSS ,response_matrixTS))
}


##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# a helper function for get_usr_chosen_design_matrix and get_usr_chosen_response

get_all_perms <- function(vec, delT_min, delT_max) {
  boundary_idx = numeric()
  init_idx = numeric()
  for (i in 1: length(vec)) {
    j=i;
    local.delT_min <- delT_min
    while ((vec[j] != 0) & (sum(vec[i:j]) <= delT_max)) {
      if(sum(vec[i:j]) >= local.delT_min){
        init_idx = c(init_idx,i);	
        boundary_idx = c(boundary_idx,j+1);
        local.delT_min <- local.delT_min + delT_min
      }
      j=j+1
    }
  }
  return (list(init_idx,boundary_idx))
}
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# input:
#	1- dMSS: design matrix steady state
#	2- dMTS: design matrix time state
#	3- rMSS: response matrix steady state
#	4- dMTS: response matrix time state
#	5- param: what final design matrix? choose from all, ts, or ss

# output:
#	resopnse and corresponding design matrices
make_final_design_and_response_matrix <- function(dMSS, dMTS, rMSS, rMTS, param) {
  if (param == 'all') {
    final_response_matrix = cbind(rMSS, rMTS) 
    final_design_matrix = cbind(dMSS, dMTS)
  } else if (param == 'ts') {
    final_response_matrix = rMTS 
    final_design_matrix = dMTS
  } else if (param == 'ss') {
    final_response_matrix = rMSS		
    final_design_matrix = dMSS
  } else {
    stop("unknown final design or response matrices read")
  }
  return (list(final_response_matrix, final_design_matrix))
}


##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# TODO remove parameter setting from this function (hack)
get_usr_chosen_dataset <- function(dset) {
  obj = list()
  all_data_files <- c("ratios.RData",
      "clusterStack.RData",   
      "colMap.RData",                 
      "tfNames.RData",
      "knockOutFilterList.RData",
      "pkoFilterList.RData",
      "knockOutCombine.RData",
      "pkoCombine.RData",
      "known_regulatory_edges.txt")
  
  if (dset == "DREAM5.1") {
    path = "input/DREAM5/net1/"
    PARAMS[["general"]][["numBoots"]] <<- 60
    PARAMS[["general"]][["tau"]] <<- 60
    PARAMS[["general"]][["pko_rerank_qtile"]] <<- .95
    PARAMS[["general"]][["processorsNumber"]] <<- 12
    PARAMS[["clr"]][["speedUp"]] <<- F
    PARAMS[["clr"]][["numGenes"]] <<- 0
  } else if (dset == "DREAM5.2") {
    path = "input/DREAM5/net2/"
    PARAMS[["general"]][["numBoots"]] <<- 80
    PARAMS[["general"]][["tau"]] <<- 40
    PARAMS[["general"]][["pko_rerank_qtile"]] <<- .95
    PARAMS[["general"]][["processorsNumber"]] <<- 12
    PARAMS[["clr"]][["speedUp"]] <<- F
    PARAMS[["clr"]][["numGenes"]] <<- 0
  } else if (dset == "DREAM5.3") {
    path = "input/DREAM5/net3/"
    PARAMS[["general"]][["numBoots"]] <<- 40
    PARAMS[["general"]][["tau"]] <<- 90
    PARAMS[["general"]][["pko_rerank_qtile"]] <<- .95
    PARAMS[["general"]][["processorsNumber"]] <<- 12
    PARAMS[["general"]][["saveInterval"]] <<- 2
    PARAMS[["clr"]][["speedUp"]] <<- T
    PARAMS[["clr"]][["numGenes"]] <<- 1000
  } else if (dset == "DREAM5.4") {
    path = "input/DREAM5/net4/"
    PARAMS[["general"]][["numBoots"]] <<- 20
    PARAMS[["general"]][["tau"]] <<- 90
    PARAMS[["general"]][["pko_rerank_qtile"]] <<- .95
    PARAMS[["general"]][["processorsNumber"]] <<- 12
    PARAMS[["general"]][["saveInterval"]] <<- 2
    PARAMS[["clr"]][["speedUp"]] <<- T
    PARAMS[["clr"]][["numGenes"]] <<- 1000
  } else if (dset == "DREAM5.2s") {
    path = "input/DREAM5/net2small/"
    PARAMS[["general"]][["numBoots"]] <<- 10
    PARAMS[["general"]][["tau"]] <<- 90
    PARAMS[["general"]][["pko_rerank_qtile"]] <<- .95
    PARAMS[["general"]][["processorsNumber"]] <<- 24
    PARAMS[["clr"]][["speedUp"]] <<- F
    PARAMS[["clr"]][["numGenes"]] <<- 0
  } else {
    return (-1) # failed to read input
  }
  PARAMS[["general"]][["delT_max"]] = PARAMS[["general"]][["tau"]] * 3
  PARAMS[["general"]][["delT_min"]] = PARAMS[["general"]][["tau"]] / 3
  
  n_loaded = 0
  # load data structures except gold standard (which may not exist)
  for (i in 1:length(all_data_files)) {
    file = all_data_files[i]
    cat("loading ", file, "\n")
    if(i == 9) {
      if(file.exists( paste(path,file, sep = "")) ){
        obj[[i]] = load_gold_standard( file=paste(path,file, sep = ""), r_names=rownames(ratios), c_names = tfNames )
      } else {
        # missing file OK for DREAM5
        # stop(file, " does not exist in dir: ", path, " bailing out...",sep='')
      }
    } else {
      if(file.exists( paste(path,file, sep = "")) ){
        load( paste(path,file, sep = "") )
      } else {
        stop(file, " does not exist in dir: ", path, " bailing out...",sep='')
      }
    }
    n_loaded = n_loaded+1
  }

  obj[[1]] = ratios
  obj[[2]] = clusterStack
  obj[[3]] = colMap
  obj[[4]] = tfNames
  obj[[5]] = knockOutFilterList
  obj[[6]] = pkoFilterList
  obj[[7]] = knockOutCombine
  obj[[8]] = pkoCombine
  return( obj )
}

##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

load_gold_standard <- function( file, r_names=NULL, c_names=NULL ) {
  x = as.matrix(read.table(file))
  if(ncol(x) == 3) {
    if(is.null(r_names) | is.null(c_names))
      stop("can't read file: ", file,". missing row names or column names.")
    y = matrix(0,length(r_names),length(c_names))
    rownames(y) = r_names
    colnames(y) = c_names
    for (i in 1:dim(x)[1]) {
      if(all(c((x[i,1] %in% rownames(y)),(x[i,2] %in% rownames(y))))){
	y[x[i,2],x[i,1]] = as.numeric(x[i,3])
      }
    }
    return(y)	
  } 
  return(x)
}

load_predictions <- function( filePath ){
  x <- as.matrix( read.table(filePath,sep="\t") )
  #switch col1 and col2, so that when usparsing rows are targers, cols are predictors...this makes
  #it consistent with all of the other matices of this type that we sue
  temp <- x[,1]
  x[,1] <- x[,2]
  x[,2] <- temp
  x[,3] <- as.double(x[,3]) #to convert from exponential notation
  return(x)
}
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
# read parameters from cmd line args (only read params specified in args.nms). 
# Bails out if arg in args.nms appear more than once or if not appear.
read.cmd.line.params.must <- function(args.nms, cmd.line.args){
	if(length(grep("--version",cmd.line.args))){
		cat("version",script.version,"\n")
		q()
	}
	args <- sapply(strsplit(cmd.line.args," "),function(i) i)
	vals <- character(length(args.nms))
	# split cmd.line to key and value pairs
	for(i in 1:length(args.nms)){
		ix <- grep(args.nms[i],args)
		if(length(ix)>1){
			stop("arg ",args.nms[i]," used more than once.  Bailing out...\n",print.error())
		} else if (length(ix)==0){
			stop("could not find --",args.nms[i],". Bailing out...\n",print.error())
		} else {
			vals[i] <- args[ix+1]
		}
	}
	return(vals)
}
read.cmd.line.params.optional <- function(args.nms, cmd.line.args){
	args <- sapply(strsplit(cmd.line.args," "),function(i) i)
	vals <- character(length(args.nms))
	# split cmd.line to key and value pairs
	for(i in 1:length(args.nms)){
		ix <- grep(args.nms[i],args)
		if(length(ix)>1){
			stop("arg ",args.nms[i]," used more than once.  Bailing out...\n",print.error())
		} else if (length(ix)==0){
			vals[i] <- NA
		} else { 
			vals[i] <- args[ix+1]
	  }
	}
	return(vals)
}


print.error <- function(){
	cat("
					DESCRIPTIION:	
					inf.R takes a ...
					
					INPUT:
					1.--data_file path to expression table\n
					2.--reg_file path to a file that specify which rows are regulators\n
					3.--meta_file path to a file that specify meta data (e.g. which experiments belong to a time series)\n
					4.--inf_max_reg maximum number of regulators considered by elastic net step\n
					5.--n_boots number of bootstraps used to estimate regulatory interactions (linear increase in run time with n_boots)\n
					6.--tau the time scale expected for transcriptional regulatory events to mannifest (only relevant if time series is used; 
					leave as default to have it estimated form data)\n
					7.--num_pred the number of top predictions to be outputed form Inferelator
					
					EXAMPLE RUN: 
					Rscript inf.R
					--data_file dataset.txt
					--reg_file regs.txt
					--meta_file meta.txt
					--inf_max_reg 30
					--n_boots 100
					--tau 15
					--num_pred 10000
					
					Please cite us if you used this script: 
					Greenfield A*, Madar A*, Ostrer H, Bonneau R, 2010 DREAM4: Combining Genetic and Dynamic Information to Identify Biological Networks and Dynamical Models. PLoS ONE 5(10): e13397. doi:10.1371/journal.pone.0013397\n
					")
}


