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

# source scripts
#source("scripts_GenePattern/init_util.R")
#source("scripts_GenePattern/common/utils.R")
#source("scripts_GenePattern/common/larsUtil.R")
#source("scripts_GenePattern/common/clr.R")
#source("scripts_GenePattern/common/bootstrapUtil.R")
#
#source("scripts_GenePattern/d5_util.R")
#source("scripts_GenePattern/prepare_data.R")

# get required packages
suppressMessages(require(multicore,quietly=TRUE))
# library(multicore)
#if require fails
#R CMD INSTALL mi
suppressMessages(require(elasticnet,quietly=TRUE))
suppressMessages(require(inline,quietly=TRUE))
# library(elasticnet)
# library(inline)
#if require fails
#install.packages("elasticnet", repos = "http://lib.stat.cmu.edu/R/CRAN")
suppressMessages(require(samr,quietly=TRUE))


# require(samr)
#code below will install samr
#first install impute, which is part of bioconductor
#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
#the install samr
#install.packages("samr", repos = "http://lib.stat.cmu.edu/R/CRAN")

######################### most useful params ##################################
PARAMS[["lars"]][["lambda"]] = c(0, 1, 100) 
# what is the maximum delta T to consider (in minutes)? (if delta T is bigger than it will be treated as steady state)
PARAMS[["general"]][["delT_max"]] = PARAMS[["general"]][["tau"]] * 3
PARAMS[["general"]][["delT_min"]] = PARAMS[["general"]][["tau"]] / 3
# how many low confidence MCZ interactions (in percentile) do you want to filter out ()i.e. remove from further consideration by tlCLR->inf)?
PARAMS[["general"]][["MCZ_fltr_prcntile"]] = .5
# how many processors to use? if use more than one need to install the multicore package
# PARAMS[["general"]][["processorsNumber"]] = 6
# save after every Nth iteration of bootstrapping
#PARAMS[["general"]][["saveInterval"]] = 10
# What n fold cross validation for inferelator should we use?
PARAMS[["lars"]][["nCv"]] = 10
# How many bins to use to calculate mutual information?
PARAMS[["clr"]][["n.bins"]] = 10
PARAMS[["clr"]][["speedUp"]] = F
PARAMS[["clr"]][["numGenes"]] = 1000
PARAMS[["general"]][["z_score_co"]] = 2
PARAMS[["general"]][["qtile_co"]] = 1
#############################################################################
# if(PARAMS[["general"]][["processorsNumber"]]>1){
	# library(multicore)
# }
PARAMS[["general"]][["inf.version"]] = "nwInf.1.2"
x = unlist(strsplit(date()," "))
PARAMS[["general"]][["date"]] = paste(x[2],x[3],x[5],x[4],sep="_")


PARAMS[["general"]][["use_t0_as_steady_state"]] = FALSE
PARAMS[["general"]][["use_mixCLR"]] = TRUE #for DREAM4
PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]] = TRUE
#making the directory name to save files into
#PARAMS[["general"]][["saveToDir"]] = 
#paste("results/distributions/",PARAMS[["general"]]$data,"_",
#gsub(":","-",PARAMS[["general"]]$date),"_nboots_",PARAMS$"general"$"numBoots",sep="")
PARAMS[["general"]][["percentCoverage"]] = 100

# are we running for 'time-series' only, 'steady state' only, or 'all' for lars and clr respectively
PARAMS[["clr"]][["what_final_design_response_matrix"]] = 'all' # choose here between ts, ss, or all
PARAMS[["lars"]][["what_final_design_response_matrix"]] = 'all' # choose here between ts, ss, or all

#x = let_usr_choose_response()
PARAMS[["clr"]][["response_matrix"]]  = "inf_1_all_intervals"
PARAMS[["lars"]][["response_matrix"]] = "inf_1_all_intervals"

#x = let_usr_choose_design_matrix()
PARAMS[["clr"]][["design_matrix"]]  = "time_delayed"
PARAMS[["lars"]][["design_matrix"]] = "time_delayed"

x = makeDataGP(PARAMS$general$d.path, PARAMS$general$f.path, PARAMS$general$tfs.path)
INPUT[["general"]][["dataset"]]            = x[["ratios"]]
INPUT[["general"]][["clusterStack"]]       = x[["cS"]]
INPUT[["general"]][["colMap"]]             = x[["cM"]]
INPUT[["general"]][["tf_names"]]           = x[["tfs"]]
INPUT[["general"]][["knockOutFilterList"]] = x[["koFilt"]]
INPUT[["general"]][["knockOutCombine"]]    = x[["koComb"]]


#get clr design matrix
#set to all_intervals
#if (PARAMS[["clr"]][["response_matrix"]] == 'inf_1_all_intervals') { 
#	x = 'all_intervals' 
#} else { 
#	x = 'consecutive' 
#}

params = c(PARAMS[["general"]][["delT_min"]], PARAMS[["general"]][["delT_max"]] ,PARAMS[["clr"]][["design_matrix"]], "all_intervals",
			  PARAMS[["general"]][["use_t0_as_steady_state"]], PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
# get clr design matrix: 1- steady_state, 2- time_series
x = get_usr_chosen_design_matrix(INPUT[["general"]][["colMap"]], INPUT[["general"]][["dataset"]],params)

INPUT[["clr"]][["design_matrix_steady_state"]] = x[[1]]
INPUT[["clr"]][["design_matrix_time_series"]] = x[[2]]

#get lars design matrix
#set to all intervals
#if (PARAMS[["lars"]][["response_matrix"]] == 'inf_1_all_intervals') { 
#	x = 'all_intervals' 
#} else { 
#	x = 'consecutive' 
#}

params = c(PARAMS[["general"]][["delT_min"]],PARAMS[["general"]][["delT_max"]],PARAMS[["lars"]][["design_matrix"]], 'all_intervals',
			  PARAMS[["general"]][["use_t0_as_steady_state"]],PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
# get lars design matrix: 1- steady_state, 2- time_series
x = get_usr_chosen_design_matrix(INPUT[["general"]][["colMap"]], as.matrix(as.data.frame(INPUT[["general"]][["dataset"]])[INPUT[["general"]][["tf_names"]],]),params)

INPUT[["lars"]][["design_matrix_steady_state"]] = x[[1]]
INPUT[["lars"]][["design_matrix_time_series"]] = x[[2]]

#get clr response matrix
params =  c(PARAMS[["general"]][["delT_min"]],PARAMS[["general"]][["delT_max"]],PARAMS[["clr"]][["response_matrix"]], PARAMS[["general"]][["tau"]],
				PARAMS[["general"]][["use_t0_as_steady_state"]],PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
if( is.null(INPUT[["general"]][["redExp"]])){
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]], INPUT[["general"]][["dataset"]], params)
}else{
	cat("getting clr response for biclusters \n")
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]],INPUT[["general"]][["redExp"]], params)
}

INPUT[["clr"]][["response_matrix_steady_state"]] = x[[1]]
INPUT[["clr"]][["response_matrix_time_series"]] = x[[2]]

#get lars response matrix
params =  c(PARAMS[["general"]][["delT_min"]],PARAMS[["general"]][["delT_max"]],PARAMS[["lars"]][["response_matrix"]], PARAMS[["general"]][["tau"]],
				PARAMS[["general"]][["use_t0_as_steady_state"]],PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
if( is.null(INPUT[["general"]][["redExp"]])){
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]], INPUT[["general"]][["dataset"]], params)
}else{
	cat("getting lars response for biclusters \n")
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]], INPUT[["general"]][["redExp"]], params)
}
						 
INPUT[["lars"]][["response_matrix_steady_state"]] = x[[1]]
INPUT[["lars"]][["response_matrix_time_series"]] = x[[2]]

# make final design/response matrices for clr
if(any(is.na(INPUT[["clr"]][["response_matrix_steady_state"]]))){
	PARAMS[["clr"]][["what_final_design_response_matrix"]] = "ts"
}else if(any(is.na(INPUT[["clr"]][["response_matrix_time_series"]]))){
	PARAMS[["clr"]][["what_final_design_response_matrix"]] = "ss"
}

x = make_final_design_and_response_matrix(INPUT[["clr"]][["design_matrix_steady_state"]] , 
												  INPUT[["clr"]][["design_matrix_time_series"]] , 
												  INPUT[["clr"]][["response_matrix_steady_state"]], 
												  INPUT[["clr"]][["response_matrix_time_series"]], 
												  PARAMS[["clr"]][["what_final_design_response_matrix"]]) 
												  
INPUT[["clr"]][["response_matrix"]] = x[[1]]
INPUT[["clr"]][["design_matrix"]] = x[[2]]

# make final design/response matrices for lars
if(any(is.na(INPUT[["lars"]][["response_matrix_steady_state"]]))){
	PARAMS[["lars"]][["what_final_design_response_matrix"]] = "ts"
}else if(any(is.na(INPUT[["lars"]][["response_matrix_time_series"]]))){
	PARAMS[["lars"]][["what_final_design_response_matrix"]] = "ss"
}

x = make_final_design_and_response_matrix(INPUT[["lars"]][["design_matrix_steady_state"]] , 
												  INPUT[["lars"]][["design_matrix_time_series"]] , 
												  INPUT[["lars"]][["response_matrix_steady_state"]], 
												  INPUT[["lars"]][["response_matrix_time_series"]], 
												  PARAMS[["lars"]][["what_final_design_response_matrix"]]) 

INPUT[["lars"]][["response_matrix"]] = x[[1]]
INPUT[["lars"]][["design_matrix"]] = x[[2]]

# remove helper variable
rm(x,params)
