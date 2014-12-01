##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

## Feb 2012 Dream5 pipeline (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
##  		     "Alex Greenfield" <ag1868@nyu.edu> 
## NYU - Center for Genomics and Systems Biology

##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

#Rscript ./scripts_GenePattern/runInf.R --data_file input/DREAM5/net3/net3_expression_data_1250.tsv --reg_file input/DREAM5/net3/net3_transcription_factors.tsv --meta_file input/DREAM5/net3/net3_chip_features.tsv --inf_max_reg 30 --n_boots 3 --tau 15 --num_pred 10000 --path_to_scripts ./scripts_GenePattern/
#Rscript ./scripts_GenePattern/runInf.R --data_file input/DREAM5/net3/net3_expression_data_1250.tsv --reg_file input/DREAM5/net3/net3_transcription_factors.tsv --meta_file input/DREAM5/net3/net3_chip_features.tsv --path_to_scripts ./scripts_GenePattern/
#Rscript ./scripts_GenePattern/runInf.R --data_file input/DREAM5/net3/net3_expression_data_1250.tsv --reg_file input/DREAM5/net3/net3_transcription_factors.tsv --inf_max_reg 30 --n_boots 3 --tau 15 --num_pred 10000 --path_to_scripts ./scripts_GenePattern/

rm(list = ls())
test.cmd = F # global parameter for us to mimic command line input
verbose = T

# run.local = F
script.version = "1.5.0"
############################
## help functions to read params from cmd line
# read parameters from cmd line args (only read params specified in args.nms). 
# Bails out if arg in args.nms appear more than once or if not appear.
print.error <- function(){
	cat("
					DESCRIPTIION:	
					the inferelator takes the following inputs.  
					inputs 1 and 2 are required, the rest are optional. 
					the output is a ranked list of predicted regulatory interactions.
					
					INPUT:
					1.--data_file path to expression table\n
					2.--reg_file path to a file that specify which rows are regulators\n
					3.--path_to_scripts the path to were helper scripts reside
					4.--meta_file path to a file that specify meta data (e.g. which experiments belong to a time series)\n
					5.--inf_max_reg maximum number of regulators considered by elastic net step\n
					6.--n_boots number of bootstraps used to estimate regulatory interactions (linear increase in run time with n_boots)\n
					7.--tau the time scale expected for transcriptional regulatory events to mannifest (only relevant if time series is used; 
					leave as default to have it estimated form data)\n
					8.--num_pred the number of top predictions to be outputed form Inferelator\n
					9.--num_processors the number of processorsNumber you want to use (run time drops linearly with core number)
					
					EXAMPLE RUN: 
					Rscript runInf.R
					--data_file dataset.txt
					--reg_file regs.txt
					--path_to_scripts scripts_GenePattern/
					--meta_file meta.txt
					--inf_max_reg 30
					--n_boots 100
					--tau 15
					--num_pred 10000
					--num_processors 2
					
					Please cite us if you used this script: 
					Greenfield A*, Madar A*, Ostrer H, Bonneau R, 2010 DREAM4: Combining Genetic and Dynamic Information to Identify Biological Networks and Dynamical Models. PLoS ONE 5(10): e13397. doi:10.1371/journal.pone.0013397\n
					")
}
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
			stop("could not find ",args.nms[i],". Bailing out...\n",print.error())
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
		} else if (length(ix)==0){ # if --param was not written in cmd line
			vals[i] <- NA
		} else if(( (ix+1) <= length(args) ) & ( length(grep("--",args[ix+1])) == 0) ){
			# if --param was written in cmd line AND not followed by another --param
		 	vals[i] <- args[ix+1]
		} else { # otherwise
			vals[i] <- NA
		}
	}
	return(vals)
}
############################

# retrieve args
#code for testing!!

if(test.cmd==T){								
	# cmd.args <- c("--path_to_scripts ./scripts_GenePattern --data_file ~/Desktop/DREAM5/net3/net3_expression_data_1250.tsv 
	# 			   --reg_file ~/Desktop/DREAM5/net3/net3_transcription_factors.tsv 
	#  					--n_boots 3 --inf_max_reg 30 --tau 40 --num_pred 100000 --num_processors 2")
				cmd.args <- c("--path_to_scripts ./scripts_GenePattern --data_file ~/Desktop/DREAM5/net1/test_expression_data.tsv --reg_file ~/Desktop/DREAM5/net1/test_transcription_factors.tsv 
							   --n_boots 100 --inf_max_reg 30 --tau 40 --num_pred 100000 --num_processors 6" )
} else {
	cmd.args <- commandArgs(trailingOnly = T);      
}

if(verbose)
	cat("\nfinshed reading args: ",cmd.args,"\n")

# args.nms.must <- c(	
# 		"--data_file",   #1
# 		"--reg_file",  	 #2
# 		"--path_to_scripts" #3
# )
args.nms.must <- c(	
		"--data_file",   #1
		"--path_to_scripts" #2
)

# n.start.numeric <- 2
# args.nms.optional <- c(	
# 		"--meta_file", 	 #1
# 		"--inf_max_reg", #2
# 		"--n_boots", 		 #3
# 		"--tau",			   #4	
# 		"--num_pred",     #5
# 		"--num_processors" #6
# )
# n.start.numeric <- 3
args.nms.optional <- c(
		"--reg_file",  	 #1	
		"--meta_file", 	 #2
		"--inf_max_reg", #3
		"--n_boots", 		 #4
		"--tau",			   #5	
		"--num_pred",     #6
		"--num_processors" #7
)


# stop("AM")
# get parameters
vals.must <- read.cmd.line.params.must(args.nms = args.nms.must, cmd.line.args = cmd.args)
if(verbose){
	cat("\nfinshed reading vals: \n")
	cat("\nvals.must", unlist(vals.must), "\n")
}

vals.optional <- read.cmd.line.params.optional(args.nms = args.nms.optional, cmd.line.args = cmd.args)
if(verbose){
	cat("\nvals.optional:", unlist(vals.optional),"\n")
}

INPUT              = list()
INPUT[["general"]] = list()
INPUT[["clr"]]     = list()
INPUT[["lars"]]    = list()

PARAMS              = list()
PARAMS[["general"]] = list()
PARAMS[["clr"]]     = list()
PARAMS[["lars"]]    = list()
PARAMS[["output"]]  = list()

PARAMS[["general"]][["d.path"]] <- vals.must[1]
# PARAMS[["general"]][["tfs.path"]] <- vals.must[2]
PARAMS[["general"]][["scripts.path"]] <- vals.must[2]

if(is.na(vals.optional[1])){
	PARAMS[["general"]][["tfs.path"]] <- NA
} else {
	PARAMS[["general"]][["tfs.path"]] <- vals.optional[1]
}

if(is.na(vals.optional[2])){
	PARAMS[["general"]][["f.path"]] <- NA
} else {
	PARAMS[["general"]][["f.path"]] <- vals.optional[2]	
}

if(is.na(vals.optional[3])){
	PARAMS[["lars"]][["max_single_preds"]] <- 25
}else{
	PARAMS[["lars"]][["max_single_preds"]] <- as.numeric(vals.optional[3])	
}

if(is.na(vals.optional[4])){
	PARAMS[["general"]][["numBoots"]] <- 10
}else{
	PARAMS[["general"]][["numBoots"]] <- as.numeric(vals.optional[4])	
}


if(is.na(vals.optional[5])){
	PARAMS[["general"]][["tau"]] <- 10
}else{
	PARAMS[["general"]][["tau"]] <- as.numeric(vals.optional[5])	
}

if(is.na(vals.optional[6]) || (vals.optional[6] == -1)){
	PARAMS[["general"]][["num.inters.out"]] <- Inf
}else{
	PARAMS[["general"]][["num.inters.out"]] <- as.numeric(vals.optional[6])
}

if(is.na(vals.optional[7])){
	PARAMS[["general"]][["processorsNumber"]] <- 2
}else{
	PARAMS[["general"]][["processorsNumber"]] <- as.numeric(vals.optional[7])
}

numeric.params <- c("numBoots", "tau", "num.inters.out")
# check if numeric params are indeed numeric
for(i in numeric.params){
	if(is.na(as.numeric(PARAMS$general[[paste(i, sep = "")]]))){
		stop("arg ",args.nms.optional[i]," is not numeric.  Bailing out...\n", print.error())
	}
}

if(verbose){
	cat("the input paramse after reading are: \n")
	str(PARAMS[["general"]])
}

source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/init_util.R"))
source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/utils.R"))
source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/larsUtil.R"))
source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/bootstrapUtil.R"))
source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/d5_util.R"))
source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/prepare_data.R"))

if(verbose)
	cat("\nFormatting Data to do Inference Run")

source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/init.R"))
source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/clr.R")) #we are sourcing CLR here bceause a needed library is loaded in init.R

if(verbose)
	cat("\nDone formatting data, executing Inferelator Pipeline\n")

source(paste(sep="",PARAMS[["general"]][["scripts.path"]],"/main.R"))

q(status=0)