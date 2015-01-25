#include <fstream>
#include <iostream>
#include <cstring>
#include <getopt.h>

using namespace std;
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "LatticeStructure.H"

#include "Vertex.H"
#include "Graph.H"

#include "FactorGraph.H"
#include "FactorManager.H"
#include "PotentialManager.H"
#include "MetaMove.H"
#include "MotifManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"
#include "MetaLearner.H"

#include "Framework.H"

Framework::Framework()
{
	epsThreshold=-1;
}

Framework::~Framework()
{
}

//We will use getopt here
//The options are 
//-m modelname
//-o outputdir
//-e epsilon to control the number of standard deviations above random
//-k maxfactorsize
//-s number of samples for approximate information estimation
//-x k for which we approximate information
//-n cnt of the top candidate MBs to save

Error::ErrorCode
Framework::init(int argc, char** argv)
{
	int c;
	int condCnt = 1;
	//char outFilePrefix[256];
	bool cvDefault = true;
	bool hDefault = true;
	bool pDefault = true;
	bool kDefault = true;
	bool rDefault = true;
	bool oDefault = true;
	bool wDefault = true;
	bool clusterDefault=true;
	bool regDefault=true;
	int optret='-';
	opterr=1;
	int oldoptind=optind;
 
	while(optret=getopt(argc,argv,"o:k:d:v:l:p:r:c:h:f:")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return Error::UNKNOWN;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind ==2)
		{
			my_optarg=argv[oldoptind+1];	
		}
		else
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 'd':
			{
				/*
				// check extension - should be .tsv
				char* ext = strrchr(optarg, '.');
				if(ext != NULL)
				{
					if(strncmp(ext + 1, "tsv", 3) != 0)
					{
						cerr << "expression data should be *.tsv "  << endl;
						return Error::DATAFILE_ERR;
					}
				}
				*/
	
				Error::ErrorCode eCode = varManager.readVariables(optarg);
				if(eCode != Error::SUCCESS)
				{
					cerr << Error::getErrorString(eCode) << endl;
					return eCode;
				}

				evManager.setVariableManager(&varManager);

				eCode = evManager.loadEvidenceFromFile_Continuous(optarg);
				if(eCode != Error::SUCCESS)
				{
					cerr << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				metaLearner.setGlobalEvidenceManager(&evManager);
				metaLearner.setVariableManager(&varManager);
			
				//strncpy(outFilePrefix, optarg, strlen(optarg) - 4);
				//outFilePrefix[strlen(optarg)-4] = '\0';

				break;
			}
			case 'o': // output file - predictions
			{
				oDefault=false;
				//metaLearner.setOutputPredictionName(optarg);
				metaLearner.setOutputDirName(optarg);	
				break;
			}
			/*
			case 'w': // output file - model
			{
				wDefault=false;
				metaLearner.setOutputModuleName(optarg);
				break;
			}
			*/
			case 'k':
			{
				int aSize = atoi(optarg);
				metaLearner.setMaxFactorSize(1);
				metaLearner.setMaxFactorSize_Approx(aSize);
				kDefault=false;
				break;
			}
			case 'v':
			{
				cvDefault=false;
				cvCnt = atoi(optarg);
				if(cvCnt <= 0) 
				{
					cerr << "Cross validation count should be greater than zero.\n";
					return Error::UNKNOWN;
				}
				break;
			}
			case 'l':
			{
				regDefault=false;
				metaLearner.setRestrictedList(optarg);
				break;
			}
			case 'p':
			{
				pDefault=false;
				metaLearner.setBeta1(atof(optarg));
				break;
			}
			case 'r':
			{
				rDefault=false;
				metaLearner.setBeta_Motif(atof(optarg));
				break;
			}
			case 'c':
			{
				clusterDefault=false;
				metaLearner.readModuleMembership(optarg);
				break;
			}
			case 'h':	
			{
				hDefault = false;
				metaLearner.setClusteringThreshold(atof(optarg));
				break;
			}
			case 'f':
			{
				metaLearner.setSpecificFold(atoi(optarg));
				break;
			}
			case '?':
			{
				/* getopt_long already printed an error message. */
				cerr << "Option error " << optopt << endl;
				return Error::UNKNOWN;
			}
			default:
			{
				cerr << "Unhandled option " << c << endl;
				return Error::UNKNOWN;
			}
		}
		oldoptind=optind;
	}

	if(regDefault)
	{
		cerr << "Please input a file of regulators. (option -l)\n";
		return Error::UNKNOWN;
	}
	if(oDefault)
	{
		cerr << "Please specify the name of output directory. (option -o)\n";
		return Error::UNKNOWN;
	}
	/*
	if(wDefault)
	{
		char moduleFName[128];
		sprintf(moduleFName,"%s_MERLIN_modules.txt", outFilePrefix);
		metaLearner.setOutputModuleName(moduleFName);
	}*/
	if(clusterDefault)
	{
		cout << "Setting to default clustering" << endl;
		metaLearner.setDefaultModuleMembership();
	}
	if(cvDefault)
	{
		cvCnt=1;	
	}
	if(hDefault)
	{
		metaLearner.setClusteringThreshold(0.6);
	}
	if(pDefault)
	{
		metaLearner.setBeta1(-5);
	}
	if(kDefault)
	{
		metaLearner.setMaxFactorSize(1);
		metaLearner.setMaxFactorSize_Approx(300);
	}
	if(rDefault)
	{
		metaLearner.setBeta_Motif(4);
	}

	metaLearner.initPartitions(condCnt);
	return Error::SUCCESS;
}

int 
Framework::start()
{
	//metaLearner.start();
	metaLearner.doCrossValidation(cvCnt);
        return 0;
}


int
main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout <<"factorGraphInf " <<  endl
			<<"-d gene_expression_file " << endl
			<< "-k maxfactorsize (default size_of_dataset)" << endl
			 << "-v cross_validation_cnt (default 1)" << endl
			 << "-l restricted_regulator_fname" << endl
			 << "-p sparsity_prior (default -5)" << endl
			 << "-r module_prior (default 4)" << endl
			 << "-o outputdirectory" << endl
			<< "-c clusterassignment (default random_partitioning) "<< endl
			<< "-h hierarchical_clustering_threshold (default 0.6)"<< endl
			<< "-f specificfold_torun (default is -1)" << endl ;
		return 0;
	}
	Framework fw;
 	if(fw.init(argc,argv)!=Error::SUCCESS)
	{
		return -1;
	}
	fw.start();
	return 0;
}

