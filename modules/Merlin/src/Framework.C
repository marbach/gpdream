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
 
	while (1)
	{
		static struct option long_options[] =
			{
				{"data", required_argument, 0, 'm'},
                                {"reg",  required_argument, 0, 'l'},
                                {"cut",  required_argument, 0, 'z'},
                                {0, 0, 0, 0}
			};

		/* getopt_long stores the option index here. */
                int option_index = 0;

                c = getopt_long (argc, argv, "m:k:d:u:l:p:q:r:s:t:c:h:f:z:o:w:v:", long_options, &option_index);

                /* Detect the end of the options. */
                if (c == -1)
                        break;

		switch(c)
		{
			case 'z': // --cut
			{
				break;
			}
			case 'm':
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
			case 'd':
			{
				double convThreshold = atof(optarg);
				metaLearner.setConvergenceThreshold(convThreshold);
				break;
			}
			case 'u':
			{
				double lambda = atof(optarg);
				cout << "Set Lambda" << lambda << endl;
				metaLearner.setLambda(lambda);
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
			case 'q':
			{
				metaLearner.setBeta_ChIP(atof(optarg));
				break;
			}
			case 'r':
			{
				rDefault=false;
				metaLearner.setBeta_Motif(atof(optarg));
				break;
			}
			case 's':
			{
				metaLearner.setChIPGraph(optarg);
				break;
			}
			case 't':
			{
				metaLearner.setMotifGraph(optarg);
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
	}

	if(clusterDefault)
	{
		cerr << "Please specify initial clustering assingment. (option -c)\n";
		return Error::UNKNOWN;
	}
	if(regDefault)
	{
		cerr << "Please input a file of regulators. (option --reg)\n";
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
			<<"--data gene_expression_file" << endl
			<< "-k maxfactorsize " << endl
			 << "-u lambda" << endl
			 << "-t convergence_threshold" << endl
			 << "-v cross_validation_cnt" << endl
			 << "--reg restrictedfname" << endl
			 << "-p beta1" << endl
			 << "-q beta_chip" << endl
			 << "-r beta_motif" << endl
			 << "-s chipgraph" << endl
			 << "-t motifgraph" << endl
			<< "-c genelist"<< endl
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

