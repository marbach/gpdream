#include <iostream>
#include <cstring>
#include <math.h>

#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager()
{
	ludecomp=NULL;
	perm=NULL;
	randgen=gsl_rng_alloc(gsl_rng_default);
	lambda=0;
	randomData=false;
}

PotentialManager::~PotentialManager()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
	for(map<int,Potential*>::iterator pIter=potBuffer.begin();pIter!=potBuffer.end();pIter++)
	{
		delete pIter->second;
	}
	potBuffer.clear();
	gsl_rng_free(randgen);
}

int 
PotentialManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int 
PotentialManager::setRestrictedNeighborSet(map<int,Variable*>& rSet)
{
	for(map<int,Variable*>::iterator vIter=rSet.begin();vIter!=rSet.end();vIter++)
	{
		restrictedNeighborSet[vIter->first]=vIter->second;
	}
	return 0;
}

int
PotentialManager::setOutputDir(const char* aDirName)
{
	strcpy(outputDir,aDirName);
	return 0;
}

int
PotentialManager::setLambda(double lVal)
{
	lambda=lVal;
	return 0;
}

int
PotentialManager::setRandom(bool flag)
{
	randomData=flag;
	return 0;
}

int
PotentialManager::initPooled()
{
	INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	estimateAllMeanCov(randomData,globalMean,globalCovar,trainEvidSet,NULL,NULL);
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::init()
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean.txt",outputDir);
	sprintf(sdFName,"%s/gauss_std.txt",outputDir);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
		estimateAllMeanCov(randomData,globalMean,globalCovar,trainEvidSet,mFName,sdFName);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}


int
PotentialManager::init(int f)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_%d.txt",outputDir,f);
	sprintf(sdFName,"%s/gauss_std_%d.txt",outputDir,f);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
		//estimateAllMeanCov(false,globalMean,globalCovar,trainEvidSet,mFName,sdFName);
		estimateAllMeanCov(randomData,globalMean,globalCovar,trainEvidSet);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}


int
PotentialManager::initValidationSet(int validSetSize)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_v%d.txt",outputDir,validSetSize);
	sprintf(sdFName,"%s/gauss_std_v%d.txt",outputDir,validSetSize);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		INTINTMAP& validationSet=evMgr->getValidationSet();
		estimateAllMeanCov(randomData,globalMean,globalCovar,validationSet,mFName,sdFName);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}


int
PotentialManager::initValidationSet(int validSetSize,int dId)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_v%d.txt",outputDir,validSetSize);
	sprintf(sdFName,"%s/gauss_std_v%d.txt",outputDir,validSetSize);
	INTINTMAP& validationSet=evMgr->getValidationSet();
	estimateAllMeanCov(randomData,globalMean,globalCovar,validationSet,mFName,sdFName,dId);
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::reset()
{
	globalMean.clear();
	for(map<int,INTDBLMAP*>::iterator cIter=globalCovar.begin();cIter!=globalCovar.end();cIter++)
	{
		cIter->second->clear();
		delete cIter->second;
	}
	globalCovar.clear();
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
		ludecomp=NULL;
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
		perm=NULL;
	}
	for(map<int,Potential*>::iterator pIter=potBuffer.begin();pIter!=potBuffer.end();pIter++)
	{
		delete pIter->second;
	}
	potBuffer.clear();
	return 0;
}

int
PotentialManager::resetCache()
{
	globalMean_Cache.clear();
	for(map<int,INTDBLMAP*>::iterator cIter=globalCovar_Cache.begin();cIter!=globalCovar_Cache.end();cIter++)
	{
		cIter->second->clear();
		delete cIter->second;
	}
	globalCovar_Cache.clear();
	return 0;
}

int
PotentialManager::initRandom()
{
	INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	//INTINTMAP& trainEvidSet=evMgr->getTestSet();
	estimateAllMeanCov(true,globalMean_Rand,globalCovar_Rand,trainEvidSet,NULL,NULL);
	return 0;
}

int
PotentialManager::estimateAllMeanCov(bool random, INTDBLMAP& gMean, map<int,INTDBLMAP*>& gCovar,INTINTMAP& trainEvidSet, const char* mFName, const char* sdFName,int leaveOutData)
{
	ofstream mFile;
	ofstream sdFile;

	if(!random)
	{
		if((mFName!=NULL) && (sdFName!=NULL))
		{
			mFile.open(mFName);
			sdFile.open(sdFName);
		}
	}

	int evidCnt=trainEvidSet.size();
	if(leaveOutData!=-1)
	{
		evidCnt=evidCnt-1;
	}
	//First get the mean and then the variance
	int dId=0;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		if(dId==leaveOutData)
		{	
			dId++;
			continue;
		}
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			if(gMean.find(vId)==gMean.end())
			{
				gMean[vId]=val;
			}
			else
			{
				gMean[vId]=gMean[vId]+val;
			}
		}
		dId++;	
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=gMean.begin();idIter!=gMean.end();idIter++)
	{
		if(idIter->first==176)
		{
			//cout <<"Stop here: Variable " << idIter->first << " mean " << idIter->second << endl;
		}
		idIter->second=idIter->second/(double) evidCnt;
		if(!random)
		{
			if(mFile.good())
			{
				mFile<<idIter->first<<"\t" << idIter->second<< endl;
			}
			globalMean_Cache[idIter->first]=idIter->second;
		}
	}
	int covPair=0;
	//Now the variance
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double vval=evid->getEvidVal();
			double vmean=gMean[vId];
			INTDBLMAP* vcov=NULL;
			if(gCovar.find(vId)==gCovar.end())
			{
				vcov=new INTDBLMAP;
				gCovar[vId]=vcov;
			}
			else
			{
				vcov=gCovar[vId];
			}
			for(EMAP_ITER uIter=vIter;uIter!=evidMap->end();uIter++)
			{
				int uId=uIter->first;
				//Don't compute covariance of vId uId pairs that both are not in the restrictedNeighborSet, when
				//the restrictedNeighborSet is empty
			/*	if((!random) && (vId!=uId) && (restrictedNeighborSet.size()>0))
				{
					if((restrictedNeighborSet.find(vId)==restrictedNeighborSet.end()) && (restrictedNeighborSet.find(uId)==restrictedNeighborSet.end()))
					{
						continue;
					}
				}*/
				Evidence* evid1=uIter->second;
				double uval=evid1->getEvidVal();
				double umean=gMean[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(gCovar.find(uId)==gCovar.end())
				{
					ucov=new INTDBLMAP;
					gCovar[uId]=ucov;
				}
				else
				{
					ucov=gCovar[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					covPair++;
					(*vcov)[uId]=diffprod;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						(*ucov)[vId]=diffprod;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}
	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		INTDBLMAP* var_unnorm=new INTDBLMAP;

		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				//vIter->second=2*vIter->second/((double)(gCovar.size()-1));
				//vIter->second=2*vIter->second/((double)(evidCnt-1));
				(*var_unnorm)[vIter->first]=vIter->second;
				vIter->second=(0.001+vIter->second)/((double)(evidCnt-1));
				double variance=vIter->second;
				if(idIter->first==176)
				{
				//	cout <<"Stop here: Variable " << idIter->first << " variance " << idIter->second << endl;
				}
			}
			else
			{
				(*var_unnorm)[vIter->first]=vIter->second;
				vIter->second=vIter->second/((double)(evidCnt-1));
				//vIter->second=vIter->second/((double)(gCovar.size()-1));
				//vIter->second=0;
			}
			if(!random)
			{
				if(sdFile.good())
				{
					sdFile<<idIter->first<<"\t" << vIter->first <<"\t" << vIter->second << endl;
				}
				INTDBLMAP* var_cache=NULL;
				if(globalCovar_Cache.find(idIter->first)==globalCovar_Cache.end())
				{
					var_cache=new INTDBLMAP;
					globalCovar_Cache[idIter->first]=var_cache;
				}
				else
				{
					var_cache=globalCovar_Cache[idIter->first];
				}
				(*var_cache)[vIter->first]=vIter->second;

			}
		}
	}
	if(!random)
	{
		if(mFile.good())
		{
			mFile.close();
		}
		if(sdFile.good())
		{
			sdFile.close();
		}
	}	
	return 0;
}


int
PotentialManager::estimateAllMeanCov(bool random, INTDBLMAP& gMean, map<int,INTDBLMAP*>& gCovar,INTINTMAP& trainEvidSet)
{
	int evidCnt=trainEvidSet.size();
	//First get the mean and then the variance
	int dId=0;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			if(gMean.find(vId)==gMean.end())
			{
				gMean[vId]=val;
			}
			else
			{
				gMean[vId]=gMean[vId]+val;
			}
		}
		dId++;	
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=gMean.begin();idIter!=gMean.end();idIter++)
	{
		idIter->second=idIter->second/(double) evidCnt;
		if(!random)
		{
			globalMean_Cache[idIter->first]=idIter->second;
		}
		INTDBLMAP* vcov=new INTDBLMAP;
		gCovar[idIter->first]=vcov;
	}
	return 0;
	int covPair=0;
	//Now the variance
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double vval=evid->getEvidVal();
			double vmean=gMean[vId];
			INTDBLMAP* vcov=NULL;
			if(gCovar.find(vId)==gCovar.end())
			{
				vcov=new INTDBLMAP;
				gCovar[vId]=vcov;
			}
			else
			{
				vcov=gCovar[vId];
			}
			for(EMAP_ITER uIter=vIter;uIter!=evidMap->end();uIter++)
			{
				int uId=uIter->first;
				//Don't compute covariance of vId uId pairs that both are not in the restrictedNeighborSet, when
				//the restrictedNeighborSet is empty
			/*	if((!random) && (vId!=uId) && (restrictedNeighborSet.size()>0))
				{
					if((restrictedNeighborSet.find(vId)==restrictedNeighborSet.end()) && (restrictedNeighborSet.find(uId)==restrictedNeighborSet.end()))
					{
						continue;
					}
				}*/
				Evidence* evid1=uIter->second;
				double uval=evid1->getEvidVal();
				double umean=gMean[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(gCovar.find(uId)==gCovar.end())
				{
					ucov=new INTDBLMAP;
					gCovar[uId]=ucov;
				}
				else
				{
					ucov=gCovar[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					covPair++;
					(*vcov)[uId]=diffprod;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						(*ucov)[vId]=diffprod;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}
	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				//vIter->second=2*vIter->second/((double)(gCovar.size()-1));
				//vIter->second=2*vIter->second/((double)(evidCnt-1));
				vIter->second=(0.001+vIter->second)/((double)(evidCnt-1));
				double variance=vIter->second;
				if(idIter->first==176)
				{
				//	cout <<"Stop here: Variable " << idIter->first << " variance " << idIter->second << endl;
				}
			}
			else
			{
				vIter->second=vIter->second/((double)(evidCnt-1));
				//vIter->second=vIter->second/((double)(gCovar.size()-1));
				//vIter->second=0;
			}
			if(!random)
			{
				INTDBLMAP* var_cache=NULL;
				if(globalCovar_Cache.find(idIter->first)==globalCovar_Cache.end())
				{
					var_cache=new INTDBLMAP;
					globalCovar_Cache[idIter->first]=var_cache;
				}
				else
				{
					var_cache=globalCovar_Cache[idIter->first];
				}
				(*var_cache)[vIter->first]=vIter->second;

			}
		}
	}
	return 0;
}


int
PotentialManager::estimateCovariance(bool random,INTDBLMAP* vcov, int uId, int vId)
{
	INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	int evidCnt=trainEvidSet.size();
	INTDBLMAP* ucov=globalCovar[uId];
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		Evidence* evid=(*evidMap)[vId];
		double vval=evid->getEvidVal();
		double vmean=globalMean[vId];
		Evidence* evid1=(*evidMap)[uId];
		double uval=evid1->getEvidVal();
		double umean=globalMean[uId];
		double diffprod=(vval-vmean)*(uval-umean);
		if(vcov->find(uId)==vcov->end())
		{
			(*vcov)[uId]=diffprod;
		}
		else
		{
			(*vcov)[uId]=(*vcov)[uId]+diffprod;
		}
		if(uId!=vId)
		{
			if(ucov->find(vId)==ucov->end())
			{
				(*ucov)[vId]=diffprod;
			}
			else
			{
				(*ucov)[vId]=(*ucov)[vId]+diffprod;
			}
		}
	}
	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	if(uId==vId)
	{
		double ssd=(*ucov)[uId];
		(*ucov)[uId]=(0.001+ssd)/((double)(evidCnt-1));
	}
	else
	{
		double ssduv=(*ucov)[vId];
		(*ucov)[vId]=ssduv/((double)(evidCnt-1));
		(*vcov)[uId]=ssduv/((double)(evidCnt-1));
	}
	return 0;
}

//Need to fill in the covariances which are empty
int
PotentialManager::estimateNewCov(INTINTMAP& trainEvidSet, INTINTMAP* varSet, Potential* aPot)
{
	int evidCnt=trainEvidSet.size();
	//First get the mean and then the variance
	int dId=0;
	int covPair=0;
	/*INTDBLMAP meanVect;
	map<int,INTDBLMAP*> covMat;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(eIter->first);
		for(INTINTMAP_ITER vIter=varSet->begin();vIter!=varSet->end();vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=(*evidMap)[vIter->first];
			double vval=evid->getEvidVal();
			if(meanVect.find(vId)==meanVect.end())
			{
				meanVect[vId]=vval;
			}
			else
			{
				meanVect[vId]=meanVect[vId]+vval;
			}
		}
	}
	for(INTDBLMAP_ITER mIter=meanVect.begin();mIter!=meanVect.end();mIter++)
	{
		mIter->second=mIter->second/evidCnt;
	}*/
	/*for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(eIter->first);
		for(INTINTMAP_ITER vIter=varSet->begin();vIter!=varSet->end();vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=(*evidMap)[vIter->first];
			double vval=evid->getEvidVal();
			double vmean=meanVect[vId];
			INTDBLMAP* vcov=NULL;
			if(covMat.find(vId)==covMat.end())
			{
				vcov=new INTDBLMAP;
				covMat[vId]=vcov;
			}
			else
			{
				vcov=covMat[vId];
			}
			for(INTINTMAP_ITER uIter=vIter;uIter!=varSet->end();uIter++)
			{
				int uId=uIter->first;
				Evidence* evid1=(*evidMap)[uIter->first];
				double uval=evid1->getEvidVal();
				double umean=meanVect[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(covMat.find(uId)==covMat.end())
				{
					ucov=new INTDBLMAP;
					covMat[uId]=ucov;
				}
				else
				{
					ucov=covMat[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					covPair++;
					(*vcov)[uId]=diffprod;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						(*ucov)[vId]=diffprod;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}*/
	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	/*for(INTINTMAP_ITER vIter=varSet->begin();vIter!=varSet->end();vIter++)
	{
		INTDBLMAP* var=covMat[vIter->first];
		INTINTMAP_ITER uIter=vIter;
		for(INTINTMAP_ITER uIter=varSet->begin();uIter!=varSet->end();uIter++)
		{
			double cov=(*var)[uIter->first]/((double)(evidCnt-1));
			(*var)[uIter->first]=cov;
		}
	}*/
	VSET& potVars=aPot->getAssocVariables();
	for(VSET_ITER vIter=potVars.begin();vIter!=potVars.end(); vIter++)
	{
		if(globalMean.find(vIter->first)==globalMean.end())
		{
			cerr <<"No var with id " << vIter->first << endl;
			exit(-1);
		}
		double mean=globalMean_Cache[vIter->first];
		INTDBLMAP* covar=globalCovar_Cache[vIter->first];
		aPot->updateMean(vIter->first,mean);
		for(VSET_ITER uIter=vIter;uIter!=potVars.end();uIter++)
		{
			if(covar->find(uIter->first)==covar->end())
			{
				cerr <<"No var " << uIter->first << " in covariance of " << vIter->first << endl;
				exit(-1);
			}
			double cval=(*covar)[uIter->first];
			aPot->updateCovariance(vIter->first,uIter->first,cval);
			aPot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}
	aPot->makeValidJPD(ludecomp,perm);
	return 0;
}


//Need to fill in the covariances which are empty
int
PotentialManager::estimateNewCov_EM(map<int,EvidenceManager*> &evMgrSet, Potential* apot, int cid, map<int,map<int,INTDBLMAP*>*>& gammasubset)
{
	int dId=0;
	int covPair=0;
	VSET& varSet=apot->getAssocVariables();
	INTDBLMAP evidCnts;
	INTDBLMAP meanVect;
	map<int,INTDBLMAP*> covMat;
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		//First get the mean and then the variance
		map<int,INTDBLMAP*>* gammaCond=gammasubset[evIter->first];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			INTDBLMAP* gMap=(*gammaCond)[eIter->first];
			double gval=(*gMap)[cid];
			for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end(); vIter++)
			{
				int vId=vIter->first;
				if(evidCnts.find(vId)==evidCnts.end())
				{
					evidCnts[vId]=gval;
				}
				else
				{
					evidCnts[vId]=evidCnts[vId]+gval;
				}
				Evidence* evid=(*evidMap)[vIter->first];
				double val=evid->getEvidVal();
				val=val*gval;
				if(meanVect.find(vId)==meanVect.end())
				{
					meanVect[vId]=val;
				}
				else
				{
					meanVect[vId]=meanVect[vId]+val;
				}
			}
		}
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=meanVect.begin();idIter!=meanVect.end();idIter++)
	{
		double evidCnt=evidCnts[idIter->first];
		idIter->second=idIter->second/evidCnt;
	}
	//Now the variance
	//Now estimate the variance
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		map<int,INTDBLMAP*>* gammaCond=gammasubset[evIter->first];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[eIter->first];
			double gval=(*gMap)[cid];
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
			{
				int vId=vIter->first;
				Evidence* evid=(*evidMap)[vIter->first];
				double vval=evid->getEvidVal();
				double vmean=meanVect[vId];
				INTDBLMAP* vcov=NULL;
				if(covMat.find(vId)==covMat.end())
				{
					vcov=new INTDBLMAP;
					covMat[vId]=vcov;
				}
				else
				{
					vcov=covMat[vId];
				}
				for(VSET_ITER uIter=vIter;uIter!=varSet.end();uIter++)
				{
					int uId=uIter->first;
					Evidence* evid1=(*evidMap)[uIter->first];
					double uval=evid1->getEvidVal();
					double umean=meanVect[uId];
					double diffprod=gval*(vval-vmean)*(uval-umean);
					INTDBLMAP* ucov=NULL;
					if(covMat.find(uId)==covMat.end())
					{
						ucov=new INTDBLMAP;
						covMat[uId]=ucov;
					}
					else
					{
						ucov=covMat[uId];
					}
					if(vcov->find(uId)==vcov->end())
					{
						covPair++;
						(*vcov)[uId]=diffprod;
					}
					else
					{
						(*vcov)[uId]=(*vcov)[uId]+diffprod;
					}
					if(uId!=vId)
					{
						if(ucov->find(vId)==ucov->end())
						{
							(*ucov)[vId]=diffprod;
						}
						else
						{
							(*ucov)[vId]=(*ucov)[vId]+diffprod;
						}
					}
				}
			}

		}
	}
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		INTDBLMAP* var=covMat[vIter->first];
		for(VSET_ITER uIter=varSet.begin();uIter!=varSet.end();uIter++)
		{
			if(uIter==vIter)
			{
				double cov=(0.001+(*var)[uIter->first])/evidCnts[vIter->first];
				(*var)[uIter->first]=cov;
			}
			else
			{
				double cov=(*var)[uIter->first]/evidCnts[vIter->first];
				(*var)[uIter->first]=cov;
			}
		}
	}
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end(); vIter++)
	{
		double mean=meanVect[vIter->first];
		INTDBLMAP* covar=covMat[vIter->first];
		apot->updateMean(vIter->first,mean);
		for(VSET_ITER uIter=vIter;uIter!=varSet.end();uIter++)
		{
			if(covar->find(uIter->first)==covar->end())
			{
				cerr <<"No var " << uIter->first << " in covariance of " << vIter->first << endl;
				exit(-1);
			}
			double cval=(*covar)[uIter->first];
			apot->updateCovariance(vIter->first,uIter->first,cval);
			apot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}
	apot->makeValidJPD(ludecomp,perm);
	meanVect.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=covMat.begin();vIter!=covMat.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	covMat.clear();
	return 0;
}


int
PotentialManager::estimateAllMeanCov_EM(map<int,EvidenceManager*>& evMgrSet,int cid,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	INTDBLMAP evidCnts;
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		//First get the mean and then the variance
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[evIter->first];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[eIter->first];
			for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
			{
				int vId=vIter->first;
				INTDBLMAP* gMap=(*gammaSet)[vId];
				double gval=(*gMap)[cid];
				if(evidCnts.find(vId)==evidCnts.end())
				{
					evidCnts[vId]=gval;
				}
				else
				{
					evidCnts[vId]=evidCnts[vId]+gval;
				}
				Evidence* evid=vIter->second;
				double val=evid->getEvidVal();
				val=val*gval;
				if(globalMean.find(vId)==globalMean.end())
				{
					globalMean[vId]=val;
				}
				else
				{
					globalMean[vId]=globalMean[vId]+val;
				}
			}
		}
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=globalMean.begin();idIter!=globalMean.end();idIter++)
	{
		double evidCnt=evidCnts[idIter->first];
		idIter->second=idIter->second/evidCnt;
	}
	int covPair=0;
	//Now the variance
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[evIter->first];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[eIter->first];
			for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
			{
				int vId=vIter->first;
				INTDBLMAP* gMap=(*gammaSet)[vId];
				double gval=(*gMap)[cid];
				Evidence* evid=vIter->second;
				double vval=evid->getEvidVal();
				double vmean=globalMean[vId];
				INTDBLMAP* vcov=NULL;
				if(globalCovar.find(vId)==globalCovar.end())
				{
					vcov=new INTDBLMAP;
					globalCovar[vId]=vcov;
				}
				else
				{
					vcov=globalCovar[vId];
				}
				for(EMAP_ITER uIter=vIter;uIter!=evidMap->end();uIter++)
				{
					int uId=uIter->first;
					Evidence* evid1=uIter->second;
					double uval=evid1->getEvidVal();
					double umean=globalMean[uId];
					double diffprod=(vval-vmean)*(uval-umean)*gval;
					INTDBLMAP* ucov=NULL;
					if(globalCovar.find(uId)==globalCovar.end())
					{
						ucov=new INTDBLMAP;
						globalCovar[uId]=ucov;
					}
					else
					{
						ucov=globalCovar[uId];
					}
					if(vcov->find(uId)==vcov->end())
					{
						covPair++;
						(*vcov)[uId]=diffprod;
					}
					else
					{
						(*vcov)[uId]=(*vcov)[uId]+diffprod;
					}
					if(uId!=vId)
					{
						if(ucov->find(vId)==ucov->end())
						{
							(*ucov)[vId]=diffprod;
						}
						else
						{
							(*ucov)[vId]=(*ucov)[vId]+diffprod;
						}
					}
				}
			}

		}
	}
	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=globalCovar.begin();idIter!=globalCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				//vIter->second=2*vIter->second/evidCnt;
				vIter->second=(0.001+vIter->second)/evidCnts[vIter->first];
			}
			else
			{
				vIter->second=vIter->second/evidCnts[vIter->first];
			}
		}
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}


int
PotentialManager::estimateAllMeanCov_EM(map<int,EvidenceManager*>& evMgrSet,int cid,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas, int componentID, INTINTMAP* varSet)
{
	INTDBLMAP evidCnts;
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		//First get the mean and then the variance
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[evIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			for(INTINTMAP_ITER vIter=varSet->begin();vIter!=varSet->end(); vIter++)
			{
				int vId=vIter->first;
				INTDBLMAP* gMap=(*gammaSet)[eIter->first];
				double gval=(*gMap)[cid];
				if(evidCnts.find(vId)==evidCnts.end())
				{
					evidCnts[vId]=gval;
				}
				else
				{
					evidCnts[vId]=evidCnts[vId]+gval;
				}
				Evidence* evid=(*evidMap)[vId];
				double val=evid->getEvidVal();
				val=val*gval;
				if(globalMean.find(vId)==globalMean.end())
				{
					globalMean[vId]=val;
				}
				else
				{
					globalMean[vId]=globalMean[vId]+val;
				}
			}
		}
	}
	//Now estimate the mean
	for(INTINTMAP_ITER vIter=varSet->begin();vIter!=varSet->end();vIter++)
	{
		double evidCnt=evidCnts[vIter->first];
		double mean=globalMean[vIter->first];
		mean=mean/evidCnt;
		globalMean[vIter->first]=mean;
	}
	int covPair=0;
	//Now the variance
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[evIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			for(INTINTMAP_ITER vIter=varSet->begin();vIter!=varSet->end(); vIter++)
			{
				int vId=vIter->first;
				INTDBLMAP* gMap=(*gammaSet)[eIter->first];
				double gval=(*gMap)[cid];
				Evidence* evid=(*evidMap)[vIter->first];
				double vval=evid->getEvidVal();
				double vmean=globalMean[vId];
				INTDBLMAP* vcov=NULL;
				if(globalCovar.find(vId)==globalCovar.end())
				{
					vcov=new INTDBLMAP;
					globalCovar[vId]=vcov;
				}
				else
				{
					vcov=globalCovar[vId];
				}
				for(INTINTMAP_ITER uIter=vIter;uIter!=varSet->end();uIter++)
				{
					int uId=uIter->first;
					Evidence* evid1=(*evidMap)[uIter->first];
					double uval=evid1->getEvidVal();
					double umean=globalMean[uId];
					double diffprod=(vval-vmean)*(uval-umean)*gval;
					INTDBLMAP* ucov=NULL;
					if(globalCovar.find(uId)==globalCovar.end())
					{
						ucov=new INTDBLMAP;
						globalCovar[uId]=ucov;
					}
					else
					{
						ucov=globalCovar[uId];
					}
					if(vcov->find(uId)==vcov->end())
					{
						covPair++;
						(*vcov)[uId]=diffprod;
					}
					else
					{
						(*vcov)[uId]=(*vcov)[uId]+diffprod;
					}
					if(uId!=vId)
					{
						if(ucov->find(vId)==ucov->end())
						{
							(*ucov)[vId]=diffprod;
						}
						else
						{
							(*ucov)[vId]=(*ucov)[vId]+diffprod;
						}
					}
				}
			}

		}
	}
	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	for(INTINTMAP_ITER vIter=varSet->begin();vIter!=varSet->end(); vIter++)
	{
		INTDBLMAP* var=globalCovar[vIter->first];
		for(INTINTMAP_ITER uIter=varSet->begin();uIter!=varSet->end();uIter++)
		{
			double cv=(*var)[uIter->first];
			if(vIter->first==uIter->first)
			{
				cv=(0.001+cv)/evidCnts[uIter->first];
				(*var)[uIter->first]=cv;
			}
			else
			{
				cv=cv/evidCnts[uIter->first];
				(*var)[uIter->first]=cv;
			}
		}
	}
	if(ludecomp==NULL)
	{
		ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	}
	if(perm==NULL)
	{
		perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	}
	return 0;
}


int
PotentialManager::readAllMeanCov(const char* mFName, const char* sdFName)
{

	ifstream mFile(mFName);
	ifstream sdFile(sdFName);
	char buffer[1024];
	while(mFile.good())
	{
		mFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId;
		double mean=0;
		while(tok!=NULL)
		{	
			if(tokCnt==0)
			{
				vId=atoi(tok);	
			}
			else if(tokCnt==1)
			{
				mean=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		globalMean[vId]=mean;
		globalMean_Cache[vId]=mean;
	}
	mFile.close();
	int lineNo=0;
	while(sdFile.good())
	{
		sdFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId=0;
		int uId=0;
		double covariance=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				uId=atoi(tok);
			}
			else if(tokCnt==1)
			{
				vId=atoi(tok);
			}
			else if(tokCnt==2)
			{
				covariance=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		INTDBLMAP* ucov=NULL;
		INTDBLMAP* ucov_cache=NULL;
		if(globalCovar.find(uId)==globalCovar.end())
		{
			ucov=new INTDBLMAP;
			ucov_cache=new INTDBLMAP;
			globalCovar[uId]=ucov;
			globalCovar_Cache[uId]=ucov_cache;
		}
		else
		{
			ucov=globalCovar[uId];
			ucov_cache=globalCovar_Cache[uId];
		}
		(*ucov)[vId]=covariance;
		(*ucov_cache)[vId]=covariance;
		if(uId==0 && vId==0)
		{
			cout << "Found uId=0 vId=0 covar="<< covariance << " at lineno " << lineNo  << endl;
		}
		lineNo++;
	}
	sdFile.close();
	return 0;
}


Error::ErrorCode
PotentialManager::populatePotentialsSlimFactors(map<int,SlimFactor*>& factorSet,VSET& varSet)
{
	//The set of flags to keep status of the potentials that have been calculated
	map<int,bool> doneFlag;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		doneFlag[fIter->first]=false;
	}
	int popFId=0;
	for(map<int,SlimFactor*>::reverse_iterator rIter=factorSet.rbegin();rIter!=factorSet.rend();rIter++)
	{
		//If we have computed the potential for this flag move one
		if(doneFlag[rIter->first])
		{
			popFId++;
			continue;
		}
		SlimFactor* sFactor=rIter->second;
		if(sFactor->fId==176)
		{
			cout <<"Stop here " << endl;
		}
		//Otherwise create the potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc,false);
		aPotFunc->calculateJointEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		if(sFactor->jointEntropy<0)
		{
		//	sFactor->jointEntropy=0;
		//	cout <<"Negative entropy for " << sFactor->fId << endl;
		}
		doneFlag[rIter->first]=true;
		delete aPotFunc;
		if(popFId%100000==0)
		{
			cout <<"Done with " << factorSet.size()-popFId << " factors " << endl;
		}
		popFId++;
	}
	return Error::SUCCESS;
}


int 
PotentialManager::estimateMarginalEntropies(map<int,SlimFactor*>& slimFactors,VSET& varSet,bool random)
{	
	for(map<int,SlimFactor*>::iterator aIter=slimFactors.begin();aIter!=slimFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		Potential* aPotFunc=new Potential;
		Variable* aVar=varSet[sFactor->vIds[0]];
		aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc,random);
		aPotFunc->calculateJointEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		delete aPotFunc;
	}
	return 0;
}

//Estimate the random information for all factors of a particular size. Assume that
//randInfo is already allocated
Error::ErrorCode
PotentialManager::estimateRandomInfo(map<int,SlimFactor*>& factorSet, 
		VSET& varSet, vector<double>& randInfo, int fSize)
{
	int rInd=0;
	for(map<int,SlimFactor*>::iterator rIter=factorSet.begin();rIter!=factorSet.end();rIter++)
	{
		SlimFactor* sFactor=rIter->second;
		if(sFactor->vCnt<fSize)
		{
			continue;
		}
		else if(sFactor->vCnt>fSize)
		{
			break;
		}
		//Otherwise create a potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc,true);
		aPotFunc->calculateJointEntropy();
		double rInfo=(-1)*aPotFunc->getJointEntropy();
		for(int j=0;j<sFactor->vCnt;j++)
		{
			SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
			rInfo=rInfo+subFactor->jointEntropy;
		}
		randInfo.push_back(rInfo);
		rInd++;
		delete aPotFunc;
	}
	return Error::SUCCESS;

}


int 
PotentialManager::populateFactor(map<int,SlimFactor*>& factorSet,VSET& varSet,SlimFactor* sFactor,bool random)
{
	Potential* aPotFunc=new Potential;
	//string fullConfStr;
	char confStr[CONSTR_LEN];
	for(int j=0;j<sFactor->vCnt;j++)
	{
		Variable* aVar=varSet[sFactor->vIds[j]];
		if(j==sFactor->vCnt-1)
		{
			
			aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		}
		else
		{
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		//sprintf(confStr,"-%d",sFactor->vIds[j]);
		//fullConfStr.append(confStr);
	}
	aPotFunc->potZeroInit();
	populatePotential(aPotFunc,random);
	aPotFunc->calculateJointEntropy();
	//jointEntropies[fullConfStr]=aPotFunc->getJointEntropy();
	double rInfo=(-1)*aPotFunc->getJointEntropy();
	for(int j=0;j<sFactor->vCnt;j++)
	{
		SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
		rInfo=rInfo+subFactor->jointEntropy;
	}
	sFactor->mutualInfo=rInfo;
	
	delete aPotFunc;
	return 0;
}


int 
PotentialManager::populateFactor_Buffer(map<int,SlimFactor*>& factorSet,VSET& varSet,SlimFactor* sFactor,bool random)
{
	Potential* aPotFunc=NULL;
	if(potBuffer.find(sFactor->vCnt)==potBuffer.end())
	{
		aPotFunc=new Potential;
		aPotFunc->initMatrix(sFactor->vCnt);
		potBuffer[sFactor->vCnt]=aPotFunc;
	}
	else
	{
		aPotFunc=potBuffer[sFactor->vCnt];
		aPotFunc->resetVarSet();
	}
	//string fullConfStr;
	char confStr[CONSTR_LEN];
	for(int j=0;j<sFactor->vCnt;j++)
	{
		Variable* aVar=varSet[sFactor->vIds[j]];
		if(j==sFactor->vCnt-1)
		{
			
			aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		}
		else
		{
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		//sprintf(confStr,"-%d",sFactor->vIds[j]);
		//fullConfStr.append(confStr);
	}
	aPotFunc->potZeroInit_MeanOnly();
	populatePotential(aPotFunc,random);
	aPotFunc->calculateJointEntropy();
	//jointEntropies[fullConfStr]=aPotFunc->getJointEntropy();
	double rInfo=(-1)*aPotFunc->getJointEntropy();
	for(int j=0;j<sFactor->vCnt;j++)
	{
		SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
		rInfo=rInfo+subFactor->jointEntropy;
	}
	sFactor->mutualInfo=rInfo;
	
	//delete aPotFunc;
	return 0;
}

int
PotentialManager::populatePotential(Potential* aPot, bool random)
{
	VSET& potVars=aPot->getAssocVariables();
	for(VSET_ITER vIter=potVars.begin();vIter!=potVars.end(); vIter++)
	{
		double mean=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
	/*	if(random)
		{
			if(globalMean_Rand.find(vIter->first)==globalMean_Rand.end())
			{
				cerr <<"No var with id " << vIter->first << endl;
				exit(-1);
			}
			mean=globalMean_Rand[vIter->first];
			covar=globalCovar_Rand[vIter->first];
		}
		else
		{*/
			if(globalMean.find(vIter->first)==globalMean.end())
			{
				cerr <<"No var with id " << vIter->first << endl;
				exit(-1);
			}
			mean=globalMean[vIter->first];
			covar=globalCovar[vIter->first];
		//}
		aPot->updateMean(vIter->first,mean);
		for(VSET_ITER uIter=vIter;uIter!=potVars.end();uIter++)
		{
			if(covar->find(uIter->first)==covar->end())
			{
				if(uIter->first==9 && vIter->first==9)
				{
					cout <<"Estimating the covariance of " << uIter->first << endl;
				}
				estimateCovariance(random,covar,uIter->first,vIter->first);
				//cerr <<"No var " << uIter->first << " in covariance of " << vIter->first << endl;
				//exit(-1);
			}
			double cval=(*covar)[uIter->first];
			aPot->updateCovariance(vIter->first,uIter->first,cval);
			aPot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}

	//aPot->makeValidJPD();
	aPot->makeValidJPD(ludecomp,perm);
	return 0;
}

int
PotentialManager::populatePotential_EM(Potential* apot,int vId, int cid, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas, bool newB)
{
	double vmean=globalMean[vId];
	INTDBLMAP* vcovars=globalCovar[vId];
	double vvar=(*vcovars)[vId];
	double b=apot->getCondBias();
	if(newB)
	{
		//b=gsl_ran_gaussian(randgen,sqrt(vvar));
		//b=sqrt(vvar);
		b=vmean;
	}
	//Now estimate the weight
	estimatePotWeight(apot,vId,cid,b,evSet,gammas);
	estimatePotCondVarBias(apot,vId,cid,b,evSet,gammas,cid);
	estimatePotWeight(apot,vId,cid,apot->getCondBias(),evSet,gammas);

	return 0;
}

int
PotentialManager::estimatePotWeight(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTINTMAP& vIDMatIDMap=apot->getVarMatrixIndexMap();
	map<int,INTDBLMAP*> sumys;
	INTDBLMAP sumxy;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[vId];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double xpart=(*evidSet)[vId]->getEvidVal()-b;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				if(sumxy.find(vIter->first)==sumxy.end())
				{	
					sumxy[vIter->first]=gammaval*xpart*vval;
				}
				else
				{
					sumxy[vIter->first]=sumxy[vIter->first]+(gammaval*xpart*vval);
				}
				VSET_ITER uIter=vIter;
				for(;uIter!=potVarSet.end();uIter++)
				{
					if(uIter->first==vId)
					{
						continue;
					}
					double uval=(*evidSet)[uIter->first]->getEvidVal();
					double uvval=vval*uval*gammaval;
					uvval=uvval+lambda;
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(vIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[vIter->first]=vrow;
					}
					else
					{
						vrow=sumys[vIter->first];
					}	
					if(vrow->find(uIter->first)==vrow->end())
					{
						(*vrow)[uIter->first]=uvval;
					}
					else
					{
						(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
					}
					if(uIter==vIter)
					{
						continue;
					}
					if(sumys.find(uIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[uIter->first]=urow;
					}
					else
					{
						urow=sumys[uIter->first];
					}
					if(urow->find(vIter->first)==urow->end())
					{
						(*urow)[vIter->first]=uvval;
					}
					else
					{
						(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
					}
				}
			}
		}
	}
	Matrix* covMatrix=new Matrix(potVarSet.size()-1,potVarSet.size()-1);
	Matrix* xyMatrix=new Matrix(1,potVarSet.size()-1);
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowvid=vIDMatIDMap[vId];
		int rowpos=vIDMatIDMap[vIter->first];
		if(rowpos>rowvid)
		{
			rowpos--;
		}
		double xyval=sumxy[vIter->first];
		xyMatrix->setValue(xyval,0,rowpos);
		for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
		{
			int colvid=vIDMatIDMap[vId];
			int colpos=vIDMatIDMap[rIter->first];
			if(colpos>colvid)
			{
				colpos--;
			}
			covMatrix->setValue(rIter->second,rowpos,colpos);
			covMatrix->setValue(rIter->second,colpos,rowpos);
		}
	}
	Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
	Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
	apot->setCondWeight(wtMatrix);
	delete covMatrix;
	delete xyMatrix;
	delete inv;
	delete wtMatrix;
	sumxy.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	sumys.clear();
	return 0;

}


int
PotentialManager::estimatePotWeight(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTINTMAP& vIDMatIDMap=apot->getVarMatrixIndexMap();
	map<int,INTDBLMAP*> sumys;
	INTDBLMAP sumxy;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double xpart=(*evidSet)[vId]->getEvidVal()-b;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				if(sumxy.find(vIter->first)==sumxy.end())
				{	
					sumxy[vIter->first]=gammaval*xpart*vval;
				}
				else
				{
					sumxy[vIter->first]=sumxy[vIter->first]+(gammaval*xpart*vval);
				}
				VSET_ITER uIter=vIter;
				for(;uIter!=potVarSet.end();uIter++)
				{
					if(uIter->first==vId)
					{
						continue;
					}
					double uval=(*evidSet)[uIter->first]->getEvidVal();
					double uvval=vval*uval*gammaval;
					uvval=uvval+lambda;
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(vIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[vIter->first]=vrow;
					}
					else
					{
						vrow=sumys[vIter->first];
					}	
					if(vrow->find(uIter->first)==vrow->end())
					{
						(*vrow)[uIter->first]=uvval;
					}
					else
					{
						(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
					}
					if(uIter==vIter)
					{
						continue;
					}
					if(sumys.find(uIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[uIter->first]=urow;
					}
					else
					{
						urow=sumys[uIter->first];
					}
					if(urow->find(vIter->first)==urow->end())
					{
						(*urow)[vIter->first]=uvval;
					}
					else
					{
						(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
					}
				}
			}
		}
	}
	Matrix* covMatrix=new Matrix(potVarSet.size()-1,potVarSet.size()-1);
	Matrix* xyMatrix=new Matrix(1,potVarSet.size()-1);
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowvid=vIDMatIDMap[vId];
		int rowpos=vIDMatIDMap[vIter->first];
		if(rowpos>rowvid)
		{
			rowpos--;
		}
		double xyval=sumxy[vIter->first];
		xyMatrix->setValue(xyval,0,rowpos);
		for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
		{
			int colvid=vIDMatIDMap[vId];
			int colpos=vIDMatIDMap[rIter->first];
			if(colpos>colvid)
			{
				colpos--;
			}
			covMatrix->setValue(rIter->second,rowpos,colpos);
			covMatrix->setValue(rIter->second,colpos,rowpos);
		}
	}
	Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
	Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
	apot->setCondWeight(wtMatrix);
	delete covMatrix;
	delete xyMatrix;
	delete inv;
	delete wtMatrix;
	sumxy.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	sumys.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_TiedParam(map<int,Potential*>& potSet, int uId, int vId, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas,int componentID)
{
	INTDBLMAP numpart;
	INTDBLMAP denpart;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				double b=apot->getCondBias();
				VSET& potVarSet=apot->getAssocVariables();
				double xpart=(*evidSet)[uId]->getEvidVal()-b;
				double gammaval=(*gMap)[pIter->first];
				INTDBLMAP& oldWeight=apot->getCondWeight();
				double ypart=0;
				double vval=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==uId)
					{
						continue;
					}
					if(vIter->first==vId)
					{
						vval=(*evidSet)[vIter->first]->getEvidVal();
						continue;
					}
					double aval=(*evidSet)[vIter->first]->getEvidVal();
					double wt=oldWeight[vIter->first];
					ypart=ypart+(aval*wt);
				}
				xpart=gammaval*(xpart-ypart);
				if(numpart.find(pIter->first)==numpart.end())
				{
					numpart[pIter->first]=xpart*vval;
				}
				else
				{
					numpart[pIter->first]=numpart[pIter->first]+(xpart*vval);
				}
				double denval=vval*vval*gammaval;
				if(denpart.find(pIter->first)==denpart.end())
				{
					denpart[pIter->first]=denval;
				}
				else
				{
					denpart[pIter->first]=denpart[pIter->first]+denval;
				}
			}
		}
	}
	double numerator=0;
	double denominator=0;
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		double condvar=pot->getCondVariance();
		double npart=numpart[pIter->first];
		npart=npart/condvar;
		numerator=numerator+npart;
		double dpart=denpart[pIter->first];
		dpart=dpart/condvar;
		denominator=denominator+dpart;
	}
	double wt=numerator/denominator;

	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		pot->setCondWeightFor(vId,wt);
	}
	numpart.clear();
	denpart.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_TiedParam(map<int,Potential*>& potSet, int uId, int vId, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	INTDBLMAP numpart;
	INTDBLMAP denpart;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				double b=apot->getCondBias();
				VSET& potVarSet=apot->getAssocVariables();
				double xpart=(*evidSet)[uId]->getEvidVal()-b;
				double gammaval=(*gMap)[pIter->first];
				INTDBLMAP& oldWeight=apot->getCondWeight();
				double ypart=0;
				double vval=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==uId)
					{
						continue;
					}
					if(vIter->first==vId)
					{
						vval=(*evidSet)[vIter->first]->getEvidVal();
						continue;
					}
					double aval=(*evidSet)[vIter->first]->getEvidVal();
					double wt=oldWeight[vIter->first];
					ypart=ypart+(aval*wt);
				}
				xpart=gammaval*(xpart-ypart);
				if(numpart.find(pIter->first)==numpart.end())
				{
					numpart[pIter->first]=xpart*vval;
				}
				else
				{
					numpart[pIter->first]=numpart[pIter->first]+(xpart*vval);
				}
				double denval=vval*vval*gammaval;
				if(denpart.find(pIter->first)==denpart.end())
				{
					denpart[pIter->first]=denval;
				}
				else
				{
					denpart[pIter->first]=denpart[pIter->first]+denval;
				}
			}
		}
	}
	double numerator=0;
	double denominator=0;
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		double condvar=pot->getCondVariance();
		double npart=numpart[pIter->first];
		npart=npart/condvar;
		numerator=numerator+npart;
		double dpart=denpart[pIter->first];
		dpart=dpart/condvar;
		denominator=denominator+dpart;
	}
	double wt=numerator/denominator;

	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		pot->setCondWeightFor(vId,wt);
	}
	numpart.clear();
	denpart.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_TiedParamSet(map<int,Potential*>& potSet, int xId, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	INTDBLMAP sumxy;
	map<int,INTDBLMAP*> sumys;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			INTDBLMAP exclParts;
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				double b=apot->getCondBias();
				VSET& potVarSet=apot->getAssocVariables();
				INTINTMAP sharedMBVars=apot->getSharedMBVars();
				INTDBLMAP& oldWeight=apot->getCondWeight();
				double exclpart=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==xId)
					{
						continue;
					}
					if(sharedMBVars.find(vIter->first)!=sharedMBVars.end())
					{
						continue;
					}
					if(oldWeight.find(vIter->first)==oldWeight.end())
					{
						cerr <<"No weight for  " << vIter->first << " in cond of " << xId << endl;
						exit(-1);
					}
					Evidence* evid=(*evidSet)[vIter->first];
					double aval=oldWeight[vIter->first]*evid->getEvidVal();
					exclpart=exclpart+aval;
				}
				double xpart=(*evidSet)[xId]->getEvidVal()-b-exclpart;
				double c=apot->getCondVariance();
				double gval=(*gMap)[pIter->first];
				xpart=(gval*xpart)/c;
				exclParts[pIter->first]=xpart;
			}
			Potential* apot=potSet.begin()->second;
			INTINTMAP& sharedMBVars=apot->getSharedMBVars();
			for(INTINTMAP_ITER sIter=sharedMBVars.begin();sIter!=sharedMBVars.end();sIter++)
			{
				double aval=(*evidSet)[sIter->first]->getEvidVal();
				double asum=0;
				for(INTDBLMAP_ITER pIter=exclParts.begin();pIter!=exclParts.end();pIter++)
				{
					asum=asum+(aval*pIter->second);
				}
				if(sumxy.find(sIter->first)==sumxy.end())
				{
					sumxy[sIter->first]=asum;
				}
				else
				{
					sumxy[sIter->first]=sumxy[sIter->first]+asum;
				}
				INTINTMAP_ITER qIter=sIter;
				for(;qIter!=sharedMBVars.end();qIter++)
				{
					double bval=(*evidSet)[qIter->first]->getEvidVal();
					double uvval=0;
					for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
					{
						Potential* apot=pIter->second;
						double gammaval=(*gMap)[pIter->first];
						double c=apot->getCondVariance();
						double dval=(gammaval*aval*bval)/c;
						uvval=uvval+dval;
					}
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(sIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[sIter->first]=vrow;
					}
					else
					{
						vrow=sumys[sIter->first];
					}	
					if(vrow->find(qIter->first)==vrow->end())
					{
						(*vrow)[qIter->first]=uvval;
					}
					else
					{
						(*vrow)[qIter->first]=(*vrow)[qIter->first]+uvval;
					}
					if(sIter==qIter)
					{
						continue;
					}
					if(sumys.find(qIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[qIter->first]=urow;
					}
					else
					{
						urow=sumys[qIter->first];
					}
					if(urow->find(qIter->first)==urow->end())
					{
						(*urow)[sIter->first]=uvval;
					}
					else
					{
						(*urow)[sIter->first]=(*urow)[sIter->first]+uvval;
					}

				}
			}
			exclParts.clear();
		}
	}
	Potential* apot=potSet.begin()->second;
	int sharedNeighbors=apot->getSharedMBVars().size();
	Matrix* covMatrix=new Matrix(sharedNeighbors,sharedNeighbors);
	Matrix* xyMatrix=new Matrix(1,sharedNeighbors);
	map<int,int> vIDMatIDMap;
	int rowvid=0;
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowpos=-1;
		if(vIDMatIDMap.find(vIter->first)==vIDMatIDMap.end())
		{
			vIDMatIDMap[vIter->first]=rowvid;
			rowvid++;
		}
		rowpos=vIDMatIDMap[vIter->first];
		double xyval=sumxy[vIter->first];
		xyMatrix->setValue(xyval,0,rowpos);
		for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
		{
			int colpos=-1;
			if(vIDMatIDMap.find(rIter->first)==vIDMatIDMap.end())
			{
				vIDMatIDMap[rIter->first]=rowvid;
				rowvid++;
			}
			colpos=vIDMatIDMap[rIter->first];
			covMatrix->setValue(rIter->second,rowpos,colpos);
			covMatrix->setValue(rIter->second,colpos,rowpos);
		}
	}
	Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
	Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* lpot=pIter->second;
		INTDBLMAP& potWt=lpot->getCondWeight();
		for(INTINTMAP_ITER vIter=vIDMatIDMap.begin();vIter!=vIDMatIDMap.end();vIter++)
		{
			int pos=vIter->second;
			double wval=wtMatrix->getValue(0,pos);
			lpot->setCondWeightFor(vIter->first,wval);
		}
	}
	vIDMatIDMap.clear();
	delete covMatrix;
	delete xyMatrix;
	delete inv;
	delete wtMatrix;
		
	sumxy.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	sumys.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_ExclParam(Potential* pot, int uId, int vId, int cid, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas, int componentID)
{
	double numerator=0;
	double denominator=0;
	double b=pot->getCondBias();
	VSET& potVarSet=pot->getAssocVariables();
	INTDBLMAP& oldWeight=pot->getCondWeight();
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xpart=(*evidSet)[uId]->getEvidVal()-b;
			double ypart=0;
			double vval=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==uId)
				{
					continue;
				}
				if(vIter->first==vId)
				{
					vval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double aval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=oldWeight[vIter->first];
				ypart=ypart+(aval*wt);
			}
			xpart=gammaval*(xpart-ypart);
			numerator=numerator+(xpart*vval);
			double denval=vval*vval*gammaval;
			denominator=denominator+denval;
		}
	}
	double wt=numerator/denominator;
	pot->setCondWeightFor(vId,wt);
	return 0;
}

int
PotentialManager::estimatePotWeight_ExclParam(Potential* pot, int uId, int vId, int cid, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	double numerator=0;
	double denominator=0;
	double b=pot->getCondBias();
	VSET& potVarSet=pot->getAssocVariables();
	INTDBLMAP& oldWeight=pot->getCondWeight();
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xpart=(*evidSet)[uId]->getEvidVal()-b;
			double ypart=0;
			double vval=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==uId)
				{
					continue;
				}
				if(vIter->first==vId)
				{
					vval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double aval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=oldWeight[vIter->first];
				ypart=ypart+(aval*wt);
			}
			xpart=gammaval*(xpart-ypart);
			numerator=numerator+(xpart*vval);
			double denval=vval*vval*gammaval;
			denominator=denominator+denval;
		}
	}
	double wt=numerator/denominator;
	pot->setCondWeightFor(vId,wt);
	return 0;
}


int 
PotentialManager::estimatePotWeight_ExclParamSet(map<int,Potential*>& potSet,int xId, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	map<int,map<int,INTDBLMAP*>*> sumys_set;
	map<int,INTDBLMAP*> sumxy_set;
	Potential* apot=potSet.begin()->second;
	INTINTMAP& sharedMBVars=potSet.begin()->second->getSharedMBVars();
	INTDBLMAP& wt=apot->getCondWeight();
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaSet=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			//First get the shared part
			//The parameters of the shared part should be the same
			double sharedpart=0;
			for(INTINTMAP_ITER vIter=sharedMBVars.begin();vIter!=sharedMBVars.end();vIter++)
			{
				if(wt.find(vIter->first)==wt.end())
				{
					cerr <<"Error! No weight value for " << vIter->first << " in cond of " << xId << endl;
					exit(-1);
				}
				Evidence* evid=(*evidSet)[vIter->first];
				double aval=wt[vIter->first]*evid->getEvidVal();
				sharedpart=sharedpart+aval;
			}
			//Now the specific parts
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				VSET& potVarSet=apot->getAssocVariables();
				double gammaval=(*gMap)[pIter->first];
				double b=apot->getCondBias();
				double xpart=(*evidSet)[xId]->getEvidVal()-b-sharedpart;
				INTDBLMAP* sumxy=NULL;
				map<int,INTDBLMAP*>* sumys=NULL;
				if(sumxy_set.find(pIter->first)==sumxy_set.end())
				{
					sumxy=new INTDBLMAP;
					sumxy_set[pIter->first]=sumxy;
					sumys=new map<int,INTDBLMAP*>;
					sumys_set[pIter->first]=sumys;
				}
				else
				{	
					sumxy=sumxy_set[pIter->first];
					sumys=sumys_set[pIter->first];
				}
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==xId)
					{
						continue;
					}
					if(sharedMBVars.find(vIter->first)!=sharedMBVars.end())
					{
						continue;
					}
					double vval=(*evidSet)[vIter->first]->getEvidVal();
					if(sumxy->find(vIter->first)==sumxy->end())
					{	
						(*sumxy)[vIter->first]=gammaval*xpart*vval;
					}
					else
					{
						(*sumxy)[vIter->first]=(*sumxy)[vIter->first]+(gammaval*xpart*vval);
					}
					VSET_ITER uIter=vIter;
					for(;uIter!=potVarSet.end();uIter++)
					{
						if(uIter->first==xId)
						{
							continue;
						}
						if(sharedMBVars.find(uIter->first)!=sharedMBVars.end())
						{
							continue;
						}
						double uval=(*evidSet)[uIter->first]->getEvidVal();
						double uvval=vval*uval*gammaval;
						INTDBLMAP* vrow=NULL;
						INTDBLMAP* urow=NULL;
						if(sumys->find(vIter->first)==sumys->end())
						{
							vrow=new INTDBLMAP;
							(*sumys)[vIter->first]=vrow;
						}
						else
						{
							vrow=(*sumys)[vIter->first];
						}	
						if(vrow->find(uIter->first)==vrow->end())
						{
							(*vrow)[uIter->first]=uvval;
						}
						else
						{
							(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
						}
						if(uIter==vIter)
						{
							continue;
						}
						if(sumys->find(uIter->first)==sumys->end())
						{
							urow=new INTDBLMAP;
							(*sumys)[uIter->first]=urow;
						}
						else
						{
							urow=(*sumys)[uIter->first];
						}
						if(urow->find(vIter->first)==urow->end())
						{
							(*urow)[vIter->first]=uvval;
						}
						else
						{
							(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
						}
					}
				}
			}
		}
	}
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* apot=pIter->second;
		VSET& potVarSet=apot->getAssocVariables();
		int unsharedNeighbors=potVarSet.size()-sharedMBVars.size()-1;
		if(unsharedNeighbors==0)
		{
			continue;
		}
		Matrix* covMatrix=new Matrix(unsharedNeighbors,unsharedNeighbors);
		Matrix* xyMatrix=new Matrix(1,unsharedNeighbors);
		map<int,INTDBLMAP*>* sumys=sumys_set[pIter->first];
		INTDBLMAP* sumxy=sumxy_set[pIter->first];
		map<int,int> vIDMatIDMap;
		int rowvid=0;
		for(map<int,INTDBLMAP*>::iterator vIter=sumys->begin();vIter!=sumys->end();vIter++)
		{
			INTDBLMAP* vrow=vIter->second;
			int rowpos=-1;
			if(vIDMatIDMap.find(vIter->first)==vIDMatIDMap.end())
			{
				vIDMatIDMap[vIter->first]=rowvid;
				rowvid++;
			}
			rowpos=vIDMatIDMap[vIter->first];
			double xyval=(*sumxy)[vIter->first];
			xyMatrix->setValue(xyval,0,rowpos);
			for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
			{
				int colpos=-1;
				if(vIDMatIDMap.find(rIter->first)==vIDMatIDMap.end())
				{
					vIDMatIDMap[rIter->first]=rowvid;
					rowvid++;
				}
				colpos=vIDMatIDMap[rIter->first];
				covMatrix->setValue(rIter->second,rowpos,colpos);
				covMatrix->setValue(rIter->second,colpos,rowpos);
			}
		}
		Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
		Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
		INTDBLMAP& potWt=apot->getCondWeight();
		for(INTINTMAP_ITER vIter=vIDMatIDMap.begin();vIter!=vIDMatIDMap.end();vIter++)
		{
			int pos=vIter->second;
			double wval=wtMatrix->getValue(0,pos);
			apot->setCondWeightFor(vIter->first,wval);
		}
		vIDMatIDMap.clear();
		delete covMatrix;
		delete xyMatrix;
		delete inv;
		delete wtMatrix;
		
		sumxy->clear();
		delete sumxy;
		for(map<int,INTDBLMAP*>::iterator vIter=sumys->begin();vIter!=sumys->end();vIter++)
		{
			vIter->second->clear();
			delete vIter->second;
		}
		sumys->clear();
		delete sumys;
	}
	sumxy_set.clear();
	sumys_set.clear();
	return 0;
}
/*

int
PotentialManager::estimatePotWeight(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTINTMAP& vIDMatIDMap=apot->getVarMatrixIndexMap();
	map<int,INTDBLMAP*> sumys;
	INTDBLMAP sumxy;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double xpart=(*evidSet)[vId]->getEvidVal()-b;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				if(sumxy.find(vIter->first)==sumxy.end())
				{	
					sumxy[vIter->first]=gammaval*xpart*vval;
				}
				else
				{
					sumxy[vIter->first]=sumxy[vIter->first]+(gammaval*xpart*vval);
				}
				VSET_ITER uIter=vIter;
				for(;uIter!=potVarSet.end();uIter++)
				{
					if(uIter->first==vId)
					{
						continue;
					}
					double uval=(*evidSet)[uIter->first]->getEvidVal();
					double uvval=vval*uval*gammaval;
					uvval=uvval+lambda;
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(vIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[vIter->first]=vrow;
					}
					else
					{
						vrow=sumys[vIter->first];
					}	
					if(vrow->find(uIter->first)==vrow->end())
					{
						(*vrow)[uIter->first]=uvval;
					}
					else
					{
						(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
					}
					if(uIter==vIter)
					{
						continue;
					}
					if(sumys.find(uIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[uIter->first]=urow;
					}
					else
					{
						urow=sumys[uIter->first];
					}
					if(urow->find(vIter->first)==urow->end())
					{
						(*urow)[vIter->first]=uvval;
					}
					else
					{
						(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
					}
				}
			}
		}
	}
	Matrix* covMatrix=new Matrix(potVarSet.size()-1,potVarSet.size()-1);
	Matrix* xyMatrix=new Matrix(1,potVarSet.size()-1);
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowvid=vIDMatIDMap[vId];
		int rowpos=vIDMatIDMap[vIter->first];
		if(rowpos>rowvid)
		{
			rowpos--;
		}
	}
	return 0;
}*/


int
PotentialManager::estimatePotCondVarBias(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas, int componentID)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double varnum=0;
	double varden=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			varden=varden+gammaval;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}

	double condBias=biasnum/varden;
	b=condBias;
	varnum=0;

	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
		}
	}
	double condCov=varnum/(varden*(1+lambda));
	if(lambda>0)
	{
		condCov=condCov/2;
	}
	apot->setCondVariance(condCov);
	apot->setUnnormCondVariance(varnum);
	apot->setCondBias(condBias);
	return 0;
}


int
PotentialManager::estimatePotCondVarBias(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double varnum=0;
	double varden=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			varden=varden+gammaval;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}
	double condBias=biasnum/varden;
	b=condBias;
	varnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
		}
	}
	double condCov=varnum/(varden*(1+lambda));
	if(lambda>0)
	{
		condCov=condCov/2;
	}
	apot->setCondVariance(condCov);
	apot->setUnnormCondVariance(varnum);
	apot->setCondBias(condBias);
	return 0;
}


int
PotentialManager::estimatePotTiedVar(map<int,Potential*>& potSet,int vId, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	double varnum=0;
	double varden=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			for(INTDBLMAP_ITER cIter=gMap->begin();cIter!=gMap->end();cIter++)
			{
				Potential* apot=potSet[cIter->first];
				VSET& potVarSet=apot->getAssocVariables();
				INTDBLMAP& condWt=apot->getCondWeight();
				double b=apot->getCondBias();
				double gammaval=(*gMap)[cIter->first];
				double xval=0;
				varden=varden+gammaval;
				EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
				double yval=0;
				double shrinkcorr=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==vId)
					{
						xval=(*evidSet)[vIter->first]->getEvidVal();
						continue;
					}
					double vval=(*evidSet)[vIter->first]->getEvidVal();
					double wt=condWt[vIter->first];
					yval=yval+(vval*wt);
				}
				varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)));
			}
		}
	}
	double condCov=varnum/varden;
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* apot=pIter->second;
		apot->setCondVariance(condCov);
		apot->setUnnormCondVariance(varnum);
	}
	return 0;
}


int
PotentialManager::estimatePotCondBias(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double biasden=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			double gammaval=(*gMap)[cid];
			double xval=0;
			biasden=biasden+gammaval;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
			}
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}
	double condBias=biasnum/biasden;
	apot->setCondBias(condBias);
	return 0;
}

int
PotentialManager::estimatePotCondBiasThenVar(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double den=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			double gammaval=(*gMap)[cid];
			den=den+gammaval;
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
			}
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}


	double condBias=biasnum/den;
	double varnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			double gammaval=(*gMap)[cid];
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-condBias)*(xval-yval-condBias)) + (lambda*shrinkcorr));
		}
	}
	double condCov=varnum/(den*(1+lambda));
	if(lambda>0)
	{
		condCov=condCov/2;
	}
	apot->setCondVariance(condCov);
	apot->setUnnormCondVariance(varnum);
	apot->setCondBias(condBias);
	return 0;
}


double
PotentialManager::getPseudoLikelihood(SlimFactor* sFactor,VSET& varSet, bool train)
{
	Potential* aPotFunc=new Potential;
	Variable* aVar=varSet[sFactor->fId];
	aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
	for(INTINTMAP_ITER aIter=sFactor->mergedMB.begin();aIter!=sFactor->mergedMB.end();aIter++)
	{
		Variable* aVar=varSet[aIter->first];
		aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
	}
	aPotFunc->potZeroInit();
	populatePotential(aPotFunc,false);
	//This function creates a submatrix of the covariance matrix and inverts it
	aPotFunc->initMBCovMean();
	INTINTMAP* dataSet=NULL;
	if(train)
	{
		dataSet=&(evMgr->getTrainingSet());
	}
	else
	{
		dataSet=&(evMgr->getTestSet());
	}
	INTDBLMAP subData;
	double pll=0;
	int thresholded=0;
	for(INTINTMAP_ITER dIter=dataSet->begin();dIter!=dataSet->end();dIter++)
	{
		EMAP* evidMap=NULL;
		evidMap=evMgr->getEvidenceAt(dIter->first);
		Evidence* evid=(*evidMap)[sFactor->fId];
		double val=evid->getEvidVal();
		subData[sFactor->fId]=val;
		for(INTINTMAP_ITER vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=(*evidMap)[vIter->first];
			double val=evid->getEvidVal();
			subData[vId]=val;
		}
		double cll=aPotFunc->getCondPotValueFor(subData);

		if(cll<1e-50)
		{
			cll=1e-50;
			thresholded++;
		}
		pll=pll+log(cll);
	}
	subData.clear();
	if(thresholded>0)
	{
	//	cout <<"Thresholded " << thresholded << " datapoints to 1e-50" << endl;
	}
	delete aPotFunc;
	return pll;
}


double 
PotentialManager::getGaussianLikelihood(map<int,SlimFactor*>& factorSet,VSET& varSet, bool train)
{
	Potential* aPotFunc=new Potential;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		Variable* aVar=varSet[fIter->first];
		if(fIter==factorSet.begin())
		{
			aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		}
		else
		{
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
	}
	aPotFunc->potZeroInit();
	//Use the graph structure to update the particular elements of the covariance and mean of this potential
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;

		double mean=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
		if(globalMean.find(fIter->first)==globalMean.end())
		{
			cerr <<"No var with id " << fIter->first << endl;
			exit(-1);
		}
		mean=globalMean[fIter->first];
		covar=globalCovar[fIter->first];
		aPotFunc->updateMean(fIter->first,mean);
		double vval=(*covar)[fIter->first];
		aPotFunc->updateCovariance(fIter->first,fIter->first,vval);
		for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
		{
			if(covar->find(mIter->first)==covar->end())
			{
				cerr <<"No var " << mIter->first << " in covariance of " << fIter->first << endl;
				exit(-1);
			}
			double cval=(*covar)[mIter->first];
			aPotFunc->updateCovariance(fIter->first,mIter->first,cval);
			aPotFunc->updateCovariance(mIter->first,fIter->first,cval);
		}
	}

	aPotFunc->makeValidJPD();
	INTDBLMAP datapt;
	double gll=0;
	int thresholded=0;
	INTINTMAP* dataSet;
	if(train)
	{
		dataSet=(&evMgr->getTrainingSet());
	}
	else
	{
		dataSet=(&evMgr->getTestSet());
	}
	for(INTINTMAP_ITER dIter=dataSet->begin();dIter!=dataSet->end();dIter++)
	{
		EMAP* evidMap=NULL;
		evidMap=evMgr->getEvidenceAt(dIter->first);
		for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
		{
			int vId=fIter->first;
			Evidence* evid=(*evidMap)[vId];
			double val=evid->getEvidVal();
			datapt[vId]=val;
		}
		double jll=aPotFunc->getJointPotValueFor(datapt);

		if(jll<1e-50)
		{
			jll=1e-50;
			thresholded++;
		}
		gll=gll+log(jll);
	}
	delete aPotFunc;
	return gll;
}


double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}


double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet,map<int,int>& visitedVertices )
{
	double dll=0;
	cout <<"Not implemented" << endl;
	return dll;
}


int
PotentialManager::estimateConditionalPotential(SlimFactor* sFactor,VSET& varSet,Potential** pot, STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int
PotentialManager::populatePotential(Potential* pot,STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}


int 
PotentialManager::estimateCanonicalPotential_Abbeel(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}



int 
PotentialManager::estimateCanonicalPotential_Approximate(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int
PotentialManager::resetPotFuncs()
{
	for(map<int,Potential*>::iterator pIter=potFuncs.begin();pIter!=potFuncs.end();pIter++)
	{
		delete pIter->second;
	}
	potFuncs.clear();
	return 0;
}



int 
PotentialManager::estimateCanonicalPotential_Joint(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

Potential*
PotentialManager::getPotential(int fId)
{
	if(potFuncs.find(fId)==potFuncs.end())
	{
		return NULL;
	}
	return potFuncs[fId];
}


double 
PotentialManager::getConditionalEntropy(int vId,INTINTMAP& fVars,VSET& varSet)
{
	double condEntropy=0;
	//string fullConfStr;
	//string partConfStr;
	char confStr[CONSTR_LEN];
	/*int varCnt=0;
	for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
	{
		sprintf(confStr,"-%d",aIter->first);
		fullConfStr.append(confStr);
		if(aIter->first!=vId)
		{
			partConfStr.append(confStr);
			varCnt++;
		}
	}*/
	double fullJointEntropy=0;
	/*if(jointEntropies.find(fullConfStr)!=jointEntropies.end())
	{
		fullJointEntropy=jointEntropies[fullConfStr];
	}
	else
	{*/
		Potential* potFunc=new Potential;
		for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
		{
			if(aIter==fVars.begin())
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::FACTOR);
			}
			else
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::MARKOV_BNKT);
			}
		}
		potFunc->potZeroInit();
		populatePotential(potFunc,false);
		potFunc->calculateJointEntropy();
		fullJointEntropy=potFunc->getJointEntropy();
		//jointEntropies[fullConfStr]=fullJointEntropy;
		delete potFunc;
	//}
	if(fVars.size()==1)
	{
		/*if(jointEntropies.size()>=20000)
		{
			jointEntropies.clear();
		}*/
		return fullJointEntropy;
	}
	double partJointEntropy=0;
	/*if(jointEntropies.find(partConfStr)!=jointEntropies.end())
	{
		partJointEntropy=jointEntropies[partConfStr];
	}
	else
	{
		Potential* potFunc=new Potential;*/
		potFunc=new Potential;
		bool setFactorVar=false;
		for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
		{
			if(aIter->first==vId)
			{
				continue;
			}
			if(!setFactorVar)
			{
				setFactorVar=true;
				potFunc->setAssocVariable(varSet[aIter->first],Potential::FACTOR);
			}
			else
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::MARKOV_BNKT);
			}
		}
		STRDBLMAP counts;
		potFunc->potZeroInit();
		populatePotential(potFunc,false);
		potFunc->calculateJointEntropy();
		partJointEntropy=potFunc->getJointEntropy();
		//jointEntropies[partConfStr]=partJointEntropy;
		delete potFunc;
	//}
	condEntropy=fullJointEntropy-partJointEntropy;
	//fullConfStr.clear();
	//partConfStr.clear();
	/*if(jointEntropies.size()>=20000)
	{
		jointEntropies.clear();
	}*/
	return condEntropy;
}


double 
PotentialManager::getSampleLikelihood(map<int,SlimFactor*>& factorSet, VSET& varSet, INTINTMAP* sample)
{
	double sampleLL=0;
	cout <<"Not implemented " <<endl;
	return sampleLL;
}


int 
PotentialManager::getVariableSample(INTINTMAP& jointConf,VSET& varSet,int vId,SlimFactor* sFactor, gsl_rng* r)
{
	Potential* pot=NULL;
	if(potFuncs.find(vId)==potFuncs.end())
	{
		pot=new Potential;
		STRDBLMAP counts;
		estimateConditionalPotential(sFactor,varSet,&pot,counts);
		potFuncs[vId]=pot;
	}
	else
	{
		pot=potFuncs[vId];
	}
	//int sample=pot->generateSample(jointConf,vId,r);
	int sample=-1;
	return sample;
}

int
PotentialManager::clearJointEntropies()
{
	jointEntropies.clear();
	return 0;
}

double
PotentialManager::estimateCanonicalValue(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* condPot)
{
	double pVal=0;
	cout <<"Not implemented " << endl;
	return 0;
}


double
PotentialManager::estimateCanonicalValue_Joint(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* jointPot)
{
	cout <<"Not implemented " << endl;
	return 0;
}
