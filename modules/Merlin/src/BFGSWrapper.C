#include <iostream>
#include <fstream>
using namespace std;

#include "gsl/gsl_randist.h"
#include "Error.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"


BFGSWrapper::BFGSWrapper()
{
	initx=NULL;
	optimizer=NULL;
	optimizer_func=NULL;
}

BFGSWrapper::~BFGSWrapper()
{
	if(initx!=NULL)
	{
		gsl_vector_free(initx);
	}
	if(optimizer!=NULL)
	{
		gsl_multimin_fdfminimizer_free(optimizer);
	}
}

int
BFGSWrapper::setParamCnt(int aCnt)
{
	paramCnt=aCnt;
	return 0;
}


int
BFGSWrapper::setFeatureCnt(int aCnt)
{
	featureCnt=aCnt;
	return 0;
}

int
BFGSWrapper::setStepSize(double asize)
{
	stepSize=asize;
	return 0;
}

int
BFGSWrapper::setTolerance(double tol)
{
	tolerance=tol;
	return 0;
}

int
BFGSWrapper::setC(double cval)
{
	C=cval;
	return 0;
}

gsl_vector*
BFGSWrapper::getInitVector()
{
	return initx;
}

int
BFGSWrapper::setInitVector(gsl_vector* x)
{
	initx=x;
	return 0;
}

map<int,map<int,double>*>& 
BFGSWrapper::getData()
{
	return data;
}

map<int,double>&
BFGSWrapper::getThetas()
{
	return thetas;
}

int 
BFGSWrapper::initializeMinimizer()
{
	bwdata.data=&data;
	bwdata.thetas=&thetas;
	bwdata.C=C;
	const gsl_multimin_fdfminimizer_type* T=gsl_multimin_fdfminimizer_vector_bfgs2;
	//const gsl_multimin_fdfminimizer_type* T=gsl_multimin_fdfminimizer_conjugate_fr;
	optimizer=gsl_multimin_fdfminimizer_alloc(T,paramCnt);
	optimizer_func=new gsl_multimin_function_fdf;
	optimizer_func->n=paramCnt;
	optimizer_func->f=&likelihood_function;
	optimizer_func->df=&likelihood_gradient;
	optimizer_func->fdf=&likelihood_function_gradient;
	optimizer_func->params=(void*)(&bwdata);
	//gsl_vector* initx=getStartingPoint();
	gsl_multimin_fdfminimizer_set(optimizer,optimizer_func,initx,stepSize,tolerance);
	return 0;
}

int
BFGSWrapper::reinitializeMinimizer()
{
	//gsl_multimin_fdfminimizer_restart(optimizer);
	gsl_multimin_fdfminimizer_set(optimizer,optimizer_func,optimizer->x,stepSize,tolerance);
	return 0;
}

int
BFGSWrapper::optimize()
{
	int iter=0;
	int status;
	do
	{
		iter++;
		status=gsl_multimin_fdfminimizer_iterate(optimizer);
		cout <<"Curr Multimin_Iter " << iter << " function_val " << optimizer->f << endl;
		if(status)
		{
			cout <<"Breaking: gsl_errno " << status << endl;
			break;
		}
		status=gsl_multimin_test_gradient(optimizer->gradient,1e-3);
		if(status==GSL_SUCCESS)
		{
			cout <<"Minima found "<< iter << endl;
		}
	} while(status==GSL_CONTINUE && iter<50);
	gsl_vector* wt=optimizer->x;
	for(int i=0;i<paramCnt;i++)
	{
	      		cout << i << " "<< gsl_vector_get(wt,i) << endl;
	}
	cout <<"function val " << optimizer->f << endl;
	return status;
}


gsl_vector* 
BFGSWrapper::getParams()
{
	gsl_vector* params=optimizer->x;
	return params;
}

double 
BFGSWrapper::getOptimalFval()
{
	double f=optimizer->f;
	f=-1*f;
	return f;
}

double 
likelihood_function(const gsl_vector* x, void* params)
{
	//cout <<"In likelihood function " << endl;
	double ans=0;
	double den=0;
	//The vector looks as follows
	//The first k dimensions have the gammas
	//The remaining dimensions have the motifs
	double gvals[MAXCOMP];
	BFGSWrapperData bwdata=*((BFGSWrapperData*) params);
	map<int,map<int,double>*>& data=*(bwdata.data);
	map<int,double>&thetas =*(bwdata.thetas);
	double C=bwdata.C;
	//cout <<"Analyzing data of size " << data.size() << " with gammas " << gammas.size() << endl;
	for(map<int,map<int,double>*>::iterator dIter=data.begin();dIter!=data.end();dIter++)
	{
		map<int,double>* dataPt=dIter->second;
		double theta_ij=thetas[dIter->first];
		double theta_sq=theta_ij*theta_ij;
		int featureCnt=dataPt->size();
		//cout <<"Analyzing datapoint " << dIter->first << endl;
		double sum=0;
		for(map<int,double>::iterator fIter=dataPt->begin();fIter!=dataPt->end();fIter++)
		{
			double val=gsl_vector_get(x,fIter->first);
			double dataPtVal=fIter->second;
			sum=sum+(val*dataPtVal);
			if(isnan(val))
			{
				cout <<"val is nan " << val << " for featureid " << fIter->first << endl;
			}
		}
		//ans=ans+(theta_sq/(2*sum))+(C*sum)+(0.5*log(2*3.14));
		ans=ans+(theta_sq/(2*sum))+(C*sum);
		if(isnan(ans))
		{
			cout <<"ans is nan  datapt:" << dIter->first<< endl;
		}
	}
//	cout <<"Current fval " << ans << endl;
	return ans;
}

//Store the gradient of the function in g using the value of x
void 
likelihood_gradient(const gsl_vector* x, void* params,gsl_vector* g)
{
	//cout <<"In likelihood_gradient function " << endl;
	//Iterate over the complete data to get the common parts first
	map<int,double> sum1;
	BFGSWrapperData bwdata=*((BFGSWrapperData*) params);
	map<int,map<int,double>*>& data=*(bwdata.data);
	map<int,double>&thetas=*(bwdata.thetas);
	double C=bwdata.C;
	for(map<int,map<int,double>*>::iterator dIter=data.begin();dIter!=data.end();dIter++)
	{
		map<int,double>* dataPt=dIter->second;
		int featureCnt=dataPt->size();
		double numsum=0;
		double densum=0;
		double theta_ij=thetas[dIter->first];
		double theta_sq=theta_ij*theta_ij;
		double sum=0;
		for(map<int,double>::iterator fIter=dataPt->begin();fIter!=dataPt->end();fIter++)
		{
			double wx=gsl_vector_get(x,fIter->first);
			double dval=fIter->second;
			sum=sum+(wx*dval);
		}
		double aval=C-(theta_sq/(2*sum*sum));
		for(map<int,double>::iterator fIter=dataPt->begin();fIter!=dataPt->end();fIter++)
		{
			double dval=fIter->second;
			double scval=dval*aval;
			if(sum1.find(fIter->first)==sum1.end())
			{
				sum1[fIter->first]=scval;
			}
			else
			{
				sum1[fIter->first]=sum1[fIter->first]+scval;
			}
		}
	}
	for(map<int,double>::iterator fIter=sum1.begin();fIter!=sum1.end();fIter++)
	{
		double val=fIter->second;
		gsl_vector_set(g,fIter->first,val);
	}
	sum1.clear();
		
}

void
likelihood_function_gradient(const gsl_vector* x, void* params, double* f, gsl_vector* g)
{
	*f=likelihood_function(x,params);
	likelihood_gradient(x,params,g);
}



//Get the starting point for the parameters that we are optimizing
gsl_vector* 
BFGSWrapper::getStartingPoint()
{
	gsl_vector* initx=gsl_vector_alloc(paramCnt);
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	for(int p=0;p<paramCnt;p++)
	{
		double w_p=gsl_ran_flat(r,0,1);
		gsl_vector_set(initx,p,w_p);
	}
	return initx;
}
