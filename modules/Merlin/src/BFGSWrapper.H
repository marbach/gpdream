#ifndef _BFGS_WRAPPER_
#define _BFGS_WRAPPER_

#include <string>
#include <map>
using namespace std;
#include "gsl/gsl_multimin.h"

#define MAXCOMP 40


class BFGSWrapper
{
	public:
		BFGSWrapper();
		~BFGSWrapper();
		int setParamCnt(int aCnt);
		int setFeatureCnt(int aCnt);
		int setStepSize(double);
		int setTolerance(double);
		int setC(double);
		gsl_vector* getInitVector();
		int setInitVector(gsl_vector*);
		//Somebody populate the data and gammas
		map<int,map<int,double>*>& getData();
		map<int,double>& getThetas();
		int initializeMinimizer();
		int reinitializeMinimizer();
		int optimize();
		gsl_vector* getParams();
		double getOptimalFval();
	private:
		gsl_vector* getStartingPoint();
		int paramCnt;
		int featureCnt;
		double stepSize;
		double tolerance;
		gsl_vector* initx;
		map<int,double> thetas;
		map<int,map<int,double>*> data;
		double optimalFval;
		BFGSWrapperData bwdata;
		gsl_multimin_fdfminimizer* optimizer;
		gsl_multimin_function_fdf* optimizer_func;
		double C;
};
#endif

extern double likelihood_function(const gsl_vector* x, void* params);
extern void likelihood_gradient(const gsl_vector* v, void* params,gsl_vector* df);
extern void likelihood_function_gradient(const gsl_vector* x, void* params, double* f, gsl_vector* g);
