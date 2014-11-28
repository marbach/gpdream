#include <iostream>
using namespace std;

#include <math.h>
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "HyperGeomPval.H"


HyperGeomPval::HyperGeomPval()
{
}

HyperGeomPval::~HyperGeomPval()
{

}

/*
double 
HyperGeomPval::getOverRepPval(int x, int n, double p)
{
	double pvalue=0;
	double initX=x;
	while(initX<=n)
	{
		double occProb=pow(p,initX);
		double nonOccProb=pow(1-p,n-initX);
		if((occProb==0) || (nonOccProb==0))
		{
			pvalue=1;
			break;
		}
		double big=0;
		double small=0;
		//need to compute n!/((initX!)*(n-initX)!)
		if(initX >= (n-initX))
		{
			big=initX;
			small=n-initX;
		}
		else
		{
			big=n-initX;
			small=initX;
		}

		double currNat=n;
		double temp1=occProb*nonOccProb;
		temp1=log(temp1);
		while(currNat>big)
		{
			//cout<< temp1 << endl;
			temp1=temp1+log(currNat);
			currNat--;
		}
		double temp2=1;
		temp2=log(temp2);
		double currNat1=small;
		while(currNat1>0)
		{
			temp2=temp2+log(currNat1);
			currNat1--;
		}
		double nPerms=exp(temp1-temp2);
		pvalue=pvalue+nPerms;
		initX++;
	}
	return pvalue;
}*/


double 
HyperGeomPval::getOverRepPval(int t, int k, int n1, int n2)
{
	double pvalue=0;
	int initX=k;
	while(initX<=t)
	{
		double temp1=gsl_ran_hypergeometric_pdf(initX,n1,n2,t);
		pvalue=temp1+pvalue;
		initX++;
	}
	return pvalue;
}

double
HyperGeomPval::getUnderRepPval(int t, int k, int n1, int n2)
{
	double pvalue=0;
	int initX=0;
	while(initX<=k)
	{
		double temp1=gsl_ran_hypergeometric_pdf(initX,n1,n2,t);
		pvalue=pvalue+temp1;
		initX++;
	}
	return pvalue;
}

/*
double
HyperGeomPval::getUnderRepPval(int x, int n, double p)
{
	double pvalue=0;
	double initX=0;
	while(initX<=x)
	{
		double occProb=pow(p,initX);
		double nonOccProb=pow(1-p,n-initX);
		if((occProb==0) || (nonOccProb==0))
		{
			pvalue=1;
			break;
		}
		double big=0;
		double small=0;
		//need to compute n!/((initX!)*(n-initX)!)
		if(initX >= (n-initX))
		{
			big=initX;
			small=n-initX;
		}
		else
		{
			big=n-initX;
			small=initX;
		}

		double currNat=n;
		double temp1=occProb*nonOccProb;
		temp1=log(temp1);
		while(currNat>big)
		{
			temp1=temp1 + log(currNat);
			currNat--;
		}
		double temp2=1;
		temp2=log(temp2);
		double currNat1=small;
		while(currNat1>0)
		{
			temp2=temp2 + log(currNat1);
			currNat1--;
		}
		double nPerms=exp(temp1-temp2);
		pvalue=pvalue+nPerms;
		initX++;
	}
	return pvalue;
}*/
