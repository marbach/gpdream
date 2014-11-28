#ifndef _POTENTIAL
#define _POTENTIAL
//A potential over continuous variables needs only store the mean and variance parameters.
//However to compute the pdf, we need to store the inverse of the covariance matrix and the
//determinant.
#include <fstream>
#define SMALL_PROB 1e-60
#define PI 3.1472
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#include "Matrix.H"
class Evidence;
class Variable;
class Potential
{
	public:
		Potential();
		~Potential();
		
		int resetVarSet();
		int initMatrix(int);
		//This enum is used to decide what role the Variable
		//plays in the potential. A Variable can be either from a Markov Blanket
		//or from the factor. This is required while converting the joint into a conditional
		//of the factor variables give Markov blanket variables
		
		typedef enum
		{
			FACTOR,
			MARKOV_BNKT
		}VariableRole;
		
		int setAssocVariable(Variable*,VariableRole);
		VSET& getAssocVariables();
		INTINTMAP& getVarMatrixIndexMap(); 
		int updateFactorVariable(int);
		int setSharedMBVar(int);
		INTINTMAP& getSharedMBVars();

		//Initialize potential 
		int potZeroInit();
		int potZeroInit_MeanOnly();
		int potCustomInit(double,double);

		double getJointPotValueFor(INTDBLMAP&);
		//double getCondPotValueFor(int vId, INTDBLMAP&);
		double getJointPotValueFor(STRINTMAP&);
		double getJointPotValueForConf(string&);
		
		//The result of marginalization is put in to the second
		//parameter
		int copyMe(Potential**);
		//Needed when we are updating the joint potential
		int updateMean(int,double);
		int updateCovariance(int,int,double);
		int dumpPotential(ostream&);
		int showVariables(ostream&);
		int makeValidJPD();
		int makeValidJPD(gsl_matrix*,gsl_permutation*);

		//Since both the joint and the conditional are in the same potential 
		//we can easily compute the entropy
		int calculateEntropy();
		int calculateJointEntropy();
		double getEntropy();
		double getJointEntropy();
		double generateSample(INTDBLMAP&, int vId,gsl_rng*);
		double generateSample(INTDBLMAP&, int vId,gsl_rng*,double);
		double predictSample(INTDBLMAP&, int vId);

		int initMBCovMean();
		double getCondPotValueFor(INTDBLMAP&);
		double getCondPotValueFor(map<int,Evidence*>*);
		double getCondPotValueFor(map<int,Evidence*>*,int);
		double getCondPotValueFor(map<int,Evidence*>* evidMap,int newMBVar,map<int,double>& partialMean,int evidID);
		int setCondVariance(double);
		double getCondVariance();
		int setUnnormCondVariance(double);
		double getUnnormCondVariance();
		int setCondBias(double);
		double getCondBias();
		int setCondWeight(Matrix*);
		int setCondWeight(INTDBLMAP&);
		int setCondWeightFor(int,double);
		INTDBLMAP& getCondWeight();
		int showParams();
		double predictSample(map<int,Evidence*>*);
		double predictSample(map<int,Evidence*>*,int&);
		INTDBLMAP& getAllPartialMeans();
		double computeLL_Tracetrick(int sampleSize);
	protected:
		VSET varSet;
		INTINTMAP factorVariables;
		INTINTMAP markovBlnktVariables;
		INTDBLMAP meanWrap;
		Matrix* covariance;
		Matrix* inverse;
		Matrix* mean;
		double determinant;
		double normFactor;
		double potEntropy;
		double jointEntropy;
		//Each variable is mapped to an index in the covariance or mean matrix
		//This map stores that information
		int matrixInd;
		INTINTMAP vIDMatIndMap;
		INTINTMAP matIndvIDMap;

		double mbcondVar;
		double mbcondVar_Unnorm;
		double mbcondMean_Part;
		INTDBLMAP mbcondMean_Vect;
		INTINTMAP sharedMBVars;
		INTDBLMAP cachePartialMeans;
};
#endif
