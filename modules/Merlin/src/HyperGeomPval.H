#ifndef _HYPERGEOM_PVAL
#define _HYPERGEOM_PVAL
class HyperGeomPval
{
	public:
		HyperGeomPval();
		~HyperGeomPval();
		//Here we ask, what is the probability of getting k out x elements
		//to belong to category A, given there n1 elements in category A
		//and n2 elements is category B.
		double getOverRepPval(int x, int k, int n1, int n2);
		double getUnderRepPval(int x, int k, int n1, int n2);
	
};

#endif
