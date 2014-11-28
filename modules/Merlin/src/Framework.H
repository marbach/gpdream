/*The class which brings it all together*/
#ifndef _FRAMEWORK
#define _FRAMEWORK

//All the managers need forward declarations
class MetaLearner;
class VariableManager;
class EvidenceManager;
class Framework
{
	public:
		Framework();
		~Framework();
		Error::ErrorCode init(int, char**);
		int start();

	private:
		MetaLearner metaLearner;
		VariableManager varManager;
		EvidenceManager evManager;
		EvidenceManager holdoutEvManager;
		double epsThreshold;
		char actionStr[256];
		int cvCnt;
		map<int,int> idVidMap;
};
#endif
