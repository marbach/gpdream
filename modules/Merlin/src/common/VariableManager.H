#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include <string>
#include <map>

using namespace std;

class Error;
class Variable;
class VariableManager
{
	public:
		VariableManager();
		~VariableManager();
		Error::ErrorCode readVariables(const char*);
		int getVarID(const char*);
		bool isValid(int,int);
		map<int,Variable*>& getVariableSet();
		Variable* getVariableAt(int);
	private:
		map<string,int> varNameIDMap;
		map<int,Variable*> variableSet;
		int globalVarId;
};
#endif
