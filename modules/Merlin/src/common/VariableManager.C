#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"


VariableManager::VariableManager()
{
}

VariableManager::~VariableManager()
{
}

//Reads the schema of the variables

Error::ErrorCode
VariableManager::readVariables(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[400000];
	int nodeCnt=0;

	if(inFile.good())
	{
		inFile.getline(buffer,400000);

		if(strlen(buffer)<=0)
		{
			cout <<"Bad format" << endl;
			return Error::VARSCHEMA_ERR;
		}

		char* tok=strtok(buffer,"\t");
		int tokCnt=0;

		while(tok!=NULL)
		{
			Variable* var=new Variable;
			var->setID(tokCnt);
			var->setName(tok);
			variableSet[tokCnt]=var;
			
			string varKey(tok);
			varNameIDMap[varKey]=tokCnt;
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
	}

	inFile.close();

	cout <<"Read information about " << variableSet.size() << " variables " << endl;

	return Error::SUCCESS;
}

int
VariableManager::getVarID(const char* varName)
{
	string varKey(varName);
	if(varNameIDMap.find(varKey)==varNameIDMap.end())
	{
		return -1;
	}
	int vId=varNameIDMap[varKey];
	return vId;
}

bool 
VariableManager::isValid(int varID,int varVal)
{
	Variable* rVar=variableSet[varID];
	return rVar->isValidValue(varVal);
}

map<int,Variable*>& 
VariableManager::getVariableSet()
{
	return variableSet;
}


Variable* 
VariableManager::getVariableAt(int vId)
{
	if(variableSet.find(vId)==variableSet.end())
	{
		cout << "Illegal variable id " << vId << endl;
		return NULL;
	}
	return variableSet[vId];
}
