#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
using namespace std;
#include "MotifManager.H"



MotifManager::MotifManager()
{
}

MotifManager::~MotifManager()
{
}

int 
MotifManager::readMotifs(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	string buffstr;
	int bufflen=0;
	int lineCnt=0;
	int motifid=0;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		if(lineCnt<=0)
		{
			lineCnt++;
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=buffstr.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		char* tok=buffer;
		int tokCnt=0;
		string motifName;
		string tfName;
		string orfName;
		string motifInterp;
		double hitcnt=0;
		while(tok!=NULL)
		{
			char* end=strchr(tok,'\t');
			if(end!=NULL)
			{
				*end='\0';
			}
			if(tokCnt==0)
			{
				motifName.append(tok);	
				tfName.append(tok);
			}
			else if(tokCnt==2)
			{
				char inter_noSpace[1024];
				int id=0;
				while(tok[id]!='\0')
				{
					if(tok[id]==' ')
					{
						inter_noSpace[id]='_';
					}
					else
					{
						inter_noSpace[id]=tok[id];
					}
					id++;
				}
				inter_noSpace[id]='\0';
				motifInterp.append(inter_noSpace);
			}
			else if(tokCnt==4)
			{
				hitcnt=atof(tok);
			}
			else if(tokCnt==6)
			{
				orfName.append(tok);
			}
			if(end!=NULL)
			{
				tok=end+1;
			}
			else
			{
				tok=end;
			}
			tokCnt++;
		}
		int currmid=-1;
		if(motifIDMap.find(motifName)==motifIDMap.end())
		{
			motifIDMap[motifName]=motifid;
			idmotifMap[motifid]=motifName;
			motifIDInterpret[motifid]=motifInterp;
			currmid=motifid;
			motifid++;
		}
		else
		{
			currmid=motifIDMap[motifName];
		}
		INTDBLMAP* mProfile=NULL;
		if(motifProfileSet.find(orfName)==motifProfileSet.end())
		{
			mProfile=new INTDBLMAP;
			motifProfileSet[orfName]=mProfile;
		}
		else
		{
			mProfile=motifProfileSet[orfName];
		}
		(*mProfile)[currmid]=hitcnt;
	}
	inFile.close();
	return 0;
}

int
MotifManager::readMotifTFMap(const char* tfMotifFName)
{
	ifstream inFile(tfMotifFName);
	char buffer[1024];
	string motifName;
	STRINTMAP* tfNames=NULL;
	while(inFile.good())
	{
		motifName.clear();
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				motifName.append(tok);
			//	cout << "Analyzing " << motifName << endl;
				if(strcmp(motifName.c_str(),"M[56]")==0)
				{
					cout << "Stop here" << endl;
				}
				tfNames=new STRINTMAP;
			}
			else if(tokCnt>1)
			{
				string atf(tok);
				(*tfNames)[atf]=0;
				STRINTMAP* motifsForTF=NULL;
				if(tfMotifMap.find(atf)==tfMotifMap.end())
				{
					motifsForTF=new STRINTMAP;
					tfMotifMap[atf]=motifsForTF;
				}
				else
				{
					motifsForTF=tfMotifMap[atf];
				}
				(*motifsForTF)[motifName]=0;
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		motifTFMap[motifName]=tfNames;
	}
	inFile.close();
	return 0;
}

INTDBLMAP* 
MotifManager::getMotifProfile(const char* orfName)
{
	string orfKey(orfName);
	if(motifProfileSet.find(orfKey)==motifProfileSet.end())
	{
		return NULL;
	}
	return motifProfileSet[orfKey];
}


map<string,INTDBLMAP*>& 
MotifManager::getMotifProfileSet()
{
	return motifProfileSet;
}

map<int,string>&
MotifManager::getMotifIDMap()
{
	return idmotifMap;
}

map<string,int>&
MotifManager::getMotifNameIDMap()
{
	return motifIDMap;
}

map<int,string>&
MotifManager::getMotifIDInterpretation()
{
	return motifIDInterpret;
}




double 
MotifManager::getMotifVal(const string& genename,int mId)
{
	double instanceVal=-1;
	if(motifProfileSet.find(genename)==motifProfileSet.end())
	{
		return instanceVal;
	}
	INTDBLMAP* motifProf=motifProfileSet[genename];
	instanceVal=0;
	if(motifProf->find(mId)!=motifProf->end())
	{
		instanceVal=(*motifProf)[mId];
	//	instanceVal=1;
	}
	return instanceVal;
}


map<string,STRINTMAP*>& 
MotifManager::getMotifTFMap()
{
	return motifTFMap;
}

map<string,STRINTMAP*>&
MotifManager::getTFMotifMap()
{
	return tfMotifMap;
}
