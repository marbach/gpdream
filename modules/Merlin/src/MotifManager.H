#ifndef _MOTIFREADER_
#define _MOTIFREADER_

#include <map>
#include <string>
#include "CommonTypes.H"
class MotifManager
{
	public:
	MotifManager();
	~MotifManager();
	int readMotifs(const char*);
	int readMotifTFMap(const char*);
	INTDBLMAP* getMotifProfile(const char*);
	map<string,INTDBLMAP*>& getMotifProfileSet();
	map<int,string>& getMotifIDMap();
	map<string,int>& getMotifNameIDMap();
	map<int,string>& getMotifIDInterpretation();
	double getMotifVal(const string&,int);
	map<string,STRINTMAP*>& getMotifTFMap();
	map<string,STRINTMAP*>& getTFMotifMap();
	private:
	map<string,INTDBLMAP*> motifProfileSet;
	map<string,int> motifIDMap;
	map<int,string> idmotifMap;
	map<int,string> motifIDInterpret;
	map<string,STRINTMAP*> motifTFMap;
	map<string,STRINTMAP*> tfMotifMap;
};
#endif
