#include <map>
using namespace std;

typedef struct
{
	map<int,map<int,double>*>* data;
	map<int,double>* thetas;
	double C;
}BFGSWrapperData;
