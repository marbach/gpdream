#ifndef _HIERARCHICAL_CLUSTER_NODE
#define _HIERARCHICAL_CLUSTER_NODE
#include <map>
#include <string>
#include <vector>
using namespace std;
class HierarchicalClusterNode;
class HierarchicalClusterNode
{
	public:
	HierarchicalClusterNode();
	~HierarchicalClusterNode();
	
	map<int,double> attrib;
	map<string,double> distToNeighbors;
	map<string,double> distToNeighbors_CC;
	vector<double> expr;
	HierarchicalClusterNode* left;
	HierarchicalClusterNode* right;
	HierarchicalClusterNode* parent;
	int id;
	string nodeName;
	int status;
	int size;
};
#endif
