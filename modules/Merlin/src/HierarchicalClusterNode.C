#include <iostream>
#include "HierarchicalClusterNode.H"

HierarchicalClusterNode::HierarchicalClusterNode()
{
	left=NULL;
	right=NULL;
	parent=NULL;
	status=0;
	size=1;
}

HierarchicalClusterNode::~HierarchicalClusterNode()
{
	//cout <<"Called the destructor of " << nodeName << endl;
	distToNeighbors.clear();
	distToNeighbors_CC.clear();
	expr.clear();
	attrib.clear();
}
