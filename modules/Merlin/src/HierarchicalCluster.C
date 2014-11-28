#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>
#include "Error.H"
#include "Distance.H"
#include "Variable.H"
#include "VariableManager.H"
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
int sortfunc(const void* first, const void* second);
int* sortingind=NULL;
double* sorteddist=NULL;


HierarchicalCluster::HierarchicalCluster()
{
	root=NULL;
}

HierarchicalCluster::~HierarchicalCluster()
{
}

int
HierarchicalCluster::setVariableManager(VariableManager* p)
{
	vMgr=p;
	return 0;
}

int
HierarchicalCluster::setOutputDir(const char* aFName)
{
	strcpy(outputDir,aFName);
	return 0;
}

int 
HierarchicalCluster::cluster(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double t)
{
	threshold=t;
	map<int,HierarchicalClusterNode*> attribs_integer;
	globalNodeID=0;
	for(map<string,HierarchicalClusterNode*>::iterator nIter=attribs.begin();nIter!=attribs.end();nIter++)
	{
		attribs_integer[globalNodeID]=nIter->second;
		nIter->second->id=globalNodeID;
		globalNodeID++;
	}
	estimatePairwiseDist(attribs_integer);
	for(map<string,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		backup[aIter->first]=aIter->second;
	}
	
	bool keepMerging=true;
	
	while(keepMerging && attribs_integer.size()>1)
	{
		struct timeval begintime;
		struct timeval endtime;
		gettimeofday(&begintime,NULL);
		keepMerging=mergePairs_LazyDelete(attribs_integer);
		gettimeofday(&endtime,NULL);
		//cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
	}
	generateModules(attribs_integer,modules,backup);
	calculatePercentVarianceExplained(modules,backup);
	//calculateSilhouetteIndex(modules,backup);
	for(map<string,HierarchicalClusterNode*>::iterator aIter=backup.begin();aIter!=backup.end();aIter++)
	{
		map<string,HierarchicalClusterNode*>::iterator cIter=attribs.find(aIter->first);
		if(cIter!=attribs.end())
		{
			attribs.erase(cIter);
		}
		delete aIter->second;
	}
	for(map<string,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		delete aIter->second;
	}
	for(map<int,HierarchicalClusterNode*>::iterator aIter=attribs_integer.begin();aIter!=attribs_integer.end();aIter++)
	{
	//	delete aIter->second;
	}
	attribs_integer.clear();
	delete [] visited;
	for(int i=0;i<treenodecnt;i++)
	{
		delete [] distvalues[i];
	}
	delete [] distvalues;
	attribs.clear();
	backup.clear();
	while(!myqueue.empty())
	{
		HierarchicalCluster::Pair* p=myqueue.top();
		myqueue.pop();
		delete p;
	}
	return 0;
}

int
HierarchicalCluster::estimatePairwiseDist(map<int,HierarchicalClusterNode*>& currNodeSet)
{
	Distance d;
	int size=currNodeSet.size();
	size=size*(size-1)/2;
	//The total number of nodes that can be there in a hierarchical cluster is 2n-1
	treenodecnt=(currNodeSet.size()*2) - 1;
	distvalues=new double*[treenodecnt];
	visited=new int[treenodecnt];
	for(int i=0;i<treenodecnt;i++)
	{
		distvalues[i]=new double[treenodecnt];
		for (int j=0;j<treenodecnt;j++)
		{
			distvalues[i][j]=-1000;
		}
	}
	for(int i=0;i<treenodecnt;i++)
	{
		visited[i]=0;
	}
	int globalPairID=0;
	//for(map<int,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	for(int i=0;i<currNodeSet.size();i++)
	{
		for(int j=i+1;j<currNodeSet.size();j++)
		{
			HierarchicalClusterNode* hcNode1=currNodeSet[i];
			HierarchicalClusterNode* hcNode2=currNodeSet[j];
			double rdist=0;
			double ccdist=d.computeCC(hcNode1->expr,hcNode2->expr);
			double sharedSign=0;
			double den1=0;
			double den2=0;
			ccdist=0.5*(1-ccdist);
			for(map<int,double>::iterator aIter=hcNode1->attrib.begin();aIter!=hcNode1->attrib.end();aIter++)
			{
				den1=den1+fabs(aIter->second);
				if(hcNode2->attrib.find(aIter->first)!=hcNode2->attrib.end())
				{	
					if((aIter->second*hcNode2->attrib[aIter->first])>=0)
					{
						sharedSign=sharedSign+(((fabs(aIter->second)+fabs(hcNode2->attrib[aIter->first])))/2.0);
					}
				}
			}
			for(map<int,double>::iterator aIter=hcNode2->attrib.begin();aIter!=hcNode2->attrib.end();aIter++)
			{
				den2=den2+fabs(aIter->second);
			}
			rdist=1- (((double)sharedSign)/((double)(den1+den2-sharedSign)));
			double dist=(ccdist+rdist)/2;
			distvalues[i][j]=dist;
			distvalues[j][i]=dist;
			HierarchicalCluster::Pair* p=new HierarchicalCluster::Pair;
			p->node1=i;
			p->node2=j;
			p->value=dist;
			myqueue.push(p);
			globalPairID++;
		}
	}
	return 0;
}


int
HierarchicalCluster::mergePairs_LazyDelete(map<int,HierarchicalClusterNode*>& currNodeSet)
{
	double maxSum=0;
	double minDist=100000;
	if(myqueue.empty())
	{
		return false;
	}
	HierarchicalCluster::Pair* p=myqueue.top();
	if(p->value>=threshold)
	{
		return false;
	}
	//Keep popping until we reach a pair whose both members have not been visited
	while(!myqueue.empty() && (visited[p->node1]==1 || visited[p->node2]==1))
	{
		delete p;
		myqueue.pop();
		p=myqueue.top();
	}
	minDist=p->value;
	
	visited[p->node1]=1;
	visited[p->node2]=1;
	if(currNodeSet.find(p->node1)==currNodeSet.end())
	{
		cerr << "node1 " << p->node1 << " not found in pair " << p->node1<<"-" <<p->node2 << endl;
		exit(-1);
	}
	if(currNodeSet.find(p->node2)==currNodeSet.end())
	{
		cerr << "node2 " << p->node2 << " not found in pair " << p->node1 <<"-" << p->node2<<  endl;
		exit(-1);
	}
	HierarchicalClusterNode* c1=currNodeSet[p->node1];
	HierarchicalClusterNode* c2=currNodeSet[p->node2];
	//cout << "Merging " << c1->nodeName << " " << c2->nodeName  << " dist=" << minDist << endl;
	HierarchicalClusterNode* c12=new HierarchicalClusterNode;
	c12->left=c1;
	c12->right=c2;
	c1->parent=c12;
	c2->parent=c12;
	c12->nodeName.append(c1->nodeName);
	c12->nodeName.append("-");
	c12->nodeName.append(c2->nodeName);
	map<int,HierarchicalClusterNode*>::iterator hIter1=currNodeSet.find(c1->id);
	map<int,HierarchicalClusterNode*>::iterator hIter2=currNodeSet.find(c2->id);
	currNodeSet.erase(hIter1);
	currNodeSet.erase(hIter2);
	double* dist_n1=distvalues[c1->id];
	double* dist_n2=distvalues[c2->id];
	double* dist_n12=distvalues[globalNodeID];
	for(map<int,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	{
		if(visited[nIter->first]==1)
		{
			continue;
		}
		double d1=dist_n1[nIter->first];
		double d2=dist_n2[nIter->first];
		double dkm_rdist=((c1->size*d1) + (c2->size*d2))/((double)(c1->size+ c2->size));
		double dist=dkm_rdist;
		dist_n12[nIter->first]=dist;
		double* dist_other=distvalues[nIter->first];
		dist_other[globalNodeID]=dist;
		if(globalPairID==14412)
		{	
			cout <<"Stop here " << endl;
		}
		HierarchicalCluster::Pair* newp=new HierarchicalCluster::Pair;
		newp->value=dist;
		if(globalNodeID<nIter->first)
		{
			newp->node1=globalNodeID;
			newp->node2=nIter->first;
		}
		else
		{
			newp->node1=nIter->first;
			newp->node2=globalNodeID;
		}
		myqueue.push(newp);
		globalPairID++;
	}
	c12->id=globalNodeID;
	currNodeSet[globalNodeID]=c12;
	globalNodeID++;
	c12->size=c1->size+c2->size;
	
	if(c1->left!=NULL || c1->right!=NULL)
	{
		backup[c1->nodeName]=c1;
	}
	if(c2->left!=NULL || c2->right!=NULL)
	{
		backup[c2->nodeName]=c2;
	}
	return true;
}

int
HierarchicalCluster::generateModules(map<int,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib)
{
	int moduleCnt=modules.size();
	/*
	char outFName[1024];
	sprintf(outFName,"%s/clusters_pw.txt",outputDir);
	ofstream oFile(outFName);
	*/
	for(map<int,HierarchicalClusterNode*>::iterator cIter=currNodeSet.begin();cIter!=currNodeSet.end();cIter++)
	{
		HierarchicalClusterNode* node=cIter->second;
		map<string,int>* moduleMembers=new map<string,int>;
		modules[moduleCnt]=moduleMembers;
		populateMembers(moduleMembers,node);

		/*
		for(map<string,int>::iterator aIter=moduleMembers->begin();aIter!=moduleMembers->end();aIter++)
		{
			HierarchicalClusterNode* childNode=origAttrib[aIter->first];
			oFile << aIter->first <<"||5\tModule||5\t" << moduleCnt <<"|1|"<< moduleCnt << endl;
			for(int i=0;i<childNode->expr.size();i++)
			{
				oFile <<aIter->first <<"||5\tExp"<< i << "||5\t" << childNode->expr[i] << "|2"<< endl; 
			}
			for(map<int,double>::iterator mIter=childNode->attrib.begin();mIter!=childNode->attrib.end();mIter++)
			{
				oFile << aIter->first <<"||5\t" << (mIter->first) << "||5\t1|2" << endl;
			}
		}
		oFile <<"|Spacer||"<<moduleCnt<<"|-" << endl;
		*/
		cout <<"Module: " << moduleCnt << "\tSize="<< moduleMembers->size() << endl;
	
		moduleCnt=moduleCnt+1;
	}
	//oFile << endl;
	return 0;
}


int
HierarchicalCluster::calculatePercentVarianceExplained(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib)
{
	double s_total=0;
	double s_err=0;
	map<int,double> globalMean;
	for(map<int,map<string,int>*>::iterator mIter=modules.begin();mIter!=modules.end();mIter++)
	{
		map<int,double> localMean;
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				if(localMean.find(i)==localMean.end())
				{
					localMean[i]=n->expr[i];
				}	
				else	
				{
					localMean[i]=localMean[i]+n->expr[i];
				}
			}
		}
		for(map<int,double>::iterator dIter=localMean.begin();dIter!=localMean.end();dIter++)
		{
			if(globalMean.find(dIter->first)==globalMean.end())
			{
				globalMean[dIter->first]=dIter->second;
			}
			else
			{
				globalMean[dIter->first]=globalMean[dIter->first]+dIter->second;
			}
			dIter->second=dIter->second/((double)geneset->size());
		}	
		double s_err_m=0;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				double diff=n->expr[i]-localMean[i];
				s_err_m=s_err_m+(diff*diff);
			}
		}
		s_err=s_err+s_err_m;
		localMean.clear();
	}
	for(map<int,double>::iterator dIter=globalMean.begin();dIter!=globalMean.end();dIter++)
	{
		dIter->second=dIter->second/((double)origAttrib.size());
	}
	for(map<int,map<string,int>*>::iterator mIter=modules.begin();mIter!=modules.end();mIter++)
	{
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				double diff=n->expr[i]-globalMean[i];
				s_total=s_total+(diff*diff);
			}
		}
	}
	double pcv=1.0-(s_err/s_total);
	cout <<"Percent variance explained " << pcv << endl;
	globalMean.clear();
	return 0;
}

int
HierarchicalCluster::calculateSilhouetteIndex(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib)
{
	map<string,double> silhouette;
	int positive_s=0;
	double totals=0;
	for(map<int,map<string,int>*>::iterator mIter=modules.begin(); mIter!=modules.end();mIter++)
	{
		double module_s=0;
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			double a=0;
			for(map<string,int>::iterator hIter=geneset->begin();hIter!=geneset->end();hIter++)
			{
				if(gIter==hIter)
				{
					continue;
				}
				a=a+n->distToNeighbors_CC[hIter->first];
			}
			if(geneset->size()>1)
			{
				a=a/(geneset->size()-1);
			}
			double minb=100;
			for(map<int,map<string,int>*>::iterator nIter=modules.begin();nIter!=modules.end();nIter++)
			{
				double b=0;
				map<string,int>* geneset2=nIter->second;
				if(mIter==nIter)
				{
					continue;
				}
				for(map<string,int>::iterator hIter=geneset2->begin();hIter!=geneset2->end();hIter++)
				{
					b=b+n->distToNeighbors_CC[hIter->first];
				}
				b=b/geneset2->size();
				if(b<minb)
				{
					minb=b;
				}
			}
			double s=minb-a;
			if(s<0)
			{
				s=s/a;
			}
			else
			{
				s=s/minb;
				positive_s=positive_s+1;
			}
			silhouette[gIter->first]=s;
			module_s=module_s+s;
		}
		totals=totals+module_s;
		cout <<"Silhoutte for module " << mIter->first <<" " << module_s/geneset->size() << endl;
	}
	cout <<"Silhoutte Avg. " << totals/origAttrib.size() << " total positive " << positive_s << " of total " << origAttrib.size()<< endl;
	return 0;
}


int
HierarchicalCluster::populateMembers(map<string,int>* members,HierarchicalClusterNode* node)
{
	if(node->left==NULL && node->right==NULL)
	{
		(*members)[node->nodeName]=0;
	}
	else
	{
		if(node->left!=NULL)
		{
			populateMembers(members,node->left);
		}
		if(node->right!=NULL)
		{
			populateMembers(members,node->right);
		}
	}
	return 0;
}


HierarchicalClusterNode*
HierarchicalCluster::getRoot()
{
	if(root==NULL)
	{
		if(backup.size()==0)
		{
			cerr <<"No nodes you fool!!" << endl;
			exit(-1);
		}
		findRoot(backup.begin()->second);
	}
	return root;
}

int
HierarchicalCluster::findRoot(HierarchicalClusterNode* n)
{
	if(n->parent==NULL)
	{
		root=n;
	}
	else
	{
		findRoot(n->parent);
	}
	return 0;
}

int 
sortfunc(const void* first, const void* second)
{
	int ind1=*((int*)first);	
	int ind2=*((int*)second);
	double pval1=sorteddist[ind1];
	double pval2=sorteddist[ind2];
	int compstat=0;
	if(pval1>pval2)
	{
		compstat=-1;
	}
	else if(pval1<pval2)
	{
		compstat=1;
	}
	return compstat;
}
