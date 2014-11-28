#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <queue>
#include "Heap.H"

Heap::Heap()
{
	root=NULL;
	left=NULL;
	right=NULL;
	parent=NULL;
}

Heap::~Heap()
{
}

Heap*
Heap::insertToHeapNoHeapify(string& n1, string& n2, double d)
{
	if((strcmp(n1.c_str(),"Pdlim4")==0) && (strcmp(n2.c_str(),"Trex1")==0))
	{
		cout <<"Stop here! "<< endl;
		showHeap();
	}
	
	Heap* n=new Heap;
	n->p.node1.append(n1.c_str());
	n->p.node2.append(n2.c_str());
	n->p.dist=d;
	if(root==NULL)
	{
		root=n;
		return n;
	}
	int depth=0;
	//Heap* potParent=findNewPosition(root,n,depth);
	Heap* potParent=findNewPosition_iterative(n);
	if(potParent==NULL)
	{
		root->parent=n;
		n->left=root;
		root=n;
		return n;
	}
	if(potParent->left==NULL)
	{
		potParent->left=n;
		n->parent=potParent;
	}
	else if(potParent->right==NULL)
	{
		potParent->right=n;
		n->parent=potParent;
	}
	else
	{
		Heap* leftchild=potParent->left;
		Heap* rightchild=potParent->right;
		if(leftchild->p.dist>=n->p.dist)
		{
			n->left=leftchild;
			n->parent=potParent;
			potParent->left=n;
			leftchild->parent=n;
		}
		else if(rightchild->p.dist>=n->p.dist)
		{
			n->right=rightchild;
			n->parent=potParent;
			potParent->right=n;
			rightchild->parent=n;
		}
		else
		{
			cerr <<"Something strange happened! Did not find a node with empty children" << endl;
			exit(-1);
		}
	}
	if(root->parent!=NULL)
	{
		cout <<"Some shit happened at " << root->p.node1 << " " << root->p.node2 << " root's parent is garbage"  << endl;
	}
	return n;	
}

Heap*
Heap::findNewPosition(Heap* currPos, Heap* findMe,int& depth)
{
	if(currPos->p.dist>findMe->p.dist)
	{
		return currPos->parent;
	}
	else if(currPos->left!=NULL && currPos->right!=NULL)
	{	
		int ldepth=depth+1;
		int rdepth=depth+1;
		Heap* fromLeft=findNewPosition(currPos->left,findMe,ldepth);
		//return fromLeft;
		Heap* fromRight=findNewPosition(currPos->right,findMe,rdepth);
		if(ldepth<rdepth)
		{
			return fromLeft;
		}
		else
		{
			return fromRight;
		}
	}
	return currPos;
}



Heap*
Heap::findNewPosition_iterative(Heap* findMe)
{
	Heap* currPos=root;
	bool notFound=true;
	int fromleft=0;
	int fromright=0;
	queue<Heap*> nodeQueue;
	nodeQueue.push(root);
	while(notFound && !nodeQueue.empty())
	{
		currPos=nodeQueue.front();
		nodeQueue.pop();
		if(currPos->p.dist>findMe->p.dist)
		{
			currPos=currPos->parent;
			notFound=false;
			break;
		}
		else if((currPos->left!=NULL)&& (currPos->right!=NULL))
		{
			nodeQueue.push(currPos->left);
			nodeQueue.push(currPos->right);	
		}
		else
		{
			notFound=false;
		}
	}
	while(!nodeQueue.empty())
	{
		nodeQueue.pop();
	}
	return currPos;
}

int 
Heap::deleteFromHeap_getLeaf(Heap* removeMe)
{
	Heap* leaf=getLeaf(removeMe);
	//Otherwise replace the removeMe by the leaf node
	if((strcmp(removeMe->p.node1.c_str(),"ADH2-BAR1-SPO11-FUS1-SIN3-RME1-SWI1")==0) && (strcmp(removeMe->p.node2.c_str(),"ALPHA1-HXT7-MFALPHA2-SNF2_SWI1-ARG5-HXT6-MCM1-STA2-STA1-SAG1-STE6-PHO11-PHO5")==0))
	{
		cout << "Stop here " << endl;
	}
	if(root->parent!=NULL)
	{
		cout <<"Some shit happened at " << root->p.node1 << " " << root->p.node2 << " root's parent is garbage"  << endl;
	}
	//cout <<"Deleting " << removeMe->p.node1 <<" " << removeMe->p.node2 << endl;
	if(leaf==removeMe)
	{
		if(leaf->parent!=NULL)
		{
			if(leaf->parent->left==leaf)
			{
				leaf->parent->left=NULL;
			}
			else if(leaf->parent->right==leaf)
			{
				leaf->parent->right=NULL;
			}
			removeMe->parent=NULL;
		}
		if(removeMe==root)
		{
			root=NULL;
		}
		delete removeMe;
		return 0;
	}
	//Now connect the current children of removeMe to leaf
	if(leaf!=removeMe->left)
	{
		leaf->left=removeMe->left;
		if(removeMe->left!=NULL)
		{
			removeMe->left->parent=leaf;
		}
	}
	if(leaf!=removeMe->right)
	{
		leaf->right=removeMe->right;
		if(removeMe->right!=NULL)
		{
			removeMe->right->parent=leaf;
		}
	}
	//next disconnect the leaf from it's.
	if(leaf->parent->left==leaf)
	{
		leaf->parent->left=NULL;	
		leaf->parent=NULL;
	}
	else if(leaf->parent->right==leaf)
	{
		leaf->parent->right=NULL;
		leaf->parent=NULL;
	}
	//Finally if removeMe had a parent then leaf's parent should be updated
	if(removeMe->parent!=NULL)
	{
		leaf->parent=removeMe->parent;
		if(removeMe->parent->left==removeMe)
		{
			removeMe->parent->left=leaf;
		}
		else if(removeMe->parent->right==removeMe)
		{
			removeMe->parent->right=leaf;
		}
		else
		{
			cerr <<"Dangling pointer for " << removeMe->p.node1<<"-" << removeMe->p.node2 << endl;	
			exit(-1);
		}
	}
	else
	{
		root=leaf;
		root->parent=NULL;
	}
	removeMe->parent=NULL;
	delete removeMe;
	heapifyDown(leaf);
	if(root->parent!=NULL)
	{
		cout <<"Some shit happened at " << root->p.node1 << " " << root->p.node2 << " root's parent is garbage"  << endl;
	}
	return 0;
}

Heap*
Heap::getLeaf(Heap* currPos)
{
	Heap* leaf=NULL;
	if(currPos->left==NULL && currPos->right==NULL)
	{
		return currPos;
	}
	if(currPos->left!=NULL)
	{
		leaf=getLeaf(currPos->left);
		return leaf;
	}
	if(currPos->right!=NULL)
	{
		leaf=getLeaf(currPos->right);
		return leaf;
	}
	return leaf;
}

int
Heap::heapifyDown(Heap* currPos)
{
	if(currPos->left==NULL && currPos->right==NULL)
	{
		return 0;
	}
	Heap* minChild=NULL;
	if(currPos->left!=NULL && currPos->right!=NULL)
	{
		minChild=currPos->right;
		if(currPos->left->p.dist<currPos->right->p.dist)
		{
			minChild=currPos->left;
		}
	}
	else if(currPos->left!=NULL)
	{
		minChild=currPos->left;
	}
	else if(currPos->right!=NULL)
	{
		minChild=currPos->right;
	}
	if(minChild->p.dist<currPos->p.dist)
	{
		//Need to shift
		if(currPos->parent!=NULL)
		{
			if(currPos->parent->left==currPos)
			{
				currPos->parent->left=minChild;
			}
			else if(currPos->parent->right==currPos)
			{
				currPos->parent->right=minChild;
			}
		}	
		minChild->parent=currPos->parent;
		currPos->parent=minChild;
		if(minChild->parent==NULL)
		{
			root=minChild;
		}
		if(currPos->left==minChild)
		{
			Heap* oldleft=minChild->left;
			minChild->left=currPos;
			currPos->left=oldleft;
			if(oldleft!=NULL)
			{
				oldleft->parent=currPos;
			}
		}
		else if(currPos->right==minChild)
		{
			Heap* oldright=minChild->right;
			minChild->right=currPos;
			currPos->right=oldright;
			if(oldright!=NULL)
			{
				oldright->parent=currPos;
			}
		}
		heapifyDown(currPos);				
	}
	return 0;
}

bool
Heap::checkHeap()
{
	bool chk=checkHeap(root);
	return chk;
}

bool
Heap::checkHeap(Heap* n)
{
	if(n==NULL)
	{
		return true;
	}
	if(n->left==NULL & n->right==NULL)
	{
		return true;
	}	
	if(n->left!=NULL && n->right!=NULL)
	{
		if(n->left->p.dist<n->p.dist)
		{
			cout <<"Heap violated at left child of " << n->p.node1<<" " << n->p.node2 << endl;
			return false;
		}
		if(n->right->p.dist<n->p.dist)
		{
			cout <<"Heap violated at right child of " << n->p.node1<<" " << n->p.node2 << endl;
			return false;
		}
	}
	else if(n->left!=NULL)
	{
		if(n->left->p.dist<n->p.dist)
		{
			cout <<"Heap violated at left child of " << n->p.node1<<" " << n->p.node2 << endl;
			return false;
		}
	}
	else if(n->right!=NULL)
	{
		if(n->right->p.dist<n->p.dist)
		{
			cout <<"Heap violated at right child of " << n->p.node1<<" " << n->p.node2 << endl;
			return false;
		}
	}
	return true;
}

bool
Heap::checkPointers()
{
	bool chk=checkPointers(root);
	return chk;
}

bool
Heap::checkPointers(Heap* n)
{
	if(n==NULL)
	{
		return true;
	}
	if(n->left==NULL && n->right==NULL)
	{
		return true;
	}
	if(n->left!=NULL)
	{
		if(n->left->parent!=n)
		{
			cerr <<"Left child pointer mismatch at " << n->p.node1 << " " << n->p.node2 << endl;
			exit(-1);
		}
	}
	if(n->right!=NULL)
	{
		if(n->right->parent!=n)
		{
			cerr <<"Left child pointer mismatch at " << n->p.node1 << " " << n->p.node2 << endl;
			exit(-1);
		}
	}
	return true;
}

Heap::Pair*
Heap::getMin()
{
	return &root->p;
}

bool
Heap::empty()
{
	if(root==NULL)
	{
		return true;
	}
	return false;
}

int
Heap::showHeap()
{
	showHeap(root);
	return 0;
}

int
Heap::showHeap(Heap* h)
{
	if(h==NULL)
	{
		return 0;
	}
	cout << h->p.node1<<"-" << h->p.node2 << " " << h->p.dist << endl;
	if(h->left!=NULL)
	{
		showHeap(h->left);
	}	
	if(h->right!=NULL)
	{
		showHeap(h->right);
	}
	return 0;
}

Heap*
Heap::getRoot()
{
	return root;
}
