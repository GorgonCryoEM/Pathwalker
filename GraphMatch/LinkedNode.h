/*
Storage structure for a node in the graph traversal.

Author: Sasakthi S. Abeysinghe
Date  : 01/23/2006
*/

#ifndef LINKEDNODE_H
#define LINKEDNODE_H

#include <vector>
#include "GlobalConstants.h"
#include "LinkedNodeStub.h"

using namespace std;


class LinkedNode : public LinkedNodeStub {
public:
	unsigned long long m1Bitmap;
	unsigned long long m2Bitmap;
	char missingNodesUsed;
	char depth;
	double cost;
	double costGStar;
public:
	LinkedNode(LinkedNode * olderNode);
	LinkedNode();
	LinkedNode(LinkedNode * olderNode, LinkedNodeStub * olderStub, int n2Node, int dummyHelixCount);
	LinkedNode(LinkedNode * olderNode, LinkedNodeStub * olderStub, int n1Node, int n2Node, int dummyHelixCount);
	~LinkedNode();
	void PrintNodeConcise(int rank, bool endOfLine = true, bool printCostBreakdown = false);
	unsigned long long GetN1Bitmap();
	unsigned long long GetN2Bitmap();
	bool operator==(LinkedNode &other);
	bool IsUserSpecifiedSolution();
	static void AddNodeToBitmap(unsigned long long & bitmap, int node);
	static void RemoveNodeFromBitmap(unsigned long long & bitmap, int node);
	static bool IsNodeInBitmap(unsigned long long bitmap, int node);
	static int RemoveSmallestNode(unsigned long long & bitmap);
};

LinkedNode::~LinkedNode() {
}

LinkedNode::LinkedNode() {
	cost = 0;
	n1Node = 0;
	n2Node = 0;
	parentNode = NULL;
	costGStar = 0;
	missingNodesUsed = 0;
	m1Bitmap = 0;
	m2Bitmap = 0;
	depth = 0;
}

LinkedNode::LinkedNode(LinkedNode * olderNode) {
	n1Node = olderNode->n1Node;
	n2Node = olderNode->n2Node;
	parentNode = olderNode->parentNode;
	m1Bitmap = olderNode->m1Bitmap;
	m2Bitmap = olderNode->m2Bitmap;
	cost = olderNode->cost;
	costGStar = olderNode->costGStar;
	missingNodesUsed = olderNode->missingNodesUsed;
	depth = olderNode->depth;
}


LinkedNode::LinkedNode(LinkedNode * olderNode, LinkedNodeStub * olderStub, int n2Node, int dummyHelixCount) {
	this->n1Node = olderNode->n1Node + dummyHelixCount + 1;
	this->n2Node = n2Node;
	this->depth = olderNode->depth + dummyHelixCount + 1;
	this->parentNode = olderStub;
	this->m1Bitmap = olderNode->m1Bitmap;
	this->m2Bitmap = olderNode->m2Bitmap;
	for(int i = olderNode->n1Node + 1; i <= this->n1Node; i++) {
		LinkedNode::RemoveNodeFromBitmap(this->m1Bitmap, i);
	}
	LinkedNode::RemoveNodeFromBitmap(this->m2Bitmap, this->n2Node);
	this->missingNodesUsed = olderNode->missingNodesUsed + dummyHelixCount;
}

LinkedNode::LinkedNode(LinkedNode * olderNode, LinkedNodeStub * olderStub, int n1Node, int n2Node, int dummyHelixCount) {
	this->n1Node = n1Node + dummyHelixCount + 1;
	this->n2Node = n2Node;
	this->depth = olderNode->depth + dummyHelixCount + 1;
	this->parentNode = olderStub;
	this->m1Bitmap = olderNode->m1Bitmap;
	this->m2Bitmap = olderNode->m2Bitmap;
	for(int i = olderNode->n1Node + 1; i <= this->n1Node; i++) {
		LinkedNode::RemoveNodeFromBitmap(this->m1Bitmap, i);
	}
	LinkedNode::RemoveNodeFromBitmap(this->m2Bitmap, this->n2Node);
	this->missingNodesUsed = olderNode->missingNodesUsed + dummyHelixCount;
}


void LinkedNode::PrintNodeConcise(int rank, bool endOfLine, bool printCostBreakdown) {
	bool used[MAX_NODES];
	int n1[MAX_NODES];
	int n2[MAX_NODES];
	int top = 0;
	for(int i = 0; i < MAX_NODES; i++) {
		used[i] = false;
	}

	LinkedNodeStub * currentNode = this;
	bool continueLoop = true;
	while(continueLoop) {
		if(currentNode->parentNode == NULL) {
			 break;
		}
		n1[top] = currentNode->n1Node;
		n2[top] = currentNode->n2Node;
		used[currentNode->n1Node] = true;
		top++;
		currentNode = currentNode->parentNode;		
	}

	for(int i = 1; i <= this->depth; i++) {
		if(!used[i]) {
			n1[top] = i;
			n2[top] = -1;
			top++;
		}
	}

	int minIndex;
	int temp;
	for(int i = 0; i < top - 1; i++) {
		minIndex = i;
		for(int j = i+1; j < top; j++) {
			if(n1[minIndex] > n1[j]) {
				minIndex = j;
			}
		}
		temp = n1[minIndex];
		n1[minIndex] = n1[i];
		n1[i] = temp;

		temp = n2[minIndex];
		n2[minIndex] = n2[i];
		n2[i] = temp;
	}

	if(IsUserSpecifiedSolution()) {
		printf("*");
	} else {
		printf(" ");
	}

	if(rank != -1) {
		printf("%d)", rank);
	}
	printf("\t");
	for(int i = 0; i < top; i++) {
		printf("%d ", n2[i]);
	}
	if(printCostBreakdown) {
		printf(" - %f = %f + %f", cost, costGStar, cost - costGStar);
	} else {
		printf(" - %f", cost);
	}
	if(endOfLine) {
		printf("\n");
	}
}

bool LinkedNode::IsUserSpecifiedSolution() {
	bool used[MAX_NODES];
	int n1[MAX_NODES];
	int n2[MAX_NODES];
	int top = 0;
	for(int i = 0; i < MAX_NODES; i++) {
		used[i] = false;
	}

	LinkedNodeStub * currentNode = this;
	bool continueLoop = true;
	while(continueLoop) {
		if(currentNode->parentNode == NULL) {
			 break;
		}
		n1[top] = currentNode->n1Node;
		n2[top] = currentNode->n2Node;
		used[currentNode->n1Node] = true;
		top++;
		currentNode = currentNode->parentNode;		
	}

	for(int i = 1; i <= this->depth; i++) {
		if(!used[i]) {
			n1[top] = i;
			n2[top] = -1;
			top++;
		}
	}

	int minIndex;
	int temp;
	for(int i = 0; i < top - 1; i++) {
		minIndex = i;
		for(int j = i+1; j < top; j++) {
			if(n1[minIndex] > n1[j]) {
				minIndex = j;
			}
		}
		temp = n1[minIndex];
		n1[minIndex] = n1[i];
		n1[i] = temp;

		temp = n2[minIndex];
		n2[minIndex] = n2[i];
		n2[i] = temp;
	}

	bool isSolution = true;
	for(int i = 0; i < top; i++) {
		isSolution = (isSolution && (n2[i] == SOLUTION[i]));
	}
	return isSolution;
}

unsigned long long LinkedNode::GetN1Bitmap() {
	unsigned long long n1Bitmap = 0;
	for(int i = 1; i <= this->n1Node; i++) {
		AddNodeToBitmap(n1Bitmap, i);
	}
	return n1Bitmap;
}

unsigned long long LinkedNode::GetN2Bitmap() {
	unsigned long long n2Bitmap = 0;
	LinkedNodeStub * currentNode = this;
	bool continueLoop = true;
	while(continueLoop) {
		if(currentNode->parentNode == NULL) {
			 break;
		}
		LinkedNode::AddNodeToBitmap(n2Bitmap, currentNode->n2Node);
		currentNode = currentNode->parentNode;
	}
	return n2Bitmap;
}

bool LinkedNode::operator==(LinkedNode &other) {
	return (n1Node == other.n1Node) && (n2Node == other.n2Node) && (m1Bitmap == other.m1Bitmap) && (m2Bitmap == other.m2Bitmap) && (parentNode == other.parentNode);
}

void LinkedNode::AddNodeToBitmap(unsigned long long & bitmap, int node) {
	bitmap = bitmap | ((unsigned long long)1 << node);
}

void LinkedNode::RemoveNodeFromBitmap(unsigned long long & bitmap, int node) {
	bitmap = bitmap - ((unsigned long long)1 << node);
}

bool LinkedNode::IsNodeInBitmap(unsigned long long bitmap, int node) {
	unsigned long long bitvalue = ((unsigned long long)1 << node);
	return ((bitmap & bitvalue) == bitvalue);
}


int LinkedNode::RemoveSmallestNode(unsigned long long & bitmap) {
	for(int i = 1; i <= MAX_NODES; i++) {
		if (IsNodeInBitmap(bitmap, i)) {
			RemoveNodeFromBitmap(bitmap, i);
			return i;
		}
	}
	assert(false);
	return 0;
}
#endif