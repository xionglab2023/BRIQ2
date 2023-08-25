/*
 * RNAGraph.h
 *
 *  Created on: 2020��11��24��
 *      Author: pengx
 */

#ifndef MOTIF_RNAGRAPH_H_
#define MOTIF_RNAGRAPH_H_

#include "motif/RGNode.h"
#include "motif/RGEdge.h"
#include "model/StructureModel.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>

namespace NSPmotif {

using namespace std;
using namespace NSPmodel;

class SubRNAGraph{
public:

	int seqLen;
	bool* connectToDownstream;
	bool* inGraph;
	bool* choosen;
	int B;
	int F;
	int* fragIndex;


	SubRNAGraph(int seqLen, bool* connectToDownstream, bool* choosen){
		this->seqLen = seqLen;
		this->inGraph = new bool[seqLen];
		for(int i=0;i<seqLen;i++){
			inGraph[i] = false;
		}
		B = 0;
		F = 0;
		this->fragIndex = new int[seqLen];
		for(int i=0;i<seqLen;i++){
			fragIndex[i] = -1;
		}
		this->connectToDownstream = new bool[seqLen];
		this->choosen = new bool[seqLen];
		for(int i=0;i<seqLen;i++){
			this->connectToDownstream[i] = connectToDownstream[i];
			this->choosen[i] = choosen[i];
		}
	}

	SubRNAGraph(vector<int> idList, int seqLen, bool* connectToDownstream, bool* choosen){
		this->seqLen = seqLen;
		this->inGraph = new bool[seqLen];
		for(int i=0;i<seqLen;i++){
			inGraph[i] = false;
		}
		for(int i=0;i<idList.size();i++){
			inGraph[idList[i]] = true;
		}

		B = 0;
		F = 0;
		this->fragIndex = new int[seqLen];
		for(int i=0;i<seqLen;i++){
			fragIndex[i] = -1;
		}
		this->connectToDownstream = new bool[seqLen];
		this->choosen = new bool[seqLen];
		for(int i=0;i<seqLen;i++){
			this->connectToDownstream[i] = connectToDownstream[i];
			this->choosen[i] = choosen[i];
		}
	}

	SubRNAGraph* copy(){
		SubRNAGraph* sg = new SubRNAGraph(seqLen, connectToDownstream, choosen);
		for(int i=0;i<seqLen;i++){
			sg->inGraph[i] = this->inGraph[i];
		}
		sg->updateInfo();
		return sg;
	}

	void updateInfo(){

		this->B = 0;
		this->F = 0;

		for(int i=0;i<seqLen;i++){
			this->fragIndex[i] = -1;
		}

		int currentFID = -1;
		bool inFrag = false;
		for(int i=0;i<seqLen;i++){
			if(!inFrag && inGraph[i]){
				F++;
				B++;
				currentFID++;
				fragIndex[i] = currentFID;
				inFrag = true;
			}
			else if(inFrag && connectToDownstream[i-1] && inGraph[i]){
				B++;
				fragIndex[i] = currentFID;
			}
			else if(inFrag && !connectToDownstream[i-1] && inGraph[i]){
				F++;
				B++;
				currentFID ++;
				fragIndex[i] = currentFID;
			}
			else if(inFrag && !inGraph[i]){
				inFrag = false;
			}
		}
	}

	void printFragmentString(){
		for(int i=0;i<seqLen;i++){
			if(fragIndex[i] == -1)
				cout << '.';
			else {
				char c = 'A'+fragIndex[i];
				cout << c;
			}
		}
		cout << endl;
	}

	int getGraphLength(){
		if(B == 0)
			updateInfo();
		return B;
	}

	vector<SubRNAGraph*> generateAllVariant(){
		vector<SubRNAGraph*> varList;

		SubRNAGraph* sg1 = removeGap();
		if(getGraphLength() != sg1->getGraphLength()){
			varList.push_back(sg1);

			SubRNAGraph* sg11 = sg1->removeSingleBase();
			if(sg1->getGraphLength() != sg11->getGraphLength()){
				varList.push_back(sg11);
			}
			else
				delete sg11;

			SubRNAGraph* sg12 = sg1->removeSingleAndDoubleBase();
			if(sg11->getGraphLength() != sg12->getGraphLength()){
				varList.push_back(sg12);
			}
			else
				delete sg12;
		}
		else
			delete sg1;

		SubRNAGraph* sg2 = removeSingleAndDoubleGap();
		if(sg1->getGraphLength() != sg2->getGraphLength())
		{
			varList.push_back(sg2);
			SubRNAGraph* sg21 = sg2->removeSingleBase();
			if(sg2->getGraphLength() != sg21->getGraphLength())
				varList.push_back(sg21);
			else
				delete sg21;
			SubRNAGraph* sg22 = sg2->removeSingleAndDoubleBase();
			if(sg22->getGraphLength() != sg21->getGraphLength())
				varList.push_back(sg22);
			else
				delete sg22;
		}
		else
			delete sg2;

		SubRNAGraph* sg3 = removeSingleBase();
		if(getGraphLength() != sg3->getGraphLength()){
			varList.push_back(sg3);
			SubRNAGraph* sg31 = sg3->removeGap();
			if(sg31->getGraphLength() != sg3->getGraphLength())
				varList.push_back(sg31);
			else
				delete sg31;
			SubRNAGraph* sg32 = sg3->removeSingleAndDoubleGap();
			if(sg32->getGraphLength() != sg31->getGraphLength())
				varList.push_back(sg32);
			else
				delete sg32;
		}
		else
			delete sg3;

		SubRNAGraph* sg4 = removeSingleAndDoubleBase();
		if(sg3->getGraphLength() != sg4->getGraphLength()){
			varList.push_back(sg4);
			SubRNAGraph* sg41 = sg4->removeGap();
			if(sg41->getGraphLength() != sg4->getGraphLength())
				varList.push_back(sg41);
			else
				delete sg41;
			SubRNAGraph* sg42 = sg4->removeSingleAndDoubleGap();
			if(sg42->getGraphLength() != sg41->getGraphLength())
				varList.push_back(sg42);
			else
				delete sg42;
		}
		else
			delete sg4;
		return varList;
	}

	SubRNAGraph* removeGap(){
		SubRNAGraph* sg = new SubRNAGraph(seqLen, connectToDownstream, choosen);
		for(int i=0;i<seqLen;i++){
			sg->inGraph[i] = this->inGraph[i];
		}
		for(int i=0;i<seqLen-2;i++){
			if(sg->inGraph[i] && sg->inGraph[i+2] && !sg->inGraph[i+1] && connectToDownstream[i] && connectToDownstream[i+1] && !choosen[i+1]){
				sg->inGraph[i+1] = true;
			}
		}
		sg->updateInfo();
		return sg;
	}

	SubRNAGraph* removeSingleAndDoubleGap(){
		SubRNAGraph* sg = new SubRNAGraph(seqLen, connectToDownstream, choosen);
		for(int i=0;i<seqLen;i++){
			sg->inGraph[i] = this->inGraph[i];
		}
		for(int i=0;i<seqLen-2;i++){
			if(sg->inGraph[i] && sg->inGraph[i+2] && !sg->inGraph[i+1]){
				sg->inGraph[i+1] = true;
			}
		}
		for(int i=0;i<seqLen-3;i++){
			if(sg->inGraph[i] && sg->inGraph[i+3] && !sg->inGraph[i+1] && !sg->inGraph[i+2] && connectToDownstream[i] && connectToDownstream[i+1] && connectToDownstream[i+2] && !choosen[i+1] && !choosen[i+2]){
				sg->inGraph[i+1] = true;
				sg->inGraph[i+2] = true;
			}
		}
		sg->updateInfo();
		return sg;
	}

	SubRNAGraph* removeSingleBase(){

		if(B == 0)
			updateInfo();

		SubRNAGraph* sg = new SubRNAGraph(seqLen, connectToDownstream, choosen);
		for(int i=0;i<seqLen;i++){
			sg->inGraph[i] = this->inGraph[i];
		}

		for(int f=0;f<F;f++){
			int nf = 0;
			for(int i=0;i<seqLen;i++){
				if(fragIndex[i] == f) nf++;
			}
			if(nf == 1){
				for(int i=0;i<seqLen;i++){
					if(fragIndex[i] == f) sg->inGraph[i] = false;
				}
			}
		}
		sg->updateInfo();
		return sg;
	}

	SubRNAGraph* removeSingleAndDoubleBase(){
		if(B == 0)
			updateInfo();

		SubRNAGraph* sg = new SubRNAGraph(seqLen, connectToDownstream, choosen);
		for(int i=0;i<seqLen;i++){
			sg->inGraph[i] = this->inGraph[i];
		}

		for(int f=0;f<F;f++){
			int nf = 0;
			for(int i=0;i<seqLen;i++){
				if(fragIndex[i] == f) nf++;
			}
			if(nf <= 2){
				for(int i=0;i<seqLen;i++){
					if(fragIndex[i] == f) sg->inGraph[i] = false;
				}
			}
		}
		sg->updateInfo();
		return sg;
	}

	virtual ~SubRNAGraph();
};

class HelixGraph{
public:
	int beginIndex;
	int length; //helix length

	int simpleHelixBegin;
	int simpleHelixLen;

	double* interactionsOutsideHelix; //helix length
	vector<int> contactBaseList;

	HelixGraph(int beginIndex, int len){
		this->beginIndex = beginIndex;
		this->length = len;
		this->interactionsOutsideHelix = new double[len];
		this->simpleHelixBegin = 0;
		this->simpleHelixLen = 0;
	}

	bool isIndependantHelix(){
		double e1 = 0;
		for(int i=1;i<length-1;i++){
			e1 += interactionsOutsideHelix[i];
		}

		if(e1 > -0.5){
			this->simpleHelixBegin = beginIndex + 1;
			this->simpleHelixLen = length-2;
			return true;
		}
		else {
			int currentLen = 0;
			int currentBegin = 0;
			bool inSimpleHelix = false;
			for(int i=1;i<length-1;i++){
				if(!inSimpleHelix && interactionsOutsideHelix[i] > -0.3){
					currentLen = 1;
					currentBegin = beginIndex + i;
					inSimpleHelix = true;
				}
				else if(inSimpleHelix && interactionsOutsideHelix[i] > -0.3){
					currentLen++;
				}
				else if(inSimpleHelix){
					if(currentLen > simpleHelixLen){
						simpleHelixLen = currentLen;
						simpleHelixBegin = currentBegin;
						inSimpleHelix = false;
					}
					else{
						inSimpleHelix = false;
					}
				}
			}
			if(inSimpleHelix && interactionsOutsideHelix[length-2] > -0.3){
				if(currentLen > simpleHelixLen){
					simpleHelixLen = currentLen;
					simpleHelixBegin = currentBegin;
				}
			}
			return false;
		}
	}



	void printInfo(){
		isIndependantHelix();
		cout << "helix begin at: " <<  beginIndex << " length: " << length << endl;
		cout << "simple helix begin at: " << simpleHelixBegin << " length " << simpleHelixLen << endl;
		for(int i=0;i<length;i++){
			printf("%-4d %8.3f\n", beginIndex+i, interactionsOutsideHelix[i]);
		}

		for(int i=0;i<contactBaseList.size();i++){
			cout << contactBaseList[i] << endl;
		}
		if(isIndependantHelix())
			cout << "independant" << endl;
		else
			cout << "contact to other region" << endl;
	}

	virtual ~HelixGraph();
};

class RNAGraph {
public:
	int seqLen;
	string seq;
	string ss;
	RGNode** allNodes;
	bool* connectToDownstream;
	int* pairingIndex;
	bool* inHelix;

	bool* choosen;

	double* contactEnergyOfNonNeighbor;
	double* allEnergy;

	vector<RGEdge*> egList;
	map<int, int> baseIndexToEgIndex;

	vector<HelixGraph*> helixList;

	int* initialGraphIndex;
	int initialClusterNumber;
	int* subGraphIndex;


	vector<SubRNAGraph*> subGraphs;

	RNAGraph(const string& ssFile, const string& eneFile);

	void splitHelix(int i);
	void removeEdgeOnNode(int i);
	void addEdge();

	void removeEdge(int i, int j);
	void findHelix();

	void initialAssign(int nodeID, int graphID);
	void assign(int nodeID, int graphID, double eneCutoff);
	void updateSubGraph(double eneCutoff, int initialClusterID);

	double getEdgeEnergy(int i, int j);

	double getLengthFactor(int B){

		if(B < 4) return 0;
		double u = (B-40)/3.0;
		return 1.0/(1.0+exp(-u));

	}
	double getGraphScore(SubRNAGraph* sg);

	int getSep(int i, int j){
		if(i == j) return 0;
		if(j < i) return getSep(j, i);
		int sep = j-i;
		if(sep == 1 && !connectToDownstream[i]) return 5;
		if(sep == 2 && (!connectToDownstream[i] || !connectToDownstream[i+1])) return 5;
		return sep;
	}

	void printConnect();
	void printSubIndex(){
		for(int i=0;i<seqLen;i++){
			cout << subGraphIndex[i] << endl;
		}
	}

	void printNodeEnergy(){
		for(int i=0;i<seqLen;i++){
			printf("e: %-3d %8.3f\n", i, contactEnergyOfNonNeighbor[i]);
		}
	}

	void printSubGraph() {
		for(int i=0;i<subGraphs.size();i++){
			subGraphs[i]->updateInfo();

			printf("subGraph: %d B=%2d F=%2d Score=%8.3f\n", i, subGraphs[i]->B, subGraphs[i]->F, getGraphScore(subGraphs[i]));
			cout << ss << endl;
			subGraphs[i]->printFragmentString();
		}
	}

	void printHelix(){
		for(int i=0;i<this->helixList.size();i++){
			this->helixList[i]->printInfo();
		}
	}



	void processLoop();

	void processHelix();

	void generateInitialClusters();
	void printInitialCluster(int id);

	SubRNAGraph* findBestSubGraph(int initialClusterID);

	vector<SubRNAGraph*> findAllGraphs();

	void deleteSubGraph(SubRNAGraph* sg);

	virtual ~RNAGraph();
};

} /* namespace NSPmotif */

#endif /* MOTIF_RNAGRAPH_H_ */