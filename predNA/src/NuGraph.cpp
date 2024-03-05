/*
 * NuGraph.cpp
 *
 *  Created on: 2023��11��15��
 *      Author: nuc
 */

#include <predNA/NuGraph.h>

namespace NSPpredNA {

NuNode::NuNode(int id, int baseType, LocalFrame& cs1, BaseRotamer* baseRot, RiboseRotamer* riboRot, AtomLib* atLib){

	this->seqID = id;
	this->baseType = baseType;
	this->connectToNeighbor = false;

	this->baseConf = new BaseConformer(baseRot, cs1);
	this->baseConfTmp = new BaseConformer(baseRot, cs1);

	this->riboseConf = new RiboseConformer(riboRot, cs1);
	this->riboseConfTmp = new RiboseConformer(riboRot, cs1);

	this->phoConf = new PhosphateConformer();
	this->phoConfTmp = new PhosphateConformer();

	this->baseConfCG = NULL;
	this->baseConfCGTmp = NULL;
	this->riboseConfCG = NULL;
	this->riboseConfCGTmp = NULL;

	this->graph = NULL;
	this->ene = 0;
	this->eneTmp = 0;
	this->eneCG = 0;
	this->eneCGTmp = 0;

	this->bbcg = 0.0;
	this->bbcgTmp = 0.0;

	this->samplingFreq = 1.0;
}

NuNode::NuNode(int id, int baseType, LocalFrame& cs1, BaseRotamerCG* baseRotCG, RiboseRotamerCG* riboRotCG, AtomLib* atLib){


	this->seqID = id;
	this->baseType = baseType;
	this->connectToNeighbor = false;


	this->baseConf = NULL;
	this->baseConfTmp = NULL;

	this->riboseConf = NULL;
	this->riboseConfTmp = NULL;

	this->phoConf = new PhosphateConformer();
	this->phoConfTmp = new PhosphateConformer();

	this->baseConfCG = new BaseConformerCG(baseRotCG, cs1);
	this->baseConfCGTmp = new BaseConformerCG(baseRotCG, cs1);
	this->riboseConfCG = new RiboseConformerCG(riboRotCG, cs1);
	this->riboseConfCGTmp = new RiboseConformerCG(riboRotCG, cs1);

	this->graph = NULL;
	this->ene = 0;
	this->eneTmp = 0;
	this->eneCG = 0;
	this->eneCGTmp = 0;

	this->bbcg = 0.0;
	this->bbcgTmp = 0.0;

	this->samplingFreq = 1.0;
}

NuNode::NuNode(int id, int baseType, LocalFrame& cs1, BaseRotamer* baseRot, BaseRotamerCG* baseRotCG,  RiboseRotamer* riboRot, RiboseRotamerCG* riboRotCG, AtomLib* atLib){


	this->seqID = id;
	this->baseType = baseType;
	this->connectToNeighbor = false;


	this->baseConf = new BaseConformer(baseRot, cs1);
	this->baseConfTmp = new BaseConformer(baseRot, cs1);

	this->riboseConf = new RiboseConformer(riboRot, cs1);
	this->riboseConfTmp = new RiboseConformer(riboRot, cs1);

	this->phoConf = new PhosphateConformer();
	this->phoConfTmp = new PhosphateConformer();

	this->baseConfCG = new BaseConformerCG(baseRotCG, cs1);
	this->baseConfCGTmp = new BaseConformerCG(baseRotCG, cs1);

	this->riboseConfCG = new RiboseConformerCG(riboRotCG, cs1);
	this->riboseConfCGTmp = new RiboseConformerCG(riboRotCG, cs1);

	this->graph = NULL;
	this->ene = 0;
	this->eneTmp = 0;
	this->eneCG = 0;
	this->eneCGTmp = 0;

	this->bbcg = 0.0;
	this->bbcgTmp = 0.0;

	this->samplingFreq = 1.0;
}

NuNode::~NuNode(){

	if(this->baseConf != NULL)
		delete this->baseConf;
	if(this->baseConfTmp != NULL)
		delete this->baseConfTmp;
	if(this->riboseConf != NULL)
		delete this->riboseConf;
	if(this->riboseConfTmp != NULL)
		delete this->riboseConfTmp;
	if(this->phoConf != NULL)
		delete this->phoConf;
	if(this->phoConfTmp != NULL)
		delete this->phoConfTmp;


	if(this->baseConfCG != NULL)
		delete this->baseConfCG;
	if(this->baseConfCGTmp != NULL)
		delete this->baseConfCGTmp;
	if(this->riboseConfCG != NULL)
		delete this->riboseConfCG;
	if(this->riboseConfCGTmp != NULL)
		delete this->riboseConfCGTmp;
}

void NuNode::updateNodeInformation(NuTree* tree){

	this->graph = tree->graph;

	int i,j;
	this->neighborList.clear();
	this->connectionBreakPoints.clear();
	this->baseGroupA.clear();
	this->riboGroupA.clear();
	this->phoGroupA.clear();
	this->phoGroupC.clear();

	this->connectToNeighbor = tree->graph->connectToDownstream[seqID];

 	for(i=0;i<graph->seqLen;i++){
		if(i == seqID) continue;
		if(tree->adjMtx[seqID*tree->graph->seqLen+i])
			this->neighborList.push_back(i);
	}

	int phoGroup[graph->seqLen];
	for(i=0;i<graph->seqLen;i++){
		phoGroup[i] = 0;
	}

	if(this->seqID>0 &&  tree->graph->connectToDownstream[this->seqID-1]) {
		this->connectionBreakPoints.push_back(this->seqID-1);
		phoGroup[seqID-1] = 2;
	}

	if(tree->graph->connectToDownstream[this->seqID]) {
		this->connectionBreakPoints.push_back(this->seqID);
		phoGroup[seqID] = 2;
	}


	for(i=0;i<graph->seqLen;i++){
		if(graph->masked[i]) continue;

		if(i != this->seqID) {
			baseGroupA.push_back(graph->allNodes[i]);
			riboGroupA.push_back(graph->allNodes[i]);
		}

		if(phoGroup[i] == 0 && graph->connectToDownstream[i])
			phoGroupA.push_back(graph->allNodes[i]);
		else if(graph->connectToDownstream[i])
			phoGroupC.push_back(graph->allNodes[i]);
	}

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->et->pb->buildPhosphate(graph->allNodes[j]->riboseConf, graph->allNodes[j+1]->riboseConf, graph->allNodes[j]->phoConf);
		graph->allNodes[j]->phoConfTmp->copyValueFrom(graph->allNodes[j]->phoConf);
	}

	this->ene = riboseConf->rot->energy;
	if(connectToNeighbor)
		this->ene += phoConf->ene;

	this->eneTmp = this->ene;
}

void NuNode::updateNodeInformationCG(NuTree* tree){

	this->graph = tree->graph;
	int i,j;
	this->neighborList.clear();
	this->connectionBreakPoints.clear();
	this->baseGroupA.clear();
	this->riboGroupA.clear();
	this->phoGroupA.clear();
	this->phoGroupC.clear();

	this->connectToNeighbor = tree->graph->connectToDownstream[seqID];

 	for(i=0;i<graph->seqLen;i++){
		if(i == seqID) continue;
		if(tree->adjMtx[seqID*tree->graph->seqLen+i])
			this->neighborList.push_back(i);
	}

	int phoGroup[graph->seqLen];
	for(i=0;i<graph->seqLen;i++){
		phoGroup[i] = 0;
	}

	if(this->seqID>0 &&  tree->graph->connectToDownstream[this->seqID-1]) {
		this->connectionBreakPoints.push_back(this->seqID-1);
		phoGroup[seqID-1] = 2;
	}

	if(tree->graph->connectToDownstream[this->seqID]) {
		this->connectionBreakPoints.push_back(this->seqID);
		phoGroup[seqID] = 2;
	}


	for(i=0;i<graph->seqLen;i++){
		if(graph->masked[i]) continue;
		if(i != this->seqID){
			baseGroupA.push_back(graph->allNodes[i]);
			riboGroupA.push_back(graph->allNodes[i]);
		}

		if(phoGroup[i] == 0 && graph->connectToDownstream[i])
			phoGroupA.push_back(graph->allNodes[i]);
		else if(graph->connectToDownstream[i])
			phoGroupC.push_back(graph->allNodes[i]);
	}

	this->eneCG = riboseConfCG->rot->energy;

	if(connectToNeighbor) {
		this->eneCG += nuConnectionEnergyCG(this->riboseConfCG, graph->allNodes[seqID+1]->riboseConfCG, graph->et);
	}
	this->eneCGTmp = this->eneCG;
}

void NuNode::printNodeInfo(){
	cout << "nodeID: " << seqID << endl;

	cout << "baseGroupA: " << endl;
	for(int i=0;i<baseGroupA.size();i++){
		cout << baseGroupA[i]->seqID << " ";
	}
	cout << endl;

	cout << "riboGroupA: " << endl;
	for(int i=0;i<riboGroupA.size();i++){
		cout << riboGroupA[i]->seqID << " ";
	}
	cout << endl;

	cout << "connection breaks: ";
	for(int i=0;i<connectionBreakPoints.size();i++){
		cout << " " << connectionBreakPoints[i];
	}
	cout << endl;

	cout << "energy: " << this->ene <<  " energyCG: " << this->eneCG << endl;
	cout << endl;
}

void NuNode::updateRiboseRotamer(RiboseRotamer* rot){
	int len = graph->seqLen;
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;

	this->riboseConfTmp->updateRotamer(rot);
	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->et->pb->buildPhosphate(graph->allNodes[j]->riboseConfTmp, graph->allNodes[j+1]->riboseConfTmp, graph->allNodes[j]->phoConfTmp);
		graph->allNodes[j]->eneTmp = graph->allNodes[j]->riboseConfTmp->rot->energy + graph->allNodes[j]->phoConfTmp->ene;

	}
	this->eneTmp = rot->energy + phoConfTmp->ene;

	/*
	 * base-ribose energy
	 */

	for(i=0;i<baseGroupA.size();i++){
		sep = graph->sepTable[baseGroupA[i]->seqID*len+seqID];
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + baseGroupA[i]->seqID];

		egA->pairEneTmp[1] = nuBaseRiboseEnergy(baseGroupA[i]->baseConfTmp, this->riboseConfTmp, sep, graph->et);
		egB->pairEneTmp[3] = egA->pairEneTmp[1];
	}

	/*
	 * base-pho
	 */
	for(i=0;i<baseGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			sep = graph->sepTable[baseGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			if(sep == 0) continue;
			egA = graph->allEdges[baseGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + baseGroupA[i]->seqID];

			egA->pairEneTmp[2] = nuBasePhoEnergy(baseGroupA[i]->baseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et);
			egB->pairEneTmp[6] = egA->pairEneTmp[2];
		}
	}

	for(j=0;j<phoGroupC.size();j++){
		sep = graph->sepTable[seqID*len + phoGroupC[j]->seqID];
		if(sep == 0) continue;
		egA = graph->allEdges[seqID*len + phoGroupC[j]->seqID];
		egB = graph->allEdges[phoGroupC[j]->seqID*len + seqID];

		egA->pairEneTmp[2] = nuBasePhoEnergy(baseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et);
		egB->pairEneTmp[6] = egA->pairEneTmp[2];
	}


	/*
	 * ribose-ribose
	 */

	for(i=0;i<riboGroupA.size();i++){

		sep = graph->sepTable[riboGroupA[i]->seqID*len+seqID];

		egA = graph->allEdges[riboGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + riboGroupA[i]->seqID];

		egA->pairEneTmp[4] = nuRiboseRiboseEnergy(riboGroupA[i]->riboseConfTmp, this->riboseConfTmp, sep, graph->et);
		egB->pairEneTmp[4] = egA->pairEneTmp[4];
	}

	/*
	 * riboseA-phoC
	 */
	for(i=0;i<riboGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			sep = graph->sepTable[riboGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			if(sep == 0) continue;
			egA = graph->allEdges[riboGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + riboGroupA[i]->seqID];

			egA->pairEneTmp[5] = nuRibosePhoEnergy(riboGroupA[i]->riboseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et);
			egB->pairEneTmp[7] = egA->pairEneTmp[5];
		}
	}

	/*
	 * riboseC-phoA
	 */
	for(i=0;i<phoGroupA.size();i++){
		sep = graph->sepTable[seqID*len+phoGroupA[i]->seqID];
		egA = graph->allEdges[seqID*len + phoGroupA[i]->seqID];
		egB = graph->allEdges[phoGroupA[i]->seqID*len + seqID];
		egA->pairEneTmp[5] = nuRibosePhoEnergy(riboseConfTmp, phoGroupA[i]->phoConfTmp, sep, graph->et);
		egB->pairEneTmp[7] = egA->pairEneTmp[5];
	}

	/*
	 * pho-pho
	 */

	for(i=0;i<phoGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			sep = graph->sepTable[phoGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			egA = graph->allEdges[phoGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len+phoGroupA[i]->seqID];
			egA->pairEneTmp[8] = nuPhoPhoEnergy(phoGroupA[i]->phoConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et);
			egB->pairEneTmp[8] = egA->pairEneTmp[8];
		}
	}

	for(i=0;i<phoGroupC.size();i++){
		for(j=i+1;j<phoGroupC.size();j++){
			sep = graph->sepTable[phoGroupC[i]->seqID*len+phoGroupC[j]->seqID];
			egA = graph->allEdges[phoGroupC[i]->seqID*len+phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len+phoGroupC[i]->seqID];
			egA->pairEneTmp[8] = nuPhoPhoEnergy(phoGroupC[i]->phoConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et);
			egB->pairEneTmp[8] = egA->pairEneTmp[8];
		}
	}
}

void NuNode::acceptRotMutation(){
	int i,j, sep;
	int len = graph->seqLen;
	NuEdge* egA;
	NuEdge* egB;
	this->riboseConf->copyValueFrom(riboseConfTmp);
	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->phoConf->copyValueFrom(graph->allNodes[j]->phoConfTmp);
		graph->allNodes[j]->ene = graph->allNodes[j]->eneTmp;
	}

	this->ene = this->eneTmp;

	/*
	 * base-ribose energy
	 */

	for(i=0;i<baseGroupA.size();i++){
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + baseGroupA[i]->seqID];

		egA->pairEne[1] = egA->pairEneTmp[1];
		egB->pairEne[3] = egB->pairEneTmp[3];
	}


	/*
	 * base-pho
	 */
	for(i=0;i<baseGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[baseGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + baseGroupA[i]->seqID];

			egA->pairEne[2] = egA->pairEneTmp[2];
			egB->pairEne[6] = egB->pairEneTmp[6];
		}
	}

	
	for(j=0;j<phoGroupC.size();j++){
		sep = graph->sepTable[seqID*len + phoGroupC[j]->seqID];
		if(sep == 0) continue;
		egA = graph->allEdges[seqID*len + phoGroupC[j]->seqID];
		egB = graph->allEdges[phoGroupC[j]->seqID*len + seqID];

		egA->pairEne[2] = egA->pairEneTmp[2];
		egB->pairEne[6] = egB->pairEneTmp[6];
	}

	/*
	 * ribose-ribose
	 */

	for(i=0;i<riboGroupA.size();i++){
		
		egA = graph->allEdges[riboGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + riboGroupA[i]->seqID];

		egA->pairEne[4] = egA->pairEneTmp[4];
		egB->pairEne[4] = egB->pairEneTmp[4];
	}

	/*
	 * riboseA-phoC
	 */
	for(i=0;i<riboGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			sep = graph->sepTable[riboGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			if(sep == 0) continue;
			egA = graph->allEdges[riboGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + riboGroupA[i]->seqID];

			egA->pairEne[5] = egA->pairEneTmp[5];
			egB->pairEne[7] = egB->pairEneTmp[7];
		}
	}

	/*
	 * riboseC-pho
	 */
	for(i=0;i<phoGroupA.size();i++){
		egA = graph->allEdges[seqID*len + phoGroupA[i]->seqID];
		egB = graph->allEdges[phoGroupA[i]->seqID*len + seqID];

		egA->pairEne[5] = egA->pairEneTmp[5];
		egB->pairEne[7] = egB->pairEneTmp[7];
	}


	/*
	 * pho-pho
	 */

	for(i=0;i<phoGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[phoGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len+phoGroupA[i]->seqID];
			egA->pairEne[8] = egA->pairEneTmp[8];
			egB->pairEne[8] = egB->pairEneTmp[8];
		}
	}

	for(i=0;i<phoGroupC.size();i++){
		for(j=i+1;j<phoGroupC.size();j++){
			egA = graph->allEdges[phoGroupC[i]->seqID*len+phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len+phoGroupC[i]->seqID];
			egA->pairEne[8] = egA->pairEneTmp[8];
			egB->pairEne[8] = egB->pairEneTmp[8];
		}
	}
}

void NuNode::clearRotMutation() {
	int i,j, sep;
	int len = graph->seqLen;
	NuEdge* egA;
	NuEdge* egB;
	this->riboseConfTmp->copyValueFrom(riboseConf);
	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->phoConfTmp->copyValueFrom(graph->allNodes[j]->phoConf);
	}


	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->eneTmp = graph->allNodes[j]->ene;
	}
	this->eneTmp = this->ene;

	/*
	 * base-ribose energy
	 */

	for(i=0;i<baseGroupA.size();i++){
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + baseGroupA[i]->seqID];

		egA->pairEneTmp[1] = egA->pairEne[1];
		egB->pairEneTmp[3] = egB->pairEne[3];
	}


	/*
	 * base-pho
	 */
	for(i=0;i<baseGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[baseGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + baseGroupA[i]->seqID];

			egA->pairEneTmp[2] = egA->pairEne[2];
			egB->pairEneTmp[6] = egB->pairEne[6];
		}
	}

	for(j=0;j<phoGroupC.size();j++){
		sep = graph->sepTable[seqID*len + phoGroupC[j]->seqID];
		if(sep == 0) continue;
		egA = graph->allEdges[seqID*len + phoGroupC[j]->seqID];
		egB = graph->allEdges[phoGroupC[j]->seqID*len + seqID];
		
		egA->pairEneTmp[2] = egA->pairEne[2];
		egB->pairEneTmp[6] = egB->pairEne[6];
	}

	/*
	 * ribose-ribose
	 */

	for(i=0;i<riboGroupA.size();i++){
		
		egA = graph->allEdges[riboGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + riboGroupA[i]->seqID];

		egA->pairEneTmp[4] = egA->pairEne[4];
		egB->pairEneTmp[4] = egB->pairEne[4];
	}

	/*
	 * riboseA-phoC
	 */
	for(i=0;i<riboGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			sep = graph->sepTable[riboGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			if(sep == 0) continue;
			egA = graph->allEdges[riboGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + riboGroupA[i]->seqID];

			egA->pairEneTmp[5] = egA->pairEne[5];
			egB->pairEneTmp[7] = egB->pairEne[7];
		}
	}

	/*
	 * riboseC-pho
	 */
	for(i=0;i<phoGroupA.size();i++){
		egA = graph->allEdges[seqID*len + phoGroupA[i]->seqID];
		egB = graph->allEdges[phoGroupA[i]->seqID*len + seqID];

		egA->pairEneTmp[5] = egA->pairEne[5];
		egB->pairEneTmp[7] = egB->pairEne[7];
	}

	/*
	 * pho-pho
	 */

	for(i=0;i<phoGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[phoGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len+phoGroupA[i]->seqID];

			egA->pairEneTmp[8] = egA->pairEne[8];
			egB->pairEneTmp[8] = egB->pairEne[8];
		}
	}

	for(i=0;i<phoGroupC.size();i++){
		for(j=i+1;j<phoGroupC.size();j++){
			egA = graph->allEdges[phoGroupC[i]->seqID*len+phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len+phoGroupC[i]->seqID];

			egA->pairEneTmp[8] = egA->pairEne[8];
			egB->pairEneTmp[8] = egB->pairEne[8];
		}
	}
}

void NuNode::updateCoordinate(LocalFrame& cs){
	this->baseConfTmp->updateCoords(cs);
	this->riboseConfTmp->updateLocalFrame(cs);
	if(connectToNeighbor)
		this->phoConfTmp->updateLocalFrame(riboseConfTmp->cs2);
}

void NuNode::acceptCoordMove(){
	this->baseConf->copyValueFrom(baseConfTmp);
	this->riboseConf->copyValueFrom(riboseConfTmp);
	if(connectToNeighbor)
		this->phoConf->copyValueFrom(phoConfTmp);
}

void NuNode::clearCoordMove() {
	this->baseConfTmp->copyValueFrom(baseConf);
	this->riboseConfTmp->copyValueFrom(riboseConf);
	if(connectToNeighbor)
		this->phoConfTmp->copyValueFrom(phoConf);
}

double NuNode::rotMutEnergy(){
	int i,j, sep;
	int len = graph->seqLen;
	NuEdge* egA;

	double mutE = 0.0;
	mutE = riboseConfTmp->rot->energy - riboseConf->rot->energy;

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		mutE += graph->allNodes[j]->phoConfTmp->ene - graph->allNodes[j]->phoConf->ene;
	}

	/*
	 * base-ribose energy
	 */

	for(i=0;i<baseGroupA.size();i++){
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		mutE += egA->pairEneTmp[1] - egA->pairEne[1];
	}

	/*
	 * base-pho
	 */
	for(i=0;i<baseGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[baseGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			mutE += egA->pairEneTmp[2] - egA->pairEne[2];
		}
	}

	for(j=0;j<phoGroupC.size();j++){
		sep = graph->sepTable[seqID*len + phoGroupC[j]->seqID];
		if(sep == 0) continue;
		egA = graph->allEdges[seqID*len + phoGroupC[j]->seqID];
		mutE += egA->pairEneTmp[2] - egA->pairEne[2];
	}

	/*
	 * ribose-ribose
	 */

	for(i=0;i<riboGroupA.size();i++){

		egA = graph->allEdges[riboGroupA[i]->seqID*len + seqID];
		mutE += egA->pairEneTmp[4] - egA->pairEne[4];
	}

	/*
	 * riboseA-phoC
	 */
	for(i=0;i<riboGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[riboGroupA[i]->seqID*len + phoGroupC[j]->seqID];
			mutE += egA->pairEneTmp[5] - egA->pairEne[5];
		}
	}

	/*
	 * riboseC-pho
	 */
	for(i=0;i<phoGroupA.size();i++){
		egA = graph->allEdges[seqID*len + phoGroupA[i]->seqID];
		mutE += egA->pairEneTmp[5] - egA->pairEne[5];
	}

	/*
	 * pho-pho
	 */

	for(i=0;i<phoGroupA.size();i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[phoGroupA[i]->seqID*len+phoGroupC[j]->seqID];
			mutE += egA->pairEneTmp[8] - egA->pairEne[8];
		}
	}

	for(i=0;i<phoGroupC.size();i++){
		for(j=i+1;j<phoGroupC.size();j++){
			egA = graph->allEdges[phoGroupC[i]->seqID*len+phoGroupC[j]->seqID];
			mutE += egA->pairEneTmp[8] - egA->pairEne[8];
		}
	}

	return mutE;
}

bool NuNode::checkEnergy(){
	double e1 = this->riboseConf->rot->energy + this->phoConf->ene;
	double e2 = this->riboseConfTmp->rot->energy + this->phoConfTmp->ene;
	bool tag = true;

	if(abs(e1 - ene) > 0.001){
		printf("Node: %3d ERROR ene %7.3f %7.3f\n",seqID, e1, ene);
		tag = false;
	}

	if(abs(e2 - eneTmp) > 0.001){
		printf("NodeTmp: %3d ERROR ene %7.3f %7.3f\n",seqID, e2, eneTmp);
		tag = false;
	}

	return tag;
}

void NuNode::updateRiboseRotamerCG(RiboseRotamerCG* rot){
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;

	int len = graph->seqLen;
	this->riboseConfCGTmp->updateRotamer(rot);

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->bbcgTmp = nuConnectionEnergyCG(graph->allNodes[j]->riboseConfCGTmp, graph->allNodes[j+1]->riboseConfCGTmp, graph->et);
		graph->allNodes[j]->eneCGTmp = graph->allNodes[j]->riboseConfCGTmp->rot->energy + graph->allNodes[j]->bbcgTmp;
	}

	this->eneCGTmp = rot->energy + this->bbcgTmp;

	/*
	 * base-ribose energy
	 * ribose-ribose energy
	 */

	for(i=0;i<baseGroupA.size();i++){

		sep = graph->sepTable[baseGroupA[i]->seqID*len+seqID];
		sepR = graph->sepTable[seqID*len+baseGroupA[i]->seqID];
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + baseGroupA[i]->seqID];

		egA->pairEneCGTmp[1] = nuBaseRiboseEnergyCG(baseGroupA[i]->baseConfCGTmp, this->riboseConfCGTmp, sep, graph->et);
		egB->pairEneCGTmp[2] = egA->pairEneCGTmp[1];
		egA->pairEneCGTmp[3] = nuRiboseRiboseEnergyCG(baseGroupA[i]->riboseConfCGTmp, this->riboseConfCGTmp, sep, graph->et);
		egB->pairEneCGTmp[3] = egA->pairEneCGTmp[3];
	}

}

void NuNode::acceptRotMutationCG(){

	this->riboseConfCG->copyValueFrom(this->riboseConfCGTmp);
	this->eneCG = this->eneCGTmp;

	int len = graph->seqLen;
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->bbcg = graph->allNodes[j]->bbcgTmp;
		graph->allNodes[j]->eneCG = graph->allNodes[j]->eneCGTmp;
	}


	for(i=0;i<baseGroupA.size();i++){
		sep = graph->sepTable[baseGroupA[i]->seqID*len+seqID];
		sepR = graph->sepTable[seqID*len+baseGroupA[i]->seqID];
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + baseGroupA[i]->seqID];

		egA->pairEneCG[1] = egA->pairEneCGTmp[1];
		egB->pairEneCG[2] = egB->pairEneCGTmp[2];
		egA->pairEneCG[3] = egA->pairEneCGTmp[3];
		egB->pairEneCG[3] = egB->pairEneCGTmp[3];
	}

}

void NuNode::clearRotMutationCG(){
	this->riboseConfCGTmp->copyValueFrom(this->riboseConfCG);
	this->eneCGTmp = this->eneCG;

	int len = graph->seqLen;
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->bbcgTmp = graph->allNodes[j]->bbcg;
		graph->allNodes[j]->eneCGTmp = graph->allNodes[j]->eneCG;
	}

	for(i=0;i<baseGroupA.size();i++){
		sep = graph->sepTable[baseGroupA[i]->seqID*len+seqID];
		sepR = graph->sepTable[seqID*len+baseGroupA[i]->seqID];
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + baseGroupA[i]->seqID];

		egA->pairEneCGTmp[1] = egA->pairEneCG[1]; //BR
		egB->pairEneCGTmp[2] = egB->pairEneCG[2]; //RB
		egA->pairEneCGTmp[3] = egA->pairEneCG[3]; //RR
		egB->pairEneCGTmp[3] = egB->pairEneCG[3]; //RR
	}
}

void NuNode::updateCoordinateCG(LocalFrame& cs){
	this->baseConfCGTmp->updateCoords(cs);
	this->riboseConfCGTmp->updateLocalFrame(cs);
}

void NuNode::acceptCoordMoveCG(){
	this->baseConfCG->copyValueFrom(baseConfCGTmp);
	this->riboseConfCG->copyValueFrom(riboseConfCGTmp);
}

void NuNode::clearCoordMoveCG(){
	this->baseConfCGTmp->copyValueFrom(baseConfCG);
	this->riboseConfCGTmp->copyValueFrom(riboseConfCG);
}

double NuNode::rotMutEnergyCG(){
	
	int i,j,sep;
	NuEdge* egA;
	int len = graph->seqLen;

	double mutE = this->riboseConfCGTmp->rot->energy - this->riboseConfCG->rot->energy;
	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		mutE += graph->allNodes[j]->bbcgTmp - graph->allNodes[j]->bbcg;
	}

	for(i=0;i<baseGroupA.size();i++){
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		mutE += egA->pairEneCGTmp[1] - egA->pairEneCG[1];
		mutE += egA->pairEneCGTmp[3] - egA->pairEneCG[3];
	}
	return mutE;
}

bool NuNode::checkEnergyCG(){
	double e3 = this->riboseConfCG->rot->energy;
	bool tag = true;

	if(connectToNeighbor)
		e3 += nuConnectionEnergyCG(this->riboseConfCG, graph->allNodes[seqID+1]->riboseConfCG, graph->et);

	double e4 = this->riboseConfCGTmp->rot->energy;
	if(connectToNeighbor)
		e4 += nuConnectionEnergyCG(this->riboseConfCGTmp, graph->allNodes[seqID+1]->riboseConfCGTmp, graph->et);

	if(abs(e3 - eneCG) > 0.001){
		printf("NodeCG %3d ERROR ene %7.3f  rotE %7.3f eneCG %7.3f\n",this->seqID, e3, this->riboseConfCG->rot->energy,  eneCG);
		tag = false;
	}

	if(abs(e4 - eneCGTmp) > 0.001){
		printf("NodeCGTmp %3d ERROR ene %7.3f  rotE %7.3f eneCG %7.3f\n",this->seqID, e4, this->riboseConfCGTmp->rot->energy, eneCGTmp);
		tag = false;
	}
	return tag;
}

vector<Atom*> NuNode::toAtomList(AtomLib* atLib){
	vector<Atom*> list;
	vector<string>* names = atLib->getRnaSidechainAtoms(this->baseType);
	vector<XYZ> tList;
	for(int i=0;i<baseConf->rot->atomNum;i++){
		tList.push_back(baseConf->coords[i]);
	}
    vector<Atom*> atomList;

    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names->at(i), tList[i]));
    }

    atomList.push_back(new Atom("C1'", riboseConf->coords[0]));
    atomList.push_back(new Atom("C2'", riboseConf->coords[1]));
    atomList.push_back(new Atom("C3'", riboseConf->coords[2]));
    atomList.push_back(new Atom("C4'", riboseConf->coords[3]));
    atomList.push_back(new Atom("O4'", riboseConf->coords[4]));
    atomList.push_back(new Atom("O3'", riboseConf->coords[5]));
    atomList.push_back(new Atom("C5'", riboseConf->coords[6]));
    if(this->baseType < 4) {
    	atomList.push_back(new Atom("O2'", riboseConf->coords[7]));
    }
    return atomList;
}

vector<Atom*> NuNode::toAtomListCG(AtomLib* atLib){
    vector<Atom*> atomList;
	atomList.push_back(new Atom("C1", baseConfCG->coords[0]));
	atomList.push_back(new Atom("C2", baseConfCG->coords[1]));
	atomList.push_back(new Atom("C3", baseConfCG->coords[2]));
	atomList.push_back(new Atom("C1'", riboseConfCG->coords[0]));
	atomList.push_back(new Atom("C1'", riboseConfCG->coords[1]));
	atomList.push_back(new Atom("C1'", riboseConfCG->coords[2]));
	return atomList;
}

vector<Atom*> NuNode::toAtomListWithPho(AtomLib* atLib){
	vector<Atom*> list;
	vector<string>* names = atLib->getRnaSidechainAtoms(this->baseType);
	vector<XYZ> tList;
	for(int i=0;i<baseConf->rot->atomNum;i++){
		tList.push_back(baseConf->coords[i]);
	}
    vector<Atom*> atomList;


    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names->at(i), tList[i]));
    }

    atomList.push_back(new Atom("C1'", riboseConf->coords[0]));
    atomList.push_back(new Atom("C2'", riboseConf->coords[1]));
    atomList.push_back(new Atom("C3'", riboseConf->coords[2]));
    atomList.push_back(new Atom("C4'", riboseConf->coords[3]));
    atomList.push_back(new Atom("O4'", riboseConf->coords[4]));
    atomList.push_back(new Atom("O3'", riboseConf->coords[5]));
    atomList.push_back(new Atom("C5'", riboseConf->coords[6]));
    if(this->baseType < 4) {
    	atomList.push_back(new Atom("O2'", riboseConf->coords[7]));
    }

    if(connectToNeighbor) {
        atomList.push_back(new Atom("P", phoConf->coords[0]));
        atomList.push_back(new Atom("O5'", phoConf->coords[1]));
        atomList.push_back(new Atom("OP1", phoConf->coords[2]));
        atomList.push_back(new Atom("OP2", phoConf->coords[3]));
    }

    return atomList;
}

NuEdge::NuEdge(NuNode* nodeA, NuNode* nodeB, NuGraph* graph){
	this->graph = graph;
	this->indexA = nodeA->seqID;
	this->indexB = nodeB->seqID;
	this->nodeA = nodeA;
	this->nodeB = nodeB;
	this->cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->cmTmp = this->cm;
	this->sep = graph->sepTable[indexA*graph->seqLen+indexB];

	int typeA = nodeA->baseType%4;
	int typeB = nodeB->baseType%4;

	this->ei = new EdgeInformation(sep, typeA, typeB, graph->pairLib);
	this->moveSet = new MixedNuPairCluster(sep, typeA*4+typeB, graph->moveLib);
	this->moveSet->updateEdgeInformation(this->ei);
	this->weight = this->ei->weight;

	this->weightRand = weight;

	this->pairLib = graph->pairLib;

	for(int i=0;i<9;i++){
		this->pairEne[i] = 0.0;
		this->pairEneTmp[i] = 0.0;
	}

	for(int i=0;i<4;i++){
		this->pairEneCG[i] = 0.0;
		this->pairEneCGTmp[i] = 0.0;
	}

	this->samplingFreq = 1.0;
}

NuEdge::NuEdge(NuNode* nodeA, NuNode* nodeB, int sep, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib){
	this->graph = NULL;
	this->indexA = nodeA->seqID;
	this->indexB = nodeB->seqID;
	this->nodeA = nodeA;
	this->nodeB = nodeB;
	this->cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->cmTmp = this->cm;
	this->sep = sep;
	int typeA = nodeA->baseType%4;
	int typeB = nodeB->baseType%4;

	this->ei = new EdgeInformation(sep, typeA, typeB, pairLib);
	this->moveSet = new MixedNuPairCluster(sep, typeA*4+typeB, moveLib);
	this->moveSet->updateEdgeInformation(this->ei);
	this->weight = this->ei->weight;
	this->weightRand = weight;

	this->pairLib = pairLib;

	for(int i=0;i<9;i++){
		this->pairEne[i] = 0.0;
		this->pairEneTmp[i] = 0.0;
	}

	for(int i=0;i<4;i++){
		this->pairEneCG[i] = 0.0;
		this->pairEneCGTmp[i] = 0.0;
	}
	this->samplingFreq = this->ei->validClusterNum*1.0;

}

NuEdge::~NuEdge(){
	delete this->ei;
	delete this->moveSet;
}

void NuEdge::initNearNativeMoveSet(){

	BaseDistanceMatrix dm(nodeA->baseConf->cs1, nodeB->baseConf->cs1);
	vector<int> neighborClusters;
	vector<double> distanceToClusterCenters;

	pairLib->getNeighborClusters(dm, nodeA->baseType, nodeB->baseType, sep, neighborClusters, distanceToClusterCenters);

	vector<double> pList;

	if(neighborClusters.size() > 0) {
		double pSum = 0;
		for(int i=0;i<distanceToClusterCenters.size();i++) {
			double d = distanceToClusterCenters[i];
			pSum += exp(-d*2);
		}

		for(int i=0;i<distanceToClusterCenters.size();i++) {
			double d = distanceToClusterCenters[i];
			pList.push_back(exp(-d*2)/pSum);
		}
	}

	this->ei->setClusterList(neighborClusters, pList, pairLib);
	this->moveSet->updateEdgeInformation(ei);
	this->weight = this->ei->weight;
	this->weightRand = this->weight;
	this->samplingFreq = this->ei->validClusterNum;
}

void NuEdge::fixNaiveMove(){
	CsMove cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->ei->setFixed(cm);
	this->moveSet->fixNativeMove(cm);
	this->weight = this->ei->weight;
	this->weightRand = this->weight;
	this->samplingFreq = 0.0;
}

void NuEdge::updateEdgeInfo(NuTree* tree){

	int i,j,k;
	if(tree->adjMtx[indexA*tree->graph->seqLen+indexB]) {
		//BFS algorithm
		this->nodeListA.clear();
		this->edgeListA.clear();
		this->nodeListB.clear();
		this->edgeListB.clear();
		this->connectionBreakPoints.clear();
		this->phoGroupA.clear();
		this->phoGroupB.clear();
		this->phoGroupC.clear();

		queue<NuNode*> qA;
		queue<NuNode*> qB;
		NuNode* vn;
		NuNode* vw;

		bool visited[graph->seqLen];
		for(i=0;i<graph->seqLen;i++){
			visited[i] = false;
		}

		qA.push(nodeA);
		qB.push(nodeB);

		nodeListA.push_back(nodeA);
		nodeListB.push_back(nodeB);

		visited[nodeA->seqID] = true;
		visited[nodeB->seqID] = true;

		while(!qA.empty()){
			vn = qA.front();
			qA.pop();

			for(i=0;i<vn->neighborList.size();i++){
				if(visited[vn->neighborList[i]]) continue;

				vw = graph->allNodes[vn->neighborList[i]];
				qA.push(vw);
				nodeListA.push_back(vw);
				edgeListA.push_back(graph->allEdges[vn->seqID*graph->seqLen+vw->seqID]);
				visited[vw->seqID] = true;
			}
		}

		while(!qB.empty()){
			vn = qB.front();
			qB.pop();
			for(i=0;i<vn->neighborList.size();i++){
				if(visited[vn->neighborList[i]]) continue;
				vw = graph->allNodes[vn->neighborList[i]];
				qB.push(vw);
				nodeListB.push_back(vw);
				edgeListB.push_back(graph->allEdges[vn->seqID*graph->seqLen+vw->seqID]);
				visited[vw->seqID] = true;
			}
		}

		int nodeGroup[graph->seqLen];
		for(i=0;i<nodeListA.size();i++){
			nodeGroup[nodeListA[i]->seqID] = 0;
		}
		for(i=0;i<nodeListB.size();i++){
			nodeGroup[nodeListB[i]->seqID] = 1;
		}

		for(i=0;i<graph->seqLen-1;i++){
			if(graph->connectToDownstream[i] && nodeGroup[i] != nodeGroup[i+1])
				connectionBreakPoints.push_back(i);
		}

		for(i=0;i<connectionBreakPoints.size();i++){
			nodeGroup[connectionBreakPoints[i]] = 2;
		}

		for(i=0;i<graph->seqLen;i++){
			if(nodeGroup[i] == 0 && graph->connectToDownstream[i])
				this->phoGroupA.push_back(graph->allNodes[i]);
			else if(nodeGroup[i] == 1 && graph->connectToDownstream[i])
				this->phoGroupB.push_back(graph->allNodes[i]);
			else if(nodeGroup[i] == 2 && graph->connectToDownstream[i])
				this->phoGroupC.push_back(graph->allNodes[i]);
		}
	}

	int sepR = -sep;
	if(abs(sepR) > 1) sepR = 2;

	if(nodeA->seqID < nodeB->seqID)
		pairEne[0] = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, graph->et);
	else
		pairEne[0] = nuBaseBaseEnergy(nodeB->baseConf, nodeA->baseConf, sepR, graph->et);

	pairEne[1] = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, graph->et);
	pairEne[2] = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, graph->et);

	pairEne[3] = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, graph->et);
	
	pairEne[4] = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, graph->et);
	
	pairEne[5] = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, graph->et);
	pairEne[6] = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, graph->et);
	pairEne[7] = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, graph->et);
	pairEne[8] = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, graph->et);

	for(i=0;i<9;i++){
		pairEneTmp[i] = pairEne[i];
	}
}

void NuEdge::updateEdgeInfoCG(NuTree* tree){

	int i,j,k;
	if(tree->adjMtx[indexA*tree->graph->seqLen+indexB]) {
		//BFS algorithm
		this->nodeListA.clear();
		this->edgeListA.clear();
		this->nodeListB.clear();
		this->edgeListB.clear();
		this->connectionBreakPoints.clear();
		this->phoGroupA.clear();
		this->phoGroupB.clear();
		this->phoGroupC.clear();

		queue<NuNode*> qA;
		queue<NuNode*> qB;
		NuNode* vn;
		NuNode* vw;

		bool visited[graph->seqLen];
		for(i=0;i<graph->seqLen;i++){
			visited[i] = false;
		}

		qA.push(nodeA);
		qB.push(nodeB);

		nodeListA.push_back(nodeA);
		nodeListB.push_back(nodeB);

		visited[nodeA->seqID] = true;
		visited[nodeB->seqID] = true;

		while(!qA.empty()){
			vn = qA.front();
			qA.pop();

			for(i=0;i<vn->neighborList.size();i++){
				if(visited[vn->neighborList[i]]) continue;

				vw = graph->allNodes[vn->neighborList[i]];
				qA.push(vw);
				nodeListA.push_back(vw);
				edgeListA.push_back(graph->allEdges[vn->seqID*graph->seqLen+vw->seqID]);
				visited[vw->seqID] = true;
			}
		}

		while(!qB.empty()){
			vn = qB.front();
			qB.pop();
			for(i=0;i<vn->neighborList.size();i++){
				if(visited[vn->neighborList[i]]) continue;
				vw = graph->allNodes[vn->neighborList[i]];
				qB.push(vw);
				nodeListB.push_back(vw);
				edgeListB.push_back(graph->allEdges[vn->seqID*graph->seqLen+vw->seqID]);
				visited[vw->seqID] = true;
			}
		}

		int nodeGroup[graph->seqLen];
		for(i=0;i<nodeListA.size();i++){
			nodeGroup[nodeListA[i]->seqID] = 0;
		}
		for(i=0;i<nodeListB.size();i++){
			nodeGroup[nodeListB[i]->seqID] = 1;
		}

		for(i=0;i<graph->seqLen-1;i++){
			if(graph->connectToDownstream[i] && nodeGroup[i] != nodeGroup[i+1])
				connectionBreakPoints.push_back(i);
		}

		for(i=0;i<connectionBreakPoints.size();i++){
			nodeGroup[connectionBreakPoints[i]] = 2;
		}

		for(i=0;i<graph->seqLen;i++){
			if(nodeGroup[i] == 0 && graph->connectToDownstream[i])
				this->phoGroupA.push_back(graph->allNodes[i]);
			else if(nodeGroup[i] == 1 && graph->connectToDownstream[i])
				this->phoGroupB.push_back(graph->allNodes[i]);
			else if(nodeGroup[i] == 2 && graph->connectToDownstream[i])
				this->phoGroupC.push_back(graph->allNodes[i]);
		}
	}

	int sepR = -sep;
	if(abs(sepR) > 1) sepR = 2;

	if(nodeA->seqID < nodeB->seqID)
		pairEneCG[0] = nuBaseBaseEnergyCG(nodeA->baseConfCG, nodeB->baseConfCG, sep, graph->et);
	else
		pairEneCG[0] = nuBaseBaseEnergyCG(nodeB->baseConfCG, nodeA->baseConfCG, sepR, graph->et);

	pairEneCG[1] = nuBaseRiboseEnergyCG(nodeA->baseConfCG, nodeB->riboseConfCG, sep, graph->et);
	pairEneCG[2] = nuBaseRiboseEnergyCG(nodeB->baseConfCG, nodeA->riboseConfCG, sepR, graph->et);
	pairEneCG[3] = nuRiboseRiboseEnergyCG(nodeA->riboseConfCG, nodeB->riboseConfCG, sep, graph->et);

	for(i=0;i<4;i++){
		pairEneCGTmp[i] = pairEneCG[i];
	}
}

void NuEdge::updateCsMove(CsMove& cm){
	int i,j,sep, sepR;
	this->cmTmp = cm;
	LocalFrame cs1 = nodeA->baseConfTmp->cs1 + cm;
	nodeB->updateCoordinate(cs1);
	for(i=0;i<edgeListB.size();i++){
		cs1 = edgeListB[i]->nodeA->baseConfTmp->cs1 + edgeListB[i]->cmTmp;
		edgeListB[i]->nodeB->updateCoordinate(cs1);
	}

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->et->pb->buildPhosphate(graph->allNodes[j]->riboseConfTmp, graph->allNodes[j+1]->riboseConfTmp, graph->allNodes[j]->phoConfTmp);
		graph->allNodes[j]->eneTmp = graph->allNodes[j]->riboseConfTmp->rot->energy + graph->allNodes[j]->phoConfTmp->ene;
	}

	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;
	NuEdge* eB;


	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nA->seqID < nB->seqID)
				eA->pairEneTmp[0] = nuBaseBaseEnergy(nA->baseConfTmp, nB->baseConfTmp, sep, graph->et);
			else
				eA->pairEneTmp[0] = nuBaseBaseEnergy(nB->baseConfTmp, nA->baseConfTmp, sepR, graph->et);

			eA->pairEneTmp[1] = nuBaseRiboseEnergy(nA->baseConfTmp, nB->riboseConfTmp, sep, graph->et);
			eA->pairEneTmp[3] = nuBaseRiboseEnergy(nB->baseConfTmp, nA->riboseConfTmp, sepR, graph->et);

			if(nA->seqID < nB->seqID)
				eA->pairEneTmp[4] = nuRiboseRiboseEnergy(nA->riboseConfTmp, nB->riboseConfTmp, sep, graph->et);
			else
				eA->pairEneTmp[4] = nuRiboseRiboseEnergy(nB->riboseConfTmp, nA->riboseConfTmp, sepR, graph->et);

			eB->pairEneTmp[0] = eA->pairEneTmp[0];
			eB->pairEneTmp[1] = eA->pairEneTmp[3];
			eB->pairEneTmp[3] = eA->pairEneTmp[1];
			eB->pairEneTmp[4] = eA->pairEneTmp[4];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListB
	 */

	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et);
			}

			eB->pairEneTmp[6] = eA->pairEneTmp[2];
			eB->pairEneTmp[7] = eA->pairEneTmp[5];

		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListC
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et);
			}
			eB->pairEneTmp[6] = eA->pairEneTmp[2];
			eB->pairEneTmp[7] = eA->pairEneTmp[5];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListA
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupA.size();j++){
			nB = phoGroupA[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et);
			}
			eB->pairEneTmp[6] = eA->pairEneTmp[2];
			eB->pairEneTmp[7] = eA->pairEneTmp[5];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListC
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et);
			}
			eB->pairEneTmp[6] = eA->pairEneTmp[2];
			eB->pairEneTmp[7] = eA->pairEneTmp[5];
		}
	}

	/*
	 * interaction energy between phoListA and phoListB
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nA->seqID < nB->seqID)
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et);

			eB->pairEneTmp[8] = eA->pairEneTmp[8];
		}
	}

	/*
	 * interaction energy between phoListA and phoListC
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nA->seqID < nB->seqID)
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et);

			eB->pairEneTmp[8] = eA->pairEneTmp[8];
		}
	}

	/*
	 * interaction energy between phoListB and phoListC
	 */
	for(i=0;i<phoGroupB.size();i++){
		nA = phoGroupB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];
			if(nA->seqID < nB->seqID)
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et);
			eB->pairEneTmp[8] = eA->pairEneTmp[8];
		}
	}

	/*
	 * interaction energy between  phoListC and phoListC
	 */
	for(i=0;i<phoGroupC.size();i++){
		nA = phoGroupC[i];
		for(j=i+1;j<phoGroupC.size();j++){
			nB = phoGroupC[j];


			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];


			if(nA->seqID < nB->seqID)
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et);

			eB->pairEneTmp[8] = eA->pairEneTmp[8];
		}
	}

}

double NuEdge::mutEnergy(){
	int i,j;
	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;
	double mutE = 0.0;

	for(i=0;i<phoGroupC.size();i++){
		mutE += phoGroupC[i]->eneTmp - phoGroupC[i]->ene;
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];

			mutE += eA->pairEneTmp[0] - eA->pairEne[0];
			mutE += eA->pairEneTmp[1] - eA->pairEne[1];
			mutE += eA->pairEneTmp[3] - eA->pairEne[3];
			mutE += eA->pairEneTmp[4] - eA->pairEne[4];

		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListB
	 */

	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];

			if(nB->connectToNeighbor) {
				mutE += eA->pairEneTmp[2] - eA->pairEne[2];
				mutE += eA->pairEneTmp[5] - eA->pairEne[5];

			}
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListC
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];

			if(nB->connectToNeighbor) {
				mutE += eA->pairEneTmp[2] - eA->pairEne[2];
				mutE += eA->pairEneTmp[5] - eA->pairEne[5];


			}
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListA
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupA.size();j++){
			nB = phoGroupA[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];

			if(nB->connectToNeighbor) {
				mutE += eA->pairEneTmp[2] - eA->pairEne[2];
				mutE += eA->pairEneTmp[5] - eA->pairEne[5];
			}
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListC
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];

			if(nB->connectToNeighbor) {
				mutE += eA->pairEneTmp[2] - eA->pairEne[2];
				mutE += eA->pairEneTmp[5] - eA->pairEne[5];

			}

		}
	}

	/*
	 * interaction energy between phoListA and phoListB
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			mutE += eA->pairEneTmp[8] - eA->pairEne[8];

		}
	}
	/*
	 * interaction energy between phoListA and phoListC
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			mutE += eA->pairEneTmp[8] - eA->pairEne[8];
		}
	}
	/*
	 * interaction energy between phoListB and phoListC
	 */
	for(i=0;i<phoGroupB.size();i++){
		nA = phoGroupB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			mutE += eA->pairEneTmp[8] - eA->pairEne[8];

		}
	}

	/*
	 * interaction energy between  phoListC and phoListC
	 */
	for(i=0;i<phoGroupC.size();i++){
		nA = phoGroupC[i];
		for(j=i+1;j<phoGroupC.size();j++){
			nB = phoGroupC[j];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			mutE += eA->pairEneTmp[8] - eA->pairEne[8];

		}
	}

	return mutE;
}

void NuEdge::acceptMutation(){
	int i,j;
	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;
	NuEdge* eB;
	int sep, sepR;

	this->cm = this->cmTmp;
	for(i=0;i<nodeListB.size();i++){
		nodeListB[i]->acceptCoordMove();
	}

	for(i=0;i<phoGroupB.size();i++){
		phoGroupB[i]->phoConf->copyValueFrom(phoGroupB[i]->phoConfTmp);
	}

	for(i=0;i<phoGroupC.size();i++){
		phoGroupC[i]->ene = phoGroupC[i]->eneTmp;
		phoGroupC[i]->phoConf->copyValueFrom(phoGroupC[i]->phoConfTmp);
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			eA->pairEne[0] = eA->pairEneTmp[0];
			eA->pairEne[1] = eA->pairEneTmp[1];
			eA->pairEne[3] = eA->pairEneTmp[3];
			eA->pairEne[4] = eA->pairEneTmp[4];

			eB->pairEne[0] = eB->pairEneTmp[0];
			eB->pairEne[1] = eB->pairEneTmp[1];
			eB->pairEne[3] = eB->pairEneTmp[3];
			eB->pairEne[4] = eB->pairEneTmp[4];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListB
	 */

	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEne[2] = eA->pairEneTmp[2];
				eA->pairEne[5] = eA->pairEneTmp[5];
			}

			eB->pairEne[6] = eB->pairEneTmp[6];
			eB->pairEne[7] = eB->pairEneTmp[7];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListC
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEne[2] = eA->pairEneTmp[2];
				eA->pairEne[5] = eA->pairEneTmp[5];
			}

			eB->pairEne[6] = eB->pairEneTmp[6];
			eB->pairEne[7] = eB->pairEneTmp[7];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListA
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupA.size();j++){
			nB = phoGroupA[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEne[2] = eA->pairEneTmp[2];
				eA->pairEne[5] = eA->pairEneTmp[5];
			}

			eB->pairEne[6] = eB->pairEneTmp[6];
			eB->pairEne[7] = eB->pairEneTmp[7];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListC
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEne[2] = eA->pairEneTmp[2];
				eA->pairEne[5] = eA->pairEneTmp[5];
			}

			eB->pairEne[6] = eB->pairEneTmp[6];
			eB->pairEne[7] = eB->pairEneTmp[7];
		}
	}

	/*
	 * interaction energy between phoListA and phoListB
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];
			eA->pairEne[8] = eA->pairEneTmp[8];
			eB->pairEne[8] = eB->pairEneTmp[8];
		}
	}
	/*
	 * interaction energy between phoListA and phoListC
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];
			eA->pairEne[8] = eA->pairEneTmp[8];
			eB->pairEne[8] = eB->pairEneTmp[8];
		}
	}
	/*
	 * interaction energy between phoListB and phoListC
	 */
	for(i=0;i<phoGroupB.size();i++){
		nA = phoGroupB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];
			eA->pairEne[8] = eA->pairEneTmp[8];
			eB->pairEne[8] = eB->pairEneTmp[8];
		}
	}

	/*
	 * interaction energy between  phoListC and phoListC
	 */
	for(i=0;i<phoGroupC.size();i++){
		nA = phoGroupC[i];
		for(j=i+1;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			eA->pairEne[8] = eA->pairEneTmp[8];
			eB->pairEne[8] = eB->pairEneTmp[8];
		}
	}
}

void NuEdge::clearMutation(){
	this->cmTmp = this->cm;
	int i,j;
	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;
	NuEdge* eB;
	int sep, sepR;

	for(i=0;i<nodeListB.size();i++){
		nodeListB[i]->clearCoordMove();
	}

	for(i=0;i<phoGroupB.size();i++){
		phoGroupB[i]->phoConfTmp->copyValueFrom(phoGroupB[i]->phoConf);
	}

	for(i=0;i<phoGroupC.size();i++){
		phoGroupC[i]->eneTmp = phoGroupC[i]->ene;
		phoGroupC[i]->phoConfTmp->copyValueFrom(phoGroupC[i]->phoConf);
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			eA->pairEneTmp[0] = eA->pairEne[0];
			eA->pairEneTmp[1] = eA->pairEne[1];
			eA->pairEneTmp[3] = eA->pairEne[3];
			eA->pairEneTmp[4] = eA->pairEne[4];

			eB->pairEneTmp[0] = eB->pairEne[0];
			eB->pairEneTmp[1] = eB->pairEne[1];
			eB->pairEneTmp[3] = eB->pairEne[3];
			eB->pairEneTmp[4] = eB->pairEne[4];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListB
	 */

	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = eA->pairEne[2];
				eA->pairEneTmp[5] = eA->pairEne[5];
			}

			eB->pairEneTmp[6] = eB->pairEne[6];
			eB->pairEneTmp[7] = eB->pairEne[7];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and phoListC
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = eA->pairEne[2];
				eA->pairEneTmp[5] = eA->pairEne[5];
			}

			eB->pairEneTmp[6] = eB->pairEne[6];
			eB->pairEneTmp[7] = eB->pairEne[7];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListA
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupA.size();j++){
			nB = phoGroupA[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = eA->pairEne[2];
				eA->pairEneTmp[5] = eA->pairEne[5];
			}

			eB->pairEneTmp[6] = eB->pairEne[6];
			eB->pairEneTmp[7] = eB->pairEne[7];
		}
	}

	/*
	 * interaction energy between (base, ribose) of nodeListB and phoListC
	 */
	for(i=0;i<nodeListB.size();i++){
		nA = nodeListB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nB->connectToNeighbor) {
				eA->pairEneTmp[2] = eA->pairEne[2];
				eA->pairEneTmp[5] = eA->pairEne[5];
			}

			eB->pairEneTmp[6] = eB->pairEne[6];
			eB->pairEneTmp[7] = eB->pairEne[7];
		}
	}

	/*
	 * interaction energy between phoListA and phoListB
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupB.size();j++){
			nB = phoGroupB[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];
			eA->pairEneTmp[8] = eA->pairEne[8];
			eB->pairEneTmp[8] = eB->pairEne[8];
		}
	}
	/*
	 * interaction energy between phoListA and phoListC
	 */
	for(i=0;i<phoGroupA.size();i++){
		nA = phoGroupA[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];
			eA->pairEneTmp[8] = eA->pairEne[8];
			eB->pairEneTmp[8] = eB->pairEne[8];
		}
	}
	/*
	 * interaction energy between phoListB and phoListC
	 */
	for(i=0;i<phoGroupB.size();i++){
		nA = phoGroupB[i];
		for(j=0;j<phoGroupC.size();j++){
			nB = phoGroupC[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];
			eA->pairEneTmp[8] = eA->pairEne[8];
			eB->pairEneTmp[8] = eB->pairEne[8];
		}
	}

	/*
	 * interaction energy between  phoListC and phoListC
	 */
	for(i=0;i<phoGroupC.size();i++){
		nA = phoGroupC[i];
		for(j=i+1;j<phoGroupC.size();j++){
			nB = phoGroupC[j];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			eA->pairEneTmp[8] = eA->pairEne[8];
			eB->pairEneTmp[8] = eB->pairEne[8];
		}
	}
}

void NuEdge::updateCsMoveCG(CsMove& cm) {

	int i,j,sep, sepR;

	this->cmTmp = cm;
	LocalFrame cs1 = nodeA->baseConfCGTmp->cs1 + cm;
	nodeB->updateCoordinateCG(cs1);

	for(i=0;i<edgeListB.size();i++){
		
		cs1 = edgeListB[i]->nodeA->baseConfCGTmp->cs1 + edgeListB[i]->cmTmp;
		edgeListB[i]->nodeB->updateCoordinateCG(cs1);
	}

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->bbcgTmp = nuConnectionEnergyCG(graph->allNodes[j]->riboseConfCGTmp, graph->allNodes[j+1]->riboseConfCGTmp, graph->et);
		graph->allNodes[j]->eneCGTmp = graph->allNodes[j]->riboseConfCGTmp->rot->energy + graph->allNodes[j]->bbcgTmp;
	}

	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;
	NuEdge* eB;

	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			if(nA->seqID < nB->seqID)
				eA->pairEneCGTmp[0] = nuBaseBaseEnergyCG(nA->baseConfCGTmp, nB->baseConfCGTmp, sep, graph->et);
			else
				eA->pairEneCGTmp[0] = nuBaseBaseEnergyCG(nB->baseConfCGTmp, nA->baseConfCGTmp, sepR, graph->et);

			eA->pairEneCGTmp[1] = nuBaseRiboseEnergyCG(nA->baseConfCGTmp, nB->riboseConfCGTmp, sep, graph->et);
			eA->pairEneCGTmp[2] = nuBaseRiboseEnergyCG(nB->baseConfCGTmp, nA->riboseConfCGTmp, sepR, graph->et);

			if(nA->seqID < nB->seqID)
				eA->pairEneCGTmp[3] = nuRiboseRiboseEnergyCG(nA->riboseConfCGTmp, nB->riboseConfCGTmp, sep, graph->et);
			else
				eA->pairEneCGTmp[3] = nuRiboseRiboseEnergyCG(nB->riboseConfCGTmp, nA->riboseConfCGTmp, sepR, graph->et);

			eB->pairEneCGTmp[0] = eA->pairEneCGTmp[0];
			eB->pairEneCGTmp[1] = eA->pairEneCGTmp[2];
			eB->pairEneCGTmp[2] = eA->pairEneCGTmp[1];
			eB->pairEneCGTmp[3] = eA->pairEneCGTmp[3];
		}
	}
}

double NuEdge::mutEnergyCG(){
	int i,j;
	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;

	double mutE = 0.0;


	for(i=0;i<phoGroupC.size();i++){
		mutE += phoGroupC[i]->eneCGTmp - phoGroupC[i]->eneCG;
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];
			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			mutE += eA->pairEneCGTmp[0] - eA->pairEneCG[0];
			mutE += eA->pairEneCGTmp[1] - eA->pairEneCG[1];
			mutE += eA->pairEneCGTmp[2] - eA->pairEneCG[2];
			mutE += eA->pairEneCGTmp[3] - eA->pairEneCG[3];
		}
	}
	return mutE;
}

void NuEdge::acceptMutationCG(){
	int i,j;
	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;
	NuEdge* eB;
	int sep, sepR;

	this->cm = this->cmTmp;
	for(i=0;i<nodeListB.size();i++){
		nodeListB[i]->acceptCoordMoveCG();
	}

	for(i=0;i<phoGroupC.size();i++){
		phoGroupC[i]->eneCG = phoGroupC[i]->eneCGTmp;
		phoGroupC[i]->bbcg = phoGroupC[i]->bbcgTmp;
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];
			sep = graph->sepTable[nA->seqID*graph->seqLen+nB->seqID];
			sepR = graph->sepTable[nB->seqID*graph->seqLen+nA->seqID];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			eA->pairEneCG[0] = eA->pairEneCGTmp[0];
			eA->pairEneCG[1] = eA->pairEneCGTmp[1];
			eA->pairEneCG[2] = eA->pairEneCGTmp[2];
			eA->pairEneCG[3] = eA->pairEneCGTmp[3];

			eB->pairEneCG[0] = eB->pairEneCGTmp[0];
			eB->pairEneCG[1] = eB->pairEneCGTmp[1];
			eB->pairEneCG[2] = eB->pairEneCGTmp[2];
			eB->pairEneCG[3] = eB->pairEneCGTmp[3];
		}
	}
}

void NuEdge::clearMutationCG(){
	int i,j;
	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;
	NuEdge* eB;
	int sep, sepR;

	this->cmTmp = this->cm;
	for(i=0;i<nodeListB.size();i++){
		nodeListB[i]->clearCoordMoveCG();
	}

	for(i=0;i<phoGroupC.size();i++){
		phoGroupC[i]->eneCGTmp = phoGroupC[i]->eneCG;
		phoGroupC[i]->bbcgTmp = phoGroupC[i]->bbcg;
	}

	/*
	 * interaction energy between (base, ribose) of nodeListA and (base, ribose) of nodeListB
	 */
	for(i=0;i<nodeListA.size();i++){
		nA = nodeListA[i];
		for(j=0;j<nodeListB.size();j++){
			nB = nodeListB[j];

			eA = graph->allEdges[nA->seqID*graph->seqLen+nB->seqID];
			eB = graph->allEdges[nB->seqID*graph->seqLen+nA->seqID];

			eA->pairEneCGTmp[0] = eA->pairEneCG[0];
			eA->pairEneCGTmp[1] = eA->pairEneCG[1];
			eA->pairEneCGTmp[2] = eA->pairEneCG[2];
			eA->pairEneCGTmp[3] = eA->pairEneCG[3];

			eB->pairEneCGTmp[0] = eB->pairEneCG[0];
			eB->pairEneCGTmp[1] = eB->pairEneCG[1];
			eB->pairEneCGTmp[2] = eB->pairEneCG[2];
			eB->pairEneCGTmp[3] = eB->pairEneCG[3];
		}
	}
}

bool NuEdge::checkEnergy(){

	bool tag = true;
	int sepR = -sep;
	if(abs(sepR) > 1)
		sepR = 2;

	double e0;
	if(nodeA->seqID < nodeB->seqID)
		e0 = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, graph->et);
	else
		e0 = nuBaseBaseEnergy(nodeB->baseConf, nodeA->baseConf, sepR, graph->et);

	double e1 = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, graph->et);
	double e2 = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, graph->et);
	double e3 = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, graph->et);
	double e4 = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, graph->et);
	double e5 = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, graph->et);
	double e6 = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, graph->et);
	double e7 = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, graph->et);
	double e8 = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, graph->et);

	double e0Tmp;
	if(nodeA->seqID < nodeB->seqID)
		e0Tmp = nuBaseBaseEnergy(nodeA->baseConfTmp, nodeB->baseConfTmp, sep, graph->et);
	else
		e0Tmp = nuBaseBaseEnergy(nodeB->baseConfTmp, nodeA->baseConfTmp, sepR, graph->et);

	double e1Tmp = nuBaseRiboseEnergy(nodeA->baseConfTmp, nodeB->riboseConfTmp, sep, graph->et);
	double e2Tmp = nuBasePhoEnergy(nodeA->baseConfTmp, nodeB->phoConfTmp, sep, graph->et);
	double e3Tmp = nuBaseRiboseEnergy(nodeB->baseConfTmp, nodeA->riboseConfTmp, sepR, graph->et);
	double e4Tmp = nuRiboseRiboseEnergy(nodeA->riboseConfTmp, nodeB->riboseConfTmp, sep, graph->et);
	double e5Tmp = nuRibosePhoEnergy(nodeA->riboseConfTmp, nodeB->phoConfTmp, sep, graph->et);
	double e6Tmp = nuBasePhoEnergy(nodeB->baseConfTmp, nodeA->phoConfTmp, sepR, graph->et);
	double e7Tmp = nuRibosePhoEnergy(nodeB->riboseConfTmp, nodeA->phoConfTmp, sepR, graph->et);
	double e8Tmp = nuPhoPhoEnergy(nodeA->phoConfTmp, nodeB->phoConfTmp, sep, graph->et);

	if(abs(e0 - this->pairEne[0]) > 0.001){
		printf("Edge %3d-%3d ERROR BB %7.3f %7.3f\n", this->indexA, this->indexB, e0, this->pairEne[0]);
		tag = false;
	}
	if(abs(e1 - this->pairEne[1]) > 0.001){
		printf("Edge %3d-%3d ERROR BR %7.3f %7.3f\n", this->indexA, this->indexB,e1, this->pairEne[1]);
		tag = false;
	}
	if(abs(e2 - this->pairEne[2]) > 0.001){
		printf("Edge %3d-%3d ERROR BP %7.3f %7.3f\n", this->indexA, this->indexB,e2, this->pairEne[2]);
		tag = false;
	}
	if(abs(e3 - this->pairEne[3]) > 0.001){
		printf("Edge %3d-%3d ERROR RB %7.3f %7.3f\n", this->indexA, this->indexB,e3, this->pairEne[3]);
		tag = false;
	}
	if(abs(e4 - this->pairEne[4]) > 0.001){
		printf("Edge %3d-%3d ERROR RR %7.3f %7.3f\n", this->indexA, this->indexB,e4, this->pairEne[4]);
		tag = false;
	}
	if(abs(e5 - this->pairEne[5]) > 0.001){
		printf("Edge %3d-%3d ERROR RP %7.3f %7.3f\n",this->indexA, this->indexB, e5, this->pairEne[5]);
		tag = false;
	}
	if(abs(e6 - this->pairEne[6]) > 0.001){
		printf("Edge %3d-%3d ERROR PB %7.3f %7.3f\n", this->indexA, this->indexB,e6, this->pairEne[6]);
		tag = false;
	}
	if(abs(e7 - this->pairEne[7]) > 0.001){
		printf("Edge %3d-%3d ERROR PR %7.3f %7.3f\n", this->indexA, this->indexB,e7, this->pairEne[7]);
		tag = false;
	}
	if(abs(e8 - this->pairEne[8]) > 0.001){
		printf("Edge %3d-%3d ERROR PP %7.3f %7.3f\n", this->indexA, this->indexB,e8, this->pairEne[8]);
		tag = false;
	}

	if(abs(e0Tmp - this->pairEneTmp[0]) > 0.001){
		printf("Edge %3d-%3d ERROR BBTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e0Tmp, this->pairEneTmp[0], this->pairEne[0]);
		tag = false;
	}
	if(abs(e1Tmp - this->pairEneTmp[1]) > 0.001){
		printf("Edge %3d-%3d ERROR BRTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e1Tmp, this->pairEneTmp[1], this->pairEne[1]);
		tag = false;
	}
	if(abs(e2Tmp - this->pairEneTmp[2]) > 0.001){
		printf("Edge %3d-%3d ERROR BPTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e2Tmp, this->pairEneTmp[2], this->pairEne[2]);
		tag = false;
	}
	if(abs(e3Tmp - this->pairEneTmp[3]) > 0.001){
		printf("Edge %3d-%3d ERROR RBTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e3Tmp, this->pairEneTmp[3], this->pairEne[3]);
		tag = false;
	}
	if(abs(e4Tmp - this->pairEneTmp[4]) > 0.001){
		printf("Edge %3d-%3d ERROR RRTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e4Tmp, this->pairEneTmp[4], this->pairEne[4]);
		tag = false;
	}
	if(abs(e5Tmp - this->pairEneTmp[5]) > 0.001){
		printf("Edge %3d-%3d ERROR RPTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e5Tmp, this->pairEneTmp[5], this->pairEne[5]);
		tag = false;
	}
	if(abs(e6Tmp - this->pairEneTmp[6]) > 0.001){
		printf("Edge %3d-%3d ERROR PBTmp %7.3f %7.3f %7.3f\n",this->indexA, this->indexB, e6Tmp, this->pairEneTmp[6], this->pairEne[6]);
		tag = false;
	}
	if(abs(e7Tmp - this->pairEneTmp[7]) > 0.001){
		printf("Edge %3d-%3d ERROR PRTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e7Tmp, this->pairEneTmp[7], this->pairEne[7]);
		tag = false;
	}
	if(abs(e8Tmp - this->pairEneTmp[8]) > 0.001){
		printf("Edge %3d-%3d ERROR PPTmp %7.3f %7.3f %7.3f\n", this->indexA, this->indexB,e8Tmp, this->pairEneTmp[8], this->pairEne[8]);
		tag = false;
	}
	return tag;
}

bool NuEdge::checkReversePair(){
	NuEdge* eB = graph->allEdges[indexB*graph->seqLen+indexA];

	bool tag = true;
	if(abs(pairEne[0] - eB->pairEne[0]) > 0.001 || abs(pairEneTmp[0] - eB->pairEneTmp[0]) > 0.001)   {
		printf("rev error: %2d %2d BB %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[0], eB->pairEne[0],pairEneTmp[0],eB->pairEneTmp[0]);
		tag = false;
	}

	if(abs(pairEne[1] - eB->pairEne[3]) > 0.001 || abs(pairEneTmp[1] - eB->pairEneTmp[3]) > 0.001)   {
		printf("rev error: %2d %2d BR %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[1], eB->pairEne[3],pairEneTmp[1],eB->pairEneTmp[3]);
		tag = false;
	}

	if(abs(pairEne[2] - eB->pairEne[6]) > 0.001 || abs(pairEneTmp[2] - eB->pairEneTmp[6]) > 0.001)   {
		printf("rev error: %2d %2d BP %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[2], eB->pairEne[6],pairEneTmp[2],eB->pairEneTmp[6]);
		tag = false;
	}

	if(abs(pairEne[3] - eB->pairEne[1]) > 0.001 || abs(pairEneTmp[3] - eB->pairEneTmp[1]) > 0.001)   {
		printf("rev error: %2d %2d RB %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[3], eB->pairEne[1],pairEneTmp[3],eB->pairEneTmp[1]);
		tag = false;
	}

	if(abs(pairEne[4] - eB->pairEne[4]) > 0.001 || abs(pairEneTmp[4] - eB->pairEneTmp[4]) > 0.001)   {
		cout << "SEP: " << this->sep << " " << eB->sep << endl;
		printf("rev error: %2d %2d RR %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[4], eB->pairEne[4],pairEneTmp[4],eB->pairEneTmp[4]);
		tag = false;
	}

	if(abs(pairEne[5] - eB->pairEne[7]) > 0.001 || abs(pairEneTmp[5] - eB->pairEneTmp[7]) > 0.001)   {
		printf("rev error: %2d %2d RP %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[5], eB->pairEne[7],pairEneTmp[5],eB->pairEneTmp[7]);
		tag = false;
	}

	if(abs(pairEne[6] - eB->pairEne[2]) > 0.001 || abs(pairEneTmp[6] - eB->pairEneTmp[2]) > 0.001)   {
		printf("rev error: %2d %2d PB %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[6], eB->pairEne[2],pairEneTmp[6],eB->pairEneTmp[2]);
		tag = false;
	}

	if(abs(pairEne[7] - eB->pairEne[5]) > 0.001 || abs(pairEneTmp[7] - eB->pairEneTmp[5]) > 0.001)   {
		printf("rev error: %2d %2d PR %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[7], eB->pairEne[5],pairEneTmp[7],eB->pairEneTmp[5]);
		tag = false;
	}

	if(abs(pairEne[8] - eB->pairEne[8]) > 0.001 || abs(pairEneTmp[8] - eB->pairEneTmp[8]) > 0.001)   {
		printf("rev error: %2d %2d RR %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEne[8], eB->pairEne[8],pairEneTmp[8],eB->pairEneTmp[8]);
		tag = false;
	}
	return tag;

}

bool NuEdge::checkReversePairCG(){
	NuEdge* eB = graph->allEdges[indexB*graph->seqLen+indexA];

	bool tag = true;
	if(abs(pairEneCG[0] - eB->pairEneCG[0]) > 0.001 || abs(pairEneCGTmp[0] - eB->pairEneCGTmp[0]) > 0.001)   {
		printf("rev error: %2d %2d BB %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEneCG[0], eB->pairEneCG[0],pairEneCGTmp[0],eB->pairEneCGTmp[0]);
		tag = false;
	}

	if(abs(pairEneCG[1] - eB->pairEneCG[2]) > 0.001 || abs(pairEneCGTmp[1] - eB->pairEneCGTmp[2]) > 0.001)   {
		printf("rev error: %2d %2d BR %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEneCG[1], eB->pairEneCG[2],pairEneCGTmp[1],eB->pairEneCGTmp[2]);
		tag = false;
	}

	if(abs(pairEneCG[2] - eB->pairEneCG[1]) > 0.001 || abs(pairEneCGTmp[2] - eB->pairEneCGTmp[1]) > 0.001)   {
		printf("rev error: %2d %2d RB %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEneCG[2], eB->pairEneCG[1],pairEneCGTmp[2],eB->pairEneCGTmp[1]);
		tag = false;
	}

	if(abs(pairEneCG[3] - eB->pairEneCG[3]) > 0.001 || abs(pairEneCGTmp[3] - eB->pairEneCGTmp[3]) > 0.001)   {
		printf("rev error: %2d %2d RR %7.3f %7.3f %7.3f %7.3f\n", indexA, indexB, pairEneCG[3], eB->pairEneCG[3],pairEneCGTmp[3],eB->pairEneCGTmp[3]);
		tag = false;
	}
	return tag;

}

bool NuEdge::checkEnergyCG(){

	bool tag = true;
	int sepR = -sep;
	if(abs(sepR) > 1)
		sepR = 2;

	double e0;
	if(nodeA->seqID < nodeB->seqID)
		e0 = nuBaseBaseEnergyCG(nodeA->baseConfCG, nodeB->baseConfCG, sep, graph->et);
	else
		e0 = nuBaseBaseEnergyCG(nodeB->baseConfCG, nodeA->baseConfCG, sepR, graph->et);

	double e1 = nuBaseRiboseEnergyCG(nodeA->baseConfCG, nodeB->riboseConfCG, sep, graph->et);
	double e2 = nuBaseRiboseEnergyCG(nodeB->baseConfCG, nodeA->riboseConfCG, sepR, graph->et);
	double e3 = nuRiboseRiboseEnergyCG(nodeA->riboseConfCG, nodeB->riboseConfCG, sep, graph->et);

	double e0Tmp;
	if(nodeA->seqID < nodeB->seqID)
		e0Tmp = nuBaseBaseEnergyCG(nodeA->baseConfCGTmp, nodeB->baseConfCGTmp, sep, graph->et);
	else
		e0Tmp = nuBaseBaseEnergyCG(nodeB->baseConfCGTmp, nodeA->baseConfCGTmp, sepR, graph->et);

	double e1Tmp = nuBaseRiboseEnergyCG(nodeA->baseConfCGTmp, nodeB->riboseConfCGTmp, sep, graph->et);
	double e2Tmp = nuBaseRiboseEnergyCG(nodeB->baseConfCGTmp, nodeA->riboseConfCGTmp, sepR, graph->et);
	double e3Tmp = nuRiboseRiboseEnergyCG(nodeA->riboseConfCGTmp, nodeB->riboseConfCGTmp, sep, graph->et);

	if(abs(e0 - this->pairEneCG[0]) > 0.001){
		printf("ERROR cgBB %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID, e0, this->pairEneCG[0]);
		tag = false;
	}
	if(abs(e1 - this->pairEneCG[1]) > 0.001){
		printf("ERROR cgBR %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID,e1, this->pairEneCG[1]);
		tag = false;
	}
	if(abs(e2 - this->pairEneCG[2]) > 0.001){
		printf("ERROR cgRB %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID,e2, this->pairEneCG[2]);
		tag = false;
	}
	if(abs(e3 - this->pairEneCG[3]) > 0.001){
		printf("ERROR cgRR %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID,e3, this->pairEneCG[3]);
		tag = false;
	}

	if(abs(e0Tmp - this->pairEneCGTmp[0]) > 0.001){
		printf("ERROR cgBBTmp %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID,e0, this->pairEneCGTmp[0]);
		tag = false;
	}
	if(abs(e1Tmp - this->pairEneCGTmp[1]) > 0.001){
		printf("ERROR cgBRTmp %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID,e1, this->pairEneCGTmp[1]);
		tag = false;
	}
	if(abs(e2Tmp - this->pairEneCGTmp[2]) > 0.001){
		printf("ERROR cgRBTmp %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID,e2, this->pairEneCGTmp[2]);
		tag = false;
	}
	if(abs(e3Tmp - this->pairEneCGTmp[3]) > 0.001){
		printf("ERROR cgRRTmp %d %d %7.3f %7.3f\n", nodeA->seqID, nodeB->seqID,e3, this->pairEneCGTmp[3]);
		tag = false;
	}
	return tag;
}

void NuEdge::printPartition(){
	int i,j,k;
	cout << "sep: " << this->sep << endl;
	cout << "sample freq: " << this->nodeA->seqID << " " << this->nodeB->seqID << " " << this->samplingFreq << endl;
	cout << "nodeListA: ";
	for(i=0;i<nodeListA.size();i++){
		cout << nodeListA[i]->seqID << " ";
	}
	cout << endl;

	cout << "nodeListB: ";
	for(i=0;i<nodeListB.size();i++){
		cout << nodeListB[i]->seqID << " ";
	}
	cout << endl;

	cout << "edgeListA: ";
	for(i=0;i<edgeListA.size();i++){
		cout << edgeListA[i]->indexA << "-" << edgeListA[i]->indexB << " ";
	}
	cout << endl;

	cout << "edgeListB: ";
	for(i=0;i<edgeListB.size();i++){
		cout << edgeListB[i]->indexA << "-" << edgeListB[i]->indexB << " ";
	}
	cout << endl;

	cout << "phoListA: " ;
	for(i=0;i<phoGroupA.size();i++){
		cout << phoGroupA[i]->seqID << " ";
	}
	cout << endl;

	cout << "phoListB: ";
	for(i=0;i<phoGroupB.size();i++){
		cout << phoGroupB[i]->seqID << " ";
	}
	cout << endl;

	cout << "phoListC: ";
	for(i=0;i<phoGroupC.size();i++){
		cout << phoGroupC[i]->seqID << " ";
	}
	cout << endl;
}

NuTree::NuTree(NuGraph* graph){
	this->graph = graph;
	this->adjMtx = new bool[graph->seqLen*graph->seqLen];
	this->masked = new bool[graph->seqLen];
	for(int i=0;i<graph->seqLen;i++){
		this->masked[i] = graph->masked[i];
	}
	this->poolSize = 100000;
	this->sampFreqNode = 0.5;
	this->sampFreqEdge = 0.5;
	this->totalSamp = 0.0;
}

NuTree::~NuTree(){
	delete [] adjMtx;
	delete [] masked;
}

void NuTree::updateNodeInfo(){
	int i,j;
	for(i=0;i<graph->seqLen;i++){
		if(masked[i]) continue;
		graph->allNodes[i]->updateNodeInformation(this);
	}
}

void NuTree::updateNodeInfoCG(){
	int i,j;
	for(i=0;i<graph->seqLen;i++){
		if(masked[i]) continue;
		graph->allNodes[i]->updateNodeInformationCG(this);
	}
}

void NuTree::printNodeInfo(){
	cout << "node sampling freq: " << endl;
	for(int i=0;i<graph->seqLen;i++){
		cout << i << " " << graph->allNodes[i]->samplingFreq << endl;
		//cout << graph->allNodes[i]->connectToNeighbor << endl;
	}
}

void NuTree::updateEdgeInfo(){
	vector<NuEdge*> tmpGeList;
	int i,j;

	for(i=0;i<graph->seqLen;i++){
		if(masked[i]) continue;
		for(j=0;j<graph->seqLen;j++){
			if(masked[j]) continue;
			if(i==j) continue;
			graph->allEdges[i*graph->seqLen+j]->updateEdgeInfo(this);
		}
	}

	for(i=0;i<geList.size();i++){
		if(masked[geList[i]->indexA]) continue;
		if(masked[geList[i]->indexB]) continue;
		tmpGeList.push_back(geList[i]);
	}

	geList.clear();

	for(i=0;i<tmpGeList.size();i++){
		NuEdge* eA = tmpGeList[i];
		NuEdge* eB = graph->allEdges[eA->indexB*graph->seqLen + eA->indexA];
		eA->updateEdgeInfo(this);
		eB->updateEdgeInfo(this);

		if(eA->nodeListA.size() < eA->nodeListB.size()){
			geList.push_back(eB);
		}
		else {
			geList.push_back(eA);
		}
	}
}

void NuTree::updateEdgeInfoCG(){


	vector<NuEdge*> tmpGeList;
	int i,j;

	for(i=0;i<graph->seqLen;i++){
		for(j=0;j<graph->seqLen;j++){
			if(i==j) continue;
			graph->allEdges[i*graph->seqLen+j]->updateEdgeInfoCG(this);
		}
	}


	for(i=0;i<geList.size();i++){
		tmpGeList.push_back(geList[i]);
	}
	geList.clear();

	for(i=0;i<tmpGeList.size();i++){
		NuEdge* eA = tmpGeList[i];
		NuEdge* eB = graph->allEdges[eA->indexB*graph->seqLen + eA->indexA];
		eA->updateEdgeInfoCG(this);
		eB->updateEdgeInfoCG(this);

		if(eA->nodeListA.size() < eA->nodeListB.size()){
			geList.push_back(eB);
		}
		else {
			geList.push_back(eA);
		}
	}
}

void NuTree::updateSamplingInfo(){

	this->sampFreqNode = 0.0;
	this->sampFreqEdge = 0.0;
	int i,j,start, end;
	for(i=0;i<graph->seqLen;i++){
		this->sampFreqNode += graph->allNodes[i]->samplingFreq;
	}
	for(i=0;i<geList.size();i++){
		this->sampFreqEdge += geList[i]->samplingFreq;
	}

	double pAdd = 0;

	/*
	 * node rand pool
	 */
	for(i=0;i<poolSize;i++){
		this->randPoolNode[i] = 0;
	}

	for(i=0;i<graph->seqLen;i++){
		start = (int)(pAdd*poolSize);
		pAdd += graph->allNodes[i]->samplingFreq/sampFreqNode;
		end = (int)(pAdd*poolSize);
		for(j=start;j<end;j++){
			randPoolNode[j] = i;
		}
	}

	/*
	 * edge rand pool
	 */
	for(i=0;i<poolSize;i++){
		this->randPoolEdge[i] = 0;
	}

	pAdd = 0;
	for(i=0;i<geList.size();i++){
		start = (int)(pAdd*poolSize);
		int fq = geList[i]->samplingFreq;
		pAdd += 1.0*fq/sampFreqEdge;
		end = (int)(pAdd*poolSize);
	for(j=start;j<end;j++){
			randPoolEdge[j] = i;
		}
	}

	double sum = sampFreqNode + sampFreqEdge;
	this->sampFreqNode = (this->sampFreqNode * graph->et->para->kNodeFreq) / (this->sampFreqNode * graph->et->para->kNodeFreq + sampFreqEdge);
	this->sampFreqEdge = this->sampFreqEdge / (this->sampFreqNode * graph->et->para->kNodeFreq + sampFreqEdge);
	this->totalSamp = sum;
}

void NuTree::randomInit(){
	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamer* randRot;
	CsMove randMove;

	for(int i=0;i<graph->seqLen;i++){
		if(graph->allNodes[i]->samplingFreq == 0) continue;
		randNode = graph->allNodes[i];
		randRot = graph->rotLib->riboseRotLib->getRandomRotamer(randNode->baseType);
		randNode->updateRiboseRotamer(randRot);
		randNode->acceptRotMutation();
	}
	for(int i=0;i<geList.size();i++){
		randEdge = geList[i];
		if(randEdge->samplingFreq == 0) continue;

		randMove = randEdge->moveSet->getRandomMove();

		/*
		BaseDistanceMatrix dm0(randEdge->cm);
		BaseDistanceMatrix dm1(randMove);
		int pairtype = graph->pairLib->getPairType(dm1, randEdge->nodeA->baseType, randEdge->nodeB->baseType, randEdge->sep);
		cout << "pairType: " << pairtype << endl;

		double dist = dm0.distanceTo(dm1);
		printf("rand init edge %d-%d, distance = %5.3f\n", geList[i]->indexA, geList[i]->indexB, dist);
		*/

		randEdge->updateCsMove(randMove);
		randEdge->acceptMutation();
	}
}

void NuTree::randomInitCG(){
	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamerCG* randRot;
	CsMove randMove;

	for(int i=0;i<graph->seqLen;i++){
		if(graph->allNodes[i]->samplingFreq == 0) continue;
		randNode = graph->allNodes[i];
		randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
		randNode->updateRiboseRotamerCG(randRot);
		randNode->acceptRotMutationCG();
	}

	for(int i=0;i<geList.size();i++){
		if(randEdge->samplingFreq == 0) continue;
		randEdge = geList[i];

		BaseDistanceMatrix dm0(randEdge->cm);
		randMove = randEdge->moveSet->getRandomMove();
		BaseDistanceMatrix dm1(randMove);
		int pairtype = graph->pairLib->getPairType(dm1, randEdge->nodeA->baseType, randEdge->nodeB->baseType, randEdge->sep);

		//double dist = dm0.distanceTo(dm1);
		//printf("rand init edge %d-%d, distance = %5.3f\n", geList[i]->indexA, geList[i]->indexB, dist);

		randEdge->updateCsMoveCG(randMove);
		randEdge->acceptMutationCG();
	}
}

void NuTree::printEdges(){
	for(int i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		printf("%-3d %-3d %7.3f\n", e->indexA, e->indexB, e->weight);
		e->moveSet->printMoveSetInfo();
	}
}

void NuTree::printEdgeInfo(const string& output){
	ofstream out;
	out.open(output, ios::out);
	char xx[200];
	for(int i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		sprintf(xx, "%-3d %-3d %7.3f %8.3f\n", e->indexA, e->indexB, e->weight, e->samplingFreq);
		out << string(xx);
	}
	out.close();
}

void NuTree::printEdgeInfo(){
	
	for(int i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		printf("%-3d %-3d %7.3f %8.3f\n", e->indexA, e->indexB, e->weight, e->samplingFreq);

	}
}

graphInfo* NuTree::runAtomicMC(){

	bool debug = false;
	randomInit();

	double T0 = this->graph->et->para->T0;
	double T1 = this->graph->et->para->T1;
	double T2 = this->graph->et->para->T2;
	double T3 = this->graph->et->para->T3;

	int stepNum1 = (int)(this->totalSamp * graph->et->para->kStepNum1);
	int stepNum2 = (int)(this->totalSamp * graph->et->para->kStepNum2);
	int stepNum3 = (int)(this->totalSamp * graph->et->para->kStepNum3);
	
	double anneal = 0.95;

	double curEne = graph->totalEnergy();
	double lastEne = curEne;

	int i,j,k, randPos, nAc, eAc, nTot, eTot;
	double T, randP, mutE;

	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamer* randRot;
	CsMove randMove;

	int len = graph->seqLen;

	if(debug) {
		cout << "check init energy: " << endl;
		graph->checkEnergy();
	}

	for(T=T0;T>T1;T=T*anneal){

		nAc = 0;
		eAc = 0;
		nTot = 0;
		eTot = 0;

		for(k=0;k<stepNum1;k++){
			randP = rand()*1.0/RAND_MAX;
			if(debug) {
				cout << "randP: " << randP << endl;
			}
			if(randP < sampFreqNode){
				/*
				 * rotamer mut
				 */
				if(debug) {
					cout << "rot mut: " << endl;
				}
				nTot ++;
				randPos = randPoolNode[rand()%poolSize];
				if(debug) {
					cout << "rot mut: " << "pos " << randPos << endl;
				}
				randNode = graph->allNodes[randPos];
				randRot = graph->rotLib->riboseRotLib->getRandomRotamer(randNode->baseType);
				if(debug) {
					cout << "update ribose rotamer" << endl;
				}
				randNode->updateRiboseRotamer(randRot);

				if(debug) {
					cout << "rot mut energy" << endl;
				}
				mutE = randNode->rotMutEnergy();

				if(debug) {
					cout << "rot mut, pos: " << randPos << " mutE: " << mutE << endl;
					graph->checkEnergy();
				}

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randNode->acceptRotMutation();
					curEne += mutE;
					nAc++;
				}
				else {
					randNode->clearRotMutation();
				}

				if(debug) {
					cout << "rot mut, pos: " << randPos << " AC/RJ " << endl;
					graph->checkEnergy();
				}
			}
			else {
				/*
				 * edge mut
				 */

				if(debug){
					cout << "edge mut: " << endl;
				}
				eTot ++;
				randPos = randPoolEdge[rand()%poolSize];
				if(debug) {
					cout << "edge mut: pos: " << randPos << endl;
				}
				
				randEdge = geList[randPos];
				randMove = randEdge->moveSet->getRandomMove();
					
				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before update cs" << endl;
					graph->checkEnergy();
					cout << "move: " << endl;
					randMove.print();
				}

				randEdge->updateCsMove(randMove);

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before mutE" << endl;
					graph->checkEnergy();
				}

				mutE = randEdge->mutEnergy();

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " mutE: " << mutE << endl;
					graph->checkEnergy();
				}

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randEdge->acceptMutation();
					curEne += mutE;
					eAc++;
				}
				else {
					randEdge->clearMutation();
				}

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " AC/RJ" << endl;
					graph->checkEnergy();
				}
			}
		}

		double totEne = graph->totalEnergy();
		graphInfo* gi = graph->getGraphInfo();
		double rms = gi->rmsd(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
	}

	cout << "fixed subClusterID move" << endl;
	for(T=T1;T>T2;T=T*anneal){

		nAc = 0;
		eAc = 0;
		nTot = 0;
		eTot = 0;

		for(k=0;k<stepNum2;k++){
			randP = rand()*1.0/RAND_MAX;

			if(randP < sampFreqNode){
		
				nTot ++;
				randPos = randPoolNode[rand()%poolSize];
				randNode = graph->allNodes[randPos];
				randRot = graph->rotLib->riboseRotLib->getRandomRotamer(randNode->baseType);
				randNode->updateRiboseRotamer(randRot);
				mutE = randNode->rotMutEnergy();

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randNode->acceptRotMutation();
					curEne += mutE;
					nAc++;
				}
				else {
					randNode->clearRotMutation();
				}
			}
			else {

				eTot ++;
				randPos = randPoolEdge[rand()%poolSize];
				randEdge = geList[randPos];
				randMove = randEdge->moveSet->getRandomMoveWithFixedSubCluster(randEdge->cm);
				randEdge->updateCsMove(randMove);
				mutE = randEdge->mutEnergy();

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randEdge->acceptMutation();
					curEne += mutE;
					eAc++;
				}
				else {
					randEdge->clearMutation();
				}
			}
		}

		double totEne = graph->totalEnergy();
		graphInfo* gi = graph->getGraphInfo();
		double rms = gi->rmsd(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
	}

	cout << "fixed SP1000ID move" << endl;
	for(T=T2;T>T3;T=T*anneal){

		nAc = 0;
		eAc = 0;
		nTot = 0;
		eTot = 0;

		for(k=0;k<stepNum3;k++){
			randP = rand()*1.0/RAND_MAX;

			if(randP < sampFreqNode){

				nTot ++;
				randPos = randPoolNode[rand()%poolSize];
				randNode = graph->allNodes[randPos];
				randRot = graph->rotLib->riboseRotLib->getRandomRotamer(randNode->baseType);
				randNode->updateRiboseRotamer(randRot);
				mutE = randNode->rotMutEnergy();

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randNode->acceptRotMutation();
					curEne += mutE;
					nAc++;
				}
				else {
					randNode->clearRotMutation();
				}
			}
			else {
				eTot ++;
				randPos = randPoolEdge[rand()%poolSize];
				randEdge = geList[randPos];
				randMove = randEdge->moveSet->getRandomMoveWithFixedSP1000Index(randEdge->cm);
				randEdge->updateCsMove(randMove);

				mutE = randEdge->mutEnergy();

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randEdge->acceptMutation();
					curEne += mutE;
					eAc++;
				}
				else {
					randEdge->clearMutation();
				}
			}
		}

		double totEne = graph->totalEnergy();
		graphInfo* gi = graph->getGraphInfo();
		double rms = gi->rmsd(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
	}
	
	double totEne = graph->totalEnergy();
	graphInfo* gi = graph->getGraphInfo();
	double rms = gi->rmsd(this->graph->initInfo);
	gi->setRMS(rms);
		
	return gi;

}

void NuTree::runCoarseGrainedMC(const string& output){
	
	bool debug = false;
	randomInitCG();

	int stepNum = 100000;

	double T0 = 5.0;
	double T1 = 0.01;
	double anneal = 0.95;

	double curEne = graph->totalEnergyCG();
	double lastEne = curEne;

	int i,j,k, randPos, nAc, eAc, nTot, eTot;
	double T, randP, mutE;

	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamerCG* randRot;
	CsMove randMove;

	int len = graph->seqLen;

	int count = 0;

	for(T=T0;T>T1;T=T*anneal){
		nAc = 0;
		eAc = 0;
		nTot = 0;
		eTot = 0;

		for(k=0;k<stepNum;k++){
			count++;

			randP = rand()*1.0/RAND_MAX;
			if(randP < sampFreqNode){
				/*
				 * rotamer mut
				 */
				nTot ++;
				randPos = randPoolNode[rand()%poolSize];

				randNode = graph->allNodes[randPos];
				randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
				randNode->updateRiboseRotamerCG(randRot);
				mutE = randNode->rotMutEnergyCG();

				if(debug) {
					cout << "rot mut, pos: " << randPos << " mutE: " << mutE << endl;
					graph->checkEnergyCG();
					cout << "finish check" << endl;

					double mutE2 = graph->totalEnergyCGTmp() - graph->totalEnergyCG();
					cout << "mutE: " << mutE << " " << "mutE2: " << mutE2 << endl;
				}

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randNode->acceptRotMutationCG();
					curEne += mutE;
					nAc++;
					if(debug) {
						cout << "rot mut, pos: " << randPos << " AC " << endl;
						graph->checkEnergyCG();
						cout << "finish check" << endl;
					}
				}
				else {
					randNode->clearRotMutationCG();
					if(debug) {
						cout << "rot mut, pos: " << randPos << " RJ " << endl;
						graph->checkEnergyCG();
						
						cout << "finish check" << endl;
					}					
				}


			}
			else {
				/*
				 * edge mut
				 */
				eTot ++;
				randPos = randPoolEdge[rand()%poolSize];
				randEdge = geList[randPos];
				randMove = randEdge->moveSet->getRandomMove();

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before update cs" << endl;
					graph->checkEnergyCG();
					cout << "finish check" << endl;
				}

				randEdge->updateCsMoveCG(randMove);

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before mutE" << endl;
					graph->checkEnergyCG();
					cout << "finish check" << endl;
				}

				mutE = randEdge->mutEnergyCG();

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " mutE: " << mutE << endl;
					graph->checkEnergyCG();
					cout << "finish check" << endl;
					double mutE2 = graph->totalEnergyCGTmp() - graph->totalEnergyCG();
					cout << "mutE: " << mutE << " " << "mutE2: " << mutE2 << endl;

				}

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randEdge->acceptMutationCG();
					curEne += mutE;
					eAc++;
					if(debug) {
						cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " AC" << endl;
						graph->checkEnergyCG();
						cout << "finish check" << endl;
					}					
				}
				else {
					randEdge->clearMutationCG();
					if(debug) {
						cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " RJ" << endl;
						graph->checkEnergyCG();
						cout << "finish check" << endl;
					}
				}
			}
		}

		double totEne = graph->totalEnergyCG();
		graphInfo* gi = graph->getGraphInfoCG();
		double rms = gi->rmsdCG(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
	}

	cout << "total step num: " << count << endl;

	graphInfo* gi = graph->getGraphInfoCG();
	gi->setRMS(gi->rmsdCG(this->graph->initInfo));
	gi->printPDBCG(output);
	delete gi;
}

graphInfo::graphInfo(int seqLen, int* seq, bool* con, bool* fixed, NuNode** nodes, double ene, AtomLib* atLib, int mode){


	this->seqLen = seqLen;
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->fixed = new bool[seqLen];
	this->nodes = new NuNode*[seqLen];
	this->ene = ene;

	for(int i=0;i<seqLen;i++){
		this->seq[i] = seq[i];
		this->connectToDownstream[i] = con[i];
		this->fixed[i] = fixed[i];
	}

	for(int i=0;i<seqLen;i++){

		NuNode* n = nodes[i];

		NuNode* node;
		if(mode == 0){
			node = new NuNode(n->seqID, n->baseType,  n->baseConf->cs1, n->baseConf->rot, n->riboseConf->rot, atLib);
		}
		else {
			node = new NuNode(n->seqID, n->baseType,  n->baseConf->cs1, n->baseConfCG->rot, n->riboseConfCG->rot, atLib);
		}

		if(node->baseConfCG != NULL)
			node->baseConfCG->copyValueFrom(n->baseConfCG);
		
		if(node->riboseConfCG != NULL) {
			node->riboseConfCG->copyValueFrom(n->riboseConfCG);
		}

		if(node->phoConf != NULL) {
			node->phoConf->copyValueFrom(n->phoConf);
		}
	
		if(node->phoConfTmp != NULL)
			node->phoConfTmp->copyValueFrom(n->phoConfTmp);
	
		node->connectToNeighbor = connectToDownstream[i];
		this->nodes[i] = node;
	}

	this->atLib = atLib;
	this->rms = 0.0;
}

graphInfo::~graphInfo(){
	delete [] seq;
	delete [] connectToDownstream;
	delete [] fixed;
	for(int i=0;i<seqLen;i++)
		delete nodes[i];
	delete [] nodes;
}

double graphInfo::rmsd(graphInfo* other){
	vector<XYZ> tList1;
	vector<XYZ> tList2;

	int seqID = 0;
	for(unsigned int i=0;i<this->seqLen;i++) {
		if(fixed[i]) continue;
		vector<Atom*> aList = nodes[i]->toAtomList(this->atLib);
		for(int j=0;j<aList.size();j++)
			tList1.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	for(unsigned int i=0;i<this->seqLen;i++) {
		if(fixed[i]) continue;
		vector<Atom*> aList = other->nodes[i]->toAtomList(this->atLib);
		for(int j=0;j<aList.size();j++)
			tList2.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}
	return NSPgeometry::simpleRMSD(tList1, tList2);
}

double graphInfo::rmsdCG(graphInfo* other){
	vector<XYZ> tList1;
	vector<XYZ> tList2;

	int seqID = 0;
	for(unsigned int i=0;i<this->seqLen;i++) {
		if(fixed[i]) continue;
		vector<Atom*> aList = nodes[i]->toAtomListCG(this->atLib);

		for(int j=0;j<aList.size();j++)
			tList1.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	for(unsigned int i=0;i<this->seqLen;i++) {
		if(fixed[i]) continue;
		vector<Atom*> aList = other->nodes[i]->toAtomListCG(this->atLib);
		for(int j=0;j<aList.size();j++)
			tList2.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	return NSPgeometry::simpleRMSD(tList1, tList2);
}

void graphInfo::printPDB(const string& outputFile){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	int seqID = 0;
	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {
		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomListWithPho(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
		atomNum += aList.size();
	}

	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of << "ene: " << this->ene << endl;
	of << "rms: " << this->rms << endl;
	of.close();

	vector<RNABase*> baseList = rc.getBaseList();
	for(int i=0;i<baseList.size();i++){
		vector<Atom*>* aList = baseList[i]->getAtomList();
		for(int j=0;j<aList->size();j++)	{
			delete aList->at(j);
		}
	}
}

void graphInfo::printPDBCG(const string& outputFile){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	int seqID = 0;
	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {

		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomListCG(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
		atomNum += aList.size();
	}

	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of << "ene: " << this->ene << endl;
	of << "rms: " << this->rms << endl;
	of.close();

	vector<RNABase*> baseList = rc.getBaseList();
	for(int i=0;i<baseList.size();i++){
		vector<Atom*>* aList = baseList[i]->getAtomList();
		for(int j=0;j<aList->size();j++)	{
			delete aList->at(j);
		}
	}
}

void graphInfo::printAlignedPDB(graphInfo* alignTarget, const string& outputFile){
	vector<XYZ> tList1;
	vector<XYZ> tList2;

	int seqID = 0;
	for(unsigned int i=0;i<this->seqLen;i++) {
		vector<Atom*> aList = nodes[i]->toAtomList(this->atLib);
		for(int j=0;j<aList.size();j++)
			tList1.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	for(unsigned int i=0;i<this->seqLen;i++) {
		vector<Atom*> aList = alignTarget->nodes[i]->toAtomList(this->atLib);
		for(int j=0;j<aList.size();j++)
			tList2.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	XYZ Acog = getCOG(tList1);
	XYZ Bcog = getCOG(tList2);


	vector<XYZ> listA;
	vector<XYZ> listB;
	for(int i=0;i<tList1.size();i++){
		XYZ a = tList1[i] - Acog;
		XYZ b = tList2[i] - Bcog;
		listA.push_back(a);
		listB.push_back(b);
	}
	TransForm tf = buildRotation(listA, listB);


	RNAChain rc;
	string s = "AUGC";
	char ss[20];

	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {

		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomListWithPho(atLib);
		for(int j=0;j<aList.size();j++){
			Atom* a = aList[j];
			XYZ newCoord = tf.transform(a->coord);
			newCoord = newCoord + Bcog;
			a->setCoord(newCoord);
			base->addAtom(aList[j]);
		}
		rc.addBase(base);
		atomNum += aList.size();
	}



	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of << "ene: " << this->ene << endl;
	of << "rms: " << this->rms << endl;
	of.close();

	vector<RNABase*> baseList = rc.getBaseList();
	for(int i=0;i<baseList.size();i++){
		vector<Atom*>* aList = baseList[i]->getAtomList();
		for(int j=0;j<aList->size();j++)	{
			delete aList->at(j);
		}
	}
}

NuGraph::NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib, RnaEnergyTable* et){

	this->pairLib = pairLib;
	this->rotLib = rotLib;
	this->atLib = atLib;
	this->moveLib = moveLib;
	this->et = et;
	this->initInfo = NULL;

}

NuGraph::NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib){
	this->pairLib = pairLib;
	this->rotLib = rotLib;
	this->atLib = atLib;
	this->moveLib = NULL;
	this->et = NULL;
	this->initInfo = NULL;
	initForMST(inputFile);
}

NuGraph::~NuGraph() {

	delete [] seq;
	delete [] wcPairPosID;
	delete [] connectToDownstream;
	delete [] masked;
	delete [] fixed;
	delete [] sepTable;

	for(int i=0;i<seqLen;i++){
		delete allNodes[i];
	}

	for(int i=0;i<seqLen*seqLen;i++){
		delete allEdges[i];
	}

	delete [] allNodes;
	delete [] allEdges;

	for(int i=0;i<seqLen;i++){
		delete initRiboseRotList[i];
		delete initBaseRotList[i];
	}
	if(initInfo != NULL)
		delete initInfo;
}

void NuGraph::init(const string& task, const string& pdbFile, const string& baseSeq, const string& baseSec, const string& cst, const string& chainBreak){
	int i,j,k;

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 */



	/*
	vector<string> templates = input.getMultiValues("template");
	vector<string> templatesAlignA = input.getMultiValues("alignNat");
	vector<string> templatesAlignB = input.getMultiValues("alignTmp");
	vector<string> templatesType = input.getMultiValues("tempType");
	*/

	seqLen = baseSeq.length();

	cout << "sequence length: " <<  seqLen << endl;

	this->seq = new int[seqLen];
	this->wcPairPosID = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->sepTable = new int[seqLen*seqLen];
	this->allNodes = new NuNode*[seqLen];
	this->allEdges = new NuEdge*[seqLen*seqLen];
	this->masked = new bool[seqLen];
	this->fixed = new bool[seqLen];


	cout << "read chain break" << endl;

	/*
	 * read chain break points
	 */
	for(i=0;i<seqLen;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[seqLen-1] = false;

	vector<string> spt;
	splitString(chainBreak, " ", &spt);
	for(i=0;i<spt.size();i++){
		int pos = atoi(spt[i].c_str());
		if(pos >= baseSeq.length() || pos < 0){
			cout << "invalid chain break position: " << pos << endl;
			exit(0);
		}
		connectToDownstream[atoi(spt[i].c_str())] = false;
	}

	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			int ij = i*seqLen+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else sepTable[ij] = 2;
		}
	}

	cout << "init nodes" << endl;

	/*
	 * init NuNodes
	 */
	RNAPDB pdb(pdbFile, "pdbid");
	vector<RNABase*> baseList = pdb.getBaseList();

	if(baseList.size() != baseSeq.length()) {
		cout << "pdb size: " << baseList.size() << " not equal to sequence length: " << baseSeq.length() << endl;
		exit(0);
	}
	else if(baseSeq.length() != baseSec.length()) {
		cout << "seq length not equal to sec length" << endl;
		exit(0);
	}
	else if(cst.length() != baseSeq.length()){
		cout << "invalid cst " << endl;
		cout << cst << endl;
		exit(1);
	}

	for(i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		this->wcPairPosID[i] = -1;

		RiboseRotamer* rot = new RiboseRotamer();

		if(!baseList[i]->backboneComplete()){
			rot->copyValueFrom(rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt));
		}
		else {
			rot->copyValueFrom(rotLib->riboseRotLib->getNearestRotamer(baseList[i]));
		}

		BaseRotamer* baseRot = new BaseRotamer(seq[i], atLib);
		BaseRotamerCG* baseRotCG = new BaseRotamerCG(seq[i], atLib);
		RiboseRotamerCG* riboRotCG = new RiboseRotamerCG(rot);

		this->initBaseRotList.push_back(baseRot);
		this->initRiboseRotList.push_back(rot);
		this->initBaseRotCGList.push_back(baseRotCG);
		this->initRiboseRotCGList.push_back(riboRotCG);
	}


	for(i=0;i<seqLen;i++){
		LocalFrame cs1 = baseList[i]->getCoordSystem();
		this->allNodes[i] = new NuNode(i, baseList[i]->baseTypeInt, cs1, initBaseRotList[i], initBaseRotCGList[i], initRiboseRotList[i], initRiboseRotCGList[i], atLib);

		if(i==0 && connectToDownstream[i] == false){
			this->allNodes[i]->samplingFreq = 0;
		}
		else if(i==seqLen-1 && connectToDownstream[i-1] == false){
			this->allNodes[i]->samplingFreq = 0;
		}
		else if(i>0 && i<seqLen-1 && connectToDownstream[i-1] == false && connectToDownstream[i] == false){
			this->allNodes[i]->samplingFreq = 0;
		}
		cout << "pos: " << i << " sampFreq: " << this->allNodes[i]->samplingFreq << endl;
		this->allNodes[i]->connectToNeighbor = connectToDownstream[i];
		this->allNodes[i]->graph = this;
	}


	/*
	 * parse secondary structure information
	 */

	char ss[seqLen];
	for(int i=0;i<seqLen;i++){
		ss[i] = baseSec[i];
	}

	map<char,char> brackets;
	brackets[')'] = '(';
	brackets[']'] = '[';
	brackets['}'] = '{';
	brackets['>'] = '<';
	brackets['a'] = 'A';
	brackets['b'] = 'B';
	brackets['c'] = 'C';
	brackets['d'] = 'D';
	brackets['e'] = 'E';
	brackets['f'] = 'F';
	brackets['g'] = 'G';
	brackets['h'] = 'H';
	brackets['i'] = 'I';
	brackets['j'] = 'J';
	map<char,char>::iterator it;

	for(i=0;i<seqLen;i++) {
		char c = ss[i];
		it = brackets.find(c);
		if(it == brackets.end()) continue;
		int preIndex = -1;
		for(j=i-1;j>=0;j--) {
			if(ss[j] == it->second) {
				preIndex = j;
				break;
			}
		}
		if(preIndex < 0) {
			cout << "invalid ssSeq: " << baseSec << endl;
			exit(1);
		}
		ss[i] = '.';
		ss[preIndex] = '.';
		wcPairPosID[i] = preIndex;
		wcPairPosID[preIndex] = i;
	}

	cout << "init edges" << endl;
	/*
	 * init NuEdges
	 */
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			this->allEdges[i*seqLen+j] = new NuEdge(allNodes[i], allNodes[j], this);
			this->allEdges[i*seqLen+j]->graph = this;
			this->allEdges[i*seqLen+j]->weight = 0.0;
			if(task == "refinement"){
				this->allEdges[i*seqLen+j]->initNearNativeMoveSet();
			}
			else if(task == "analysis") {
				//this->allEdges[i*seqLen+j]->weight = xxx; to be modified
			}
		}
	}

	cout << "init constraints" << endl;
	/*
	 * init constaints
	 */

	string cstString = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	for(i=0;i<seqLen;i++){
		char c = cst[i];
		for(j=0;j<cstString.length();j++){
			if(c == cstString[j]) {
				for(k=i+1;k<seqLen;k++){
					char d = cst[k];
					if(c == d){
						this->allEdges[i*seqLen+k]->weight = -999.9;
						this->allEdges[k*seqLen+i]->weight = -999.9;
						this->allEdges[i*seqLen+k]->fixNaiveMove();
						this->allEdges[k*seqLen+i]->fixNaiveMove();
					}
				}
			}
		}
	}

	cout << "init masked residue" << endl;
	/*
	 * init masked residues
	 */
	for(i=0;i<seqLen;i++) {
		char c = cst[i];
		if(c == '-') {
			this->masked[i] = true;
			this->fixed[i] = true;
			this->allNodes[i]->samplingFreq = 0.0;
			if(i > 0)
				this->connectToDownstream[i-1] = false;
		}
		else if(c == 'F') {
			this->masked[i] = false;
			this->fixed[i] = true;
		}
		else {
			this->masked[i] = false;
			this->fixed[i] = false;
		}
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			this->geList.push_back(allEdges[i*seqLen+j]);
		}
	}
	cout << "finish init" << endl;
}

void NuGraph::initPho(){
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i]){
			et->pb->buildPhosphate(allNodes[i]->riboseConf, allNodes[i]->riboseConf, allNodes[i]->phoConf);
			allNodes[i]->phoConfTmp->copyValueFrom(allNodes[i]->phoConf);
		}
	}

	
}

void NuGraph::initPho(PO3Builder* pb) {
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i]){
			pb->buildPhosphate(allNodes[i]->riboseConf, allNodes[i]->riboseConf, allNodes[i]->phoConf);
			allNodes[i]->phoConfTmp->copyValueFrom(allNodes[i]->phoConf);
		}
	}
}


void NuGraph::initForMC(const string& inputFile){

	NSPtools::InputParser input(inputFile);
	input.printOptions();

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 */

	string task = input.getValue("task");
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string cst = input.getValue("cst");
	string chainBreak = input.getValue("break");

	init(task, pdbFile, baseSeq, baseSec, cst, chainBreak);
	
	initPho();
	this->initInfo = new graphInfo(seqLen, seq, connectToDownstream, masked, allNodes, 0.0, atLib, 0);	
}

void NuGraph::initForCGMC(const string& inputFile){

	NSPtools::InputParser input(inputFile);
	input.printOptions();

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 */

	string task = input.getValue("task");
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string cst = input.getValue("cst");
	string chainBreak = input.getValue("break");

	init(task, pdbFile, baseSeq, baseSec, cst, chainBreak);
	this->initInfo = new graphInfo(seqLen, seq, connectToDownstream, masked, allNodes, 0.0, atLib, 1);	
}

void NuGraph::initForMST(const string& inputFile){
	NSPtools::InputParser input(inputFile);
	input.printOptions();

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 *      
	 */

	string task = input.getValue("task");
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string cst = input.getValue("cst");
	string chainBreak = input.getValue("break");

	init(task, pdbFile, baseSeq, baseSec, cst, chainBreak);
}

void NuGraph::initForSingleResiduePrediction(const string& inputFile, int pos){
	NSPtools::InputParser input(inputFile);
	input.printOptions();

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 */

	string task = input.getValue("task");
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");

	string cst = input.getValue("cst");
	string chainBreak = input.getValue("break");

	init(task, pdbFile, baseSeq, baseSec, cst, chainBreak);
	initPho();
	this->initInfo = new graphInfo(seqLen, seq, connectToDownstream, masked, allNodes, 0.0, atLib, 0);	

	set<int> maskedPositions;
	NuNode* nodeTarget = allNodes[pos];
	vector<XYZ> coordTarget;
	for(int i=0;i<nodeTarget->baseConf->rot->atomNum;i++){
		coordTarget.push_back(nodeTarget->baseConf->coords[i]);
	}
	for(int i=0;i<nodeTarget->riboseConf->rot->atomNum;i++){
		coordTarget.push_back(nodeTarget->riboseConf->coords[i]);
	}
	for(int i=0;i<4;i++){
		coordTarget.push_back(nodeTarget->phoConf->coords[i]);
	}

	double d;
	int i,j,k;
	for(i=0;i<seqLen;i++){
		if(i==pos) continue;
		NuNode* nodeA = allNodes[i];
		vector<XYZ> coordA;
		for(j=0;j<nodeA->baseConf->rot->atomNum;j++){
			coordA.push_back(nodeA->baseConf->coords[j]);
		}
		for(j=0;j<nodeA->riboseConf->rot->atomNum;j++){
			coordA.push_back(nodeA->riboseConf->coords[j]);
		}
		for(j=0;j<4;j++){
			coordA.push_back(nodeA->phoConf->coords[j]);
		}
		double minD = 99.9;
		for(int j=0;j<coordTarget.size();j++){
			for(int k=0;k<coordA.size();k++){
				d = coordTarget[j].distance(coordA[k]);
				if(d < minD){
					minD = d;
				}
			}
		}
		if(minD > 10.0) {
			maskedPositions.insert(i);
		}
	}

	int seqLen = baseSeq.length();
	char xx[seqLen+1];
	for(i=0;i<seqLen;i++)
		xx[i] = '0';
	xx[seqLen] = '\0';
	for(i=0;i<seqLen;i++){
		if(maskedPositions.contains(i))
			xx[i] = '-';
		else if(i==pos)
			xx[i] = '0';
		else 
			xx[i] = 'F';
	}
	cst = string(xx);

	cout << "single residue prediction cst: " << endl;
	cout << cst << endl;

	string cstString = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	for(i=0;i<seqLen;i++){
		char c = cst[i];
		for(j=0;j<cstString.length();j++){
			if(c == cstString[j]) {
				for(k=i+1;k<seqLen;k++){
					char d = cst[k];
					if(c == d){
						this->allEdges[i*seqLen+k]->weight = -999.9;
						this->allEdges[k*seqLen+i]->weight = -999.9;
						this->allEdges[i*seqLen+k]->fixNaiveMove();
						this->allEdges[k*seqLen+i]->fixNaiveMove();
					}
				}
			}
		}
	}


	for(i=0;i<seqLen;i++) {
		char c = cst[i];
		if(c == '-') {
			this->masked[i] = true;
			this->fixed[i] = true;
			this->allNodes[i]->samplingFreq = 0.0;
			if(i > 0)
				this->connectToDownstream[i-1] = false;
		}
		else if(c >= 'A' && c <= 'Z') {
			this->masked[i] = false;
			this->fixed[i] = true;
			if(this->allNodes[i]->samplingFreq != 0)
				this->allNodes[i]->samplingFreq = 0.5;
		}
		else {
			this->masked[i] = false;
			this->fixed[i] = false;
		}
	}



}

bool MST_cmp_weight(NuEdge* e1, NuEdge* e2){
	return (e1->weightRand < e2->weightRand);
}

int MST_find(int* parent, int f){
	while(parent[f] >= 0){
		f = parent[f];
	}
	return f;
}

void NuGraph::initRandWeight(){

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			allEdges[i*seqLen+j]->weightRand = allEdges[i*seqLen+j]->weight - rand()*1.0/RAND_MAX;
		}
	}
}

void NuGraph::MST_kruskal(NuTree* outputTree){
	int i,j,n,m,a,b;
	int parent[seqLen];
	for(i=0;i<seqLen;i++){
		parent[i] = -1;
	}

	int ne = geList.size();

	NuEdge* sortedGeList[ne];

	for(i=0;i<ne;i++){
		sortedGeList[i] = geList[i];
	}

	sort(sortedGeList, sortedGeList+ne, MST_cmp_weight);

	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			outputTree->adjMtx[i*seqLen+j] = false;
		}
	}
	outputTree->geList.clear();

	for(i=0;i<ne;i++){
		n = MST_find(parent, sortedGeList[i]->indexA);
		m = MST_find(parent, sortedGeList[i]->indexB);

		if(n != m) {
			a = sortedGeList[i]->indexA;
			b = sortedGeList[i]->indexB;
			parent[n] = m;
			outputTree->adjMtx[a*seqLen+b] = true;
			outputTree->adjMtx[b*seqLen+a] = true;
			if(a<b)
				outputTree->geList.push_back(allEdges[a*seqLen+b]);
			else
				outputTree->geList.push_back(allEdges[b*seqLen+a]);
		}
	}
}

void NuGraph::printAllEdge(){
	for(int i=0;i<seqLen;i++){
		for(int j=i+1;j<seqLen;j++){
			NuEdge* e = allEdges[i*seqLen+j];
			printf("%-3d %-3d %7.3f\n", i, j, e->weight);
		}
	}
}

void NuGraph::checkEnergy(){
	for(int i=0;i<seqLen;i++){
		if(masked[i]) continue;
		allNodes[i]->checkEnergy();
	}
	for(int i=0;i<seqLen;i++){
		if(masked[i]) continue;
		for(int j=0;j<seqLen;j++){
			if(masked[j]) continue;
			if(i==j) continue;
			allEdges[i*seqLen+j]->checkEnergy();
			allEdges[i*seqLen+j]->checkReversePair();
		}
	}
}

void NuGraph::checkEnergyCG(){
	for(int i=0;i<seqLen;i++){
		if(masked[i]) continue;
		allNodes[i]->checkEnergyCG();
	}
	for(int i=0;i<seqLen;i++){
		if(masked[i]) continue;
		for(int j=0;j<seqLen;j++){
			if(masked[j]) continue;
			if(i==j) continue;
			allEdges[i*seqLen+j]->checkEnergyCG();
			allEdges[i*seqLen+j]->checkReversePairCG();
		}
	}
}

double NuGraph::totalEnergy(){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){

		if(masked[i]) continue;
		ene += allNodes[i]->riboseConf->rot->energy;
		ene += allNodes[i]->phoConf->ene;
	}

	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		for(j=i+1;j<seqLen;j++){
			if(masked[j]) continue;
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergy(allNodes[i]->baseConf, allNodes[j]->baseConf, sep, et);
			ene += nuBaseRiboseEnergy(allNodes[i]->baseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuBaseRiboseEnergy(allNodes[j]->baseConf, allNodes[i]->riboseConf, sepR, et);
			ene += nuBasePhoEnergy(allNodes[i]->baseConf, allNodes[j]->phoConf, sep, et);
			ene += nuBasePhoEnergy(allNodes[j]->baseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuRiboseRiboseEnergy(allNodes[i]->riboseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuRibosePhoEnergy(allNodes[i]->riboseConf, allNodes[j]->phoConf, sep, et);
			ene += nuRibosePhoEnergy(allNodes[j]->riboseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuPhoPhoEnergy(allNodes[i]->phoConf, allNodes[j]->phoConf, sep, et);
		}
	}
	return ene;
}

double NuGraph::totalEnergyCG(){
		double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		ene += allNodes[i]->riboseConfCG->rot->energy;
		ene += allNodes[i]->bbcg;
	}

	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		for(j=i+1;j<seqLen;j++){
			if(masked[j]) continue;
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergyCG(allNodes[i]->baseConfCG, allNodes[j]->baseConfCG, sep, et);
			ene += nuBaseRiboseEnergyCG(allNodes[i]->baseConfCG, allNodes[j]->riboseConfCG, sep, et);
			ene += nuBaseRiboseEnergyCG(allNodes[j]->baseConfCG, allNodes[i]->riboseConfCG, sepR, et);
			ene += nuRiboseRiboseEnergyCG(allNodes[i]->riboseConfCG, allNodes[j]->riboseConfCG, sep, et);
		}
	}
	return ene;
}

double NuGraph::totalEnergyCGTmp(){
		double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		ene += allNodes[i]->riboseConfCGTmp->rot->energy;
		ene += allNodes[i]->bbcgTmp;
	}

	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		for(j=i+1;j<seqLen;j++){
			if(masked[j]) continue;
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergyCG(allNodes[i]->baseConfCGTmp, allNodes[j]->baseConfCGTmp, sep, et);
			ene += nuBaseRiboseEnergyCG(allNodes[i]->baseConfCGTmp, allNodes[j]->riboseConfCGTmp, sep, et);
			ene += nuBaseRiboseEnergyCG(allNodes[j]->baseConfCGTmp, allNodes[i]->riboseConfCGTmp, sepR, et);
			ene += nuRiboseRiboseEnergyCG(allNodes[i]->riboseConfCGTmp, allNodes[j]->riboseConfCGTmp, sep, et);
		}
	}
	return ene;
}

double NuGraph::totalEnergyTmp(){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		ene += allNodes[i]->riboseConfTmp->rot->energy;
		ene += allNodes[i]->phoConfTmp->ene;
	}

	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		for(j=i+1;j<seqLen;j++){
			if(masked[j]) continue;
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergy(allNodes[i]->baseConfTmp, allNodes[j]->baseConfTmp, sep, et);
			ene += nuBaseRiboseEnergy(allNodes[i]->baseConfTmp, allNodes[j]->riboseConfTmp, sep, et);
			ene += nuBaseRiboseEnergy(allNodes[j]->baseConfTmp, allNodes[i]->riboseConfTmp, sepR, et);
			ene += nuBasePhoEnergy(allNodes[i]->baseConfTmp, allNodes[j]->phoConfTmp, sep, et);
			ene += nuBasePhoEnergy(allNodes[j]->baseConfTmp, allNodes[i]->phoConfTmp, sepR, et);
			ene += nuRiboseRiboseEnergy(allNodes[i]->riboseConfTmp, allNodes[j]->riboseConfTmp, sep, et);
			ene += nuRibosePhoEnergy(allNodes[i]->riboseConfTmp, allNodes[j]->phoConfTmp, sep, et);
			ene += nuRibosePhoEnergy(allNodes[j]->riboseConfTmp, allNodes[i]->phoConfTmp, sepR, et);
			ene += nuPhoPhoEnergy(allNodes[i]->phoConfTmp, allNodes[j]->phoConfTmp, sep, et);
		}
	}
	return ene;
}

double NuGraph::totalEnergy2(){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		ene += allNodes[i]->ene;
	}

	for(i=0;i<seqLen;i++){
		if(masked[i]) continue;
		for(j=i+1;j<seqLen;j++){
			if(masked[j]) continue;
			for(k=0;k<9;k++)
				ene += allEdges[i*seqLen+j]->pairEne[k];
		}
	}
	return ene;
}

void NuGraph::printEnergy(){
	for(int i=0;i<seqLen;i++){
		printf("node: %2d ene: %7.3f eneTmp: %7.3f\n", i, allNodes[i]->ene, allNodes[i]->eneTmp);
	}
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			NuEdge* eg = this->allEdges[i*seqLen+j];
			printf("edge: %2d %2d sep: %2d BB %7.3f BR %7.3f BP %7.3f RB %7.3f RR %7.3f RP %7.3f PB %7.3f PR %7.3f PP %7.3f\n", i, j, eg->sep, eg->pairEne[0], eg->pairEne[1], eg->pairEne[2],eg->pairEne[3],eg->pairEne[4],eg->pairEne[5],eg->pairEne[6],eg->pairEne[7],eg->pairEne[8]);
			printf("eTmp: %2d %2d sep: %2d BB %7.3f BR %7.3f BP %7.3f RB %7.3f RR %7.3f RP %7.3f PB %7.3f PR %7.3f PP %7.3f\n", i, j, eg->sep, eg->pairEneTmp[0], eg->pairEneTmp[1], eg->pairEneTmp[2],eg->pairEneTmp[3],eg->pairEneTmp[4],eg->pairEneTmp[5],eg->pairEneTmp[6],eg->pairEneTmp[7],eg->pairEneTmp[8]);
		}
	}
}

void NuGraph::printEnergyCG(){
	for(int i=0;i<seqLen;i++){
		printf("node: %2d ene: %7.3f eneTmp: %7.3f\n", i, allNodes[i]->eneCG, allNodes[i]->eneCGTmp);
	}
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			NuEdge* eg = this->allEdges[i*seqLen+j];
			printf("edge: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", i, j, eg->sep, eg->pairEneCG[0], eg->pairEneCG[1], eg->pairEneCG[2],eg->pairEneCG[3]);
			printf("eTmp: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", i, j, eg->sep, eg->pairEneCGTmp[0], eg->pairEneCGTmp[1], eg->pairEneCGTmp[2],eg->pairEneCGTmp[3]);
		}
	}
}

graphInfo* NuGraph::getGraphInfo(){
	graphInfo* gi = new graphInfo(seqLen, seq, connectToDownstream, masked, allNodes, totalEnergy(), atLib, 0);
	return gi;
}

graphInfo* NuGraph::getGraphInfoCG(){
	graphInfo* gi = new graphInfo(seqLen, seq, connectToDownstream, masked, allNodes, totalEnergyCG(), atLib, 1);
	return gi;
}

} /* namespace NSPpredNA */
