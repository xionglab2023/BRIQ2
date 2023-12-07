/*
 * NuGraph.cpp
 *
 *  Created on: 2023Äê11ÔÂ15ÈÕ
 *      Author: nuc
 */

#include <predNA/NuGraph.h>

namespace NSPpredNA {

NuNode::NuNode(int id, int baseType, LocalFrame& cs1, RiboseRotamer* riboRot, AtomLib* atLib){


	this->seqID = id;
	this->baseType = baseType;
	this->connectToNeighbor = false;

	BaseRotamer* baseRot = new BaseRotamer(baseType, atLib);
	BaseRotamer* baseRotTmp = new BaseRotamer(baseType, atLib);

	this->baseConf = new BaseConformer(baseRot, cs1);
	this->baseConfTmp = new BaseConformer(baseRotTmp, cs1);

	this->riboseConf = new RiboseConformer(riboRot, cs1);
	this->riboseConfTmp = new RiboseConformer(riboRot, cs1);

	this->phoConf = new PhosphateConformer();
	this->phoConfTmp = new PhosphateConformer();

	BaseRotamerCG* baseRotCG = new BaseRotamerCG(baseType, atLib);
	BaseRotamerCG* baseRotCGTmp = new BaseRotamerCG(baseType, atLib);

	this->baseConfCG = new BaseConformerCG(baseRotCG, cs1);
	this->baseConfCGTmp = new BaseConformerCG(baseRotCGTmp, cs1);

	RiboseRotamerCG* riboRotCG = new RiboseRotamerCG(riboRot);
	RiboseRotamerCG* riboRotCGTmp = new RiboseRotamerCG(riboRot);

	this->riboseConfCG = new RiboseConformerCG(riboRotCG, cs1);
	this->riboseConfCGTmp = new RiboseConformerCG(riboRotCGTmp, cs1);

	this->graph = NULL;
	this->ene = 0;
	this->eneTmp = 0;
	this->eneCG = 0;
	this->eneCGTmp = 0;

	this->samplingFreq = 1.0;
}

NuNode::~NuNode(){

	delete this->baseConf->rot;
	delete this->baseConfTmp->rot;

	delete this->baseConf;
	delete this->baseConfTmp;
	delete this->riboseConf;
	delete this->riboseConfTmp;
	delete this->phoConf;
	delete this->phoConfTmp;

	delete this->baseConfCG->rot;
	delete this->baseConfCGTmp->rot;
	delete this->riboseConfCG->rot;
	delete this->riboseConfCGTmp->rot;
	delete this->baseConfCG;
	delete this->baseConfCGTmp;
	delete this->riboseConfCG;
	delete this->riboseConfCGTmp;
}

void NuNode::updateNodeInformation(NuTree* tree){

	this->graph = tree->graph;

	int i,j;
	this->neighborList.clear();
	this->connectionBreakPoints.clear();
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

	this->eneCG = riboseConfCG->rot->energy;
	if(connectToNeighbor)
		this->eneCG += nuConnectionEnergyCG(graph->allNodes[connectionBreakPoints[seqID]]->riboseConfCG, graph->allNodes[connectionBreakPoints[seqID]+1]->riboseConfCG, graph->et);


	this->eneTmp = this->ene;
	this->eneCGTmp = this->eneCG;
}

void NuNode::printNodeInfo(){
	cout << "nodeID: " << seqID << endl;
	cout << "neighbors: ";
	for(int i=0;i<neighborList.size();i++){
		cout << " " << neighborList[i];
	}

	cout << endl;
	cout << "connection breaks: ";
	for(int i=0;i<connectionBreakPoints.size();i++){
		cout << " " << connectionBreakPoints[i];
	}
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

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		sep = graph->sepTable[i*len+seqID];

		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->pairEneTmp[1] = nuBaseRiboseEnergy(graph->allNodes[i]->baseConfTmp, this->riboseConfTmp, sep, graph->et);
		egB->pairEneTmp[3] = egA->pairEneTmp[1];
	}


	/*
	 * base-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){
			sep = graph->sepTable[i*len+phoGroupC[j]->seqID];
			if(sep == 0) continue;
			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + i];

			egA->pairEneTmp[2] = nuBasePhoEnergy(graph->allNodes[i]->baseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et);
			egB->pairEneTmp[6] = egA->pairEneTmp[2];
		}
	}

	/*
	 * ribose-ribose
	 */

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		sep = graph->sepTable[i*len+seqID];

		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->pairEneTmp[4] = nuRiboseRiboseEnergy(graph->allNodes[i]->riboseConfTmp, this->riboseConfTmp, sep, graph->et);
		egB->pairEneTmp[4] = egA->pairEneTmp[4];
	}

	/*
	 * ribose-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){
			sep = graph->sepTable[i*len+phoGroupC[j]->seqID];
			if(sep == 0) continue;
			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + i];

			egA->pairEneTmp[5] = nuRibosePhoEnergy(graph->allNodes[i]->riboseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et);
			egB->pairEneTmp[7] = egA->pairEneTmp[5];
		}
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
	}


	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->ene = graph->allNodes[j]->eneTmp;
	}
	this->ene = this->eneTmp;

	/*
	 * base-ribose energy
	 */

	for(i=0;i<len;i++){
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->pairEne[1] = egA->pairEneTmp[1];
		egB->pairEne[3] = egB->pairEneTmp[3];
	}


	/*
	 * base-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + i];

			egA->pairEne[2] = egA->pairEneTmp[2];
			egB->pairEne[6] = egB->pairEneTmp[6];
		}
	}

	/*
	 * ribose-ribose
	 */

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->pairEne[4] = egA->pairEneTmp[4];
		egB->pairEne[4] = egB->pairEneTmp[4];
	}

	/*
	 * ribose-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){

			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + i];

			egA->pairEne[5] = egA->pairEneTmp[5];
			egB->pairEne[7] = egB->pairEneTmp[7];

		}
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

	for(i=0;i<len;i++){
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->pairEneTmp[1] = egA->pairEne[1];
		egB->pairEneTmp[3] = egB->pairEne[3];
	}


	/*
	 * base-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + i];

			egA->pairEneTmp[2] = egA->pairEne[2];
			egB->pairEneTmp[6] = egB->pairEne[6];
		}
	}

	/*
	 * ribose-ribose
	 */

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->pairEneTmp[4] = egA->pairEne[4];
		egB->pairEneTmp[4] = egB->pairEne[4];
	}

	/*
	 * ribose-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){

			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len + i];

			egA->pairEneTmp[5] = egA->pairEne[5];
			egB->pairEneTmp[7] = egB->pairEne[7];

		}
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

	for(i=0;i<len;i++){
		egA = graph->allEdges[i*len + seqID];
		mutE += egA->pairEneTmp[1] - egA->pairEne[1];
	}


	/*
	 * base-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){
			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];

			mutE += egA->pairEneTmp[2] - egA->pairEne[2];
		}
	}

	/*
	 * ribose-ribose
	 */

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		egA = graph->allEdges[i*len + seqID];
		mutE += egA->pairEneTmp[4] - egA->pairEne[4];
	}

	/*
	 * ribose-pho
	 */
	for(i=0;i<len;i++){
		for(j=0;j<phoGroupC.size();j++){

			egA = graph->allEdges[i*len + phoGroupC[j]->seqID];
			mutE += egA->pairEneTmp[5] - egA->pairEne[5];

		}
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

bool NuNode::checkEnergyCG(){
	double e3 = this->riboseConfCG->rot->energy;
	bool tag = true;


	if(connectToNeighbor)
		e3 += nuConnectionEnergyCG(this->riboseConfCG, graph->allNodes[seqID+1]->riboseConfCG, graph->et);
	double e4 = this->riboseConfCGTmp->rot->energy;
	if(connectToNeighbor)
		e4 += nuConnectionEnergyCG(this->riboseConfCGTmp, graph->allNodes[seqID+1]->riboseConfCGTmp, graph->et);

	if(abs(e3 - eneCG) > 0.001){
		printf("ERROR ene %7.3f %7.3f\n",e3, eneCG);
		tag = false;
	}

	if(abs(e4 - eneTmp) > 0.001){
		printf("ERROR ene %7.3f %7.3f\n",e4, eneTmp);
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

	this->eneCGTmp = rot->energy;

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		this->eneCGTmp += nuConnectionEnergyCG(graph->allNodes[connectionBreakPoints[i]]->riboseConfCG, graph->allNodes[connectionBreakPoints[i]+1]->riboseConfCG, graph->et);
	}

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		sep = graph->sepTable[i*len+seqID];
		sepR = graph->sepTable[seqID*len+i];
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->eneCGTmp[1] = nuBaseRiboseEnergyCG(graph->allNodes[i]->baseConfCGTmp, this->riboseConfCGTmp, sep, graph->et);
		egB->eneCGTmp[2] = egA->eneCGTmp[1];
		egA->eneCGTmp[3] = nuRiboseRiboseEnergyCG(graph->allNodes[i]->riboseConfCGTmp, this->riboseConfCGTmp, sep, graph->et);
		egB->eneCGTmp[3] = egA->eneCGTmp[3];
	}

}

void NuNode::acceptRotMutationCG(){
	this->riboseConfCG->copyValueFrom(this->riboseConfCGTmp);
	this->eneCG = this->eneCGTmp;

	int len = graph->seqLen;
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		sep = graph->sepTable[i*len+seqID];
		sepR = graph->sepTable[seqID*len+i];
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->eneCG[1] = egA->eneCGTmp[1];
		egB->eneCG[2] = egB->eneCGTmp[2];
		egA->eneCG[3] = egA->eneCGTmp[3];
		egB->eneCG[3] = egB->eneCGTmp[3];
	}
}

void NuNode::clearRotMutationCG(){
	this->riboseConfCGTmp->copyValueFrom(this->riboseConfCG);
	this->eneCGTmp = this->eneCG;
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;
	int len = graph->seqLen;

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		sep = graph->sepTable[i*len+seqID];
		sepR = graph->sepTable[seqID*len+i];
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		egA->eneCGTmp[1] = egA->eneCG[1];
		egB->eneCGTmp[2] = egB->eneCG[2];
		egA->eneCGTmp[3] = egA->eneCG[3];
		egB->eneCGTmp[3] = egB->eneCG[3];
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

double NuNode::mutEnergyCG(){
	double mutE = this->eneCGTmp - this->eneCG;
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;
	int len = graph->seqLen;

	for(i=0;i<len;i++){
		if(i == seqID) continue;
		sep = graph->sepTable[i*len+seqID];
		sepR = graph->sepTable[seqID*len+i];
		egA = graph->allEdges[i*len + seqID];
		egB = graph->allEdges[seqID*len + i];

		mutE += egA->eneCGTmp[1] - egA->eneCG[1];
		mutE += egA->eneCGTmp[3] = egA->eneCG[3];
	}
	return mutE;
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
		this->eneCG[i] = 0.0;
		this->eneCGTmp[i] = 0.0;
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
		this->eneCG[i] = 0.0;
		this->eneCGTmp[i] = 0.0;
	}

	this->samplingFreq = 1.0;
}

NuEdge::~NuEdge(){
	delete this->ei;
	delete this->moveSet;
}

void NuEdge::initNativeMoveSet(){

	BaseDistanceMatrix dm(nodeA->baseConf->cs1, nodeB->baseConf->cs1);
	int clusterID = pairLib->getPairType(dm, nodeA->baseType, nodeB->baseType, sep);
	this->ei->setUniqueCluster(clusterID, pairLib);
	this->moveSet->updateEdgeInformation(ei);
	this->weight = this->ei->weight;
	this->weightRand = this->weight;
}

void NuEdge::fixNaiveMove(){
	CsMove cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->ei->setFixed(cm);
	this->moveSet->fixNativeMove(cm);
	this->weight = this->ei->weight;
	this->weightRand = this->weight;
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

	pairEne[0] = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, graph->et);
	pairEne[1] = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, graph->et);
	pairEne[2] = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, graph->et);

	pairEne[3] = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, graph->et);
	pairEne[4] = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, graph->et);
	pairEne[5] = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, graph->et);
	pairEne[6] = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, graph->et);
	pairEne[7] = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, graph->et);
	pairEne[8] = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, graph->et);


	eneCG[0] = nuBaseBaseEnergyCG(nodeA->baseConfCG, nodeB->baseConfCG, sep, graph->et);
	eneCG[1] = nuBaseRiboseEnergyCG(nodeA->baseConfCG, nodeB->riboseConfCG, sep, graph->et);
	eneCG[2] = nuBaseRiboseEnergyCG(nodeB->baseConfCG, nodeA->riboseConfCG, sepR, graph->et);
	eneCG[3] = nuRiboseRiboseEnergyCG(nodeA->riboseConfCG, nodeB->riboseConfCG, sep, graph->et);

	for(i=0;i<9;i++){
		pairEneTmp[i] = pairEne[i];
	}
	for(i=0;i<4;i++){
		eneCGTmp[i] = eneCG[i];
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

			eA->pairEneTmp[0] = nuBaseBaseEnergy(nA->baseConfTmp, nB->baseConfTmp, sep, graph->et);
			eA->pairEneTmp[1] = nuBaseRiboseEnergy(nA->baseConfTmp, nB->riboseConfTmp, sep, graph->et);
			eA->pairEneTmp[3] = nuBaseRiboseEnergy(nB->baseConfTmp, nA->riboseConfTmp, sepR, graph->et);
			eA->pairEneTmp[4] = nuRiboseRiboseEnergy(nA->riboseConfTmp, nB->riboseConfTmp, sep, graph->et);

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
			eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
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
			eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
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
			eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
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


			eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et);
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
	this->cmTmp = cm;
	LocalFrame cs1 = nodeA->baseConfCGTmp->cs1 + cm;
	nodeB->updateCoordinateCG(cs1);
	for(int i=0;i<edgeListB.size();i++){
		cs1 = edgeListB[i]->nodeA->baseConfCGTmp->cs1 + edgeListB[i]->cmTmp;
		edgeListB[i]->nodeB->updateCoordinateCG(cs1);
	}
}

double NuEdge::mutEnergyCG(){
	return 0.0;
}

void NuEdge::acceptMutationCG(){
	this->cm = this->cmTmp;
	for(int i=0;i<nodeListB.size();i++){
		nodeListB[i]->acceptCoordMoveCG();
	}
}

void NuEdge::clearMutationCG(){
	this->cmTmp = this->cm;
	for(int i=0;i<nodeListB.size();i++){
		nodeListB[i]->clearCoordMoveCG();
	}
}

bool NuEdge::checkEnergy(){

	bool tag = true;
	int sepR = -sep;
	if(abs(sepR) > 1)
		sepR = 2;

	double e0 = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, graph->et);
	double e1 = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, graph->et);
	double e2 = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, graph->et);
	double e3 = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, graph->et);
	double e4 = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, graph->et);
	double e5 = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, graph->et);
	double e6 = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, graph->et);
	double e7 = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, graph->et);
	double e8 = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, graph->et);

	double e0Tmp = nuBaseBaseEnergy(nodeA->baseConfTmp, nodeB->baseConfTmp, sep, graph->et);
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

bool NuEdge::checkEnergyCG(){

	bool tag = true;
	int sepR = -sep;
	if(abs(sepR) > 1)
		sepR = 2;

	double e0 = nuBaseBaseEnergyCG(nodeA->baseConfCG, nodeB->baseConfCG, sep, graph->et);
	double e1 = nuBaseRiboseEnergyCG(nodeA->baseConfCG, nodeB->riboseConfCG, sep, graph->et);
	double e2 = nuBaseRiboseEnergyCG(nodeB->baseConfCG, nodeA->riboseConfCG, sepR, graph->et);
	double e3 = nuRiboseRiboseEnergyCG(nodeA->riboseConfCG, nodeB->riboseConfCG, sep, graph->et);

	double e0Tmp = nuBaseBaseEnergyCG(nodeA->baseConfCGTmp, nodeB->baseConfCGTmp, sep, graph->et);
	double e1Tmp = nuBaseRiboseEnergyCG(nodeA->baseConfCGTmp, nodeB->riboseConfCGTmp, sep, graph->et);
	double e2Tmp = nuBaseRiboseEnergyCG(nodeB->baseConfCGTmp, nodeA->riboseConfCGTmp, sepR, graph->et);
	double e3Tmp = nuRiboseRiboseEnergyCG(nodeA->riboseConfCGTmp, nodeB->riboseConfCGTmp, sep, graph->et);

	if(abs(e0 - this->eneCG[0]) > 0.001){
		printf("ERROR cgBB %7.3f %7.3f\n", e0, this->eneCG[0]);
		tag = false;
	}
	if(abs(e1 - this->eneCG[1]) > 0.001){
		printf("ERROR cgBR %7.3f %7.3f\n", e1, this->eneCG[1]);
		tag = false;
	}
	if(abs(e2 - this->eneCG[2]) > 0.001){
		printf("ERROR cgRB %7.3f %7.3f\n", e2, this->eneCG[2]);
		tag = false;
	}
	if(abs(e3 - this->eneCG[3]) > 0.001){
		printf("ERROR cgRR %7.3f %7.3f\n", e3, this->eneCG[3]);
		tag = false;
	}

	if(abs(e0Tmp - this->eneCGTmp[0]) > 0.001){
		printf("ERROR cgBBTmp %7.3f %7.3f\n", e0, this->eneCGTmp[0]);
		tag = false;
	}
	if(abs(e1Tmp - this->eneCGTmp[1]) > 0.001){
		printf("ERROR cgBRTmp %7.3f %7.3f\n", e1, this->eneCGTmp[1]);
		tag = false;
	}
	if(abs(e2Tmp - this->eneCGTmp[2]) > 0.001){
		printf("ERROR cgRBTmp %7.3f %7.3f\n", e2, this->eneCGTmp[2]);
		tag = false;
	}
	if(abs(e3Tmp - this->eneCGTmp[3]) > 0.001){
		printf("ERROR cgRRTmp %7.3f %7.3f\n", e3, this->eneCGTmp[3]);
		tag = false;
	}
	return tag;
}

void NuEdge::printPartition(){
	int i,j,k;
	cout << "sep: " << this->sep << endl;
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
	this->poolSize = 100000;
	this->sampFreqNode = 1.0;
	this->sampFreqEdge = 1.0;
}

NuTree::~NuTree(){
	delete [] adjMtx;
}

void NuTree::updateNodeInfo(){
	int i,j;
	for(i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->updateNodeInformation(this);
	}
}

void NuTree::updateEdgeInfo(){


	vector<NuEdge*> tmpGeList;
	int i,j;

	for(i=0;i<graph->seqLen;i++){
		for(j=0;j<graph->seqLen;j++){
			if(i==j) continue;
			graph->allEdges[i*graph->seqLen+j]->updateEdgeInfo(this);
		}
	}


	for(i=0;i<geList.size();i++){
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
		pAdd += graph->geList[i]->samplingFreq/sampFreqEdge;
		end = (int)(pAdd*poolSize);
		for(j=start;j<end;j++){
			randPoolEdge[j] = i;
		}
	}

	double sum = sampFreqNode + sampFreqEdge;
	this->sampFreqNode = this->sampFreqNode/sum;
	this->sampFreqEdge = this->sampFreqEdge/sum;

}


void NuTree::printEdges(){
	for(int i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		printf("%-3d %-3d %7.3f\n", e->indexA, e->indexB, e->weight);
	}
}

void NuTree::runAtomicMC(){

	int stepNum = 10;

	double T0 = 5.0;
	double T1 = 0.001;
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

	for(T=T0;T>T1;T=T*anneal){

		nAc = 0;
		eAc = 0;
		nTot = 0;
		eTot = 0;

		for(k=0;k<stepNum;k++){
			randP = rand()*1.0/RAND_MAX;
			if(randP < sampFreqNode){
				/*
				 * rotamer mut
				 */
				nTot ++;

				randPos = randPoolNode[rand()%poolSize];
				randNode = graph->allNodes[randPos];
				randRot = graph->rotLib->riboseRotLib->getRandomRotamer(randNode->baseType);
				randNode->updateRiboseRotamer(randRot);
				mutE = randNode->rotMutEnergy();

				double mutE2 = graph->totalEnergyTmp() - graph->totalEnergy();


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
				randMove = randEdge->moveSet->getRandomMove();
				randEdge->updateCsMove(randMove);
				mutE = randEdge->mutEnergy();
				double mutE2 = graph->totalEnergyTmp() - graph->totalEnergy();

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
		double totEne2 = graph->totalEnergy2();

		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f totEne2=%8.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, totEne2);
	}

}

void NuTree::runCoarseGrainedMC(){

}



NuGraph::NuGraph(const string& inputFile){

	this->pairLib = new BasePairLib();
	this->rotLib = new RotamerLib();
	this->atLib = new AtomLib();
	this->moveLib = new NuPairMoveSetLibrary();
	this->et = new RnaEnergyTable();

	init(inputFile);
}

NuGraph::~NuGraph() {
	delete this->pairLib;
	delete this->rotLib;
	delete this->atLib;
	delete this->moveLib;

	delete [] seq;
	delete [] wcPairPosID;
	delete [] connectToDownstream;
	delete [] sepTable;
	for(int i=0;i<seqLen;i++){
		delete allNodes[i];
	}
	for(int i=0;i<seqLen*seqLen;i++){
		delete allEdges[i];
	}
	delete [] allNodes;
	delete [] allEdges;
	delete [] connectToDownstream;
	for(int i=0;i<seqLen;i++){
		delete initRiboseRotList[i];
	}

}

void NuGraph::init(const string& inputFile){
	int i,j,k;

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
	vector<string> templates = input.getMultiValues("template");
	vector<string> templatesAlignA = input.getMultiValues("alignNat");
	vector<string> templatesAlignB = input.getMultiValues("alignTmp");
	vector<string> templatesType = input.getMultiValues("tempType");

	seqLen = baseSeq.length();

	cout << "sequence length: " <<  seqLen << endl;

	this->seq = new int[seqLen];
	this->wcPairPosID = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->sepTable = new int[seqLen*seqLen];
	this->allNodes = new NuNode*[seqLen];
	this->allEdges = new NuEdge*[seqLen*seqLen];

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

		this->initRiboseRotList.push_back(rot);
	}


	for(i=0;i<seqLen;i++){
		LocalFrame cs1 = baseList[i]->getCoordSystem();
		this->allNodes[i] = new NuNode(i, baseList[i]->baseTypeInt, cs1, initRiboseRotList[i], atLib);
		this->allNodes[i]->graph = this;
	}

	cout << "parse secondary structure" << endl;
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
				this->allEdges[i*seqLen+j]->initNativeMoveSet();
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
					}
				}
			}
		}
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			this->geList.push_back(allEdges[i*seqLen+j]);
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
		allNodes[i]->checkEnergy();
	}
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			if(i==j) continue;
			allEdges[i*seqLen+j]->checkEnergy();
			allEdges[i*seqLen+j]->checkReversePair();
		}
	}
}

double NuGraph::totalEnergy(){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		ene += allNodes[i]->riboseConf->rot->energy;
		ene += allNodes[i]->phoConf->ene;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
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

double NuGraph::totalEnergyTmp(){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		ene += allNodes[i]->riboseConfTmp->rot->energy;
		ene += allNodes[i]->phoConfTmp->ene;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
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
		ene += allNodes[i]->ene;

	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
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


} /* namespace NSPpredNA */
