/*
 * NuGraph.cpp
 *
 *  Created on: 2023
 *      Author: nuc
 */

#include <predNA/NuGraph.h>
#include "model/AssignRNASasa.h"

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

void NuNode::updateEnergy(double clashRescale, double connectRescale){
	this->ene = riboseConf->rot->energy;
	if(connectToNeighbor)
		this->ene += phoConf->ene * connectRescale;

	this->eneTmp = this->ene;
}

int NuNode::hbondNumTo(NuNode* other, AtomLib* atLib){
	int hbNum  = 0;
	if(this->baseConf->cs1.origin_.squaredDistance(other->baseConf->cs1.origin_) > 225)
		return 0;
	//base - base
	for(int i=0;i<this->baseConf->rot->polarAtomNum;i++){
		int uniqueA = this->baseConf->rot->polarAtomUniqueID[i];
		XYZ coordA = this->baseConf->coords[this->baseConf->rot->polarAtomIndex[i]];
		
		for(int j=0;j<other->baseConf->rot->polarAtomNum;j++){
			int uniqueB = other->baseConf->rot->polarAtomUniqueID[j];
			XYZ coordB = other->baseConf->coords[other->baseConf->rot->polarAtomIndex[j]];		
			double d = coordA.distance(coordB);
			if(d < 2.0 || d > 3.4) continue;
			int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
			int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

			int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
			int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

			if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
				hbNum++;
		}
	}
	//base - ribose
	for(int i=0;i<this->baseConf->rot->polarAtomNum;i++){
		int uniqueA = this->baseConf->rot->polarAtomUniqueID[i];
		XYZ coordA = this->baseConf->coords[this->baseConf->rot->polarAtomIndex[i]];
		int uniqueB = 177; //O2'
		XYZ coordB = other->riboseConf->o2Polar.origin_;
		double d = coordA.distance(coordB);
		if(d < 2.0 || d > 3.4) continue;
		int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
		int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

		int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
		int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

		if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
			hbNum++;
	}
	//base - pho
	for(int i=0;i<this->baseConf->rot->polarAtomNum;i++){
		int uniqueA = this->baseConf->rot->polarAtomUniqueID[i];
		XYZ coordA = this->baseConf->coords[this->baseConf->rot->polarAtomIndex[i]];
		int uniqueB = 168; //OP1
		{
			XYZ coordB = other->phoConf->op1Polar.origin_;
			double d = coordA.distance(coordB);
			if(d> 2.0 && d < 3.4){
				int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
				int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

				int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
				int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

				if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
					hbNum++;
			}
		}

		{
			XYZ coordB = other->phoConf->op2Polar.origin_;
			double d = coordA.distance(coordB);
			if(d> 2.0 && d < 3.4){
				int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
				int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

				int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
				int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

				if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
					hbNum++;
			}
		}

	}
	//ribose - base 
	for(int i=0;i<other->baseConf->rot->polarAtomNum;i++){
		int uniqueA = other->baseConf->rot->polarAtomUniqueID[i];
		XYZ coordA = other->baseConf->coords[other->baseConf->rot->polarAtomIndex[i]];
		int uniqueB = 177; //O2'
		XYZ coordB = this->riboseConf->o2Polar.origin_;
		double d = coordA.distance(coordB);
		if(d < 2.0 || d > 3.4) continue;
		int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
		int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

		int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
		int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

		if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
			hbNum++;
	}

	//ribose - ribose
	{
		int uniqueA = 177;
		int uniqueB = 177;
		XYZ coordA = this->riboseConf->o2Polar.origin_;
		XYZ coordB = other->riboseConf->o2Polar.origin_;
		double d = coordA.distance(coordB);
		if(d> 2.0 && d < 3.4){
			int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
			int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

			int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
			int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

			if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
				hbNum++;
		}
	}

	//ribose - pho
	{
		int uniqueA = 177;
		int uniqueB = 168;
		{
			XYZ coordA = this->riboseConf->o2Polar.origin_;
			XYZ coordB = other->phoConf->op1Polar.origin_;
			double d = coordA.distance(coordB);
			if(d> 2.0 && d < 3.4){
				int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
				int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

				int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
				int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

				if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
					hbNum++;
			}
		}

		{
			XYZ coordA = this->riboseConf->o2Polar.origin_;
			XYZ coordB = other->phoConf->op2Polar.origin_;
			double d = coordA.distance(coordB);
			if(d> 2.0 && d < 3.4){
				int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
				int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

				int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
				int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

				if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
					hbNum++;
			}
		}

	}
	//pho - base
	
	for(int i=0;i<other->baseConf->rot->polarAtomNum;i++){
			int uniqueA = other->baseConf->rot->polarAtomUniqueID[i];
			XYZ coordA = other->baseConf->coords[other->baseConf->rot->polarAtomIndex[i]];
			int uniqueB = 168; //OP1
			{
				XYZ coordB = this->phoConf->op1Polar.origin_;
				double d = coordA.distance(coordB);
				if(d> 2.0 && d < 3.4){
					int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
					int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

					int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
					int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

					if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
						hbNum++;
				}
			}

			{
				XYZ coordB = this->phoConf->op2Polar.origin_;
				double d = coordA.distance(coordB);
				if(d> 2.0 && d < 3.4){
					int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
					int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

					int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
					int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

					if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
						hbNum++;
				}
			}

	}
	
	//pho - ribose 
	{
		int uniqueA = 177;
		int uniqueB = 168;
		{
			XYZ coordA = other->riboseConf->o2Polar.origin_;
			XYZ coordB = this->phoConf->op1Polar.origin_;
				double d = coordA.distance(coordB);
				if(d> 2.0 && d < 3.4){
					int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
					int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

					int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
					int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

					if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
						hbNum++;
				}
		}

		{
			XYZ coordA = other->riboseConf->o2Polar.origin_;
			XYZ coordB = this->phoConf->op2Polar.origin_;
				double d = coordA.distance(coordB);
				if(d> 2.0 && d < 3.4){
					int donorIndexA = atLib->uniqueIDToDonorID[uniqueA];
					int acceptorIndexA = atLib->uniqueIDToAcceptorID[uniqueA];

					int donorIndexB = atLib->uniqueIDToDonorID[uniqueB];
					int acceptorIndexB = atLib->uniqueIDToAcceptorID[uniqueB];

					if(donorIndexA >= 0 && acceptorIndexB >=0 || donorIndexB >=0 && acceptorIndexA >=0 )
						hbNum++;
				}
		}

	}
	return hbNum;
}

void NuNode::updateEnergyCG(double clashRescale, double connectRescale){
	this->eneCG = riboseConfCG->rot->energy;

	if(connectToNeighbor) {
		this->bbcg = nuConnectionEnergyCG(this->riboseConfCG, graph->allNodes[seqID+1]->riboseConfCG, graph->et, connectRescale);
		this->eneCG += this->bbcg;
	}
	this->bbcgTmp = this->bbcg;
	this->eneCGTmp = this->eneCG;
}

void NuNode::updateNodeInformation(NuTree* tree, double clashRescale, double connectRescale){

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
		this->ene += phoConf->ene * connectRescale;

	this->eneTmp = this->ene;
}

void NuNode::updateNodeInformationCG(NuTree* tree, double clashRescale, double connectRescale){

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
		this->eneCG += nuConnectionEnergyCG(this->riboseConfCG, graph->allNodes[seqID+1]->riboseConfCG, graph->et, connectRescale) ;
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

void NuNode::updateRiboseRotamer(RiboseRotamer* rot, double clashRescale, double connectRescale){
	int len = graph->seqLen;
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;

	this->riboseConfTmp->updateRotamer(rot);
	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->et->pb->buildPhosphate(graph->allNodes[j]->riboseConfTmp, graph->allNodes[j+1]->riboseConfTmp, graph->allNodes[j]->phoConfTmp);
		graph->allNodes[j]->eneTmp = graph->allNodes[j]->riboseConfTmp->rot->energy + graph->allNodes[j]->phoConfTmp->ene * connectRescale;

	}
	this->eneTmp = rot->energy + phoConfTmp->ene * connectRescale;

	/*
	 * base-ribose energy
	 */

	for(i=0;i<baseGroupA.size();i++){
		sep = graph->sepTable[baseGroupA[i]->seqID*len+seqID];
		egA = graph->allEdges[baseGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + baseGroupA[i]->seqID];

		egA->pairEneTmp[1] = nuBaseRiboseEnergy(baseGroupA[i]->baseConfTmp, this->riboseConfTmp, sep, graph->et, clashRescale);
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

			egA->pairEneTmp[2] = nuBasePhoEnergy(baseGroupA[i]->baseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et, clashRescale);
			egB->pairEneTmp[6] = egA->pairEneTmp[2];
		}
	}

	for(j=0;j<phoGroupC.size();j++){
		sep = graph->sepTable[seqID*len + phoGroupC[j]->seqID];
		if(sep == 0) continue;
		egA = graph->allEdges[seqID*len + phoGroupC[j]->seqID];
		egB = graph->allEdges[phoGroupC[j]->seqID*len + seqID];

		egA->pairEneTmp[2] = nuBasePhoEnergy(baseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et, clashRescale);
		egB->pairEneTmp[6] = egA->pairEneTmp[2];
	}


	/*
	 * ribose-ribose
	 */

	for(i=0;i<riboGroupA.size();i++){

		sep = graph->sepTable[riboGroupA[i]->seqID*len+seqID];

		egA = graph->allEdges[riboGroupA[i]->seqID*len + seqID];
		egB = graph->allEdges[seqID*len + riboGroupA[i]->seqID];

		egA->pairEneTmp[4] = nuRiboseRiboseEnergy(riboGroupA[i]->riboseConfTmp, this->riboseConfTmp, sep, graph->et, clashRescale);
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

			egA->pairEneTmp[5] = nuRibosePhoEnergy(riboGroupA[i]->riboseConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et, clashRescale);
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
		egA->pairEneTmp[5] = nuRibosePhoEnergy(riboseConfTmp, phoGroupA[i]->phoConfTmp, sep, graph->et, clashRescale);
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
			egA->pairEneTmp[8] = nuPhoPhoEnergy(phoGroupA[i]->phoConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et, clashRescale);
			egB->pairEneTmp[8] = egA->pairEneTmp[8];
		}
	}

	for(i=0;i<phoGroupC.size();i++){
		for(j=i+1;j<phoGroupC.size();j++){
			sep = graph->sepTable[phoGroupC[i]->seqID*len+phoGroupC[j]->seqID];
			egA = graph->allEdges[phoGroupC[i]->seqID*len+phoGroupC[j]->seqID];
			egB = graph->allEdges[phoGroupC[j]->seqID*len+phoGroupC[i]->seqID];
			egA->pairEneTmp[8] = nuPhoPhoEnergy(phoGroupC[i]->phoConfTmp, phoGroupC[j]->phoConfTmp, sep, graph->et, clashRescale);
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

double NuNode::rotMutEnergy(double connectRescale){
	int i,j, sep;
	int len = graph->seqLen;
	NuEdge* egA;

	double mutE = 0.0;
	mutE = riboseConfTmp->rot->energy - riboseConf->rot->energy;

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		mutE += (graph->allNodes[j]->phoConfTmp->ene - graph->allNodes[j]->phoConf->ene)*connectRescale;
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

bool NuNode::checkEnergy(double clashRescale, double connectRescale){
	double e1 = this->riboseConf->rot->energy + this->phoConf->ene*connectRescale;
	double e2 = this->riboseConfTmp->rot->energy + this->phoConfTmp->ene*connectRescale;
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

void NuNode::updateRiboseRotamerCG(RiboseRotamerCG* rot, double clashRescale, double connectRescale){
	int i,j,sep, sepR;
	NuEdge* egA;
	NuEdge* egB;

	int len = graph->seqLen;
	this->riboseConfCGTmp->updateRotamer(rot);

	for(i=0;i<connectionBreakPoints.size();i++){
		j = connectionBreakPoints[i];
		graph->allNodes[j]->bbcgTmp = nuConnectionEnergyCG(graph->allNodes[j]->riboseConfCGTmp, graph->allNodes[j+1]->riboseConfCGTmp, graph->et,connectRescale);
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

		egA->pairEneCGTmp[1] = nuBaseRiboseEnergyCG(baseGroupA[i]->baseConfCGTmp, this->riboseConfCGTmp, sep, graph->et, clashRescale);
		egB->pairEneCGTmp[2] = egA->pairEneCGTmp[1];
		egA->pairEneCGTmp[3] = nuRiboseRiboseEnergyCG(baseGroupA[i]->riboseConfCGTmp, this->riboseConfCGTmp, sep, graph->et, clashRescale);
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

bool NuNode::checkEnergyCG(double clashRescale, double connectRescale){
	double e3 = this->riboseConfCG->rot->energy;
	bool tag = true;

	if(connectToNeighbor)
		e3 += nuConnectionEnergyCG(this->riboseConfCG, graph->allNodes[seqID+1]->riboseConfCG, graph->et, connectRescale);

	double e4 = this->riboseConfCGTmp->rot->energy;
	if(connectToNeighbor)
		e4 += nuConnectionEnergyCG(this->riboseConfCGTmp, graph->allNodes[seqID+1]->riboseConfCGTmp, graph->et, connectRescale);

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
	vector<string> names;
	atLib->getRnaSidechainAtoms(this->baseType, names);
	vector<XYZ> tList;
	for(int i=0;i<baseConf->rot->atomNum;i++){
		tList.push_back(baseConf->coords[i]);
	}
    vector<Atom*> atomList;

    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names.at(i), tList[i]));
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

vector<Atom*> NuNode::toPhoAtomList(AtomLib* atLib){
    vector<Atom*> atomList;
    {
        atomList.push_back(new Atom("P", phoConf->coords[0]));
        atomList.push_back(new Atom("O5'", phoConf->coords[1]));
        atomList.push_back(new Atom("OP1", phoConf->coords[2]));
        atomList.push_back(new Atom("OP2", phoConf->coords[3]));
    }

    return atomList;
}

vector<Atom*> NuNode::toAtomListOnlyBase(AtomLib* atLib){
	vector<Atom*> list;
	vector<string> names;
	atLib->getRnaSidechainAtoms(this->baseType, names);
	vector<XYZ> tList;
	for(int i=0;i<baseConf->rot->atomNum;i++){
		tList.push_back(baseConf->coords[i]);
	}
    vector<Atom*> atomList;

    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names.at(i), tList[i]));
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
	vector<string> names;
	atLib->getRnaSidechainAtoms(this->baseType, names);
	vector<XYZ> tList;
	for(int i=0;i<baseConf->rot->atomNum;i++){
		tList.push_back(baseConf->coords[i]);
	}
    vector<Atom*> atomList;


    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names.at(i), tList[i]));
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
	this->moveSet = NULL;


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
	this->moveSet = NULL;
	//this->moveSet = new MixedNuPairCluster(sep, typeA*4+typeB, moveLib);
	//this->moveSet->updateEdgeInformation(this->ei);

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
	if(this->moveSet != NULL)
		delete this->moveSet;
}

void NuEdge::initNearNativeMoveSet(double distanceCutoff){

	BaseDistanceMatrix dm(nodeA->baseConf->cs1, nodeB->baseConf->cs1);
	initMoveSet(dm, distanceCutoff);

}

void NuEdge::initMoveSet(BaseDistanceMatrix& dm, double distanceCutoff){
	vector<int> neighborClusters;
	vector<double> distanceToClusterCenters;
	vector<double> pList;

	if(distanceCutoff == 0.0){
		pairLib->getNearestCluster(dm, nodeA->baseType, nodeB->baseType, sep, neighborClusters, distanceToClusterCenters);
	}
	else {
		pairLib->getNeighborClusters(dm, nodeA->baseType, nodeB->baseType, sep, neighborClusters, distanceToClusterCenters, distanceCutoff);
	}

	if(neighborClusters.size() > 0) {
		double pSum = 0;
		for(int i=0;i<distanceToClusterCenters.size();i++) {
			double d = distanceToClusterCenters[i];
			
			pSum += exp(-d*2);
		}

		for(int i=0;i<distanceToClusterCenters.size();i++) {
			double d = distanceToClusterCenters[i];
			//printf("cluster: %4d distanceToCluster: %5.3f probability: %7.5f\n", neighborClusters[i], d, exp(-d*2)/pSum);
			pList.push_back(exp(-d*2)/pSum);
		}
	}

	this->ei->setClusterList(neighborClusters, pList, pairLib);
	//this->moveSet->updateEdgeInformation(ei);
	
	/*
	if(this->ei->weight < 0) {
		this->ei->print();
		this->moveSet->printMoveSetInfo();
	}
	*/

	this->weight = this->ei->weight;
	this->weightRand = this->weight;
	this->samplingFreq = this->ei->validClusterNum;
	if(this->samplingFreq > 20) 
		this->samplingFreq = 20;
}

void NuEdge::fixNaiveMove(){
	CsMove cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->ei->setFixed(cm);
	//this->moveSet->fixNativeMove(cm);
	this->weight = this->ei->weight;
	this->weightRand = this->weight;
	this->samplingFreq = 0.0;
}

void NuEdge::updateEnergy(double clashRescale, double connectRescale){
	int sepR = -sep;
	if(abs(sepR) > 1) sepR = 2;

	if(nodeA->seqID < nodeB->seqID)
		pairEne[0] = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, graph->et, clashRescale);
	else
		pairEne[0] = nuBaseBaseEnergy(nodeB->baseConf, nodeA->baseConf, sepR, graph->et, clashRescale);

	pairEne[1] = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, graph->et, clashRescale);
	pairEne[2] = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, graph->et, clashRescale);

	pairEne[3] = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, graph->et, clashRescale);
	
	pairEne[4] = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, graph->et, clashRescale);
	
	pairEne[5] = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, graph->et, clashRescale);
	pairEne[6] = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, graph->et, clashRescale);
	pairEne[7] = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, graph->et, clashRescale);
	pairEne[8] = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, graph->et, clashRescale);

	for(int i=0;i<9;i++){
		pairEneTmp[i] = pairEne[i];
	}	
}

void NuEdge::updateEnergyCG(double clashRescale, double connectRescale){
	int sepR = -sep;
	if(abs(sepR) > 1) sepR = 2;

	if(nodeA->seqID < nodeB->seqID)
		pairEneCG[0] = nuBaseBaseEnergyCG(nodeA->baseConfCG, nodeB->baseConfCG, sep, graph->et, clashRescale);
	else
		pairEneCG[0] = nuBaseBaseEnergyCG(nodeB->baseConfCG, nodeA->baseConfCG, sepR, graph->et, clashRescale);

	pairEneCG[1] = nuBaseRiboseEnergyCG(nodeA->baseConfCG, nodeB->riboseConfCG, sep, graph->et, clashRescale);
	pairEneCG[2] = nuBaseRiboseEnergyCG(nodeB->baseConfCG, nodeA->riboseConfCG, sepR, graph->et, clashRescale);
	pairEneCG[3] = nuRiboseRiboseEnergyCG(nodeA->riboseConfCG, nodeB->riboseConfCG, sep, graph->et, clashRescale);

	for(int i=0;i<4;i++){
		pairEneCGTmp[i] = pairEneCG[i];
	}
}

void NuEdge::updateEdgeInfo(NuTree* tree ,double clashRescale, double connectRescale){

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
		pairEne[0] = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, graph->et, clashRescale);
	else
		pairEne[0] = nuBaseBaseEnergy(nodeB->baseConf, nodeA->baseConf, sepR, graph->et, clashRescale);

	pairEne[1] = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, graph->et, clashRescale);
	pairEne[2] = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, graph->et, clashRescale);

	pairEne[3] = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, graph->et, clashRescale);
	
	pairEne[4] = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, graph->et, clashRescale);
	
	pairEne[5] = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, graph->et, clashRescale);
	pairEne[6] = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, graph->et, clashRescale);
	pairEne[7] = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, graph->et, clashRescale);
	pairEne[8] = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, graph->et, clashRescale);

	for(i=0;i<9;i++){
		pairEneTmp[i] = pairEne[i];
	}
}

void NuEdge::updateEdgeInfoCG(NuTree* tree,double clashRescale, double connectRescale){

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
		pairEneCG[0] = nuBaseBaseEnergyCG(nodeA->baseConfCG, nodeB->baseConfCG, sep, graph->et, clashRescale);
	else
		pairEneCG[0] = nuBaseBaseEnergyCG(nodeB->baseConfCG, nodeA->baseConfCG, sepR, graph->et, clashRescale);

	pairEneCG[1] = nuBaseRiboseEnergyCG(nodeA->baseConfCG, nodeB->riboseConfCG, sep, graph->et, clashRescale);
	pairEneCG[2] = nuBaseRiboseEnergyCG(nodeB->baseConfCG, nodeA->riboseConfCG, sepR, graph->et, clashRescale);
	pairEneCG[3] = nuRiboseRiboseEnergyCG(nodeA->riboseConfCG, nodeB->riboseConfCG, sep, graph->et, clashRescale);

	for(i=0;i<4;i++){
		pairEneCGTmp[i] = pairEneCG[i];
	}
}

void NuEdge::updateCsMove(CsMove& cm, double clashRescale, double connectRescale){
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
		graph->allNodes[j]->eneTmp = graph->allNodes[j]->riboseConfTmp->rot->energy + graph->allNodes[j]->phoConfTmp->ene * connectRescale;
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
				eA->pairEneTmp[0] = nuBaseBaseEnergy(nA->baseConfTmp, nB->baseConfTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneTmp[0] = nuBaseBaseEnergy(nB->baseConfTmp, nA->baseConfTmp, sepR, graph->et, clashRescale);

			eA->pairEneTmp[1] = nuBaseRiboseEnergy(nA->baseConfTmp, nB->riboseConfTmp, sep, graph->et, clashRescale);
			eA->pairEneTmp[3] = nuBaseRiboseEnergy(nB->baseConfTmp, nA->riboseConfTmp, sepR, graph->et, clashRescale);

			if(nA->seqID < nB->seqID)
				eA->pairEneTmp[4] = nuRiboseRiboseEnergy(nA->riboseConfTmp, nB->riboseConfTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneTmp[4] = nuRiboseRiboseEnergy(nB->riboseConfTmp, nA->riboseConfTmp, sepR, graph->et, clashRescale);

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
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
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
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
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
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
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
				eA->pairEneTmp[2] = nuBasePhoEnergy(nA->baseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
				eA->pairEneTmp[5] = nuRibosePhoEnergy(nA->riboseConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
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
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et, clashRescale);

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
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et, clashRescale);

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
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et, clashRescale);
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
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nA->phoConfTmp, nB->phoConfTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneTmp[8] = nuPhoPhoEnergy(nB->phoConfTmp, nA->phoConfTmp, sepR, graph->et, clashRescale);

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

void NuEdge::updateCsMoveCG(CsMove& cm, double clashRescale, double connectRescale) {

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
		graph->allNodes[j]->bbcgTmp = nuConnectionEnergyCG(graph->allNodes[j]->riboseConfCGTmp, graph->allNodes[j+1]->riboseConfCGTmp, graph->et, connectRescale);
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
				eA->pairEneCGTmp[0] = nuBaseBaseEnergyCG(nA->baseConfCGTmp, nB->baseConfCGTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneCGTmp[0] = nuBaseBaseEnergyCG(nB->baseConfCGTmp, nA->baseConfCGTmp, sepR, graph->et, clashRescale);

			eA->pairEneCGTmp[1] = nuBaseRiboseEnergyCG(nA->baseConfCGTmp, nB->riboseConfCGTmp, sep, graph->et, clashRescale);
			eA->pairEneCGTmp[2] = nuBaseRiboseEnergyCG(nB->baseConfCGTmp, nA->riboseConfCGTmp, sepR, graph->et, clashRescale);

			if(nA->seqID < nB->seqID)
				eA->pairEneCGTmp[3] = nuRiboseRiboseEnergyCG(nA->riboseConfCGTmp, nB->riboseConfCGTmp, sep, graph->et, clashRescale);
			else
				eA->pairEneCGTmp[3] = nuRiboseRiboseEnergyCG(nB->riboseConfCGTmp, nA->riboseConfCGTmp, sepR, graph->et, clashRescale);

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

void NuEdge::printMutEnergyCG(){
	int i,j;
	NuNode* nA;
	NuNode* nB;
	NuEdge* eA;

	double mutE = 0.0;

	for(i=0;i<phoGroupC.size();i++){
		mutE += phoGroupC[i]->eneCGTmp - phoGroupC[i]->eneCG;
		cout << phoGroupC[i]->seqID << " " << phoGroupC[i]->eneCG << " " << phoGroupC[i]->eneCGTmp << endl;
		
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

			printf("edge: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", nodeListA[i]->seqID, nodeListA[j]->seqID, eA->sep, eA->pairEneCG[0], eA->pairEneCG[1], eA->pairEneCG[2],eA->pairEneCG[3]);
			printf("eTmp: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", nodeListA[i]->seqID, nodeListA[j]->seqID, eA->sep, eA->pairEneCGTmp[0], eA->pairEneCGTmp[1], eA->pairEneCGTmp[2],eA->pairEneCGTmp[3]);
		}
	}
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

bool NuEdge::checkEnergy(double clashRescale, double connectRescale){

	bool tag = true;
	int sepR = -sep;
	if(abs(sepR) > 1)
		sepR = 2;

	double e0;
	if(nodeA->seqID < nodeB->seqID)
		e0 = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, graph->et, clashRescale);
	else
		e0 = nuBaseBaseEnergy(nodeB->baseConf, nodeA->baseConf, sepR, graph->et, clashRescale);

	double e1 = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, graph->et, clashRescale);
	double e2 = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, graph->et, clashRescale);
	double e3 = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, graph->et, clashRescale);
	double e4 = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, graph->et, clashRescale);
	double e5 = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, graph->et, clashRescale);
	double e6 = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, graph->et, clashRescale);
	double e7 = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, graph->et, clashRescale);
	double e8 = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, graph->et, clashRescale);

	double e0Tmp;
	if(nodeA->seqID < nodeB->seqID)
		e0Tmp = nuBaseBaseEnergy(nodeA->baseConfTmp, nodeB->baseConfTmp, sep, graph->et, clashRescale);
	else
		e0Tmp = nuBaseBaseEnergy(nodeB->baseConfTmp, nodeA->baseConfTmp, sepR, graph->et, clashRescale);

	double e1Tmp = nuBaseRiboseEnergy(nodeA->baseConfTmp, nodeB->riboseConfTmp, sep, graph->et, clashRescale);
	double e2Tmp = nuBasePhoEnergy(nodeA->baseConfTmp, nodeB->phoConfTmp, sep, graph->et, clashRescale);
	double e3Tmp = nuBaseRiboseEnergy(nodeB->baseConfTmp, nodeA->riboseConfTmp, sepR, graph->et, clashRescale);
	double e4Tmp = nuRiboseRiboseEnergy(nodeA->riboseConfTmp, nodeB->riboseConfTmp, sep, graph->et, clashRescale);
	double e5Tmp = nuRibosePhoEnergy(nodeA->riboseConfTmp, nodeB->phoConfTmp, sep, graph->et, clashRescale);
	double e6Tmp = nuBasePhoEnergy(nodeB->baseConfTmp, nodeA->phoConfTmp, sepR, graph->et, clashRescale);
	double e7Tmp = nuRibosePhoEnergy(nodeB->riboseConfTmp, nodeA->phoConfTmp, sepR, graph->et, clashRescale);
	double e8Tmp = nuPhoPhoEnergy(nodeA->phoConfTmp, nodeB->phoConfTmp, sep, graph->et, clashRescale);

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

bool NuEdge::checkEnergyCG(double clashRescale, double connectRescale){

	bool tag = true;
	int sepR = -sep;
	if(abs(sepR) > 1)
		sepR = 2;

	double e0;
	if(nodeA->seqID < nodeB->seqID)
		e0 = nuBaseBaseEnergyCG(nodeA->baseConfCG, nodeB->baseConfCG, sep, graph->et, clashRescale);
	else
		e0 = nuBaseBaseEnergyCG(nodeB->baseConfCG, nodeA->baseConfCG, sepR, graph->et, clashRescale);

	double e1 = nuBaseRiboseEnergyCG(nodeA->baseConfCG, nodeB->riboseConfCG, sep, graph->et, clashRescale);
	double e2 = nuBaseRiboseEnergyCG(nodeB->baseConfCG, nodeA->riboseConfCG, sepR, graph->et, clashRescale);
	double e3 = nuRiboseRiboseEnergyCG(nodeA->riboseConfCG, nodeB->riboseConfCG, sep, graph->et, clashRescale);

	double e0Tmp;
	if(nodeA->seqID < nodeB->seqID)
		e0Tmp = nuBaseBaseEnergyCG(nodeA->baseConfCGTmp, nodeB->baseConfCGTmp, sep, graph->et, clashRescale);
	else
		e0Tmp = nuBaseBaseEnergyCG(nodeB->baseConfCGTmp, nodeA->baseConfCGTmp, sepR, graph->et, clashRescale);

	double e1Tmp = nuBaseRiboseEnergyCG(nodeA->baseConfCGTmp, nodeB->riboseConfCGTmp, sep, graph->et, clashRescale);
	double e2Tmp = nuBaseRiboseEnergyCG(nodeB->baseConfCGTmp, nodeA->riboseConfCGTmp, sepR, graph->et, clashRescale);
	double e3Tmp = nuRiboseRiboseEnergyCG(nodeA->riboseConfCGTmp, nodeB->riboseConfCGTmp, sep, graph->et, clashRescale);

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

bool NuEdge::isWC() {
	return graph->wcPairPosID[indexA] == indexB;
};  // check if this is currently a Watson-Crick edge

NuTree::NuTree(NuGraph* graph){
	this->graph = graph;
	this->adjMtx = new bool[graph->seqLen*graph->seqLen];
	this->poolSize = 100000;
	this->sampFreqNode = 0.5;
	this->sampFreqEdge = 0.5;
	this->totalSamp = 0.0;
}

NuTree::~NuTree(){
	delete [] adjMtx;
}

void NuTree::updateNodeInfo(double clashRescale, double connectRescale){
	int i,j;
	for(i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->updateNodeInformation(this ,  clashRescale,  connectRescale);
	}
}

void NuTree::updateNodeInfoCG(double clashRescale, double connectRescale){
	int i,j;
	for(i=0;i<graph->seqLen;i++){
		graph->allNodes[i]->updateNodeInformationCG(this,  clashRescale,  connectRescale);
	}
}

void NuTree::printNodeInfo(){
	cout << "node sampling freq: " << endl;
	for(int i=0;i<graph->seqLen;i++){
		cout << i << " " << graph->allNodes[i]->samplingFreq << endl;
		//cout << graph->allNodes[i]->connectToNeighbor << endl;
	}
}

void NuTree::updateEdgeInfo(double clashRescale, double connectRescale){
	cout << "geListSize: " << geList.size() << endl;
	vector<NuEdge*> tmpGeList;
	int i,j;

	for(i=0;i<graph->seqLen;i++){
		for(j=0;j<graph->seqLen;j++){
			if(i==j) continue;
			graph->allEdges[i*graph->seqLen+j]->updateEdgeInfo(this, clashRescale, connectRescale);
		}
	}

	for(i=0;i<geList.size();i++){
		tmpGeList.push_back(geList[i]);
	}

	geList.clear();

	for(i=0;i<tmpGeList.size();i++){
		NuEdge* eA = tmpGeList[i];
		NuEdge* eB = graph->allEdges[eA->indexB*graph->seqLen + eA->indexA];
		eA->updateEdgeInfo(this, clashRescale, connectRescale);
		eB->updateEdgeInfo(this, clashRescale, connectRescale);

		if(eA->nodeListA.size() < eA->nodeListB.size()){
			geList.push_back(eB);
		}
		else {
			geList.push_back(eA);
		}
	}

	for(i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		if(e->moveSet == NULL) {
			e->moveSet = new MixedNuPairCluster(e->sep, e->nodeA->baseType*4+e->nodeB->baseType, graph->moveLib);
		}
		else {
			cout << e->moveSet->pairType << endl;
		}
		e->moveSet->updateEdgeInformation(e->ei);
		printf("edge: %3d-%-3d %12s clusterNum: %4d\n", e->indexA, e->indexB, e->ei->ssSepKey.c_str(), e->ei->validClusterNum);
	}

}

void NuTree::updateEdgeInfoCG(double clashRescale, double connectRescale){


	vector<NuEdge*> tmpGeList;
	int i,j;

	for(i=0;i<graph->seqLen;i++){
		for(j=0;j<graph->seqLen;j++){
			if(i==j) continue;
			graph->allEdges[i*graph->seqLen+j]->updateEdgeInfoCG(this, clashRescale, connectRescale);
		}
	}


	for(i=0;i<geList.size();i++){
		tmpGeList.push_back(geList[i]);
	}
	geList.clear();

	for(i=0;i<tmpGeList.size();i++){
		NuEdge* eA = tmpGeList[i];
		NuEdge* eB = graph->allEdges[eA->indexB*graph->seqLen + eA->indexA];
		eA->updateEdgeInfoCG(this, clashRescale, connectRescale);
		eB->updateEdgeInfoCG(this, clashRescale, connectRescale);

		if(eA->nodeListA.size() < eA->nodeListB.size()){
			geList.push_back(eB);
		}
		else {
			geList.push_back(eA);
		}
	}

	for(i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		if(e->moveSet == NULL) {
			e->moveSet = new MixedNuPairCluster(e->sep, e->nodeA->baseType*4+e->nodeB->baseType, graph->moveLib);
		}
		e->moveSet->updateEdgeInformation(e->ei);
		printf("edge: %3d-%-3d %12s clusterNum: %4d\n", e->indexA, e->indexB, e->ei->ssSepKey.c_str(), e->ei->validClusterNum);
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

	if(sampFreqNode == 0 && sampFreqEdge == 0) {
		cout << "sampFreqNode == 0 && sampFreqEdge == 0" << endl;
		exit(0);
	}

	double pAdd = 0;
	int fgAdd = 0;
	

	/*
	 * node rand pool
	 */

	if(sampFreqNode > 0) {
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

		if(randPoolNode[poolSize-1] == 0){
			randPoolNode[poolSize-1] = randPoolNode[poolSize-2];
		}
	}

	/*
	 * edge rand pool
	 */
	for(i=0;i<poolSize;i++){
		this->randPoolEdge[i] = 0;
	}

	if(sampFreqEdge > 0) {
 		pAdd = 0;
		for(i=0;i<geList.size();i++){
			start = (int)(pAdd*poolSize);
			pAdd += 1.0*geList[i]->samplingFreq/sampFreqEdge;
			end = (int)(pAdd*poolSize);
			
			for(j=start;j<end;j++){
				randPoolEdge[j] = i;
			}
		}

		if(randPoolEdge[poolSize-1] == 0)
			randPoolEdge[poolSize-1] = randPoolEdge[poolSize-2];
	}


	double sum = sampFreqNode + sampFreqEdge;
	this->sampFreqNode = (this->sampFreqNode * graph->et->para->kNodeFreq) / (this->sampFreqNode * graph->et->para->kNodeFreq + sampFreqEdge);
	this->sampFreqEdge = this->sampFreqEdge / (this->sampFreqNode * graph->et->para->kNodeFreq + sampFreqEdge);
	if(this->sampFreqNode < 0.1)
	{
		this->sampFreqNode = 0.1;
		this->sampFreqEdge = 0.9;
	}	

	this->totalSamp = sum;
}

void NuTree::randomInit(double clashRescale, double connectRescale){
	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamer* randRot;

	CsMove oldMove;
	CsMove randMove;


	for(int i=0;i<graph->seqLen;i++){
		if(graph->allNodes[i]->samplingFreq == 0) continue;
		randNode = graph->allNodes[i];
		randRot = graph->rotLib->riboseRotLib->getRandomRotamer(randNode->baseType);
		randNode->updateRiboseRotamer(randRot, clashRescale, connectRescale);
		randNode->acceptRotMutation();
	}
	for(int i=0;i<geList.size();i++){

		randEdge = geList[i];
		//cout << "rand edge: " << randEdge->indexA << " " << randEdge->indexB << endl;
		//cout << "edge cluster num: " << randEdge->moveSet->clusterIDList.size() << endl;


		//cout << "rand edge sampling freq: " << randEdge->samplingFreq << endl;
		if(randEdge->samplingFreq == 0) continue;

		//oldMove = randEdge->cm;
		//BaseDistanceMatrix dmOld(oldMove);

		//cout << "rand move" << endl;
		randMove = randEdge->moveSet->getRandomMove();
		
		//BaseDistanceMatrix dmRand(randMove);
		//double moveDistance = dmOld.distanceTo(dmRand);
		//double minDistance = randEdge->moveSet->minDistanceToMoveSet(oldMove, graph->oi);

		//printf("rand move distance: %6.3f minDistaneToMoveLib: %6.3f\n", moveDistance, minDistance);
	
		randEdge->updateCsMove(randMove, clashRescale, connectRescale);

		randEdge->acceptMutation();

	}

}

void NuTree::randomInitCG(double clashRescale, double connectRescale){
	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamerCG* randRot;
	CsMove randMove;

	cout << "init rand node" << endl;
	for(int i=0;i<graph->seqLen;i++){
		if(graph->allNodes[i]->samplingFreq == 0) continue;
		randNode = graph->allNodes[i];
		randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
		randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
		randNode->acceptRotMutationCG();
	}

	cout << "init rand edge: " << geList.size() << endl;
	for(int i=0;i<geList.size();i++){
		
		randEdge = geList[i];
		if(randEdge->samplingFreq == 0) continue;

		//cout << "get rand move: " << i  << endl;
		randMove = randEdge->moveSet->getRandomMove();
		
		randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale );
		randEdge->acceptMutationCG();

		//cout << "rand cm clusterID: " << randEdge->cm.clusterID << endl;
	}

}

void NuTree::printEdges(){

	for(int i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		
		int clusterID = e->cm.clusterID;
		EdgeInformation* ei = e->ei;

		double p = ei->pContact;
		string ss = ei->ssSepKey;
		double w1 = e->weight;
		double w2 = e->weightRand;

		printf("%-3d %-3d %7.3f %7.3f %5.3f %10s %4d\n", e->indexA, e->indexB, e->weight, e->weightRand, p, ss.c_str(), clusterID);
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
	cout << "print edge info: " << endl;
	for(int i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		printf("%-3d %-3d %7.3f %8.3f %s\n", e->indexA, e->indexB, e->weight, e->samplingFreq, e->ei->ssSepKey.c_str());

	}
}

void NuTree::printTreeInfomation(const string& output) {
	ofstream out;
	out.open(output.c_str(), ios::out|ios::app);
	string augc = "AUGCatgc";
	string seq = "";
	string cnt = "";
	for(int i=0;i<graph->seqLen;i++){
		seq.append(augc.substr(graph->seq[i], 1));
		if(graph->connectToDownstream[i])
			cnt.append("-");
		else
			cnt.append("|");
	}

	out << "model " << seq << " " << cnt << " " << graph->getTotalEnergyForModelSelection() << endl;

	for(int i=0;i<graph->seqLen;i++){
		out <<  graph->allNodes[i]->riboseConf->rot->rotType << " ";
	}
	out << endl;

	for(int i=0;i<graph->seqLen;i++){
		if(graph->connectToDownstream[i])
			out <<  graph->allNodes[i]->phoConf->rot->dihed1 << " ";
		else 
			out << "0 ";
	}
	out << endl;	

	for(int i=0;i<graph->seqLen;i++){
		if(graph->connectToDownstream[i])
			out  <<  graph->allNodes[i]->phoConf->rot->dihed2 << " ";
		else 
			out << "0 " ;
	}
	out << endl;

	for(int i=1;i<graph->seqLen;i++){
		CsMove cm = graph->allNodes[i]->baseConf->cs1 - graph->allNodes[0]->baseConf->cs1;
		out << cm.toString() << endl;

	}

	out.close();

}

graphInfo* NuTree::runAtomicMC(){

	double clashRescale = 1.0;
	double connectRescale = 1.0;

	bool debug = false;

	XYZ n1n2Move;

	if(graph->et->para->withRandomInit)
		randomInit(clashRescale, connectRescale);

	double T0 = this->graph->et->para->T0;
	double T1 = this->graph->et->para->T1;
	double T2 = this->graph->et->para->T2;
	double T3 = this->graph->et->para->T3;

	int stepNum1 = (int)(this->totalSamp * graph->et->para->kStepNum1);
	int stepNum2 = (int)(this->totalSamp * graph->et->para->kStepNum2);
	int stepNum3 = (int)(this->totalSamp * graph->et->para->kStepNum3);
	
	double anneal = 0.95;

	if(debug) {
		stepNum1 = 100;
		stepNum2 = 100;
		stepNum3 = 100;
	}

	if(debug) {
		cout << "calculate total energy" << endl;
	}	
	
	clashRescale = graph->et->para->clashRescale;
	connectRescale = graph->et->para->connectRescale;

	graph->updateEnergy(clashRescale, connectRescale);
	double lamda = 1.07;

	
	double curEne = graph->totalEnergy(clashRescale, connectRescale);


	if(debug) {
		cout << "total energy: " << curEne << endl;
	}

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
		graph->checkEnergy(clashRescale, connectRescale);
	}

	
	cout << "start MC: " << endl;
	cout << "total step num: " << stepNum1+stepNum2+stepNum3 << endl;

	cout << "clash rescale: " << clashRescale << " connect rescale: " << connectRescale << endl;

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
				randNode->updateRiboseRotamer(randRot, clashRescale, connectRescale);

				if(debug) {
					cout << "rot mut energy" << endl;
				}
				mutE = randNode->rotMutEnergy(connectRescale);

				if(debug) {
					cout << "rot mut, pos: " << randPos << " mutE: " << mutE << endl;
					graph->checkEnergy(clashRescale, connectRescale);
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
					graph->checkEnergy(clashRescale, connectRescale);
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
					cout << "edge mut: edge: " << geList[randPos]->indexA << " " << geList[randPos]->indexB << endl;
				}
				
				randEdge = geList[randPos];
				//cout << randEdge->indexA << " " << randEdge->indexB << endl;
				randMove = randEdge->moveSet->getRandomMove();

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before update cs" << endl;
					graph->checkEnergy(clashRescale, connectRescale);
					cout << "move: " << endl;
					randMove.print();
				}
				randEdge->updateCsMove(randMove, clashRescale, connectRescale);

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before mutE" << endl;
					graph->checkEnergy(clashRescale, connectRescale);
				}

				//mutE = randEdge->mutEnergy(connectRescale);
				mutE = randEdge->mutEnergy();

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " mutE: " << mutE << endl;
					graph->checkEnergy(clashRescale, connectRescale);
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
					graph->checkEnergy(clashRescale, connectRescale);
				}
				
			}

			if(debug) {
				printf("step %d curEne: %8.3f eTot: %d totEne: %8.3f\n", k, curEne, eTot, graph->totalEnergy(clashRescale, connectRescale));
			}
		}
		
		double totEne = graph->totalEnergy(clashRescale, connectRescale);
		
		graphInfo* gi = graph->getGraphInfo();
		
		double rms = gi->rmsd(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f clashRescale: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms, clashRescale);

		clashRescale = clashRescale*lamda;
		connectRescale = connectRescale*lamda;
		if(clashRescale > 1.0)
			clashRescale = 1.0;
		if(connectRescale > 1.0)
			connectRescale = 1.0;
		curEne = graph->totalEnergy(clashRescale, connectRescale);
		graph->updateEnergy(clashRescale, connectRescale);
		
	}

	clashRescale = 1.0;
	connectRescale = 1.0;
	curEne = graph->totalEnergy(clashRescale, connectRescale);
	graph->updateEnergy(clashRescale, connectRescale);
	cout << "fixed subClusterID move" << endl;
	cout << "clash rescale: " << clashRescale << " connect rescale: " << connectRescale << endl;
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
				randNode->updateRiboseRotamer(randRot, clashRescale, connectRescale);
				mutE = randNode->rotMutEnergy(connectRescale);

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
				if(randEdge->cm.clusterID < 0) continue;
				
				randMove = randEdge->moveSet->getRandomMoveWithFixedSubCluster(randEdge->cm);
				randEdge->updateCsMove(randMove, clashRescale, connectRescale);
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

		double totEne = graph->totalEnergy(clashRescale, connectRescale);
		graphInfo* gi = graph->getGraphInfo();
		double rms = gi->rmsd(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
	}

	cout << "fixed SP1000ID move" << endl;
	cout << "clash rescale: " << clashRescale << " connect rescale: " << connectRescale << endl;
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
				randNode->updateRiboseRotamer(randRot, clashRescale, connectRescale);
				mutE = randNode->rotMutEnergy(connectRescale);

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
				randEdge->updateCsMove(randMove, clashRescale, connectRescale);

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
	
		double totEne = graph->totalEnergy(clashRescale, connectRescale);
		graphInfo* gi = graph->getGraphInfo();
		double rms = gi->rmsd(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
	}
	
	graphInfo* gi = graph->getGraphInfo();
	double rms = gi->rmsd(this->graph->initInfo);
	gi->setRMS(rms);
	return gi;
}

void NuTree::runCoarseGrainedMC(const string& output){
	
	int clashRescale = 1.0;
	int connectRescale = 1.0;

	bool debug = false;

	randomInitCG(clashRescale, connectRescale);

	printEdgeInfo();

	int stepNum = 100000;

	double T0 = 5.0;
	double T1 = 0.01;
	double anneal = 0.95;

	double curEne = graph->totalEnergyCG(clashRescale, connectRescale);
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

			if(debug) {
				cout << "step: " << k << " randP: " << randP << endl;
				cout << sampFreqNode << endl;
			}

			if(randP < sampFreqNode){
				/*
				 * rotamer mut
				 */
				nTot ++;
				randPos = randPoolNode[rand()%poolSize];

				randNode = graph->allNodes[randPos];
				randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
				randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
				mutE = randNode->rotMutEnergyCG();

				if(debug) {
					cout << "rot mut, pos: " << randPos << " mutE: " << mutE << endl;
					graph->checkEnergyCG(clashRescale, connectRescale);
					cout << "finish check" << endl;

					double mutE2 = graph->totalEnergyCGTmp(clashRescale, connectRescale) - graph->totalEnergyCG(clashRescale, connectRescale);
					cout << "mutE: " << mutE << " " << "mutE2: " << mutE2 << endl;
				}

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randNode->acceptRotMutationCG();
					curEne += mutE;
					nAc++;
					if(debug) {
						cout << "rot mut, pos: " << randPos << " AC " << endl;
						graph->checkEnergyCG(clashRescale, connectRescale);
						cout << "finish check" << endl;
					}
				}
				else {
					randNode->clearRotMutationCG();
					if(debug) {
						cout << "rot mut, pos: " << randPos << " RJ " << endl;
						graph->checkEnergyCG(clashRescale, connectRescale);
						
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
					graph->checkEnergyCG(clashRescale, connectRescale);
					cout << "finish check" << endl;
				}

				randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale);

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before mutE" << endl;
					graph->checkEnergyCG(clashRescale, connectRescale);
					cout << "finish check" << endl;
				}

				mutE = randEdge->mutEnergyCG();

				if(debug) {
					cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " mutE: " << mutE << endl;
					graph->checkEnergyCG(clashRescale, connectRescale);
					cout << "finish check" << endl;
					double mutE2 = graph->totalEnergyCGTmp(clashRescale, connectRescale) - graph->totalEnergyCG(clashRescale, connectRescale);
					cout << "mutE: " << mutE << " " << "mutE2: " << mutE2 << endl;

				}

				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					randEdge->acceptMutationCG();
					curEne += mutE;
					eAc++;
					if(debug) {
						cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " AC" << endl;
						graph->checkEnergyCG(clashRescale, connectRescale);
						cout << "finish check" << endl;
					}					
				}
				else {
					randEdge->clearMutationCG();
					if(debug) {
						cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " RJ" << endl;
						graph->checkEnergyCG(clashRescale, connectRescale);
						cout << "finish check" << endl;
					}
				}
			}
		}

		double totEne = graph->totalEnergyCG(clashRescale, connectRescale);
		graphInfo* gi = graph->getGraphInfoCG();
		double rms = gi->rmsdCG(this->graph->initInfo);
		delete gi;
		printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
	}

	cout << "total step num: " << count << endl;
	cout << "cg to all atom" << endl;
	graph->cgToAllAtom();

	cout << "get graph info" << endl;
	graphInfo* gi = graph->getGraphInfoCG();
	double rms = gi->rmsdCG(this->graph->initInfo);
	delete gi;

	cout << "pring pdb" << endl;
	gi = graph->getGraphInfo(0.0);
	gi->printPDB(output);
	delete gi;
	
}

graphInfo::graphInfo(int seqLen, int* seq, bool* con, bool* fixed, NuNode** nodes, double ene, AtomLib* atLib, int mode){
	this->seqLen = seqLen;
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->fixed = new bool[seqLen];
	this->nodes = new NuNode*[seqLen];
	this->ene = ene;
	this->sepTable = new int[seqLen*seqLen];

	for(int i=0;i<seqLen;i++){
		this->seq[i] = seq[i];
		this->connectToDownstream[i] = con[i];
		this->fixed[i] = fixed[i];
	}


	for(int i=0;i<seqLen;i++){
		NuNode* n = nodes[i];
		NuNode* node;
		if(mode == 0){ //all atom
			node = new NuNode(n->seqID, n->baseType,  n->baseConf->cs1, n->baseConf->rot, n->riboseConf->rot, atLib);
		}
		else { //coarse grained
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

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			int ij = i*seqLen+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else sepTable[ij] = 2;
		}
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
	delete [] sepTable;
}

double graphInfo::rmsd(graphInfo* other){
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

		vector<Atom*> aList = other->nodes[i]->toAtomList(this->atLib);

		for(int j=0;j<aList.size();j++)
			tList2.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}
	return NSPgeometry::simpleRMSD(tList1, tList2);
}

double graphInfo::rmsd(graphInfo* other, int pos){
	vector<XYZ> tList1;
	vector<XYZ> tList2;

			
	vector<Atom*> aList = nodes[pos]->toAtomListOnlyBase(this->atLib);
	for(int j=0;j<aList.size();j++)
		tList1.push_back(aList[j]->coord);
	for(int k=0;k<aList.size();k++){
		delete aList[k];
	}

	aList = other->nodes[pos]->toAtomListOnlyBase(this->atLib);
	for(int j=0;j<aList.size();j++)
		tList2.push_back(aList[j]->coord);
	for(int k=0;k<aList.size();k++){
		delete aList[k];
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
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
		atomNum += aList.size();
	}

	for(int i=0;i<this->seqLen-1;i++){
		if(!connectToDownstream[i]) continue;
		vector<Atom*> aList = nodes[i]->toPhoAtomList(atLib);
		RNABase* base = rc.getBaseList()[i+1];
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
	}

	vector<RNABase*> baseList = rc.getBaseList();

	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);

	string augc = "AUGCatgc";
	of << "seq: ";
	for(int i=0;i<this->seqLen;i++){
		of << augc[this->nodes[i]->baseType];
	}
	of << endl;

	of << "cnt: ";
	for(int i=0;i<this->seqLen;i++){
		if(connectToDownstream[i])
			of << "-";
		else 
			of << "|";
	}
	of << endl;

	of << "ene: " << this->ene << endl;
	of << "rms: " << this->rms << endl;
	of.close();


	for(int i=0;i<baseList.size();i++){
		vector<Atom*>* aList = baseList[i]->getAtomList();
		for(int j=0;j<aList->size();j++)	{
			delete aList->at(j);
		}
	}
}

void graphInfo::printPDBWithPairMtx(const string& outputFile, BasePairLib* bpLib){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	int seqID = 0;
	int atomNum = 0;

	for(int i=0;i<this->seqLen;i++) {
		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
		atomNum += aList.size();
	}

	for(int i=0;i<this->seqLen-1;i++){
		if(!connectToDownstream[i]) continue;
		vector<Atom*> aList = nodes[i]->toPhoAtomList(atLib);
		RNABase* base = rc.getBaseList()[i+1];
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
	}

	vector<RNABase*> baseList = rc.getBaseList();

	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);

	string augc = "AUGCatgc";
	of << "seq: ";
	for(int i=0;i<this->seqLen;i++){
		of << augc[this->nodes[i]->baseType];
	}
	of << endl;

	of << "cnt: ";
	for(int i=0;i<this->seqLen;i++){
		if(connectToDownstream[i])
			of << "-";
		else 
			of << "|";
	}
	of << endl;

	of << "mtx: " << this->seqLen << endl;

	int clusterIDList[this->seqLen][this->seqLen];
	for(int i=0;i<this->seqLen;i++){
		for(int j=0;j<this->seqLen;j++){
			clusterIDList[i][j] = -1;
		}
	}
	for(int i=0;i<this->seqLen;i++){
		for(int j=i+1;j<this->seqLen;j++){
			if(i==j) continue;
			if(!baseList[i]->contactTo(baseList[j], 5.0)) continue;
			int clusterID = -1;
			BaseDistanceMatrix dm(*baseList[i], *baseList[j]);
			if(j==i+1 && connectToDownstream[i]){
				clusterID = bpLib->getPairType(dm, this->nodes[i]->baseType, this->nodes[j]->baseType, 1);
			}
			else {
				clusterID = bpLib->getPairType(dm, this->nodes[i]->baseType, this->nodes[j]->baseType, 2);
			}
			
			
			clusterIDList[i][j] = clusterID;
			clusterIDList[j][i] = clusterID;
		}
	}

	for(int i=0;i<this->seqLen;i++){
		for(int j=0;j<this->seqLen;j++){
			sprintf(ss, "%-3d ", clusterIDList[i][j]);
			of << string(ss);
		}
		of << endl;
	}

	of << "ene: " << this->ene << endl;
	of << "rms: " << this->rms << endl;
	of.close();


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

void graphInfo::printDetailEnergy(const string& outputFile, BasePairLib* bpLib, AtomLib* atLib, RnaEnergyTable* et){
	ofstream out;
	out.open(outputFile.c_str(), ios::out);
	NuNode* nodeA;
	NuNode* nodeB;
	int sep, sepR;

	char xx[200];

	double totEnergy = 0.0;


	for(int i=0;i<seqLen;i++){
		nodeA = this->nodes[i];
		double eRot = nodeA->riboseConf->rot->energy;
		double ePho = 0.0;
		if(connectToDownstream[i]){
			ePho = nodeA->phoConf->ene;
		}

		totEnergy += eRot + ePho;

		sprintf(xx, "ROT %d %8.3f", nodeA->baseConf->rot->baseType, eRot);
		out << string(xx) << endl;
		sprintf(xx, "PHO %d %8.3f", nodeA->baseConf->rot->baseType, ePho);
		out << string(xx) << endl;

		double e1,e2;
		double BB, BR, BP, RB, RR, RP, PB, PR, PP;
		for(int j=i+1;j<seqLen;j++){
			nodeB = this->nodes[j];
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			e1 = nuBaseBaseEnergyPrintDetail(nodeA->baseConf, nodeB->baseConf, sep, et, bpLib, out);
			
			BB = e1;
			totEnergy += e1;

			e2 = nuBaseBaseEnergy(nodeA->baseConf, nodeB->baseConf, sep, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "base base error" << endl;
			}

			e1 = nuBaseRiboseEnergyPrintDetail(nodeA->baseConf, nodeB->riboseConf, sep, et, atLib, out);
			totEnergy += e1;
			BR = e1;
			e2 = nuBaseRiboseEnergy(nodeA->baseConf, nodeB->riboseConf, sep, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "base ribose error" << endl;
			}

			e1 = nuBaseRiboseEnergyPrintDetail(nodeB->baseConf, nodeA->riboseConf, sepR, et, atLib, out);
			totEnergy += e1;
			RB = e1;
			e2 = nuBaseRiboseEnergy(nodeB->baseConf, nodeA->riboseConf, sepR, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "ribose base error" << endl;
			}

			e1 = nuBasePhoEnergyPrintDetail(nodeA->baseConf, nodeB->phoConf, sep, et, atLib, out);
			totEnergy += e1;
			BP = e1;
			e2 = nuBasePhoEnergy(nodeA->baseConf, nodeB->phoConf, sep, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "base pho error" << endl;
			}


			e1 = nuBasePhoEnergyPrintDetail(nodeB->baseConf, nodeA->phoConf, sepR, et, atLib, out);
			totEnergy += e1;
			PB = e1;
			e2 = nuBasePhoEnergy(nodeB->baseConf, nodeA->phoConf, sepR, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "pho base error" << endl;
			}

			e1 = nuRiboseRiboseEnergyPrintDetail(nodeA->riboseConf, nodeB->riboseConf, sep, et, atLib, out);
			totEnergy += e1;
			RR = e1;
			e2 = nuRiboseRiboseEnergy(nodeA->riboseConf, nodeB->riboseConf, sep, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "ribose ribose error" << endl;
			}

			e1 = nuRibosePhoEnergyPrintDetail(nodeA->riboseConf, nodeB->phoConf, sep, et, atLib, out);
			totEnergy += e1;
			RP = e1;
			e2 = nuRibosePhoEnergy(nodeA->riboseConf, nodeB->phoConf, sep, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "ribose pho error" << endl;
			}

			e1 = nuRibosePhoEnergyPrintDetail(nodeB->riboseConf, nodeA->phoConf, sepR, et, atLib, out);
			totEnergy += e1;
			PR = e1;
			e2 = nuRibosePhoEnergy(nodeB->riboseConf, nodeA->phoConf, sepR, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "pho ribose error" << endl;
			}

			e1 = nuPhoPhoEnergyPrintDetail(nodeA->phoConf, nodeB->phoConf, sep, et, atLib, out);
			//out << "PP "<< i << " " << j << " " <<  e1 << endl;
			totEnergy += e1;
			PP = e1;
			e2 = nuPhoPhoEnergy(nodeA->phoConf, nodeB->phoConf, sep, et, 1.0);
			if(abs(e1-e2) > 0.001){
				cout << "pho pho error" << endl;
			}
			//printf("edge: %2d %2d sep: %2d BB %7.3f BR %7.3f BP %7.3f RB %7.3f RR %7.3f RP %7.3f PB %7.3f PR %7.3f PP %7.3f\n", i, j,sep, BB, BR, BP, RB, RR, RP, PB, PR, PP);
		}
	}

	cout << "total energy: " << totEnergy << endl;
	out.close();

}

NuGraph::NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et){

	this->pairLib = pairLib;
	this->rotLib = rotLib;
	this->atLib = atLib;
	this->moveLib = moveLib;
	this->et = et;
	this->eiLib = eiLib;
	this->initInfo = NULL;

}

NuGraph::NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib){
	this->pairLib = pairLib;
	this->rotLib = rotLib;
	this->atLib = atLib;
	this->moveLib = NULL;
	this->eiLib = NULL;
	this->et = NULL;
	this->initInfo = NULL;
	initForMST(inputFile);
}

NuGraph::NuGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, RnaEnergyTable* et, int initMode){
	this->pairLib = pairLib;
	this->rotLib = rotLib;
	this->atLib = atLib;
	this->moveLib = NULL;
	this->eiLib = NULL;
	this->et = et;
	this->initInfo = NULL;
	switch (initMode) {
		case 1:
			initForMC(inputFile);
			break;
		default:
			cerr << "Unknown initMode " << initMode << endl;
	}
}

NuGraph::NuGraph(RNAPDB* pdb, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib,  RnaEnergyTable* et, const string& cnt){
	this->pairLib = pairLib;
	this->rotLib = rotLib;
	this->atLib = atLib;
	this->moveLib = moveLib;
	this->et = et;
	this->eiLib = eiLib;

	AssignRNASS* ar = new AssignRNASS(pdb, atLib);

	string task = "refinement";
	vector<RNABase*> baseList = pdb->getBaseList();
	string baseSeq = ar->seq;
	string baseSec = ar->ssSeq;
	delete ar;

	int i,j,k;
	vector<string> spt;
	this->oi = new OrientationIndex();

	seqLen = baseSeq.length();
	if(baseList.size() != seqLen){
		cout << "invalid pdb" << endl;
		exit(0);
	}

	this->seq = new int[seqLen];
	this->wcPairPosID = new int[seqLen];
	this->stemIndex = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->sepTable = new int[seqLen*seqLen];
	this->allNodes = new NuNode*[seqLen];
	this->allEdges = new NuEdge*[seqLen*seqLen];
	this->fixed = new bool[seqLen];


	/*
	 * read chain break points
	 */
	for(i=0;i<seqLen;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[seqLen-1] = false;


	if(cnt.length() != seqLen){
		for(i=0;i<seqLen-1;i++){
			if(baseList[i]->connectToNeighbor(baseList[i+1]))
				connectToDownstream[i] = true;
			else {
				connectToDownstream[i] = false;
			}
		}
	}
	else {
		for(i=0;i<seqLen-1;i++){
			if(cnt[i] == '-')
				connectToDownstream[i] = true;
			else 
				connectToDownstream[i] = false;
		}
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

	/*
	 * init NuNodes
	 */

	if(baseList.size() != baseSeq.length()) {
		cout << "pdb size: " << baseList.size() << " not equal to sequence length: " << baseSeq.length() << endl;
		exit(0);
	}
	else if(baseSeq.length() != baseSec.length()) {
		cout << "seq length not equal to sec length" << endl;
		exit(0);
	}

	string augc = "AUGCatgc";
	string pdbSeq = "";
	for(i=0;i<seqLen;i++){
		if(baseList[i]->baseTypeInt < 0 || baseList[i]->baseTypeInt > 3){
			cout << "invalid base type: " << baseList[i]->baseType << endl;
			exit(0);
		}
		this->seq[i] = baseList[i]->baseTypeInt;
		pdbSeq = pdbSeq + augc.substr(seq[i], 1);

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

	//check pdb sequence
	if(baseSeq != pdbSeq) {
		cout << "pdb sequence is inconsistent with input sequence: " << endl;
		cout << pdbSeq << endl;
		cout << baseSeq << endl;
		exit(0);
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
		this->allNodes[i]->connectToNeighbor = connectToDownstream[i];
		this->allNodes[i]->graph = this;
	}

	/*
	 * parse secondary structure information
	 */
	cout << "parse secondary structure " << endl;

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

	int currentStemID = 1;
	int curStemLen = 0;
	for(int i=0;i<seqLen;i++){
		this->stemIndex[i] = 0;
	}
	for(int i=0;i<seqLen-1;i++) {
		if(wcPairPosID[i] < 0 ) continue;
		if(wcPairPosID[i] > i && wcPairPosID[i+1] == wcPairPosID[i]-1 && connectToDownstream[i] && connectToDownstream[wcPairPosID[i+1]]) {
			stemIndex[i] = currentStemID;
			stemIndex[i+1] = currentStemID;
			stemIndex[wcPairPosID[i]] = currentStemID;
			stemIndex[wcPairPosID[i+1]] = currentStemID;
			curStemLen ++;
		}
		else if(curStemLen > 0) {
			curStemLen = 0;
			currentStemID++;
		}
	}

	cout << "stems: " << endl;
	for(i=0;i<seqLen;i++){
		cout << stemIndex[i];
	}
	cout << endl;


	cout << "init edges" << endl;
	/*
	 * init NuEdges for refinement
	 */

	char xs[200];
	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			this->allEdges[i*seqLen+j] = new NuEdge(allNodes[i], allNodes[j], this);

			this->allEdges[i*seqLen+j]->graph = this;
			this->allEdges[i*seqLen+j]->weight = 0.0;

			sprintf(xs, "%c%c", augc[this->allNodes[i]->baseType], augc[this->allNodes[j]->baseType]);
			string pairType = string(xs);

			if(i==j) continue;

			string ssSepType = "UNK";
			int sep = this->allEdges[i*seqLen+j]->sep;
			int pairID1 = this->wcPairPosID[i];
			int pairID2 = this->wcPairPosID[j];
			this->allEdges[i*seqLen+j]->initNearNativeMoveSet(0.0);
		}
	}
	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			this->geList.push_back(allEdges[i*seqLen+j]);
		}
	}
	this->initInfo = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, 0.0, atLib, 0);	
}


NuGraph::~NuGraph() {

	if (oi != NULL){
		delete oi;
		oi = NULL;
	}
	if (seq != nullptr){
		delete [] seq;
		seq = nullptr;
	}
	if (wcPairPosID != nullptr){
		delete [] wcPairPosID;
		wcPairPosID = nullptr;
	}
	if(stemIndex != nullptr){
		delete [] stemIndex;
		stemIndex = nullptr;
	}
	if(connectToDownstream != nullptr){
		delete [] connectToDownstream;
		connectToDownstream = nullptr;
	}
	if (fixed != nullptr){
		delete [] fixed;	
		fixed = nullptr;
	}
	if (sepTable != nullptr){
		delete [] sepTable;
		sepTable = nullptr;
	}
	if (allNodes){
		cout << seqLen << endl;
		for(int i=0;i<seqLen;i++){	
			if (allNodes[i]!= NULL){
				delete allNodes[i];
				allNodes[i] = NULL;
			}
		}
		delete [] allNodes;
		allNodes = nullptr;
	}

	if (allEdges) {
		for(int i=0;i<seqLen*seqLen;i++){
			if (allEdges[i]!= NULL){
				delete allEdges[i];
				allEdges[i] = NULL;
			}
		}
		delete [] allEdges;
		allEdges = nullptr;
	}

	for(int i=0;i<seqLen;i++){
		delete initRiboseRotList[i];
		delete initBaseRotList[i];
	}
	if(initInfo != NULL){
		delete initInfo;
	}
}

void NuGraph::init(const string& task, RNAPDB* pdb, const string& baseSeq, const string& baseSec, const string& csn, const string& cst, const string& cnt, const string& contactKey, vector<string>& ctList){
	int i,j,k;
	vector<string> spt;

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 *      singleBasePredict: single base prediction
	 *      analysis: without energy calculation
	 */

	cout << "task: " << task << endl;  

	/*
	vector<string> templates = input.getMultiValues("template");
	vector<string> templatesAlignA = input.getMultiValues("alignNat");
	vector<string> templatesAlignB = input.getMultiValues("alignTmp");
	vector<string> templatesType = input.getMultiValues("tempType");
	*/

	this->oi = new OrientationIndex();

	seqLen = baseSeq.length();


	this->seq = new int[seqLen];
	this->wcPairPosID = new int[seqLen];
	this->stemIndex = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->sepTable = new int[seqLen*seqLen];
	this->allNodes = new NuNode*[seqLen];
	this->allEdges = new NuEdge*[seqLen*seqLen];
	this->fixed = new bool[seqLen];


	/*
	 * read chain break points
	 */
	for(i=0;i<seqLen;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[seqLen-1] = false;

	if(cnt.length() != seqLen) {
		cout << "invalid cnt string: " << cnt << endl;
		cout << "cnt length: " << cnt.length() << " " << seqLen << endl;
		exit(0);
	}
	for(i=0;i<seqLen-1;i++){
		if(cnt[i] == '|'){
			connectToDownstream[i] = false;
		}
		else if(cnt[i] == '-'){
			connectToDownstream[i] = true;
		}
		else {
			cout << "invalid cnt string: " << cnt << endl;
			cout << "cnt example: ------|-------" << endl;
			exit(0);
		}
	}

	cout << "sepTable" << endl;

	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			int ij = i*seqLen+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else sepTable[ij] = 2;
		}
	}

	/*
	 * init NuNodes
	 */
	cout << "init nodes" << endl;

	vector<RNABase*> baseList = pdb->getBaseList();

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
	else if(csn.length() != baseSeq.length()){
		cout << "invalid csn " << endl;
		cout << csn << endl;
		exit(1);
	}

	string augc = "AUGCatgc";
	string pdbSeq = "";
	for(i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		pdbSeq = pdbSeq + augc.substr(seq[i], 1);

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

	//check pdb sequence
	if(baseSeq != pdbSeq) {
		cout << "pdb sequence is inconsistent with input sequence: " << endl;
		cout << pdbSeq << endl;
		cout << baseSeq << endl;
		exit(0);
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
		this->allNodes[i]->connectToNeighbor = connectToDownstream[i];
		this->allNodes[i]->graph = this;
	}


	/*
	 * parse secondary structure information
	 */
	cout << "parse secondary structure " << endl;

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

	int currentStemID = 1;
	int curStemLen = 0;
	for(int i=0;i<seqLen;i++){
		this->stemIndex[i] = 0;
	}
	for(int i=0;i<seqLen-1;i++) {
		if(wcPairPosID[i] < 0 ) continue;
		if(wcPairPosID[i] > i && wcPairPosID[i+1] == wcPairPosID[i]-1 && connectToDownstream[i] && connectToDownstream[wcPairPosID[i+1]]) 
		{
			stemIndex[i] = currentStemID;
			stemIndex[i+1] = currentStemID;
			stemIndex[wcPairPosID[i]] = currentStemID;
			stemIndex[wcPairPosID[i+1]] = currentStemID;
			curStemLen ++;
		}
		else if(wcPairPosID[i] > i){
			stemIndex[i] = currentStemID;
			stemIndex[wcPairPosID[i]] = currentStemID;
			curStemLen ++;
		}
		else if(curStemLen > 0) {
			curStemLen = 0;
			currentStemID++;
		}
	}

	cout << "stems: " << endl;
	for(i=0;i<seqLen;i++){
		cout << stemIndex[i];
	}
	cout << endl;


	cout << "init edges" << endl;
	/*
	 * init NuEdges
	 */


	vector<int> h3List;
	vector<int> h5List;
	vector<int> h3IDList;
	vector<int> h5IDList;
	for(i=0;i<seqLen;i++) {
			
			if(stemIndex[i] == 0) {
				int h3Sep = -1;
				int h3Id = -1;
				int h5Sep = -1;
				int h5Id = -1;
				
				for(int j=i-1;j>=0;j--) {
					if(!connectToDownstream[j]) break;
					if(stemIndex[j] > 0) {
						h3Sep = (i-j);
						h3Id = j;
						break;
					}
				}
				
				for(int j=i;j<seqLen;j++) {
					if(j>0 && !connectToDownstream[j-1]) break;
					if(stemIndex[j] >0) {
						h5Sep = (j-i);
						h5Id = j;
						break;
					}
				}
				
				if(h3Sep > 5){
					h3Sep = -1;
					h3Id = -1;
				}
				if(h5Sep > 5){
					h5Sep = -1;
					h5Id = -1;
				}
				
				h3List.push_back(h3Sep);
				h3IDList.push_back(h3Id);
				h5List.push_back(h5Sep);
				h5IDList.push_back(h5Id);
			}
			else if(i< seqLen-1 && stemIndex[i]>0 && stemIndex[i+1]==0) {
				int h3Sep = 0;
				int h3Id = i;
				int h5Sep = -1;
				int h5Id = -1;
				h3List.push_back(h3Sep);
				h3IDList.push_back(h3Id);
				h5List.push_back(h5Sep);
				h5IDList.push_back(h5Id);
			}
			else if(i>0 && stemIndex[i]>0 && stemIndex[i-1]==0) {
				int h3Sep = -1;
				int h3Id = -1;
				int h5Sep = 0;
				int h5Id = i;
				h3List.push_back(h3Sep);
				h3IDList.push_back(h3Id);
				h5List.push_back(h5Sep);
				h5IDList.push_back(h5Id);
			}
			else {
				h3List.push_back(-1);
				h3IDList.push_back(-1);
				h5List.push_back(-1);
				h5IDList.push_back(-1);
			}
	}

	cout << "h3List: ";
	for(int i=0;i<seqLen;i++){
		cout << h3List[i] << " ";
	}
	cout << endl;

	cout << "h3IDList: ";
	for(int i=0;i<seqLen;i++){
		cout << h3IDList[i] << " ";
	}
	cout << endl;

	cout << "h5List: ";
	for(int i=0;i<seqLen;i++){
		cout << h5List[i] << " ";
	}
	cout << endl;

	cout << "h5IDList: ";
	for(int i=0;i<seqLen;i++){
		cout << h5IDList[i] << " ";
	}
	cout << endl;

	char xs[200];
	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			if(task == "predict"){
				this->allEdges[i*seqLen+j] = new NuEdge(allNodes[i], allNodes[j], this->sepTable[i*seqLen+j], this->pairLib, this->moveLib);
			}
			else {
				this->allEdges[i*seqLen+j] = new NuEdge(allNodes[i], allNodes[j], this);
			}

			this->allEdges[i*seqLen+j]->graph = this;
			this->allEdges[i*seqLen+j]->weight = 0.0;

			sprintf(xs, "%c%c", augc[this->allNodes[i]->baseType], augc[this->allNodes[j]->baseType]);
			string pairType = string(xs);

			if(i==j) continue;

			string ssSepType = "UNK";
			int sep = this->allEdges[i*seqLen+j]->sep;
			int pairID1 = this->wcPairPosID[i];
			int pairID2 = this->wcPairPosID[j];

			if(sep == 1) {
				ssSepType = "nb-"+pairType+"-";
				if(pairID1 >=0 && pairID2 >=0 && stemIndex[i] == stemIndex[j])
					ssSepType = ssSepType + "HH";
				else if(pairID1 >=0 && pairID2 >=0)
					ssSepType = ssSepType + "Hh";
				else if(pairID1 >= 0)
					ssSepType = ssSepType + "HL";
				else if(pairID2 >= 0)
					ssSepType = ssSepType + "LH";
				else 
					ssSepType = ssSepType + "LL";
			}
			else if(sep == -1) {
				ssSepType = "rnb-"+pairType+"-";
				if(pairID1 >=0 && pairID2 >=0 && stemIndex[i] == stemIndex[j])
					ssSepType = ssSepType + "HH";
				else if(pairID1 >=0 && pairID2 >=0)
					ssSepType = ssSepType + "hH";
				else if(pairID1 >= 0)
					ssSepType = ssSepType + "HL";
				else if(pairID2 >= 0)
					ssSepType = ssSepType + "LH";
				else 
					ssSepType = ssSepType + "LL";				
			}
			else if(sep == 2){
				
				int h3SepA = h3List[i];
				int h3IDA = h3IDList[i];
				int h5SepA = h5List[i];
				int h5IDA = h5IDList[i];

				int h3SepB = h3List[j];
				int h3IDB = h3IDList[j];

				int h5SepB = h5List[j];
				int h5IDB = h5IDList[j];

				if(wcPairPosID[i] == j){ //WC pair
					sprintf(xs, "wc-%s",pairType.c_str());
					ssSepType = string(xs);
				}
				else if(h3IDA >=0 && wcPairPosID[h3IDA] == h5IDB && h3IDB >=0 && wcPairPosID[h3IDB] == h5IDA){
					if(h3SepA + h5SepB < h3SepB + h5SepA){
						sprintf(xs, "nnb-%s-3P5P-%d-%d",pairType.c_str(), h3SepA, h5SepB);
						ssSepType = string(xs);
					}
					else {
						sprintf(xs, "nnb-%s-5P3P-%d-%d",pairType.c_str(), h5SepA, h3SepB);
						ssSepType = string(xs);
					}
				}
				else if(h3IDA >=0 && wcPairPosID[h3IDA] == h5IDB){
					sprintf(xs, "nnb-%s-3P5P-%d-%d",pairType.c_str(), h3SepA, h5SepB);
					ssSepType = string(xs);
				}
				else if(h3IDB >=0 && wcPairPosID[h3IDB] == h5IDA){
					sprintf(xs, "nnb-%s-5P3P-%d-%d",pairType.c_str(), h5SepA, h3SepB);
					ssSepType = string(xs);
				}
				else {
					ssSepType = "UNK";
				}
			}
			
			
			if(task == "refinement"){
				//cout << "init near native move set: "  << endl;
				this->allEdges[i*seqLen+j]->initNearNativeMoveSet();

				//if(this->allEdges[i*seqLen+j]->weight < 0) {
				//	double minD = this->allEdges[i*seqLen+j]->moveSet->minDistanceToMoveSet(this->allEdges[i*seqLen+j]->cm, oi);
				//	double maxD = this->allEdges[i*seqLen+j]->moveSet->maxDistanceToMoveSet(this->allEdges[i*seqLen+j]->cm, oi);
				//	cout << "minD: " << minD << " maxD: " << maxD << endl;
				//}
			}
			else if(task == "singleBasePredict"){
				this->allEdges[i*seqLen+j]->initNearNativeMoveSet(2.0);
			}
			else if(task == "predict"){
				if(eiLib == NULL) {
					cout << "edgeInformation lib not initialized" << endl;
					exit(0);
				}
				if(ssSepType != "UNK") {


					this->allEdges[i*seqLen+j]->ei->setToLibPCluster(ssSepType, eiLib);
					//this->allEdges[i*seqLen+j]->moveSet->updateEdgeInformation(this->allEdges[i*seqLen+j]->ei);
					this->allEdges[i*seqLen+j]->samplingFreq = this->allEdges[i*seqLen+j]->ei->validClusterNum*1.0;
					if(this->allEdges[i*seqLen+j]->ei->validClusterNum > 20) {
						this->allEdges[i*seqLen+j]->samplingFreq = 20.0;
					}
					this->allEdges[i*seqLen+j]->weight = this->allEdges[i*seqLen+j]->ei->weight;
					this->allEdges[i*seqLen+j]->weightRand = this->allEdges[i*seqLen+j]->weight;

					cout << "init edge: " << i << "-" << j << " " << ssSepType <<  " clusterNum: " << this->allEdges[i*seqLen+j]->ei->validClusterNum << endl;
					if(ssSepType == "nnb-GG-3P5P-2-3")
						this->allEdges[i*seqLen+j]->ei->print();
				}
			}
			else if(task == "analysis") {
				NuNode* nodeA = this->allEdges[i*seqLen+j]->nodeA;
				NuNode* nodeB = this->allEdges[i*seqLen+j]->nodeB;
				BaseDistanceMatrix dm(nodeA->baseConf->cs1, nodeB->baseConf->cs1);
				int clusterID = pairLib->getPairType(dm, nodeA->baseType, nodeB->baseType, this->sepTable[i*seqLen+j]);
				this->allEdges[i*seqLen+j]->ei->setUniqueCluster(clusterID, pairLib);
				this->allEdges[i*seqLen+j]->weight = this->allEdges[i*seqLen+j]->ei->weight;
				this->allEdges[i*seqLen+j]->weightRand = this->allEdges[i*seqLen+j]->weight;
			}
		}
	}

	/*
	 * init from contact key:
	 */
	cout << "init from key" << endl;
	if(contactKey.length() > 0){
		int n = contactKey.length()/6;
		int a,b,c,d,e,f, sep;
		int indexA, indexB, clusterID, typeA, typeB;
		BasePairLib* bpLib = new BasePairLib("stat");
		cout << "init from contactKey: " << endl;

		for(i=0;i<seqLen;i++){
			for(j=0;j<seqLen;j++){
				this->allEdges[i*seqLen+j]->weightRand = this->allEdges[i*seqLen+j]->weightRand+10.0;
			}
		}

		for(i=0;i<n;i++){
			a = contactKey[i*6] - '!';
			b = contactKey[i*6+1] - '!';
			c = contactKey[i*6+2] - '!';
		 	d = contactKey[i*6+3] - '!';
		 	e = contactKey[i*6+4] - '!';
		 	f = contactKey[i*6+5] - '!';
			indexA = a*90+b;
			indexB = c*90+d;
			clusterID = e*90+f;
			sep = this->sepTable[indexA*seqLen+indexB];
			typeA = this->allNodes[indexA]->baseType;
			typeB = this->allNodes[indexB]->baseType;
			BaseDistanceMatrix dm;
			if(sep == 1)
				dm = bpLib->nbDMClusterCenters[typeA*4+typeB][clusterID];
			else if(sep == -1)
				dm = bpLib->revNbDMClusterCenters[typeA*4+typeB][clusterID];
			else 
				dm = bpLib->nnbDMClusterCenters[typeA*4+typeB][clusterID];

			this->allEdges[indexA*seqLen+indexB]->initMoveSet(dm, 1.6);

			dm = dm.reverse();
			this->allEdges[indexB*seqLen+indexA]->initMoveSet(dm, 1.6);

		}
		delete bpLib;
	}
	
	/*
	 * init from ct 
	 */

	cout << "init from CT: " << ctList.size() << endl;

	for(i=0;i<ctList.size();i++){
        splitString(ctList[i], " ", &spt);
        if(spt.size() != 3 && spt.size() != 4){
            cout << "invalid connect line" << endl;
            cout << ctList[i] << endl;
            exit(0);
    	}
        int idA = atoi(spt[0].c_str());
        int idB = atoi(spt[1].c_str());
        string ctType = spt[2];

		if(ctType == "nat") {
			this->allEdges[idA*seqLen+idB]->weight = -999.9;
			this->allEdges[idB*seqLen+idA]->weight = -999.9;
			this->allEdges[idA*seqLen+idB]->fixNaiveMove();
			this->allEdges[idB*seqLen+idA]->fixNaiveMove();
		}
		else {
			this->allEdges[idA*seqLen+idB]->ei->setToLibPCluster(ctType, eiLib);
			//this->allEdges[idA*seqLen+idB]->moveSet->updateEdgeInformation(this->allEdges[idA*seqLen+idB]->ei);
			this->allEdges[idB*seqLen+idA]->ei->setToLibReversePCluster(ctType, eiLib);
			//this->allEdges[idB*seqLen+idA]->moveSet->updateEdgeInformation(this->allEdges[idB*seqLen+idA]->ei);
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
						cout << "fix move: " << i << " " << k << endl;
						cout << "fix move: " << k << " " << i << endl;
						this->allEdges[i*seqLen+k]->fixNaiveMove();
						this->allEdges[k*seqLen+i]->fixNaiveMove();
					}
				}
			}
		}
	}

	for(i=0;i<seqLen;i++){
		if(csn[i] == 'F')
			this->allNodes[i]->samplingFreq = 0;
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
			et->pb->buildPhosphate(allNodes[i]->riboseConf, allNodes[i+1]->riboseConf, allNodes[i]->phoConf);
			allNodes[i]->phoConfTmp->copyValueFrom(allNodes[i]->phoConf);
		}
	}
}

void NuGraph::initPho(PO3Builder* pb) {
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i]){
			pb->buildPhosphate(allNodes[i]->riboseConf, allNodes[i+1]->riboseConf, allNodes[i]->phoConf);
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
	string csn = input.getValue("csn");
	string key = input.getValue("key");
	vector<string> ctList = input.getMultiValues("ct");

	if(cst.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			cst = cst + "0";
		}
	}

	if(csn.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			csn = csn + "0";
		}
	}

	string cnt = input.getValue("cnt");
	if(cnt.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			cnt = cnt + "-";
		}
	}

	RNAPDB* pdb = new RNAPDB(pdbFile);

	init(task, pdb, baseSeq, baseSec, csn, cst, cnt, key, ctList);
	

	initPho();
	this->initInfo = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, 0.0, atLib, 0);	

	delete pdb;
}

void NuGraph::initForCGMC(const string& inputFile){

	NSPtools::InputParser input(inputFile);
	input.printOptions();

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 */

	string task = "predict";
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string cst = input.getValue("cst");
	string csn = input.getValue("csn");
	string key = input.getValue("key");
	vector<string> ctList = input.getMultiValues("ct");
	
	if(cst.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			cst = cst + "0";
		}
	}

	if(csn.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			csn = csn + "0";
		}
	}

	string cnt = input.getValue("cnt");
	if(cnt.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			cnt = cnt + "-";
		}
	}

	RNAPDB* pdb = new RNAPDB(pdbFile);
	init(task, pdb, baseSeq, baseSec, csn, cst, cnt, key, ctList);
	this->initInfo = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, 0.0, atLib, 1);	
	delete pdb;
}

void NuGraph::initForMST(const string& inputFile){
	NSPtools::InputParser input(inputFile);
	input.printOptions();

	string task = "analysis";
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string cst = input.getValue("cst");
	string csn = input.getValue("csn");
	string key = input.getValue("key");
	vector<string> ctList = input.getMultiValues("ct");
	if(csn.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			csn = csn + "0";
		}
	}

	string cnt = input.getValue("cnt");
	RNAPDB* pdb = new RNAPDB(pdbFile);
	init(task, pdb, baseSeq, baseSec, csn, cst, cnt, key, ctList);
	delete pdb;
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
	string csn = input.getValue("csn");
	vector<string> ctList = input.getMultiValues("ct");
	if(csn.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			csn = csn + "0";
		}
	}

	string key = input.getValue("key");

	string cnt = input.getValue("cnt");

	RNAPDB* pdb = new RNAPDB(pdbFile);
	init(task, pdb, baseSeq, baseSec, csn, cst, cnt, key, ctList);

	initPho();
	this->initInfo = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, 0.0, atLib, 0);	
	delete pdb;

	NuNode* nodeTarget = allNodes[pos];
	double d;
	int i,j,k;

	int seqLen = baseSeq.length();
	char xx[seqLen+1];
	for(i=0;i<seqLen;i++)
		xx[i] = '0';
	xx[seqLen] = '\0';

	for(i=0;i<=pos-3;i++){
		if(allNodes[i]->samplingFreq > 0.0)
			allNodes[i]->samplingFreq = 0.2;
	}
	for(i=pos+3;i<seqLen;i++){
		if(allNodes[i]->samplingFreq > 0.0)
			allNodes[i]->samplingFreq = 0.2;
	}

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

	for(i=0;i<seqLen;i++){
		if(csn[i] == 'F')
			this->allNodes[i]->samplingFreq = 0;
	}

	for(i=0;i<seqLen;i++) {
		char c = cst[i];
		if(c >= 'A' && c <= 'Z') {
			this->fixed[i] = true;
			if(this->allNodes[i]->samplingFreq != 0)
				this->allNodes[i]->samplingFreq = 0.3;
		}
		else {
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
		for(int j=i+1;j<seqLen;j++){
			//cout << "edge: " << i << " " << j << " : " << allEdges[i*seqLen+j]->ei->getRandomWeight() << endl;
			allEdges[i*seqLen+j]->weightRand = allEdges[i*seqLen+j]->ei->getRandomWeight();
			allEdges[j*seqLen+i]->weightRand = allEdges[i*seqLen+j]->weightRand;
		}
	}
}

void NuGraph::initNearestNativeEdge(){
	for(int i=0;i<seqLen;i++){
		for(int j=i+1;j<seqLen;j++){
			allEdges[i*seqLen+j]->initNearNativeMoveSet();
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

void NuGraph::updateEnergy(double clashRescale, double connectRescale){
	for(int i=0;i<seqLen;i++){
		allNodes[i]->updateEnergy(clashRescale, connectRescale);
	}

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			allEdges[i*seqLen+j]->updateEnergy(clashRescale, connectRescale);
		}
	}
}

void NuGraph::updateEnergyCG(double clashRescale, double connectRescale){
	for(int i=0;i<seqLen;i++){
		allNodes[i]->updateEnergyCG(clashRescale, connectRescale);
	}

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			allEdges[i*seqLen+j]->updateEnergyCG(clashRescale, connectRescale);
		}
	}
}

string NuGraph::toContactMapHashKeyCG(){

	string key;
	int i,j, k, sep, clusterID;
	char a,b;
	double ene;

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = this->sepTable[i*seqLen+j];
			double wt = 0.0;

			if(this->allEdges[i*seqLen+j]->pairEneCG[0] < 0){
				BaseDistanceMatrix dm(this->allNodes[i]->baseConfCG->cs1, allNodes[j]->baseConfCG->cs1);
				clusterID = pairLib->getPairType(dm, allNodes[i]->baseType, allNodes[j]->baseType, sep);
				if(clusterID < 0) continue;

				a = i/90 + '!';
				b = i%90 + '!';
				key.push_back(a);
				key.push_back(b);

				a = j/90 + '!';
				b = j%90 + '!';
				key.push_back(a);
				key.push_back(b);

				a = clusterID/90 + '!';
				b = clusterID%90 + '!';
				key.push_back(a);
				key.push_back(b);
			}		
		}
	}
	return key;
}

void NuGraph::keyToContactMatrix(const string& key){
	int n = key.length()/6;
	int a,b,c,d,e,f;
	int indexA, indexB, clusterID;
	for(int i=0;i<n;i++){
		 a = key[i*6] - '!';
		 b = key[i*6+1] - '!';
		 c = key[i*6+2] - '!';
		 d = key[i*6+3] - '!';
		 e = key[i*6+4] - '!';
		 f = key[i*6+5] - '!';
		indexA = a*90+b;
		indexB = c*90+d;
		clusterID = e*90+f;
		cout << indexA << " " << indexB << " " << clusterID << endl;
	}
}

void NuGraph::keyToMatrixFile(const string& key, const string& outfile){
	
	int mtxNat[seqLen][seqLen];
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			LocalFrame csA = this->allNodes[i]->baseConf->cs1;
			LocalFrame csB = this->allNodes[j]->baseConf->cs1;

			BaseDistanceMatrix dm(csA, csB);
			int sep = this->sepTable[i*seqLen+j];
			int clusterID = pairLib->getPairType(dm, allNodes[i]->baseType, allNodes[j]->baseType, sep);
			mtxNat[i][j] = clusterID;
		}
	}
	
	
	int mtx[seqLen][seqLen];
	char compareMtx[seqLen][seqLen];

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			mtx[i][j] = -1;
			compareMtx[i][j] = '-';
		}
	}

	int n = key.length()/6;
	int a,b,c,d,e,f;
	int indexA, indexB, clusterID;
	for(int i=0;i<n;i++){
		 a = key[i*6] - '!';
		 b = key[i*6+1] - '!';
		 c = key[i*6+2] - '!';
		 d = key[i*6+3] - '!';
		 e = key[i*6+4] - '!';
		 f = key[i*6+5] - '!';
		indexA = a*90+b;
		indexB = c*90+d;
		clusterID = e*90+f;
		mtx[indexA][indexB] = clusterID;

		if(mtx[indexA][indexB] != mtxNat[indexA][indexB]){
			compareMtx[indexA][indexB] = 'W';
		}
	}

	ofstream out;
	out.open(outfile.c_str(), ios::out);
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			out << compareMtx[i][j] << " ";
		}
		out << endl;
	}
	out.close();

}

void NuGraph::printContactMatrix(const string& outfile){
	ofstream out;
	out.open(outfile.c_str(), ios::out);

	char xx[20];
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			LocalFrame csA = this->allNodes[i]->baseConf->cs1;
			LocalFrame csB = this->allNodes[j]->baseConf->cs1;

			BaseDistanceMatrix dm(csA, csB);
			int sep = this->sepTable[i*seqLen+j];
			int clusterID = pairLib->getPairType(dm, allNodes[i]->baseType, allNodes[j]->baseType, sep);
			sprintf(xx, "%-4d ", clusterID);
			out << string(xx);
		}
		out << endl;
	}
	out.close();
}

double NuGraph::keyAccuracy(const string& key){
	//NuGraph init for analysis
	int n = key.length()/6;
	int a,b,c,d,e,f;
	int indexA, indexB, clusterID;
	double matchedNum = 0;
	double totalNum = seqLen - 1.0;

	int fixedNum = 0;

	for(int i=0;i<n;i++){
		 a = key[i*6] - '!';
		 b = key[i*6+1] - '!';
		 c = key[i*6+2] - '!';
		 d = key[i*6+3] - '!';
		 e = key[i*6+4] - '!';
		 f = key[i*6+5] - '!';
		indexA = a*90+b;
		indexB = c*90+d;
		clusterID = e*90+f;
		//cout << key.substr(i*6, 6) << " " << indexA << " " << indexB << " " << clusterID << " " << e << " " << f << endl;
		
		LocalFrame csA = this->allNodes[indexA]->baseConf->cs1;
		LocalFrame csB = this->allNodes[indexB]->baseConf->cs1;
		BaseDistanceMatrix dm(csA, csB);
		int sep = this->sepTable[indexA*seqLen+indexB];
		if(this->allEdges[indexA*seqLen+indexB]->ei->fixed)
			fixedNum++;

		int clusterIDNat = pairLib->getPairType(dm, allNodes[indexA]->baseType, allNodes[indexB]->baseType, sep);

		BaseDistanceMatrix dmNat = dm;
		BaseDistanceMatrix dmPred;

		if(sep == 1) {
			dmPred = pairLib->nbDMClusterCenters[allNodes[indexA]->baseType*4+allNodes[indexB]->baseType][clusterID];
		}
		else if(sep == -1) {
			dmPred = pairLib->revNbDMClusterCenters[allNodes[indexA]->baseType*4+allNodes[indexB]->baseType][clusterID];
		}
		else {
			dmPred = pairLib->nnbDMClusterCenters[allNodes[indexA]->baseType*4+allNodes[indexB]->baseType][clusterID];
		}

		double dist = dmNat.distanceTo(dmPred);
		//printf("base pair: %2d %2d clusterNat: %3d clusterPred: %3d distance: %6.3f\n", indexA, indexB, clusterIDNat, clusterID, dist);

		if(dist < 1.0) {
			matchedNum += 1.0;
		}
		else if(dist < 1.5){
			matchedNum += (1.5 - dist)*2;
		}
	}
	return (matchedNum-fixedNum)*1.0/(totalNum-fixedNum);
}

void NuGraph::checkEnergy(double clashRescale, double connectRescale){
	for(int i=0;i<seqLen;i++){
		allNodes[i]->checkEnergy(clashRescale, connectRescale);
	}
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			if(i==j) continue;
			allEdges[i*seqLen+j]->checkEnergy(clashRescale, connectRescale);
			allEdges[i*seqLen+j]->checkReversePair();
		}
	}
}

void NuGraph::checkEnergyCG(double clashRescale, double connectRescale){
	for(int i=0;i<seqLen;i++){
		allNodes[i]->checkEnergyCG(clashRescale, connectRescale);
	}
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			if(i==j) continue;
			allEdges[i*seqLen+j]->checkEnergyCG(clashRescale, connectRescale);
			allEdges[i*seqLen+j]->checkReversePairCG();
		}
	}
}

double NuGraph::totalEnergy(double clashRescale, double connectRescale){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){

		ene += allNodes[i]->riboseConf->rot->energy;
		ene += allNodes[i]->phoConf->ene * connectRescale;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergy(allNodes[i]->baseConf, allNodes[j]->baseConf, sep, et, clashRescale);
			ene += nuBaseRiboseEnergy(allNodes[i]->baseConf, allNodes[j]->riboseConf, sep, et, clashRescale);
			ene += nuBaseRiboseEnergy(allNodes[j]->baseConf, allNodes[i]->riboseConf, sepR, et, clashRescale);
			ene += nuBasePhoEnergy(allNodes[i]->baseConf, allNodes[j]->phoConf, sep, et, clashRescale);
			ene += nuBasePhoEnergy(allNodes[j]->baseConf, allNodes[i]->phoConf, sepR, et, clashRescale);
			ene += nuRiboseRiboseEnergy(allNodes[i]->riboseConf, allNodes[j]->riboseConf, sep, et, clashRescale);
			ene += nuRibosePhoEnergy(allNodes[i]->riboseConf, allNodes[j]->phoConf, sep, et, clashRescale);
			ene += nuRibosePhoEnergy(allNodes[j]->riboseConf, allNodes[i]->phoConf, sepR, et, clashRescale);
			ene += nuPhoPhoEnergy(allNodes[i]->phoConf, allNodes[j]->phoConf, sep, et, clashRescale);
		}
	}
	return ene;
}

double NuGraph::nbEnergy(double clashRescale, double connectRescale) {
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){

		ene += allNodes[i]->riboseConf->rot->energy;
		ene += allNodes[i]->phoConf->ene * clashRescale;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			if(sep == 1)
				ene += nuBaseBaseEnergy(allNodes[i]->baseConf, allNodes[j]->baseConf, sep, et, clashRescale);
		}
	}
	return ene;	
}

double NuGraph::nnbEnergy(double clashRescale, double connectRescale) {
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){

		ene += allNodes[i]->riboseConf->rot->energy;
		ene += allNodes[i]->phoConf->ene* clashRescale;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			if(sep > 1)
				ene += nuBaseBaseEnergy(allNodes[i]->baseConf, allNodes[j]->baseConf, sep, et, clashRescale);
		}
	}
	return ene;	
}

double NuGraph::totalEnergyCG(double clashRescale, double connectRescale){
		double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		ene += allNodes[i]->riboseConfCG->rot->energy;
		ene += allNodes[i]->bbcg;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergyCG(allNodes[i]->baseConfCG, allNodes[j]->baseConfCG, sep, et, clashRescale);
			ene += nuBaseRiboseEnergyCG(allNodes[i]->baseConfCG, allNodes[j]->riboseConfCG, sep, et, clashRescale);
			ene += nuBaseRiboseEnergyCG(allNodes[j]->baseConfCG, allNodes[i]->riboseConfCG, sepR, et, clashRescale);
			ene += nuRiboseRiboseEnergyCG(allNodes[i]->riboseConfCG, allNodes[j]->riboseConfCG, sep, et, clashRescale);
		}
	}
	return ene;
}

double NuGraph::totalEnergyCGTmp(double clashRescale, double connectRescale){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		ene += allNodes[i]->riboseConfCGTmp->rot->energy;
		ene += allNodes[i]->bbcgTmp;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergyCG(allNodes[i]->baseConfCGTmp, allNodes[j]->baseConfCGTmp, sep, et, clashRescale);
			ene += nuBaseRiboseEnergyCG(allNodes[i]->baseConfCGTmp, allNodes[j]->riboseConfCGTmp, sep, et, clashRescale);
			ene += nuBaseRiboseEnergyCG(allNodes[j]->baseConfCGTmp, allNodes[i]->riboseConfCGTmp, sepR, et, clashRescale);
			ene += nuRiboseRiboseEnergyCG(allNodes[i]->riboseConfCGTmp, allNodes[j]->riboseConfCGTmp, sep, et, clashRescale);
		}
	}
	return ene;
}

double NuGraph::totalEnergyTmp(double clashRescale, double connectRescale){
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
			ene += nuBaseBaseEnergy(allNodes[i]->baseConfTmp, allNodes[j]->baseConfTmp, sep, et, clashRescale);
			ene += nuBaseRiboseEnergy(allNodes[i]->baseConfTmp, allNodes[j]->riboseConfTmp, sep, et, clashRescale);
			ene += nuBaseRiboseEnergy(allNodes[j]->baseConfTmp, allNodes[i]->riboseConfTmp, sepR, et, clashRescale);
			ene += nuBasePhoEnergy(allNodes[i]->baseConfTmp, allNodes[j]->phoConfTmp, sep, et, clashRescale);
			ene += nuBasePhoEnergy(allNodes[j]->baseConfTmp, allNodes[i]->phoConfTmp, sepR, et, clashRescale);
			ene += nuRiboseRiboseEnergy(allNodes[i]->riboseConfTmp, allNodes[j]->riboseConfTmp, sep, et, clashRescale);
			ene += nuRibosePhoEnergy(allNodes[i]->riboseConfTmp, allNodes[j]->phoConfTmp, sep, et, clashRescale);
			ene += nuRibosePhoEnergy(allNodes[j]->riboseConfTmp, allNodes[i]->phoConfTmp, sepR, et, clashRescale);
			ene += nuPhoPhoEnergy(allNodes[i]->phoConfTmp, allNodes[j]->phoConfTmp, sep, et, clashRescale);
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
		for(int j=i+1;j<seqLen;j++){
			NuEdge* eg = this->allEdges[i*seqLen+j];
			printf("edge: %2d %2d sep: %2d BB %7.3f BR %7.3f BP %7.3f RB %7.3f RR %7.3f RP %7.3f PB %7.3f PR %7.3f PP %7.3f\n", i, j, eg->sep, eg->pairEne[0], eg->pairEne[1], eg->pairEne[2],eg->pairEne[3],eg->pairEne[4],eg->pairEne[5],eg->pairEne[6],eg->pairEne[7],eg->pairEne[8]);
			//printf("eTmp: %2d %2d sep: %2d BB %7.3f BR %7.3f BP %7.3f RB %7.3f RR %7.3f RP %7.3f PB %7.3f PR %7.3f PP %7.3f\n", i, j, eg->sep, eg->pairEneTmp[0], eg->pairEneTmp[1], eg->pairEneTmp[2],eg->pairEneTmp[3],eg->pairEneTmp[4],eg->pairEneTmp[5],eg->pairEneTmp[6],eg->pairEneTmp[7],eg->pairEneTmp[8]);
		}
	}
}

void NuGraph::printEnergyCG(double clashRescale){

	double delta = 0.0;
	int sep, sepR;

	for(int i=0;i<seqLen;i++){
		//if(allNodes[i]->eneCG == allNodes[i]->eneCGTmp) continue;
		printf("node: %2d ene: %7.3f eneTmp: %7.3f\n", i, allNodes[i]->eneCG, allNodes[i]->eneCGTmp);
		delta += allNodes[i]->eneCGTmp - allNodes[i]->eneCG;
		printf("NODE: %2d ene: %7.3f eneTmp: %7.3f\n", i, allNodes[i]->riboseConfCG->rot->energy+allNodes[i]->bbcg, allNodes[i]->riboseConfCGTmp->rot->energy+allNodes[i]->bbcgTmp);
	}
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			NuEdge* eg = this->allEdges[i*seqLen+j];
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];

			if(eg->pairEneCG[0] == eg->pairEneCGTmp[0] && eg->pairEneCG[1] == eg->pairEneCGTmp[1] && eg->pairEneCG[2] == eg->pairEneCGTmp[2] && eg->pairEneCG[3] == eg->pairEneCGTmp[3]) continue;
			
			if(i < j) {
				delta += eg->pairEneCGTmp[0] - eg->pairEneCG[0];
				delta += eg->pairEneCGTmp[1] - eg->pairEneCG[1];
				delta += eg->pairEneCGTmp[2] - eg->pairEneCG[2];
				delta += eg->pairEneCGTmp[3] - eg->pairEneCG[3];

			double e1 = nuBaseBaseEnergyCG(allNodes[i]->baseConfCG, allNodes[j]->baseConfCG, sep, et, clashRescale);
			double e2 = nuBaseRiboseEnergyCG(allNodes[i]->baseConfCG, allNodes[j]->riboseConfCG, sep, et, clashRescale);
			double e3 = nuBaseRiboseEnergyCG(allNodes[j]->baseConfCG, allNodes[i]->riboseConfCG, sepR, et, clashRescale);
			double e4 = nuRiboseRiboseEnergyCG(allNodes[i]->riboseConfCG, allNodes[j]->riboseConfCG, sep, et, clashRescale);
			
			double f1 = nuBaseBaseEnergyCG(allNodes[i]->baseConfCGTmp, allNodes[j]->baseConfCGTmp, sep, et, clashRescale);
			double f2 = nuBaseRiboseEnergyCG(allNodes[i]->baseConfCGTmp, allNodes[j]->riboseConfCGTmp, sep, et, clashRescale);
			double f3 = nuBaseRiboseEnergyCG(allNodes[j]->baseConfCGTmp, allNodes[i]->riboseConfCGTmp, sepR, et, clashRescale);
			double f4 = nuRiboseRiboseEnergyCG(allNodes[i]->riboseConfCGTmp, allNodes[j]->riboseConfCGTmp, sep, et, clashRescale);


				printf("edge: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", i, j, eg->sep, eg->pairEneCG[0], eg->pairEneCG[1], eg->pairEneCG[2],eg->pairEneCG[3]);
				printf("eTmp: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", i, j, eg->sep, eg->pairEneCGTmp[0], eg->pairEneCGTmp[1], eg->pairEneCGTmp[2],eg->pairEneCGTmp[3]);
				printf("EDGE: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", i, j, eg->sep, e1, e2, e3, e4);
				printf("ETMP: %2d %2d sep: %2d BB %7.3f BR %7.3f RB %7.3f RR %7.3f\n", i, j, eg->sep, f1, f2, f3, f4);
			}

		}
	}
	printf("delta: %7.3f\n", delta);

}

void NuGraph::cgToAllAtom(){

	for(int i=0;i<seqLen;i++){
		this->allNodes[i]->baseConf->updateCoords(allNodes[i]->baseConfCG->cs1);
		RiboseRotamer* rot = rotLib->riboseRotLib->rotCGToRot(allNodes[i]->riboseConfCG->rot);
		this->allNodes[i]->riboseConf->updateLocalFrameAndRotamer(allNodes[i]->baseConfCG->cs1, rot);
		if(connectToDownstream[i]) {
			et->pb->buildPhosphate(allNodes[i]->riboseConf, allNodes[i+1]->riboseConf, allNodes[i]->phoConf);
		}
	}
}

double NuGraph::getTotalEnergyForModelSelection(){
	double ene = 0.0;
	int i,j,k, sep, sepR;
	for(i=0;i<seqLen;i++){
		ene += allNodes[i]->phoConf->ene * 0.7;
	}

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene += nuBaseBaseEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->baseConf, sep, et);
			ene += nuBaseRiboseEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuBaseRiboseEnergyForSelection(allNodes[j]->baseConf, allNodes[i]->riboseConf, sepR, et);
			ene += nuBasePhoEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->phoConf, sep, et);
			ene += nuBasePhoEnergyForSelection(allNodes[j]->baseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuRiboseRiboseEnergyForSelection(allNodes[i]->riboseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuRibosePhoEnergyForSelection(allNodes[i]->riboseConf, allNodes[j]->phoConf, sep, et);
			ene += nuRibosePhoEnergyForSelection(allNodes[j]->riboseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuPhoPhoEnergyForSelection(allNodes[i]->phoConf, allNodes[j]->phoConf, sep, et);
		}
	}

	BaseSasaPoints* bs = new BaseSasaPoints();
	vector<NuNode*> nodes;
	for(i=0;i<seqLen;i++){
		nodes.push_back(allNodes[i]);
	}
	AssignRNASasa* as = new AssignRNASasa(bs, nodes);
	for(i=0;i<this->seqLen;i++){
		int expNum = as->getExposeNum(i);
		ene += -0.12*sqrt(1.0*expNum);
	}

	return ene;	
}

void NuGraph::printBaseEnergyList(const string& outfile){
	ofstream out;
	out.open(outfile.c_str(), ios::out | ios::app);
	int i,j,sep, sepR;

	double* eneList = new double[seqLen];
	for(i=0;i<seqLen;i++){
		eneList[i] = 0.0;
	}

	double ene;

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene = nuBaseBaseEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->baseConf, sep, et);
			ene += nuBaseRiboseEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuBaseRiboseEnergyForSelection(allNodes[j]->baseConf, allNodes[i]->riboseConf, sepR, et);
			ene += nuBasePhoEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->phoConf, sep, et);
			ene += nuBasePhoEnergyForSelection(allNodes[j]->baseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuRiboseRiboseEnergyForSelection(allNodes[i]->riboseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuRibosePhoEnergyForSelection(allNodes[i]->riboseConf, allNodes[j]->phoConf, sep, et);
			ene += nuRibosePhoEnergyForSelection(allNodes[j]->riboseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuPhoPhoEnergyForSelection(allNodes[i]->phoConf, allNodes[j]->phoConf, sep, et);
			
			eneList[i] += 0.5*ene;
			eneList[j] += 0.5*ene;

			cout << "pos: " << i << " " << j << " " << ene << endl;
		}
	}

	string augc = "AUGCatgc";
	char xx[200];
	vector<double> eList;
	for(i=0;i<seqLen;i++){
		sprintf(xx, "eBase: %-3d %c %7.3f", i, augc[this->seq[i]], eneList[i]);
		out << string(xx) << endl;
	}
	out.close();
}

void NuGraph::printPairwiseEnergy(const string& outfile){
	ofstream out;
	out.open(outfile.c_str(), ios::out | ios::app);
	int i,j,sep, sepR;

	double* eneMatrix = new double[seqLen*seqLen];
	for(i=0;i<seqLen*seqLen;i++){
		eneMatrix[i] = 0.0;
	}

	double ene;

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			sep = sepTable[i*seqLen+j];
			sepR = sepTable[j*seqLen+i];
			ene = nuBaseBaseEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->baseConf, sep, et);
			ene += nuBaseRiboseEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuBaseRiboseEnergyForSelection(allNodes[j]->baseConf, allNodes[i]->riboseConf, sepR, et);
			ene += nuBasePhoEnergyForSelection(allNodes[i]->baseConf, allNodes[j]->phoConf, sep, et);
			ene += nuBasePhoEnergyForSelection(allNodes[j]->baseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuRiboseRiboseEnergyForSelection(allNodes[i]->riboseConf, allNodes[j]->riboseConf, sep, et);
			ene += nuRibosePhoEnergyForSelection(allNodes[i]->riboseConf, allNodes[j]->phoConf, sep, et);
			ene += nuRibosePhoEnergyForSelection(allNodes[j]->riboseConf, allNodes[i]->phoConf, sepR, et);
			ene += nuPhoPhoEnergyForSelection(allNodes[i]->phoConf, allNodes[j]->phoConf, sep, et);
			
			eneMatrix[i*seqLen+j] = ene;
			eneMatrix[j*seqLen+i] = ene;
		}
	}

	string augc = "AUGCatgc";
	char xx[200];

	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			if(abs(eneMatrix[i*seqLen+j]) < 0.01) continue;
			sprintf(xx, "ePair: %-4d %c %-4d %c %7.3f", i, augc[this->seq[i]], j, augc[this->seq[j]], eneMatrix[i*seqLen+j]);
			out << string(xx) << endl;
		}
		
	}
	out.close();	
}

void NuGraph::generateSubGraph(const string& inputFile, int corePos, int* subGraphPosList, int* fixedPositions, NuGraph* subGraph, vector<int>& outsubGraphPosList, vector<vector<int>>& Bclusters, vector<vector<int>>& Cclusters) {
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
	string csn = input.getValue("csn");
	string key = input.getValue("key");
	vector<string> ctList = input.getMultiValues("ct");

	if(cst.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			cst = cst + "0";
		}
	}

	if(csn.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			csn = csn + "0";
		}
	}

	string cnt = input.getValue("cnt");
	if(cnt.length() == 0){
		for(int i=0;i<baseSeq.length();i++){
			cnt = cnt + "-";
		}
	}
	cout << "corePos:seq[" << corePos << "]:" << seq[corePos] << endl;
	int subGraphPosCount = 0;  // 用于跟踪 subGraphPosList 的索引
	int fixedPosCount = 0;  // 用于跟踪 fixedPositions 的索引
	// 第一步：遍历所有序列节点，找到与核心节点平方距离小于 900.0 的节点
	for (int pos = 0; pos < seqLen; pos++) {	
		if (squareDistance(allNodes[pos]->baseConf->coords[0], allNodes[corePos]->baseConf->coords[0]) < 400.0) {
			subGraphPosList[subGraphPosCount] = pos;
			subGraphPosCount++;
		}	
	}
	cout << "A Array:" << endl;
	for (int i = 0; i < subGraphPosCount; i++) {
		cout << subGraphPosList[i] << " ";
	}
	cout << endl;

	for (int idx = 0; idx < subGraphPosCount; idx++) {
		int subNode = subGraphPosList[idx];
		// 遍历所有 seqLen 中的节点，寻找与 subNode 最近原子平方距离小于 64 的节点
		for (int neighbor = 0; neighbor < seqLen; neighbor++) {
			if (neighbor != subNode) {
				// 找到最小的平方距离
				double minSDistance = minSquareDistance(neighbor, subNode); 

				if (minSDistance < 64.0) {
					bool alreadyInSubGraph = false;

					// 检查节点 neighbor 是否已经在 subGraphPosList 中
					for (int i = 0; i < subGraphPosCount; i++) {
						if (subGraphPosList[i] == neighbor) {
							alreadyInSubGraph = true;
							break;
						}
					}

					// 如果节点 neighbor 不在 subGraphPosList 中，则将其加入 fixedPositions
					if (!alreadyInSubGraph) {
						fixedPositions[fixedPosCount] = neighbor;
						fixedPosCount++;
						// cout << "neighbor:" << neighbor << "->subNode:" << subNode << endl;
					}
				}
			}
		}
	}

	// 使用 set 来去除重复的元素
    set<int> uniquePositions(fixedPositions, fixedPositions + fixedPosCount);
    // 将 set 中的元素拷贝回原始数组
    int newfixedPosCount = 0;
    for (const int& pos : uniquePositions) {
        fixedPositions[newfixedPosCount] = pos;
        newfixedPosCount++;
    }
    // 更新 fixedPosCount 为新数组的大小
    fixedPosCount = newfixedPosCount;

	cout << "B Array: " << endl;
	for (int i = 0; i < fixedPosCount; i++) {
		cout << fixedPositions[i] << " ";
	}
	cout << endl;

	// 第三步：将 fixedPositions 中的元素加入到 subGraphPosList 中, 并对 subGraphPosList 中的元素重新排序
	for (int i = 0; i < fixedPosCount; i++) {
		subGraphPosList[subGraphPosCount] = fixedPositions[i];
		subGraphPosCount++;
	}
    sort(subGraphPosList, subGraphPosList + subGraphPosCount); // 使用 sort 进行排序

	cout << "A+B Array: " << endl;
	for (int i = 0; i < subGraphPosCount; i++) {
		cout << subGraphPosList[i] << " ";
	}
	cout << endl;

	// 集合C：所有不在 subGraphPosList 中的节点
	// 使用 unordered_set 存储 subGraphPosList 中的元素以提高查找效率
    unordered_set<int> subGraphSet(subGraphPosList, subGraphPosList + subGraphPosCount);
    // 计算集合C的大小
    int CSize = 0;
    for (int i = 0; i < seqLen; i++) {
        if (subGraphSet.find(i) == subGraphSet.end()) {
            CSize++;
        }
    }
    // 创建新的 int 数组
    int* CArray = new int[CSize];
    int Cindex = 0;
    // 填充集合C
    for (int i = 0; i < seqLen; i++) {
        if (subGraphSet.find(i) == subGraphSet.end()) {
            CArray[Cindex++] = i;
        }
    }
    // 输出集合C
    cout << "C Array: ";
    for (int i = 0; i < CSize; i++) {
        cout << CArray[i] << " ";
    }
    cout << endl;

	// 集合B+C 判断连通性
	int BCSize = fixedPosCount + CSize;
	int* BCArray = new int[BCSize];
	for (int i = 0; i < fixedPosCount; i++) {
		BCArray[i] = fixedPositions[i];	
	}
	for (int i = fixedPosCount; i < BCSize; i++) {
		BCArray[i] = CArray[i - fixedPosCount];
	}
    sort(BCArray, BCArray + BCSize); // 使用 sort 进行排序

	// 输出集合B+C
    cout << "B+C Array: ";
    for (int i = 0; i < BCSize; i++) {
        cout << BCArray[i] << " ";
    }
    cout << endl;

	unordered_set<int> BSet(fixedPositions, fixedPositions + fixedPosCount);
    unordered_set<int> CSet(CArray, CArray + CSize);

	// 对 BCArray 中的元素按联通性进行聚类（局部变量 clusters）
	vector<vector<int>> clusters;
	clusterContacts(BCArray, BCSize, clusters);
	int clusterIndex = 1;
	for (const auto& cluster : clusters) {
		cout << "Cluster " << clusterIndex << ".size: " << cluster.size() << endl;
		cout << "Cluster " << clusterIndex << ": ";
		vector<int> Bcluster, Ccluster;
        for (int node : cluster) {
            cout << node << " ";
            if (BSet.count(node)) {
                Bcluster.push_back(node);
            } else if (CSet.count(node)) {
                Ccluster.push_back(node);
            }
        }
        Bclusters.push_back(Bcluster);
        Cclusters.push_back(Ccluster);
		cout << endl;
		clusterIndex++;
	}

	//cst
	// 初始化 subcst，默认值为 '0'
	string subcst(subGraphPosCount, '0');
	// 标记每个聚类，使用字母 'A', 'B', 'C' 等
	char clusterLabel = 'A';
	for (const auto& cluster : clusters) {
		// 遍历当前聚类中的每个 fixedPositions 元素(B+C中的元素)
		for (int fixedPos : cluster) {
			// 查找 fixedPos 在 subGraphPosList 中的索引
			for (int i = 0; i < subGraphPosCount; i++) {
				if (subGraphPosList[i] == fixedPos) {
					// 将 subcnt 中相应位置标记为当前的字母
					subcst[i] = clusterLabel;
				}
			}
		}
		// 每处理完一个聚类后，更新到下一个字母
		clusterLabel++;
		// 如果字母超过 'Z'，可以视情况重用字母或处理异常
		if (clusterLabel > 'Z') {
			cout << "Error: Too many clusters to label with single letters." << endl;
			break;
		}
	}

	// 释放内存
    delete[] CArray, BCArray;

    // seq
	string augc = "AUGCatgc";
	string subseq = "";
	for (int i = 0; i < subGraphPosCount; i++) {
        subseq = subseq + augc.substr(this->seq[subGraphPosList[i]],1);  // 复制元素
    }
	string _subseq = "";	
	for (int i = 0; i < subGraphPosCount; i++) {
        _subseq = _subseq + baseSeq[subGraphPosList[i]];  // 复制元素
    }
	// cout << "_Subseq: " << _subseq << endl;

	if (subseq.compare(_subseq) == 0){
		// cout << "Subseq is equal to _Subseq" << endl;
	}else{
		cout << "Subseq is not equal to _Subseq" << endl;
		exit(0);
	}

	// sec
	string sec(this->seqLen, '.');
	sec = indexToBractString(this->wcPairPosID, this->seqLen);// 由wcPairPosID转化为sec
    // cout << "sec: ";
    // for (int i = 0; i < seqLen; i++) {
    //     cout << sec[i];
    // }
    // cout << endl;
	//检查 wcPairPosID 然后转回sec 是否与输入的sec一致；
	if (sec.compare(baseSec) == 0){
		// cout << "sec is equal to baseSec" << endl;
	}else{
		cout << "sec is not equal to baseSec" << endl;
		exit(0);
	}
	string _subsec(subGraphPosCount, '.');// 直接复制而得的sec
    for (int i = 0; i < subGraphPosCount; i++) {
        _subsec[i] = sec[subGraphPosList[i]];  
    }
	cout << "_Subsec: " << _subsec << endl;

	// 检查 subGraphPosList 中的配对关系，重新生成sec
	int* wcPairPosIDCopy = new int[seqLen];
	// cout << "wcPairPosIDCopy: " ;
	for (int i = 0; i < seqLen; i++) {
		wcPairPosIDCopy[i] = this->wcPairPosID[i];
		// cout << wcPairPosIDCopy[i] << " ";
	}
	// cout << endl;
	int* subwcPairPosID = new int[subGraphPosCount];
	cout << "BsubwcPairPosID: " ;
	for (int i = 0; i < subGraphPosCount; i++) {
        subwcPairPosID[i] = wcPairPosIDCopy[subGraphPosList[i]];  // 复制元素
		cout << subwcPairPosID[i] << " ";
    }
	cout << endl;

    // 建立位置到索引的映射，以便快速查找原始配对位置是否在 subGraphPosList 中 如87位与71位配对 则创建映射 wcPairPosIDCopy[index] 87-> index 71
    unordered_map<int, int> posToIndex;
    for (int i = 0; i < subGraphPosCount; i++) {
        int index = wcPairPosIDCopy[subGraphPosList[i]];
        if(index != -1 && wcPairPosIDCopy[index] != -1) {// 不记录不配对位置        
            posToIndex[wcPairPosIDCopy[index]] = index;
            // cout << "posToIndex[" << wcPairPosIDCopy[index] << "] = " << index << endl;
        }
    }
    unordered_map<int, int> posToIndex_2;// 记录配对位置的映射 如 87位与71位配对 87位在subwcPairPosID中的28 71位在subwcPairPosID中的48 则创建映射 28 -> 48
    // 遍历 subwcPairPosID，调整配对位置
    for (int i = 0; i < subGraphPosCount; i++) {
        int originalPairPos = subwcPairPosID[i];
        if (originalPairPos == -1) {
            continue;  // 没有配对，跳过
        }
        // 查找配对位置在 subGraphPosList 中的索引
        if (posToIndex.find(originalPairPos) != posToIndex.end()){
            int pairedIndex = posToIndex[originalPairPos];
            // cout << "originalPairPos = " << originalPairPos << ", pairedIndex = " << pairedIndex << endl;
            for (int j = 0; j < subGraphPosCount; j++) {
                if (subwcPairPosID[j] == pairedIndex) {
                    posToIndex_2[i] = j;// 记录配对位置
                    // cout << "posToIndex_2 = [" << i << "] = " << j << endl;
                    break;
                }
            }
        } else {
            // 若配对位置不在 subGraphPosList 中，将其设置为 -1
            subwcPairPosID[i] = -1;
        }
    }
    // 重新给 subwcPairPosID 编号
    for (const auto& pair : posToIndex_2) {
        int keyIndex = pair.first;
        int valueIndex = pair.second;
        if (keyIndex < valueIndex && keyIndex < subGraphPosCount){
            subwcPairPosID[keyIndex] = valueIndex;       
            subwcPairPosID[valueIndex] = keyIndex;
        }
    }
    // 输出更新后的 subwcPairPosID
    cout << "UsubwcPairPosID: ";
    for (int i = 0; i < subGraphPosCount; i++) {
        cout << subwcPairPosID[i] << " ";
    }
    cout << endl;

	string subsec(subGraphPosCount, '.');
	subsec = indexToBractString(subwcPairPosID, subGraphPosCount);
	cout << "Subsec: " << subsec << endl;
	// 释放内存
    delete[] subwcPairPosID,wcPairPosIDCopy;

	// csn
	for (int i = 0; i < fixedPosCount; i++) {
		csn[fixedPositions[i]] = 'F';
	}
	string subcsn(subGraphPosCount, '0');
    for (int i = 0; i < subGraphPosCount; i++) {
        subcsn[i] = csn[subGraphPosList[i]];  // 复制元素
    }

	//cnt
	string _cnt(this->seqLen, '-');
	for (int i = 0; i < seqLen; i++) {
		if (this->connectToDownstream[i]) {
			_cnt[i] = '-';  // 如果是 true，则用 '-'
		} else {
			_cnt[i] = '|';  // 否则用 '.'
		}
	}
	string subcnt = "";
	for (size_t i = 0; i < subGraphPosCount; i++) {
		int currentPos = subGraphPosList[i];
		subcnt += _cnt[currentPos];  // 将 subGraphPosList 中的字符加入 subcnt

		// 检查是否遇到不连续片段的末端
		if (i < subGraphPosCount - 1 && subGraphPosList[i + 1] != currentPos + 1) {
			subcnt.back() = '|';  // 替换当前段末端的字符为 '|'
		}
	}
	// 最后一个不连续片段的末端字符设为 '|'
	subcnt.back() = '|';

	cout << "Subcst: " << subcst << endl;
	cout << "Subseq: " << subseq << endl;
	cout << "Subcsn: " << subcsn << endl;
	cout << "Subcnt: " << subcnt << endl;

	//pdb
 	// 调用函数，提取符合条件的 ATOM 行并生成新文件

	RNAPDB* subpdb;
	vector<NuNode*> subNodeList;
	for(size_t i=0;i<subGraphPosCount;i++){
		subNodeList.push_back(this->allNodes[subGraphPosList[i]]);
	}
	nodeListToPDBWithoutPho(subNodeList, subpdb);

    //extractPDBAtoms(pdbFile, "output_subGraph.pdb", subGraphPosList, subGraphPosCount);
    //cout << "New PDB file generated: output_subGraph.pdb" << endl;

	subGraph->init(task, subpdb, subseq, subsec, subcsn, subcst, subcnt, key, ctList);
	
	subGraph->initPho();
	subGraph->initInfo = new graphInfo(subGraph->seqLen, subGraph->seq, subGraph->connectToDownstream, subGraph->fixed, subGraph->allNodes, 0.0, atLib, 0);

    // 将元素复制到 vector
    for (int i = 0; i < subGraphPosCount; i++) {
        outsubGraphPosList.push_back(subGraphPosList[i]);  // 将元素添加到 vector 中
    }
}

void NuGraph::nodeListToPDBWithoutPho(vector<NuNode*> nodeList, RNAPDB* outpdb){
	RNAChain* rc;
	string s = "AUGCatgc";
	char ss[20];
	int seqID = 0;

	for(int i=0;i<nodeList.size();i++) {
		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[nodeList[i]->baseType]);
		vector<Atom*> aList = nodeList[i]->toAtomList(this->atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc->addBase(base);
	}
	outpdb->chains.clear();
	outpdb->addChain(rc);
	
}

// 函数：从 PDB 文件中提取符合 subGraphPosList 的 ATOM 行，并生成新的 PDB 文件
void NuGraph::extractPDBAtoms(const string& inputFile, const string& outputFile, int* subGraphPosList, int subGraphPosCount) {
    ifstream infile(inputFile);
    ofstream outfile(outputFile);
    string line;

    // 数组用于保存旧编号与新编号的映射
    int oldNumbers[1000];  // 假设最多1000个不同的旧编号
    int newNumbers[1000];  // 对应的新编号
    int mapSize = 0;       // 当前映射表的大小
    int nextNucleotideNumber = 1;  // 新编号从1开始

    // 原子序号计数器
    int atomSerialNumber = 1;

    // 检查输入文件是否成功打开
    if (!infile.is_open()) {
        cerr << "Error: Unable to open input file " << inputFile << endl;
        return;
    }
    // 检查输出文件是否成功打开
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file " << outputFile << endl;
        return;
    }

    // 遍历 PDB 文件的每一行
    while (getline(infile, line)) {
        if (line.substr(0, 4) == "ATOM") {
            // 获取核苷酸编号 (如 "69")，一般在行的第23到26列
            int nucleotideNumber = stoi(line.substr(22, 4));// 核苷酸编号是从1开始的 而subGraphPosList是从0开始的 所以判断条件需要给subGraphPosList[i]加1

            // 检查该编号是否在 subGraphPosList 中
            bool foundInSubGraph = false;
            for (int i = 0; i < subGraphPosCount; i++) {
                if (nucleotideNumber == subGraphPosList[i] + 1) {
                    foundInSubGraph = true;
                    break;
                }
            }

            if (foundInSubGraph) {
                // 查找旧编号是否已经有对应的新编号
                int newNucleotideNumber = -1;
                for (int i = 0; i < mapSize; i++) {
                    if (oldNumbers[i] == nucleotideNumber) {
                        newNucleotideNumber = newNumbers[i];
                        break;
                    }
                }

                // 如果没有找到，则分配一个新的编号
                if (newNucleotideNumber == -1) {
                    newNucleotideNumber = nextNucleotideNumber;
                    oldNumbers[mapSize] = nucleotideNumber;
                    newNumbers[mapSize] = newNucleotideNumber;
                    mapSize++;
                    nextNucleotideNumber++;
                }

                // 创建新的行，替换核苷酸编号和原子序号
                string formattedLine = "ATOM  ";  // ATOM字段，前6列
                string atomSerialStr = to_string(atomSerialNumber);
                while (atomSerialStr.length() < 5) {
                    atomSerialStr = " " + atomSerialStr;  // 补齐原子序号长度到5位
                }
                formattedLine += atomSerialStr + line.substr(11, 11); // 填充原子序号和其它字段

                // 将新编号转为字符串，手动补齐到4位
                string newNumStr = to_string(newNucleotideNumber);
                while (newNumStr.length() < 4) {
                    newNumStr = " " + newNumStr;  // 前面补空格
                }

                // 拼接新的编号和剩余的行内容
                formattedLine += newNumStr + line.substr(26);

                // 输出到文件
                outfile << formattedLine << endl;

                // 更新原子序号
                atomSerialNumber++;
            }
        }
    }

    // 关闭文件
    infile.close();
    outfile.close();
}

string NuGraph::indexToBractString(int* index, int indexCount){

	char ss[indexCount+1];
	for(int i=0;i<indexCount;i++){
		ss[i] = '.';
	}

	ss[indexCount] = '\0';// 末尾加上字符串终止符
	for(int i=0;i<indexCount;i++){
		int pi = index[i];//遍历每个碱基的位置 i，从 index 数组中读取其配对的碱基位置 pi。
		if(pi < i) continue;//如果 pi 小于 i，说明当前碱基已经被处理过（因为配对关系是对称的），跳过。
		if(pi == -1) continue;//如果 pi 等于 -1，说明该碱基没有配对，也跳过。

		//确定碱基 i 和 pi 之间哪种括号类型是未使用的（即在此区间内平衡的括号类型）。
		//使用六种不同的括号对 (), [], {}, <>, Aa, Bb, Cc，并通过计数每种类型的开括号和闭括号来判断当前类型是否可用。
		int idA = i;
		int idB = pi;
		int a = 0;
		int b = 0;
		int c = 0;
		int d = 0;
		int e = 0;
		int f = 0;
		int g = 0;
		for(int k=idA+1;k<idB;k++){
			char cc = ss[k];
			if(cc == '(') a++;
			if(cc == ')') a--;
			if(cc == '[') b++;
			if(cc == ']') b--;
			if(cc == '{') c++;
			if(cc == '}') c--;
			if(cc == '<') d++;
			if(cc == '>') d--;
			if(cc == 'A') e++;
			if(cc == 'a') e--;
			if(cc == 'B') f++;
			if(cc == 'b') f--;
			if(cc == 'C') g++;
			if(cc == 'c') g--;
		}
		//根据检查结果，分配未被使用的括号类型给当前的配对碱基 i 和 pi。
		if(a==0) {
			ss[idA] = '(';
			ss[idB] = ')';
		}
		else if(b == 0){
			ss[idA] = '[';
			ss[idB] = ']';
		}
		else if(c == 0){
			ss[idA] = '{';
			ss[idB] = '}';
		}
		else if(d == 0){
			ss[idA] = '<';
			ss[idB] = '>';
		}
		else if(e == 0){
			ss[idA] = 'A';
			ss[idB] = 'a';
		}
		else if(f == 0){
			ss[idA] = 'B';
			ss[idB] = 'b';
		}
		else if(g == 0){
			ss[idA] = 'C';
			ss[idB] = 'c';
		}
	}
	//最后，返回生成的括号字符串表示 RNA 二级结构。
	return string(ss);
}

void NuGraph::clusterContacts(int* array, int array_size, vector<vector<int>>& clusters) {
	vector<bool> visited(array_size, false);  // 记录每个节点是否已经被访问

	// 遍历 array 中的每个节点
	for (int i = 0; i < array_size; i++) {
		if (!visited[i]) {
			vector<int> cluster;  // 当前聚类
			dfsCluster(array, array_size, i, visited, cluster);
			clusters.push_back(cluster);  // 将当前聚类加入聚类列表
		}
	}
}

void NuGraph::dfsCluster(int* array, int array_size, int currentIdx, vector<bool>& visited, vector<int>& cluster) {
	visited[currentIdx] = true;  // 标记当前节点为已访问
	cluster.push_back(array[currentIdx]);  // 将当前节点加入到聚类中

	int currentNode = array[currentIdx];  // 获取当前节点在 seqLen 中的索引

	// 向下查找连接的节点
    for (int i = 0; i < array_size; i++) {
        if (!visited[i]) {
            int neighborNode = array[i];

			double minSDistance = minSquareDistance(currentNode, neighborNode);  // 获取当前节点与邻居节点之间的最小距离
            // 检查两个节点是否在8Å范围内连接
            if (minSDistance < 64.0) {
                dfsCluster(array, array_size, i, visited, cluster);
            }
        }
    }
}

void NuGraph::mergeClusters(vector<vector<int>>& clusters, double threshold) {
	for (size_t i = 0; i < clusters.size(); i++) {
		for (size_t j = i + 1; j < clusters.size(); j++) {
			// 检查两个聚类之间的最近节点距离
			double minDistance = findMinDistanceBetweenClusters(clusters[i], clusters[j]);
			if (minDistance < threshold) {
				// 如果距离小于阈值，将这两个聚类合并
				clusters[i].insert(clusters[i].end(), clusters[j].begin(), clusters[j].end());
				clusters.erase(clusters.begin() + j);  // 删除已合并的聚类
				j--;  // 调整 j 的索引，因为集合已经改变
			}
		}
	}
}

double NuGraph::findMinDistanceBetweenClusters(const vector<int>& clusterA, const vector<int>& clusterB) {
	double minDistance = numeric_limits<double>::max();  // 初始化为最大值

	for (int nodeA : clusterA) {
		for (int nodeB : clusterB) {
			double distance = squareDistance(allNodes[nodeA]->baseConf->coords[0], allNodes[nodeB]->baseConf->coords[0]);
			if (distance < minDistance) {
				minDistance = distance;
			}
		}
	}

	return sqrt(minDistance);  // 返回最小距离
}

double NuGraph::minSquareDistance(int Node1,int Node2) {
	vector<XYZ> Node1Coords, Node2Coords;

	// 收集 neighbor 的所有原子坐标
	for(int i=0;i<allNodes[Node1]->baseConf->rot->atomNum;i++){
		Node1Coords.push_back(allNodes[Node1]->baseConf->coords[i]);
	}
	for(int i=0;i<allNodes[Node2]->baseConf->rot->atomNum;i++){
		Node2Coords.push_back(allNodes[Node2]->baseConf->coords[i]);
	}
	for(int i=0;i<8;i++){
		Node1Coords.push_back(allNodes[Node1]->riboseConf->coords[i]);
		Node2Coords.push_back(allNodes[Node2]->riboseConf->coords[i]);
	}
	for(int i=0;i<4;i++){
		Node1Coords.push_back(allNodes[Node1]->phoConf->coords[i]);
		Node2Coords.push_back(allNodes[Node2]->phoConf->coords[i]);
	}

	// 找到最小的平方距离
	double minSquareDistance = 1e9; // 初始为一个很大的值
	for (const auto& coord1 : Node1Coords) {
		for (const auto& coord2 : Node2Coords) {
			double dist = squareDistance(coord1, coord2);
			minSquareDistance = min(minSquareDistance, dist);
		}
	}
	Node1Coords.clear();
	Node2Coords.clear();
	return minSquareDistance;
}

void NuGraph::refineSubGraph(const string& inputFile, RotamerLib* rotLib, AtomLib* atLib, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib, EdgeInformationLib* eiLib, RnaEnergyTable* et){

	double clashRescale = 1.0;
	double connectRescale = 1.0;

	int* subGraphPosList = new int[this->seqLen]();
	int* fixedPositions = new int[this->seqLen]();
	vector<int> count(seqLen, 0);  // 初始化每个元素的计数
	vector<int> outsubGraphPosList;
	vector<vector<int>> Bclusters;
	vector<vector<int>> Cclusters;
	// int corePos = 100;

	while (true) {
        int corePos = rand() % (this->seqLen);  // 随机选择一个位置作为 corePos
		corePos = 321;
		NuGraph* subGraph = new NuGraph(inputFile, rotLib, atLib, pairLib, moveLib, eiLib, et);
		this->generateSubGraph(inputFile, corePos, subGraphPosList, fixedPositions, subGraph, outsubGraphPosList, Bclusters, Cclusters);
		
		vector<pair<int, int>> ranges; // 存储连续序列的首尾
		if (outsubGraphPosList.empty()) {
			cout << "outsubGraphPosList is empty" << endl;
			exit(0);
		}
		int start = outsubGraphPosList[0];
		for (size_t i = 1; i < outsubGraphPosList.size(); i++) {
			// 检查当前元素与前一个元素的差值
			if (outsubGraphPosList[i] != outsubGraphPosList[i - 1] + 1) {
				// 找到一个连续序列的结束
				ranges.emplace_back(start, outsubGraphPosList[i - 1]);
				start = outsubGraphPosList[i]; // 更新开始位置
			}
		}
		// 添加最后一个连续序列
		ranges.emplace_back(start, outsubGraphPosList.back());
		// 输出结果
		for (const auto& range : ranges) {
			cout << "Continuous sequence: " << range.first << " to " << range.second << endl;
		}

		//在run MC之前/更新节点坐标前记录cm关系
		vector<vector<CsMove>> csMoveLists;
		for(size_t i = 0; i < Bclusters.size(); i++){
			vector<CsMove> csMoveList;
			int Bpos = Bclusters[i][0];
			for(size_t j = 0; j < Cclusters[i].size(); j++) {
				int Cpos = Cclusters[i][j];
				//if(Bpos < Cpos) {
					csMoveList.push_back(this->allEdges[Bpos*this->seqLen + Cpos]->cm);
				//}
			}
			csMoveLists.push_back(csMoveList);
		}

		subGraph->initRandWeight();
		NuTree* subtree = new NuTree(subGraph);
		subGraph->MST_kruskal(subtree);
		subtree->updateNodeInfo(1.0, 1.0);
		subtree->updateEdgeInfo(1.0, 1.0);
		subtree->updateSamplingInfo();

		cout << "run MC" << endl;
		graphInfo* subgi = subtree->runAtomicMC();
		const string& outputFile1 = inputFile + ".subgraph.pdb";
		subgi->printPDB(outputFile1);

		vector<int> subconnectionBreakPoints;
		int indexcBP = 0;
		for (const auto& range : ranges) {
			subconnectionBreakPoints.push_back(range.first);
			indexcBP++;
			cout << "connectionBreakPoints[" << indexcBP << "] = " <<  range.first << ' ';
		}
		cout << endl;

		for (size_t i = 0; i < outsubGraphPosList.size(); i++) {
			// 检查 outsubGraphPosList[i] 是否不在 subconnectionBreakPoints 中
			if (find(subconnectionBreakPoints.begin(), subconnectionBreakPoints.end(), outsubGraphPosList[i]) == subconnectionBreakPoints.end()) {
				this->allNodes[outsubGraphPosList[i]] = subGraph->allNodes[i];
				cout << "node " << outsubGraphPosList[i] << " updated" << endl;
			}
		}

		for(size_t i = 0; i < Bclusters.size(); i++){
			int Bpos = Bclusters[i][0];
			cout << "Bpos:" << Bpos << endl;
			for(size_t j = 0; j < Cclusters[i].size(); j++) {
				int Cpos = Cclusters[i][j];
				cout << "Cpos:" << Cpos << " ";
				LocalFrame cs1 = this->allEdges[Bpos*this->seqLen + Cpos]->nodeA->baseConfTmp->cs1 + csMoveLists[i][j];
				this->allEdges[Bpos*this->seqLen + Cpos]->nodeB->updateCoordinate(cs1);
				this->allEdges[Bpos*this->seqLen + Cpos]->nodeB->acceptCoordMove();
			}
			cout << endl;
		}


		for(size_t i = 0;i < subconnectionBreakPoints.size();i++){
			if(subconnectionBreakPoints[i] > 0){//第一个位置怎么办？
				int j = subconnectionBreakPoints[i]-1;
				cout << "subconnectionBreakPoints[" << i << "] = " <<  subconnectionBreakPoints[i] << ' ';
				this->et->pb->buildPhosphate(this->allNodes[j]->riboseConfTmp, this->allNodes[j+1]->riboseConfTmp, this->allNodes[j]->phoConfTmp);
				this->allNodes[j+1]->phoConf->copyValueFrom(this->allNodes[j+1]->phoConfTmp);
			}
		}
		cout << "update NuEdge" << endl;

		for(size_t i = 0; i < this->seqLen; i++){
			for(size_t j = 0; j < this->seqLen; j++){
				this->allEdges[i*this->seqLen +j] = new NuEdge(this->allNodes[i], this->allNodes[j], this);
			}
		}

		cout << "clean ranges" << endl;
		ranges.clear();
		csMoveLists.clear();
		subconnectionBreakPoints.clear();
		Bclusters.clear();
		Cclusters.clear();
		
		cout << "clean subtree" << endl;
		delete subgi;//delete gi;
		delete subtree;//delete tree;

		delete [] subGraph->allNodes;
		subGraph->allNodes = nullptr;
		delete [] subGraph->allEdges;
		subGraph->allEdges = nullptr;
		delete subGraph;// 由于子图的allNodes和allEdges和主图中subGraphPosList中的元素对应的allNodes和allEdges指向相同的内存，
		// 所以这一步先不释放内存，直接将对应的指针置空，等待主图释放内存。
		cout << "clean subGraph" << endl;
		break;

		

        // 更新 count 数组，记录 subGraphPosList 中的元素出现次数
        for (int pos : outsubGraphPosList) {
            count[pos]++;
        }
        // 检查是否所有元素都出现至少三次
        bool allElementsMet = true;
        for (int i = 0; i < this->seqLen; ++i) {
            if (count[i] < 3) {
                allElementsMet = false;
                break;
            }
        }
        // 如果所有元素都至少出现三次，结束循环
        if (allElementsMet) {
            break;
        }
		outsubGraphPosList.clear();

    }
	delete [] subGraphPosList;
    delete [] fixedPositions;
	count.clear();
	cout << "clean count" << endl;
}

graphInfo* NuGraph::getGraphInfo(){
	graphInfo* gi = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, getTotalEnergyForModelSelection(), atLib, 0);
	return gi;
}

graphInfo* NuGraph::getGraphInfoCG(){
	graphInfo* gi = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, totalEnergyCG(1.0, 1.0), atLib, 1);
	return gi;
}

graphInfo* NuGraph::getGraphInfo(double ene){
	graphInfo* gi = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, ene, atLib, 0);
	return gi;
}

graphInfo* NuGraph::getGraphInfoCG(double ene){
	graphInfo* gi = new graphInfo(seqLen, seq, connectToDownstream, fixed, allNodes, ene, atLib, 1);
	return gi;
}

} /* namespace NSPpredNA */
