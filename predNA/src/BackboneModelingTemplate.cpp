/*
 * BackboneModelingTemplate.cpp
 *
 */

#include "predNA/BackboneModelingTemplate.h"

namespace NSPpredna {

BackboneModelingTemplate::BackboneModelingTemplate(const string& inputPDB, ForceFieldPara* para) {

	cout << "init backbone modeling template: " << endl;
	RNAPDB pdb(inputPDB, "xxxx");
	vector<RNABase*> baseList = pdb.getBaseList();
	

	for(int i=0;i<baseList.size();i++){
		baseList[i]->baseSeqID = i;
	}

	for(int i=0;i<baseList.size();i++){
		mutNodeList.push_back(i);
	}

	this->seqLen = baseList.size();
	this->nodes = new BRNode*[seqLen];
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];



	this->para = para;

	this->et = new RnaEnergyTableSimple(para);
	this->rotLib = new RotamerLib(para);
	this->pb = new PO3Builder(para);

	for(int i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		RiboseRotamer* rot;
		if(!baseList[i]->backboneComplete())
			rot = rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
		else
			rot = new RiboseRotamer(baseList[i]);

		this->nodes[i] = new BRNode(baseList[i], rotLib);
		nodes[i]->riboseConf->updateRotamer(rot);
		nodes[i]->riboseConfTmp->updateRotamer(rot);
	}




	this->sepTable = new int[seqLen*seqLen];
	this->allBaseRiboseE = new double[seqLen*seqLen];
	this->tmpBaseRiboseE = new double[seqLen*seqLen];
	this->allBasePhoE = new double[seqLen*seqLen];
	this->tmpBasePhoE = new double[seqLen*seqLen];
	this->allRiboseRiboseE = new double[seqLen*seqLen];
	this->tmpRiboseRiboseE = new double[seqLen*seqLen];
	this->allRibosePhoE = new double[seqLen*seqLen];
	this->tmpRibosePhoE = new double[seqLen*seqLen];
	this->allPhoPhoE = new double[seqLen*seqLen];
	this->tmpPhoPhoE = new double[seqLen*seqLen];
	this->allRotE = new double[seqLen];
	this->tmpRotE = new double[seqLen];
	this->allRcE = new double[seqLen];
	this->tmpRcE = new double[seqLen];


	for(int i=0;i<seqLen-1;i++){
		if(baseList[i]->connectToNeighbor(baseList[i+1])){
			this->connectToDownstream[i] = true;
			nodes[i]->connectToNeighbor = true;
		}
		else {
			this->connectToDownstream[i] = false;
			nodes[i]->connectToNeighbor = false;
		}
	}

	connectToDownstream[seqLen-1] = false;
	nodes[seqLen-1]->connectToNeighbor = false;


	vector<RiboseRotamer*> rotList;
	for(int i=0;i<seqLen;i++){
		RiboseRotamer* rot = nodes[i]->riboseConf->rot;
		rotList.push_back(rot);
		LocalFrame cs = nodes[i]->baseConf->cs1;


		if(i>0 && i < seqLen-1 && connectToDownstream[i-1] && connectToDownstream[i]) {
			for(int j=0;j<rot->atomNum;j++){
				this->initBackboneAtomList.push_back(nodes[i]->riboseConf->coords[j]);
			}
			this->initBackboneAtomList.push_back(baseList[i+1]->getAtom("P")->coord);
			this->initBackboneAtomList.push_back(baseList[i+1]->getAtom("O5'")->coord);
		}

	}

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			int ij = i*seqLen+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else if(j == i+2 && connectToDownstream[i] && connectToDownstream[i+1]) sepTable[ij] = 2;
			else if(j == i-2 && connectToDownstream[j] && connectToDownstream[j+1]) sepTable[ij] = -2;
			else sepTable[ij] = 3;
		}
	}

	for(int i=0;i<seqLen;i++){
		nodes[i]->updateChildInfo(nodes, seqLen);

	}

	for(int i=0;i<seqLen;i++){

		buildPho(i, false);
		acceptTmpRotamer(nodes[i], false);
	}

	updateEnergy();


}

BackboneModelingTemplate::BackboneModelingTemplate(RNAPDB* pdb, const string& cnt, const string& csn, ForceFieldPara* para){
	vector<RNABase*> baseList = pdb->getBaseList();
	
	for(int i=0;i<baseList.size();i++){
		baseList[i]->baseSeqID = i;
	}

	if(baseList.size() != cnt.length()){
		cout << "base list: " << baseList.size() << endl;
		cout << "cnt length: " << cnt.length() << endl;
		exit(0);
	}

	this->seqLen = baseList.size();
	this->nodes = new BRNode*[seqLen];
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];




	this->para = para;

	this->et = new RnaEnergyTableSimple(para);
	this->rotLib = new RotamerLib(para);
	this->pb = new PO3Builder(para);

	for(int i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		RiboseRotamer* rot;
		if(!baseList[i]->backboneComplete())
			rot = rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
		else
			rot = new RiboseRotamer(baseList[i]);

		this->nodes[i] = new BRNode(baseList[i], rotLib);
		nodes[i]->riboseConf->updateRotamer(rot);
		nodes[i]->riboseConfTmp->updateRotamer(rot);
	}




	this->sepTable = new int[seqLen*seqLen];
	this->allBaseRiboseE = new double[seqLen*seqLen];
	this->tmpBaseRiboseE = new double[seqLen*seqLen];
	this->allBasePhoE = new double[seqLen*seqLen];
	this->tmpBasePhoE = new double[seqLen*seqLen];
	this->allRiboseRiboseE = new double[seqLen*seqLen];
	this->tmpRiboseRiboseE = new double[seqLen*seqLen];
	this->allRibosePhoE = new double[seqLen*seqLen];
	this->tmpRibosePhoE = new double[seqLen*seqLen];
	this->allPhoPhoE = new double[seqLen*seqLen];
	this->tmpPhoPhoE = new double[seqLen*seqLen];
	this->allRotE = new double[seqLen];
	this->tmpRotE = new double[seqLen];
	this->allRcE = new double[seqLen];
	this->tmpRcE = new double[seqLen];

	for(int i=0;i<baseList.size();i++){
		if(cnt[i] == '-') {
			connectToDownstream[i] = true;
			nodes[i]->connectToNeighbor = true;
		}
		else {
			connectToDownstream[i] = false;
			nodes[i]->connectToNeighbor = false;
		}
	}

	bool* selectNodeList = new bool[baseList.size()];
	for(int i=0;i<baseList.size();i++){
		selectNodeList[i] = false;
	}
	for(int i=0;i<baseList.size()-1;i++){
		if(csn[i] != 'F'){
			for(int j=i-2;j<i+3;j++){
				if(j >=0 && j < seqLen){
					selectNodeList[i] = true;
				}
			}
		}
	}
	for(int i=0;i<seqLen;i++){
		if(selectNodeList[i])
			mutNodeList.push_back(i);
	}

	cout << "backbone modeling node num: " << mutNodeList.size() << endl;

	vector<RiboseRotamer*> rotList;
	for(int i=0;i<seqLen;i++){
		RiboseRotamer* rot = nodes[i]->riboseConf->rot;
		rotList.push_back(rot);
		LocalFrame cs = nodes[i]->baseConf->cs1;


		if(i>0 && i < seqLen-1 && connectToDownstream[i-1] && connectToDownstream[i]) {
			for(int j=0;j<rot->atomNum;j++){
				this->initBackboneAtomList.push_back(nodes[i]->riboseConf->coords[j]);
			}
			this->initBackboneAtomList.push_back(baseList[i+1]->getAtom("P")->coord);
			this->initBackboneAtomList.push_back(baseList[i+1]->getAtom("O5'")->coord);
		}

	}

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			int ij = i*seqLen+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else if(j == i+2 && connectToDownstream[i] && connectToDownstream[i+1]) sepTable[ij] = 2;
			else if(j == i-2 && connectToDownstream[j] && connectToDownstream[j+1]) sepTable[ij] = -2;
			else sepTable[ij] = 3;
		}
	}

	for(int i=0;i<seqLen;i++){
		nodes[i]->updateChildInfo(nodes, seqLen);

	}

	for(int i=0;i<seqLen;i++){

		buildPho(i, false);
		acceptTmpRotamer(nodes[i], false);
	}

	updateEnergy();
}

double BackboneModelingTemplate::buildPho(int seqID, bool verbose) {
	if(seqID < 0) return 0;
	if(seqID >= seqLen-1) return 0;
	if(!connectToDownstream[seqID]) return 0;
	if(verbose) {
		cout << "update pho: " << seqID << endl;
	}
	BRNode* nodeA = nodes[seqID];
	BRNode* nodeB = nodes[seqID+1];
	this->pb->buildPhosphate(nodeA->riboseConfTmp, nodeB->riboseConfTmp, nodeA->phoConfTmp);

	return nodeA->phoConfTmp->ene;
}


void BackboneModelingTemplate::updateTmpRotamer(BRNode* node, RiboseRotamer* rot, bool verbose){

	node->riboseConfTmp->updateRotamer(rot);

	int chainBreakNum = node->phoGroupC.size();
	int indexA, indexB;

	BRNode* nodeA;
	BRNode* nodeB;
	for(int i=0;i<chainBreakNum;i++){
		if(verbose){
			cout << "connection: " << indexA << " " << indexB << endl;
		}
		indexA = node->phoGroupC[i];
		indexB = indexA+1;
		nodeA = nodes[indexA];
		nodeB = nodes[indexB];
		pb->buildPhosphate(nodeA->riboseConfTmp, nodeB->riboseConfTmp, nodeA->phoConfTmp);
	}
}

void BackboneModelingTemplate::clearTmpRotamer(BRNode* nodeSelect, bool verbose){
	int i,j;

	nodeSelect->riboseConfTmp->copyValueFrom(nodeSelect->riboseConf);

	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoConfTmp->copyValueFrom(node->phoConf);
	}
}

void BackboneModelingTemplate::acceptTmpRotamer(BRNode* nodeSelect, bool verbose){
	int i,j;

	nodeSelect->riboseConf->copyValueFrom(nodeSelect->riboseConfTmp);
	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoConf->copyValueFrom(node->phoConfTmp);
	}
}

double BackboneModelingTemplate::mutEnergy(BRNode* node, bool verbose){

	double tot = 0;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	double mutE = 0;

	//int baseGroupANum = node->baseGroupA.size();


	int riboseGroupANum = node->riboseGroupA.size();
	int riboseGroupCNum = node->riboseGroupC.size();

	int phoGroupANum = node->phoGroupA.size();
	int phoGroupCNum = node->phoGroupC.size();

	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		mutE += (nodes[indexX]->riboseConfTmp->rot->energy - nodes[indexX]->riboseConf->rot->energy);

	}

	for(i=0;i<phoGroupCNum;i++){
		indexX = node->phoGroupC[i];
		mutE += (nodes[indexX]->phoConfTmp->ene - nodes[indexX]->phoConf->ene);
	}

	/*
	 * base-ribose energy
	 */

	//A-C
	for(i=0;i<seqLen;i++){
		indexX = i;
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
				mutE += (getBaseRiboseEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getBaseRiboseEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	/*
	 * base-pho energy
	 */

	//A-C
	for(i=0;i<seqLen;i++){
		indexX = i;
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
				mutE += (getBasePhoEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getBasePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	/*
	 * ribose-ribose energy
	 */

	//A-C
	for(i=0;i<riboseGroupANum;i++){
		indexX = node->riboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
			//	printf("mut: %-2d %-2d %2d %7.3f %7.3f %7.3f\n", indexX, indexY, sepTable[pi], getRiboseRiboseEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose), getRiboseRiboseEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose), getRiboseRiboseEnergyBM(nodeY, nodeX, -sepTable[pi], et, verbose));

				mutE += (getRiboseRiboseEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getRiboseRiboseEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}


	/*
	 * ribose-pho energy
	 */

	//A-C
	for(i=0;i<riboseGroupANum;i++){
		indexX = node->riboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
				mutE += (getRibosePhoEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getRibosePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	//C-A
	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = node->phoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
				mutE += (getRibosePhoEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getRibosePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	//C-C
	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
				mutE += (getRibosePhoEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getRibosePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	/*
	 * pho-pho energy
	 */

	//A-C
	for(i=0;i<phoGroupANum;i++){
		indexX = node->phoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
				mutE += (getPhoPhoEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getPhoPhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	//C-C
	for(i=0;i<phoGroupCNum;i++){
		indexX = node->phoGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			if(squareDistance(nodeX->baseConf->cs1.origin_, nodeY->baseConf->cs1.origin_) < 144){
				mutE += (getPhoPhoEnergyTmpBM(nodeX, nodeY, sepTable[pi], et, verbose) - getPhoPhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	return mutE;
}

BRTreeInfo* BackboneModelingTemplate::toTreeInfo(){
	vector<BRNode*> flexNodes;
	for(int i=0;i<seqLen;i++){
		flexNodes.push_back(nodes[i]);
	}
	return new BRTreeInfo(seqLen, seq, connectToDownstream, nodes, flexNodes, 0.0, rotLib);
}

double BackboneModelingTemplate::totEnergy(bool verbose){
	double tot = 0;
	int i,j;

	//cout << "total energy:" << endl;
	//ribose-ribose, pho-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];

			//printf("%-2d %-2d %8.3f\n", i, j, getRiboseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose));
			tot += getRiboseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			tot += getPhoPhoEnergyBM(nodeA, nodeB, sep, et, verbose);
		}
	}


	//base-ribose, base-pho, ribose-pho

	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			tot += getBaseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			tot += getBasePhoEnergyBM(nodeA, nodeB, sep,  et, verbose);
			tot += getRibosePhoEnergyBM(nodeA, nodeB, sep, et, verbose);
		}
	}



	for(int i=0;i<seqLen;i++){
		if(verbose){
			cout << "rotEnergy: " << i << " " << nodes[i]->riboseConf->rot->energy << endl;
		}
		tot += nodes[i]->riboseConf->rot->energy;
	}


	for(int i=0;i<seqLen;i++){
		if(verbose){
			cout << "phoEnergy: " << i << " " <<  nodes[i]->phoConf->ene << endl;
		}
		tot += nodes[i]->phoConf->ene;
	}

	return tot;
}

void BackboneModelingTemplate::printEnergyDetail(){
	double e = 0;

	for(int i=0;i<seqLen;i++){
		e += nodes[i]->riboseConf->rot->energy;
		e += nodes[i]->phoConf->ene;
//		printf("Node: %2d RotamerE: %9.3f PhoE: %9.3f\n", i, nodes[i]->rot->energy, nodes[i]->pho.e);
	}
	bool verbose = false;

	for(int i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(int j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int pi = i*seqLen+j;
			int pj = j*seqLen+i;
			int sepi = sepTable[pi];
			int sepj = sepTable[pj];
			{
				double tot = 0;
				double bb = 0;
				double br = getBaseRiboseEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double rb = getBaseRiboseEnergyBM(nodeB, nodeA, sepj, et, verbose);
				double bp = getBasePhoEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double pb = getBasePhoEnergyBM(nodeB, nodeA, sepj, et, verbose);
				double rr = getRiboseRiboseEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double rp = getRibosePhoEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double pr = getRibosePhoEnergyBM(nodeB, nodeA, sepj, et, verbose);
				double pp = getPhoPhoEnergyBM(nodeA, nodeB, sepi, et, verbose);
				tot = bb + br + rb + bp + pb + rr + rp + pr + pp;
				if(abs(tot) > 0.001)
					printf("Pair: %2d-%-2d bb %7.3f br %7.3f rb %7.3f bp %7.3f pb %7.3f rr %7.3f rp %7.3f pr %7.3f pp %7.3f\n",i,j, bb, br, rb, bp, pb, rr, rp, pr, pp);
				e += tot;
			}

			{
				double tot = 0;
				double bb = 0;
				double br = getBaseRiboseEnergyTmpBM(nodeA, nodeB, sepi, et, verbose);
				double rb = getBaseRiboseEnergyTmpBM(nodeB, nodeA, sepj, et, verbose);
				double bp = getBasePhoEnergyTmpBM(nodeA, nodeB, sepi, et, verbose);
				double pb = getBasePhoEnergyTmpBM(nodeB, nodeA, sepj, et, verbose);
				double rr = getRiboseRiboseEnergyTmpBM(nodeA, nodeB, sepi, et, verbose);
				double rp = getRibosePhoEnergyTmpBM(nodeA, nodeB, sepi, et, verbose);
				double pr = getRibosePhoEnergyTmpBM(nodeB, nodeA, sepj, et, verbose);
				double pp = getPhoPhoEnergyTmpBM(nodeA, nodeB, sepi, et, verbose);
				tot = bb + br + rb + bp + pb + rr + rp + pr + pp;
				if(abs(tot) > 0.001)
					printf("TmPr: %2d-%-2d bb %7.3f br %7.3f rb %7.3f bp %7.3f pb %7.3f rr %7.3f rp %7.3f pr %7.3f pp %7.3f\n",i,j, bb, br, rb, bp, pb, rr, rp, pr, pp);
			}

		}
	}

	printf("total energy: %9.3f\n",  e);
	cout << endl;
}


double BackboneModelingTemplate::fragMC(){

	int stepNum = 2000;
	double curEne = totEnergy(false);
	double lastEne = curEne;
	double mutE;
	srand(time(0));
	int randIndex, randType;
	BRNode* node;
	RiboseRotamer* rotMut;

	double minEne = 9999999999.9;
	double rmsd = 999.9;
	int fragLen = 7;
	for(int startIndex=0;startIndex<seqLen-fragLen+4;startIndex+=3){
		cout << startIndex << endl;
		for(double T=20.0;T>0.05;T=T*0.8){
			for(int k=0;k<stepNum;k++){
				randIndex = rand()%fragLen + startIndex;
				if(randIndex >= seqLen) continue;
				node = nodes[randIndex];
				randType = rand()%2;
				rotMut = rotLib->riboseRotLib->getRandomRotamer(node->baseType);
				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}
	}

	stepNum = 3000;
	for(int startIndex=0;startIndex<seqLen-fragLen+4;startIndex+=3){
		double T0 = 10.0;
		double T1 = 0.01;
		cout << startIndex << endl;
		for(double T=T0;T>T1;T=T*0.9){
			for(int k=0;k<stepNum;k++){

				randIndex = rand()%fragLen + startIndex;
				if(randIndex >= seqLen) continue;
				node = nodes[randIndex];
				randType = rand()%2;
				rotMut = rotLib->riboseRotLib->getRandomRotamer(node->baseType);

				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;

						printf("R1: ene: %8.3f rmsd: %6.4f\n", curEne,rmsd);
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}
	}

	for(int startIndex=0;startIndex<seqLen-fragLen+4;startIndex+=3){
		double T0 = 10.0;
		double T1 = 0.01;
		cout << startIndex << endl;
		for(double T=T0;T>T1;T=T*0.9){
			for(int k=0;k<stepNum;k++){

				randIndex = rand()%fragLen + startIndex;
				if(randIndex >= seqLen) continue;
				node = nodes[randIndex];
				randType = rand()%2;
				rotMut = rotLib->riboseRotLib->getRandomRotamer(node->baseType);

				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
						printf("R2: ene: %8.3f rmsd: %6.4f\n", curEne,rmsd);
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}
	}

	for(int i=0;i<seqLen;i++){
		node = nodes[i];
		for(int j=0;j<200;j++){
			int type1 = node->riboseConf->rot->rotTypeLv1;
			rotMut = rotLib->riboseRotLib->rotLib[node->baseType][type1*200+j];
			updateTmpRotamer(node, rotMut, false);
			mutE = mutEnergy(node, false);
			if(mutE < 0) {
				curEne += mutE;
				acceptTmpRotamer(node, false);
				if(curEne < minEne){
					minEne = curEne;
				}
			}
			else{
				clearTmpRotamer(node, false);
			}
		}
	}
	return rmsd;
}

void BackboneModelingTemplate::debug(){
	int stepNum = 100;
	double curEne = totEnergy(false);
	double lastEne = curEne;
	double mutE;
	srand(time(0));
	int randIndex, randType;
	BRNode* node;
	RiboseRotamer* rotMut;
	double T = 2.0;
	double minEne = 99999.9;

	double e = totEnergy(true);
	cout << "curE: " << curEne << " totE: " << e << endl;

	for(int k=0;k<1000;k++){
		cout << "step: " << k << endl;
		randIndex = rand()%seqLen;
		node = nodes[randIndex];
		randType = rand()%2;
		rotMut = rotLib->riboseRotLib->getRandomRotamer(node->baseType);

		printEnergyDetail();

		updateTmpRotamer(node, rotMut,false);
		mutE = mutEnergy(node,false);


		if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
			cout << "accept: " << mutE << endl;
			curEne += mutE;
			acceptTmpRotamer(node, false);
			if(curEne < minEne){
				minEne = curEne;
			}
		}
		else{
			cout <<"reject" << mutE <<  endl;
			clearTmpRotamer(node, false);
		}

		printEnergyDetail();

		/*
		checkTotalEnergy();
		checkTotalEnergyTmp();
		*/
		double totE = totEnergy(false);
		if(abs(curEne - totE) > 0.001){
			cout << "error" << endl;
		}
		printf(" T=%8.5f  addEne: %9.3f totEne: %9.3f rms: %5.3f\n", T, curEne, totE, rms());
	}
}

double BackboneModelingTemplate::runMC(){

	int stepNum = 1000*seqLen;
	double curEne = totEnergy(false);
	double mutE;
	srand(time(0));
	int randIndex, randType;
	BRNode* node;
	RiboseRotamer* rotMut;

	int mutNodeNum = this->mutNodeList.size();

	double minEne = 9999999999.9;
	double minERMS = 99.9;

	for(int n=0;n<1;n++){
		for(double T=10.0;T>0.01;T=T*0.9){
			int acNum = 0;
			for(int k=0;k<stepNum;k++){
				randIndex = rand()%mutNodeNum;
				node = nodes[mutNodeList[randIndex]];
				randType = rand()%2;
				rotMut = rotLib->riboseRotLib->getRandomRotamer(node->baseType);
				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					acNum ++;
					if(curEne < minEne){
						minEne = curEne;
						minERMS = rms();
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
			double totE = totEnergy(false);
			printf("R: %d T=%8.5f accept=%5d addEne: %9.3f totEne: %9.3f minEne: %9.3f rms: %5.3f\n", n, T, acNum, curEne, totE,minEne, rms());
		}

		for(int i=0;i<mutNodeNum;i++){
			node = nodes[mutNodeList[i]];
			for(int j=0;j<1500;j++){
			//	int type1 = node->riboseConf->rot->rotTypeLv1;
				rotMut = rotLib->riboseRotLib->rotLib[node->baseType][j];
				updateTmpRotamer(node, rotMut,false);
				mutE = mutEnergy(node,false);
				if(mutE < 0) {
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
						minERMS = rms();
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}

		for(int i=0;i<mutNodeNum;i++){
			node = nodes[mutNodeList[i]];
			for(int j=0;j<1500;j++){
			//	int type1 = node->riboseConf->rot->rotTypeLv1;
				rotMut = rotLib->riboseRotLib->rotLib[node->baseType][j];
				updateTmpRotamer(node, rotMut,false);
				mutE = mutEnergy(node,false);
				if(mutE < 0) {
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
						minERMS = rms();
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}
	}

	return minERMS;
}

void BackboneModelingTemplate::checkTotalEnergy(){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e, e1, e2;



	//base-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBaseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allBaseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allBaseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//base-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBasePhoEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allBasePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allBasePhoE[i*seqLen+j], e);
			}
		}
	}

	//ribose-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRiboseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allRiboseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allRiboseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//ribose-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRibosePhoEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allRibosePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allRibosePhoE[i*seqLen+j], e);
			}
		}
	}

	//pho-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getPhoPhoEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allPhoPhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: pho: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allPhoPhoE[i*seqLen+j], e);
			}
		}
	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->riboseConf->rot->energy;
		if(abs(allRotE[i] - e) > 0.001)
			cout << "DIFF rot: " << i << " " << allRotE[i] << " " << e << endl;

	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->phoConf->ene;
		if(abs(allRcE[i] -e) > 0.001)
			cout << "DIFF connect: " << i << " " << allRcE[i] << " " << e << endl;
	}
}

void BackboneModelingTemplate::checkTotalEnergyTmp(){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e, e1, e2;

	//base-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBaseRiboseEnergyTmpBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpBaseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpBaseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//base-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBasePhoEnergyTmpBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpBasePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpBasePhoE[i*seqLen+j], e);
			}
		}
	}

	//ribose-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRiboseRiboseEnergyTmpBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpRiboseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpRiboseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//ribose-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRibosePhoEnergyTmpBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpRibosePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpRibosePhoE[i*seqLen+j], e);
			}
		}
	}

	//pho-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getPhoPhoEnergyTmpBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpPhoPhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: pho: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpPhoPhoE[i*seqLen+j], e);
			}
		}
	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->riboseConfTmp->rot->energy;
		if(abs(tmpRotE[i] - e) > 0.001)
			cout << "DIFF rot: " << i << " " << tmpRotE[i] << " " << e << endl;
	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->phoConfTmp->ene;
		if(abs(tmpRcE[i] -e) > 0.001)
			cout << "DIFF connect: " << i << " " << tmpRcE[i] << " " << e << endl;
	}
}

BackboneModelingTemplate::~BackboneModelingTemplate() {
	for(int i=0;i<seqLen;i++){
		delete nodes[i];
	}
	delete [] nodes;
	delete [] seq;
	delete [] connectToDownstream;
	delete pb;
	delete et;
	delete rotLib;

	delete [] allBaseRiboseE;
	delete [] tmpBaseRiboseE;
	delete [] allBasePhoE;
	delete [] tmpBasePhoE;
	delete [] allRiboseRiboseE;
	delete [] tmpRiboseRiboseE;
	delete [] allRibosePhoE;
	delete [] tmpRibosePhoE;
	delete [] allPhoPhoE;
	delete [] tmpPhoPhoE;
	delete [] allRotE;
	delete [] tmpRotE;
	delete [] allRcE;
	delete [] tmpRcE;
}

} /* namespace NSPforcefield */
