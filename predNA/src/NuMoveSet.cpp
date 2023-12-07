/*
 * NuMoveSet.cpp
 *
 *  Created on: 2023Äê11ÔÂ23ÈÕ
 *      Author: nuc
 */

#include <predNA/NuMoveSet.h>

namespace NSPpredNA {

IndividualNuPairMoveSet::IndividualNuPairMoveSet(int sep, int pairType, int clusterID, OrientationIndex* oi) {
	// TODO Auto-generated constructor stub
	this->sep = sep;
	this->pairType = pairType;
	this->clusterID = clusterID;

	string path = NSPdataio::datapath()+"pairMove2/";
	string augc = "AUGC";
	int typeA = pairType/4;
	int typeB = pairType%4;
	string fileName;

	char xx[20];
	sprintf(xx, "%c%c%d", augc[typeA], augc[typeB], clusterID);

	if(sep == 1 || sep == -1)
		fileName = path+"nb/"+string(xx)+".move";
	else
		fileName = path+"nnb/"+string(xx)+".move";

	ifstream file;
	file.open(fileName.c_str(), ios::in);
	int moveID, binID;
	double p;
	while(file >> moveID >> p >> binID){
		if(sep > 0)
			moveIndexList[binID].push_back(moveID);
		else if(sep == -1){
			CsMove move = oi->index1000ToCsMove(moveID);
			CsMove revMove = move.reverse();
			moveIndexList[binID].push_back(oi->moveToIndex1000(revMove));
		}
	}
	file.close();
}

CsMove IndividualNuPairMoveSet::getRandomMove(OrientationIndex* oi){
	int selectBin = rand()% 20;
	int selectID = rand()% moveIndexList[selectBin].size();
	CsMove cm =  oi->index1000ToCsMove(moveIndexList[selectBin][selectID]);
	return cm;
}

IndividualNuPairMoveSet::~IndividualNuPairMoveSet() {
	// TODO Auto-generated destructor stub
}


NuPairMoveSetLibrary::NuPairMoveSetLibrary(){

	BasePairLib* bpLib = new BasePairLib();
	this->oi = new OrientationIndex();


	for(int i=0;i<16;i++){
		int clusterNum = bpLib->nbBasePairNum[i];
		for(int j=0;j<clusterNum;j++){
			nbMoveList[i].push_back(new IndividualNuPairMoveSet(1, i, j, oi));
			revNbMoveList[i].push_back(new IndividualNuPairMoveSet(-1, i, j, oi));
		}

		clusterNum = bpLib->nnbBasePairNum[i];
		for(int j=0;j<clusterNum;j++){
			nnbMoveList[i].push_back(new IndividualNuPairMoveSet(2, i, j, oi));
		}
	}


	delete bpLib;
}

NuPairMoveSetLibrary::~NuPairMoveSetLibrary(){
	for(int i=0;i<16;i++){
		for(int j=0;j<nbMoveList[i].size();j++){
			delete nbMoveList[i][j];
		}

		for(int j=0;j<revNbMoveList[i].size();j++){
			delete revNbMoveList[i][j];
		}

		for(int j=0;j<nnbMoveList[i].size();j++){
			delete nnbMoveList[i][j];
		}
	}
}

MixedNuPairCluster::MixedNuPairCluster(int sep, int pairType, NuPairMoveSetLibrary* lib){
	this->sep = sep;
	this->pairType = pairType;
	this->moveLib = lib;
	this->fixedNativeCM = false;
}

void MixedNuPairCluster::updateEdgeInformation(EdgeInformation* ei){
	double psum = 0;

	clusterIDList.clear();
	clusterPList.clear();

	for(int i=0;i<ei->totalClusterNum;i++){
		double p = ei->pCluster[i];
		if(p >0){
			clusterIDList.push_back(i);
			clusterPList.push_back(p);
			psum += p;
		}
	}

	if(psum == 0) {
		cout << "edge information error: total probability is zero!" << endl;
		exit(0);
	}

	double pAdd = 0;
	for(int i=0;i<10000;i++){
		randPool[i] = 0;
	}

	for(int i=0;i<clusterIDList.size();i++){
		int start = (int)(pAdd*10000);
		pAdd += clusterPList[i]/psum;
		int end = (int)(pAdd*10000);

		if(end > 10000){
			cout << "pAdd " << pAdd << endl;
			exit(0);
		}

		for(int j=start;j<end;j++){
			randPool[j] = clusterIDList[i];
		}
	}
}

CsMove MixedNuPairCluster::getRandomMove(){
	int cluster = randPool[rand()%10000];
	if(sep == 1)
		return moveLib->nbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);
	else if(sep == -1)
		return moveLib->revNbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);
	else
		return moveLib->nnbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);

}

MixedNuPairCluster::~MixedNuPairCluster(){

}


} /* namespace NSPpredNA */
