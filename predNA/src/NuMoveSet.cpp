/*
 * NuMoveSet.cpp
 *
 *  Created on: 2023��11��23��
 *      Author: nuc
 */

#include <predNA/NuMoveSet.h>
#include "tools/ThreadPool.h"
#include "dataio/dirOperations.h"
#include <string.h>
#ifdef TIMING
#include <time.h>
#endif

namespace NSPpredNA {

using namespace NSPdataio;
using namespace NSPthread;

IndividualNuPairMoveSet::IndividualNuPairMoveSet(int sep, int pairType, int clusterID, OrientationIndex* oi, BinaryBook* bb) {
	// TODO Auto-generated constructor stub
	this->sep = sep;
	this->pairType = pairType;
	this->clusterID = clusterID;

	string augc = "AUGC";
	int typeA = pairType/4;
	int typeB = pairType%4;
	string fileName;

	char xx[20];
	sprintf(xx, "%c%c%d", augc[typeA], augc[typeB], clusterID);

	char rxx[20];
	sprintf(rxx, "%c%c%d", augc[typeB], augc[typeA], clusterID);

	char yy[20];
	sprintf(yy, "%d", clusterID);

	if(!bb){
		string path = NSPdataio::datapath()+"pairMove4/";

		if(pairType == -1) //neighbor non-contact base pair
		{
			fileName = path+"nb/nonContact/nc" + string(yy)+".move";
		}
		else if(sep == 1)
			fileName = path+"nb/contact/"+string(xx)+".move";
		else if(sep == -1)
			fileName = path+"nb/contact/"+string(rxx)+".move";
		else
			fileName = path+"nnb/"+string(xx)+".move";

		ifstream file;
		file.open(fileName.c_str(), ios::in);
		int moveID, binID;
		double p;
		while(file >> moveID >> p >> binID){
			CsMove m1 = oi->index1000ToCsMove(moveID);
			double d = m1.oriMove.length();
			if(d < 0.001){
				printf("BAD MOVE: %s %d\n", fileName, moveID);
			}

			if(sep > 0)
				moveIndexList[binID].push_back(moveID);
			else if(sep == -1){
				CsMove revMove = m1.reverse();
				moveIndexList[binID].push_back(oi->moveToIndex1000(revMove));
			}
		}
		file.close();
	} else {
		string tableName = string(xx) + ".move";
		BinaryTable* btab = bb->tables_map.at(tableName);
		auto& moveIDCol = get<BinaryColumn<int>>(*(btab->cols[0]));
		auto& pCol = get<BinaryColumn<double>>(*(btab->cols[1]));
		auto& binIDCol = get<BinaryColumn<int>>(*(btab->cols[2]));
		for(int i=0;i<btab->nRow;i++) {
			int moveID = moveIDCol[i];
			double p = pCol[i];
			int binID = binIDCol[i];

			CsMove m1 = oi->index1000ToCsMove(moveID);
			double d = m1.oriMove.length();
			if(d < 0.001){
				printf("%s %d\n", tableName, moveID);
			}

			if(sep > 0)
				moveIndexList[binID].emplace_back(moveID);
			else if(sep == -1){
				CsMove revMove = m1.reverse();
				moveIndexList[binID].emplace_back(oi->moveToIndex1000(revMove));
			}
		}
	}
}

int IndividualNuPairMoveSet::dump(ostream& outs) {
	char z0[4] = {'\0'};
	outs.write(reinterpret_cast<char*>(&sep), sizeof(int));
	outs.write(reinterpret_cast<char*>(&pairType), sizeof(int));
	outs.write(reinterpret_cast<char*>(&clusterID), sizeof(int));
	outs.write(reinterpret_cast<char*>(z0), sizeof(int));
	int vecSize[50];
	for(int i=0;i<50;i++) {
		vecSize[i] = moveIndexList[i].size();
	}
	outs.write(reinterpret_cast<char *>(vecSize), 50*sizeof(int));
	for(int i=0;i<50;i++) {
		outs.write(reinterpret_cast<char*>(&moveIndexList[i][0]), vecSize[i]*sizeof(int));
		if(vecSize[i]%2) {
			outs.write(reinterpret_cast<char*>(z0), sizeof(int));  // align Bytes
		}
	}
	return EXIT_SUCCESS;
}

int IndividualNuPairMoveSet::load(istream& ins) {
	int tint[4];
	ins.read(reinterpret_cast<char*>(tint), sizeof(int)*4);
	sep = tint[0];
	pairType = tint[1];
	clusterID = tint[2];
	int vecSize[50];
	ins.read(reinterpret_cast<char*>(vecSize), sizeof(int)*50);
	for(int i=0;i<50;i++) {
		moveIndexList[i].resize(vecSize[i]);
		if(vecSize[i]%2) {
			int memSize = vecSize[i]+1;
			int tint2[memSize];
			ins.read(reinterpret_cast<char*>(tint2), memSize*sizeof(int));
			memcpy(&moveIndexList[i][0], tint2, vecSize[i]*sizeof(int));
		} else {
			ins.read(reinterpret_cast<char*>(&moveIndexList[i][0]), vecSize[i]*sizeof(int));
		}
	}
	return EXIT_SUCCESS;
}

CsMove IndividualNuPairMoveSet::getRandomMove(OrientationIndex* oi){

	int selectBin = rand()% 50;
	int selectID = rand()% moveIndexList[selectBin].size();
	
	CsMove cm =  oi->index1000ToCsMoveWithRandPerturbation(moveIndexList[selectBin][selectID]);
	cm.subClusterID = selectBin;
	return cm;
}

CsMove IndividualNuPairMoveSet::getRandomMoveWithFixedSubCluster(OrientationIndex* oi, int subClusterID){

	if(subClusterID < 0 || subClusterID > 49)
		subClusterID = rand()%50;
	int selectID = rand()% moveIndexList[subClusterID].size();
	
	CsMove cm =  oi->index1000ToCsMoveWithRandPerturbation(moveIndexList[subClusterID][selectID]);
	cm.subClusterID = subClusterID;
	return cm;
}

IndividualNuPairMoveSet::~IndividualNuPairMoveSet() {
	// TODO Auto-generated destructor stub
}


int NuPairMoveSetLibrary::dump() {
	// serialized dump method
	string outpath = datapath() + "../binaryCache";
	string fileName = "NuPairMoveSetLibrary";
	if(! makeDirs(outpath)) {
		throw("[Error]Unable to create " + outpath);
	}
	ofstream outs;
	outs.open(outpath + "/" + fileName, ios::out|ios::binary);
	if(!outs.is_open()) {
		throw("[Error]Fail to open " + outpath + "/" + fileName);
	}
	oi->dump(outs);
	outs.write(reinterpret_cast<char*>(nbContactClusterNum), sizeof(int)*16);
	outs.write(reinterpret_cast<char*>(revNbContactClusterNum), sizeof(int)*16);
	int nnbClustNum[16];
	for(int i=0;i<16; i++) {
		nnbClustNum[i] = nnbMoveList[i].size();
	}
	outs.write(reinterpret_cast<char*>(nnbClustNum), sizeof(int)*16);
	for(int i=0; i<16; i++) {
		for(int j=0;j<nbContactClusterNum[i];j++) {
			nbMoveList[i][j]->dump(outs);
		}
		for(int j=0;j<revNbContactClusterNum[i];j++) {
			revNbMoveList[i][j]->dump(outs);
		}
		for(int j=0;j<nnbClustNum[i];j++) {
			nnbMoveList[i][j]->dump(outs);
		}		
	}
	#ifdef DEBUG
		cout << "nbNonContactClusterNum=" << nbNonContactClusterNum << endl;
	#endif
	outs.write(reinterpret_cast<char*>(&nbNonContactClusterNum), sizeof(int));
	char z0[4] = {'\0'};
	outs.write(reinterpret_cast<char*>(z0), sizeof(int));

	for(int i=0;i<nbNonContactClusterNum;i++) {
		nbNonContactMoveList[i]->dump(outs);
		revNbNonContactMoveList[i]->dump(outs);
	}

	// outs.write(reinterpret_cast<char*>(this), sizeof(*this));
	// outs.write(reinterpret_cast<char*>(this->oi), sizeof(*(this->oi)));
	// for(int i=0; i<16; i++) {
	// 	int lv = this->nbMoveList[i].size();
	// 	// outs.write(reinterpret_cast<char*>(&lv), sizeof(int));
	// 	for(int j=0; j<lv; j++) {
	// 		outs.write(reinterpret_cast<char*>(this->nbMoveList[i][j]), sizeof(IndividualNuPairMoveSet));
	// 	}
	// 	int lv = this->revNbMoveList[i].size();
	// 	// outs.write(reinterpret_cast<char*>(&lv), sizeof(int));
	// 	for(int j=0; j<lv; j++) {
	// 		outs.write(reinterpret_cast<char*>(this->revNbMoveList[i][j]), sizeof(IndividualNuPairMoveSet));
	// 	}
	// 	int lv = this->nnbMoveList[i].size();
	// 	// outs.write(reinterpret_cast<char*>(&lv), sizeof(int));
	// 	for(int j=0; j<lv; j++) {
	// 		outs.write(reinterpret_cast<char*>(this->nnbMoveList[i][j]), sizeof(IndividualNuPairMoveSet));
	// 	}
	// }
	outs.close();
	return EXIT_SUCCESS;
}

int NuPairMoveSetLibrary::load() {
	ifstream ins;
	string fileName = datapath() + "../binaryCache/NuPairMoveSetLibrary";
	ins.open(fileName, ios::in|ios::binary);
	if(!ins.is_open()) {
		throw("[Error]Fail to open " + fileName);
	}
	#ifdef DEBUG
		cout << "reading OrientationIndex" << endl;
	#endif
	oi = new OrientationIndex(true);
	oi->load(ins);
	#ifdef DEBUG
		cout <<"reading clusterNums" << endl;
	#endif
	ins.read(reinterpret_cast<char*>(nbContactClusterNum), sizeof(int)*16);
	ins.read(reinterpret_cast<char*>(revNbContactClusterNum), sizeof(int)*16);
	int nnbClustNum[16];
	ins.read(reinterpret_cast<char*>(nnbClustNum), sizeof(int)*16);
	#ifdef DEBUG
		cout << "reading MoveList set 1" << endl;
	#endif
	for(int i=0; i<16; i++) {
		nbMoveList[i].resize(nbContactClusterNum[i]);
		for(int j=0;j<nbContactClusterNum[i];j++) {
			nbMoveList[i][j] = new IndividualNuPairMoveSet();
			nbMoveList[i][j]->load(ins);
		}
		revNbMoveList[i].resize(revNbContactClusterNum[i]);
		for(int j=0;j<revNbContactClusterNum[i];j++) {
			revNbMoveList[i][j] = new IndividualNuPairMoveSet();
			revNbMoveList[i][j]->load(ins);
		}
		nnbMoveList[i].resize(nnbClustNum[i]);
		for(int j=0;j<nnbClustNum[i];j++) {
			nnbMoveList[i][j] = new IndividualNuPairMoveSet();
			nnbMoveList[i][j]->load(ins);
		}
	}
	#ifdef DEBUG
		cout << "reading nbNonContactClusterNum" << endl;
	#endif
	ins.read(reinterpret_cast<char*>(&nbNonContactClusterNum), sizeof(int));
	int z0;
	ins.read(reinterpret_cast<char*>(&z0), sizeof(int));
	#ifdef DEBUG
		cout << "nbNonContactClusterNum=" << nbNonContactClusterNum << endl;
		cout << "reading NonContactMoveList" << endl;
	#endif
	nbNonContactMoveList.resize(nbNonContactClusterNum);
	revNbNonContactMoveList.resize(nbNonContactClusterNum);
	for(int i=0; i<nbNonContactClusterNum; i++) {
		nbNonContactMoveList[i] = new IndividualNuPairMoveSet();
		nbNonContactMoveList[i]->load(ins);
		revNbNonContactMoveList[i] = new IndividualNuPairMoveSet();
		revNbNonContactMoveList[i]->load(ins);
	}
	
	ins.close();	
	return EXIT_SUCCESS;
}

int runIndividualNuPairMoveSetMT(int sep, int pairType, int clusterID, 
OrientationIndex* oi, BinaryBook* bb, vector<IndividualNuPairMoveSet*> * rvec) {
	try {
		rvec->at(clusterID) = new IndividualNuPairMoveSet(sep, pairType, clusterID, oi, bb);
		return EXIT_SUCCESS;
	} catch(exception& e) {
		string errInfo = e.what() + '\n';
		cerr << errInfo;
		return EXIT_FAILURE;
	}
}

NuPairMoveSetLibrary::NuPairMoveSetLibrary(bool withBinary, int binaryMode){



	if(withBinary && binaryMode == 1) {
		// Read from binaryCache
		this->load();
	} else {
		#ifdef TIMING
		clock_t start, end;
		start = clock(); 
		#endif
		BasePairLib* bpLib = new BasePairLib();
		#ifdef TIMING
		double bpLibTime = (double) (clock()-start)/CLOCKS_PER_SEC;
		cout << "bpLib initialization time: " << bpLibTime << " seconds" << endl;
		#endif
		this->oi = new OrientationIndex();
		#ifdef TIMING
		double oiTime = (double) (clock()-start)/CLOCKS_PER_SEC - bpLibTime;
		cout << "oi initialization time: " << oiTime << " seconds, total " <<
		oiTime + bpLibTime << " seconds" << endl;
		#endif

		if(!withBinary) {
			/* Read from TXT*/
			for(int i=0;i<16;i++){
				int typeA = i/4;
				int typeB = i%4;
				int reverseType = typeB*4+typeA;

				int clusterNum = bpLib->nbContactBasePairNum[i];
				int reverseClusterNum = bpLib->nbContactBasePairNum[reverseType];

				for(int j=0;j<clusterNum;j++){
					nbMoveList[i].push_back(new IndividualNuPairMoveSet(1, i, j, oi));
					
				}

				for(int j=0;j<reverseClusterNum;j++){
					revNbMoveList[i].push_back(new IndividualNuPairMoveSet(-1, i, j, oi));
				}

				this->nbContactClusterNum[i] = clusterNum;
				this->revNbContactClusterNum[i] = reverseClusterNum;

				clusterNum = bpLib->nnbBasePairNum[i];
				for(int j=0;j<clusterNum;j++){
					nnbMoveList[i].push_back(new IndividualNuPairMoveSet(2, i, j, oi));
				}
			}

			for(int i=0;i<bpLib->nbNonContactBasePairNum;i++){
				this->nbNonContactMoveList.push_back(new IndividualNuPairMoveSet(1, -1, i, oi));
				this->revNbNonContactMoveList.push_back(new IndividualNuPairMoveSet(-1, -1, i, oi));
			}
			this->nbNonContactClusterNum = bpLib->nbNonContactBasePairNum;

		} else if(binaryMode == 2) {
			/* Read from BinaryTable */
			ifstream ins;
			string FileNb = NSPdataio::datapath()+"../binaryData/pairMove3/nb";
			string FileNnb = NSPdataio::datapath()+"../binaryData/pairMove3/nnb";
 		    ins.open(FileNb,ios::in | ios::binary);
			if(!ins.is_open()) {
            	throw("Unable to open " + FileNb);
        	}
 		    auto* bbNb = new BinaryBook;
 		    bbNb->read(ins);
			ins.close();
			ins.open(FileNnb,ios::in | ios::binary);
			if(!ins.is_open()) {
            	throw("Unable to open " + FileNnb);
        	}
 		    auto* bbNnb = new BinaryBook;
 		    bbNnb->read(ins);
			ins.close();

			char *ntstr=std::getenv("OMP_NUM_THREADS");
			int nt = *ntstr ? atoi(ntstr) : 1;
			shared_ptr<ThreadPool> thrPool(new ThreadPool(nt));
			size_t jid = 0;
			for(int i=0;i<16;i++){
				int clusterNum = bpLib->nbBasePairNum[i];
				nbMoveList[i].resize(clusterNum);
				revNbMoveList[i].resize(clusterNum);
				for(int j=0;j<clusterNum;j++){
					shared_ptr<IntFuncTask> request(new IntFuncTask);
					request->asynBind(runIndividualNuPairMoveSetMT, 1, i, j, oi, bbNb, &nbMoveList[i]);
 	               jid++;
 	               thrPool->addTask(request);
					shared_ptr<IntFuncTask> request2(new IntFuncTask);
					request2->asynBind(runIndividualNuPairMoveSetMT, -1, i, j, oi, bbNb, &revNbMoveList[i]);
 	               jid++;
 	               thrPool->addTask(request2);
					// nbMoveList[i].emplace_back(new IndividualNuPairMoveSet(1, i, j, oi, bbNb));
					// revNbMoveList[i].emplace_back(new IndividualNuPairMoveSet(-1, i, j, oi, bbNb));
				}

				clusterNum = bpLib->nnbBasePairNum[i];
				nnbMoveList[i].resize(clusterNum);
				for(int j=0;j<clusterNum;j++){
					shared_ptr<IntFuncTask> request(new IntFuncTask);
					request->asynBind(runIndividualNuPairMoveSetMT, 2, i, j, oi, bbNnb, &nnbMoveList[i]);
 	               jid++;
 	               thrPool->addTask(request);
					// nnbMoveList[i].emplace_back(new IndividualNuPairMoveSet(2, i, j, oi, bbNnb));
				}
			}
 	       while(true) {
 	           sleep(1);
 	           if(thrPool->getTaskCount() == 0) {
 	               break;
 	           }
 	       }
			delete bbNb;
			delete bbNnb;
		}
		#ifdef TIMING
		double nuPairMoveSetTime = (double) (clock()-start)/CLOCKS_PER_SEC - bpLibTime - oiTime;
		cout << "NuPairMoveSet initialization time: " << nuPairMoveSetTime << " seconds, total " <<
		oiTime + bpLibTime + nuPairMoveSetTime << " seconds" << endl;
		#endif

		delete bpLib;
	}
}

void NuPairMoveSetLibrary::printMoveLibInfo(){
	cout << "neighbor move: " << endl;
	for(int i=0;i<16;i++){
		int imNum = nbMoveList[i].size();
		cout << "cluster num: " << i << " " << imNum << endl;
		for(int j=0;j<imNum;j++){
			int moveNum = 0;
			for(int k=0;k<50;k++){
				moveNum += nbMoveList[i][j]->moveIndexList[k].size();
			}
			cout << "cluster: " << j << " moveNum: " << moveNum << endl;
		}
	}

	cout << "reverse neighbor move: " << endl;
	for(int i=0;i<16;i++){
		int imNum = revNbMoveList[i].size();
		cout << "cluster num: " << i << " " << imNum << endl;
		for(int j=0;j<imNum;j++){
			int moveNum = 0;
			for(int k=0;k<50;k++){
				moveNum += revNbMoveList[i][j]->moveIndexList[k].size();
			}
			cout << "cluster: " << j << " moveNum: " << moveNum << endl;
		}
	}

	cout << "non-contact neighbor move: " << endl;
	int num = nbNonContactMoveList.size();
	cout << "cluster num: " << num << endl;
	for(int j=0;j<num;j++){
		int moveNum = 0;
		for(int k=0;k<50;k++){
			moveNum += nbNonContactMoveList[j]->moveIndexList[k].size();
		}
		cout << "cluster: " << j << " moveNum: " << moveNum << endl;
	}

	cout << "non-contact reverse neighbor move: " << endl;
	num =  revNbNonContactMoveList.size();
	cout << "cluster num: " << num << endl;
	for(int j=0;j<num;j++){
		int moveNum = 0;
		for(int k=0;k<50;k++){
			moveNum += revNbNonContactMoveList[j]->moveIndexList[k].size();
		}
		cout << "cluster: " << j << " moveNum: " << moveNum << endl;
	}

	cout << "non-neighbor move:" << endl;
	for(int i=0;i<16;i++){
		int imNum = nnbMoveList[i].size();
		cout << "cluster num: " << i << " " << imNum << endl;
		for(int j=0;j<imNum;j++){
			int moveNum = 0;
			for(int k=0;k<50;k++){
				moveNum += nnbMoveList[i][j]->moveIndexList[k].size();
			}
			cout << "cluster: " << j << " moveNum: " << moveNum << endl;
		}
	}	
}

IndividualNuPairMoveSet* NuPairMoveSetLibrary::getMoveSet(int type, int clusterID, int sep){
	if(type < 0 || type >= 16){
		cout << "invalid move set base type: " << type << endl;
		exit(0);
	}

	if(sep==1){
		if(clusterID < 0) {
			cout << "invalid move set clusterID: " << clusterID << " pairType: " << type << endl;
			exit(0);
		}
		else if(clusterID < this->nbContactClusterNum[type]){
			return this->nbMoveList[type][clusterID];
		}
		else if(clusterID < this->nbContactClusterNum[type] + this->nbNonContactClusterNum) {
			return this->nbNonContactMoveList[clusterID - this->nbContactClusterNum[type]];
		}
		else {
			cout << "invalid move set clusterID: " << clusterID << " pairType: " << type << endl;
			exit(0);
		}
	}
	else if(sep == -1){


		if(clusterID < 0) {
			cout << "invalid move set clusterID: " << clusterID << " pairType: " << type << endl;
			exit(0);
		}
		else if(clusterID < this->revNbContactClusterNum[type]){
			return this->revNbMoveList[type][clusterID];
		}
		else if(clusterID < this->revNbContactClusterNum[type] + this->nbNonContactClusterNum) {
			return this->revNbNonContactMoveList[clusterID - this->revNbContactClusterNum[type]];
		}
		else {
			cout << "invalid move set clusterID: " << clusterID << " pairType: " << type << endl;
			exit(0);
		}		
	}
	else{
		return this->nnbMoveList[type][clusterID];
	}
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
	delete oi;
}

MixedNuPairCluster::MixedNuPairCluster(int sep, int pairType, NuPairMoveSetLibrary* lib){
	this->sep = sep;
	this->pairType = pairType;
	this->moveLib = lib;
	this->fixedNativeCM = false;
	this->type = "unknown";
}

void MixedNuPairCluster::updateEdgeInformation(EdgeInformation* ei){

	double psum = 0;
	this->type = ei->ssSepKey;
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
	for(int i=0;i<100000;i++){
		randPool[i] = 0;
	}

	for(int i=0;i<clusterIDList.size();i++){
		int start = (int)(pAdd*100000);
		pAdd += clusterPList[i]/psum;
		int end = (int)(pAdd*100000);

		if(end > 100000){
			cout << "pAdd " << pAdd << endl;
			exit(0);
		}

		for(int j=start;j<end;j++){
			randPool[j] = clusterIDList[i];
		}
	}
}

CsMove MixedNuPairCluster::getRandomMove(){

	if(fixedNativeCM){
		return natCM;
	}

	int cluster = randPool[rand()%100000];
	IndividualNuPairMoveSet* selectMoveSet = moveLib->getMoveSet(pairType, cluster, sep);

	CsMove cm = moveLib->getMoveSet(pairType, cluster, sep)->getRandomMove(moveLib->oi);

	cm.clusterID = cluster;
	return cm;
}

CsMove MixedNuPairCluster::getRandomMoveWithFixedCluster(CsMove& move){
	if(fixedNativeCM){
		return natCM;
	}

	int cluster = move.clusterID;
	if(cluster < 0)
		cluster = randPool[rand()%100000];
	

	CsMove cm = moveLib->getMoveSet(pairType, cluster, sep)->getRandomMove(moveLib->oi);
	cm.clusterID = cluster;
	return cm;
}

CsMove MixedNuPairCluster::getRandomMoveWithFixedSubCluster(CsMove& move){
	if(fixedNativeCM){
		return natCM;
	}

	int cluster = move.clusterID;
	if(cluster < 0)
		cluster = randPool[rand()%100000];
	int subCluster = move.subClusterID;
	
	CsMove cm = moveLib->getMoveSet(pairType, cluster, sep)->getRandomMoveWithFixedSubCluster(moveLib->oi, subCluster);
	cm.clusterID = cluster;
	return cm;
	
}

CsMove MixedNuPairCluster::getRandomMoveWithFixedSP1000Index(CsMove& move) {
	return moveLib->oi->fixIndex1000WithRandPerturbation(move);
}

void MixedNuPairCluster::printMoveSetInfo(){

	cout << "move set info: " << endl;
	cout << this->type << endl;
	cout << this->pairType << endl;
	cout << this->clusterIDList.size() << endl;

	if(this->type == "all") {
		cout << "all" << endl;
		return;
	}

	set<int> clusterSet;
	for(int i=0;i<100000;i++){
		clusterSet.insert(randPool[i]);
	}

	cout << "cluster num: " << clusterSet.size() << endl;

	set<int>::iterator it;
	for(it=clusterSet.begin();it!=clusterSet.end();it++){
		cout << *it << " ";
	}
	cout << endl;

}

MixedNuPairCluster::~MixedNuPairCluster(){

}


} /* namespace NSPpredNA */
