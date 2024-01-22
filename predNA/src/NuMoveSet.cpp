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

	if(!bb){
		string path = NSPdataio::datapath()+"pairMove3/";

		if(sep == 1 || sep == -1)
			fileName = path+"nb/"+string(xx)+".move";
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
				printf("%s %d\n", fileName, moveID);
			}


			if(sep > 0)
				moveIndexList[binID].push_back(moveID);
			else if(sep == -1){
				CsMove move = oi->index1000ToCsMove(moveID);
				CsMove revMove = move.reverse();
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
				// CsMove move = oi->index1000ToCsMove(moveID);
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
	int vecSize[20];
	for(int i=0;i<20;i++) {
		vecSize[i] = moveIndexList[i].size();
	}
	outs.write(reinterpret_cast<char *>(vecSize), 20*sizeof(int));
	for(int i=0;i<20;i++) {
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
	int vecSize[20];
	ins.read(reinterpret_cast<char*>(vecSize), sizeof(int)*20);
	for(int i=0;i<20;i++) {
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
	int selectBin = rand()% 20;
	int selectID = rand()% moveIndexList[selectBin].size();
	CsMove cm =  oi->index1000ToCsMoveWithRandPerturbation(moveIndexList[selectBin][selectID]);
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
	for(int i=0; i<16; i++) {
		int lnb = nbMoveList[i].size();
		outs.write(reinterpret_cast<char*>(&lnb), sizeof(int));
		int lnnb = nnbMoveList[i].size();
		outs.write(reinterpret_cast<char*>(&lnnb), sizeof(int));
		for(int j=0;j<lnb;j++) {
			nbMoveList[i][j]->dump(outs);
		}
		for(int j=0;j<lnb;j++) {
			revNbMoveList[i][j]->dump(outs);
		}
		for(int j=0;j<lnnb;j++) {
			nnbMoveList[i][j]->dump(outs);
		}
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
	oi = new OrientationIndex(true);
	oi->load(ins);
	for(int i=0; i<16; i++) {
		int lnb, lnnb, tmp[2];
		ins.read(reinterpret_cast<char*>(&tmp), 2*sizeof(int));
		lnb = tmp[0];
		lnnb = tmp[1];
		nbMoveList[i].reserve(lnb);
		for(int j=0;j<lnb;j++) {
			nbMoveList[i][j] = new IndividualNuPairMoveSet();
			nbMoveList[i][j]->load(ins);
		}
		revNbMoveList[i].reserve(lnb);
		for(int j=0;j<lnb;j++) {
			revNbMoveList[i][j] = new IndividualNuPairMoveSet();
			revNbMoveList[i][j]->load(ins);
		}
		nnbMoveList[i].reserve(lnnb);
		for(int j=0;j<lnnb;j++) {
			nnbMoveList[i][j] = new IndividualNuPairMoveSet();
			nnbMoveList[i][j]->load(ins);
		}
	}

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
			for(int k=0;k<20;k++){
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
			for(int k=0;k<20;k++){
				moveNum += revNbMoveList[i][j]->moveIndexList[k].size();
			}
			cout << "cluster: " << j << " moveNum: " << moveNum << endl;
		}
	}

	cout << "non-neighbor move:" << endl;
	for(int i=0;i<16;i++){
		int imNum = nnbMoveList[i].size();
		cout << "cluster num: " << i << " " << imNum << endl;
		for(int j=0;j<imNum;j++){
			int moveNum = 0;
			for(int k=0;k<20;k++){
				moveNum += nnbMoveList[i][j]->moveIndexList[k].size();
			}
			cout << "cluster: " << j << " moveNum: " << moveNum << endl;
		}
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
	this->type = ei->moveType;

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

	if(sep == 1)
		return moveLib->nbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);
	else if(sep == -1) {
		return moveLib->revNbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);
	}
	else
		return moveLib->nnbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);

}

void MixedNuPairCluster::printMoveSetInfo(){

	if(this->type == "all") {
		cout << "all" << endl;
		return;
	}

	set<int> clusterSet;
	for(int i=0;i<100000;i++){
		clusterSet.insert(randPool[i]);
	}
	set<int>::iterator it;
	for(it=clusterSet.begin();it!=clusterSet.end();it++){
		cout << *it << " ";
	}
	cout << endl;

}

MixedNuPairCluster::~MixedNuPairCluster(){

}


} /* namespace NSPpredNA */
