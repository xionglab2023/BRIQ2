/*
 * NuMoveSet.cpp
 *
 *  Created on: 2023��11��23��
 *      Author: nuc
 */

#include <predNA/NuMoveSet.h>
#include "tools/ThreadPool.h"
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
		string path = NSPdataio::datapath()+"pairMove2/";

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

CsMove IndividualNuPairMoveSet::getRandomMove(OrientationIndex* oi){
	int selectBin = rand()% 20;
	int selectID = rand()% moveIndexList[selectBin].size();
	CsMove cm =  oi->index1000ToCsMove(moveIndexList[selectBin][selectID]);
	return cm;
}

IndividualNuPairMoveSet::~IndividualNuPairMoveSet() {
	// TODO Auto-generated destructor stub
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

NuPairMoveSetLibrary::NuPairMoveSetLibrary(bool withBinary){

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
	} else {
		/* Read from Binary */
		ifstream ins;
		string FileNb = NSPdataio::datapath()+"../binaryData/pairMove2/nb";
		string FileNnb = NSPdataio::datapath()+"../binaryData/pairMove2/nnb";
 	    ins.open(FileNb,ios::in | ios::binary);
 	    auto* bbNb = new BinaryBook;
 	    bbNb->read(ins);
		ins.close();
		ins.open(FileNnb,ios::in | ios::binary);
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

	if(fixedNativeCM){
		return natCM;
	}

	int cluster = randPool[rand()%10000];

	if(sep == 1)
		return moveLib->nbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);
	else if(sep == -1)
		return moveLib->revNbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);
	else
		return moveLib->nnbMoveList[pairType][cluster]->getRandomMove(moveLib->oi);

}

void MixedNuPairCluster::printMoveSetInfo(){
	set<int> clusterSet;
	for(int i=0;i<10000;i++){
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
