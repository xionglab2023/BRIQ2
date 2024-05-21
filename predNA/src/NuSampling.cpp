/*
 * NuSampling.cpp
 *
 *  Created on: 2024,3,6
 *      Author: pengx
 */

#include <predNA/NuSampling.h>

namespace NSPpredNA {

nuContactMatrix::nuContactMatrix(int seqLen){
    this->clusterIDMtx = new int[seqLen*seqLen];
}

nuContactMatrix::~nuContactMatrix(){
    delete [] clusterIDMtx;
}

singleContactInfo::singleContactInfo(int seqLen, int* seq, bool* connectToDownstream, bool* fixed, nuContactMatrix* scm){
    this->seqLen = seqLen;
    this->seq = new int[seqLen];
    this->connectToDownstream = new bool[seqLen];
    this->fixed = new bool[seqLen];

    for(int i=0;i<seqLen;i++){
        this->seq[i] = seq[i];
        this->connectToDownstream[i] = connectToDownstream[i];
        this->fixed[i] = fixed[i];
    }
    this->scm = scm;
}

singleContactInfo::~singleContactInfo(){
    delete [] seq;
    delete [] connectToDownstream;
    delete [] fixed;
}

mixedContactInfo::mixedContactInfo(int seqLen, int* seq, bool* connectToDownstream, bool* fixed){
    this->seqLen = seqLen;
    this->seq = new int[seqLen];
    this->connectToDownstream = new bool[seqLen];
    this->fixed = new bool[seqLen];

    for(int i=0;i<seqLen;i++){
        this->seq[i] = seq[i];
        this->connectToDownstream[i] = connectToDownstream[i];
        this->fixed[i] = fixed[i];
    }
}

mixedContactInfo::~mixedContactInfo(){
    delete [] seq;
    delete [] connectToDownstream;
    delete [] fixed;
}

double mixedContactInfo::totalConfEntropy(){
    return 0.0;
}

double mixedContactInfo::relativeConfEntropy(singleContactInfo* natConf){
    return 0.0;
}

NuSampling::NuSampling(NuGraph* graph, NuTree* tree){
    this->graph = graph;
    this->tree = tree;

    this->poolSize = tree->poolSize;
    this->sampFreqEdge = tree->sampFreqEdge;
    this->sampFreqNode = tree->sampFreqNode;
    this->totalSamp = tree->totalSamp;
    for(int i=0;i<poolSize;i++){
        randPoolNode[i] = tree->randPoolNode[i];
        randPoolEdge[i] = tree->randPoolEdge[i];
    }
}

void NuSampling::runCoarseGrainedMC(map<string,double>& results, const string& outFile){
	bool debug = false;

	double T0 = 10.0;
	double T1 = 0.1;
	double T2 = 0.01;
	double T3 = 0.001;

	int stepNum1 = (int)(this->totalSamp * graph->et->para->kStepNum1CG);
	int stepNum2 = (int)(this->totalSamp * graph->et->para->kStepNum2CG);
	int stepNum3 = (int)(this->totalSamp * graph->et->para->kStepNum3CG);

	cout << "stepNum: " << stepNum1+stepNum2+stepNum3 << endl;
	
	double anneal = 0.95;

	double curEne = graph->totalEnergyCG(1.0, 1.0);
	double lastEne = curEne;
	int i,j,k, randPos, nAc, eAc, nTot, eTot;
	double T, randP, mutE;
	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamerCG* randRot;                     
	CsMove randMove;

	double initClashRescale = 0.03;
	double initConnectRescale = 0.1;
	double lamda = 1.07;

	int len = graph->seqLen;
	int count = 0;
	int roundNum = 10;
	map<string, double>::iterator it;
	char xx[200];
	for(int round=0;round<roundNum;round++) {
		double clashRescale = initClashRescale;
		double connectRescale = initConnectRescale;

		cout << "init rand weight" << endl;
		graph->initRandWeight();
		cout << "init mst" << endl;
		graph->MST_kruskal(tree);
		cout << "update node info" << endl;
		tree->updateNodeInfoCG(clashRescale, connectRescale);
		cout << "update edge info" << endl;
		tree->updateEdgeInfoCG(clashRescale, connectRescale);
		cout << "update sampling info" << endl;
		tree->updateSamplingInfo();
		cout << "rand init cg" << endl;
		tree->randomInitCG(clashRescale, connectRescale);
		cout << "print edge" << endl;
		tree->printEdges();
		curEne = graph->totalEnergyCG(clashRescale, connectRescale);
		lastEne = curEne;

		sprintf(xx, "%d", round);
		string outFile_round = outFile.substr(0, outFile.length()-4) + "-" + string(xx) + ".pdb";

		for(T=T0;T>T1;T=T*anneal){
			nAc = 0;
			eAc = 0;
			nTot = 0;
			eTot = 0;

			for(k=0;k<stepNum1;k++){
				count++;

				randP = rand()*1.0/RAND_MAX;
				if(randP < tree->sampFreqNode){
					/*
					 * rotamer mut
					 */
					nTot ++;
					randPos = randPoolNode[rand()%poolSize];

					randNode = graph->allNodes[randPos];
					randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
					randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
					mutE = randNode->rotMutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randNode->acceptRotMutationCG();
						curEne += mutE;
						nAc++;
					}
					else {
						randNode->clearRotMutationCG();					
					}
				}
				else {
					/*
					 * edge mut
					 */
					eTot ++;
					randPos = randPoolEdge[rand()%poolSize];
					randEdge = tree->geList[randPos];
					randMove = randEdge->moveSet->getRandomMove();

					randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale);

					mutE = randEdge->mutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randEdge->acceptMutationCG();
						curEne += mutE;
						eAc++;					
					}
					else {
						randEdge->clearMutationCG();
					}
				}
			}
			double totEne = graph->totalEnergyCG(clashRescale, connectRescale);
			graphInfo* gi = graph->getGraphInfoCG();
			double rms = gi->rmsdCG(this->graph->initInfo);
			delete gi;
			printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f clashRescale: %5.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms, clashRescale);

			clashRescale = clashRescale * lamda;
			connectRescale = connectRescale * lamda;
			if(clashRescale > 1.0) clashRescale = 1.0;
			if(connectRescale > 1.0) connectRescale = 1.0;
			curEne = graph->totalEnergyCG(clashRescale, connectRescale);
			graph->updateEnergyCG(clashRescale, connectRescale);
		}


		clashRescale = 1.0;
		connectRescale = 1.0;
		curEne = graph->totalEnergyCG(clashRescale, connectRescale);
		graph->updateEnergyCG(clashRescale, connectRescale);

		for(T=T1;T>T2;T=T*anneal){
			nAc = 0;
			eAc = 0;
			nTot = 0;
			eTot = 0;

			for(k=0;k<stepNum2;k++){
				count++;

				randP = rand()*1.0/RAND_MAX;
				if(randP < tree->sampFreqNode){
					/*
					 * rotamer mut
					 */
					nTot ++;
					randPos = randPoolNode[rand()%poolSize];

					randNode = graph->allNodes[randPos];
					randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
					randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
					mutE = randNode->rotMutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randNode->acceptRotMutationCG();
						curEne += mutE;
						nAc++;
					}
					else {
						randNode->clearRotMutationCG();					
					}
				}
				else {
					/*
					 * edge mut
					 */
					eTot ++;
					//cout << "rand pos" << endl;
					randPos = randPoolEdge[rand()%poolSize];
					//cout << "rand edge" << endl;
					randEdge = tree->geList[randPos];
					//if(randEdge->cm.clusterID < 0)
					//   cout << "rand move: " << randEdge->nodeA->baseType << " " << randEdge->nodeB->baseType << " " << randEdge->sep << " " << randEdge->cm.clusterID  << endl;
					randMove = randEdge->moveSet->getRandomMoveWithFixedSubCluster(randEdge->cm);
					//cout << "update edge" << endl;
					randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale);

					mutE = randEdge->mutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randEdge->acceptMutationCG();
						curEne += mutE;
						eAc++;					
					}
					else {
						randEdge->clearMutationCG();
					}
				}
			}
			double totEne = graph->totalEnergyCG(clashRescale, connectRescale);
			graphInfo* gi = graph->getGraphInfoCG();
			double rms = gi->rmsdCG(this->graph->initInfo);
			delete gi;
			printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
		}

		for(T=T2;T>T3;T=T*anneal){
			nAc = 0;
			eAc = 0;
			nTot = 0;
			eTot = 0;

			for(k=0;k<stepNum3;k++){
				count++;

				randP = rand()*1.0/RAND_MAX;
				if(randP < tree->sampFreqNode){
					/*
					 * rotamer mut
					 */
					nTot ++;
					randPos = randPoolNode[rand()%poolSize];

					randNode = graph->allNodes[randPos];
					randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
					randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
					mutE = randNode->rotMutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randNode->acceptRotMutationCG();
						curEne += mutE;
						nAc++;
					}
					else {
						randNode->clearRotMutationCG();					
					}
				}
				else {
					/*
					 * edge mut
					 */
					eTot ++;
					randPos = randPoolEdge[rand()%poolSize];
					randEdge = tree->geList[randPos];
					randMove = randEdge->moveSet->getRandomMoveWithFixedSP1000Index(randEdge->cm);
					randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale);
					mutE = randEdge->mutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randEdge->acceptMutationCG();
						curEne += mutE;
						eAc++;					
					}
					else {
						randEdge->clearMutationCG();
					}
				}
			}
			double totEne = graph->totalEnergyCG(clashRescale, connectRescale);
			graphInfo* gi = graph->getGraphInfoCG();
			double rms = gi->rmsdCG(this->graph->initInfo);
			delete gi;
			printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
		}

		cout << "to cg model" << endl;
		graph->cgToAllAtom();
		graphInfo* gi = graph->getGraphInfo(0.0);
		gi->printPDB(outFile_round);
		delete gi;

		string key = graph->toContactMapHashKeyCG();
		cout << key << endl;
		graph->keyToContactMatrix(key);

		it = results.find(key);
		if(it != results.end()){
			if(curEne < results[key])
				results[key] = curEne;
		}
		else {
			results[key] = curEne;
		}
	}

	cout << "mc result key num: " << results.size() << endl;
}

void NuSampling::runCoarseGrainedMC(map<string,double>& results, int roundNum){

	bool debug = false;

	double T0 = 10.0;
	double T1 = 0.1;
	double T2 = 0.01;
	double T3 = 0.001;

	int stepNum1 = (int)(this->totalSamp * graph->et->para->kStepNum1CG);
	int stepNum2 = (int)(this->totalSamp * graph->et->para->kStepNum2CG);
	int stepNum3 = (int)(this->totalSamp * graph->et->para->kStepNum3CG);

	if(debug) {
		stepNum1 = 10;
		stepNum2 = 10;
		stepNum3 =10;
	}


	double anneal = 0.95;

	double curEne = graph->totalEnergyCG(1.0, 1.0);
	double lastEne = curEne;
	int i,j,k, randPos, nAc, eAc, nTot, eTot;
	double T, randP, mutE;
	NuNode* randNode;
	NuEdge* randEdge;
	RiboseRotamerCG* randRot;                     
	CsMove randMove;

	double initClashRescale = 0.03;
	double initConnectRescale = 0.1;
	double lamda = 1.07;


	int len = graph->seqLen;
	int count = 0;

	map<string, double>::iterator it;
	char xx[200];
	for(int round=0;round<roundNum;round++) {
		double clashRescale = initClashRescale;
		double connectRescale = initConnectRescale;

		graph->updateEnergyCG(clashRescale, connectRescale);
		curEne = graph->totalEnergyCG(clashRescale, connectRescale);
		
		graph->initRandWeight();
		graph->MST_kruskal(tree);
		tree->updateNodeInfoCG(clashRescale, connectRescale);
		tree->updateEdgeInfoCG(clashRescale, connectRescale);
		tree->updateSamplingInfo();
		tree->randomInitCG(clashRescale, connectRescale);
		tree->printEdges();
		
		graph->updateEnergyCG(clashRescale, connectRescale);
		curEne = graph->totalEnergyCG(clashRescale, connectRescale);
		

		lastEne = curEne;
		
		for(T=T0;T>T1;T=T*anneal){
			nAc = 0;
			eAc = 0;
			nTot = 0;
			eTot = 0;

			for(k=0;k<stepNum1;k++){
				count++;

				randP = rand()*1.0/RAND_MAX;
				if(randP < tree->sampFreqNode){
					/*
					 * rotamer mut
					 */
					nTot ++;
					randPos = randPoolNode[rand()%poolSize];

					randNode = graph->allNodes[randPos];
					randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
					randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
					mutE = randNode->rotMutEnergyCG();

					double oldEne = curEne;
				

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randNode->acceptRotMutationCG();
						curEne += mutE;
						nAc++;

						if(abs(curEne - graph->totalEnergyCG(clashRescale, connectRescale)) > 0.01) {
							cout << "rotamer mut: " << randPos << " k=" << k << " oldEne: " << oldEne << " mutE: " << mutE << " totalE: " << graph->totalEnergyCG(clashRescale, connectRescale) << endl;
							exit(0);
						}
					}
					else {
						randNode->clearRotMutationCG();					
					}
				}
				else {
					/*
					 * edge mut
					 */
					eTot ++;
					randPos = randPoolEdge[rand()%poolSize];
					randEdge = tree->geList[randPos];

					if(debug) {
						cout << "edge mut: " << randEdge->indexA << "-" << randEdge->indexB << endl;
					}

					randMove = randEdge->moveSet->getRandomMove();

					if(debug) {
						cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before update cs" << endl;
						graph->checkEnergyCG(clashRescale, connectRescale);
						cout << "move: " << endl;
						randMove.print();
					}

					randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale);
					
					
					mutE = randEdge->mutEnergyCG();

					if(debug) {
						cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " before mutE" << endl;
						graph->checkEnergyCG(clashRescale, connectRescale);
						graph->printEnergyCG(clashRescale);
					}

					double oldEne = curEne;

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randEdge->acceptMutationCG();
						curEne += mutE;
						eAc++;
					}
					else {
						randEdge->clearMutationCG();
					}

					if(debug) {
						cout << "edge mut, edge: " << randEdge->indexA << " " << randEdge->indexB << " AC/RJ" << endl;
						graph->checkEnergyCG(clashRescale, connectRescale);
					}
				}
			}
			double totEne = graph->totalEnergyCG(clashRescale, connectRescale);
			graphInfo* gi = graph->getGraphInfoCG();
			double rms = gi->rmsdCG(this->graph->initInfo);
			delete gi;
			printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f clashRescale: %5.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms, clashRescale);

			clashRescale = clashRescale * lamda;
			connectRescale = connectRescale * lamda;

			if(clashRescale > 1.0) clashRescale = 1.0;
			if(connectRescale > 1.0) connectRescale = 1.0;

			
			graph->updateEnergyCG(clashRescale, connectRescale);
			curEne = graph->totalEnergyCG(clashRescale, connectRescale);
		}


		clashRescale = 1.0;
		connectRescale = 1.0;
		curEne = graph->totalEnergyCG(clashRescale, connectRescale);
		graph->updateEnergyCG(clashRescale, connectRescale);

		for(T=T1;T>T2;T=T*anneal){
			nAc = 0;
			eAc = 0;
			nTot = 0;
			eTot = 0;

			for(k=0;k<stepNum2;k++){
				count++;

				randP = rand()*1.0/RAND_MAX;
				if(randP < tree->sampFreqNode){
					/*
					 * rotamer mut
					 */
					nTot ++;
					randPos = randPoolNode[rand()%poolSize];

					randNode = graph->allNodes[randPos];
					randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
					randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
					mutE = randNode->rotMutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randNode->acceptRotMutationCG();
						curEne += mutE;
						nAc++;
					}
					else {
						randNode->clearRotMutationCG();					
					}
				}
				else {
					/*
					 * edge mut
					 */
					eTot ++;
					//cout << "rand pos" << endl;
					randPos = randPoolEdge[rand()%poolSize];
					//cout << "rand edge" << endl;
					randEdge = tree->geList[randPos];
					//if(randEdge->cm.clusterID < 0)
					//   cout << "rand move: " << randEdge->nodeA->baseType << " " << randEdge->nodeB->baseType << " " << randEdge->sep << " " << randEdge->cm.clusterID  << endl;
					randMove = randEdge->moveSet->getRandomMoveWithFixedSubCluster(randEdge->cm);
					//cout << "update edge" << endl;
					randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale);

					mutE = randEdge->mutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randEdge->acceptMutationCG();
						curEne += mutE;
						eAc++;					
					}
					else {
						randEdge->clearMutationCG();
					}
				}
			}
			double totEne = graph->totalEnergyCG(clashRescale, connectRescale);
			graphInfo* gi = graph->getGraphInfoCG();
			double rms = gi->rmsdCG(this->graph->initInfo);
			delete gi;
			printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
		}

		for(T=T2;T>T3;T=T*anneal){
			nAc = 0;
			eAc = 0;
			nTot = 0;
			eTot = 0;

			for(k=0;k<stepNum3;k++){
				count++;

				randP = rand()*1.0/RAND_MAX;
				if(randP < tree->sampFreqNode){
					/*
					 * rotamer mut
					 */
					nTot ++;
					randPos = randPoolNode[rand()%poolSize];

					randNode = graph->allNodes[randPos];
					randRot = graph->rotLib->riboseRotLib->getRandomRotamerCG(randNode->baseType);
					randNode->updateRiboseRotamerCG(randRot, clashRescale, connectRescale);
					mutE = randNode->rotMutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randNode->acceptRotMutationCG();
						curEne += mutE;
						nAc++;
					}
					else {
						randNode->clearRotMutationCG();					
					}
				}
				else {
					/*
					 * edge mut
					 */
					eTot ++;
					randPos = randPoolEdge[rand()%poolSize];
					randEdge = tree->geList[randPos];
					randMove = randEdge->moveSet->getRandomMoveWithFixedSP1000Index(randEdge->cm);
					randEdge->updateCsMoveCG(randMove, clashRescale, connectRescale);
					mutE = randEdge->mutEnergyCG();

					if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
						randEdge->acceptMutationCG();
						curEne += mutE;
						eAc++;					
					}
					else {
						randEdge->clearMutationCG();
					}
				}
			}
			double totEne = graph->totalEnergyCG(clashRescale, connectRescale);
			graphInfo* gi = graph->getGraphInfoCG();
			double rms = gi->rmsdCG(this->graph->initInfo);
			delete gi;
			printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
		}

		cout << "to cg model" << endl;
		graph->cgToAllAtom();
		string key = graph->toContactMapHashKeyCG();
		cout << key << endl;
		graph->keyToContactMatrix(key);

		it = results.find(key);
		if(it != results.end()){
			if(curEne < results[key])
				results[key] = curEne;
		}
		else {
			results[key] = curEne;
		}
	}

	cout << "mc result key num: " << results.size() << endl;
}

}