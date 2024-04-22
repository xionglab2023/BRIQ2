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

	int stepNum1 = (int)(this->totalSamp * 5000);
	int stepNum2 = (int)(this->totalSamp * 500);
	int stepNum3 = (int)(this->totalSamp * 500);

	cout << "stepNum: " << stepNum1+stepNum2+stepNum3 << endl;
	
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
	int roundNum = 10;
	map<string, double>::iterator it;
	char xx[200];
	for(int round=0;round<roundNum;round++) {

		graph->initRandWeight();
		graph->MST_kruskal(tree);
		tree->randomInitCG();
		tree->updateNodeInfoCG();
		tree->updateEdgeInfoCG();
		tree->updateSamplingInfo();
		tree->printEdges();
		curEne = graph->totalEnergyCG();
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
					randNode->updateRiboseRotamerCG(randRot);
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

					randEdge->updateCsMoveCG(randMove);

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
			double totEne = graph->totalEnergyCG();
			graphInfo* gi = graph->getGraphInfoCG();
			double rms = gi->rmsdCG(this->graph->initInfo);
			delete gi;
			printf("T=%7.4f nTot=%7d pN=%6.4f eTot=%7d pE=%6.4f curE=%8.3f totEne=%8.3f rms: %6.3f\n", T, nTot, nAc*1.0/nTot, eTot, eAc*1.0/eTot, curEne, totEne, rms);
		}

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
					randNode->updateRiboseRotamerCG(randRot);
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
					randMove = randEdge->moveSet->getRandomMoveWithFixedSubCluster(randEdge->cm);

					randEdge->updateCsMoveCG(randMove);

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
			double totEne = graph->totalEnergyCG();
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
					randNode->updateRiboseRotamerCG(randRot);
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
					randEdge->updateCsMoveCG(randMove);
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
			double totEne = graph->totalEnergyCG();
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