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

void NuSampling::runCoarseGrainedMC(mixedContactInfo* out){
	bool debug = false;
	//randomInitCG();

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
				randEdge = tree->geList[randPos];
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

    
}

}