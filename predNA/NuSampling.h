#ifndef PREDNA_NUSAMPLING_H_
#define PREDNA_NUSAMPLING_H_

#include <vector>
#include <string>
#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include "predNA/NuGraph.h"

namespace NSPpredNA {

class nuContactMatrix{
    
public:
    int* clusterIDMtx; //seqLen*seqLen
    
    nuContactMatrix(int seqLen);

    virtual ~nuContactMatrix();
};

class singleContactInfo{

public:

    int seqLen;
	int* seq;
	bool* connectToDownstream;
	bool* fixed;
    nuContactMatrix* scm;

    singleContactInfo(int seqLen, int* seq, bool* connectToDownstream, bool* fixed, nuContactMatrix* scm);
    virtual ~singleContactInfo();
};

class mixedContactInfo{

public:

    int seqLen;
	int* seq;
	bool* connectToDownstream;
	bool* fixed;

    vector<nuContactMatrix*> scmList;

    mixedContactInfo(int seqLen, int* seq, bool* connectToDownstream, bool* fixed);
    double totalConfEntropy();
    double relativeConfEntropy(singleContactInfo* natConf);
    virtual ~mixedContactInfo();
};

class NuSampling{

public:

    NuGraph* graph;
    NuTree* tree;

    int poolSize = 100000;
	double sampFreqNode;
	int randPoolNode[100000];
	double sampFreqEdge;
	int randPoolEdge[100000];
	double totalSamp;

    NuSampling(NuGraph* graph, NuTree* tree);

	void runCoarseGrainedMC(map<string, double>& results, int roundNum);
    void runCoarseGrainedMC(const string& outFile, int roundNum);
};


}

#endif