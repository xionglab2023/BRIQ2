/*
 * NuMoveSet.h
 *
 *  Created on: 2023��11��23��
 *      Author: nuc
 */

#ifndef PREDNA_NUMOVESET_H_
#define PREDNA_NUMOVESET_H_
// #define TIMING

#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "dataio/datapaths.h"
#include "dataio/binaryTable.h"
#include "geometry/CsMove.h"
#include "model/BasePairLib.h"
#include "predNA/EdgeInformation.h"
#include "geometry/OrientationIndex.h"

namespace NSPpredNA {

using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPdataio;

class IndividualNuPairMoveSet {
public:

	int sep; //1: nb   -1: revNb  2: nnb
	int pairType; //0~15: AA, AU, AG, AC, UA, ..., CC
	int clusterID;

	/*
	 * Ust SP1000 move index, distance 50 bins, dihedral angle 40 bins, sphere 1000*1000
	 * move set were divided into 20 groups, the total probability of each group is 5%
	 */
	vector<int> moveIndexList[20];

	/**
	 * @brief Construct an empty IndividualNuPairMoveSet object 
	 * in advance of object loading
	 * 
	 */
	IndividualNuPairMoveSet() {
	};

	IndividualNuPairMoveSet(int sep, int pairType, int clusterID, OrientationIndex* oi, BinaryBook* bb=nullptr);
	/**
	 * @brief dump object to binary cache
	 * 
	 * @return int 
	 */
	int dump();

	/**
	 * @brief load object from binary cache
	 * 
	 * @return int 
	 */
	int load();

	CsMove getRandomMove(OrientationIndex* oi);

	virtual ~IndividualNuPairMoveSet();
};

class NuPairMoveSetLibrary {
public:
	vector<IndividualNuPairMoveSet*> nbMoveList[16];
	vector<IndividualNuPairMoveSet*> revNbMoveList[16];
	vector<IndividualNuPairMoveSet*> nnbMoveList[16];
	OrientationIndex* oi;

	/**
	 * @brief Construct a new NuPairMoveSetLibrary object. 
	 * 
	 * @param withBinary Bool, if true, read from binary; if false, read from txt parameter tables.
	 * @param binaryMode Int, if 1, read dumped object; if 2, read BinaryTable
	 */
	NuPairMoveSetLibrary(bool withBinary=true, int binaryMode=2);
	int dump();
	int load();
	virtual ~NuPairMoveSetLibrary();
};


class MixedNuPairCluster {
public:

	int sep; //1: nb   -1: revNb  2: nnb
	int pairType; //0~15: AA, AU, AG, AC, UA, ..., CC
	vector<int> clusterIDList;
	vector<double> clusterPList;
	NuPairMoveSetLibrary* moveLib;

	bool fixedNativeCM;
	CsMove natCM;

	int randPool[10000];

	MixedNuPairCluster(int sep, int pairType, NuPairMoveSetLibrary* lib);
	void updateEdgeInformation(EdgeInformation* ei);
	void fixNativeMove(CsMove& cm){
		this->natCM = cm;
		this->fixedNativeCM = true;
	}
	CsMove getRandomMove();

	void printMoveSetInfo();

	virtual ~MixedNuPairCluster();
};


} /* namespace NSPpredNA */

#endif /* PREDNA_NUMOVESET_H_ */
