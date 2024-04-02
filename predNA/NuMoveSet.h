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
#include <unistd.h>
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
	int pairType; //0~15: AA, AU, AG, AC, UA, ..., CC; pairType = -1 means non-contact base pair cluster
	int clusterID; //non-contact base pair clusterID : 0~1050

	/*
	 * Ust SP1000 move index, distance 50 bins, dihedral angle 40 bins, sphere 1000*1000
	 * move set were divided into 50 groups, the total probability of each group is 2%
	 */
	vector<int> moveIndexList[50];

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
	 * @param outs: ostream object to dump to 
	 * @return int 
	 */
	int dump(ostream& outs);

	/**
	 * @brief load object from binary cache
	 * 
	 * @return int 
	 */
	int load(istream& ins);


	CsMove getRandomMove(OrientationIndex* oi);

	CsMove getRandomMoveWithFixedSubCluster(OrientationIndex* oi, int subClusterID);
	
	void printInfo(){
		cout << "sep: " << sep << endl;
		cout << "pair: " << pairType << endl;
		cout << "clusterID: " << clusterID << endl;
		for(int i=0;i<20;i++){
			cout << i << " " << moveIndexList[i].size() << endl;
		}
	}


	virtual ~IndividualNuPairMoveSet();
};

class NuPairMoveSetLibrary {
private:
	vector<IndividualNuPairMoveSet*> nbMoveList[16];
	vector<IndividualNuPairMoveSet*> revNbMoveList[16];

	int nbContactClusterNum[16];
	int revNbContactClusterNum[16];

	int nbNonContactClusterNum;

	vector<IndividualNuPairMoveSet*> nbNonContactMoveList;
	vector<IndividualNuPairMoveSet*> revNbNonContactMoveList;

	vector<IndividualNuPairMoveSet*> nnbMoveList[16];


	/**
	 * @brief Construct a new NuPairMoveSetLibrary object. 
	 * 
	 * @param withBinary Bool, if true, read from binary; if false, read from txt parameter tables.
	 * @param binaryMode Int, if 1, read dumped BinaryCache; if 2, read BinaryTable
	 */
	

public:
	OrientationIndex* oi;
	NuPairMoveSetLibrary(bool withBinary=true, int binaryMode=1);
	int dump();
	int load();
	void printMoveLibInfo();
	IndividualNuPairMoveSet* getMoveSet(int type, int clusterID, int sep);
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
	string type;

	int randPool[100000];

	MixedNuPairCluster(int sep, int pairType, NuPairMoveSetLibrary* lib);
	void updateEdgeInformation(EdgeInformation* ei);
	void fixNativeMove(CsMove& cm){
		this->natCM = cm;
		this->fixedNativeCM = true;
	}

	CsMove getRandomMove();
	CsMove getRandomMoveWithFixedCluster(CsMove& move);
	CsMove getRandomMoveWithFixedSubCluster(CsMove& move);
	CsMove getRandomMoveWithFixedSP1000Index(CsMove& move);

	void printMoveSetInfo();

	virtual ~MixedNuPairCluster();
};


} /* namespace NSPpredNA */

#endif /* PREDNA_NUMOVESET_H_ */
