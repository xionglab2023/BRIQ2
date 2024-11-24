/*
 * BRFoldingTree.h
 *
 */

#ifndef predNA_BRFOLDINGTREE_H_
#define predNA_BRFOLDINGTREE_H_

#include "model/RNABaseName.h"
#include "model/RiboseRotamerLib.h"
#include "forcefield/PO3Builder.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include <iostream>
#include <set>

#include "forcefield/RnaEnergyTable.h"
#include "predNA/BaseMoveLibrary.h"
#include "predNA/BRConnection.h"
#include "predNA/BRNode.h"
#include "predNA/FragmentLibrary.h"
#include "predNA/MoveMutator.h"

namespace NSPpredNA {

using namespace std;
using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPgeometry;

class BRTreeInfo {
public:
	int seqLen;
	int* seq;
	bool* connectToDownstream;
	BRNode** nodes;
	double ene;
	vector<int> freeNodeIDs;


	BRTreeInfo(int seqLen, int* seq, bool* con, BRNode** nodes, double ene, RotamerLib* rotLib);

	BRTreeInfo(int seqLen, int* seq, bool* con, BRNode** nodes, vector<BRNode*>& flexibleNodes, double ene, RotamerLib* rotLib);
	//BRTreeInfo(const string& pdbFile);

	double rmsd(BRTreeInfo* other);
	void printPDB(ofstream& of, int modelID);
	void printPDB(const string& outputFile);
	void printTmpPDB(const string& outputFile);
	virtual ~BRTreeInfo();
};


class BRFoldingTree {
public:
	int seqLen;
	int* seq;
	int* wcPairPosID;
	int* nwcPairPosID;

	string group;
	bool* fixed;
	bool* connectToDownstream;
	bool* loopRevPoints;
	bool* chainBreakPoints;

	bool* nodeConnectMatrix;

	vector<int> freeNodeIDs; //nodes for rmsd calculation, single base move and riboReverseMove could be applied to these nodes

	vector<int> bulgeList;
	vector<int> fixedList;
	vector<int> agList;

	BRNode* pseudoNode;

	RiboseRotamer** initRotList;

	vector<XYZ> constraintCoordList;
	vector<BaseDistanceMatrix> constraintDMList;

	BRNode** nodes;
	BRNode* rootNode;

	vector<BRConnection*> fixedConnectionList; //move fixed or child fixed
	vector<BRConnection*> flexibleConnectionList; //could be changed during sampling
	vector<BRNode*> flexibleNodes; //base position could be changed during sampling
	vector<BRNode*> riboFlexibleNodes; //ribose rotamer could be changed during sampling

	BRTreeInfo* initTreeInfo;
	FragmentLibrary* fragLib;
	RotamerLib* rotLib;

	int* sepTable;

	double* allBaseClashE;
	double* tmpBaseClashE;
	double* allBaseBaseE;
	double* tmpBaseBaseE;
	double* allBaseRiboseE;
	double* tmpBaseRiboseE;
	double* allBasePhoE;
	double* tmpBasePhoE;
	double* allRiboseRiboseE;
	double* tmpRiboseRiboseE;
	double* allRibosePhoE;
	double* tmpRibosePhoE;
	double* allPhoPhoE;
	double* tmpPhoPhoE;
	double* allRotE;
	double* tmpRotE;
	double* allRcE;
	double* tmpRcE;
	double* allConstraint;
	double* tmpConstraint;
	double* allPairConstraint;
	double* tmpPairConstraint;

	double* baseConstraintFactor;
	double* basePairConstraintFactor;

	RnaEnergyTable* et;
	BRFoldingTree(const string& inputFile); //only build connections
	BRFoldingTree(const string& inputFile, RnaEnergyTable* et);
	void randInit();
	void initFromKey(const string& keyInfo);
	void printBaseConnection(){
		for(int i=0;i<seqLen;i++){
			for(int j=0;j<seqLen;j++){
				if(nodeConnectMatrix[i*seqLen+j])
					cout << "1 ";
				else
					cout << "0 ";
			}
			cout << endl;
		}
	}

	void buildFromInput(vector<string>& ctList);
	void buildFixedNodes();
	void buildGroups();
	void buildAGPairs();
	void buildBasePairOfFixedNodes();
	void buildPairs();
	void buildNeighbor();
	void buildBulged13();
	void buildFrom2(BRNode* node);

	void findMidPoint(){
		for(int i=0;i<seqLen;i++){
			loopRevPoints[i] = false;
			chainBreakPoints[i] = false;
		}
		int pairID[this->seqLen];
		for(int i=0;i<seqLen;i++){
			pairID[i] = this->wcPairPosID[i];
		}

		for(int i=0;i<seqLen;i++){
			if(this->nwcPairPosID[i] >=0) {
				if(pairID[i] >= 0 || pairID[nwcPairPosID[i]] >=0)
					continue;
				pairID[i] = nwcPairPosID[i];
				pairID[nwcPairPosID[i]] = i;
			}
		}


		int loopStart = -1;
		int loopEnd = -1;
		int loopRegion[seqLen];
		for(int i=0;i<seqLen;i++){
			if(fixed[i])
				loopRegion[i] = 2; //fixed region
			else if(group[i] != '-')
				loopRegion[i] = 2; //group region
			else if(wcPairPosID[i] == -1 && nwcPairPosID[i] == -1)
				loopRegion[i] = 0; //loop region
			else if((wcPairPosID[i] > -1 && wcPairPosID[i] > i) || (nwcPairPosID[i] > -1 && nwcPairPosID[i] > i))
				loopRegion[i] = 1; //right bracket
			else
				loopRegion[i] = -1; // left bracket
		}

		for(int i=0;i<seqLen-1;i++){
			if(loopRegion[i] !=0  && loopRegion[i+1] == 0){
				loopStart = i+1;
				loopEnd = i+1;
			}

			if(loopStart > -1 && loopRegion[i] == 0 && loopRegion[i+1] != 0){
				loopEnd = i;
				if((loopEnd-loopStart) >2 && nodes[i+1]->father != NULL){
					int breakPoint = loopStart + rand()%(loopEnd-loopStart+1);
					chainBreakPoints[breakPoint] = true;
					cout << "chain break: " << breakPoint << endl;
					for(int j=breakPoint;j<=loopEnd;j++){
						loopRevPoints[j] = true;
						cout << "loopRev" << j << endl;
					}
				}

				loopStart = -1;
				loopEnd = -1;
			}
		}
	}

	void buildFrom3(BRNode* node);
	void buildReverseNodes();

	void keyInfoToRMS(const string& keyFile, const string& outFile);
	void updateCoordinate(BRNode* node);
	void updatePhoGroups();
	void updateEnergies();

	string toCtKey();
	string toCtDetailString();

	pair<double,double> ctMutEnergy(BRConnection* selectConnect, double breakCTWT, double connectWT, double clashWT,  bool verbose);
	void debugCtMutEnergy(BRConnection* selectConnect, double breakCTWT, double connectWT, double clashWT,  bool verbose);

	vector<double> getF2MutEnergyDetail(BRConnection* selectConnect, double breakCTWT,double connectWT, double clashWT,  bool verbose);

	pair<double,double> f2MutEnergy(BRConnection* selectConnect, double breakCTWT,double connectWT, double clashWT, bool verbose);
	void debugF2MutEnergy(BRConnection* selectConnect, double breakCTWT,double connectWT, double clashWT,  bool verbose);

	pair<double,double> f3MutEnergy(BRConnection* ct1,  double breakCTWT,double connectWT, double clashWT, bool verbose);

	pair<double,double> singleBaseMutEnergy(BRNode* node, double breakCTWT,double connectWT, double clashWT, bool verbose);
	void debugSingleBaseMutEnergy(BRNode* node, double breakCTWT,double connectWT, double clashWT, bool verbose);

	pair<double,double> riboseRotamerMutEnergy(BRNode* node, double breakCTWT,double connectWT, double clashWT, bool verbose);
	void debugRiboseRotamerMutEnergy(BRNode* node, double breakCTWT,double connectWT, double clashWT, bool verbose);

	double totalEnergy(double breakCTWT, double connectWT,double clashWT, bool verbose);
	double totalTmpEnergy(double breakCTWT, double connectWT,double clashWT, bool verbose);

	void printDetailEnergy();
	void printDetailEnergy(ofstream& of);


	void updateCtChildTmpCs(BRConnection* ct, CsMove& cm, bool verbose);
	void clearCtChildTmpCs(BRConnection* ct, bool verbose);
	void acceptCtChildTmpCs(BRConnection* ct, bool verbose);


	void updateF2ChildCs(BRConnection* ct, F2Fragment* frag, bool verbose);
	void updateF2ChildTmpCs(BRConnection* ct, F2Fragment* frag, bool verbose);
	void clearF2ChildTmpCs(BRConnection* ct, bool verbose);
	void acceptF2ChildTmpCs(BRConnection* ct, bool verbose);

	/*
	void updateF3ChildTmpCs(BRConnection* ct1, F3Fragment* frag, bool verbose);
	void clearF3ChildTmpCs(BRConnection* ct1, bool verbose);
	void acceptF3ChildTmpCs(BRConnection* ct1,  bool verbose);
	*/

	void updateRiboseRotamerTmp(BRNode* node, RiboseRotamer* rot, bool verbose);
	void clearRiboseRotamerTmp(BRNode* node, bool verbose);
	void acceptRiboseRotamerTmp(BRNode* node, bool verbose);

	void updateReverseRiboseRotamerTmp(BRNode* node, RiboseRotamer* rot, bool verbose);
	void clearReverseRiboseRotamerTmp(BRNode* node, bool verbose);
	void acceptReverseRiboseRotamerTmp(BRNode* node, bool verbose);

	void updateSingleBaseCoordTmp(BRNode* node, CsMove& move, bool verbose);
	void clearSingleBaseCoordTmp(BRNode* node, bool verbose);
	void acceptSingleBaseCoordTmp(BRNode* node, bool verbose);

	void trackCoordinateChangeCt(BRConnection* selectConnect);
	void trackCoordinateChangeF2(BRConnection* selectConnect);
	void trackCoordinateChangeSingleBase(BRNode* node);
	void trackCoordinateChangeRotamer(BRNode* node);

	void checkConnection();
	void checkNode();
	void checkRibose();
	void checkPho();

	void checkTmpNode();
	void energyChange();
	void checkTotalEnergy(double shift);
	void checkTmpTotalEnergy(double shift);

	BRTreeInfo* getTreeInfo();

	void printTree(int index);
	int printConnections();

	virtual ~BRFoldingTree();
};

inline double getBaseBaseEnergy(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){

	if(nodeA->seqID > nodeB->seqID)
		return getBaseBaseEnergy(nodeB, nodeA, abs(sep), et, verbose);

	double bpEnergy = 0.0;

	double minDD, dd;
	int i,j, nA, nB;
	if(squareDistance(nodeA->baseConf->coords[0], nodeB->baseConf->coords[0]) < 225.0) {
		minDD = 999999.9;
		nA = nodeA->baseConf->rot->atomNum;
		nB = nodeB->baseConf->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(nodeA->baseConf->coords[i], nodeB->baseConf->coords[j]);
				if(dd < minDD){
					minDD = dd;
				}
			}
		}

		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->baseConf->cs1, nodeB->baseConf->cs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}
	if(verbose){
		printf("minDD: %8.3f\n", minDD);
		printf("base %d base %d bpEnergy: %7.3f\n",nodeA->seqID, nodeB->seqID, bpEnergy);
	}
	return bpEnergy;
}

inline double baseBaseClash(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){

	if(nodeA->seqID > nodeB->seqID)
		return baseBaseClash(nodeB, nodeA, abs(sep), et, verbose);
	double clashEnergy = 0.0;
	double dd;
	int i,j, nA, nB;
	if(sep == 0) return 0;
	if(squareDistance(nodeA->baseConf->coords[0], nodeB->baseConf->coords[0]) < 256.0) {
		nA = nodeA->baseConf->rot->atomNum;
		nB = nodeB->baseConf->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(nodeA->baseConf->coords[i], nodeB->baseConf->coords[j]);
				if(dd < 16.0) {
					clashEnergy += et->acET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, dd, sep);
				}
			}
		}
	}
	return clashEnergy;
}

inline double getBaseBaseEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){

	if(nodeA->seqID > nodeB->seqID)
		return getBaseBaseEnergyTmp(nodeB, nodeA, abs(sep), et, verbose);
	double bpEnergy = 0.0;
	double minDD, dd;
	int i,j, nA, nB;
	if(squareDistance(nodeA->baseConfTmp->coords[0], nodeB->baseConfTmp->coords[0]) < 225.0) {
		minDD = 999999.9;
		nA = nodeA->baseConfTmp->rot->atomNum;
		nB = nodeB->baseConfTmp->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(nodeA->baseConfTmp->coords[i], nodeB->baseConfTmp->coords[j]);
				if(dd < minDD){
					minDD = dd;
				}
			}
		}
		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->baseConfTmp->cs1, nodeB->baseConfTmp->cs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}
	if(verbose){
		printf("minDD: %8.3f\n", minDD);
		printf("base %d base %d bpTmpEnergy: %7.3f\n",nodeA->seqID, nodeB->seqID, bpEnergy);
	}
	return bpEnergy;
}

inline double baseBaseClashTmp(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){

	if(nodeA->seqID > nodeB->seqID)
		return baseBaseClashTmp(nodeB, nodeA, abs(sep), et, verbose);
	double clashEnergy = 0.0;
	double dd;
	int i,j, nA, nB;
	if(sep == 0) return 0;
	if(squareDistance(nodeA->baseConfTmp->coords[0], nodeB->baseConfTmp->coords[0]) < 256.0) {
		nA = nodeA->baseConfTmp->rot->atomNum;
		nB = nodeB->baseConfTmp->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(nodeA->baseConfTmp->coords[i], nodeB->baseConfTmp->coords[j]);
				if(dd < 16.0) {
					clashEnergy += et->acET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, dd, sep);
				}
			}
		}
	}
	return clashEnergy;
}

inline double getBaseRiboseEnergy(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	double baseOxyEnergy = 0.0;
	double hbondEnergy = 0.0;
	double clashEnergy = 0.0;

	int i,j;
	double dd;
	LocalFrame csA, csB;

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(nodeA->baseConf->coords[0], nodeB->riboseConf->coords[2]) < 144.0){
		if(abs(sep) > 0) {

			baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 0, global2local(nodeA->baseConf->cs1, nodeB->riboseConf->coords[5]), sep); //O3'
			baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 1, global2local(nodeA->baseConf->cs1, nodeB->riboseConf->coords[4]), sep); //O4'

			for(i=0;i<nodeA->baseConf->rot->polarAtomNum;i++){
				hbondEnergy += et->hbET->getEnergy(nodeA->baseConf->rot->polarAtomUniqueID[i], nodeA->baseConf->csPolar[i], O2UniqueID, nodeB->riboseConf->o2Polar);
			}

			for(i=0;i<nodeA->baseConf->rot->atomNum;i++){
				for(j=0;j<nodeB->riboseConf->rot->atomNum;j++){
					dd = squareDistance(nodeA->baseConf->coords[i], nodeB->riboseConf->coords[j]);
					if(dd < 16)
						clashEnergy += et->acET->getBaseRiboseEnergy(nodeA->baseType, i, j, dd, sep);
				}
			}
		}
	}
	if(verbose){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "hbond energy: " << hbondEnergy << endl;
		cout << "clash: " << clashEnergy << endl;

	}
	return baseOxyEnergy + hbondEnergy + clashEnergy;
}


inline double getBaseRiboseEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	double hbondEnergy = 0.0;
	int i,j;
	double dd;
	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(nodeA->baseConfTmp->coords[0], nodeB->riboseConfTmp->coords[2]) < 144.0){
		if(abs(sep) > 0) {
			baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 0, global2local(nodeA->baseConfTmp->cs1, nodeB->riboseConfTmp->coords[5]), sep); //O3'
			baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 1, global2local(nodeA->baseConfTmp->cs1, nodeB->riboseConfTmp->coords[4]), sep); //O4'

			//hbondEnergy: NodeA-PolarAtom -> NodeB-O2'
			for(i=0;i<nodeA->baseConfTmp->rot->polarAtomNum;i++){
				hbondEnergy += et->hbET->getEnergy(nodeA->baseConfTmp->rot->polarAtomUniqueID[i], nodeA->baseConfTmp->csPolar[i], O2UniqueID, nodeB->riboseConfTmp->o2Polar);
			}

			for(i=0;i<nodeA->baseConfTmp->rot->atomNum;i++){
				for(j=0;j<nodeB->riboseConfTmp->rot->atomNum;j++){
					dd = squareDistance(nodeA->baseConfTmp->coords[i], nodeB->riboseConfTmp->coords[j]);
					if(dd < 16)
						clashEnergy += et->acET->getBaseRiboseEnergy(nodeA->baseType, i, j, dd, sep);
				}
			}
		}
	}
	if(verbose){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "hbond energy: " << hbondEnergy << endl;
		cout << "clash: " << clashEnergy << endl;
	}
	return baseOxyEnergy + hbondEnergy + clashEnergy;
}

inline double getBasePhoEnergy(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	if(!nodeB->connectToNeighbor) return 0;

	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	double hbondEnergy = 0.0;
	int i,j, nA;
	double dd;

	int uniqueIDOP1 = 168;
	int uniqueIDOP2 = 169;

	if(sep == 0) return 0.0;


	//hbondEnergy: NodeA-PolarAtom -> NodeB-O2'
	for(i=0;i<nodeA->baseConf->rot->polarAtomNum;i++){
		hbondEnergy += et->hbET->getEnergy(nodeA->baseConf->rot->polarAtomUniqueID[i], nodeA->baseConf->csPolar[i], uniqueIDOP1, nodeB->phoConf->op1Polar);
		hbondEnergy += et->hbET->getEnergy(nodeA->baseConf->rot->polarAtomUniqueID[i], nodeA->baseConf->csPolar[i], uniqueIDOP2, nodeB->phoConf->op2Polar);
	}


	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 2, global2local(nodeA->baseConf->cs1, nodeB->phoConf->coords[1]), sep); //O5'
	nA = nodeA->baseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->baseConf->coords[i], nodeB->phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getBasePhoEnergy(nodeA->baseType, i, j, dd, sep);
		}
	}

	if(verbose && (baseOxyEnergy+clashEnergy)>100){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "clash: " << clashEnergy << endl;
		cout << "hbond: " << hbondEnergy << endl;
	}
	return baseOxyEnergy + hbondEnergy + clashEnergy;
}

inline double getBasePhoEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	if(!nodeB->connectToNeighbor) return 0;
	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	double hbondEnergy = 0.0;
	int i,j, nA;
	double dd;

	int uniqueIDOP1 = 168;
	int uniqueIDOP2 = 169;

	if(sep == 0) return 0.0;

	//hbondEnergy: NodeA-PolarAtom -> NodeB-O2'
	for(i=0;i<nodeA->baseConfTmp->rot->polarAtomNum;i++){
		hbondEnergy += et->hbET->getEnergy(nodeA->baseConfTmp->rot->polarAtomUniqueID[i], nodeA->baseConfTmp->csPolar[i], uniqueIDOP1, nodeB->phoConfTmp->op1Polar);
		hbondEnergy += et->hbET->getEnergy(nodeA->baseConfTmp->rot->polarAtomUniqueID[i], nodeA->baseConfTmp->csPolar[i], uniqueIDOP2, nodeB->phoConfTmp->op2Polar);
	}


	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 2, global2local(nodeA->baseConfTmp->cs1, nodeB->phoConfTmp->coords[1]), sep); //O5'
	//need to be modified
	//hbondEnergy =

	nA = nodeA->baseConfTmp->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->baseConfTmp->coords[i], nodeB->phoConfTmp->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getBasePhoEnergy(nodeA->baseType, i, j, dd, sep);
		}
	}


	if(verbose && (baseOxyEnergy+clashEnergy)>10){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "clash: " << clashEnergy << endl;
		cout << "hbond: " << hbondEnergy << endl;
	}
	return baseOxyEnergy + hbondEnergy + clashEnergy;
}

inline double getRiboseRiboseEnergy(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	if(nodeA->seqID > nodeB->seqID){
		return getRiboseRiboseEnergy(nodeB, nodeA, -sep, et, verbose);
	}
	int i,j, nA, nB;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(nodeA->riboseConf->coords[2], nodeB->riboseConf->coords[2]) < 81.0) {

		if(abs(sep) >= 2)
		{
			hbondEnergy = et->hbET->getEnergy(O2UniqueID, nodeA->riboseConf->o2Polar, O2UniqueID, nodeB->riboseConf->o2Polar);
			if(verbose){
				if(abs(hbondEnergy) > 0.1){
					printf("hbond energy: %7.3f\n", hbondEnergy);
				}
			}

		}

		double dO4O2C2AB = nodeA->riboseConf->coords[4].distance((nodeB->riboseConf->coords[7] + nodeB->riboseConf->coords[1])*0.5);
		double dO4O2C2BA = nodeB->riboseConf->coords[4].distance((nodeA->riboseConf->coords[7] + nodeA->riboseConf->coords[1])*0.5);

		if(dO4O2C2AB < 5.0)
			hbondEnergy += et->hbET->getO4O2C2Energy(dO4O2C2AB, sep);
		if(dO4O2C2BA < 5.0)
			hbondEnergy += et->hbET->getO4O2C2Energy(dO4O2C2BA, sep);

		nA = nodeA->riboseConf->rot->atomNum;
		nB = nodeB->riboseConf->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(nodeA->riboseConf->coords[i], nodeB->riboseConf->coords[j]);
				if(dd < 16)
					clashEnergy += et->acET->getRiboseRiboseEnergy(i,j,dd, sep);
				if(verbose){
					if(dd < 36){
						printf("ribose ribose energy: ");
						printf("distance: %7.3f energy: %7.3f\n", sqrt(dd), et->acET->getRiboseRiboseEnergy(i,j,dd, sep));
					}
				}
			}
		}

	}
	return clashEnergy + hbondEnergy;
}

inline double getRiboseRiboseEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	if(nodeA->seqID > nodeB->seqID){
		return getRiboseRiboseEnergyTmp(nodeB, nodeA, -sep, et, verbose);
	}
	int i,j, nA, nB;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;
	int O2UniqueID = 177; //uniqueID: A-O2'


	if(squareDistance(nodeA->riboseConfTmp->coords[2], nodeB->riboseConfTmp->coords[2]) < 81.0) {

		if(abs(sep) >= 2)
		{
			hbondEnergy = et->hbET->getEnergy(O2UniqueID, nodeA->riboseConfTmp->o2Polar, O2UniqueID, nodeB->riboseConfTmp->o2Polar);
		}

		double dO4O2C2AB = nodeA->riboseConfTmp->coords[4].distance((nodeB->riboseConfTmp->coords[7] + nodeB->riboseConfTmp->coords[1])*0.5);
		double dO4O2C2BA = nodeB->riboseConfTmp->coords[4].distance((nodeA->riboseConfTmp->coords[7] + nodeA->riboseConfTmp->coords[1])*0.5);

		if(dO4O2C2AB < 5.0)
			hbondEnergy += et->hbET->getO4O2C2Energy(dO4O2C2AB, sep);
		if(dO4O2C2BA < 5.0)
			hbondEnergy += et->hbET->getO4O2C2Energy(dO4O2C2BA, sep);

		nA = nodeA->riboseConfTmp->rot->atomNum;
		nB = nodeB->riboseConfTmp->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(nodeA->riboseConfTmp->coords[i], nodeB->riboseConfTmp->coords[j]);
				if(dd < 16)
					clashEnergy += et->acET->getRiboseRiboseEnergy(i,j,dd, sep);
				if(verbose){
					//
				}
			}
		}

	}
	return clashEnergy + hbondEnergy;
}

inline double getRibosePhoEnergy(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	if(sep == 0 || sep == -1) return 0;
	if(!nodeB->connectToNeighbor) return 0;
	int i,j, nA;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	int O2UniqueID = 177; //uniqueID: A-O2'
	int uniqueIDOP1 = 168; //uniqueID: A-OP1
	int uniqueIDOP2 = 169; //uniqueID: A-OP2

	hbondEnergy += et->hbET->getEnergy(O2UniqueID, nodeA->riboseConf->o2Polar, uniqueIDOP1, nodeB->phoConf->op1Polar);
	hbondEnergy += et->hbET->getEnergy(O2UniqueID, nodeA->riboseConf->o2Polar, uniqueIDOP2, nodeB->phoConf->op2Polar);


	nA = nodeA->riboseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->riboseConf->coords[i], nodeB->phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getRibosePhoEnergy(i,j,dd, sep);
		}
	}

	return clashEnergy + hbondEnergy;
}

inline double getRibosePhoEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	if(sep == 0 || sep == -1) return 0;
	if(!nodeB->connectToNeighbor) return 0;
	int i,j, nA;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;
	int O2UniqueID = 177; //uniqueID: A-O2'
	int uniqueIDOP1 = 168; //uniqueID: A-OP1
	int uniqueIDOP2 = 169; //uniqueID: A-OP2

	hbondEnergy += et->hbET->getEnergy(O2UniqueID, nodeA->riboseConfTmp->o2Polar, uniqueIDOP1, nodeB->phoConfTmp->op1Polar);
	hbondEnergy += et->hbET->getEnergy(O2UniqueID, nodeA->riboseConfTmp->o2Polar, uniqueIDOP2, nodeB->phoConfTmp->op2Polar);

	nA = nodeA->riboseConfTmp->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->riboseConfTmp->coords[i], nodeB->phoConfTmp->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getRibosePhoEnergy(i,j,dd, sep);
		}
	}
	return clashEnergy + hbondEnergy;
}

inline double getPhoPhoEnergy(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){

	if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;

	if(nodeA->seqID > nodeB->seqID){
		return getPhoPhoEnergy(nodeB, nodeA, -sep, et, verbose);
	}
	if(squareDistance(nodeA->phoConf->coords[0], nodeB->phoConf->coords[0]) > 36.0) return 0.0;

	int i,j;
	double clashEnergy = 0;

	double dd;
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->phoConf->coords[i], nodeB->phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getPhoPhoEnergy(i,j,dd, sep);
		}
	}

	return clashEnergy;
}

inline double getPhoPhoEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTable* et, bool verbose){
	if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;

	if(nodeA->seqID > nodeB->seqID){
		return getPhoPhoEnergyTmp(nodeB, nodeA, -sep, et, verbose);
	}
	if(squareDistance(nodeA->phoConfTmp->coords[0], nodeB->phoConfTmp->coords[0]) > 36.0) return 0.0;

	int i,j;
	double clashEnergy = 0;

	double dd;
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->phoConfTmp->coords[i], nodeB->phoConfTmp->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getPhoPhoEnergy(i,j,dd, sep);
		}
	}

	return clashEnergy;
}

}


#endif /* predNA_BRFOLDINGTREE_H_ */
