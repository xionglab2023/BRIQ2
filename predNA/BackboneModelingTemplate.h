/*
 * BackboneModelingTemplate.h
 *
 */

#ifndef predNA_BACKBONEMODELINGTEMPLATE_H_
#define predNA_BACKBONEMODELINGTEMPLATE_H_

#include <iostream>
#include <set>

#include "forcefield/RnaEnergyTableSimple.h"
#include "forcefield/PO3Builder.h"
#include "forcefield/ForceFieldPara.h"
#include "model/RNABaseName.h"
#include "model/RotamerLib.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include "geometry/Angles.h"
#include "predNA/BaseMoveLibrary.h"
#include "predNA/BRConnection.h"
#include "predNA/BRFoldingTree.h"
#include "predNA/BRNode.h"

namespace NSPpredNA {

using namespace NSPgeometry;
using namespace std;
using namespace NSPforcefield;
using namespace NSPmodel;

class BackboneModelingTemplate {
public:

	int seqLen;
	int* seq;
	bool* connectToDownstream;
	BRNode** nodes;
	RnaEnergyTableSimple* et;
	PO3Builder* pb;
	RotamerLib* rotLib;
	vector<XYZ> initBackboneAtomList;
	ForceFieldPara* para;

	vector<int> mutNodeList;

	int* sepTable;
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

	BackboneModelingTemplate(const string& inputPDB, ForceFieldPara* para);

	BackboneModelingTemplate(RNAPDB* pdb, const string& cnt, const string& csn, ForceFieldPara* para);

	double buildPho(int seqID, bool verbose); //updated method, improper dependent dihedral energy

	//double phoEnergy3(int seqID); //pho energy is backbone dependent

	void rebuildPho(bool verbose){
		for(int i=0;i<seqLen;i++){
			buildPho(i, verbose);
			this->nodes[i]->phoConf->copyValueFrom(this->nodes[i]->phoConfTmp);
		}
	}

	void randomInit();

	void updateEnergy(){
		bool verbose = false;
		int i,j, pi, pj;
		for(int i=0;i<seqLen;i++){
			for(int j=0;j<seqLen;j++){
				pi = i*seqLen+j;
				allBaseRiboseE[pi] = 0;
				allBasePhoE[pi] = 0;
				allRiboseRiboseE[pi] = 0;
				allRibosePhoE[pi] = 0;
				allPhoPhoE[pi] = 0;
			}
		}

		for(i=0;i<seqLen;i++){
			BRNode* nodeA = nodes[i];
			for(j=i+1;j<seqLen;j++){
				BRNode* nodeB = nodes[j];
				pi = nodeA->seqID*seqLen+nodeB->seqID;
				pj = nodeB->seqID*seqLen+nodeA->seqID;


				allBaseRiboseE[pi] = getBaseRiboseEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allBaseRiboseE[pj] = getBaseRiboseEnergyBM(nodeB, nodeA, sepTable[pj], et, verbose);

				allBasePhoE[pi] = getBasePhoEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allBasePhoE[pj] = getBasePhoEnergyBM(nodeB, nodeA, sepTable[pj], et, verbose);

				allRiboseRiboseE[pi] = getRiboseRiboseEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allRiboseRiboseE[pj] = allRibosePhoE[pi];

				allRibosePhoE[pi] = getRibosePhoEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allRibosePhoE[pj] = getRibosePhoEnergyBM(nodeB, nodeA, sepTable[pj], et, verbose);

				allPhoPhoE[pi] = getPhoPhoEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allPhoPhoE[pj] = allPhoPhoE[pi];
			}
		}


		for(i=0;i<seqLen;i++){
			allRotE[i] = nodes[i]->riboseConf->rot->energy;
		}

		for(i=0;i<seqLen;i++){
			allRcE[i] = nodes[i]->phoConf->ene;
		}

		for(i=0;i<seqLen;i++){
			tmpRotE[i] = allRotE[i];
			tmpRcE[i] = allRcE[i];
			for(j=0;j<seqLen;j++){
				tmpBaseRiboseE[i*seqLen+j] = allBaseRiboseE[i*seqLen+j];
				tmpBasePhoE[i*seqLen+j] = allBasePhoE[i*seqLen+j];
				tmpRiboseRiboseE[i*seqLen+j] = allRiboseRiboseE[i*seqLen+j];
				tmpRibosePhoE[i*seqLen+j] = allRibosePhoE[i*seqLen+j];
				tmpPhoPhoE[i*seqLen+j] = allPhoPhoE[i*seqLen+j];
			}
		}
	}

	double mutEnergy(BRNode* node, bool verbose);
	void updateTmpRotamer(BRNode* node, RiboseRotamer* rot, bool verbose);
	void clearTmpRotamer(BRNode* node, bool verbose);
	void acceptTmpRotamer(BRNode* node, bool verbose);

	BRTreeInfo* toTreeInfo();
	double totEnergy(bool verbose);

	void debug();

	double rms(){
		vector<XYZ> backboneCoords;
		for(int i=0;i<seqLen;i++){

			if(i>0 && i < seqLen-1 && connectToDownstream[i-1] && connectToDownstream[i]) {
				for(int j=0;j<nodes[i]->riboseConf->rot->atomNum;j++){
					backboneCoords.push_back(nodes[i]->riboseConf->coords[j]);
				}
				backboneCoords.push_back(nodes[i]->phoConf->coords[0]);
				backboneCoords.push_back(nodes[i]->phoConf->coords[1]);
			}
		}

		double rms = 0;
		for(int i=0;i<initBackboneAtomList.size();i++){
			rms += squareDistance(initBackboneAtomList[i], backboneCoords[i]);
		}
		return sqrt(rms/initBackboneAtomList.size());
	}


	void printEnergyDetail();
	double runMC();
	double fragMC();

	void checkTotalEnergy();
	void checkTotalEnergyTmp();

	inline double getBaseRiboseEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
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


	inline double getBaseRiboseEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
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

	inline double getBasePhoEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
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

	inline double getBasePhoEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
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

	inline double getRiboseRiboseEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		if(nodeA->seqID > nodeB->seqID){
			return getRiboseRiboseEnergyBM(nodeB, nodeA, -sep, et, verbose);
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

	inline double getRiboseRiboseEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		if(nodeA->seqID > nodeB->seqID){
			return getRiboseRiboseEnergyTmpBM(nodeB, nodeA, -sep, et, verbose);
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

	inline double getRibosePhoEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
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

	inline double getRibosePhoEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
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


	inline double getPhoPhoEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){

		if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;

		if(nodeA->seqID > nodeB->seqID){
			return getPhoPhoEnergyBM(nodeB, nodeA, -sep, et, verbose);
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

	inline double getPhoPhoEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;

		if(nodeA->seqID > nodeB->seqID){
			return getPhoPhoEnergyTmpBM(nodeB, nodeA, -sep, et, verbose);
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


	virtual ~BackboneModelingTemplate();
};

} /* namespace NSPforcefield */

#endif /* predNA_BACKBONEMODELINGTEMPLATE_H_ */
