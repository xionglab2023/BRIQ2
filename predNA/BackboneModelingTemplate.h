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
#include "forcefield/XPara.h"
#include "model/RNABaseName.h"
#include "model/RotamerLib.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include "geometry/Angles.h"
#include "predNA/BaseMoveLibrary.h"
#include "predNA/BRConnection.h"
#include "predNA/BRFoldingTree.h"
#include "predNA/BRNode.h"

namespace NSPpredna {

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
	XPara* para;

	vector<vector<double>> initDihedsList;
	vector<vector<double>> predDihedsList;

	vector<XYZ> initBackboneAtomList;

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

	BackboneModelingTemplate(const string& inputPDB, const string& paraFile);
	BackboneModelingTemplate(const string& inputPDB);

	void updateDiheds();

	double calPhoEnergy(double len, double xang3, double xang4, int dihed1, int dihed2, int dihed3, int dihed4, int dihed5, XYZ& p, XYZ& o5, XYZ& op1, XYZ op2, int seqID){
		double e = 0;
		double u;

		double len1 = 1.605;
		double len2 = 1.592;
		double len3 = 1.422;
		double ang1 = 120.1;
		double ang2 = 103.5;
		double ang3 = 120.7;
		double ang4 = 111.1;

		u = (len-len3)*para->kBond;
		//e += u*u;

		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang3-ang3)*para->kAng;
		//e += u*u;

		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang4-ang4)*para->kAng;
		//e += u*u;

		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += pb->eDihed1[dihed1]*para->wtDihed;
		e += pb->eDihed2[dihed2]*para->wtDihed;
		e += pb->eDihed3[dihed3]*para->wtDihed;
		e += pb->eDihed4[dihed4]*para->wtDihed;
		e += pb->eDihed5[dihed5]*para->wtDihed;

		double baseOxyEnergy = 0;
		double clashEnergy = 0;
		/*

		BRNode* node;
		for(int i=0;i<seqLen;i++){
			if(i == seqID || i == seqID + 1) continue;
			node = nodes[i];
			if(squareDistance(node->cs1.origin_, p) < 144){
				// base pho energy
				baseOxyEnergy += et->roET.getEnergy(node->baseType, 3, global2local(node->cs1, o5));
				baseOxyEnergy += et->roET.getEnergy(node->baseType, 4, global2local(node->cs1, op1));
				baseOxyEnergy += et->roET.getEnergy(node->baseType, 4, global2local(node->cs1, op2));
				for(int j=0;j<node->baseAtomNum;j++) {
					clashEnergy += et->atET.getBasePhoEnergy(node->baseType, j, 1, squareDistance(node->baseAtomCoords[j], o5));
					clashEnergy += et->atET.getBasePhoEnergy(node->baseType, j, 2, squareDistance(node->baseAtomCoords[j], op1));
					clashEnergy += et->atET.getBasePhoEnergy(node->baseType, j, 3, squareDistance(node->baseAtomCoords[j], op2));
				}
				// ribose pho energy

				for(int j=0;j<8;j++){
					clashEnergy += et->atET.getRibosePhoEnergy(j,1,squareDistance(node->riboAtomCoords[j], o5));
					clashEnergy += et->atET.getRibosePhoEnergy(j,2,squareDistance(node->riboAtomCoords[j], op1));
					clashEnergy += et->atET.getRibosePhoEnergy(j,3,squareDistance(node->riboAtomCoords[j], op2));
				}

				// pho pho energy
				for(int j=1;j<4;j++){
					clashEnergy += et->atET.getPhoPhoEnergy(j,1,squareDistance(node->pho.tList[j], o5));
					clashEnergy += et->atET.getPhoPhoEnergy(j,2,squareDistance(node->pho.tList[j], op1));
					clashEnergy += et->atET.getPhoPhoEnergy(j,3,squareDistance(node->pho.tList[j], op2));
				}
			}
		}
		*/

		return e+baseOxyEnergy+clashEnergy;
	}


	//double fastBuildPho(int seqID);

	//double phoEnergy(int seqID, bool verbose); //old method

	double phoEnergy2(int seqID, bool verbose); //updated method, improper dependent dihedral energy

	//double phoEnergy3(int seqID); //pho energy is backbone dependent

	void rebuildPho(bool verbose){
		for(int i=0;i<seqLen;i++){
			phoEnergy2(i, verbose);
			this->nodes[i]->phoConf->copyValueFrom(this->nodes[i]->phoConfTmp);
		}
	}

	void randomInit();

	void updateEnergy(){
		bool verbose = false;
		int i,j, pi, pj, sep;
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


	void printPhoEnergy();
	void printPhoTmpEnergy();
	BRTreeInfo* toTreeInfo();
	double totEnergy(bool verbose);

	void debug();

	double rms(){
		vector<XYZ> backboneCoords;
		for(int i=0;i<seqLen;i++){
			for(int j=0;j<nodes[i]->riboseConf->rot->atomNum;j++){
				backboneCoords.push_back(nodes[i]->riboseConf->coords[j]);
			}
			if(connectToDownstream[i]){
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

	void printDiheds(const string& outfile);
	void printEnergyDetail();
	void runMC();
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

		if(squareDistance(nodeA->baseConf->coords[0], nodeB->riboseConf->coords[2]) < 144.0){
			if(abs(sep) > 0) {
				baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 0, global2local(nodeA->baseConf->cs1, nodeB->riboseConf->coords[5]), sep); //O3'
				baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 1, global2local(nodeA->baseConf->cs1, nodeB->riboseConf->coords[4]), sep); //O4'

				//need to be modified
				//hbondEnergy =
				for(i=0;i<nodeA->baseConf->rot->polarAtomNum;i++){
					csA = nodeA->baseConf->csPolar[i];
					csB = nodeB->riboseConf->o2Polar;

				}

				for(i=0;i<nodeA->baseConf->rot->atomNum;i++){
					for(j=0;j<nodeB->riboseConf->rot->atomNum;j++){
						dd = squareDistance(nodeA->baseConf->coords[i], nodeB->riboseConf->coords[j]);
						if(dd < 36)
							clashEnergy += et->acET->getBaseRiboseEnergy(nodeA->baseType, i, j, dd);
					}
				}
			}
		}
		if(verbose && (baseOxyEnergy+clashEnergy)>100){
			cout << "base oxy: " << baseOxyEnergy << endl;
			cout << "clash: " << clashEnergy << endl;
		}
		return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;
	}


	inline double getBaseRiboseEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		double baseOxyEnergy = 0.0;
		double clashEnergy = 0.0;
		int i,j;
		double dd;

		if(squareDistance(nodeA->baseConfTmp->coords[0], nodeB->riboseConfTmp->coords[2]) < 144.0){
			if(abs(sep) > 0) {
				baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 0, global2local(nodeA->baseConfTmp->cs1, nodeB->baseConfTmp->coords[5]), sep); //O3'
				baseOxyEnergy = et->roET->getEnergy(nodeA->baseType, 1, global2local(nodeA->baseConfTmp->cs1, nodeB->baseConfTmp->coords[4]), sep); //O4'
				//need to be modified
				//hbondEnergy =

				for(i=0;i<nodeA->baseConfTmp->rot->atomNum;i++){
					for(j=0;j<nodeB->riboseConfTmp->rot->atomNum;j++){
						dd = squareDistance(nodeA->baseConfTmp->coords[i], nodeB->riboseConfTmp->coords[j]);
						if(dd < 36)
							clashEnergy += et->acET->getBaseRiboseEnergy(nodeA->baseType, i, j, dd);
					}
				}
			}
		}
		if(verbose && (baseOxyEnergy+clashEnergy)>100){
			cout << "base oxy: " << baseOxyEnergy << endl;
			cout << "clash: " << clashEnergy << endl;
		}
		return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;
	}

	inline double getBasePhoEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		if(!nodeB->connectToNeighbor) return 0;
		double baseOxyEnergy = 0.0;
		double clashEnergy = 0.0;
		double hbondEnergy = 0.0;
		int i,j, nA;
		double dd;

		if(sep == 0) return 0.0;

		baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 2, global2local(nodeA->baseConf->cs1, nodeB->phoConf->coords[1]), sep); //O5'
		//need to be modified
		//hbondEnergy =

		nA = nodeA->baseConf->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=1;j<4;j++){
				dd = squareDistance(nodeA->baseConf->coords[i], nodeB->phoConf->coords[j]);
				if(dd < 36)
					clashEnergy += et->acET->getBasePhoEnergy(nodeA->baseType, i, j, dd);
			}
		}

		if(verbose && (baseOxyEnergy+clashEnergy)>100){
			cout << "base oxy: " << baseOxyEnergy << endl;
			cout << "clash: " << clashEnergy << endl;
		}
		return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;
	}

	inline double getBasePhoEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		if(!nodeB->connectToNeighbor) return 0;
		double baseOxyEnergy = 0.0;
		double clashEnergy = 0.0;
		double hbondEnergy = 0.0;
		int i,j, nA;
		double dd;

		if(sep == 0) return 0.0;

		baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 2, global2local(nodeA->baseConfTmp->cs1, nodeB->baseConfTmp->coords[1]), sep); //O5'
		//need to be modified
		//hbondEnergy =

		nA = nodeA->baseConfTmp->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=1;j<4;j++){
				dd = squareDistance(nodeA->baseConfTmp->coords[i], nodeB->baseConfTmp->coords[j]);
				if(dd < 36)
					clashEnergy += et->acET->getBasePhoEnergy(nodeA->baseType, i, j, dd);
			}
		}

		if(verbose && (baseOxyEnergy+clashEnergy)>100){
			cout << "base oxy: " << baseOxyEnergy << endl;
			cout << "clash: " << clashEnergy << endl;
		}
		return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;
	}

	inline double getRiboseRiboseEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		if(nodeA->seqID > nodeB->seqID){
			return getRiboseRiboseEnergyBM(nodeB, nodeA, -sep, et, verbose);
		}
		int i,j, nA, nB;
		double clashEnergy = 0;
		double hbondEnergy = 0;
		double dd;

		if(squareDistance(nodeA->riboseConf->coords[2], nodeB->riboseConf->coords[2]) < 81.0) {

			if(abs(sep) >= 2)
			{

				//need to be modified
				//hbondEnergy =
			}

			nA = nodeA->riboseConf->rot->atomNum;
			nB = nodeB->riboseConf->rot->atomNum;
			for(i=0;i<nA;i++){
				for(j=0;j<nB;j++){
					dd = squareDistance(nodeA->riboseConf->coords[i], nodeB->riboseConf->coords[j]);
					if(dd < 36)
						clashEnergy += et->acET->getRiboseRiboseEnergy(i,j,dd, sep);
					if(verbose){
						//
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

		if(squareDistance(nodeA->riboseConfTmp->coords[2], nodeB->riboseConfTmp->coords[2]) < 81.0) {

			if(abs(sep) >= 2)
			{

				//need to be modified
				//hbondEnergy =
			}

			nA = nodeA->riboseConfTmp->rot->atomNum;
			nB = nodeB->riboseConfTmp->rot->atomNum;
			for(i=0;i<nA;i++){
				for(j=0;j<nB;j++){
					dd = squareDistance(nodeA->riboseConfTmp->coords[i], nodeB->riboseConfTmp->coords[j]);
					if(dd < 36)
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
		if(sep > 2)
		{
			//need to be modified
			//hbondEnergy =
		}

		nA = nodeA->riboseConf->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=1;j<4;j++){

				dd = squareDistance(nodeA->riboseConf->coords[i], nodeB->phoConf->coords[j]);
				if(dd < 36)
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
		if(sep > 2)
		{
			//need to be modified
			//hbondEnergy =
		}

		nA = nodeA->riboseConfTmp->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=1;j<4;j++){

				dd = squareDistance(nodeA->riboseConfTmp->coords[i], nodeB->phoConfTmp->coords[j]);
				if(dd < 36)
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
		double hbondEnergy = 0;
		if(sep > 1)
		{
			//need to be modified
			//hbondEnergy =
		}

		double dd;
		for(i=1;i<4;i++){
			for(j=1;j<4;j++){
				dd = squareDistance(nodeA->phoConf->coords[i], nodeB->phoConf->coords[j]);
				if(dd < 25)
					clashEnergy += et->acET->getPhoPhoEnergy(i,j,dd, sep);
			}
		}
		return clashEnergy+hbondEnergy;
	}

	inline double getPhoPhoEnergyTmpBM(BRNode* nodeA, BRNode* nodeB, int sep, RnaEnergyTableSimple* et, bool verbose){
		if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;

		if(nodeA->seqID > nodeB->seqID){
			return getPhoPhoEnergyTmpBM(nodeB, nodeA, -sep, et, verbose);
		}
		if(squareDistance(nodeA->phoConfTmp->coords[0], nodeB->phoConfTmp->coords[0]) > 36.0) return 0.0;

		int i,j;
		double clashEnergy = 0;
		double hbondEnergy = 0;
		if(sep > 1)
		{
			//need to be modified
			//hbondEnergy =
		}

		double dd;
		for(i=1;i<4;i++){
			for(j=1;j<4;j++){
				dd = squareDistance(nodeA->phoConfTmp->coords[i], nodeB->phoConfTmp->coords[j]);
				if(dd < 25)
					clashEnergy += et->acET->getPhoPhoEnergy(i,j,dd, sep);
			}
		}
		return clashEnergy+hbondEnergy;
	}

	virtual ~BackboneModelingTemplate();
};

} /* namespace NSPforcefield */

#endif /* predNA_BACKBONEMODELINGTEMPLATE_H_ */
