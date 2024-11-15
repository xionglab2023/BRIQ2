/*
 * NuEnergyCalculator.h
 *
 *  Created on: 2023��11��23��
 *      Author: nuc
 */

#ifndef PREDNA_NUENERGYCALCULATOR_H_
#define PREDNA_NUENERGYCALCULATOR_H_
#include <vector>
#include <string>
#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "model/StructureModel.h"
#include "model/AtomLib.h"
#include "model/RotamerLib.h"
#include "model/BaseRotamer.h"
#include "model/RiboseRotamer.h"
#include "model/PhosphateRotamer.h"
#include "model/BasePairLib.h"
#include "tools/InputParser.h"
#include "tools/StringTool.h"
#include "forcefield/RnaEnergyTable.h"

namespace NSPpredNA {

using namespace NSPmodel;
using namespace std;

inline double nuBaseBaseEnergy(BaseConformer* baseConfA, BaseConformer* baseConfB, int sep, RnaEnergyTable* et, double clashRescale){

	double bpEnergy = 0.0;
	double clashEnergy = 0.0;

	double minDD, dd;
	int i,j, nA, nB;
	if(squareDistance(baseConfA->coords[0], baseConfB->coords[0]) < 225.0) {
		minDD = 999999.9;
		nA = baseConfA->rot->atomNum;
		nB = baseConfB->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(baseConfA->coords[i], baseConfB->coords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				if(dd < 16.0) {
					clashEnergy += et->acET->getBaseBaseEnergy(baseConfA->rot->baseType, i, baseConfB->rot->baseType, j, dd, sep);
				}
			}
		}

		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, sep, sqrt(minDD));

			/*
			if(sep == 1 && bpEnergy < -0.01){
				BaseDistanceMatrix dm(baseConfA->cs1, baseConfB->cs1);
				int clusterID = et->bpLib->getNeighborPairFirstFiveClusterID(dm, baseConfA->rot->baseType, baseConfB->rot->baseType);
				bpEnergy = bpEnergy * et->para->nbPairEnergyRescale[baseConfA->rot->baseType*4+baseConfB->rot->baseType][clusterID];
			}
			else if(sep == -1 && bpEnergy < -0.01){
				BaseDistanceMatrix dm(baseConfB->cs1, baseConfA->cs1);
				int clusterID = et->bpLib->getNeighborPairFirstFiveClusterID(dm, baseConfB->rot->baseType, baseConfA->rot->baseType);
				bpEnergy = bpEnergy * et->para->nbPairEnergyRescale[baseConfB->rot->baseType*4+baseConfA->rot->baseType][clusterID];
			}
			*/
		}
	}

	return bpEnergy+clashEnergy*clashRescale;
}

inline double nuBaseRiboseEnergy(BaseConformer* baseConf, RiboseConformer* riboseConf, int sep, RnaEnergyTable* et, double clashRescale){
	double baseOxyEnergy = 0.0;
	double hbondEnergy = 0.0;
	double clashEnergy = 0.0;

	int i,j;
	double dd;
	LocalFrame csA, csB;

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(baseConf->coords[0], riboseConf->coords[2]) < 144.0){
		if(abs(sep) > 0) {

			baseOxyEnergy = et->roET->getEnergy(baseConf->rot->baseType, 0, global2local(baseConf->cs1, riboseConf->coords[5]), sep); //O3'
			baseOxyEnergy += et->roET->getEnergy(baseConf->rot->baseType, 1, global2local(baseConf->cs1, riboseConf->coords[4]), sep); //O4'

			for(i=0;i<baseConf->rot->polarAtomNum;i++){
				hbondEnergy += et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], O2UniqueID, riboseConf->o2Polar);
			}

			for(i=0;i<baseConf->rot->atomNum;i++){
				for(j=0;j<riboseConf->rot->atomNum;j++){
					dd = squareDistance(baseConf->coords[i], riboseConf->coords[j]);
					if(dd < 16)
						clashEnergy += et->acET->getBaseRiboseEnergy(baseConf->rot->baseType, i, j, dd, sep);
				}
			}
		}
	}
	if(hbondEnergy > 0){
		hbondEnergy = hbondEnergy*clashRescale;
	}
	return baseOxyEnergy + hbondEnergy + clashEnergy*clashRescale;
}

inline double nuBasePhoEnergy(BaseConformer* baseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et, double clashRescale){

	if(phoConf->rot == NULL) return 0.0;

	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	double hbondEnergy = 0.0;
	int i,j, nA;
	double dd;

	int uniqueIDOP1 = 168;
	int uniqueIDOP2 = 169;

	if(sep == 0) return 0.0;

	//hbondEnergy: NodeA-PolarAtom -> NodeB-O2'
	for(i=0;i<baseConf->rot->polarAtomNum;i++){
		hbondEnergy += et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], uniqueIDOP1, phoConf->op1Polar);
		hbondEnergy += et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], uniqueIDOP2, phoConf->op2Polar);
	}

	baseOxyEnergy += et->roET->getEnergy(baseConf->rot->baseType, 2, global2local(baseConf->cs1, phoConf->coords[1]), sep); //O5'
	nA = baseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(baseConf->coords[i], phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getBasePhoEnergy(baseConf->rot->baseType, i, j, dd, sep);
		}
	}

	if(hbondEnergy > 0){
		hbondEnergy = hbondEnergy*clashRescale;
	}
	return baseOxyEnergy + hbondEnergy + clashEnergy*clashRescale;
}

inline double nuRiboseRiboseEnergy(RiboseConformer* riboseConfA, RiboseConformer* riboseConfB, int sep, RnaEnergyTable* et, double clashRescale){
	int i,j, nA, nB;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(riboseConfA->coords[2], riboseConfB->coords[2]) < 81.0) {

		if(abs(sep) >= 2)
		{
			hbondEnergy += et->hbET->getEnergy(O2UniqueID, riboseConfA->o2Polar, O2UniqueID, riboseConfB->o2Polar);
			
		}

		
		double dO4O2C2AB = riboseConfA->coords[4].distance((riboseConfB->coords[7] + riboseConfB->coords[1])*0.5);
		double dO4O2C2BA = riboseConfB->coords[4].distance((riboseConfA->coords[7] + riboseConfA->coords[1])*0.5);

		if(dO4O2C2AB < 5.0)
			hbondEnergy += et->hbET->getO4O2C2Energy(dO4O2C2AB, sep);
		if(dO4O2C2BA < 5.0)
			hbondEnergy += et->hbET->getO4O2C2Energy(dO4O2C2BA, sep);

		nA = riboseConfA->rot->atomNum;
		nB = riboseConfB->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(riboseConfA->coords[i], riboseConfB->coords[j]);
				if(dd < 16)
					clashEnergy += et->acET->getRiboseRiboseEnergy(i,j,dd, sep);
			}
		}
	}

	if(hbondEnergy > 0)
		hbondEnergy = hbondEnergy*clashRescale;
	return clashEnergy*clashRescale + hbondEnergy;
}

inline double nuRibosePhoEnergy(RiboseConformer* riboseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et, double clashRescale){

	if(phoConf->rot == NULL) return 0.0;

	if(sep == 0 || sep == -1) return 0;

	int i,j, nA;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	int O2UniqueID = 177; //uniqueID: A-O2'
	int uniqueIDOP1 = 168; //uniqueID: A-OP1
	int uniqueIDOP2 = 169; //uniqueID: A-OP2

	hbondEnergy += et->hbET->getEnergy(O2UniqueID, riboseConf->o2Polar, uniqueIDOP1, phoConf->op1Polar);
	hbondEnergy += et->hbET->getEnergy(O2UniqueID, riboseConf->o2Polar, uniqueIDOP2, phoConf->op2Polar);

	if(hbondEnergy > 0){
		hbondEnergy = hbondEnergy*clashRescale;
	}

	nA = riboseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(riboseConf->coords[i], phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getRibosePhoEnergy(i,j,dd, sep);
		}
	}
	return clashEnergy*clashRescale + hbondEnergy;
}

inline double nuPhoPhoEnergy(PhosphateConformer* phoConfA, PhosphateConformer* phoConfB, int sep, RnaEnergyTable* et, double clashRescale){

	if(phoConfA->rot == NULL) return 0.0;
	if(phoConfB->rot == NULL) return 0.0;

	if(squareDistance(phoConfA->coords[0], phoConfB->coords[0]) > 36.0) return 0.0;

	int i,j;
	double clashEnergy = 0;
	double dd;
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(phoConfA->coords[i], phoConfB->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getPhoPhoEnergy(i,j,dd, sep);
		}
	}
	return clashEnergy*clashRescale;
}

inline double nuBaseBaseEnergyCG(BaseConformerCG* baseConfA, BaseConformerCG* baseConfB, int sep, RnaEnergyTable* et, double clashRescale){
	double bpEnergy = 0.0;
	double clashEnergy = 0.0;
	double minDD, dd;
	int i,j, nA, nB;
	if(squareDistance(baseConfA->coords[0], baseConfB->coords[0]) < 225.0) {
		minDD = 999999.9;
		nA = 3;
		nB = 3;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(baseConfA->coords[i], baseConfB->coords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				if(dd < 16.0) {
					clashEnergy += et->acET->getClashEnergyCG(baseConfA->rot->uniqueIDs[i],  baseConfB->rot->uniqueIDs[j], dd);
				}
			}
		}
		if(minDD < 20.25){
			bpEnergy = et->bpcgET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, sep, sqrt(minDD));
		}
	}
	return bpEnergy+clashEnergy*clashRescale;
}

inline double nuBaseBaseEnergyCG(BaseConformerCG* baseConfA, BaseConformerCG* baseConfB, int sep, RnaEnergyTable* et, double clashRescale, EdgeClusterRegister* ecr){
	double bpEnergy = 0.0;
	double clashEnergy = 0.0;
	double minDD, dd;
	int i,j, nA, nB;
	if(squareDistance(baseConfA->coords[0], baseConfB->coords[0]) < 225.0) {
		minDD = 999999.9;
		nA = 3;
		nB = 3;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(baseConfA->coords[i], baseConfB->coords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				if(dd < 16.0) {
					clashEnergy += et->acET->getClashEnergyCG(baseConfA->rot->uniqueIDs[i],  baseConfB->rot->uniqueIDs[j], dd);
				}
			}
		}
		if(minDD < 20.25){
			bpEnergy = et->bpcgET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, sep, sqrt(minDD), ecr);
		}
	}
	return bpEnergy+clashEnergy*clashRescale;
}

inline double nuBaseRiboseEnergyCG(BaseConformerCG* baseConf, RiboseConformerCG* riboConf, int sep, RnaEnergyTable* et, double clashRescale){

	double clashEnergy = 0.0;
	int i,j;
	double dd;

	if(squareDistance(baseConf->coords[0], riboConf->coords[0]) < 144.0){
		if(abs(sep) > 0) {

			for(i=0;i<3;i++){
				for(j=0;j<3;j++){
					dd = squareDistance(baseConf->coords[i], riboConf->coords[j]);
					if(dd < 16)
						clashEnergy += et->acET->getClashEnergyCG(baseConf->rot->uniqueIDs[i], riboConf->rot->uniqueIDs[j], dd);
				}
			}
		}
	}
	return  clashEnergy*clashRescale;
}

inline double nuRiboseRiboseEnergyCG(RiboseConformerCG* riboConfA, RiboseConformerCG* riboConfB, int sep, RnaEnergyTable* et, double clashRescale){
	double clashEnergy = 0.0;
	int i,j;
	double dd;

	if(squareDistance(riboConfA->coords[0], riboConfB->coords[0]) < 144.0){
		if(abs(sep) > 0) {

			for(i=0;i<3;i++){
				for(j=0;j<3;j++){
					dd = squareDistance(riboConfA->coords[i], riboConfB->coords[j]);
					if(dd < 16)
						clashEnergy += et->acET->getClashEnergyCG(riboConfA->rot->uniqueIDs[i], riboConfB->rot->uniqueIDs[j], dd);
				}
			}
		}
	}
	return  clashEnergy*clashRescale;
}

inline double nuConnectionEnergyCG(RiboseConformerCG* riboConfA, RiboseConformerCG* riboConfB, RnaEnergyTable* et, double connectRescale){

	XYZ c1A = riboConfA->coords[0];
	XYZ o3A = riboConfA->coords[1];
	XYZ c5B = riboConfB->coords[2];
	XYZ c1B = riboConfB->coords[0];
	double d = o3A.distance(c5B);
	double ang1 = angleX(c1A, o3A, c5B);
	double ang2 = angleX(o3A, c5B, c1B);
	double dihed = dihedral(c1A, o3A, c5B, c1B);

	return et->bbcgET->getEnergy(d, ang1, ang2, dihed)*connectRescale;
}

inline double nuBaseBaseEnergyPrintDetail(BaseConformer* baseConfA, BaseConformer* baseConfB, int sep, RnaEnergyTable* et, BasePairLib* bpLib, ofstream& out){

	double totE = 0;

	double bpEnergy = 0.0;
	double clashEnergy = 0.0;

	double minDD, dd;
	int i,j, nA, nB;

	char xx[200];
	if(squareDistance(baseConfA->coords[0], baseConfB->coords[0]) < 225.0) {
		minDD = 999999.9;
		nA = baseConfA->rot->atomNum;
		nB = baseConfB->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(baseConfA->coords[i], baseConfB->coords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				if(dd < 16.0) {
					clashEnergy += et->acET->getBaseBaseEnergy(baseConfA->rot->baseType, i, baseConfB->rot->baseType, j, dd, sep);
				}
			}
		}

		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, sep, sqrt(minDD));
			if(bpEnergy < 0) {
				BaseDistanceMatrix dm(baseConfA->cs1, baseConfB->cs1);
				int clusterID = bpLib->getPairType(dm, baseConfA->rot->baseType, baseConfB->rot->baseType, sep);
				//baseTypeA, baseTypeB, sep, clusterID, ene
				sprintf(xx, "BASEPAIR %d %d %d %4d %8.3f",  baseConfA->rot->baseType, baseConfB->rot->baseType, sep, clusterID, bpEnergy);
				out << string(xx) << endl;
				totE += bpEnergy;
			}
		}
	}

	if(clashEnergy != 0){
		//sep ene
		sprintf(xx, "CLASH %d %8.3f", sep, clashEnergy);
		out << string(xx) << endl;
	}

	totE += clashEnergy;
	return totE;

}

inline double nuBaseRiboseEnergyPrintDetail(BaseConformer* baseConf, RiboseConformer* riboseConf, int sep, RnaEnergyTable* et, AtomLib* atLib, ofstream& out){
	double baseOxyEnergy = 0.0;
	double hbondEnergy = 0.0;
	double clashEnergy = 0.0;
	double totE = 0;

	int i,j;
	double dd;
	LocalFrame csA, csB;
	char xx[200];

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(baseConf->coords[0], riboseConf->coords[2]) < 144.0){
		if(abs(sep) > 0) {

			baseOxyEnergy = et->roET->getEnergy(baseConf->rot->baseType, 0, global2local(baseConf->cs1, riboseConf->coords[5]), sep); //O3'
			if(abs(baseOxyEnergy) > 0) {
				//baseType oxyType sep energy
				sprintf(xx, "BASEOXY %d %d %d %8.3f", baseConf->rot->baseType, 0, sep, baseOxyEnergy);
				out << string(xx) << endl;
				totE += baseOxyEnergy;
			}

			baseOxyEnergy = et->roET->getEnergy(baseConf->rot->baseType, 1, global2local(baseConf->cs1, riboseConf->coords[4]), sep); //O4'
			if(abs(baseOxyEnergy) > 0) {
				//baseType oxyType sep energy
				sprintf(xx, "BASEOXY %d %d %d %8.3f", baseConf->rot->baseType, 1, sep, baseOxyEnergy);
				out << string(xx) << endl;
				totE += baseOxyEnergy;
			}

			for(i=0;i<baseConf->rot->polarAtomNum;i++){
				hbondEnergy = et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], O2UniqueID, riboseConf->o2Polar);
				if(abs(hbondEnergy) > 0){
					//donorIndex acceptorIndex energy
					sprintf(xx, "HBOND %d %d %8.3f", atLib->uniqueIDToDonorID[baseConf->rot->polarAtomUniqueID[i]], atLib->uniqueIDToAcceptorID[O2UniqueID],  hbondEnergy);
					out << string(xx) << endl;	
					totE += hbondEnergy;
				}
			}

			for(i=0;i<baseConf->rot->atomNum;i++){
				for(j=0;j<riboseConf->rot->atomNum;j++){
					dd = squareDistance(baseConf->coords[i], riboseConf->coords[j]);
					if(dd < 16)
						clashEnergy += et->acET->getBaseRiboseEnergy(baseConf->rot->baseType, i, j, dd, sep);
				}
			}
		}
	}
	if(clashEnergy != 0){
		//sep ene
		sprintf(xx, "CLASH %d %8.3f", sep, clashEnergy);
		out << string(xx) << endl;
	}

	return totE + clashEnergy;

}

inline double nuBasePhoEnergyPrintDetail(BaseConformer* baseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et, AtomLib* atLib, ofstream& out){

	if(phoConf->rot == NULL) return 0.0;
	double totE = 0;
	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	double hbondEnergy = 0.0;
	int i,j, nA;
	double dd;
	char xx[200];

	int uniqueIDOP1 = 168;
	int uniqueIDOP2 = 169;

	if(sep == 0) return 0.0;

	//hbondEnergy: NodeA-PolarAtom -> NodeB-O2'
	for(i=0;i<baseConf->rot->polarAtomNum;i++){
		hbondEnergy = et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], uniqueIDOP1, phoConf->op1Polar);
		if(abs(hbondEnergy) > 0){
			//donorIndex acceptorIndex energy
			sprintf(xx, "HBOND %d %d %8.3f", atLib->uniqueIDToDonorID[baseConf->rot->polarAtomUniqueID[i]], atLib->uniqueIDToAcceptorID[uniqueIDOP1],  hbondEnergy);
			out << string(xx) << endl;	
			totE += hbondEnergy;
		}

		hbondEnergy = et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], uniqueIDOP2, phoConf->op2Polar);
		if(abs(hbondEnergy) > 0){
			//donorIndex acceptorIndex energy
			sprintf(xx, "HBOND %d %d %8.3f", atLib->uniqueIDToDonorID[baseConf->rot->polarAtomUniqueID[i]], atLib->uniqueIDToAcceptorID[uniqueIDOP2],  hbondEnergy);
			out << string(xx) << endl;	
			totE += hbondEnergy;
		}
	}

	baseOxyEnergy = et->roET->getEnergy(baseConf->rot->baseType, 2, global2local(baseConf->cs1, phoConf->coords[1]), sep); //O5'
	if(abs(baseOxyEnergy) > 0){
		sprintf(xx, "BASEOXY %d %d %d %8.3f", baseConf->rot->baseType, 2, sep, baseOxyEnergy);
		out << string(xx) << endl;
		totE += baseOxyEnergy;
	}

	nA = baseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(baseConf->coords[i], phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getBasePhoEnergy(baseConf->rot->baseType, i, j, dd, sep);
		}
	}


	if(clashEnergy != 0){
		//sep ene
		sprintf(xx, "CLASH %d %8.3f", sep, clashEnergy);
		out << string(xx) << endl;
	}

	return totE + clashEnergy;
}

inline double nuRiboseRiboseEnergyPrintDetail(RiboseConformer* riboseConfA, RiboseConformer* riboseConfB, int sep, RnaEnergyTable* et, AtomLib* atLib, ofstream& out){
	int i,j, nA, nB;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;
	double totE = 0;

	char xx[200];

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(riboseConfA->coords[2], riboseConfB->coords[2]) < 81.0) {

		if(abs(sep) >= 2)
		{
			hbondEnergy = et->hbET->getEnergy(O2UniqueID, riboseConfA->o2Polar, O2UniqueID, riboseConfB->o2Polar);
			if(abs(hbondEnergy) > 0){
				//donorIndex acceptorIndex energy
				sprintf(xx, "HBOND %d %d %8.3f", atLib->uniqueIDToDonorID[O2UniqueID], atLib->uniqueIDToAcceptorID[O2UniqueID],  hbondEnergy);
				out << string(xx) << endl;	
				totE += hbondEnergy;
			}			
		}

		double dO4O2C2AB = riboseConfA->coords[4].distance((riboseConfB->coords[7] + riboseConfB->coords[1])*0.5);
		double dO4O2C2BA = riboseConfB->coords[4].distance((riboseConfA->coords[7] + riboseConfA->coords[1])*0.5);

		if(dO4O2C2AB < 5.0)
		{
			hbondEnergy = et->hbET->getO4O2C2Energy(dO4O2C2AB, sep);
			if(abs(hbondEnergy) > 0){
				sprintf(xx, "O2O4 %2d %8.3f", sep, hbondEnergy);
				out << string(xx) << endl;	
				totE += hbondEnergy;
			}

		}	
		if(dO4O2C2BA < 5.0) {			
			hbondEnergy = et->hbET->getO4O2C2Energy(dO4O2C2BA, sep);
			if(abs(hbondEnergy) > 0){
				sprintf(xx, "O2O4 %2d %8.3f", sep, hbondEnergy);
				out << string(xx) << endl;	
				totE += hbondEnergy;
			}
		}

		nA = riboseConfA->rot->atomNum;
		nB = riboseConfB->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(riboseConfA->coords[i], riboseConfB->coords[j]);
				if(dd < 16)
					clashEnergy += et->acET->getRiboseRiboseEnergy(i,j,dd, sep);
			}
		}

		if(clashEnergy != 0){
			//sep ene
			sprintf(xx, "CLASH %d %8.3f", sep, clashEnergy);
			out << string(xx) << endl;
		}
	}

	return totE + clashEnergy;

}

inline double nuRibosePhoEnergyPrintDetail(RiboseConformer* riboseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et, AtomLib* atLib, ofstream& out){

	if(phoConf->rot == NULL) return 0.0;

	if(sep == 0 || sep == -1) return 0.0;

	int i,j, nA;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double totE = 0;
	double dd;
	
	char xx[200];

	int O2UniqueID = 177; //uniqueID: A-O2'
	int uniqueIDOP1 = 168; //uniqueID: A-OP1
	int uniqueIDOP2 = 169; //uniqueID: A-OP2

	hbondEnergy = et->hbET->getEnergy(O2UniqueID, riboseConf->o2Polar, uniqueIDOP1, phoConf->op1Polar);
	if(abs(hbondEnergy) > 0){
		//donorIndex acceptorIndex energy
		sprintf(xx, "HBOND %d %d %8.3f", atLib->uniqueIDToDonorID[O2UniqueID], atLib->uniqueIDToAcceptorID[uniqueIDOP1],  hbondEnergy);
		out << string(xx) << endl;	
		totE += hbondEnergy;
	}	

	hbondEnergy = et->hbET->getEnergy(O2UniqueID, riboseConf->o2Polar, uniqueIDOP2, phoConf->op2Polar);
	if(abs(hbondEnergy) > 0){
		//donorIndex acceptorIndex energy
		sprintf(xx, "HBOND %d %d %8.3f", atLib->uniqueIDToDonorID[O2UniqueID], atLib->uniqueIDToAcceptorID[uniqueIDOP2],  hbondEnergy);
		out << string(xx) << endl;	
		totE += hbondEnergy;
	}



	nA = riboseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(riboseConf->coords[i], phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getRibosePhoEnergy(i,j,dd, sep);
		}
	}
	
	if(clashEnergy != 0){
			//sep ene
			sprintf(xx, "CLASH %d %8.3f", sep, clashEnergy);
			out << string(xx) << endl;
	}

	return totE + clashEnergy;

}

inline double nuPhoPhoEnergyPrintDetail(PhosphateConformer* phoConfA, PhosphateConformer* phoConfB, int sep, RnaEnergyTable* et, AtomLib* atLib, ofstream& out){

	if(phoConfA->rot == NULL) return 0.0;
	if(phoConfB->rot == NULL) return 0.0;

	if(squareDistance(phoConfA->coords[0], phoConfB->coords[0]) > 36.0) return 0.0;

	int i,j;
	double clashEnergy = 0;
	double dd;
	char xx[200];

	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(phoConfA->coords[i], phoConfB->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getPhoPhoEnergy(i,j,dd, sep);
		}
	}

	if(clashEnergy != 0){
		//sep ene
		sprintf(xx, "PHOREP %d %8.3f", sep, clashEnergy);
		out << string(xx) << endl;
	}

	return clashEnergy;
}

inline double nuBaseBaseEnergyForSelection(BaseConformer* baseConfA, BaseConformer* baseConfB, int sep, RnaEnergyTable* et){

	double bpEnergy = 0.0;
	double clashEnergy = 0.0;

	double minDD, dd;
	int i,j, nA, nB;
	if(squareDistance(baseConfA->coords[0], baseConfB->coords[0]) < 225.0) {
		minDD = 999999.9;
		nA = baseConfA->rot->atomNum;
		nB = baseConfB->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(baseConfA->coords[i], baseConfB->coords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				if(dd < 16.0) {
					clashEnergy += et->acET->getBaseBaseEnergy(baseConfA->rot->baseType, i, baseConfB->rot->baseType, j, dd, sep);
				}
			}
		}

		if(minDD < 20.25){

		
			bpEnergy = et->bpET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, sep, sqrt(minDD));
			if(sep == 1)
				bpEnergy = bpEnergy*0.7;
		}
	}

	return bpEnergy+clashEnergy;
}

inline double nuBaseRiboseEnergyForSelection(BaseConformer* baseConf, RiboseConformer* riboseConf, int sep, RnaEnergyTable* et){
	double baseOxyEnergy = 0.0;
	double hbondEnergy = 0.0;
	double clashEnergy = 0.0;

	int i,j;
	double dd;
	LocalFrame csA, csB;

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(baseConf->coords[0], riboseConf->coords[2]) < 144.0){
		if(abs(sep) > 0) {

			baseOxyEnergy = et->roET->getEnergy(baseConf->rot->baseType, 0, global2local(baseConf->cs1, riboseConf->coords[5]), sep); //O3'
			baseOxyEnergy += et->roET->getEnergy(baseConf->rot->baseType, 1, global2local(baseConf->cs1, riboseConf->coords[4]), sep); //O4'

			for(i=0;i<baseConf->rot->polarAtomNum;i++){
				hbondEnergy += et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], O2UniqueID, riboseConf->o2Polar);
			}

			for(i=0;i<baseConf->rot->atomNum;i++){
				for(j=0;j<riboseConf->rot->atomNum;j++){
					dd = squareDistance(baseConf->coords[i], riboseConf->coords[j]);
					if(dd < 16)
						clashEnergy += et->acET->getBaseRiboseEnergy(baseConf->rot->baseType, i, j, dd, sep);
				}
			}
		}
	}
	
	return baseOxyEnergy + hbondEnergy + clashEnergy;
}

inline double nuBasePhoEnergyForSelection(BaseConformer* baseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et){

	if(phoConf->rot == NULL) return 0.0;

	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	double hbondEnergy = 0.0;
	int i,j, nA;
	double dd;

	int uniqueIDOP1 = 168;
	int uniqueIDOP2 = 169;

	if(sep == 0) return 0.0;

	//hbondEnergy: NodeA-PolarAtom -> NodeB-O2'
	for(i=0;i<baseConf->rot->polarAtomNum;i++){
		hbondEnergy += et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], uniqueIDOP1, phoConf->op1Polar);
		hbondEnergy += et->hbET->getEnergy(baseConf->rot->polarAtomUniqueID[i], baseConf->csPolar[i], uniqueIDOP2, phoConf->op2Polar);
	}

	baseOxyEnergy += et->roET->getEnergy(baseConf->rot->baseType, 2, global2local(baseConf->cs1, phoConf->coords[1]), sep); //O5'
	nA = baseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(baseConf->coords[i], phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getBasePhoEnergy(baseConf->rot->baseType, i, j, dd, sep);
		}
	}

	
	return baseOxyEnergy + hbondEnergy + clashEnergy;
}

inline double nuRiboseRiboseEnergyForSelection(RiboseConformer* riboseConfA, RiboseConformer* riboseConfB, int sep, RnaEnergyTable* et){
	int i,j, nA, nB;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(riboseConfA->coords[2], riboseConfB->coords[2]) < 81.0) {

		if(abs(sep) >= 2)
		{
			hbondEnergy += et->hbET->getEnergy(O2UniqueID, riboseConfA->o2Polar, O2UniqueID, riboseConfB->o2Polar);
			
		}

		nA = riboseConfA->rot->atomNum;
		nB = riboseConfB->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(riboseConfA->coords[i], riboseConfB->coords[j]);
				if(dd < 16)
					clashEnergy += et->acET->getRiboseRiboseEnergy(i,j,dd, sep);
			}
		}
	}

	return clashEnergy + hbondEnergy;
}

inline double nuRibosePhoEnergyForSelection(RiboseConformer* riboseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et){

	if(phoConf->rot == NULL) return 0.0;

	if(sep == 0 || sep == -1) return 0;

	int i,j, nA;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	int O2UniqueID = 177; //uniqueID: A-O2'
	int uniqueIDOP1 = 168; //uniqueID: A-OP1
	int uniqueIDOP2 = 169; //uniqueID: A-OP2

	hbondEnergy += et->hbET->getEnergy(O2UniqueID, riboseConf->o2Polar, uniqueIDOP1, phoConf->op1Polar);
	hbondEnergy += et->hbET->getEnergy(O2UniqueID, riboseConf->o2Polar, uniqueIDOP2, phoConf->op2Polar);

	nA = riboseConf->rot->atomNum;
	for(i=0;i<nA;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(riboseConf->coords[i], phoConf->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getRibosePhoEnergy(i,j,dd, sep);
		}
	}
	return clashEnergy + hbondEnergy;
}

inline double nuPhoPhoEnergyForSelection(PhosphateConformer* phoConfA, PhosphateConformer* phoConfB, int sep, RnaEnergyTable* et){

	if(phoConfA->rot == NULL) return 0.0;
	if(phoConfB->rot == NULL) return 0.0;

	if(squareDistance(phoConfA->coords[0], phoConfB->coords[0]) > 36.0) return 0.0;

	int i,j;
	double clashEnergy = 0;
	double dd;
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(phoConfA->coords[i], phoConfB->coords[j]);
			if(dd < 16)
				clashEnergy += et->acET->getPhoPhoEnergy(i,j,dd, sep);
		}
	}
	return clashEnergy;
}


inline double baseLigandEnergy(){
	return 0.0;
}

inline double riboseLigandEnergy(){
	return 0.0;
}

inline double phoLigandEnergy(){
	return 0.0;
}

} /* namespace NSPpredNA */

#endif /* PREDNA_NUENERGYCALCULATOR_H_ */
