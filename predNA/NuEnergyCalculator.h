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

inline double nuBaseBaseEnergy(BaseConformer* baseConfA, BaseConformer* baseConfB, int sep, RnaEnergyTable* et){

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

			
			if(sep == 1 || sep == -1){
				/*test code, need to be modified*/
				BaseDistanceMatrix dm(baseConfA->cs1, baseConfB->cs1);
				double pCluster = et->bpLib->getPairClusterProportion(dm, baseConfA->rot->baseType, baseConfB->rot->baseType, sep);
				if(pCluster < et->para->pNbClusterCutoff){
					bpEnergy = et->bpET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, 2, sqrt(minDD));
				}
				else 
					bpEnergy = et->bpET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, sep, sqrt(minDD));
			}
			else 
				bpEnergy = et->bpET->getEnergy(baseConfA->cs1, baseConfB->cs1, baseConfA->rot->baseType, baseConfB->rot->baseType, sep, sqrt(minDD));
		}
	}

	return bpEnergy+clashEnergy;
}

inline double nuBaseRiboseEnergy(BaseConformer* baseConf, RiboseConformer* riboseConf, int sep, RnaEnergyTable* et){
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
			baseOxyEnergy = et->roET->getEnergy(baseConf->rot->baseType, 1, global2local(baseConf->cs1, riboseConf->coords[4]), sep); //O4'

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

inline double nuBasePhoEnergy(BaseConformer* baseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et){

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

inline double nuRiboseRiboseEnergy(RiboseConformer* riboseConfA, RiboseConformer* riboseConfB, int sep, RnaEnergyTable* et){
	int i,j, nA, nB;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	int O2UniqueID = 177; //uniqueID: A-O2'

	if(squareDistance(riboseConfA->coords[2], riboseConfB->coords[2]) < 81.0) {

		if(abs(sep) >= 2)
		{
			hbondEnergy = et->hbET->getEnergy(O2UniqueID, riboseConfA->o2Polar, O2UniqueID, riboseConfB->o2Polar);
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
	return clashEnergy + hbondEnergy;
}

inline double nuRibosePhoEnergy(RiboseConformer* riboseConf, PhosphateConformer* phoConf, int sep, RnaEnergyTable* et){

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

inline double nuPhoPhoEnergy(PhosphateConformer* phoConfA, PhosphateConformer* phoConfB, int sep, RnaEnergyTable* et){

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

inline double nuBaseBaseEnergyCG(BaseConformerCG* baseConfA, BaseConformerCG* baseConfB, int sep, RnaEnergyTable* et){
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
	return bpEnergy+clashEnergy;
}

inline double nuBaseRiboseEnergyCG(BaseConformerCG* baseConf, RiboseConformerCG* riboConf, int sep, RnaEnergyTable* et){

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
	return  clashEnergy;
}

inline double nuRiboseRiboseEnergyCG(RiboseConformerCG* riboConfA, RiboseConformerCG* riboConfB, int sep, RnaEnergyTable* et){
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
	return  clashEnergy;
}

inline double nuConnectionEnergyCG(RiboseConformerCG* riboConfA, RiboseConformerCG* riboConfB, RnaEnergyTable* et){

	XYZ c1A = riboConfA->coords[0];
	XYZ o3A = riboConfA->coords[1];
	XYZ c5B = riboConfB->coords[2];
	XYZ c1B = riboConfB->coords[0];
	double d = o3A.distance(c5B);
	double ang1 = angleX(c1A, o3A, c5B);
	double ang2 = angleX(o3A, c5B, c1B);
	double dihed = dihedral(c1A, o3A, c5B, c1B);

	return et->bbcgET->getEnergy(d, ang1, ang2, dihed);
}

} /* namespace NSPpredNA */

#endif /* PREDNA_NUENERGYCALCULATOR_H_ */
