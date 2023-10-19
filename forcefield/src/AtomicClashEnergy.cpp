/*
 * AtomicClashEnergy.cpp
 *
 *  Created on: 2023Äê4ÔÂ25ÈÕ
 *      Author: pengx
 */

#include "forcefield/AtomicClashEnergy.h"

namespace NSPforcefield {

AtomicClashEnergy::AtomicClashEnergy(ForceFieldPara* ffp) {
	// TODO Auto-generated constructor stub
	this->ffp = ffp;
	this->atLib = new AtomLib();

	double d0,d;
	double shift = 0.0;

	for(int i=0;i<200;i++){
		d0 = 2.0 + i*0.01 + 0.005;
		for(int j=0;j<2500;j++){
			d = sqrt(j*0.01+0.005);
			clashEnergyTableNb[i][j] = clash(d0,d,ffp->clashLamdaNb,shift);
			clashEnergyTableNnb[i][j] = clash(d0,d,ffp->clashLamdaNnb,shift);
		}
	}


	//uniqueID of base A and base DA
	for(int i=0;i<10;i++){
		baseAtomUniqueID[0][i] = 179+i;
		baseAtomUniqueID[4][i] = 263+i;
	}


	//uniqueID of base U and base DT
	for(int i=0;i<8;i++){
		baseAtomUniqueID[1][i] = 201+i;
	}
	for(int i=0;i<9;i++){
		baseAtomUniqueID[5][i] = 284+i;
	}

	//uniqueID of base G and base DG
	for(int i=0;i<11;i++){
		baseAtomUniqueID[2][i] = 221+i;
		baseAtomUniqueID[6][i] = 304+i;
	}

	//uniqueID of base C and base DC
	for(int i=0;i<8;i++){
		baseAtomUniqueID[3][i] = 244+i;
		baseAtomUniqueID[7][i] = 326+i;
	}

	//uniqueID of ribose
	riboseUniqueID[0] = 178;
	riboseUniqueID[1] = 176;
	riboseUniqueID[2] = 174;
	riboseUniqueID[3] = 172;
	riboseUniqueID[4] = 173;
	riboseUniqueID[5] = 175;
	riboseUniqueID[6] = 171;
	riboseUniqueID[7] = 177;

	//uniqueID of phosphate
	phosphateUniqueID[0] = 167;
	phosphateUniqueID[1] = 170;
	phosphateUniqueID[2] = 168;
	phosphateUniqueID[3] = 169;

	for(int i=0;i<362;i++){
		AtomProperty* ap = atLib->getAtomProperty(i);
		atomRadius[i] = ap->vdwRadius;
		isDonor[i] = ap->isHDonor;
		isAcceptor[i] = ap->isHAcceptor;
	}

}

AtomicClashEnergy::~AtomicClashEnergy() {
	// TODO Auto-generated destructor stub
	delete atLib;
}

} /* namespace NSPmodel */
