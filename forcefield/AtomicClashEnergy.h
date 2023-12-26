/*
 * AtomicClashEnergy.h
 *
 *  Created on: 2023��4��25��
 *      Author: nuc
 */

#ifndef FORCEFIELD_ATOMICCLASHENERGY_H_
#define FORCEFIELD_ATOMICCLASHENERGY_H_

#include <vector>
#include <map>
#include <fstream>
#include "geometry/CsMove.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "geometry/localframe.h"
#include "dataio/datapaths.h"
#include "forcefield/ForceFieldPara.h"
#include "tools/StringTool.h"
#include "model/AtomLib.h"
#include "model/BaseRotamer.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPmodel;
using namespace NSPgeometry;
using namespace NSPtools;

class AtomicClashEnergy {
public:

	ForceFieldPara* ffp;
	AtomLib* atLib;
	/*
	 * pre-calculated clash energy, d0 ranges from 2 angstrom to 4 angstrom, d*d ranges from 0 to 25 angstrom^2
	 */
	double clashEnergyTableNb[200][2500];
	double clashEnergyTableNnb[200][2500];

	/*
	 * Clash energy is only calculated between atoms separated by more than five chemical bonds.
	 */

	/*
	 * base atom uniqueID
	 * atomName atomUniqueID
	 */
	int baseAtomUniqueID[8][11];

	/*
	 * ribose atom uniqueID
	 * 0:  C1' -> 178
	 * 1:  C2' -> 176
	 * 2:  C3' -> 174
	 * 3:  C4' -> 172
	 * 4:  O4' -> 173
	 * 5:  O3' -> 175
	 * 6:  C5' -> 171
	 * 7:  O2' -> 177
	 */
	int riboseUniqueID[8];

	/*
	 * phosphate atom uniqueID
	 * P       0 -> 167
	 * O5'     1 -> 170
	 * OP1     2 -> 168
	 * OP2     3 -> 169
	 */
	int phosphateUniqueID[4];

	double atomRadius[362];
	bool isDonor[362];
	bool isAcceptor[362];

	AtomicClashEnergy(ForceFieldPara* ffp);

	double clash(double d0, double d, double lamda, double shift){
		double u, e, slop, e1;

		if(d < d0+shift - 0.4){
			u = -0.4*lamda;
			slop = 4*u*u*u*lamda;
			e1 = u*u*u*u;
			e = e1 + slop*(d-d0-shift+0.4);
		}
		else if(d < d0+shift) {
			u = (d-d0-shift)*lamda;
			e = u*u*u*u;
		}
		else
			e = 0;
		return e;
	}

	double getBaseBaseEnergy(int baseTypeA, int atomTypeA, int baseTypeB, int atomTypeB, double dd, int sep){
		if(dd >= 16.00) return 0;
		int uniqueIDA = baseAtomUniqueID[baseTypeA][atomTypeA];
		int uniqueIDB = baseAtomUniqueID[baseTypeB][atomTypeB];
		if(isDonor[uniqueIDA] && isAcceptor[uniqueIDB]) return 0.0;
		if(isDonor[uniqueIDB] && isAcceptor[uniqueIDA]) return 0.0;
		double d0 = atomRadius[uniqueIDA] + atomRadius[uniqueIDB];
		if(abs(sep) == 1)
			return this->clashEnergyTableNb[(int)(d0*100-200)][(int)(dd*100)];
		else
			return this->clashEnergyTableNnb[(int)(d0*100-200)][(int)(dd*100)];
	}

	double getBaseRiboseEnergy(int baseTypeA, int atomTypeA, int riboseAtomTypeB, double dd, int sep){
		if(dd >= 16.00) return 0;
		int uniqueIDA = baseAtomUniqueID[baseTypeA][atomTypeA];
		int uniqueIDB = riboseUniqueID[riboseAtomTypeB];
		if(isDonor[uniqueIDA] && isAcceptor[uniqueIDB]) return 0.0;
		if(isDonor[uniqueIDB] && isAcceptor[uniqueIDA]) return 0.0;
		double d0 = atomRadius[uniqueIDA] + atomRadius[uniqueIDB];
		if(abs(sep) == 1)
			return this->clashEnergyTableNb[(int)(d0*100-200)][(int)(dd*100)];
		else
			return this->clashEnergyTableNnb[(int)(d0*100-200)][(int)(dd*100)];
	}

	double getBasePhoEnergy(int baseTypeA, int atomTypeA, int phoAtomTypeB, double dd, int sep){
		if(dd >= 16.00) return 0;
		int uniqueIDA = baseAtomUniqueID[baseTypeA][atomTypeA];
		int uniqueIDB = phosphateUniqueID[phoAtomTypeB];
		if(isDonor[uniqueIDA] && isAcceptor[uniqueIDB]) return 0.0;
		if(isDonor[uniqueIDB] && isAcceptor[uniqueIDA]) return 0.0;
		double d0 = atomRadius[uniqueIDA] + atomRadius[uniqueIDB];
		if(abs(sep) == 1)
			return this->clashEnergyTableNb[(int)(d0*100-200)][(int)(dd*100)];
		else
			return this->clashEnergyTableNnb[(int)(d0*100-200)][(int)(dd*100)];
	}

	double getRiboseRiboseEnergy(int riboseTypeA, int riboseTypeB, double dd, int sep){
		/*
		 * sep = seqIdB - seqIdA
		 * if |sep| > 1, sep = 2
		 * We don't calculate the clash energy between O3' of ribose i and C5' of ribose i+1
		 */
		if(sep == 1 && riboseTypeA == 5 && riboseTypeB == 6) return 0.0;
		if(sep == -1 && riboseTypeA == 6 && riboseTypeB == 5) return 0.0;
		if(dd >= 16.00) return 0;
		int uniqueIDA = riboseUniqueID[riboseTypeA];
		int uniqueIDB = riboseUniqueID[riboseTypeB];
		if(isDonor[uniqueIDA] && isAcceptor[uniqueIDB]) return 0.0;
		if(isDonor[uniqueIDB] && isAcceptor[uniqueIDA]) return 0.0;
		double d0 = atomRadius[uniqueIDA] + atomRadius[uniqueIDB];
		if(abs(sep) == 1)
			return this->clashEnergyTableNb[(int)(d0*100-200)][(int)(dd*100)];
		else
			return this->clashEnergyTableNnb[(int)(d0*100-200)][(int)(dd*100)];
	}

	double getRibosePhoEnergy(int riboseTypeA, int phoTypeB, double dd, int sep){
		if(sep == -1) return 0.0;
		if(dd >= 16.00) return 0;
		int uniqueIDA = riboseUniqueID[riboseTypeA];
		int uniqueIDB = phosphateUniqueID[phoTypeB];
		if(isDonor[uniqueIDA] && isAcceptor[uniqueIDB]) return 0.0;
		if(isDonor[uniqueIDB] && isAcceptor[uniqueIDA]) return 0.0;
		double d0 = atomRadius[uniqueIDA] + atomRadius[uniqueIDB];
		if(abs(sep) == 1)
			return this->clashEnergyTableNb[(int)(d0*100-200)][(int)(dd*100)];
		else
			return this->clashEnergyTableNnb[(int)(d0*100-200)][(int)(dd*100)];
	}

	double getPhoPhoEnergy(int phoTypeA, int phoTypeB, double dd, int sep){
		if(dd >= 16.00) return 0.0;
		int uniqueIDA = phosphateUniqueID[phoTypeA];
		int uniqueIDB = phosphateUniqueID[phoTypeB];
		if(isDonor[uniqueIDA] && isAcceptor[uniqueIDB]) return 0.0;
		if(isDonor[uniqueIDB] && isAcceptor[uniqueIDA]) return 0.0;
		double d0 = atomRadius[uniqueIDA] + atomRadius[uniqueIDB];
		if(abs(sep) == 1)
			return this->clashEnergyTableNb[(int)(d0*100-200)][(int)(dd*100)];
		else
			return this->clashEnergyTableNnb[(int)(d0*100-200)][(int)(dd*100)];
	}

	double getClashEnergyCG(int uniqueIDA, int uniqueIDB, double dd){
		if(dd >= 16.00) return 0;
		double d0 = atomRadius[uniqueIDA] + atomRadius[uniqueIDB];
		if(isDonor[uniqueIDA] && isAcceptor[uniqueIDB]) d0 = 2.8;
		if(isDonor[uniqueIDB] && isAcceptor[uniqueIDA]) d0 = 2.8;
		return this->clashEnergyTableNnb[(int)(d0*100-200)][(int)(dd*100)];
	}	


	virtual ~AtomicClashEnergy();
};

} /* namespace NSPmodel */

#endif /* FORCEFIELD_ATOMICCLASHENERGY_H_ */
