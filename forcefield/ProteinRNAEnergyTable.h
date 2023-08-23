/*
 * ProteinRNAEnergyTable.h
 *
 *  Created on: 2022Äê7ÔÂ19ÈÕ
 *      Author: pengx
 */

#ifndef FORCEFIELD_PROTEINRNAENERGYTABLE_H_
#define FORCEFIELD_PROTEINRNAENERGYTABLE_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "dataio/datapaths.h"
#include "forcefield/XPara.h"
#include "tools/StringTool.h"
#include "model/AtomLib.h"
#include "model/AtomLib.h"
#include <time.h>
#include <math.h>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;
using namespace NSPforcefield;
using namespace NSPmodel;

class ProteinRNAEnergyTable {
public:

	/*
	 * polarAtoms:
	 * 0: BBN
	 * 1: BBO
	 */

	vector<double> etList[8];
	vector<bool> directDep[8];
	vector<XYZ> hbAtom[8];
	vector<double> ang1List[8];
	vector<double> ang2List[8];
	vector<double> ang3List[8];
	XPara* para;

	vector<double> baseVdwRadii[4];
	vector<double> riboseVdwRadii;
	vector<double> phoVdwRadii;

	AtomLib* atomLib;

	ProteinRNAEnergyTable();

	double clashEnergy(double d0, double d){
		double u, e, slop, e1;

		if(d < d0 - 0.4){
			u = -0.4*para->lamdaClash;
			slop = 4*u*u*u*para->lamdaClash;
			e1 = u*u*u*u;
			e = e1 + slop*(d-d0+0.4);
		}
		else if(d < d0) {
			u = (d-d0)*para->lamdaClash;
			e = u*u*u*u;
		}
		else
			e = 0;
		return e*para->wtClash;
	}

	double getBaseProClash(int baseType, int atomIndex, int proUniqueID, double d, double shift){
		return 0.0;
	}

	double getRiboseProClash(int atomIndex, int proUniqueID, double d, double shift){
		return 0.0;
	}

	double getPhoProClash(int atomIndex, int proUniqueID, double d, double shift){
		return 0.0;
	}

	double getPolarEnergy(int baseType, int polarAtomType, const XYZ& localCoreCoord, const XYZ& localSupportCoord){
		/*
		 * protein polarAtomType:
		 * 0: BBN
		 * 1: BBO
		 */

		if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
		if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
		if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
		int idX = (int)((localCoreCoord.x_+3.0)*5);
		int idY = (int)((localCoreCoord.y_+6.0)*5);
		int idZ = (int)((localCoreCoord.z_+5.0)*5);
		int type = baseType*2+polarAtomType;
		int index = idX*4000+idY*50+idZ;

		double ene = etList[type][index];
		if(ene > 0)
			return ene;
		if(!directDep[type][index])
			return ene;

		double ang = NSPgeometry::angleX(hbAtom[type][index], localCoreCoord, localSupportCoord);
		double ang1 = ang1List[type][index];
		double ang2 = ang2List[type][index];
		double ang3 = ang3List[type][index];
		if(ang < ang1 || ang > ang3)
			return 0;

		if(ang < ang2){
			double u = (ang - ang2)/(ang1 - ang2);
			return ene * (1 - u*u);
		}
		else {
			double u = (ang - ang2)/(ang2 - ang3);
			return ene * (1 - u*u);
		}
	}
	virtual ~ProteinRNAEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_PROTEINRNAENERGYTABLE_H_ */
