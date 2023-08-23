/*
 * ProAtomicEnergyTable.h
 *
 */

#ifndef FORCEFIELD_PROATOMICENERGYTABLE_H_
#define FORCEFIELD_PROATOMICENERGYTABLE_H_
#include <vector>
#include <map>
#include <fstream>
#include "geometry/CsMove.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "geometry/localframe.h"
#include "model/AtomLib.h"
#include "dataio/datapaths.h"
#include "forcefield/XPara.h"
#include "para/DesignPara.h"
#include "para/ProParameter.h"

namespace NSPforcefield {
using namespace std;
using namespace NSPpara;
using namespace NSPgeometry;
using namespace NSPmodel;

class ProAtomicEnergyTable {

public:

	double vdwCurve[1000];
	double dsEnergy[1458000];
	double distAve[54];
	double distSd[54];
	double donerAngleAve[54];
	double donerAngleSdLeft[54];
	double donerAngleSdRight[54];
	double acceptorAngleAve[54];
	double acceptorAngleSdLeft[54];
	double acceptorAngleSdRight[54];



	ProAtomicEnergyTable(DesignPara* para);
	ProAtomicEnergyTable(ProParameter* para);
	ProAtomicEnergyTable();

	double vdwEnergy(double d, double d0, double shift, double wd, double lamda);
	double vdwEnergyDesign(double d, double d0, double shift, double wd, double ldamda, double range);
	double vdwEnergyABACUS(double d, double d0, double kRep, double kAtr, double wd);

	float getHBEnergy(float d, float d0, double wd, double rescale){
		float e;
		if(d < d0)
		{
			float u = (d-d0)/d0/0.1;
			e = u*u+wd;
		}
		else
		{
			float u = (d-d0)/d0/0.15;
			if(u>4)
				e = 0.0;
			else
				e = wd*expf(-u*u);
		}
		if(e > 0)
			return e;
		else
			return e*rescale;
	}

	double getAtomEnergy(int atomIDA, const XYZ& tA, int atomIDB, const XYZ& tB, int sep, ProParameter* para);
	double getAtomEnergy(int atomIDA, const LocalFrame& csA, int atomIDB, const LocalFrame& csB, int sep, ProParameter* para);

	double getAtomEnergyDesign(int atomIDA, const LocalFrame& csA, int atomIDB, const LocalFrame& csB, int sep, float sai, DesignPara* para);

	double getAtomEnergyABACUS(int atomIDA, const LocalFrame& csA, int atomIDB, const LocalFrame& csB, int sep, float sai, DesignPara* para);
	double printDetailEnergyABACUS(int atomIDA, const LocalFrame& csA, int atomIDB, const LocalFrame& csB, int sep, float sai, DesignPara* para);


	double disulfideBondEnergy(double d, double ang1, double ang2, double dihed1, double dihed2, double dihed3, ProParameter* para);
	double disulfideBondEnergy(double d, double ang1, double ang2, double dihed1, double dihed2, double dihed3, DesignPara* para);

	virtual ~ProAtomicEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_PROATOMICENERGYTABLE_H_ */
