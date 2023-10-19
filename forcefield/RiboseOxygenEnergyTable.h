/*
 * RiboseOxygenEnergyTable.h
 *
 */

#ifndef FORCEFIELD_RIBOSEOXYGENENERGYTABLE_H_
#define FORCEFIELD_RIBOSEOXYGENENERGYTABLE_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "dataio/datapaths.h"
#include "forcefield/ForceFieldPara.h"
#include "tools/StringTool.h"
#include <time.h>
#include <math.h>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;

class RiboseOxygenEnergyTable {
public:

	/* Base type
	 * 0: A
	 * 1: U
	 * 2: G
	 * 3: C
	 * 4: a
	 * 5: t
	 * 6: g
	 * 7: c
	 * Oxygen atom type
	 * 0: O3'
	 * 1: O4'
	 * 2: O5'
	 */

	vector<double> etList[24];
	double wtList1[24];
	/*
	 * sep=-1
	 * base-O4
	 */

	vector<double> etListM1[8];
	double wtList2[8];

	RiboseOxygenEnergyTable(ForceFieldPara* para);

	double getEnergy(int baseType, int oxygenType, const XYZ& localCoreCoord, int sep){
		if(sep == 0 || sep == 1)
			return 0.0;
		else if(sep == -1){
			if(oxygenType == 1){ //O4'
				if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
				if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
				if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
				int idX = (int)((localCoreCoord.x_+3.0)*5);
				int idY = (int)((localCoreCoord.y_+6.0)*5);
				int idZ = (int)((localCoreCoord.z_+5.0)*5);
				int type = baseType;
				int index = idX*4000+idY*50+idZ;
				return etListM1[type][index];
			}
			else
				return 0.0;
		}
		else {
			if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
			if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
			if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
			int idX = (int)((localCoreCoord.x_+3.0)*5);
			int idY = (int)((localCoreCoord.y_+6.0)*5);
			int idZ = (int)((localCoreCoord.z_+5.0)*5);
			int type = baseType*3+oxygenType;
			int index = idX*4000+idY*50+idZ;

			return etList[type][index];
		}
	}

	inline double riboseOxyEnergyRescale(double e){
		if(e < 0.81)
			return e;
		else
			return 1.8*sqrt(e) - 0.81;
	}

	virtual ~RiboseOxygenEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_RIBOSEOXYGENENERGYTABLE_H_ */
