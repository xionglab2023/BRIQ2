/*
 * HbondEnergy.h
 *
 *  Created on: 2023Äê8ÔÂ11ÈÕ
 *      Author: nuc
 */

#ifndef FORCEFIELD_HBONDENERGY_H_
#define FORCEFIELD_HBONDENERGY_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "geometry/localframe.h"
#include "model/AtomLib.h"
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
using namespace NSPmodel;


class HbondEnergy {
public:

	/*
	 * donor atom list
	 * 0: O2'
	 * 1: SER-OG
	 * 2: THR-OG1
	 * 3: TYR-OH
	 * 4: HIS-ND1 HIS-NE2
	 * 5: N
	 * 6: ASN-ND2 GLN-NE2
	 * 7: ARG-NH1 ARG-NH2
	 * 8: ARG-NE
	 * 9: LYS-NZ
	 * 10: TRP-NE1
	 * 11: A-N6 DA-N6
	 * 12: U-N3 DT-N3
	 * 13: G-N1 DG-N1
	 * 14: G-N2 DG-N2
	 * 15: C-N4 DC-N4
	 */

	/*
	 * acceptor atom list
	 * 0: O2'
	 * 1: SER-OG
	 * 2: THR-OG1
	 * 3: TYR-OH
	 * 4: HIS-ND1 HIS-NE2
	 * 5: O
	 * 6: ASP-OD1 ASP-OD2
	 * 7: GLU-OE1 GLU-OE2
	 * 8: ASN-OD1 GLN-OE1
	 * 9: OP1 OP2
	 * 10: A-N7 DA-N7
	 * 11: A-N1 DA-N1
	 * 12: A-N3 DA-N3
	 * 13: U-O2 U-O4 DT-O2 DT-O4
	 * 14: G-N7 DG-N7
	 * 15: G-O6 DG-O6
	 * 16: G-N3 DG-N3
	 * 17: C-O2 DC-O2
	 * 18: C-N3 DC-N3
	 */

	//In total 304 hbond types (16*19)

	double hbEne[16][19];
	double hbDist[16][19];

	double donorSphereLogDensity[16][19][5000];
	double acceptorSphereLogDensity[16][19][5000];

	//sphere point number: 5000
	map<int, int> sphereKeyMap;
	map<int, int>::iterator it;
	map<int, int> keyIdToListID; //key id range: 0~400*400*400
	int spIndex1[753848];
	int spIndex2[753848];
	int spIndex3[753848];
	double spWt1[753848];
	double spWt2[753848];
	double spWt3[753848];

	AtomLib* atLib;
	ForceFieldPara* para;

	HbondEnergy(ForceFieldPara* para);
	double getEnergy(int uniqueIDA, LocalFrame& csA, int uniqueIDB, LocalFrame& csB);
	double getO4O2C2Energy(double distance, int sep);
	virtual ~HbondEnergy();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_HBONDENERGY_H_ */
