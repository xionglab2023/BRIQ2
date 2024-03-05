/*
 * PO3Builder.h
 *
 */

#ifndef FORCEFIELD_PO3BUILDER_H_
#define FORCEFIELD_PO3BUILDER_H_

#include <vector>
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "model/PhosphateRotamer.h"
#include "model/PhosphateRotamerLib.h"
#include "model/RiboseRotamer.h"
#include "forcefield/ForceFieldPara.h"
#include "geometry/Angles.h"
#include <fstream>
#include <vector>
#include "geometry/CsMove.h"
#include "dataio/datapaths.h"
#include "RnaAtomicEnergyTable.h"
#include "tools/StringTool.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;


class PO3Builder {
public:

	PhosphateRotamerLib* rotLib;

	vector<XYZ> pList;
	vector<XYZ> op1List;
	vector<XYZ> op2List;
	vector<XYZ> o5List;

	double len1 = 1.605;
	double len2 = 1.592;
	double len3 = 1.422;
	double ang1 = 120.1;
	double ang2 = 103.5;
	double ang3 = 120.7;
	double ang4 = 111.1;

	vector<double> eImpD1D2;
	/*
	 * six regions: regionIndexA
	 * R0: impIndexA < 10 && D2 >  250
	 * R1: impIndexA < 10 && D2 <= 250 && D1 <  210
	 * R2: impIndexA < 10 && D2 <= 250 && D1 >= 210
	 * R3: impIndexA >= 10 && D2 >  250
	 * R4: impIndexA >= 10 && D2 <= 250 && D1 <  210
	 * R5: impIndexA >= 10 && D2 <= 250 && D1 >= 210
	 */


	vector<double> eImpD4D5;
	/*
	 *six regions: regionIndexB
	 * R0: impIndexB < 10 && D5 >  240
	 * R1: impIndexB < 10 && D5 > 120 && D5 <= 240
	 * R2: impIndexB < 10 && D5 <= 120
	 * R3: impIndexB >= 10 && D5 >  240
	 * R4: impIndexB >= 10 && D5 > 120 && D5 <= 240
	 * R5: impIndexB >= 10 && D5 <= 120
	 */


	vector<double> eD2D4D3;
	/*
	 * six regions: regionIndexC
	 * R0: D2 > 250 && D3 > 240
	 * R1: D2 > 250 && D3 > 120 && D3 <= 240
	 * R2: D2 > 250 && D3 <= 120
	 * R3: D2 <= 250 && D3 > 240
	 * R4: D2 <= 250 && D3 > 120 && D3 <= 240
	 * R5: D2 <= 250 && D3 <= 120
	 */

	vector<XYZ> d1d2Lib1A; //dihed1-dihed2 library level 1, impA < 0
	vector<vector<XYZ>> d1d2Lib2A; //dihed1-dihed2 library level 2 impA < 0

	vector<XYZ> d1d2Lib1B; //dihed1-dihed2 library level 1, impA > 0
	vector<vector<XYZ>> d1d2Lib2B; //dihed1-dihed2 library level 2, impA > 0

	vector<double> lib2ErrorA;
	vector<double> lib2ErrorB;

	vector<double> eDihed1;
	vector<double> eDihed2;
	vector<double> eDihed3;
	vector<double> eDihed4;
	vector<double> eDihed5;

	ForceFieldPara* para;

	PO3Builder(ForceFieldPara* para);

	void buildPhosphate(RiboseConformer* riboConfA, RiboseConformer* riboConfB, PhosphateConformer* outPhoConf);
	double getEnergy(RiboseConformer* riboConfA, RiboseConformer* riboConfB);
	double getEnergyFast(RiboseConformer* riboConfA, RiboseConformer* riboConfB);

	virtual ~PO3Builder();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_PO3BUILDER_H_ */
