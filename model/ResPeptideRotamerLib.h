/*
 * ResPeptideRotamerLib.h
 *
 */

#ifndef MODEL_RESPEPTIDEROTAMERLIB_H_
#define MODEL_RESPEPTIDEROTAMERLIB_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "model/ResPeptideRotamer.h"
#include "dataio/datapaths.h"
#include "model/StructureModel.h"
#include "tools/StringTool.h"

namespace NSPmodel {

using namespace NSPgeometry;
using namespace std;

class ResPeptideRotamerLib {
public:

	ResPeptideRotamer* pepLib[32000]; //20*20*80
	float ene[32000][200]; //20*20*80*200

	float ppLibPsi[200]; //200 rep points for psi-phi
	float ppLibPhi[200];

	int psiphiToPPIndex[130321]; //361*361
	ResPeptideRotamerLib();

	ResPeptideRotamer* getRandomPeptide(int typeA, int typeB){
		return pepLib[typeA*1600+typeB*80+rand()%80];
	}

	double getEnergy(int typeA, int typeB, int rotID, double psi, double phi){
		if(psi < -180 || psi >180 || phi < -180 || phi > 180){
			cerr << "invalid input for peptide rotamer energy "  << psi << " " << phi << endl;
			exit(1);
		}
		int id1 = int(psi + 180);
		int id2 = int(phi + 180);
		return ene[typeA*1600 + typeB*80 + rotID][psiphiToPPIndex[id1*361+id2]];
	}


	virtual ~ResPeptideRotamerLib();
};

}
#endif /* MODEL_RESPEPTIDEROTAMERLIB_H_ */
