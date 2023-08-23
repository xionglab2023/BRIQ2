/*
 * ResBBRotamer.h
 *
 */

#ifndef MODEL_RESBBROTAMER_H_
#define MODEL_RESBBROTAMER_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/ResName.h"
#include "model/DistanceMatrixHbond.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "model/AtomLib.h"
#include "model/StructureModel.h"

namespace NSPmodel {

using namespace NSPgeometry;
using namespace std;

class ResBBRotamer {
public:

	XYZ localTList[5]; //pre_C, N, CA, C, O, coordinate in LocalFrame2

	int uniqueID[4];
	//Localframe1: N
	//Localframe2: CA
	//Localframe3: C
	//Localframe4: O

	CsMove cm21;
	CsMove cm23;
	CsMove cm12;
	CsMove cm32;
	CsMove cm13;
	CsMove cm31;
	CsMove cm24;

	int aaType;
	double phi;
	double psi;

	double ene;

	int index1W; //index in 10000 rotamer lib
	int index1K; //index in 1000 rotamer lib

	DistanceMatrixHbond dm;

	ResBBRotamer(int type, const XYZ& preC, const XYZ& N, const XYZ& CA, const XYZ& C, const XYZ& O, AtomLib* atLib){
		this->aaType = type;

		LocalFrame csX = generateLocalFrameResidueStyle(N, CA, C);
		localTList[0] = global2local(csX, preC);
		localTList[1] = global2local(csX, N);
		localTList[2] = global2local(csX, CA);
		localTList[3] = global2local(csX, C);
		localTList[4] = global2local(csX, O);

		LocalFrame cs1 = generateLocalFrameResidueStyle(localTList[0], localTList[1], localTList[2]);
		LocalFrame cs2;
		LocalFrame cs3 = generateLocalFrameResidueStyle(localTList[2], localTList[3], localTList[4]);
		LocalFrame cs4 = LocalFrame(localTList[2], localTList[3], localTList[4]); //O

		dm = DistanceMatrixHbond(cs1, cs3);

		cm21 = cs1 - cs2;
		cm23 = cs3 - cs2;
		cm12 = cs2 - cs1;
		cm32 = cs2 - cs3;
		cm13 = cs3 - cs1;
		cm31 = cs1 - cs3;
		cm24 = cs4 - cs2;

		this->phi = dihedral(localTList[0], localTList[1], localTList[2], localTList[3]);
		this->psi = dihedral(localTList[1], localTList[2], localTList[3], localTList[4]);
		psi = psi + 180;
		if(psi > 180)
			psi = psi - 360;

		this->ene = 0.0;
		this->index1K = 0;
		this->index1W = 0;

		for(int i=0;i<4;i++)
		{
			this->uniqueID[i] = atLib->aaUniqueIDs[aaType]->at(i);
		}
	}

	ResBBRotamer(Residue* resP, Residue* res, AtomLib* atLib);
	ResBBRotamer(Residue* res, AtomLib* atLib);
	ResBBRotamer(const string& line, AtomLib* atLib);
	void setAAType(int aa, AtomLib* atLib);
	double distanceTo(ResBBRotamer* other){
		return this->dm.distanceTo(other->dm);
	}
	string toString();

	virtual ~ResBBRotamer();
};

class ResBBConformer {
public:
	XYZ coords[5]; //pre_C, N, CA, C, O

	LocalFrame cs1; //pre_C, N, CA, origin at N  , polar CsN
	LocalFrame cs2; //N, CA, C, origin at CA
	LocalFrame cs3; //CA, C, O, origin at C, polar CsO
	LocalFrame cs4; //CA, C, O, origin at O

	ResBBRotamer* rot;
	int aaType;

	ResBBConformer(){
		this->aaType = 0;
		this->rot = NULL;
	}

	void init(ResBBRotamer* rot, const LocalFrame& cs2);

	void updateRotFixCs1(ResBBRotamer* rot);
	void updateRotFixCs2(ResBBRotamer* rot);
	void updateRotFixCs3(ResBBRotamer* rot);

	void updateCs1(LocalFrame& cs1);
	void updateCs2(LocalFrame& cs2);
	void updateCs3(LocalFrame& cs3);

	void copyValueFrom(ResBBConformer* other);

	virtual ~ResBBConformer();
};

}
#endif /* MODEL_RESBBROTAMER_H_ */
