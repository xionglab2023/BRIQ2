/*
 * PhosphateRotamer.h
 *
 *  Created on: 2023Äê8ÔÂ14ÈÕ
 *      Author: nuc
 */

#ifndef MODEL_PHOSPHATEROTAMER_H_
#define MODEL_PHOSPHATEROTAMER_H_

#include <vector>
#include <map>
#include <fstream>
#include "geometry/CsMove.h"
#include "geometry/localframe.h"
#include "tools/StringTool.h"
#include "model/StructureModel.h"

namespace NSPmodel {
using namespace NSPgeometry;
using namespace NSPtools;

class PhosphateRotamer {
public:

	/*
	 * atom order:
	 * p
	 * O5'
	 * OP1
	 * OP2
	 */

	XYZ localCoords[4];
	CsMove cmOP1;
	CsMove cmOP2;
	double dihed1;
	double dihed2;

	PhosphateRotamer(double d1, double d2){
		this->dihed1 = d1;
		this->dihed2 = d2;
		double len1 = 1.605;
		double len2 = 1.592;
		double ang1 = 120.1;
		double ang2 = 103.5;
		LocalFrame cs0;
		LocalFrame cs1 = cs0.csNext(len1, ang1, d1);
		LocalFrame cs2 = cs1.csNext(len2, ang2, d2);
		XYZ t0 = XYZ(0.0, 0.0, 0.0);
		localCoords[0] = cs1.origin_;
		localCoords[1] = cs2.origin_;
		localCoords[2] = local2global(cs2, XYZ(-2.062, 0.570, 1.278));
		localCoords[3] = local2global(cs2, XYZ(-2.054, 0.589, -1.275));

		LocalFrame op1(t0, localCoords[0], localCoords[2]);
		LocalFrame op2(t0, localCoords[0], localCoords[3]);

		cmOP1 = op1 - cs0;
		cmOP2 = op2 - cs0;
	}

	PhosphateRotamer(RNABase* base, RNABase* baseNext) {
		Atom* a = base->getAtom("C2'");
		Atom* b = base->getAtom("C3'");
		Atom* c = base->getAtom("O3'");
		Atom* d = baseNext->getAtom("P");
		Atom* e = baseNext->getAtom("O5'");
		double d1, d2;
		if(a == NULL || b == NULL || c == NULL || d == NULL || e == NULL)
		{
			d1 = 0.0;
			d2 = 0.0;
		}
		else {
			d1 = dihedral(a->coord, b->coord, c->coord, d->coord);
			d2 = dihedral(b->coord, c->coord, d->coord, e->coord);
		}
		this->dihed1 = d1;
		this->dihed2 = d2;
		double len1 = 1.605;
		double len2 = 1.592;
		double ang1 = 120.1;
		double ang2 = 103.5;
		LocalFrame cs0;
		LocalFrame cs1 = cs0.csNext(len1, ang1, d1);
		LocalFrame cs2 = cs1.csNext(len2, ang2, d2);
		XYZ t0 = XYZ(0.0, 0.0, 0.0);
		localCoords[0] = cs1.origin_;
		localCoords[1] = cs2.origin_;
		localCoords[2] = local2global(cs2, XYZ(-2.062, 0.570, 1.278));
		localCoords[3] = local2global(cs2, XYZ(-2.054, 0.589, -1.275));

		LocalFrame op1(t0, localCoords[0], localCoords[2]);
		LocalFrame op2(t0, localCoords[0], localCoords[3]);

		cmOP1 = op1 - cs0;
		cmOP2 = op2 - cs0;
	}

	PhosphateRotamer& operator=(const PhosphateRotamer& other){
		for(int i=0;i<4;i++){
			this->localCoords[i] = other.localCoords[i];
		}
		this->cmOP1 = other.cmOP1;
		this->cmOP2 = other.cmOP2;
		this->dihed1 = other.dihed1;
		this->dihed2 = other.dihed2;
		return *this;
	}

	bool equalTo(PhosphateRotamer* other){
		if(squareDistance(localCoords[0], other->localCoords[0]) > 0.000001) return false;
		if(squareDistance(localCoords[1], other->localCoords[1]) > 0.000001) return false;
		if(squareDistance(localCoords[2], other->localCoords[2]) > 0.000001) return false;
		if(squareDistance(localCoords[3], other->localCoords[3]) > 0.000001) return false;
		return true;
	}

	virtual ~PhosphateRotamer();
};

class PhosphateConformer{

public:
	PhosphateRotamer* rot;
	LocalFrame cs2;

	XYZ coords[4];
	LocalFrame op1Polar;
	LocalFrame op2Polar;
	double ene;

	PhosphateConformer();
	PhosphateConformer(PhosphateRotamer* rot, LocalFrame& cs2);

	void setEnergy(double e){
		this->ene = e;
	}
	void copyValueFrom(PhosphateConformer* other);
	void updateLocalFrameAndRotamer(LocalFrame& cs2, PhosphateRotamer* rot, double ene);
	void updateLocalFrame(LocalFrame& cs2);
	void updateRotamer(PhosphateRotamer* rot);
	bool equalTo(PhosphateConformer* other){
		if(squareDistance(coords[0], other->coords[0]) > 0.000001) return false;
		if(squareDistance(coords[1], other->coords[1]) > 0.000001) return false;
		if(squareDistance(coords[2], other->coords[2]) > 0.000001) return false;
		if(squareDistance(coords[3], other->coords[3]) > 0.000001) return false;
		return true;
	}
	double distanceTo(PhosphateConformer* other);
	virtual ~PhosphateConformer();

};

} /* namespace NSPmodel */

#endif /* MODEL_PHOSPHATEROTAMER_H_ */
