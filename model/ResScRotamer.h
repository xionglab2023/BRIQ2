/*
 * ResScRotamer.h
 *
 */

#ifndef MODEL_RESSCROTAMER_H_
#define MODEL_RESSCROTAMER_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/AtomLib.h"
#include "model/ResName.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "model/StructureModel.h"

namespace NSPmodel {


class ResScRotamer {
public:

	int aaType;
	int rotID;
	int atomNum;

	XYZ coordsLocal[10];
	int uniqueIDs[10];

	int polarAtomNum;
	CsMove polarCmList[3];
	int polarAtomIndex[3];


	ResScRotamer();
	ResScRotamer(const string& line, AtomLib* atLib);
	ResScRotamer(Residue* res, AtomLib* atLib);
	double distanceTo(ResScRotamer* other);

	virtual ~ResScRotamer();

};


class ResScConformer {
public:

	int aaType;
	int atomNum;
	XYZ coords[10];

	int polarAtomNum;
	LocalFrame polarCsList[3];

	ResScRotamer* rot;
	LocalFrame cs; //N, CA, C, origin at CA

	ResScConformer();
	void init(ResScRotamer* rot, const LocalFrame& cs);
	void updateRotamer(ResScRotamer* rot);
	void updateCs(const LocalFrame& cs);
	void copyValueFrom(ResScConformer* other);
	bool equals(ResScConformer* other);

	virtual ~ResScConformer();
};

}
#endif /* MODEL_RESSCROTAMER_H_ */
