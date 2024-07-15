/*
 * BaseRotamer.h
 *
 */

#ifndef MODEL_BASEROTAMER_H_
#define MODEL_BASEROTAMER_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include <vector>
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/AtomLib.h"
#include "model/ResName.h"
#include "model/StructureModel.h"
#include "tools/StringTool.h"

namespace NSPmodel {

using namespace NSPgeometry;

class BaseRotamer{
public:
	int baseType;
	int atomNum;
	XYZ coordsLocal[11];
	int uniqueIDs[11];

	int polarAtomNum;
	CsMove polarCmList[5];
	int polarAtomIndex[5];
	int polarAtomUniqueID[5];

	BaseRotamer(int baseType, AtomLib* atLib);
	virtual ~BaseRotamer();
};

class BaseRotamerCG{
public:
	int baseType;
	XYZ coordsLocal[3];
	int uniqueIDs[3];

	BaseRotamerCG(int baseType, AtomLib* atLib);

	void printInfo(){
		printf("%d \n", baseType);
		printf("%7.3f %7.3f %7.3f\n", coordsLocal[0].x_,coordsLocal[0].y_, coordsLocal[0].z_);
		printf("%7.3f %7.3f %7.3f\n", coordsLocal[1].x_,coordsLocal[1].y_, coordsLocal[1].z_);
		printf("%7.3f %7.3f %7.3f\n", coordsLocal[2].x_,coordsLocal[2].y_, coordsLocal[2].z_);
	}

	virtual ~BaseRotamerCG();
};


class BaseConformer{
public:

	BaseRotamer* rot;
	XYZ coords[11];
	LocalFrame cs1;
	LocalFrame csPolar[5];

	BaseConformer(BaseRotamer* rot, LocalFrame& cs);
	void copyValueFrom(BaseConformer* other);
	void updateCoords(LocalFrame& cs);
	double distanceTo(BaseConformer* other);
	virtual ~BaseConformer();
};

class BaseConformerCG{
public:
	BaseRotamerCG* rot;
	XYZ coords[3];
	LocalFrame cs1;

	BaseConformerCG(BaseRotamerCG* rot, LocalFrame& cs);
	void copyValueFrom(BaseConformerCG* other);
	void updateCoords(LocalFrame& cs);
	double distanceTo(BaseConformerCG* other);
	virtual ~BaseConformerCG();
};

} /* namespace NSPforcefield */

#endif /* MODEL_BASEROTAMER_H_ */
