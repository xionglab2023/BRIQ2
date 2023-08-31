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

} /* namespace NSPforcefield */

#endif /* MODEL_BASEROTAMER_H_ */
