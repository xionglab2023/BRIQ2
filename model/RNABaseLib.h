/*
 * RNABaseLib.h
 *
 */

#ifndef MODEL_RNABASELIB_H_
#define MODEL_RNABASELIB_H_

#include "model/AtomLib.h"
#include "geometry/localframe.h"
#include "geometry/RMSD.h"
#include "model/RiboseRotamer.h"
#include "model/StructureModel.h"
#include "model/PhosphateRotamer.h"

namespace NSPmodel {

using namespace NSPgeometry;

class RNABaseLib {
public:
	AtomLib atLib;
	RNABaseLib();

	vector<XYZ> getBaseSidechainCoords(int baseType);
	vector<XYZ> getPolarAtomCoords(int baseType);
	RNABase* getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs);
	RNABase* getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs, RiboseRotamer* rot);
	RNABase* getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs, RiboseRotamer* rot, PhosphateRotamer* pRot);

	RNABase* toStandardBase(RNABase* base);

	virtual ~RNABaseLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_RNABASELIB_H_ */
