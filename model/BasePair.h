/*
 * BasePair.h
 *
 */

#ifndef MODEL_BASEPAIR_H_
#define MODEL_BASEPAIR_H_

#include "model/BaseDistanceMatrix.h"
#include "model/AtomLib.h"
#include "model/StructureModel.h"

namespace NSPmodel {

class BasePair {
public:
	RNABase* baseA;
	RNABase* baseB;
	BaseDistanceMatrix dm;
	string type;
	bool isHbondPair;
	int hbNum;

	BasePair(RNABase* baseA, RNABase* baseB, AtomLib* atLib);
	bool isWCPair();
	bool isHbondedPair();
	double distanceToWCPair();

	virtual ~BasePair();
};

} /* namespace NSPmodel */

#endif /* MODEL_BASEPAIR_H_ */
