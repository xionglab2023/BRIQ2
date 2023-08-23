/*
 * AssignRNASS.h
 *
 */

#ifndef MODEL_ASSIGNRNASS_H_
#define MODEL_ASSIGNRNASS_H_

#include "model/BasePair.h"
#include "model/BaseDistanceMatrix.h"
#include "geometry/localframe.h"
#include "model/StructureModel.h"

namespace NSPmodel {

class AssignRNASS {
public:
	vector<RNABase*> baseList;
	int len;
	int* pairIndex;
	double* pairDistanceToWC;
	int* nwcPairIndex;
	int* nwcPairHbondNum;

	vector<int> breakList;
	bool* connectWithNeighbor;
	string ssSeq;
	string nwcSeq;
	string seq;

	AssignRNASS(RNAPDB* pdb, AtomLib* atLib);
	string indexToBractString(int* index);
	void printInfo(const string& outFile);
	virtual ~AssignRNASS();
};

} /* namespace NSPforcefield */

#endif /* MODEL_ASSIGNRNASS_H_ */
