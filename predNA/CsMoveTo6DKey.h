/*
 * CsMoveTo6DKey.h
 *
 */

#ifndef predNA_CSMOVETO6DKEY_H_
#define predNA_CSMOVETO6DKEY_H_

#include <vector>
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include <fstream>
#include <map>
#include "geometry/CsMove.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"


namespace NSPpredNA {
using namespace std;
using namespace NSPgeometry;

class CsMoveTo6DKey {

private:
	map<int, int> sphereKeyMap;

public:

	CsMoveTo6DKey();
	string toKey(const LocalFrame& csA, const LocalFrame& csB);

	virtual ~CsMoveTo6DKey();
};

} /* namespace NSPforcefield */

#endif /* predNA_CSMOVETO6DKEY_H_ */
