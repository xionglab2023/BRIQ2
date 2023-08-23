/*
 * DistanceMatrixRes.h
 *
 */

#ifndef MODEL_DISTANCEMATRIXRES_H_
#define MODEL_DISTANCEMATRIXRES_H_
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "tools/StringTool.h"
#include "geometry/TransMatrix.h"
#include <string>

namespace NSPmodel {

using namespace std;
using namespace NSPgeometry;

class DistanceMatrixRes {
public:
	double dm[16];
	DistanceMatrixRes();
	DistanceMatrixRes(LocalFrame& csA, LocalFrame& csB);
	DistanceMatrixRes(CsMove& cm);
	DistanceMatrixRes(const string& line);
	DistanceMatrixRes& operator=(const DistanceMatrixRes& other){
		for(int i=0;i<16;i++){
			this->dm[i] = other.dm[i];
		}
		return *this;
	}
	double distanceTo(DistanceMatrixRes& other);
	double squareDistance(DistanceMatrixRes& other);
	virtual ~DistanceMatrixRes();
};


} /* namespace NSPmodel */

#endif /* MODEL_DISTANCEMATRIXRES_H_ */
