/*
 * DistanceMatrixHbond.h
 *
 */

#ifndef MODEL_DISTANCEMATRIXHBOND_H_
#define MODEL_DISTANCEMATRIXHBOND_H_
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "tools/StringTool.h"
#include "geometry/TransMatrix.h"
#include <string>
#include "model/StructureModel.h"

namespace NSPmodel {

class DistanceMatrixHbond {
public:
	double dm[16];
	DistanceMatrixHbond();
	DistanceMatrixHbond(LocalFrame& csA, LocalFrame& csB);
	DistanceMatrixHbond(const string& line);

	DistanceMatrixHbond& operator=(const DistanceMatrixHbond& other){
		for(int i=0;i<16;i++){
			this->dm[i] = other.dm[i];
		}
		return *this;
	}

	double distanceTo(DistanceMatrixHbond& other);
	double squareDistance(DistanceMatrixHbond& other);

	virtual ~DistanceMatrixHbond();
};

} /* namespace NSPmodel */

#endif /* MODEL_DISTANCEMATRIXHBOND_H_ */
