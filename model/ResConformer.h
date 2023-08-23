/*
 * ResConformer.h
 *
 *  Created on: 2022Äê3ÔÂ3ÈÕ
 *      Author: pengx
 */

#ifndef MODEL_RESCONFORMER_H_
#define MODEL_RESCONFORMER_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "model/ResBBRotamer.h"
#include "model/ResScRotamer.h"
#include "model/StructureModel.h"

namespace NSPmodel {

using namespace NSPgeometry;
using namespace std;

class ResConformer {
public:

	int aaType;
	ResBBConformer* bbConf;
	ResScConformer* scConf;

	ResConformer();
	ResConformer(int aaType);

	void init(ResScRotamer* rotSc, ResBBRotamer* rotBb, const LocalFrame& cs2);
	void copyValueFrom(ResConformer* other);
	void updateScRotamer(ResScRotamer* rotSc);


	virtual ~ResConformer();
};

} /* namespace NSPmodel */

#endif /* MODEL_RESCONFORMER_H_ */
