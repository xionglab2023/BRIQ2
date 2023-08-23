/*
 * BaseRotamerLib.h
 *
 */

#ifndef MODEL_BASEROTAMERLIB_H_
#define MODEL_BASEROTAMERLIB_H_

#include <time.h>
#include "model/StructureModel.h"
#include "model/BaseRotamer.h"

namespace NSPmodel {


class BaseRotamerLib {
public:
	BaseRotamer* baseLib[8];
	BaseRotamerLib(AtomLib* atLib);
	void testRotamerLib(AtomLib* atLib);
	virtual ~BaseRotamerLib();
};


} /* namespace NSPforcefield */

#endif /* MODEL_BASEROTAMERLIB_H_ */
