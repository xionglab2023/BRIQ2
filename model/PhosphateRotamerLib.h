/*
 * PhosphateRotamerLib.h
 *
 *  Created on: 2023Äê8ÔÂ14ÈÕ
 *      Author: nuc
 */

#ifndef MODEL_PHOSPHATEROTAMERLIB_H_
#define MODEL_PHOSPHATEROTAMERLIB_H_


#include <time.h>
#include "model/StructureModel.h"
#include "model/PhosphateRotamer.h"


namespace NSPmodel {

class PhosphateRotamerLib {
public:

	PhosphateRotamer* prLib[360][360];

	PhosphateRotamerLib();
	virtual ~PhosphateRotamerLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_PHOSPHATEROTAMERLIB_H_ */
