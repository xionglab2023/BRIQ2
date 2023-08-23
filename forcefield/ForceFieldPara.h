/*
 * ForceFieldPara.h
 *
 *  Created on: 2023Äê4ÔÂ25ÈÕ
 *      Author: pengx
 */

#ifndef FORCEFIELD_FORCEFIELDPARA_H_
#define FORCEFIELD_FORCEFIELDPARA_H_

namespace NSPforcefield {

class ForceFieldPara {
public:

	double clashLamda;
	double clashShift;

	ForceFieldPara();
	virtual ~ForceFieldPara();
};

} /* namespace NSPmodel */

#endif /* FORCEFIELD_FORCEFIELDPARA_H_ */
