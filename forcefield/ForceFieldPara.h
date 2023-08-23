/*
 * ForceFieldPara.h
 *
 *  Created on: 2023��4��25��
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
