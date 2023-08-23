/*
 * BaseMoveLibrary.cpp
 *
 */

#include "predNA/BaseMoveLibrary.h"

namespace NSPpredna {

BaseMoveLibrary::BaseMoveLibrary() {
	// TODO Auto-generated constructor stub

	for(int i=0;i<100000;i++) {
		TransMatrix tm;
		tm = tm.randomRot(2.0);
		double rx = 0.2*(1.0*rand()/RAND_MAX - 0.5);
		double ry = 0.2*(1.0*rand()/RAND_MAX - 0.5);
		double rz = 0.2*(1.0*rand()/RAND_MAX - 0.5);
		XYZ ori(rx,ry,rz);
		this->randMoveList1.push_back(CsMove(ori,tm));
	}

	for(int i=0;i<100000;i++) {
		TransMatrix tm;
		tm = tm.randomRot(5.0);
		double rx = 0.5*(1.0*rand()/RAND_MAX - 0.5);
		double ry = 0.5*(1.0*rand()/RAND_MAX - 0.5);
		double rz = 0.5*(1.0*rand()/RAND_MAX - 0.5);
		XYZ ori(rx,ry,rz);
		this->randMoveList2.push_back(CsMove(ori,tm));
	}

	this->randMoveNum = randMoveList1.size();
}



BaseMoveLibrary::~BaseMoveLibrary() {

}

} /* namespace NSPpred */
