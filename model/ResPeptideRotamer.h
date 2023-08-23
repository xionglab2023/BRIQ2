/*
 * ResPeptideRotamer.h
 *
 */

#ifndef MODEL_RESPEPTIDEROTAMER_H_
#define MODEL_RESPEPTIDEROTAMER_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "model/DistanceMatrixHbond.h"
#include "model/StructureModel.h"

namespace NSPmodel {


class ResPeptideRotamer {
public:
	int typeA;
	int typeB;
	double psi;
	double phi;
	XYZ tList[5]; //pre-CA, pre-C, pre-O, N, CA
	CsMove cm31; //preCs3 -> curCs1, sequential connection
	CsMove cm13; //curCs1 -> preCs3, reverse connection
	DistanceMatrixHbond dm;
	int rotID;

	ResPeptideRotamer(const string& line);
	void setRotID(int id) {
		this->rotID = id;
	}
	ResPeptideRotamer& operator=(const ResPeptideRotamer& other){
		this->typeA = other.typeA;
		this->typeB = other.typeB;
		this->psi = other.psi;
		this->phi = other.phi;
		for(int i=0;i<5;i++){
			this->tList[i] = other.tList[i];
		}
		this->cm31 = other.cm31;
		this->cm13 = other.cm13;
		this->dm = other.dm;
		this->rotID = other.rotID;
		return *this;
	}

	double distanceTo(ResPeptideRotamer& other){
		return this->dm.distanceTo(other.dm);
	}

	virtual ~ResPeptideRotamer();
};

}
#endif /* MODEL_RESPEPTIDEROTAMER_H_ */
