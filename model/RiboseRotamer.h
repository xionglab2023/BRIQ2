/*
 * RiboseRotamer.h
 *
 *  Created on: 2022Äê9ÔÂ6ÈÕ
 *      Author: pengx
 */

#ifndef MODEL_RIBOSEROTAMER_H_
#define MODEL_RIBOSEROTAMER_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include <vector>

#include "model/StructureModel.h"
#include "tools/StringTool.h"

namespace NSPmodel {

using namespace NSPgeometry;

class RiboseRotamer {

	/*
	 * Atom Order for RNA: (old order)
	 * 0:  C1'
	 * 1:  C2'
	 * 2:  C3'
	 * 3:  C4'
	 * 4:  O4'
	 * 5:  O2'
	 * 6:  O3'
	 * 7:  C5'
	 *
	 * Atom Order for RNA: (new order)
	 * 0:  C1'
	 * 1:  C2'
	 * 2:  C3'
	 * 3:  C4'
	 * 4:  O4'
	 * 5:  O3'
	 * 6:  C5'
	 * 7:  O2'
	 *
	 * Atom Order for DNA:
	 * 0:  C1'
	 * 1:  C2'
	 * 2:  C3'
	 * 3:  C4'
	 * 4:  O4'
	 * 5:  O3'
	 * 6:  C5'
	 *
	 * cs2: C2'-C3'-O3'
	 * cs3: O4'-C4'-C5'
	 */


public:

	int resType;
	int rotTypeLv1;
	int rotType;

	int atomNum;
	XYZ localCoords[8]; //local coordinate in cs1
	XYZ center;
	CsMove mv12;
	CsMove mv13;
	CsMove mv31;
	CsMove mv21;

	CsMove mv1O2;

	//local frame X: generate by atom N9(N1)-C1'-C4'
	CsMove mv1X;
	CsMove mvX1;

	double chi;
	double improper;
	double energy;

	RiboseRotamer(){
		this->resType = 0;
		this->atomNum = 8;
		this->rotTypeLv1 = 0;
		this->rotType = 0;
		this->energy = 0;
		this->improper = 0;
		this->chi = 0;
	}

	RiboseRotamer(RNABase* base);
	RiboseRotamer(const string& line);


	RiboseRotamer& operator=(const RiboseRotamer& other){
		this->resType = other.resType;
		this->rotTypeLv1 = other.rotTypeLv1;
		this->rotType = other.rotType;
		this->atomNum = other.atomNum;
		for(int i=0;i<atomNum;i++){
			localCoords[i] = other.localCoords[i];
		}
		center = other.center;
		mv12 = other.mv12;
		mv13 = other.mv13;
		mv31 = other.mv31;
		mv21 = other.mv21;
		mv1O2 = other.mv1O2;
		mv1X = other.mv1X;
		mvX1 = other.mvX1;
		energy = other.energy;
		return *this;
	}

	void copyValueFrom(RiboseRotamer* other){
		this->resType = other->resType;
		this->rotTypeLv1 = other->rotTypeLv1;
		this->rotType = other->rotType;
		this->atomNum = other->atomNum;
		for(int i=0;i<atomNum;i++){
			localCoords[i] = other->localCoords[i];
		}
		center = other->center;
		mv12 = other->mv12;
		mv13 = other->mv13;
		mv31 = other->mv31;
		mv21 = other->mv21;
		mv1O2 = other->mv1O2;
		mv1X = other->mv1X;
		mvX1 = other->mvX1;
		energy = other->energy;
	}


	bool equalTo(const RiboseRotamer& other){
		if(this->resType != other.resType) return false;
		if(this->rotTypeLv1 != other.rotTypeLv1) return false;
		if(this->rotType != other.rotType) return false;
		for(int i=0;i<atomNum;i++){
			if(squareDistance(localCoords[i], other.localCoords[i]) > 0.0001) return false;
		}
		if(squareDistance(center, other.center) > 0.0001) return false;
		return true;
	}

	double distanceTo(RiboseRotamer* other){
		double d = 0;
		for(int i=0;i<atomNum;i++){
			d += localCoords[i].squaredDistance(other->localCoords[i]);
		}
		return sqrt(d/atomNum);
	}

	virtual ~RiboseRotamer();
};


class RiboseRotamerCG {

	/*
	 * Atom Order
	 * 0:  C1'
	 * 6:  O3'
	 * 7:  C5'
	 */


public:

	int resType;
	int rotTypeLv1;
	int rotType;

	XYZ localCoords[3]; //local coordinate in cs1
	double energy;

	RiboseRotamerCG(){
		this->resType = 0;
		this->rotTypeLv1 = 0;
		this->rotType = 0;
		this->energy = 0;
	}

	RiboseRotamerCG(RNABase* base);
	RiboseRotamerCG(const string& line);


	RiboseRotamerCG& operator=(const RiboseRotamerCG& other){
		this->resType = other.resType;
		this->rotTypeLv1 = other.rotTypeLv1;
		this->rotType = other.rotType;
		for(int i=0;i<3;i++){
			localCoords[i] = other.localCoords[i];
		}
		energy = other.energy;
		return *this;
	}

	void copyValueFrom(RiboseRotamerCG* other){
		this->resType = other->resType;
		this->rotTypeLv1 = other->rotTypeLv1;
		this->rotType = other->rotType;
		for(int i=0;i<3;i++){
			localCoords[i] = other->localCoords[i];
		}
		energy = other->energy;
	}


	bool equalTo(const RiboseRotamerCG& other){
		if(this->resType != other.resType) return false;
		if(this->rotTypeLv1 != other.rotTypeLv1) return false;
		if(this->rotType != other.rotType) return false;
		for(int i=0;i<3;i++){
			if(squareDistance(localCoords[i], other.localCoords[i]) > 0.0001) return false;
		}
		return true;
	}

	double distanceTo(RiboseRotamerCG* other){
		double d = 0;
		for(int i=0;i<3;i++){
			d += localCoords[i].squaredDistance(other->localCoords[i]);
		}
		return sqrt(d/3);
	}

	virtual ~RiboseRotamerCG();
};


class RiboseConformer {
public:
	RiboseRotamer* rot;

	XYZ coords[8];
	LocalFrame cs1;
	LocalFrame cs2; //C2'-C3'-O3'
	LocalFrame cs3; //O4'-C4'-C5'

	bool hasO2;

	LocalFrame o2Polar;

	RiboseConformer();
	RiboseConformer(RiboseRotamer* rot, LocalFrame& cs1);
	void copyValueFrom(RiboseConformer* other);
	void updateLocalFrame(LocalFrame& cs1);

	void updateRotamer(RiboseRotamer* rot);
	void updateRotamerCs2Fixed(RiboseRotamer* rot);
	void updateRotamerCs3Fixed(RiboseRotamer* rot);
	void updateLocalFrameAndRotamer(LocalFrame& cs1, RiboseRotamer* rot);

	double distanceTo(RiboseConformer* other);
	virtual ~RiboseConformer();
};

class RiboseConformerCG {
public:
	RiboseRotamerCG* rot;

	XYZ coords[3];
	LocalFrame cs1;

	RiboseConformerCG();
	RiboseConformerCG(RiboseRotamerCG* rot, LocalFrame& cs1);
	void copyValueFrom(RiboseConformerCG* other);
	void updateLocalFrame(LocalFrame& cs1);
	void updateRotamer(RiboseRotamerCG* rot);
	void updateLocalFrameAndRotamer(LocalFrame& cs1, RiboseRotamerCG* rot);
	double distanceTo(RiboseConformerCG* other);
	virtual ~RiboseConformerCG();
};

} /* namespace NSPmodel */

#endif /* MODEL_RIBOSEROTAMER_H_ */
