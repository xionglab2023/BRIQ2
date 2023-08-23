/*
 * RiboseRotamer.cpp
 *
 *  Created on: 2022Äê9ÔÂ6ÈÕ
 *      Author: pengx
 */

#include "model/RiboseRotamer.h"

namespace NSPmodel {

RiboseRotamer::RiboseRotamer(RNABase* base) {
	vector<string> nameList;

	if(base->baseTypeInt < 4) {
		nameList.push_back("C1'"); //0
		nameList.push_back("C2'"); //1
		nameList.push_back("C3'"); //2
		nameList.push_back("C4'"); //3
		nameList.push_back("O4'"); //4
		nameList.push_back("O3'"); //5
		nameList.push_back("C5'"); //6
		nameList.push_back("O2'"); //7
		this->atomNum = 8;
	}
	else {
		nameList.push_back("C1'");
		nameList.push_back("C2'");
		nameList.push_back("C3'");
		nameList.push_back("C4'");
		nameList.push_back("O4'");
		nameList.push_back("O3'");
		nameList.push_back("C5'");
		this->atomNum = 7;
	}

	LocalFrame cs = base->getCoordSystem();
	for(int i=0;i<atomNum;i++){
		Atom* a = base->getAtom(nameList[i]);
		if(a == NULL) {
			cout << "base backbone incomplete: " << base->baseID << " " << base->baseType << endl;
			exit(1);
		}
		localCoords[i] = cs.global2localcrd(a->coord);
		//cout << "base rotamer: " <<  tList1[i].toString() << endl;
	}

	double x=0,y=0,z=0;
	for(int i=0;i<atomNum;i++){
		x += localCoords[i].x_;
		y += localCoords[i].y_;
		z += localCoords[i].z_;
	}
	center = XYZ(x*0.125, y*0.125, z*0.125);
	this->energy = 0.0;

	XYZ N = XYZ(1.475, 0, 0);
	LocalFrame cs1;
	LocalFrame cs2 = LocalFrame(localCoords[1], localCoords[2], localCoords[5]); //C2'-C3'-O3'
	LocalFrame cs3 = LocalFrame(localCoords[4], localCoords[3], localCoords[6]); //O4'-C4'-C5'

	LocalFrame csX = LocalFrame(N, localCoords[0], localCoords[3]);
	mv12 = cs2 - cs1;
	mv13 = cs3 - cs1;
	mv31 = cs1 - cs3;
	mv21 = cs1 - cs2;

	if(atomNum == 8){
		LocalFrame csO2 = LocalFrame(localCoords[0], localCoords[1], localCoords[7]);
		mv1O2 = csO2 - cs1;
	}


	mv1X = csX - cs1;
	mvX1 = cs1 - csX;

	this->rotTypeLv1 = 0;
	this->rotType = 0;
	this->resType = base->baseTypeInt;
	double ang = dihedral(localCoords[0], localCoords[1], localCoords[3], localCoords[2]);
	if(ang > 0)
		ang = 180 - ang;
	else
		ang = -ang - 180;
	this->improper = ang;
	XYZ a  = localCoords[4];
	XYZ b = localCoords[0];
	XYZ c(1.48, 0, 0);
	XYZ d(2.318, -1.073, 0);
	this->chi = dihedral(a,b,c,d);

}

RiboseRotamer::RiboseRotamer(const string& line) {

	map<char,int> typeToInt;
	typeToInt['A'] = 0;
	typeToInt['U'] = 1;
	typeToInt['G'] = 2;
	typeToInt['C'] = 3;
	typeToInt['a'] = 4;
	typeToInt['t'] = 5;
	typeToInt['g'] = 6;
	typeToInt['c'] = 7;

	this->resType = typeToInt[line[0]];
	if(resType < 4)
		atomNum = 8;
	else
		atomNum = 7;

	vector<string> spt;
	splitString(line," ",&spt);

	localCoords[0] = XYZ();

	for(int i=1;i<atomNum;i++){
		localCoords[i] = XYZ(atof(spt[i*3-2].c_str()), atof(spt[i*3-1].c_str()), atof(spt[i*3].c_str()));
	}

	double x=0,y=0,z=0;
	for(int i=0;i<atomNum;i++){
		x += localCoords[i].x_;
		y += localCoords[i].y_;
		z += localCoords[i].z_;
	}
	center = XYZ(x*0.125, y*0.125, z*0.125);
	this->energy = atof(spt[3*atomNum-2].c_str());

	XYZ N = XYZ(1.475, 0, 0);
	LocalFrame cs1;
	LocalFrame cs2 = LocalFrame(localCoords[1], localCoords[2], localCoords[5]);
	LocalFrame cs3 = LocalFrame(localCoords[4], localCoords[3], localCoords[6]);

	LocalFrame csX = LocalFrame(N, localCoords[0], localCoords[3]);
	mv12 = cs2 - cs1;
	mv13 = cs3 - cs1;
	mv31 = cs1 - cs3;
	mv21 = cs1 - cs2;

	if(atomNum == 8){
		LocalFrame csO2 = LocalFrame(localCoords[0], localCoords[1], localCoords[7]);
		mv1O2 = csO2 - cs1;
	}


	mv1X = csX - cs1;
	mvX1 = cs1 - csX;

	this->rotTypeLv1 = 0;
	this->rotType = 0;
	double ang = dihedral(localCoords[0], localCoords[1], localCoords[3], localCoords[2]);
	if(ang > 0)
		ang = 180 - ang;
	else
		ang = -ang - 180;
	this->improper = ang;
	XYZ a  = localCoords[4];
	XYZ b = localCoords[0];
	XYZ c(1.48, 0, 0);
	XYZ d(2.318, -1.073, 0);
	this->chi = dihedral(a,b,c,d);
}

RiboseConformer::RiboseConformer(){
	this->hasO2 = false;
	this->rot = NULL;
}

RiboseConformer::RiboseConformer(RiboseRotamer* rot, LocalFrame& cs1){
	this->rot = rot;
	this->cs1 = cs1;
	this->cs2 = cs1 + rot->mv12;
	this->cs3 = cs1 + rot->mv13;
	if(rot->resType < 4)
		hasO2 = true;
	else
		hasO2 = false;

	if(hasO2){
		this->o2Polar = cs1 + rot->mv1O2;
	}
}

void RiboseConformer::copyValueFrom(RiboseConformer* other){
	this->rot = other->rot;
	for(int i=0;i<8;i++){
		this->coords[i] = other->coords[i];
	}
	cs1 = other->cs1;
	cs2 = other->cs2;
	cs3 = other->cs3;
	hasO2 = other->hasO2;
	o2Polar = other->o2Polar;

}

void RiboseConformer::updateLocalFrame(LocalFrame& cs1){
	this->cs1 = cs1;
	this->cs2 = cs1 + rot->mv12;
	this->cs3 = cs1 + rot->mv13;

	for(int i=0;i<rot->atomNum;i++){
		coords[i] = local2global(cs1, rot->localCoords[i]);
	}

	if(hasO2){
		this->o2Polar = cs1 + rot->mv1O2;
	}
}

void RiboseConformer::updateRotamer(RiboseRotamer* rot){
	this->rot = rot;
	//cs1 fixed
	this->cs2 = cs1 + rot->mv12;
	this->cs3 = cs1 + rot->mv13;

	for(int i=0;i<rot->atomNum;i++){
		coords[i] = local2global(cs1, rot->localCoords[i]);
	}

	if(rot->resType < 4)
		hasO2 = true;
	else
		hasO2 = false;

	if(hasO2){
		this->o2Polar = cs1 + rot->mv1O2;
	}
}

void RiboseConformer::updateRotamerCs2Fixed(RiboseRotamer* rot){
	this->rot = rot;
	//cs2 fixed
	this->cs1 = cs2 + rot->mv21;
	this->cs3 = cs1 + rot->mv13;

	for(int i=0;i<rot->atomNum;i++){
		coords[i] = local2global(cs1, rot->localCoords[i]);
	}

	if(rot->resType < 4)
		hasO2 = true;
	else
		hasO2 = false;

	if(hasO2){
		this->o2Polar = cs1 + rot->mv1O2;
	}
}

void RiboseConformer::updateRotamerCs3Fixed(RiboseRotamer* rot){
	this->rot = rot;
	//cs3 fixed
	this->cs1 = cs3 + rot->mv31;
	this->cs2 = cs1 + rot->mv12;

	for(int i=0;i<rot->atomNum;i++){
		coords[i] = local2global(cs1, rot->localCoords[i]);
	}

	if(rot->resType < 4)
		hasO2 = true;
	else
		hasO2 = false;

	if(hasO2){
		this->o2Polar = cs1 + rot->mv1O2;
	}
}

void RiboseConformer::updateLocalFrameAndRotamer(LocalFrame& cs1, RiboseRotamer* rot){
	this->rot = rot;
	this->cs1 = cs1;
	this->cs2 = cs1 + rot->mv12;
	this->cs3 = cs1 + rot->mv13;

	for(int i=0;i<rot->atomNum;i++){
		coords[i] = local2global(cs1, rot->localCoords[i]);
	}

	if(rot->resType < 4)
		hasO2 = true;
	else
		hasO2 = false;

	if(hasO2){
		this->o2Polar = cs1 + rot->mv1O2;
	}
}



RiboseRotamer::~RiboseRotamer() {
	// TODO Auto-generated destructor stub
}

RiboseConformer::~RiboseConformer(){

}

} /* namespace NSPmodel */
