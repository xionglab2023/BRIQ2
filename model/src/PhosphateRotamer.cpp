/*
 * PhosphateRotamer.cpp
 *
 *  Created on: 2023Äê8ÔÂ14ÈÕ
 *      Author: nuc
 */

#include <model/PhosphateRotamer.h>

namespace NSPmodel {


PhosphateRotamer::~PhosphateRotamer() {
	// TODO Auto-generated destructor stub
}

PhosphateConformer::PhosphateConformer(){
	this->rot = NULL;
	this->ene = 0.0;
}

PhosphateConformer::PhosphateConformer(PhosphateRotamer* rot, LocalFrame& cs2){
	this->rot = rot;
	this->cs2 = cs2;
	for(int i=0;i<4;i++){
		coords[i] = local2global(cs2, rot->localCoords[i]);
	}
	op1Polar = cs2 + rot->cmOP1;
	op2Polar = cs2 + rot->cmOP2;
	this->ene = 0.0;
}

void PhosphateConformer::copyValueFrom(PhosphateConformer* other){
	this->rot = other->rot;
	this->cs2 = other->cs2;
	for(int i=0;i<4;i++){
		coords[i] = other->coords[i];
	}
	op1Polar = other->op1Polar;
	op2Polar = other->op2Polar;
	this->ene = other->ene;
}

void PhosphateConformer::updateLocalFrameAndRotamer(LocalFrame& cs2, PhosphateRotamer* rot, double ene){
	this->cs2 = cs2;
	this->rot = rot;
	for(int i=0;i<4;i++){
		coords[i] = local2global(cs2, rot->localCoords[i]);
	}
	op1Polar = cs2 + rot->cmOP1;
	op2Polar = cs2 + rot->cmOP2;
	this->ene = ene;
}

void PhosphateConformer::updateLocalFrame(LocalFrame& cs2){
	this->cs2 = cs2;
	for(int i=0;i<4;i++){
		coords[i] = local2global(cs2, rot->localCoords[i]);
	}
	op1Polar = cs2 + rot->cmOP1;
	op2Polar = cs2 + rot->cmOP2;
}

void PhosphateConformer::updateRotamer(PhosphateRotamer* rot){
	this->rot = rot;
	for(int i=0;i<4;i++){
		coords[i] = local2global(cs2, rot->localCoords[i]);
	}
	op1Polar = cs2 + rot->cmOP1;
	op2Polar = cs2 + rot->cmOP2;
}

double PhosphateConformer::distanceTo(PhosphateConformer* other){
	double dd = 0.0;

	for(int i=0;i<4;i++){
		dd += squareDistance(coords[i], other->coords[i]);
	}
	return sqrt(dd/4);
}

PhosphateConformer::~PhosphateConformer(){

}

} /* namespace NSPmodel */
