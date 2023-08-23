/*
 * ResConformer.cpp
 *
 *  Created on: 2022Äê3ÔÂ3ÈÕ
 *      Author: pengx
 */

#include "model/ResConformer.h"

namespace NSPmodel {

ResConformer::ResConformer() {
	this->aaType = 0;
	this->bbConf = new ResBBConformer();
	this->scConf = new ResScConformer();
}

ResConformer::ResConformer(int aaType) {
	this->aaType = aaType;
	this->bbConf = new ResBBConformer();
	this->scConf = new ResScConformer();
}

void ResConformer::init(ResScRotamer* rotSc, ResBBRotamer* rotBB, const LocalFrame& cs2){
	this->aaType = rotSc->aaType;
	this->bbConf->init(rotBB, cs2);
	this->scConf->init(rotSc, cs2);
}

void ResConformer::copyValueFrom(ResConformer* other){
	this->aaType = other->aaType;
	this->scConf->copyValueFrom(other->scConf);
	this->bbConf->copyValueFrom(other->bbConf);
}

void ResConformer::updateScRotamer(ResScRotamer* rotSc){
	this->aaType = rotSc->aaType;
	this->scConf->updateRotamer(rotSc);
}

ResConformer::~ResConformer() {
	delete bbConf;
	delete scConf;
}

} /* namespace NSPmodel */
