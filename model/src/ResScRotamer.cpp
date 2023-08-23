/*
 * ResScRotamer.cpp
 *
 */

#include "model/ResScRotamer.h"

namespace NSPmodel {

ResScRotamer::ResScRotamer(){
	this->aaType = 5;
	this->rotID = 0;
	this->atomNum = 0;

	this->polarAtomNum = 0;

}

ResScRotamer::ResScRotamer(const string& line, AtomLib* atLib) {
	// TODO Auto-generated constructor stub
	ResName rn;
	vector<string> spt;
	splitString(line," ",&spt);
	this->aaType= rn.triToInt(spt[0]);
	this->rotID = 0;

	this->atomNum = (spt.size()-1)/3;


	for(int i=0;i<atomNum;i++){
		this->coordsLocal[i] = XYZ(atof(spt[i*3+1].c_str()), atof(spt[i*3+2].c_str()), atof(spt[i*3+3].c_str()));
	}

	vector<int> scIDs = atLib->getAASidechainUniqueIDs(aaType);
	if(scIDs.size() != atomNum){
		cout << "scID number not equal to atom number: " << scIDs.size() << " " << atomNum << endl;
		cout << line << endl;
	}

	for(int i=0;i<atomNum;i++){
		uniqueIDs[i] = scIDs[i];
	}

	this->rotID = 0;

	this->polarAtomNum = 0;

	LocalFrame cs0;
	LocalFrame cs1;
	LocalFrame cs2;
	LocalFrame cs3;

	if(spt[0] == "ASP"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 2;
		this->polarAtomIndex[1] = 3;
		cs1 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[2]);
		cs2 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[3]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(spt[0] == "GLU"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 3;
		this->polarAtomIndex[1] = 4;
		cs1 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
		cs2 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(spt[0] == "HIS"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 2;
		this->polarAtomIndex[1] = 4;
		cs1 = generateLocalFrameResidueStyle(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
		cs2 = generateLocalFrameResidueStyle(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(spt[0] == "LYS"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 4;
		cs1 = LocalFrame(coordsLocal[2], coordsLocal[3], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(spt[0] == "ASN"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 2;
		this->polarAtomIndex[1] = 3;
		cs1 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[2]);
		cs2 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[3]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(spt[0] == "GLN"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 3;
		this->polarAtomIndex[1] = 4;
		cs1 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
		cs2 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(spt[0] == "ARG"){
		this->polarAtomNum = 3;
		this->polarAtomIndex[0] = 3;
		this->polarAtomIndex[1] = 5;
		this->polarAtomIndex[2] = 6;
		cs1 = generateLocalFrameResidueStyle(coordsLocal[2], coordsLocal[3], coordsLocal[4]);
		cs2 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
		cs3 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[6]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
		this->polarCmList[2] = cs3 - cs0;
	}
	else if(spt[0] == "SER"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 1;
		XYZ t;
		cs1 = LocalFrame(t, coordsLocal[0], coordsLocal[1]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(spt[0] == "THR"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 1;
		XYZ t;
		cs1 = LocalFrame(t, coordsLocal[0], coordsLocal[1]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(spt[0] == "TRP"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 3;
		cs1 = generateLocalFrameResidueStyle(coordsLocal[2], coordsLocal[3], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(spt[0] == "TYR"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 5;
		cs1 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
		this->polarCmList[0] = cs1 - cs0;
	}
}

ResScRotamer::ResScRotamer(Residue* res, AtomLib* atLib){
	if(!res->sidechainComplete(atLib)) {
		cout << "sidechain not complate" << endl;
		exit(1);
	}

	ResName rn;

	this->aaType= res->intName;
	this->rotID = 0;

	this->atomNum = atLib->aaScAtomNames[aaType]->size();

	LocalFrame csRes = res->getCoordSystem();

	for(int i=0;i<atomNum;i++){
		this->coordsLocal[i] = global2local(csRes, res->getAtom(atLib->aaScAtomNames[aaType]->at(i))->coord);
	}

	vector<int> scIDs = atLib->getAASidechainUniqueIDs(aaType);
	if(scIDs.size() != atomNum){
		cout << "scID number not equal to atom number: " << scIDs.size() << " " << atomNum << endl;
		cout << "res: " << res->resID << endl;
	}

	for(int i=0;i<atomNum;i++){
		uniqueIDs[i] = scIDs[i];
	}

	this->rotID = 0;
	string triName = res->triName;
	this->polarAtomNum = 0;

	LocalFrame cs0;
	LocalFrame cs1;
	LocalFrame cs2;
	LocalFrame cs3;

	if(triName == "ASP"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 2;
		this->polarAtomIndex[1] = 3;
		cs1 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[2]);
		cs2 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[3]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(triName == "GLU"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 3;
		this->polarAtomIndex[1] = 4;
		cs1 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
		cs2 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(triName == "HIS"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 2;
		this->polarAtomIndex[1] = 4;
		cs1 = generateLocalFrameResidueStyle(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
		cs2 = generateLocalFrameResidueStyle(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(triName == "LYS"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 4;
		cs1 = LocalFrame(coordsLocal[2], coordsLocal[3], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(triName == "ASN"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 2;
		this->polarAtomIndex[1] = 3;
		cs1 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[2]);
		cs2 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[3]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(triName == "GLN"){
		this->polarAtomNum = 2;
		this->polarAtomIndex[0] = 3;
		this->polarAtomIndex[1] = 4;
		cs1 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
		cs2 = LocalFrame(coordsLocal[1], coordsLocal[2], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
	}
	else if(triName == "ARG"){
		this->polarAtomNum = 3;
		this->polarAtomIndex[0] = 3;
		this->polarAtomIndex[1] = 5;
		this->polarAtomIndex[2] = 6;
		cs1 = generateLocalFrameResidueStyle(coordsLocal[2], coordsLocal[3], coordsLocal[4]);
		cs2 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
		cs3 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[6]);
		this->polarCmList[0] = cs1 - cs0;
		this->polarCmList[1] = cs2 - cs0;
		this->polarCmList[2] = cs3 - cs0;
	}
	else if(triName == "SER"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 1;
		XYZ t;
		cs1 = LocalFrame(t, coordsLocal[0], coordsLocal[1]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(triName == "THR"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 1;
		XYZ t;
		cs1 = LocalFrame(t, coordsLocal[0], coordsLocal[1]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(triName == "TRP"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 3;
		cs1 = generateLocalFrameResidueStyle(coordsLocal[2], coordsLocal[3], coordsLocal[4]);
		this->polarCmList[0] = cs1 - cs0;
	}
	else if(triName == "TYR"){
		this->polarAtomNum = 1;
		this->polarAtomIndex[0] = 5;
		cs1 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
		this->polarCmList[0] = cs1 - cs0;
	}
}

double ResScRotamer::distanceTo(ResScRotamer* other){
	if(this->atomNum != other->atomNum)
		return 9.9;

	if(this->atomNum == 0)
		return 0.0;

	double dd = 0;
	for(int i=0;i<atomNum;i++){
		dd += coordsLocal[i].squaredDistance(other->coordsLocal[i]);
	}
	return sqrt(dd/atomNum);
}

ResScRotamer::~ResScRotamer() {

	// TODO Auto-generated destructor stub
}


ResScConformer::ResScConformer() {
	// TODO Auto-generated constructor stub
	this->aaType = 5;
	this->atomNum = 0;
	this->polarAtomNum = 0;
	this->rot = NULL;
}

void ResScConformer::init(ResScRotamer* rot, const LocalFrame& cs){

	this->aaType = rot->aaType;
	this->atomNum = rot->atomNum;
	this->polarAtomNum = rot->polarAtomNum;

	this->cs = cs;
	this->rot = rot;

	for(int i=0;i<atomNum;i++){
		this->coords[i] = local2global(cs, rot->coordsLocal[i]);
	}

	for(int i=0;i<polarAtomNum;i++){
		this->polarCsList[i] = cs + rot->polarCmList[i];
	}
}

void ResScConformer::updateCs(const LocalFrame& cs){
	this->cs = cs;
	for(int i=0;i<atomNum;i++){
		this->coords[i] = local2global(cs, rot->coordsLocal[i]);
	}

	for(int i=0;i<polarAtomNum;i++){
		this->polarCsList[i] = cs + rot->polarCmList[i];
	}
}

void ResScConformer::updateRotamer(ResScRotamer* rot){


	this->rot = rot;
	this->aaType = rot->aaType;
	this->atomNum = rot->atomNum;
	this->polarAtomNum = rot->polarAtomNum;
	for(int i=0;i<atomNum;i++){
		this->coords[i] = local2global(cs, rot->coordsLocal[i]);
	}

	for(int i=0;i<polarAtomNum;i++){
		this->polarCsList[i] = cs + rot->polarCmList[i];
	}
}

void ResScConformer::copyValueFrom(ResScConformer* other){

	this->aaType = other->aaType;
	this->atomNum = other->atomNum;
	this->polarAtomNum = other->polarAtomNum;
	this->rot = other->rot;
	this->cs = other->cs;

	for(int i=0;i<other->atomNum;i++){
		this->coords[i] = other->coords[i];
	}

	for(int i=0;i<other->polarAtomNum;i++){
		this->polarCsList[i] = other->polarCsList[i];
	}
}

bool ResScConformer::equals(ResScConformer* other){
	if(this->aaType != other->aaType)
		return false;
	if(this->atomNum != other->atomNum)
		return false;
	if(this->rot->rotID != other->rot->rotID)
		return false;
	for(int i=0;i<other->atomNum;i++){
		if(this->coords[i].x_ != other->coords[i].x_)
			return false;
		if(this->coords[i].y_ != other->coords[i].y_)
			return false;
		if(this->coords[i].z_ != other->coords[i].z_)
			return false;
	}
	return true;
}

ResScConformer::~ResScConformer() {

}

}
