/*
 * ResBBRotamer.cpp
 *
 */

#include "model/ResBBRotamer.h"

namespace NSPmodel {

ResBBRotamer::ResBBRotamer(Residue* resP, Residue* res, AtomLib* atLib){
	this->aaType = res->intName;

	Atom* preC = resP->getAtom("C");
	Atom* N = res->getAtom("N");
	Atom* CA = res->getAtom("CA");
	Atom* C = res->getAtom("C");
	Atom* O = res->getAtom("O");

	if(preC==NULL || N == NULL || CA == NULL || C == NULL || O == NULL)
	{
		cout << "backbone not complete: res: " << res->resID << endl;
		exit(0);
	}

	LocalFrame csX(res->getAtom("C")->getCoord(), res->getAtom("CA")->getCoord(), res->getAtom("N")->getCoord());


	XYZ tpc = preC->coord;
	double d = preC->coord.distance(N->coord);
	if(d > 1.50)
	{
		tpc = csX.coordNext(1.32, 122.01, -139.02);
	}

	LocalFrame cs = res->getCoordSystem();

	localTList[0] = global2local(cs, tpc);
	localTList[1] = global2local(cs, N->coord);
	localTList[2] = global2local(cs, CA->coord);
	localTList[3] = global2local(cs, C->coord);
	localTList[4] = global2local(cs, O->coord);

	LocalFrame cs1 = generateLocalFrameResidueStyle(localTList[0], localTList[1], localTList[2]);
	LocalFrame cs2;
	LocalFrame cs3 = generateLocalFrameResidueStyle(localTList[2], localTList[3], localTList[4]);
	LocalFrame cs4 = LocalFrame(localTList[2], localTList[3], localTList[4]); //O

	dm = DistanceMatrixHbond(cs1, cs3);

	cm21 = cs1 - cs2;
	cm23 = cs3 - cs2;
	cm12 = cs2 - cs1;
	cm32 = cs2 - cs3;
	cm13 = cs3 - cs1;
	cm31 = cs1 - cs3;
	cm24 = cs4 - cs2;

	this->phi = dihedral(localTList[0], localTList[1], localTList[2], localTList[3]);
	this->psi = dihedral(localTList[1], localTList[2], localTList[3], localTList[4]);
	psi = psi + 180;
	if(psi > 180)
		psi = psi - 360;

	this->ene = 0.0;
	this->index1K = 0;
	this->index1W = 0;

	for(int i=0;i<4;i++)
	{
		this->uniqueID[i] = atLib->aaUniqueIDs[aaType][i];
	}
}

ResBBRotamer::ResBBRotamer(Residue* res, AtomLib* atLib){
	this->aaType = res->intName;

	Atom* N = res->getAtom("N");
	Atom* CA = res->getAtom("CA");
	Atom* C = res->getAtom("C");
	Atom* O = res->getAtom("O");

	if(N == NULL || CA == NULL || C == NULL || O == NULL)
	{
		cout << "backbone not complete: res: " << res->resID << endl;
		exit(0);
	}

	LocalFrame csX(res->getAtom("C")->getCoord(), res->getAtom("CA")->getCoord(), res->getAtom("N")->getCoord());
	XYZ tpc = csX.coordNext(1.32, 122.01, -139.02);
	LocalFrame cs = res->getCoordSystem();

	localTList[0] = global2local(cs, tpc);
	localTList[1] = global2local(cs, N->coord);
	localTList[2] = global2local(cs, CA->coord);
	localTList[3] = global2local(cs, C->coord);
	localTList[4] = global2local(cs, O->coord);

	LocalFrame cs1 = generateLocalFrameResidueStyle(localTList[0], localTList[1], localTList[2]);
	LocalFrame cs2;
	LocalFrame cs3 = generateLocalFrameResidueStyle(localTList[2], localTList[3], localTList[4]);
	LocalFrame cs4 = LocalFrame(localTList[2], localTList[3], localTList[4]); //O

	dm = DistanceMatrixHbond(cs1, cs3);

	cm21 = cs1 - cs2;
	cm23 = cs3 - cs2;
	cm12 = cs2 - cs1;
	cm32 = cs2 - cs3;
	cm13 = cs3 - cs1;
	cm31 = cs1 - cs3;
	cm24 = cs4 - cs2;

	this->phi = dihedral(localTList[0], localTList[1], localTList[2], localTList[3]);
	this->psi = dihedral(localTList[1], localTList[2], localTList[3], localTList[4]);
	psi = psi + 180;
	if(psi > 180)
		psi = psi - 360;

	this->ene = 0.0;
	this->index1K = 0;
	this->index1W = 0;

	for(int i=0;i<4;i++)
	{
		this->uniqueID[i] = atLib->aaUniqueIDs[aaType][i];
	}
}

ResBBRotamer::ResBBRotamer(const string& line, AtomLib* atLib) {

	ResName rn;
	vector<string> spt;
	splitString(line," ",&spt);
	this->aaType= rn.triToInt(spt[0]);
	this->localTList[0] = XYZ(atof(spt[1].c_str()), atof(spt[2].c_str()), atof(spt[3].c_str()));
	this->localTList[1] = XYZ(atof(spt[4].c_str()), atof(spt[5].c_str()), atof(spt[6].c_str()));
	this->localTList[2] = XYZ();
	this->localTList[3] = XYZ(atof(spt[7].c_str()), atof(spt[8].c_str()), atof(spt[9].c_str()));
	this->localTList[4] = XYZ(atof(spt[10].c_str()), atof(spt[11].c_str()), atof(spt[12].c_str()));

	if(spt.size() > 13)
		this->ene = atof(spt[13].c_str());
	else
		this->ene = 0.0;

	if(spt.size() > 15){
		this->index1K = atoi(spt[14].c_str());
		this->index1W = atoi(spt[15].c_str());
	}

	LocalFrame cs1 = generateLocalFrameResidueStyle(localTList[0], localTList[1], localTList[2]); //N
	LocalFrame cs2; //CA
	LocalFrame cs3 = generateLocalFrameResidueStyle(localTList[2], localTList[3], localTList[4]); //C
	LocalFrame cs4 = LocalFrame(localTList[2], localTList[3], localTList[4]); //O

	dm = DistanceMatrixHbond(cs1, cs3);

	cm21 = cs1 - cs2;
	cm23 = cs3 - cs2;
	cm12 = cs2 - cs1;
	cm32 = cs2 - cs3;
	cm13 = cs3 - cs1;
	cm31 = cs1 - cs3;
	cm24 = cs4 - cs2;

	this->phi = dihedral(localTList[0], localTList[1], localTList[2], localTList[3]);
	this->psi = dihedral(localTList[1], localTList[2], localTList[3], localTList[4]);
	psi = psi + 180;
	if(psi > 180)
		psi = psi - 360;

	for(int i=0;i<4;i++)
	{
		this->uniqueID[i] = atLib->aaUniqueIDs[aaType][i];
	}
}

void ResBBRotamer::setAAType(int aa, AtomLib* atLib){
	this->aaType = aa;
	for(int i=0;i<4;i++)
	{
		this->uniqueID[i] = atLib->aaUniqueIDs[aaType][i];
	}
}

string ResBBRotamer::toString(){
	cout << "tostring" << endl;
	char xx[200];
	ResName rn;
	sprintf(xx, "%s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f",rn.intToTri(this->aaType).c_str(), localTList[0].x_,localTList[0].y_,localTList[0].z_,localTList[1].x_,localTList[1].y_,localTList[1].z_, localTList[3].x_,localTList[3].y_,localTList[3].z_,localTList[4].x_,localTList[4].y_,localTList[4].z_);
	return string(xx);

}

void ResBBConformer::init(ResBBRotamer* rot, const LocalFrame& cs2){
	this->cs2 = cs2;
	this->cs1 = cs2 + rot->cm21;
	this->cs3 = cs2 + rot->cm23;
	this->cs4 = cs2 + rot->cm24;

	this->aaType = rot->aaType;
	this->rot = rot;
	for(int i=0;i<5;i++){
		this->coords[i] = local2global(cs2, rot->localTList[i]);
	}
}

void ResBBConformer::updateRotFixCs1(ResBBRotamer* rot){
	this->aaType = rot->aaType;
	this->rot = rot;

	this->cs2 = this->cs1 + rot->cm12;
	this->cs3 = this->cs1 + rot->cm13;
	this->cs4 = this->cs2 + rot->cm24;

	for(int i=0;i<5;i++){
		this->coords[i] = local2global(cs2, rot->localTList[i]);
	}
}

void ResBBConformer::updateRotFixCs2(ResBBRotamer* rot){
	this->aaType = rot->aaType;
	this->rot = rot;

	this->cs1 = this->cs2 + rot->cm21;
	this->cs3 = this->cs2 + rot->cm23;
	this->cs4 = this->cs2 + rot->cm24;

	for(int i=0;i<5;i++){
		this->coords[i] = local2global(cs2, rot->localTList[i]);
	}
}

void ResBBConformer::updateRotFixCs3(ResBBRotamer* rot){
	this->aaType = rot->aaType;
	this->rot = rot;

	this->cs1 = this->cs3 + rot->cm31;
	this->cs2 = this->cs3 + rot->cm32;
	this->cs4 = this->cs2 + rot->cm24;

	for(int i=0;i<5;i++){
		this->coords[i] = local2global(cs2, rot->localTList[i]);
	}
}

void ResBBConformer::updateCs1(LocalFrame& cs1){
	this->cs1 = cs1;
	this->cs2 = this->cs1 + rot->cm12;
	this->cs3 = this->cs1 + rot->cm13;
	this->cs4 = this->cs2 + rot->cm24;

	for(int i=0;i<5;i++){
		this->coords[i] = local2global(cs2, rot->localTList[i]);
	}
}

void ResBBConformer::updateCs2(LocalFrame& cs2){
	this->cs2 = cs2;
	this->cs1= this->cs2 + rot->cm21;
	this->cs3 = this->cs2 + rot->cm23;
	this->cs4 = this->cs2 + rot->cm24;

	for(int i=0;i<5;i++){
		this->coords[i] = local2global(cs2, rot->localTList[i]);
	}
}

void ResBBConformer::updateCs3(LocalFrame& cs3){
	this->cs3 = cs3;
	this->cs1 = this->cs3 + rot->cm31;
	this->cs2 = this->cs3 + rot->cm32;
	this->cs4 = this->cs2 + rot->cm24;

	for(int i=0;i<5;i++){
		this->coords[i] = local2global(cs2, rot->localTList[i]);
	}
}

void ResBBConformer::copyValueFrom(ResBBConformer* other){
	for(int i=0;i<5;i++){
		coords[i] = other->coords[i];
	}

	cs1 = other->cs1;
	cs2 = other->cs2;
	cs3 = other->cs3;
	cs4 = other->cs4;
	this->rot = other->rot;
	this->aaType = other->aaType;
}

ResBBConformer::~ResBBConformer() {
	// TODO Auto-generated destructor stub
}

ResBBRotamer::~ResBBRotamer() {
	// TODO Auto-generated destructor stub
}

}
