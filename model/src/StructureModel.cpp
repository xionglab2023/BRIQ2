/*
 * StructureModel.cpp
 *
 */

#include <array>
#include "model/StructureModel.h"

namespace NSPmodel {

Atom::Atom(){
	this->name = "X";
	this->type = "X";
	this->resType = "";
	this->alt = ' ';
	this->coord = XYZ();
}

Atom::Atom(const string& line){
	string s = line.substr(12,4);
	this->name = trimString(s);
	if(this->name == "OXT" || this->name == "OXT1")
		this->name = "O";
	this->alt = line.at(16);
	this->type = "";
	float x = atof(line.substr(30,8).c_str());
	float y = atof(line.substr(38,8).c_str());
	float z = atof(line.substr(46,8).c_str());
	this->coord = XYZ(x,y,z);
	guessAtomType();
	s = line.substr(17,3);
	this->resType = trimString(s);
}

Atom::Atom(const string& line, int fileType){
	if(fileType == 0){
		//pdb file
		string s = line.substr(12,4);
		this->name = trimString(s);
		if(this->name == "OXT" || this->name == "OXT1")
			this->name = "O";
		this->alt = line.at(16);
		this->type = "";
		float x = atof(line.substr(30,8).c_str());
		float y = atof(line.substr(38,8).c_str());
		float z = atof(line.substr(46,8).c_str());
		this->coord = XYZ(x,y,z);
		guessAtomType();
		s = line.substr(17,3);
		this->resType = trimString(s);
	}
	else {
		//cif file
		vector<string> spt;
		splitString(line, " ", &spt);
		this->name = spt[3];
		if(this->name[0] == '"')
			this->name = name.substr(1,this->name.length()-2);

		//cout << name << " " << name.length() << endl;
		
		if(this->name == "OXT" || this->name == "OXT1")
			this->name = "O";
		this->alt = spt[9][0];
		if(alt == '?') this->alt = ' ';
		this->type = spt[2];
		float x = atof(spt[10].c_str());
		float y = atof(spt[11].c_str());
		float z = atof(spt[12].c_str());
		this->coord = XYZ(x,y,z);
		guessAtomType();
		this->resType = spt[5];
	}
}

Atom::Atom(string name, const XYZ& coord){
	this->name = name;
	this->coord = coord;
	this->type = "X";
	this->resType = "UNK";
	this->alt = ' ';
}

Atom::Atom(const Atom& other) {
	this->name = other.name;
	this->coord = other.coord;
	this->type = other.type;
	this->resType = other.resType;
	this->alt = other.alt;
}

Atom& Atom::operator=(const Atom& other)
{
	if(this == &other)
		return *this;
	this->name = other.name;
	this->type = other.type;
	this->coord = other.coord;
	this->resType = other.resType;
	return *this;
}

bool Atom::isBackboneAtom() const
{
	return name=="N" || name=="CA" || name=="C" || name=="O" || name=="OXT" || name=="OT1" || name=="OT2";
}

void Atom::guessAtomType()
{
	char c = name.at(0);
	if(name.length() == 1)
		this->type = name;
	else if(name == "FE" || name == "ZN" || name=="MG" || name=="MN" || name=="CL" || name == "SE" || name=="CU" || name=="NA")
		this->type = name;
	else if(c>='A' && c <= 'Z')
	{
	    this->type = string(1,c);

	}
	else
		this->type = string(1,name.at(1));
}

bool Atom::isIon() const
{
	if(type == "FE" || type == "ZN" || type == "MG" || type == "MN" || type == "CA" || type == "CU" || type == "NI")
		return true;
	return false;
}

string& Atom::getName()
{
	return this->name;
}

string& Atom::getType()
{
	return this->type;
}

string& Atom::getResType()
{
	return this->resType;
}



XYZ& Atom::getCoord()
{
	return this->coord;
}

void Atom::setCoord(const XYZ& coord)
{
	this->coord = coord;
}

void Atom::setCoord(const XYZ&& coord)
{
	this->coord = coord;
}

float Atom::distance(const Atom& other) const
{
	return this->coord.distance(other.coord);
}

void Atom::setAtomType(string type)
{
	this->type = type;
}

void Atom::setResType(string resType)
{
	this->resType = resType;
}

string Atom::nameString() const
{
	char s[10];
	if(this->type.length() > 1)
		sprintf(s,"%-5s",this->name.c_str());
	else
		sprintf(s," %-4s",this->name.c_str());
	string ss = string(s);
	return ss;
}

Atom::~Atom(){

}

Residue::Residue() {
	this->resID = "1";
	this->resSeqID = 0;
	this->atomNum = 0;
	this->triName = "UNK";
	this->intName = 20;
	this->chainID = "-";
	this->hasLocalFrame = false;
	this->coordSys;
	this->altLoc = ' ';
	this->hasAltConf = false;

}

Residue::Residue(string resID, string chainID, string triName)
{
	ResName rn;
	this->resID = resID;
	this->chainID = chainID;
	this->triName = triName;
	this->intName = rn.triToInt(triName);
	this->resSeqID = 0;
	this->atomNum = 0;
	this->hasLocalFrame = false;
	this->altLoc = ' ';
	this->hasAltConf = false;
}

void Residue::addAtom(Atom* a)
{
	this->atomList.push_back(a);
	if(a->isBackboneAtom())
		this->backboneAtoms.push_back(a);
	else
		this->sidechainAtoms.push_back(a);
	this->atomMap[a->getName()] = a;
	this->atomNum ++;
}

void Residue::setResSeqID(int id)
{
    this->resSeqID = id;
}

bool Residue::hasThreeCoreAtoms() const
{
	if(atomMap.find("N") == atomMap.end())
		return false;
	if(atomMap.find("CA") == atomMap.end())
		return false;
	if(atomMap.find("C") == atomMap.end())
		return false;
	return true;
}

void Residue::updateCoordSystem()
{
	if(!this->hasThreeCoreAtoms())
		return;
	Atom* N = this->atomMap.at("N");
	Atom* C = this->atomMap.at("C");
	Atom* CA = this->atomMap.at("CA");

	XYZ n = N->getCoord();
	XYZ c = C->getCoord();
	XYZ ca = CA->getCoord();

	XYZ can = ~(N->getCoord() - CA->getCoord());
	XYZ cac = ~(C->getCoord() - CA->getCoord());

	XYZ z = ~(can^cac);
	XYZ x = ~(can+cac);

	XYZ y = ~(z^x);

	this->coordSys = LocalFrame(ca, x, y, z);

}

bool Residue::hasAtom(const string& atomName) const
{
    if(atomMap.find(atomName) == atomMap.end())
		return false;
	else
        return true;
}

Atom* Residue::getAtom(const string& atomName)
{
    map<string,Atom*>::const_iterator it = atomMap.find(atomName);
    if(it == atomMap.end())
        return NULL;
    else
        return it->second;
}

vector<Atom*>* Residue::getAtomList()
{
    return &this->atomList;
}

vector<Atom*>* Residue::getBackboneAtoms()
{
    return &this->backboneAtoms;
}

vector<Atom*>* Residue::getSidechainAtoms()
{
    return &this->sidechainAtoms;
}

bool Residue::contactToNeighbor(Residue* other){
	Atom* c = getAtom("C");
	Atom* n = other->getAtom("N");
	if(c == NULL || n == NULL) return false;
	if(c->coord.distance(n->coord) < 1.6)
		return true;
	else
		return false;
}

bool Residue::sidechainComplete(AtomLib* atomLib) const
{
	vector<string> scAtoms;
	atomLib->getAminoAcidSidechainAtomNames(this->intName, scAtoms);
	if(scAtoms.size() == 0) return false;
	vector<string>::iterator it;
	for(it=scAtoms.begin();it<scAtoms.end();it++)
	{
		string s = *it;
		if(!hasAtom(s))
			return false;
	}
	return true;
}

XYZ Residue::getCbCoord()
{
	if(!hasThreeCoreAtoms())
	{
		cerr << "lack backbone atom information" << endl;
		exit(1);
	}
	if(!this->hasLocalFrame)
		updateCoordSystem();
	XYZ localCB = XYZ(-0.942, 0.009, 1.208);
	return this->coordSys.local2globalcrd(localCB);
}

LocalFrame& Residue::getCoordSystem()
{
	if(!hasThreeCoreAtoms())
	{
		cerr << "lack backbone atom information" << endl;
		exit(1);
	}
	if(!this->hasLocalFrame)
		updateCoordSystem();
	return this->coordSys;
}

string Residue::getChainID() const
{
	return this->chainID;
}

int Residue::getResSeqID() const
{
	return this->resSeqID;
}

string Residue::getResID() const
{
	return this->resID;
}

string Residue::getType() const
{
	return this->triName;
}

int Residue::printPDBFormat(ofstream& out, int startAtomID) const
{
    char c = this->resID.at(resID.length()-1);
    char s[100];
    vector<Atom*>::const_iterator it;
    int atomID = startAtomID;
    if(c >= '0' && c <= '9')
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s%3s %c%4s    %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->triName.c_str(),this->chainID,this->resID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    else
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s%3s %c%5s   %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->triName.c_str(),this->chainID,this->resID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    return atomID;
}

Residue::~Residue() {
}

RNABase::RNABase() {
	this->baseID = "1";
	this->chainID = 'A';
	this->baseType = 'N';
	this->baseTypeInt = 4;
	this->baseSeqID = -1;
	this->hasLocalFrame = false;
	this->hasAltConf = false;
	this->altLoc = '-';
}

RNABase::RNABase(const string& baseID, const string& chainID, char baseType){
	this->baseID = baseID;
	this->chainID = chainID;
	this->baseType = baseType;
	if(baseType == 'A')
		this->baseTypeInt = 0;
	else if(baseType == 'U')
		this->baseTypeInt = 1;
	else if(baseType == 'G')
		this->baseTypeInt = 2;
	else if(baseType == 'C')
		this->baseTypeInt = 3;
	else if(baseType == 'a')
		this->baseTypeInt = 4;
	else if(baseType == 't')
		this->baseTypeInt = 5;
	else if(baseType == 'g')
		this->baseTypeInt = 6;
	else if(baseType == 'c')
		this->baseTypeInt = 7;
	else
		this->baseTypeInt = -1;

	this->baseSeqID = -1;
	this->hasLocalFrame = false;
	this->hasAltConf = false;
	this->altLoc = '-';
}

RNABase::RNABase(const RNABase& other) {
	this->baseID = other.baseID;
	this->chainID = other.chainID;
	this->baseType = other.baseType;
	this->baseTypeInt = other.baseTypeInt;
	this->baseSeqID = other.baseSeqID;
	this->hasLocalFrame = other.hasLocalFrame;
	this->hasAltConf = other.hasAltConf;
	this->altLoc = other.altLoc;

	int la = other.atomList.size();
	map<Atom*, Atom*> old2newAtMap;
	for(int i=0; i<la; i++) {
		this->atomList.emplace_back(new Atom{*other.atomList[i]});
		old2newAtMap.emplace(other.atomList[i], this->atomList[i]);
	}

	la = other.backboneAtoms.size();
	for(int i=0; i<la; i++) {
		this->backboneAtoms.emplace_back(old2newAtMap.at(other.backboneAtoms[i]));
	}

	la = other.sidechainAtoms.size();
	for(int i=0; i<la; i++) {
		this->sidechainAtoms.emplace_back(old2newAtMap.at(other.sidechainAtoms[i]));
	}

	for(auto & iter : other.atomMap) {
		this->atomMap.emplace(iter.first, old2newAtMap.at(iter.second));
	}
}

RNABase::RNABase(RNABase&& other) noexcept {  // noexcept 意味着不会调用 new，delete 和 IO
	this->baseID = move(other.baseID);  // string 具有移动构造函数
	this->chainID = move(other.chainID);  // string
	this->baseType = other.baseType;
	this->baseTypeInt = other.baseTypeInt;
	this->baseSeqID = other.baseSeqID;
	this->hasLocalFrame = other.hasLocalFrame;
	this->hasAltConf = other.hasAltConf;
	this->altLoc = other.altLoc;
	this->coordSys = move(other.coordSys);  // coordSys 会调用拷贝构造?
	this->atomList = move(other.atomList);  // vector  具有移动构造函数
	this->backboneAtoms = move(other.backboneAtoms);  // vector 移动构造
	this->sidechainAtoms = move(other.sidechainAtoms);  // vector 移动构造
	this->atomMap = move(other.atomMap);  // map 移动构造
}

void RNABase::addAtom(Atom* a) {
	this->atomList.push_back(a);
	string atomName = a->name;
	if(atomName == "P" || atomName == "OP1" || atomName == "OP2") {
		this->backboneAtoms.push_back(a);
	}
	else if(atomName.length() == 3 && atomName[2] == '\'')
		this->backboneAtoms.push_back(a);
	else
		this->sidechainAtoms.push_back(a);
	atomMap[a->name] = a;
}

void RNABase::updateCoordSystem() {
	if(hasLocalFrame) return;
	if(this->baseType == 'A' || this->baseType == 'G' || this->baseType == 'a' || this->baseType == 'g') {
		if(hasAtom("C1'") && hasAtom("N9") && hasAtom("C4")) {
			XYZ c1 = this->atomMap["C1'"]->coord;
			XYZ n9 = this->atomMap["N9"]->coord;
			XYZ c4 = this->atomMap["C4"]->coord;

			XYZ x = ~(n9-c1);
			XYZ z = ~((n9-c4)^x);
			XYZ y = z^x;
			this->coordSys = LocalFrame(c1, x, y, z);
			this->hasLocalFrame = true;
		}
	}
	else if(this->baseType == 'U' || this->baseType == 'C' || this->baseType == 't' || this->baseType == 'c') {
		if(hasAtom("C1'") && hasAtom("N1") && hasAtom("C2")) {
			XYZ c1 = this->atomMap["C1'"]->coord;
			XYZ n1 = this->atomMap["N1"]->coord;
			XYZ c2 = this->atomMap["C2"]->coord;

			XYZ x = ~(n1-c1);
			XYZ z = ~((n1-c2)^x);
			XYZ y = z^x;
			this->coordSys = LocalFrame(c1, x, y, z);
			this->hasLocalFrame = true;
		}
	}
}

bool RNABase::sidechainComplete(AtomLib* atLib) const{
	vector<string> atomNames;
	atLib->getRnaSidechainAtoms(this->baseTypeInt, atomNames);
	for(unsigned int i=0;i<atomNames.size();i++) {
		if(this->atomMap.find(atomNames.at(i)) == this->atomMap.end())
			return false;
	}
	return true;
}

array<XYZ,4> RNABase::getFourPseudoAtomCoords(){
	array<XYZ,4> list;
	if(!hasLocalFrame)
		updateCoordSystem();
	if(!hasLocalFrame) {
		cerr << "base atom not complete, can't generate coordinate system!" << endl;
		return list;
	}


	XYZ a(2.158 ,  3.826 ,  1.427);
	XYZ b(-0.789 , -0.329 , -1.273);
	XYZ c(4.520 , -3.006 ,  1.586);
	XYZ d(6.018 ,  1.903 , -1.638);
	get<0>(list) = local2global(coordSys, a);
	get<1>(list) = local2global(coordSys, b);
	get<2>(list) = local2global(coordSys, c);
	get<3>(list) = local2global(coordSys, d);
	return list;
}

ConvexPolygon RNABase::getBaseConvexPolygon(LocalFrame& cs, AtomLib* atLib){
	if(!sidechainComplete(atLib)) return ConvexPolygon();
	vector<XY> points;
	if(baseTypeInt == 0 || baseTypeInt == 4){

				XYZ a = global2local(cs, getAtom("N9")->coord);
				XYZ b = global2local(cs, getAtom("C8")->coord);
				XYZ c = global2local(cs, getAtom("N7")->coord);
				XYZ d = global2local(cs, getAtom("N6")->coord);
				XYZ e = global2local(cs, getAtom("N1")->coord);
				XYZ f = global2local(cs, getAtom("C2")->coord);
				XYZ g = global2local(cs, getAtom("N3")->coord);
				
    			points.push_back(XY(a.x_, a.y_));
    			points.push_back(XY(b.x_, b.y_));
    			points.push_back(XY(c.x_, c.y_));
    			points.push_back(XY(d.x_, d.y_));
    			points.push_back(XY(e.x_, e.y_));
    			points.push_back(XY(f.x_, f.y_));
    			points.push_back(XY(g.x_, g.y_));
	}
	else if(baseTypeInt == 2 || baseTypeInt == 6){

				XYZ a = global2local(cs, getAtom("N9")->coord);
				XYZ b = global2local(cs, getAtom("C8")->coord);
				XYZ c = global2local(cs, getAtom("N7")->coord);
				XYZ d = global2local(cs, getAtom("O6")->coord);
				XYZ e = global2local(cs, getAtom("N2")->coord);
				XYZ f = global2local(cs, getAtom("N3")->coord);
				
    			points.push_back(XY(a.x_, a.y_));
    			points.push_back(XY(b.x_, b.y_));
    			points.push_back(XY(c.x_, c.y_));
    			points.push_back(XY(d.x_, d.y_));
    			points.push_back(XY(e.x_, e.y_));
    			points.push_back(XY(f.x_, f.y_));
	}
	else if(baseTypeInt == 3 || baseTypeInt == 7){

				XYZ a = global2local(cs, getAtom("N1")->coord);
				XYZ b = global2local(cs, getAtom("O2")->coord);
				XYZ c = global2local(cs, getAtom("N4")->coord);
				XYZ d = global2local(cs, getAtom("C5")->coord);
				XYZ e = global2local(cs, getAtom("C6")->coord);
				
    			points.push_back(XY(a.x_, a.y_));
    			points.push_back(XY(b.x_, b.y_));
    			points.push_back(XY(c.x_, c.y_));
    			points.push_back(XY(d.x_, d.y_));
    			points.push_back(XY(e.x_, e.y_));
	}
	else if(baseTypeInt == 1 || baseTypeInt == 5){
				XYZ a = global2local(cs, getAtom("N1")->coord);
				XYZ b = global2local(cs, getAtom("O2")->coord);
				XYZ c = global2local(cs, getAtom("O4")->coord);
				XYZ d = global2local(cs, getAtom("C5")->coord);
				XYZ e = global2local(cs, getAtom("C6")->coord);
				
    			points.push_back(XY(a.x_, a.y_));
    			points.push_back(XY(b.x_, b.y_));
    			points.push_back(XY(c.x_, c.y_));
    			points.push_back(XY(d.x_, d.y_));
    			points.push_back(XY(e.x_, e.y_));
	}

	return ConvexPolygon(points);
}

bool RNABase::isStackingTo(RNABase* other, AtomLib* atLib){

	LocalFrame csA = getCoordSystem();
	LocalFrame csB = other->getCoordSystem();
	if(csA.origin_.distance(csB.origin_) > 15.0) return false;
	double planeDistance = this->planeDistance(other);
	if(planeDistance > 4.0) return false;
	double planeAng = planeAngle(other);
	if(planeAng > 40.0) return false;


	TransMatrix tmA = csA.localToGlobalTM;
	XYZ z1 = XYZ(tmA.mtx[0][2], tmA.mtx[1][2], tmA.mtx[2][2]);
	TransMatrix tmB = csB.localToGlobalTM;
	XYZ z2 = XYZ(tmB.mtx[0][2], tmB.mtx[1][2], tmB.mtx[2][2]);
	double ang = angleX(z1, z2);


	XYZ meanZ;
	if(ang > 90.0)
		meanZ = z1 - z2;
	else
		meanZ = z1 + z2;
	
	double len = meanZ.length();
	if(len == 0)
		meanZ = XYZ(0,0,1);
	else 
		meanZ = ~meanZ;
		


	XYZ x = XYZ(tmA.mtx[0][0], tmA.mtx[1][0], tmA.mtx[2][0]);
	XYZ y = meanZ ^ x;
	XYZ origin = XYZ(0,0,0);
	LocalFrame cs = LocalFrame(origin, x, y, meanZ);

	ConvexPolygon cpA = getBaseConvexPolygon(cs, atLib);
    ConvexPolygon cpB = other->getBaseConvexPolygon(cs, atLib);

    if(cpA.overlap(cpB))
    	return true;

    return false;
}

int RNABase::printPDBFormat(ostream& out, int startAtomID) const{
    char c = this->baseID.at(baseID.length()-1);
    char s[100];
    vector<Atom*>::const_iterator it;
    int atomID = startAtomID;
    if(c >= '0' && c <= '9')
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s  %c %c%4s    %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->baseType,this->chainID[0],this->baseID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    else
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s  %c %c%5s   %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->baseType,this->chainID[0],this->baseID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    return atomID;
}

int RNABase::printCIFFormat(ostream& out, int startAtomID) const{
    char c = this->baseID.at(baseID.length()-1);
    char s[100];
    vector<Atom*>::const_iterator it;
    int atomID = startAtomID;
	char insCode;
    if(c >= '0' && c <= '9')
    {
		 insCode = '.';
    }
	else
    {
		 insCode = c;
	}
    for(it=atomList.begin();it!=atomList.end();it++)
	{
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
			string baseTypeString;
			if(baseType == 'A')
				baseTypeString = "A";
			else if(baseType == 'U')
				baseTypeString = "U";
			else if(baseType == 'G')
				baseTypeString = "G";
			else if(baseType == 'C')
				baseTypeString = "C";
			else if(baseType == 'a')
				baseTypeString = "DA";				
			else if(baseType == 't')
				baseTypeString = "DT";
			else if(baseType == 'g')
				baseTypeString = "DG";
			else if(baseType == 'c')
				baseTypeString = "DC";

			string atomName = (*it)->getName();
			cout << atomName << " " << atomName.length() << endl;
			
			if(atomName[atomName.length()-1] == '\'')
			{
				atomName = "\"" + atomName + "\"";
				cout << atomName << endl;
			}	

			sprintf(s, "ATOM   %-7d %s %-5s %c %s  %s  1 %-4s ? %8.3f %8.3f %8.3f %4.2f %6.2f ? %s %s %s %-5s 1", atomID, (*it)->type.c_str(), atomName.c_str(), insCode, baseTypeString.c_str(), chainID.c_str(), baseID.c_str(),coord[0],coord[1],coord[2], 1.00, 0.00, baseID.c_str(), baseTypeString.c_str(), chainID.c_str(), atomName.c_str());
            out << s << endl;
            atomID++;
    }

    return atomID;
}

string RNABase::print() {
	string ret = chainID + "_" + baseType + baseID;
	return ret;
}

RNABase::~RNABase() {
}

PolarAtom::PolarAtom(Residue* res, string atomName) {

	this->uniqueName = res->triName + "-" + atomName;
	this->isDonor = false;
	this->isAcceptor = false;

	if(uniqueName == "ASP-OD1" || uniqueName == "ASP-OD2" || uniqueName == "ASN-OD1")
	{
		this->support = res->getAtom("CG")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isAcceptor = true;
	}
	else if(uniqueName == "ASN-ND2")
	{
		this->support = res->getAtom("CG")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "GLU-OE1" || uniqueName == "GLU-OE2" || uniqueName == "GLN-OE1")
	{
		this->support = res->getAtom("CD")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isAcceptor = true;
	}
	else if(uniqueName == "GLN-NE2")
	{
		this->support = res->getAtom("CD")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "HIS-ND1")
	{
		XYZ n =  res->getAtom("ND1")->getCoord();
		XYZ sup1 = res->getAtom("CG")->getCoord();
		XYZ sup2 = res->getAtom("CE1")->getCoord();

		this->core = n;
		this->support = sup1 + sup2 -n;

		this->isAcceptor = true;
		this->isDonor = true;
	//	this->HCore = true;

	}
	else if(uniqueName == "HIS-NE2")
	{
		XYZ n = res->getAtom("NE2")->getCoord();

		XYZ sup1 = res->getAtom("CD2")->getCoord();
		XYZ sup2 = res->getAtom("CE1")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isAcceptor = true;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "LYS-NZ")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CE")->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "ARG-NE")
	{
		XYZ n = res->getAtom("NE")->getCoord();
		XYZ sup1 = res->getAtom("CD")->getCoord();
		XYZ sup2 = res->getAtom("CZ")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "ARG-NH1" || uniqueName == "ARG-NH2")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CZ")->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "SER-OG" || uniqueName == "THR-OG1")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CB")->getCoord();
		this->isAcceptor = true;
		this->isDonor = true;
	}
	else if(uniqueName == "TRP-NE1")
	{
		XYZ n = res->getAtom("NE1")->getCoord();
		XYZ sup1 = res->getAtom("CD1")->getCoord();
		XYZ sup2 = res->getAtom("CE2")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "TYR-OH")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CZ")->getCoord();
		this->isAcceptor = true;
		this->isDonor = true;
	}
	else if(atomName == "O" || atomName == "OXT" || atomName == "OT1" || atomName == "OT2")
	{
		this->uniqueName = res->triName + "-O";
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("C")->getCoord();
		this->isAcceptor = true;
	}
	else
	{
		cerr << "not polar atom " << uniqueName << endl;
		exit(0);
	}

	if(atomName.at(0) == 'O')
		this->vdwRadius = 1.6;
	else
		this->vdwRadius = 1.7;

}

PolarAtom::PolarAtom(RNABase* base, string atomName){

	char xx[20];
	this->isAcceptor = false;
	this->isDonor = false;
	this->core = XYZ(0, 0, 0);
	this->support = XYZ(0, 0, 0);

	sprintf(xx, "%c-%s", base->baseType, atomName.c_str());
	this->uniqueName = string(xx);

	if(atomName == "O2'" && base->getAtom("O2'") != NULL && base->getAtom("C2'") != NULL){
		XYZ c = base->getAtom("O2'")->getCoord();
		XYZ sup = base->getAtom("C2'")->getCoord();
		this->core = c;
		this->support = sup;
		this->isAcceptor = true;
		this->isDonor = true;
	}
	else if(atomName == "OP1" && base->getAtom("P") != NULL && base->getAtom("OP1") != NULL){
		XYZ c = base->getAtom("OP1")->getCoord();
		XYZ sup = base->getAtom("P")->getCoord();
		this->core = c;
		this->support = sup;
		this->isAcceptor = true;
	}
	else if(atomName == "OP2" && base->getAtom("P") != NULL && base->getAtom("OP2") != NULL){
		XYZ c = base->getAtom("OP2")->getCoord();
		XYZ sup = base->getAtom("P")->getCoord();
		this->core = c;
		this->support = sup;
		this->isAcceptor = true;
	}

	if(base->baseTypeInt == 0 || base->baseTypeInt == 4){
		if(atomName == "N1"){
			if(base->getAtom("N1") != NULL && base->getAtom("C2") != NULL && base->getAtom("C6") != NULL){
				XYZ c = base->getAtom("N1")->coord;
				XYZ sup1 = base->getAtom("C2")->coord;
				XYZ sup2 = base->getAtom("C6")->coord;
				this->core = c;
				this->support = sup1 + sup2 - c;
				this->isDonor = true;
			}

		}
		else if(atomName == "N3"){
			if(base->getAtom("N3") != NULL && base->getAtom("C2") != NULL && base->getAtom("C4") != NULL){
				XYZ c = base->getAtom("N3")->coord;
				XYZ sup1 = base->getAtom("C2")->coord;
				XYZ sup2 = base->getAtom("C4")->coord;
				this->core = c;
				this->support = sup1 + sup2 - c;
				this->isDonor = true;
			}
		}
		else if(atomName == "N6"){
			if(base->getAtom("N6") != NULL && base->getAtom("C6") != NULL){
				XYZ c = base->getAtom("N6")->coord;
				XYZ sup = base->getAtom("C6")->coord;
				this->core = c;
				this->support = sup;
				this->isDonor = true;
			}

		}
		else if(atomName == "N7"){
			if(base->getAtom("N7") != NULL && base->getAtom("C2") != NULL && base->getAtom("C8") != NULL){
				XYZ c = base->getAtom("N7")->coord;
				XYZ sup1 = base->getAtom("C5")->coord;
				XYZ sup2 = base->getAtom("C8")->coord;
				this->core = c;
				this->support = sup1 + sup2 - c;
				this->isDonor = true;
			}
		}
	}
	else if(base->baseTypeInt == 1 || base->baseTypeInt == 5){
		if(atomName == "O2" && base->getAtom("O2") != NULL && base->getAtom("C2") != NULL){
			XYZ c = base->getAtom("O2")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			this->core = c;
			this->support = sup1;
			this->isAcceptor = true;
		}
		else if(atomName == "N3" && base->getAtom("N3") != NULL && base->getAtom("C2") != NULL && base->getAtom("C4") != NULL){
			XYZ c = base->getAtom("N3")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "O4" && base->getAtom("O4") != NULL && base->getAtom("C4") != NULL){
			XYZ c = base->getAtom("O4")->coord;
			XYZ sup = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup;
			this->isAcceptor = true;
		}
	}
	else if(base->baseTypeInt == 2 || base->baseTypeInt == 6){
		if(atomName == "N1" && base->getAtom("N1") != NULL && base->getAtom("C2") != NULL && base->getAtom("C6") != NULL ){
			XYZ c = base->getAtom("N1")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C6")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "N2" && base->getAtom("N2") != NULL && base->getAtom("C2") != NULL){
			XYZ c = base->getAtom("N2")->coord;
			XYZ sup = base->getAtom("C2")->coord;
			this->core = c;
			this->support = sup;
			this->isDonor = true;
		}
		else if(atomName == "N3" && base->getAtom("N3") != NULL && base->getAtom("C2") != NULL && base->getAtom("C4") != NULL ){
			XYZ c = base->getAtom("N3")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "O6" && base->getAtom("O6") != NULL && base->getAtom("C6") != NULL){
			XYZ c = base->getAtom("O6")->coord;
			XYZ sup = base->getAtom("C6")->coord;
			this->core = c;
			this->support = sup;
			this->isAcceptor = true;
		}
		else if(atomName == "N7" && base->getAtom("N7") != NULL && base->getAtom("C5") != NULL && base->getAtom("C8") != NULL){
			XYZ c = base->getAtom("N7")->coord;
			XYZ sup1 = base->getAtom("C5")->coord;
			XYZ sup2 = base->getAtom("C8")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
	}
	else if(base->baseTypeInt == 3 || base->baseTypeInt == 7){
		if(atomName == "O2" && base->getAtom("O2") != NULL && base->getAtom("C2") != NULL){
			XYZ c = base->getAtom("O2")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			this->core = c;
			this->support = sup1;
			this->isAcceptor = true;
		}
		else if(atomName == "N3" && base->getAtom("N3") != NULL && base->getAtom("C2") != NULL && base->getAtom("C4") != NULL){
			XYZ c = base->getAtom("N3")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "N4" && base->getAtom("N4") != NULL && base->getAtom("C4") != NULL){
			XYZ c = base->getAtom("N4")->coord;
			XYZ sup = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup;
			this->isAcceptor = true;
		}
	}
	if(atomName[0] == 'N')
		this->vdwRadius = 1.7;
	else
		this->vdwRadius = 1.6;

}

PolarAtom::PolarAtom(Residue* res, string atomName, Atom* preC)
{
	if(atomName != "N" || res->triName == "PRO")
	{
		cerr << "constructor for only backbone polar N" << res->triName << "-" <<  atomName << endl;
		exit(0);
	}
	this->uniqueName = res->triName + "-N";
	XYZ n = res->getAtom("N")->getCoord();
	XYZ sup1 = res->getAtom("CA")->getCoord();
	XYZ sup2 = preC->getCoord();
	this->core = n;
	this->support = sup1 + sup2 -n;
	this->isDonor = true;
	this->isAcceptor = false;
	this->vdwRadius = 1.7;
//	this->HCore = true;
}

PolarAtom::~PolarAtom() {
}

ProteinChain::ProteinChain() {
	this->pdbID = "xxxx";
	this->chainID = '-';
	this->chainLen = 0;
}

ProteinChain::ProteinChain(string chainID) {
	this->pdbID = "xxxx";
	this->chainID = chainID;
	this->chainLen = 0;
}

ProteinChain::ProteinChain(string pdbID, string chainID){
	this->pdbID = pdbID;
	this->chainID = chainID;
	this->chainLen = 0;
}

void ProteinChain::setPDBID(string pdbID){
	this->pdbID = pdbID;
}

void ProteinChain::setChainID(string c){
	this->chainID = c;
}

string ProteinChain::getPDBID() const{
	return this->pdbID;
}

string ProteinChain::getChainID() const{
	return this->chainID;
}

int ProteinChain::getChainLength() const{
	return this->chainLen;
}

vector<Residue*>& ProteinChain::getResList(){
	return this->resList;
}

Residue* ProteinChain::getResidue(const string& resID) {
	map<string,Residue*>::const_iterator it = resMap.find(resID);
	if(it != resMap.end())
		return it->second;
	else
		return NULL;
}

void ProteinChain::addResidue(Residue* res){
	this->resList.push_back(res);
	this->resMap[res->resID] = res;
	res->setResSeqID(this->resList.size()-1);
	this->chainLen ++;
}

string ProteinChain::getSequence() const{
	string s = "";
	ResName rn = ResName();
	for(int i=0;i<chainLen;i++){
		char c = rn.triToSin(this->resList.at(i)->triName);
		s = s + c;
	}
	return s;
}

int ProteinChain::printPDBFormat(ofstream& out, int startAtomID) const{

	for(int i=0;i<resList.size();i++){
		startAtomID = resList.at(i)->printPDBFormat(out,startAtomID);
	}
	return startAtomID;
}

int ProteinChain::printPDBFormatNoHydrogen(ofstream& out, int startAtomID) const{

	return 0;
}

ProteinChain::~ProteinChain() {
}

RNAChain::RNAChain() {
	this->pdbID = "pdbx";
	this->chainID = "A";
	this->chainLen = 0;
}


RNAChain::RNAChain(const string& chainID) {
	this->pdbID = "pdbx";
	this->chainID = chainID;
	this->chainLen = 0;
}

RNAChain::RNAChain(const string& pdbID, const string& chainID) {
	this->pdbID = pdbID;
	this->chainID = chainID;
	this->chainLen = 0;
}

RNAChain::RNAChain(const RNAChain& other) {
	this->pdbID = other.pdbID;
	this->chainID = other.chainID;
	this->chainLen = other.chainLen;

	int lb = other.baseList.size();
	map<RNABase*, RNABase*> old2newBsMap;
	for(int i=0; i<lb; i++) {
		this->baseList.emplace_back(new RNABase{*other.baseList[i]});
		old2newBsMap.emplace(other.baseList[i], this->baseList[i]);
	}

	for(auto & iter : other.baseMap) {
		this->baseMap.emplace(iter.first, old2newBsMap.at(iter.second));
	}
}

RNAChain::RNAChain(RNAChain&& other) noexcept {
	this->pdbID = move(other.pdbID);
	this->chainID = move(other.chainID);
	this->chainLen = other.chainLen;
	this->baseList = move(other.baseList);
	this->baseMap = move(other.baseMap);
}

int RNAChain::printPDBFormat(ostream& out, int startAtomID) const{
	for(int i=0;i<baseList.size();i++) {
		startAtomID = baseList[i]->printPDBFormat(out, startAtomID);
	}
	return startAtomID;
}

int RNAChain::printCIFFormat(ostream& out, int startAtomID) const{

	out << "loop_" << endl;
	out << "_atom_site.group_PDB" << endl;
	out << "_atom_site.id" << endl;
	out << "_atom_site.type_symbol" << endl;
	out << "_atom_site.label_atom_id" << endl;
	out << "_atom_site.label_alt_id" << endl;
	out << "_atom_site.label_comp_id" << endl;
	out << "_atom_site.label_asym_id" << endl;
	out << "_atom_site.label_entity_id" << endl;
	out << "_atom_site.label_seq_id" << endl;
	out << "_atom_site.pdbx_PDB_ins_code" << endl;
	out << "_atom_site.Cartn_x" << endl;
	out << "_atom_site.Cartn_y" << endl;
	out << "_atom_site.Cartn_z" << endl;
	out << "_atom_site.occupancy" << endl;
	out << "_atom_site.B_iso_or_equiv" << endl;
	out << "_atom_site.pdbx_formal_charge" << endl;
	out << "_atom_site.auth_seq_id" << endl;
	out << "_atom_site.auth_comp_id" << endl;
	out << "_atom_site.auth_asym_id" << endl;
	out << "_atom_site.auth_atom_id" << endl;
	out << "_atom_site.pdbx_PDB_model_num" << endl;

	for(int i=0;i<baseList.size();i++) {
		startAtomID = baseList[i]->printCIFFormat(out, startAtomID);
	}
	return startAtomID;
}

PDB::PDB() {
	this->pdbID = "pdbx";
}

PDB::PDB(const string& pdbFile, const string& pdbID)
{
	this->pdbID = pdbID;

}

void PDB::readPDB(const string& pdbFile){
	ifstream input;
	input.open(pdbFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << pdbFile << endl;
        exit (1);
    }
    string s;
    int models = 0;
    int len;
    string lastChainID = "@";
    string lastResID = "XXX";

    string curChainID;
    string curResID;
    char altLoc;
    string resName;

    ProteinChain* curChain;
    Residue* curResidue;

    ResName rn = ResName();
	while(getline(input,s))
    {
        len = s.length();
        if(len < 6) continue;
        string prefix = s.substr(0,6);
        if(prefix == "MODEL ")
        {
            if(models > 0)
                break;
            else
                models = 1;
        }
        if(len < 54) continue;
        if(prefix != "ATOM  " && prefix != "HETATM") continue;
        resName = s.substr(17,3);
        if(resName == "HOH") continue;

        if(!rn.isStandardAminoAcid(resName) && resName != "MSE") continue;


        curChainID = s[21];
        curResID = trimString(s.substr(22,5));
        altLoc = s[16];


        Atom* a = new Atom(s);
        if(curChainID != lastChainID)
        {
             if(getChain(curChainID) == NULL){
                  curChain = new ProteinChain(curChainID);
                  this->chains.push_back(curChain);
             }
             else
            	  curChain = getChain(curChainID);

             lastChainID = curChainID;
             lastResID = "XXX";
        }

        if(curResID != lastResID)
        {
        	curResidue = new Residue(curResID,curChainID,resName);
            curChain->addResidue(curResidue);
            residues.push_back(curResidue);
            lastResID = curResID;
        }
        if(altLoc != ' ' && altLoc != 'A' && altLoc != '1') {
        	curResidue->hasAltConf = true;
        	continue;
        }

        curResidue->addAtom(a);
        curResidue->setAltLoc(altLoc);
    }

    input.close();
}

void PDB::readCIF(const string& cifFile){
	ifstream input;
	input.open(cifFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << cifFile << endl;
        exit (1);
    }

    string s;
    int models = 0;
    int len;
    string lastChainID = "@";
    string lastResID = "XXX";

    string curChainID;
    string curResID;
    int curResSeqID = 0;
    string insCode;
	string altCode;
    string rawResName;
    string resName;

    ProteinChain* curChain;
    Residue* curResidue;
    ResName rn;

    vector<string> spt;
	map<string, int> atomSiteMap;

	int atomSiteIndex = 0;
	while(getline(input,s))
    {
		//cout << s << endl;
        len = s.length();
        if(len < 5) continue;
        string prefix = s.substr(0,5);
		if(prefix == "loop_"){
			atomSiteMap.clear();
			atomSiteIndex = 0;
			/*
			atomSiteMap["label_alt_id"] = 4;
			atomSiteMap["label_comp_id"] = 5;
			atomSiteMap["label_asym_id"] = 6;
			atomSiteMap["pdbx_PDB_ins_code"] = 9;
			atomSiteMap["auth_seq_id"] = 16; //resID
			atomSiteMap["auth_comp_id"] = 17; //resName
			atomSiteMap["auth_asym_id"] = 18; //chainID
			atomSiteMap["pdbx_PDB_model_num"] = 20;
			*/
		}
		else if(prefix == "_atom"){
			atomSiteMap[s.substr(11, s.length()-12)] = atomSiteIndex;
			atomSiteIndex++;
		}

		

        if(prefix != "ATOM " && prefix != "HETAT") continue;

        Atom* a = new Atom(s, 1);
        splitString(s, " ", &spt);
	
        if(spt.size() < atomSiteMap.size()) continue;
		if(atomSiteMap.find("pdbx_PDB_model_num") != atomSiteMap.end()) {
			int modelID = atoi(spt[atomSiteMap["pdbx_PDB_model_num"]].c_str());
			if(modelID  > 1 ) 
				continue;
		}

        rawResName = spt[atomSiteMap["label_comp_id"]];
        if(rawResName == "HOH") continue;

        rawResName = trimString(rawResName);

		if(!rn.isAminoAcid(rawResName)) continue;
		resName = rawResName;

        curChainID = spt[atomSiteMap["auth_asym_id"]];
        insCode = spt[atomSiteMap["pdbx_PDB_ins_code"]];
        if(insCode == "?")
        	insCode = "";
		
		curResID = spt[atomSiteMap["auth_seq_id"]]+insCode;

		altCode = spt[atomSiteMap["label_alt_id"]];

        if(a->type == "H") continue;
        if(curChainID != lastChainID)
        {
             if(getChain(curChainID) == NULL)
             {
                  curChain = new ProteinChain(curChainID);
                  this->chains.push_back(curChain);
             }
             else
            	  curChain = getChain(curChainID);

             lastChainID = curChainID;
             lastResID = "XXX";
        }

        if(curResID != lastResID)
        {
        	curResidue = new Residue(curResID, curChainID, resName);
        	curResidue->setResSeqID(curResSeqID);
        	curResSeqID++;
        	curChain->addResidue(curResidue);
            residues.push_back(curResidue);
            lastResID = curResID;
        }

        if(altCode != "." && altCode != "A" && altCode != "1") {
        	curResidue->hasAltConf = true;
        	continue;
        }

        curResidue->addAtom(a);
        curResidue->setAltLoc(altCode[0]);
    }

    input.close();
}

PDB& PDB::operator=(const PDB& other){
	cout << "operator '=' is inhibited in class PDB" << endl;
	abort();
	return *this;
}

vector<ProteinChain*>& PDB::getChains()
{
    return chains;
}

ProteinChain* PDB::getFirstChain()
{
    if(this->chains.size() > 0)
        return this->chains.at(0);
    return NULL;
}

ProteinChain* PDB::getChain(string& c)
{
    for(unsigned int i=0;i<this->chains.size();i++)
    {
        if(this->chains.at(i)->getChainID() == c)
            return chains.at(i);
    }
    return NULL;
}

string PDB::getFirstSeq()
{
    return getFirstChain()->getSequence();
}

vector<Residue*>& PDB::getResList()
{
    return residues;
}

void PDB::printPDBFormat(ofstream& out) const
{
    int startID = 1;
    for(unsigned int i=0;i<this->chains.size();i++)
    {
        ProteinChain* pc = this->chains.at(i);
        for(int j=0;j<pc->getChainLength();j++)
        {
            Residue* res = pc->getResList().at(j);
            startID = res->printPDBFormat(out,startID);
        }
    }
}


PDB::~PDB() {
	unsigned int i,j;
	ProteinChain* p;
	for(i = 0;i<chains.size();i++)
	{
		p = chains.at(i);
		delete p;

	}

	Residue* p2;
	Atom* p3;
	for(i=0;i<residues.size();i++)
	{

		p2 = residues.at(i);

		for(j=0;j<p2->getAtomList()->size();j++)
		{
			p3 = p2->getAtomList()->at(j);
			delete p3;

		}
		delete p2;
	}
}

RNAPDB::RNAPDB() {
	this->pdbID = "pdbx";
}

RNAPDB::RNAPDB(const string& pdbFile){
	this->pdbID = "pdbx";
	if(pdbFile.substr(pdbFile.length()-3, 3) == "pdb"){
		readPDB(pdbFile);
	}
	else if(pdbFile.substr(pdbFile.length()-3, 3) == "cif"){
		readCIF(pdbFile);
	}
}

RNAPDB::RNAPDB(const string& pdbFile, const string& pdbID) {

	this->pdbID = pdbID;
	if(pdbFile.substr(pdbFile.length()-3, 3) == "pdb"){
		readPDB(pdbFile);
	}
	else if(pdbFile.substr(pdbFile.length()-3, 3) == "cif"){
		readCIF(pdbFile);
	}
}

void RNAPDB::readPDB(const string& pdbFile){
	ifstream input;
	input.open(pdbFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << pdbFile << endl;
        exit (1);
    }
    string s;
    int models = 0;
    int len;
    string lastChainID = "@";
    string lastResID = "XXX";

    string curChainID;
    string curResID;
    int curResSeqID = 0;
    char altLoc;
    string rawResName;
    string resName;

    RNAChain* curChain;
    RNABase* curResidue;
	Ligand* curLigand = NULL;
    ResName rn;
	vector<string> spt;

	while(getline(input,s))
    {
		//cout << s << endl;
		if(s.length() > 3 && s.substr(0,3) == "ene") {
			splitString(s, " " , &spt);
			this->ene = atof(spt[1].c_str());
		}

        len = s.length();
        if(len < 6) continue;
        string prefix = s.substr(0,6);
        if(prefix == "MODEL ")
        {
            if(models > 0)
                break;
            else
                models = 1;
        }
        if(len < 54) continue;
        if(prefix != "ATOM  " && prefix != "HETATM") continue;
        rawResName = s.substr(17,3);
        if(rawResName == "HOH") continue;

        rawResName = trimString(rawResName);


        curChainID = s.substr(21,1);
        curResID = trimString(s.substr(22,5));
        altLoc = s.at(16);


        Atom* a = new Atom(s);
        if(rawResName == "4SU" && a->name == "S4"){
        	a->name = "O4";
        	a->type = "O";
        }
        if(rawResName == "QUO" && a->name == "C7"){
        	a->name = "N7";
        	a->type = "N";
        }

        if(a->type == "H") continue;

		if(rn.isRNABase(rawResName)){
			resName = rn.toStandardBase(rawResName);
        	if(curChainID != lastChainID)
       		{
            	if(getChain(curChainID) == NULL)
             	{
                	curChain = new RNAChain(curChainID);
                	this->chains.push_back(curChain);
             	}
             	else
            		curChain = getChain(curChainID);

             	lastChainID = curChainID;
             	lastResID = "XXX";
        	}

        	if(curResID != lastResID)
        	{
        		curResidue = new RNABase(curResID, curChainID, resName[0]);
				//cout << curResidue->getType() << endl;
        		curResidue->setResSeqID(curResSeqID);
        		curResSeqID++;
        		curChain->addBase(curResidue);
            	baseList.push_back(curResidue);
            	lastResID = curResID;
        	}
        	if(altLoc != ' ' && altLoc != 'A' && altLoc != '1') {
        		curResidue->hasAltConf = true;
        		continue;
        	}

        	curResidue->addAtom(a);
        	curResidue->setAltLoc(altLoc);
		}
		else if(rn.isAminoAcid(rawResName)) {
			continue; //ignore protein atoms
		}
		else {
			//read ligands
			if(curChainID != lastChainID || curResID != lastResID || curLigand == NULL){
				lastChainID = curChainID;
				lastResID = curResID;
				curLigand = new Ligand();
				curLigand->setChainID(curChainID);
				curLigand->setLigandName(rawResName);
				curLigand->setResID(curResID);
				this->ligList.push_back(curLigand);
			}

			curLigand->addAtom(a);
		}
    }

    input.close();
}

void RNAPDB::readCIF(const string& cifFile){
	ifstream input;
	input.open(cifFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << cifFile << endl;
        exit (1);
    }

    string s;
    int models = 0;
    int len;
    string lastChainID = "@";
    string lastResID = "XXX";

    string curChainID = "";
    string curResID = "";
	string curLigID = "";
    int curResSeqID = 0;
    string insCode;
	string altCode;
    string rawResName;
    string resName;

    RNAChain* curChain;
    RNABase* curResidue;
    ResName rn;

	Ligand* curLigand = NULL;

    vector<string> spt;
	map<string, int> atomSiteMap;

	int atomSiteIndex = 0;
	while(getline(input,s))
    {
		//cout << s << endl;

		if(s.length() > 3 && s.substr(0,3) == "ene") {
			splitString(s, " " , &spt);
			this->ene = atof(spt[1].c_str());
		}
		
        len = s.length();
        if(len < 5) continue;
        string prefix = s.substr(0,5);
		if(prefix == "loop_"){
			atomSiteMap.clear();
			atomSiteIndex = 0;
			/*
			atomSiteMap["label_alt_id"] = 4;
			atomSiteMap["label_comp_id"] = 5;
			atomSiteMap["label_asym_id"] = 6;
			atomSiteMap["pdbx_PDB_ins_code"] = 9;
			atomSiteMap["auth_seq_id"] = 16; //resID
			atomSiteMap["auth_comp_id"] = 17; //resName
			atomSiteMap["auth_asym_id"] = 18; //chainID
			atomSiteMap["pdbx_PDB_model_num"] = 20;
			*/
		}
		else if(prefix == "_atom"){
			atomSiteMap[s.substr(11, s.length()-12)] = atomSiteIndex;
			atomSiteIndex++;
		}

		

        if(prefix != "ATOM " && prefix != "HETAT") continue;

        Atom* a = new Atom(s, 1);
        splitString(s, " ", &spt);
	
        if(spt.size() < atomSiteMap.size()) continue;
		if(atomSiteMap.find("pdbx_PDB_model_num") != atomSiteMap.end()) {
			int modelID = atoi(spt[atomSiteMap["pdbx_PDB_model_num"]].c_str());
			if(modelID  > 1 ) 
				continue;
		}

        rawResName = spt[atomSiteMap["label_comp_id"]];
        if(rawResName == "HOH") continue;



        curChainID = spt[atomSiteMap["label_asym_id"]];

      
        insCode = spt[atomSiteMap["pdbx_PDB_ins_code"]];

        if(insCode == "?")
        	insCode = "";
		
		curResID = spt[atomSiteMap["label_seq_id"]]+insCode;
		altCode = spt[atomSiteMap["label_alt_id"]];

        if(rawResName == "4SU" && a->name == "S4"){
        	a->name = "O4";
        	a->type = "O";
        }
        if(rawResName == "QUO" && a->name == "C7"){
        	a->name = "N7";
        	a->type = "N";
        }

        if(a->type == "H") continue;

        rawResName = trimString(rawResName);
        if(!rn.isRNABase(rawResName)) continue;
        
		if(rn.isRNABase(rawResName)) {
			resName = rn.toStandardBase(rawResName);
        	if(curChainID != lastChainID)
        	{
             	if(getChain(curChainID) == NULL)
             	{
                	curChain = new RNAChain(curChainID);
                	this->chains.push_back(curChain);
             	}
             	else
            		curChain = getChain(curChainID);

            	lastChainID = curChainID;
             	lastResID = "XXX";
        	}

        	if(curResID != lastResID)
        	{
        		curResidue = new RNABase(curResID, curChainID, resName[0]);
        		curResidue->setResSeqID(curResSeqID);
        		curResSeqID++;
        		curChain->addBase(curResidue);
            	baseList.push_back(curResidue);
            	lastResID = curResID;
        	}

        	if(altCode != "." && altCode != "A" && altCode != "1") {
        		curResidue->hasAltConf = true;
        		continue;
        	}

        	curResidue->addAtom(a);
        	curResidue->setAltLoc(altCode[0]);
		}
		else if(rn.isAminoAcid(rawResName)) {
			continue; //ignore protein atoms
		}
		else {
			//read ligands
			if(curChainID != lastChainID || curResID != lastResID || curLigand == NULL){
				lastChainID = curChainID;
				lastResID = curResID;
				curLigand = new Ligand();
				curLigand->setChainID(curChainID);
				curLigand->setLigandName(rawResName);
				curLigand->setResID(curResID);
				this->ligList.push_back(curLigand);
			}

			curLigand->addAtom(a);
		}
    }

    input.close();
}

void RNAPDB::printPDBFormat(ostream& out) const{
    int startID = 1;
    for(unsigned int i=0;i<this->chains.size();i++)
    {
    	RNAChain* pc = this->chains.at(i);
        for(int j=0;j<pc->getChainLength();j++)
        {
            RNABase* res = pc->getBaseList().at(j);
            startID = res->printPDBFormat(out,startID);
        }
    }
}

void RNAPDB::printCIFFormat(ostream& out) const{
    int startID = 1;
    for(unsigned int i=0;i<this->chains.size();i++)
    {
    	RNAChain* pc = this->chains.at(i);
		startID = pc->printCIFFormat(out, startID);
    }
}

RNAPDB::~RNAPDB() {
	unsigned int i,j;
	RNAChain* p;
	for(i = 0;i<chains.size();i++)
	{
		p = chains.at(i);
		delete p;

	}

	RNABase* p2;
	Atom* p3;
	for(i=0;i<baseList.size();i++)
	{

		p2 = baseList.at(i);

		for(j=0;j<p2->getAtomList()->size();j++)
		{
			p3 = p2->getAtomList()->at(j);
			delete p3;

		}
		delete p2;
	}
}


float Phipsi::distance(const Phipsi& other) const{
	float delX = this->phi - other.phi;
	float delY = this->psi - other.psi;
	if(delX > 180)
		delX = 360 - delX;
	else if(delX < -180)
		delX = 360 + delX;

	if(delY > 180)
		delY = 360 - delY;
	else if(delY < -180)
		delY = 360 + delY;
	return sqrt(delX*delX+delY*delY);
}

char Phipsi::regionAB() const{
	if(phi < 0 && psi > -100 && psi < 60)
	{
		return 'A';
	}
	return 'B';
}

Phipsi::~Phipsi() {
}

PhipsiLib::PhipsiLib(){
	this->pointNum = 200;
	this->ppList.push_back(new Phipsi( -62.71 ,  -41.57));
	this->ppList.push_back(new Phipsi( -60.44 ,  -36.89));
	this->ppList.push_back(new Phipsi( -65.79 ,  -36.65));
	this->ppList.push_back(new Phipsi( -67.17 ,  -41.65));
	this->ppList.push_back(new Phipsi( -59.91 ,  -45.30));
	this->ppList.push_back(new Phipsi( -56.70 ,  -41.11));
	this->ppList.push_back(new Phipsi( -64.61 ,  -46.65));
	this->ppList.push_back(new Phipsi( -59.68 ,  -50.44));
	this->ppList.push_back(new Phipsi( -54.61 ,  -46.94));
	this->ppList.push_back(new Phipsi( -71.59 ,  -46.40));
	this->ppList.push_back(new Phipsi( -72.03 ,  -37.93));
	this->ppList.push_back(new Phipsi( -62.47 ,  -30.49));
	this->ppList.push_back(new Phipsi( -54.86 ,  -32.60));
	this->ppList.push_back(new Phipsi( -69.62 ,  -31.26));
	this->ppList.push_back(new Phipsi( -48.60 ,  -39.19));
	this->ppList.push_back(new Phipsi( -65.64 ,  -55.36));
	this->ppList.push_back(new Phipsi( -53.82 ,  -55.75));
	this->ppList.push_back(new Phipsi( -46.67 ,  -48.77));
	this->ppList.push_back(new Phipsi( -57.50 ,  -24.42));
	this->ppList.push_back(new Phipsi( -66.32 ,  -23.90));
	this->ppList.push_back(new Phipsi( -78.56 ,  -32.35));
	this->ppList.push_back(new Phipsi( -81.49 ,  -42.96));
	this->ppList.push_back(new Phipsi( -75.23 ,  -23.39));
	this->ppList.push_back(new Phipsi( -62.73 ,  -17.29));
	this->ppList.push_back(new Phipsi( -70.93 ,  -15.55));
	this->ppList.push_back(new Phipsi( -80.44 ,  -56.71));
	this->ppList.push_back(new Phipsi( -86.95 ,  -23.82));
	this->ppList.push_back(new Phipsi( -92.42 ,  -35.65));
	this->ppList.push_back(new Phipsi( -81.19 ,  -14.33));
	this->ppList.push_back(new Phipsi( -66.67 ,   -7.74));
	this->ppList.push_back(new Phipsi( -76.53 ,   -6.37));
	this->ppList.push_back(new Phipsi( -96.72 ,  -50.52));
	this->ppList.push_back(new Phipsi( -91.95 ,  -12.06));
	this->ppList.push_back(new Phipsi(-100.84 ,  -23.37));
	this->ppList.push_back(new Phipsi( -85.93 ,   -3.37));
	this->ppList.push_back(new Phipsi(-108.71 ,  -38.01));
	this->ppList.push_back(new Phipsi( -76.47 ,    6.23));
	this->ppList.push_back(new Phipsi(-104.47 ,  -10.15));
	this->ppList.push_back(new Phipsi( -96.36 ,   -0.46));
	this->ppList.push_back(new Phipsi( -89.35 ,    6.81));
	this->ppList.push_back(new Phipsi(-116.16 ,  -21.63));
	this->ppList.push_back(new Phipsi( -38.10 ,  -59.62));
	this->ppList.push_back(new Phipsi(-107.80 ,    3.09));
	this->ppList.push_back(new Phipsi(-100.12 ,   11.45));
	this->ppList.push_back(new Phipsi(-117.52 ,   -6.73));
	this->ppList.push_back(new Phipsi(-117.12 ,  -57.58));
	this->ppList.push_back(new Phipsi( -95.66 ,  -77.23));
	this->ppList.push_back(new Phipsi( -91.15 ,   20.75));
	this->ppList.push_back(new Phipsi(-129.66 ,  -36.11));
	this->ppList.push_back(new Phipsi(-112.26 ,   14.81));
	this->ppList.push_back(new Phipsi(-133.22 ,  -11.58));
	this->ppList.push_back(new Phipsi(-124.53 ,    6.65));
	this->ppList.push_back(new Phipsi(-105.87 ,   25.91));
	this->ppList.push_back(new Phipsi( -58.02 ,  -83.01));
	this->ppList.push_back(new Phipsi(-124.36 ,   23.11));
	this->ppList.push_back(new Phipsi( -78.31 ,   41.41));
	this->ppList.push_back(new Phipsi(-142.95 ,  -64.95));
	this->ppList.push_back(new Phipsi(-143.81 ,   13.38));
	this->ppList.push_back(new Phipsi(-119.01 ,   38.30));
	this->ppList.push_back(new Phipsi(-119.93 ,  -94.85));
	this->ppList.push_back(new Phipsi(-138.92 ,   38.56));
	this->ppList.push_back(new Phipsi(-172.29 ,  -27.81));
	this->ppList.push_back(new Phipsi(  -5.37 ,  -81.80));
	this->ppList.push_back(new Phipsi(-109.17 ,  132.79));
	this->ppList.push_back(new Phipsi(-108.81 ,  123.49));
	this->ppList.push_back(new Phipsi(-117.57 ,  127.44));
	this->ppList.push_back(new Phipsi(-118.05 ,  136.97));
	this->ppList.push_back(new Phipsi(-118.00 ,  118.36));
	this->ppList.push_back(new Phipsi(-125.28 ,  131.48));
	this->ppList.push_back(new Phipsi(-127.56 ,  122.79));
	this->ppList.push_back(new Phipsi(-108.62 ,  114.24));
	this->ppList.push_back(new Phipsi(-100.17 ,  128.43));
	this->ppList.push_back(new Phipsi( -99.45 ,  118.02));
	this->ppList.push_back(new Phipsi(-110.11 ,  143.24));
	this->ppList.push_back(new Phipsi(-100.40 ,  139.45));
	this->ppList.push_back(new Phipsi(-128.65 ,  139.96));
	this->ppList.push_back(new Phipsi(-121.40 ,  146.83));
	this->ppList.push_back(new Phipsi(-134.95 ,  131.25));
	this->ppList.push_back(new Phipsi(-127.40 ,  111.21));
	this->ppList.push_back(new Phipsi(-115.32 ,  105.65));
	this->ppList.push_back(new Phipsi(-138.60 ,  119.63));
	this->ppList.push_back(new Phipsi(-138.74 ,  141.26));
	this->ppList.push_back(new Phipsi(-131.24 ,  150.02));
	this->ppList.push_back(new Phipsi(-113.22 ,  153.33));
	this->ppList.push_back(new Phipsi( -91.01 ,  134.29));
	this->ppList.push_back(new Phipsi( -90.81 ,  123.32));
	this->ppList.push_back(new Phipsi(-101.02 ,  152.94));
	this->ppList.push_back(new Phipsi(-123.72 ,  157.97));
	this->ppList.push_back(new Phipsi( -90.39 ,  145.82));
	this->ppList.push_back(new Phipsi(-100.97 ,  105.75));
	this->ppList.push_back(new Phipsi( -90.37 ,  110.74));
	this->ppList.push_back(new Phipsi( -81.99 ,  129.04));
	this->ppList.push_back(new Phipsi( -80.58 ,  139.50));
	this->ppList.push_back(new Phipsi( -81.23 ,  117.30));
	this->ppList.push_back(new Phipsi(-110.58 ,  164.44));
	this->ppList.push_back(new Phipsi( -88.14 ,  157.53));
	this->ppList.push_back(new Phipsi( -78.65 ,  149.95));
	this->ppList.push_back(new Phipsi( -73.22 ,  125.73));
	this->ppList.push_back(new Phipsi( -96.67 ,  167.33));
	this->ppList.push_back(new Phipsi( -71.22 ,  135.23));
	this->ppList.push_back(new Phipsi( -70.37 ,  144.53));
	this->ppList.push_back(new Phipsi(-134.04 ,  160.70));
	this->ppList.push_back(new Phipsi(-140.74 ,  152.22));
	this->ppList.push_back(new Phipsi(-123.92 ,  169.98));
	this->ppList.push_back(new Phipsi( -76.83 ,  160.13));
	this->ppList.push_back(new Phipsi(-146.93 ,  130.97));
	this->ppList.push_back(new Phipsi( -81.98 ,  169.43));
	this->ppList.push_back(new Phipsi( -68.82 ,  154.34));
	this->ppList.push_back(new Phipsi(-108.56 , -179.84));
	this->ppList.push_back(new Phipsi(-149.54 ,  142.95));
	this->ppList.push_back(new Phipsi(-136.72 ,  172.17));
	this->ppList.push_back(new Phipsi(-145.19 ,  163.00));
	this->ppList.push_back(new Phipsi(-151.06 ,  154.55));
	this->ppList.push_back(new Phipsi(-108.83 ,   93.20));
	this->ppList.push_back(new Phipsi( -68.43 ,  116.31));
	this->ppList.push_back(new Phipsi( -78.79 ,  103.37));
	this->ppList.push_back(new Phipsi( -63.16 ,  128.77));
	this->ppList.push_back(new Phipsi( -90.94 ,   95.27));
	this->ppList.push_back(new Phipsi( -61.74 ,  138.40));
	this->ppList.push_back(new Phipsi( -61.79 ,  147.62));
	this->ppList.push_back(new Phipsi( -68.50 ,  165.94));
	this->ppList.push_back(new Phipsi(-126.07 ,   94.90));
	this->ppList.push_back(new Phipsi( -90.20 , -176.79));
	this->ppList.push_back(new Phipsi(-141.05 ,  103.07));
	this->ppList.push_back(new Phipsi( -58.84 ,  157.36));
	this->ppList.push_back(new Phipsi( -53.54 ,  132.85));
	this->ppList.push_back(new Phipsi( -52.94 ,  143.20));
	this->ppList.push_back(new Phipsi( -74.20 ,  179.78));
	this->ppList.push_back(new Phipsi( -54.33 ,  121.37));
	this->ppList.push_back(new Phipsi(-126.43 , -170.97));
	this->ppList.push_back(new Phipsi(-155.86 ,  116.09));
	this->ppList.push_back(new Phipsi(-160.77 ,  148.35));
	this->ppList.push_back(new Phipsi(-162.21 ,  134.10));
	this->ppList.push_back(new Phipsi(-153.57 ,  171.62));
	this->ppList.push_back(new Phipsi(-159.20 ,  161.49));
	this->ppList.push_back(new Phipsi(-145.25 , -175.88));
	this->ppList.push_back(new Phipsi( -80.19 ,   82.18));
	this->ppList.push_back(new Phipsi(-122.33 ,   77.75));
	this->ppList.push_back(new Phipsi( -93.64 ,   75.65));
	this->ppList.push_back(new Phipsi(-138.24 ,   79.85));
	this->ppList.push_back(new Phipsi( -54.25 ,   96.28));
	this->ppList.push_back(new Phipsi( -40.61 ,  127.84));
	this->ppList.push_back(new Phipsi(-104.40 , -156.29));
	this->ppList.push_back(new Phipsi( -81.68 , -159.28));
	this->ppList.push_back(new Phipsi(-173.28 ,  156.71));
	this->ppList.push_back(new Phipsi(-167.83 ,  170.51));
	this->ppList.push_back(new Phipsi(-159.27 ,   91.56));
	this->ppList.push_back(new Phipsi(-162.77 , -174.54));
	this->ppList.push_back(new Phipsi(-148.38 , -155.90));
	this->ppList.push_back(new Phipsi( -81.65 ,   64.08));
	this->ppList.push_back(new Phipsi(-130.51 ,   61.59));
	this->ppList.push_back(new Phipsi(-152.79 ,   64.32));
	this->ppList.push_back(new Phipsi(-100.02 ,   51.12));
	this->ppList.push_back(new Phipsi(-122.55 , -135.41));
	this->ppList.push_back(new Phipsi( -92.38 , -124.40));
	this->ppList.push_back(new Phipsi(-157.80 , -114.40));
	this->ppList.push_back(new Phipsi(  79.97 ,   11.05));
	this->ppList.push_back(new Phipsi(  67.21 ,   13.11));
	this->ppList.push_back(new Phipsi(  75.20 ,   24.73));
	this->ppList.push_back(new Phipsi(  77.13 ,   -1.28));
	this->ppList.push_back(new Phipsi(  92.68 ,   16.58));
	this->ppList.push_back(new Phipsi(  91.44 ,    1.17));
	this->ppList.push_back(new Phipsi(  60.91 ,   25.89));
	this->ppList.push_back(new Phipsi(  67.07 ,   37.97));
	this->ppList.push_back(new Phipsi(  55.06 ,   37.75));
	this->ppList.push_back(new Phipsi(  87.84 ,  -12.94));
	this->ppList.push_back(new Phipsi( 111.79 ,    6.45));
	this->ppList.push_back(new Phipsi( 102.18 ,  -12.68));
	this->ppList.push_back(new Phipsi(  97.17 ,   46.11));
	this->ppList.push_back(new Phipsi(  58.70 ,   52.03));
	this->ppList.push_back(new Phipsi(  48.74 ,   46.96));
	this->ppList.push_back(new Phipsi(  43.41 ,   60.83));
	this->ppList.push_back(new Phipsi(  62.67 ,   78.25));
	this->ppList.push_back(new Phipsi( 104.48 ,  -30.61));
	this->ppList.push_back(new Phipsi(  73.14 ,  -44.39));
	this->ppList.push_back(new Phipsi( 131.11 ,  -17.09));
	this->ppList.push_back(new Phipsi( 160.05 ,   42.52));
	this->ppList.push_back(new Phipsi( 142.52 ,  -65.36));
	this->ppList.push_back(new Phipsi(   4.13 ,   98.81));
	this->ppList.push_back(new Phipsi( 123.61 , -160.45));
	this->ppList.push_back(new Phipsi(  92.55 , -154.26));
	this->ppList.push_back(new Phipsi( 100.14 ,  178.37));
	this->ppList.push_back(new Phipsi( 128.04 ,  168.58));
	this->ppList.push_back(new Phipsi(  76.51 , -174.33));
	this->ppList.push_back(new Phipsi(  84.21 ,  162.79));
	this->ppList.push_back(new Phipsi( 103.76 ,  146.02));
	this->ppList.push_back(new Phipsi(  67.39 , -154.11));
	this->ppList.push_back(new Phipsi( 152.15 , -168.57));
	this->ppList.push_back(new Phipsi( 114.53 , -122.60));
	this->ppList.push_back(new Phipsi( 146.81 , -140.10));
	this->ppList.push_back(new Phipsi(  77.05 , -123.32));
	this->ppList.push_back(new Phipsi(  56.48 , -134.67));
	this->ppList.push_back(new Phipsi( 159.19 ,  165.30));
	this->ppList.push_back(new Phipsi( 175.66 , -154.48));
	this->ppList.push_back(new Phipsi( 176.13 , -177.66));
	this->ppList.push_back(new Phipsi( 147.35 ,  128.99));
	this->ppList.push_back(new Phipsi(  52.69 , -117.81));
	this->ppList.push_back(new Phipsi( 106.36 ,  108.41));
	this->ppList.push_back(new Phipsi(  68.96 ,  118.77));
	this->ppList.push_back(new Phipsi(  81.03 ,  -73.88));
	creatIndexTable();
}
int PhipsiLib::findNearestPointsWithoutIndex(Phipsi* pp)
{
	float minDis = 10000.0;
	int minIndex = -1;
	float d;
	for(int i=0;i<pointNum;i++)
	{
		d = pp->distance(*(ppList.at(i)));
		//cout << d << endl;
		if(d < minDis)
		{
			minDis = d;
			minIndex = i;
		}
	}
	return minIndex;
}

void PhipsiLib::creatIndexTable()
{
	for(int i=0;i<36;i++)
	{
		for(int j=0;j<36;j++)
		{
			for(int k=0;k<20;k++)
			{
				this->indexTable[i][j][k] = -1;
			}
		}
	}

	set<unsigned int> possibleNeighbor;
	Phipsi *pp, *pp1, *pp2, *pp3, *pp4;
	set<unsigned int>::const_iterator it;

	for(int i=0;i<36;i++)
	{
		for(int j=0;j<36;j++)
		{
			possibleNeighbor.clear();
			float phi0 = i*10-180;
			float psi0 = j*10-180;
			for(unsigned int a=0;a<ppList.size();a++)
			{
				pp = ppList.at(a);
				if((pp->phi) > phi0 && (pp->phi) < phi0+10 && (pp->psi) > psi0 && (pp->psi) < psi0+10)
					possibleNeighbor.insert(a);
			}

			for(int a=0;a<11;a++)
			{
				pp1 = new Phipsi(phi0+a,psi0);
				pp2 = new Phipsi(phi0+a, psi0+10);
				pp3 = new Phipsi(phi0, psi0+a);
				pp4 = new Phipsi(phi0+10, psi0+a);
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp1));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp2));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp3));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp4));
				delete pp1;
				delete pp2;
				delete pp3;
				delete pp4;
			}

			int n=0;
			for(it = possibleNeighbor.begin();it!=possibleNeighbor.end();it++)
			{
				int k = *it;
				this->indexTable[i][j][n] = k;
				n++;
			}
		}
	}
}

Phipsi* PhipsiLib::indexToPhipsi(int id) const
{
	return this->ppList.at(id);
}

int PhipsiLib::phipsiToIndex(const Phipsi* pp) const
{
	int i = (int)((pp->phi+180)/10);
	int j = (int)((pp->psi+180)/10);
	double minD = 10000.0;
	int minIndex = -1;
	for(int k=0;k<20;k++)
	{
		int id = this->indexTable[i][j][k];
		if(id < 0)
			break;
		double dist = pp->distance(* ppList.at(id));
		if(dist < minD)
		{
			minD = dist;
			minIndex = id;
		}
	}
	return minIndex;
}

vector<pair<int,double>> PhipsiLib::neighborPhipsiIndexList(const Phipsi* pp) const
{
	int n = 3;
	int indexList[n];
	double distList[n];
	for(int i=0;i<n;i++){
		indexList[i] = -1;
		distList[i] = 999999.9;
	}

	double d;
	int k;
	for(int i=0;i<pointNum;i++){
		d = pp->distance(*ppList[i]);
		//cout << "index: " << i << " distance: " << d << endl;
		if(d == 0)
			d = 0.001;
        if(d < distList[n-1]) {
            distList[n-1] = d;
            indexList[n-1] = i;
        }
        else
            continue;

        for(int j=n-2;j>=0;j--) {
            if(distList[j+1] < distList[j]) {
                d = distList[j];
                distList[j] = distList[j+1];
                distList[j+1] = d;
                k = indexList[j];
                indexList[j] = indexList[j+1];
                indexList[j+1] = k;
            }
            else
                break;
        }
	}

	vector<pair<int,double>> result;
	for(int i=0;i<n;i++){
		pair<int,double> p(indexList[i],distList[i]);
		result.push_back(p);
	}
	return result;
}

PhipsiLib::~PhipsiLib() {
	for(int i=0;i<this->pointNum;i++)
	{
		delete this->ppList.at(i);
	}
}

ResInfo::ResInfo(int aa, char ss, double sai, int bbIndex){
	this->aaType = aa;
	this->ss = ss;
	this->sai = sai;
	this->bbIndex = bbIndex;
	int saiRep = 0;
	if(sai < 0.173) saiRep = 0;
	else if(sai < 0.258) saiRep = 1;
	else if(sai < 0.295) saiRep = 2;
	else if(sai < 0.331) saiRep = 3;
	else if(sai < 0.375) saiRep = 4;
	else if(sai < 0.425) saiRep = 5;
	else if(sai < 0.475) saiRep = 6;
	else if(sai < 0.525) saiRep = 7;
	else if(sai < 0.575) saiRep = 8;
	else if(sai < 0.625) saiRep = 9;
	else if(sai < 0.675) saiRep = 10;
	else if(sai < 0.725) saiRep = 11;
	else if(sai < 0.775) saiRep = 12;
	else if(sai < 0.825) saiRep = 13;
	else if(sai < 0.875) saiRep = 14;
	else if(sai < 0.925) saiRep = 15;
	else if(sai < 0.975) saiRep = 16;
	else saiRep = 17;

	if(ss == 'H')
	{
		this->ssType = 0;
		this->repIndex = saiRep*1000 + bbIndex;
	}
	else if(ss == 'E') {
		this->ssType = 1;
		this->repIndex = 18000 + saiRep*1000 + bbIndex;
	}
	else {
		this->ssType = 2;
		this->repIndex = 36000 + saiRep*1000 + bbIndex;
	}
}

ResInfo::ResInfo(const string& line){
	vector<string> spt;
	splitString(line, " ", &spt);
	this->aaType = atoi(spt[0].c_str());
	this->ss = spt[1][0];
	this->sai = atof(spt[2].c_str());
	this->bbIndex = atoi(spt[3].c_str());
	int saiRep = 0;
	if(sai < 0.173) saiRep = 0;
	else if(sai < 0.258) saiRep = 1;
	else if(sai < 0.295) saiRep = 2;
	else if(sai < 0.331) saiRep = 3;
	else if(sai < 0.375) saiRep = 4;
	else if(sai < 0.425) saiRep = 5;
	else if(sai < 0.475) saiRep = 6;
	else if(sai < 0.525) saiRep = 7;
	else if(sai < 0.575) saiRep = 8;
	else if(sai < 0.625) saiRep = 9;
	else if(sai < 0.675) saiRep = 10;
	else if(sai < 0.725) saiRep = 11;
	else if(sai < 0.775) saiRep = 12;
	else if(sai < 0.825) saiRep = 13;
	else if(sai < 0.875) saiRep = 14;
	else if(sai < 0.925) saiRep = 15;
	else if(sai < 0.975) saiRep = 16;
	else saiRep = 17;

	if(ss == 'H')
	{
		this->ssType = 0;
		this->repIndex = saiRep*1000 + bbIndex;
	}
	else if(ss == 'E') {
		this->ssType = 1;
		this->repIndex = 18000 + saiRep*1000 + bbIndex;
	}
	else {
		this->ssType = 2;
		this->repIndex = 36000 + saiRep*1000 + bbIndex;
	}
}

string ResInfo::toString(){
	char xx[200];
	sprintf(xx, "%-2d %c %5.3f %4d\n", aaType, ss, sai, bbIndex);
	return string(xx);
}

ResPairInfo::ResPairInfo(int indexA, int indexB, int typeA, int typeB, char ssA, char ssB, double saiA, double saiB, CsMove cm, int sep) {
	// TODO Auto-generated constructor stub
	this->indexA = indexA;
	this->indexB = indexB;
	this->typeAInt = typeA;
	this->typeBInt = typeB;
	this->ssA = ssA;
	this->ssB = ssB;
	this->saiA = saiA;
	this->saiB = saiB;

	if(sep < 0) sep = 0-sep;
	if(sep > 5) sep = 5;

	this->sep = sep;

	this->cm = cm;
	this->dm = DistanceMatrixRes(cm);
	char xx[200];
	sprintf(xx, "%c%c%d", ssA, ssB, sep);
	key = string(xx);
	int ssIDA = 0;
	int ssIDB = 0;
	if(ssA == 'H')
		ssIDA = 0;
	else if(ssA == 'E')
		ssIDA = 1;
	else if(ssA == 'C')
		ssIDA = 2;

	if(ssB == 'H')
		ssIDB = 0;
	else if(ssB == 'E')
		ssIDB = 1;
	else if(ssB == 'C')
		ssIDB = 2;
	keyID = ssIDA*15 + ssIDB * 5 + sep - 1;
	this->bbIndexA = 0;
	this->bbIndexB = 0;
	XYZ localCB = XYZ(-0.942, 0.009, 1.208);
	LocalFrame cs1;
	LocalFrame cs2 = cs1.add(cm);
	XYZ t2 = local2global(cs2, localCB);
	this->cbDistance = localCB.distance(t2);
}

ResPairInfo::ResPairInfo(const string& line){
	vector<string> spt;
	splitString(line, " ", &spt);
	this->indexA = -1;
	this->indexB = -1;
	ResName rn;
	this->typeAInt = rn.triToInt(spt[0]);
	this->typeBInt = rn.triToInt(spt[4]);
	this->ssA = spt[1].at(0);
	this->ssB = spt[5].at(0);
	this->saiA = atof(spt[2].c_str());
	this->saiB = atof(spt[6].c_str());
	this->bbIndexA = atoi(spt[3].c_str());
	this->bbIndexB = atoi(spt[7].c_str());
	string cmString = line.substr(34, 130);
	this->cm = CsMove(cmString);
	this->sep = atoi(spt[20].c_str());
	this->dm = DistanceMatrixRes(cm);
	char xx[200];
	sprintf(xx, "%c%c%d", ssA, ssB, sep);
	key = string(xx);
	int ssIDA = 0;
	int ssIDB = 0;
	if(ssA == 'H')
		ssIDA = 0;
	else if(ssA == 'E')
		ssIDA = 1;
	else if(ssA == 'C')
		ssIDA = 2;

	if(ssB == 'H')
		ssIDB = 0;
	else if(ssB == 'E')
		ssIDB = 1;
	else if(ssB == 'C')
		ssIDB = 2;
	keyID = ssIDA*15 + ssIDB * 5 + sep - 1;
	XYZ localCB = XYZ(-0.942, 0.009, 1.208);
	LocalFrame cs1;
	LocalFrame cs2 = cs1.add(cm);
	XYZ t2 = local2global(cs2, localCB);
	this->cbDistance = localCB.distance(t2);
}

double ResPairInfo::saiDistance(ResPairInfo* other){
	double d1 = this->saiA - other->saiA;
	double d2 = this->saiB - other->saiB;
	return sqrt(d1*d1 + d2*d2);
}

double ResPairInfo::distance(ResPairInfo* other){
	return this->dm.distanceTo(other->dm);
}

string ResPairInfo::toString(){
	ResName rn;
	string typeA = rn.intToTri(typeAInt);
	string typeB = rn.intToTri(typeBInt);
	char xx[200];
	sprintf(xx, "%s %c %5.3f %4d %s %c %5.3f %4d %s %d", typeA.c_str(), ssA, saiA, bbIndexA, typeB.c_str(), ssB, saiB, bbIndexB, cm.toString().c_str(), sep);
	return string(xx);
}

string ResPairInfo::toRevString(){
	ResName rn;
	string typeA = rn.intToTri(typeAInt);
	string typeB = rn.intToTri(typeBInt);
	char xx[200];
	sprintf(xx, "%s %c %5.3f %4d %s %c %5.3f %4d %s %d", typeB.c_str(), ssB, saiB, bbIndexB, typeA.c_str(), ssA, saiA, bbIndexA, cm.reverse().toString().c_str(), sep);
	return string(xx);
}

ResPairInfo::~ResPairInfo() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
