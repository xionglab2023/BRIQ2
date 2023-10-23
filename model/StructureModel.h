/*
 * StructureModel.h
 *
 */

#ifndef MODEL_STRUCTUREMODEL_H_
#define MODEL_STRUCTUREMODEL_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "geometry/Angles.h"
#include "tools/StringTool.h"
#include "model/AtomLib.h"
#include "model/DistanceMatrixRes.h"
#include <string>
#include "stdio.h"
#include <vector>
#include <map>
#include <set>

namespace NSPmodel {
using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPtools;

class Atom{
public:
	string resType;
	string name;
	string type;
	XYZ coord;
	char alt;

	Atom();
	Atom(const string& line);
	Atom(const string& line, int fileType);
	Atom(string name, const XYZ& coord);
	Atom(const Atom& other);  // 拷贝构造函数
	Atom& operator=(const Atom& other);
	bool isBackboneAtom() const;
	void guessAtomType();
	void setAtomType(string type);
	void setResType(string resType);
	bool isIon() const;
	string& getName();
	string& getType();
	XYZ& getCoord();
	string& getResType();
	void setCoord(const XYZ& coord);
	void setCoord(const XYZ&& coord);
	float distance(const Atom& other) const;
	string nameString()	 const;
	virtual ~Atom();
};

class Residue{
private:
	vector<Atom*> atomList;
	vector<Atom*> backboneAtoms;
	vector<Atom*> sidechainAtoms;
	map<string, Atom*> atomMap;
public:
	string resID;
	int resSeqID;
	string triName;
	int intName;
	char chainID;
	int atomNum;
	bool hasLocalFrame;
	bool hasAltConf;
	LocalFrame coordSys;
	char altLoc;
	Residue();
	Residue(string resID, char chainID, string triName);

	void addAtom(Atom* a);
	void setResSeqID(int id);
	void setAltLoc(char c) {this->altLoc = c;}
	bool hasThreeCoreAtoms() const;
	void updateCoordSystem();
	bool hasAtom(const string& atomName) const;
	Atom* getAtom(const string& atomName);
	vector<Atom*>* getAtomList();
	vector<Atom*>* getBackboneAtoms();
	vector<Atom*>* getSidechainAtoms();
	bool contactToNeighbor(Residue* other);
	bool sidechainComplete(AtomLib* atomLib) const;
	XYZ getCbCoord();
	LocalFrame& getCoordSystem();
	char getChainID() const;
	int getResSeqID() const;
	string getResID() const;
	string getType() const;

	int printPDBFormat(ofstream& out, int startAtomID) const;
	virtual ~Residue();
};

class RNABase {
private:
	vector<Atom*> atomList;
	vector<Atom*> backboneAtoms;
	vector<Atom*> sidechainAtoms;
	map<string, Atom*> atomMap;
public:
	string baseID;  // ID number from PDB
	int baseSeqID;
	char baseType;  // AUGCatgc
	int baseTypeInt;  // 0-7, and -1
	string chainID;
	int atomNum;
	bool hasAltConf;
	LocalFrame coordSys;
	bool hasLocalFrame;
	char altLoc;
	RNABase();
	RNABase(const string& baseID, const string& chainID, char baseType);

	/**
	 * @brief 拷贝构造函数.
	 * @note 此函数构造了新的Atom副本，但本类的析构函数为空，
	 * 因此本函数构造的对象在生命期结束前要显式释放其Atom内存。
	*/
	RNABase(const RNABase& other); 

	/**
	 * @brief 移动构造函数.
	 * @note 此函数通过移动获取了Atom内存，但本类的析构函数为空，
	 * 因此本函数构造的对象在生命期结束前要显式释放其Atom内存。
	*/
	RNABase(RNABase&& other) noexcept; 

	void addAtom(Atom* a);
	void setResSeqID(int id) {this->baseSeqID = id;}
	void setAltLoc(char c) {this->altLoc = c;}
	void updateCoordSystem();
	XYZ getCenter(){
		double x = 0;
		double y = 0;
		double z = 0;
		int n = this->sidechainAtoms.size();
		if(n == 0) return XYZ(x,y,z);
		for(int i=0;i<n;i++) {
			x += sidechainAtoms[i]->coord.x_;
			y += sidechainAtoms[i]->coord.y_;
			z += sidechainAtoms[i]->coord.z_;
		}
		return XYZ(x/n, y/n, z/n);
	}
	bool hasAtom(const string& atomName) const {
		return this->atomMap.find(atomName) != atomMap.end();
	}
	Atom* getAtom(const string& atomName) {
		if(this->atomMap.find(atomName) != atomMap.end())
			return this->atomMap[atomName];
		else
			return NULL;
	}
	vector<Atom*>* getAtomList() {return &this->atomList;}
	vector<Atom*>* getBackboneAtoms() {return &this->backboneAtoms;}
	vector<Atom*>* getSidechainAtoms() {return &this->sidechainAtoms;}
	double distanceTo(RNABase& other){
		XYZ a = getCenter();
		XYZ b = other.getCenter();
		return a.distance(b);
	}

	array<XYZ,4> getFourPseudoAtomCoords();

	XYZ getBaseNormVector() {
		XYZ n;
		if(this->baseTypeInt == 0 || this->baseTypeInt == 2){
			Atom* a = getAtom("C1'");
			Atom* b = getAtom("N9");
			Atom* c = getAtom("C4");
			if(a == NULL) return n;
			if(b == NULL) return n;
			if(c == NULL) return n;
			XYZ ba = a->coord - b->coord;
			XYZ bc = c->coord - b->coord;
			XYZ z = ~(bc^ba);
			return z;
		}
		else if(this->baseTypeInt == 1 || this->baseTypeInt == 3){
			Atom* a = getAtom("C1'");
			Atom* b = getAtom("N1");
			Atom* c = getAtom("C2");
			if(a == NULL) return n;
			if(b == NULL) return n;
			if(c == NULL) return n;
			XYZ ba = a->coord - b->coord;
			XYZ bc = c->coord - b->coord;
			XYZ z = ~(bc^ba);
			return z;
		}
		return n;
	}

	double planeAngle(RNABase* other){
		XYZ t1 = getBaseNormVector();
		XYZ t2 = other->getBaseNormVector();
		double ang = NSPgeometry::angleX(t1,t2);
		if(ang > 90)
			ang = 180 - 90;
		return ang;
	}

	double planeDistance(RNABase* other){
		LocalFrame cs1 = getCoordSystem();
		LocalFrame cs2 = other->getCoordSystem();
		XYZ center1 = getCenter();
		XYZ center2 = other->getCenter();
		XYZ localCenter1 = global2local(cs2, center1);
		XYZ localCenter2 = global2local(cs1, center2);
		double d1 = abs(localCenter1.z_);
		double d2 = abs(localCenter2.z_);
		if(d1 < d2) return d1;
		else return d2;
	}

	bool sidechainComplete(AtomLib* atLib) const;

	bool backboneComplete() {
		if(getAtom("C1'") == NULL) return false;
		if(getAtom("C2'") == NULL) return false;
		if(getAtom("C3'") == NULL) return false;
		if(getAtom("C4'") == NULL) return false;
		if(getAtom("O4'") == NULL) return false;
		if(getAtom("O3'") == NULL) return false;
		if(getAtom("C5'") == NULL) return false;

		if(baseTypeInt<4 && getAtom("O2'") == NULL) return false;

		return true;
	}

	LocalFrame getCoordSystem() {
		if(!hasLocalFrame) {
			updateCoordSystem();
		}
		return this->coordSys;
	}

	string getChainID() const {return this->chainID;}
	int getResSeqID() const {return this->baseSeqID;}
	string getResID() const {return this->baseID;}
	char getType() const {return this->baseType;}
	bool connectToNeighbor(RNABase* other) {
		Atom* a = getAtom("O3'");
		Atom* b = other->getAtom("P");
		if(a == NULL || b == NULL)
			return false;
		else if(a->distance(*b) < 2.0)
			return true;
		else
			return false;
	}

	bool contactTo(RNABase* other){
		for(int i=0;i<this->sidechainAtoms.size();i++){
			Atom* a = this->sidechainAtoms[i];
			for(int j=0;j<other->sidechainAtoms.size();j++){
				Atom* b = other->sidechainAtoms[j];
				double d = a->coord.distance(b->coord);
				if(d < 4.5)
					return true;
			}
		}
		return false;
	}

	int printPDBFormat(ofstream& out, int startAtomID) const;
	virtual ~RNABase();
};

class PolarAtom{
private:
	string uniqueName;
	XYZ support;
	XYZ core;
	bool isDonor;
	bool isAcceptor;
	float vdwRadius;
public:
	PolarAtom(Residue* res, string atomName);
	PolarAtom(RNABase* base, string atomName);
	PolarAtom(Residue* res, string atomName, Atom* preC);

	XYZ& getCore() {return this->core;}
	XYZ& getSupport() {return this->support;}
	bool hbondedTo(PolarAtom* other){
		double d = this->core.distance(other->core);
		if(d < 1.5) return false;
		if(d > this->vdwRadius + other->vdwRadius - 0.1) return false;
		double ang1 = angleX(this->support, this->core, other->core);
		double ang2 = angleX(this->core, other->core, other->support);
		if(ang1 < 100) return false;
		if(ang2 < 100) return false;

		return true;
	}

	bool isEmpty() {
		if(this->core.length() == 0 && this->core.distance(this->support) == 0) return true;
		else return false;
	}
	bool isDonerAtom() {return this->isDonor;}
	bool isAcceptorAtom() {return this->isAcceptor;}
	bool neighborTo(PolarAtom* other) {return this->core.distance(other->core) < (this->vdwRadius + other->vdwRadius + 2.0);}
	string getName() {return this->uniqueName;}
	float getVdwRadius() {return this->vdwRadius;}

	virtual ~PolarAtom();
};


class ProteinChain{
private:
	string pdbID;
	char chainID;
	int chainLen;
	vector<Residue*> resList;
	map<string,Residue*> resMap;

public:
	ProteinChain();
	ProteinChain(string pdbID, char chainID);
	ProteinChain(char chainID);

	void setPDBID(string pdbID);
	void setChainID(char c);
	string getPDBID() const;
	char getChainID() const;
	int getChainLength() const;
	vector<Residue*>& getResList();
	Residue* getResidue(const string& resID);
	void addResidue(Residue* res);
	string getSequence() const;
	int printPDBFormat(ofstream& out, int startAtomID) const;
	int printPDBFormatNoHydrogen(ofstream& out, int startAtomID) const;
	virtual ~ProteinChain();
};

class RNAChain{
private:
	string pdbID;
	string chainID;
	int chainLen;
	vector<RNABase*> baseList;
	map<string, RNABase*> baseMap;
public:
	RNAChain();
	RNAChain(const string& pdbID, const string& chainID);
	RNAChain(const string& chainID);

	/**
	 * @brief 拷贝构造函数.
	 * @note 此函数构造了新的 RNABase 和 Atom 副本，但本类的析构函数为空，
	 * 因此本函数构造的对象在生命期结束前要显式释放其 RNABase 和 Atom 内存。
	*/
	RNAChain(const RNAChain& other);

	/**
	 * @brief 移动构造函数.
	 * @note 此函数移动获取了 RNABase 和 Atom 内存，但本类的析构函数为空，
	 * 因此本函数构造的对象在生命期结束前要显式释放其 RNABase 和 Atom 内存。
	*/
	RNAChain(RNAChain&& other) noexcept;

	void setPDBID(const string& pdbID) {
		this->pdbID = pdbID;
	}
	void setChainID(const string& c) {this->chainID = c;}
	string getPDBID() {return this->pdbID;}
	string getChainID() {return this->chainID;}
	int getChainLength() {return this->chainLen;}
	vector<RNABase*>& getBaseList() {return this->baseList;}
	RNABase* getBase(const string& baseID) {
		if(baseMap.find(baseID) != baseMap.end())
			return baseMap[baseID];
		return NULL;
	}
	void addBase(RNABase* base) {
		this->chainLen++;
		this->baseList.push_back(base);
		this->baseMap[base->baseID] = base;
	}
	string getSequence() {
		char s[baseList.size() +1];
		for(int i=0;i<baseList.size();i++) {
			s[i] = baseList[i]->baseType;
		}
		s[baseList.size()] = '\0';
		return string(s);
	}
	int printPDBFormat(ofstream& out, int startAtomID) const;

};

class PDB{
private:
	string pdbID;
	vector<ProteinChain*> chains;
	vector<Residue*> residues;
public:
	PDB();
	PDB(const string& pdbFile, const string& pdbID);
	PDB& operator=(const PDB& other);
	vector<ProteinChain*>& getChains();
	ProteinChain* getFirstChain();
	ProteinChain* getChain(char c);
	string getFirstSeq();
	vector<Residue*>& getResList();
	vector<Residue*> getValidResList(){
		vector<Residue*> list;
		for(int i=0;i<residues.size();i++){
			Residue* res = residues.at(i);
			if(res->hasThreeCoreAtoms() && res->hasAtom("O"))
				list.push_back(res);
		}
		return list;
	}
	void printPDBFormat(ofstream& out) const;
	void printPDBFormatNoHydrogen(ofstream& out) const;
	string getPDBID()
	{
		return this->pdbID;
	}
	virtual ~PDB();
};

class RNAPDB{
private:
	string pdbID;
	vector<RNAChain*> chains;
	vector<RNABase*> baseList;
public:
	RNAPDB();
	RNAPDB(const string& pdbFile);
	RNAPDB(const string& pdbFile, const string& pdbID);

	void readPDB(const string& pdbFile);
	void readCIF(const string& cifFile);

	void addChain(RNAChain* rc){
		this->chains.push_back(rc);
		for(RNABase* b : rc->getBaseList()){
			this->baseList.push_back(b);
		}
	}

	vector<RNAChain*>& getChains() {return this->chains;}
	RNAChain* getFirstChain() {return this->chains[0];}
	RNAChain* getChain(const string& c) {
		for(int i=0;i<this->chains.size();i++) {
			if(this->chains[i]->getChainID() == c)
				return this->chains[i];
		}
		return NULL;
	}
	vector<RNABase*>& getBaseList() {return this->baseList;}
	vector<RNABase*> getValidBaseList(AtomLib* atLib) {
		vector<RNABase*> list;
		for(int i=0;i<baseList.size();i++){
			if(baseList[i]->sidechainComplete(atLib))
				list.push_back(baseList[i]);
		}
		return list;
	}
	void printPDBFormat(ofstream& out) const;
	virtual ~RNAPDB();
};

class Phipsi{
public:
	float phi;
	float psi;
	Phipsi() {this->phi = 0; this->psi = 0;}
	Phipsi(float phi, float psi){
		// TODO Auto-generated constructor stub
		this->phi = phi;
		this->psi = psi;
		if(phi == 360)
			this->phi = -80;
		if(psi == 360)
			this->psi = 130;

		if(this->phi >= 180)
			this->phi = this->phi - 360;
		if(this->phi < -180)
			this->phi = this->phi + 360;
		if(this->psi >= 180)
			this->psi = this->psi - 360;
		if(this->psi < -180)
			this->psi = this->psi + 360;

		if(this->phi > 180 || this->phi < -180 || this->psi > 180 || this->psi < -180)
		{
			cerr << "invalid phi psi: " << phi << " " << psi << endl;
			exit(1);
		}
	}
	float distance(const Phipsi& other) const;
	char regionAB() const;
	virtual ~Phipsi();
};

class PhipsiLib{
private:
	int pointNum;
	vector<Phipsi*> ppList;
	int indexTable[36][36][20];
	void creatIndexTable();
public:
	PhipsiLib();
	Phipsi* indexToPhipsi(int id) const;
	int phipsiToIndex(const Phipsi* pp) const;
	vector<pair<int,double>> neighborPhipsiIndexList(const Phipsi* pp) const;
	int findNearestPointsWithoutIndex(Phipsi* pp);
	virtual ~PhipsiLib();
};


class ResInfo{
public:
	int aaType;
	char ss;
	int ssType;
	double sai;
	int bbIndex; //index from 1k
	int repIndex;

	ResInfo(const string& line);
	ResInfo(int aa, char ss, double sai, int bbIndex);
	string toString();
};


class ResPairInfo {
public:
	int indexA;
	int indexB;
	int typeAInt;
	int typeBInt;
	char ssA;
	char ssB;
	double saiA;
	double saiB;
	double cbDistance;
	int bbIndexA;
	int bbIndexB;
	int sep;
	CsMove cm;
	DistanceMatrixRes dm;
	string key;
	int keyID;

	ResPairInfo(int indexA, int indexB, int typeA, int typeB, char ssA, char ssB, double saiA, double saiB, CsMove cm, int sep);
	ResPairInfo(const string& line);
	double saiDistance(ResPairInfo* other);
	double distance(ResPairInfo* other);
	string toString();
	string toRevString();

	virtual ~ResPairInfo();
};

class SaiPair{
public:
	float saiA;
	float saiB;

	SaiPair(float x, float y){
		this->saiA = x;
		this->saiB = y;
	}

	float distanceSquare(SaiPair* other){
		float dx = this->saiA - other->saiA;
		float dy = this->saiB - other->saiB;
		return dx*dx+dy*dy;
	}
};

} /* namespace NSPmodel*/




#endif /* MODEL_STRUCTUREMODEL_H_ */
