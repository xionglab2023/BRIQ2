/*
 * BRNode.cpp
 *
 */

#include "predNA/BRNode.h"

namespace NSPpredna {
BRNode& BRNode::operator =(const BRNode& other){
	this->baseType = other.baseType;
	this->seqID = other.seqID;

	this->baseConf = new BaseConformer(other.baseConf->rot, other.baseConf->cs1);
	this->baseConfTmp = new BaseConformer(other.baseConfTmp->rot, other.baseConfTmp->cs1);

	this->riboseConf = new RiboseConformer(other.riboseConf->rot, other.riboseConf->cs1);
	this->riboseConfTmp = new RiboseConformer(other.riboseConfTmp->rot, other.riboseConfTmp->cs1);

	LocalFrame cs2 = riboseConf->cs2;

	this->phoConf = new PhosphateConformer(other.phoConf->rot, other.riboseConf->cs2);
	this->phoConfTmp = new PhosphateConformer(other.phoConfTmp->rot, other.riboseConfTmp->cs2);

	this->fixed = other.fixed;
	this->father = other.father;
	this->leftChild = other.leftChild;
	this->midChild = other.midChild;
	this->rightChild = other.rightChild;
	this->reverseChild = other.reverseChild;
	this->upConnection = other.upConnection;
	this->connectToNeighbor = other.connectToNeighbor;

	return *this;
}

void BRNode::copyValueFrom(const BRNode& other) {
	this->baseType = other.baseType;
	this->seqID = other.seqID;
	this->baseConf->copyValueFrom(other.baseConf);
	this->baseConfTmp->copyValueFrom(other.baseConfTmp);
	this->riboseConf->copyValueFrom(other.riboseConf);
	this->riboseConfTmp->copyValueFrom(other.riboseConfTmp);
	this->phoConf->copyValueFrom(other.phoConf);
	this->phoConfTmp->copyValueFrom(other.phoConfTmp);

	this->fixed = other.fixed;
	this->father = other.father;
	this->leftChild = other.leftChild;
	this->midChild = other.midChild;
	this->rightChild = other.rightChild;
	this->reverseChild = other.reverseChild;
	this->upConnection = other.upConnection;
	this->connectToNeighbor = other.connectToNeighbor;
}


vector<Atom*> BRNode::toAtomList(AtomLib& atLib) {

    vector<Atom*> list;
    vector<string>* names = atLib.getRnaSidechainAtoms(this->baseType);
    vector<XYZ> tList;
    for(int i=0;i<baseConf->rot->atomNum;i++){
    	tList.push_back(baseConf->coords[i]);
    }


    names->push_back("C1'");
    names->push_back("C2'");
    names->push_back("C3'");
    names->push_back("C4'");
    names->push_back("O4'");
    names->push_back("O3'");
    names->push_back("C5'");
    if(riboseConf->hasO2)
    	names->push_back("O2'");

    int n = 7;
    if(riboseConf->hasO2)
    	n = 8;

    for(int i=0;i<n;i++){
    	tList.push_back(riboseConf->coords[i]);
    }


    vector<Atom*> atomList;
    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names->at(i), tList[i]));
    }

    return atomList;
}

vector<Atom*> BRNode::toBaseAtomList(AtomLib& atLib) {


    vector<string>* names = atLib.getRnaSidechainAtoms(this->baseType);
    vector<XYZ> tList;
    for(int i=0;i<baseConf->rot->atomNum;i++){
    	tList.push_back(baseConf->coords[i]);
    }
    vector<Atom*> atomList;
    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names->at(i), tList[i]));
    }
    return atomList;
}

vector<Atom*> BRNode::phoAtoms(){
	 vector<Atom*> atomList;
	 if(connectToNeighbor) {
	 	atomList.push_back(new Atom("P", phoConf->coords[0]));
	 	atomList.push_back(new Atom("O5'", phoConf->coords[1]));
	 	atomList.push_back(new Atom("OP1", phoConf->coords[2]));
	 	atomList.push_back(new Atom("OP2", phoConf->coords[3]));
	 }
	 return atomList;
}

vector<Atom*> BRNode::toTmpAtomList(AtomLib& atLib) {

    vector<Atom*> list;
    vector<string>* names = atLib.getRnaSidechainAtoms(this->baseType);
    vector<XYZ> tList;
    for(int i=0;i<baseConf->rot->atomNum;i++){
    	tList.push_back(baseConf->coords[i]);
    }


    names->push_back("C1'");
    names->push_back("C2'");
    names->push_back("C3'");
    names->push_back("C4'");
    names->push_back("O4'");
    names->push_back("O3'");
    names->push_back("C5'");
    if(riboseConf->hasO2)
    	names->push_back("O2'");


    int n = 7;
    if(riboseConf->hasO2)
    	n = 8;

    for(int i=0;i<n;i++){
    	tList.push_back(riboseConf->coords[i]);
    }


    vector<Atom*> atomList;
    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names->at(i), tList[i]));
    }

    if(connectToNeighbor) {
    	atomList.push_back(new Atom("P", phoConfTmp->coords[0]));
    	atomList.push_back(new Atom("O5", phoConfTmp->coords[1]));
    	atomList.push_back(new Atom("OP1", phoConfTmp->coords[2]));
    	atomList.push_back(new Atom("OP2", phoConfTmp->coords[3]));
    }
    return atomList;
}

bool BRNode::baseConsistent(){
	for(int i=0;i<baseConf->rot->atomNum;i++){
		if(squareDistance(this->baseConf->coords[i], this->baseConfTmp->coords[i]) > 0.0000001) return false;
	}
	return true;
}

bool BRNode::riboConsistent(){
	for(int i=0;i<8;i++){
		if(squareDistance(this->riboseConf->coords[i], this->riboseConfTmp->coords[i]) > 0.0000001) return false;
	}
	return true;
}

bool BRNode::rotamerConsistent(){
	for(int i=0;i<11;i++){
		if(this->riboseConf->rot->distanceTo(this->riboseConfTmp->rot) > 0.0000001) return false;
	}
	return true;
}

bool BRNode::consistent(){
	return baseConsistent() && riboConsistent();
}

bool BRNode::phoConsistent(){
	if(this->phoConf->equalTo(this->phoConfTmp)) return true;
	return false;
}

bool BRNode::phoLocalConsistent(){
	if(this->phoConf->rot->equalTo(this->phoConfTmp->rot)) return true;
	return false;
}

void BRNode::checkRotamer(){
	LocalFrame xcs2 = LocalFrame(local2global(baseConf->cs1, riboseConf->rot->localCoords[1]), local2global(baseConf->cs1, riboseConf->rot->localCoords[2]), local2global(baseConf->cs1,riboseConf->rot->localCoords[5]));
	if(!xcs2.equalTo(riboseConf->cs2)) {
		cout << "cs2 error" << endl;
		xcs2.print();
		cout << endl;
		riboseConf->cs2.print();
	}
	cout << endl;
}


void BRNode::checkTmpRotamer(){
	LocalFrame xcs2 = LocalFrame(local2global(baseConfTmp->cs1, riboseConfTmp->rot->localCoords[1]), local2global(baseConfTmp->cs1, riboseConfTmp->rot->localCoords[2]), local2global(baseConfTmp->cs1,riboseConfTmp->rot->localCoords[5]));
	if(!xcs2.equalTo(riboseConfTmp->cs2)) {
		cout << "cs2 error" << endl;
		xcs2.print();
		cout << endl;
		riboseConfTmp->cs2.print();
	}
	cout << endl;
}

BRNode::~BRNode() {

	delete this->baseConf;
	delete this->baseConfTmp;
	delete this->riboseConf;
	delete this->riboseConfTmp;
	delete this->phoConf;
	delete this->phoConfTmp;
}

} /* namespace NSPpred */
