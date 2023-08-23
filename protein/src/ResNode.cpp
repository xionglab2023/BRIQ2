/*
 * ResNode.cpp
 *
 */

#include "protein/ResNode.h"

namespace NSPprotein{

ResNode::ResNode(){
	this->aaType = 0;
	this->nodeID = 0;
	this->seqID = 0;
	this->sai = 0;
	this->ss = 'H';
	this->ssInt = 0;

	this->conf = NULL;
	this->confTmp = NULL;

	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->hbondAChild = NULL;
	this->hbondBChild = NULL;
	this->jumpChild = NULL;
	this->upConnection = NULL;
	scRotmerNeedDelete = false;
}


ResNode::ResNode(const string & line, int nodeID, ResName* rn, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib){

	this->nodeID = nodeID;
	vector<string> spt;
	splitString(line, " ", &spt);
	LocalFrame cs;
	ResBBRotamer* bbRot;
	ResScRotamer* scRot;
	if(spt.size() == 16) {
		this->aaType = rn->triToInt(line.substr(0, 3));
		if(aaType < 0 || aaType > 19){
			cout << "invalid aaType: " << line.substr(0,3) << endl;
			exit(1);
		}
		this->seqID = atoi(spt[1].c_str());
		int bbRotID = atoi(spt[2].c_str());
		bbRot = bbLib->allRotLib1w[aaType][bbRotID];
		int scRotID = atoi(spt[3].c_str());
		if(scRotID == -1)
			scRotID = 0;
		scRot = scLib->rotList[aaType][scRotID];
		double mtx[3][3];
		mtx[0][0] = atof(spt[4].c_str());
		mtx[0][1] = atof(spt[5].c_str());
		mtx[0][2] = atof(spt[6].c_str());
		mtx[1][0] = atof(spt[7].c_str());
		mtx[1][1] = atof(spt[8].c_str());
		mtx[1][2] = atof(spt[9].c_str());
		mtx[2][0] = atof(spt[10].c_str());
		mtx[2][1] = atof(spt[11].c_str());
		mtx[2][2] = atof(spt[12].c_str());
		XYZ origin = XYZ(atof(spt[13].c_str()), atof(spt[14].c_str()), atof(spt[15].c_str()));
		TransMatrix tm(mtx);
		cs = LocalFrame(origin, tm);
		this->ss = 'C';
		this->ssInt = 2;
		this->sai = 0.5;
	}
	else if(spt.size() == 18) {
		this->aaType = rn->triToInt(line.substr(0, 3));
		if(aaType < 0 || aaType > 19){
			cout << "invalid aaType: " << line.substr(0,3) << endl;
			exit(1);
		}
		this->seqID = atoi(spt[1].c_str());
		this->ss = spt[2][0];
		if(ss=='H')
			ssInt = 0;
		else if(ss == 'E')
			ssInt = 1;
		else
			ssInt = 2;
		this->sai = atof(spt[3].c_str());
		int bbRotID = atoi(spt[4].c_str());
		bbRot = bbLib->allRotLib1w[aaType][bbRotID];
		int scRotID = atoi(spt[5].c_str());
		if(scRotID == -1)
			scRotID = 0;
		scRot = scLib->rotList[aaType][scRotID];
		double mtx[3][3];
		mtx[0][0] = atof(spt[6].c_str());
		mtx[0][1] = atof(spt[7].c_str());
		mtx[0][2] = atof(spt[8].c_str());
		mtx[1][0] = atof(spt[9].c_str());
		mtx[1][1] = atof(spt[10].c_str());
		mtx[1][2] = atof(spt[11].c_str());
		mtx[2][0] = atof(spt[12].c_str());
		mtx[2][1] = atof(spt[13].c_str());
		mtx[2][2] = atof(spt[14].c_str());
		XYZ origin = XYZ(atof(spt[15].c_str()), atof(spt[16].c_str()), atof(spt[17].c_str()));
		TransMatrix tm(mtx);
		cs = LocalFrame(origin, tm);
	}
	else {
		cout << "invalid resNode line: " << line << endl;
		exit(0);
	}

	this->cs2 = cs;
	this->cs2Tmp = cs;
	this->cs1 = cs2 + bbRot->cm21;
	this->cs1Tmp = this->cs1;
	this->cs3 = cs2 + bbRot->cm23;
	this->cs3Tmp = this->cs3;

	this->conf = new ResConformer(aaType);
	this->conf->init(scRot, bbRot, cs2);
	this->confTmp = new ResConformer(aaType);
	this->confTmp->init(scRot, bbRot, cs2);

	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->hbondAChild = NULL;
	this->hbondBChild = NULL;
	this->jumpChild = NULL;

	this->upConnection = NULL;
	scRotmerNeedDelete = false;
}

ResNode::ResNode(Residue* resP, Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib){
	this->aaType = res->intName;
	this->seqID = seqID;
	this->nodeID = nodeID;
	this->ss = 'C';
	this->ssInt = 2;
	this->sai = 0.5;
	scRotmerNeedDelete = false;

	ResBBRotamer* bbRot = new ResBBRotamer(resP, res, atLib);
	int rotID = bbLib->getRotamerIndex1W(bbRot);
	delete bbRot;

	int scRotID = 0;

	this->conf = new ResConformer(aaType);
	this->confTmp = new ResConformer(aaType);

	this->cs2 = res->getCoordSystem();
	this->cs2Tmp = cs2;
	this->cs1 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm21;
	this->cs1Tmp = this->cs1;
	this->cs3 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm23;
	this->cs3Tmp = this->cs3;

	if(res->sidechainComplete(atLib))
	{
		ResScRotamer* scRot = new ResScRotamer(res, atLib);
		this->conf->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		scRotmerNeedDelete = true;
	}
	else {
		this->conf->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
	}



	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->hbondAChild = NULL;
	this->hbondBChild = NULL;
	this->jumpChild = NULL;
	this->upConnection = NULL;
}

ResNode::ResNode(Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib){
	this->aaType = res->intName;
	this->seqID = seqID;
	this->nodeID = nodeID;
	this->ss = 'C';
	this->ssInt = 2;
	this->sai = 0.5;
	scRotmerNeedDelete = false;
	ResBBRotamer* bbRot = new ResBBRotamer(res, atLib);
	int rotID = bbLib->getRotamerIndex1W(bbRot);
	delete bbRot;


	this->cs2 = res->getCoordSystem();
	this->cs2Tmp = cs2;
	this->cs1 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm21;
	this->cs1Tmp = this->cs1;
	this->cs3 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm23;
	this->cs3Tmp = this->cs3;

	int scRotID = 0;
	this->conf = new ResConformer(aaType);
	this->confTmp = new ResConformer(aaType);

	if(res->sidechainComplete(atLib))
	{

		ResScRotamer* scRot = new ResScRotamer(res, atLib);
		this->conf->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		scRotmerNeedDelete = true;
		//scRotID = scLib->getRotamerID(scRot);
		//delete scRot;
	}
	else {
		this->conf->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
	}




	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->hbondAChild = NULL;
	this->hbondBChild = NULL;
	this->jumpChild = NULL;
	this->upConnection = NULL;
}

ResNode::ResNode(Residue* resP, Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib){
	this->aaType = res->intName;
	this->seqID = seqID;
	this->nodeID = nodeID;
	this->ss = 'C';
	this->ssInt = 2;
	this->sai = 0.5;
	scRotmerNeedDelete = false;
	ResBBRotamer* bbRot = new ResBBRotamer(resP, res, atLib);
	int rotID = bbLib->getRotamerIndex1W(bbRot);
	delete bbRot;

	int scRotID = 0;

	this->conf = new ResConformer(aaType);
	this->confTmp = new ResConformer(aaType);

	this->cs2 = res->getCoordSystem();
	this->cs2Tmp = cs2;
	this->cs1 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm21;
	this->cs1Tmp = this->cs1;
	this->cs3 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm23;
	this->cs3Tmp = this->cs3;

	if(res->sidechainComplete(atLib))
	{
		ResScRotamer* scRot = new ResScRotamer(res, atLib);
		this->conf->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		scRotmerNeedDelete = true;
	}
	else {
		this->conf->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
	}



	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->hbondAChild = NULL;
	this->hbondBChild = NULL;
	this->jumpChild = NULL;
	this->upConnection = NULL;
}

ResNode::ResNode(Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib){

	this->aaType = res->intName;
	this->seqID = seqID;
	this->nodeID = nodeID;
	this->ss = 'C';
	this->ssInt = 2;
	this->sai = 0.5;
	scRotmerNeedDelete = false;

	ResBBRotamer* bbRot = new ResBBRotamer(res, atLib);
	int rotID = bbLib->getRotamerIndex1W(bbRot);
	delete bbRot;


	this->cs2 = res->getCoordSystem();
	this->cs2Tmp = cs2;
	this->cs1 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm21;
	this->cs1Tmp = this->cs1;
	this->cs3 = cs2 + bbLib->allRotLib1w[aaType][rotID]->cm23;
	this->cs3Tmp = this->cs3;


	int scRotID = 0;
	this->conf = new ResConformer(aaType);
	this->confTmp = new ResConformer(aaType);

	if(res->sidechainComplete(atLib)){
		ResScRotamer* scRot = new ResScRotamer(res, atLib);
		this->conf->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scRot, bbLib->allRotLib1w[aaType][rotID], cs2);
		scRotmerNeedDelete = true;
	}
	else {
		this->conf->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
		this->confTmp->init(scLib->rotList[aaType][scRotID], bbLib->allRotLib1w[aaType][rotID], cs2);
	}


	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->hbondAChild = NULL;
	this->hbondBChild = NULL;
	this->jumpChild = NULL;
	this->upConnection = NULL;
}

ResNode::~ResNode() {

	delete this->conf;
	delete this->confTmp;
}
}


