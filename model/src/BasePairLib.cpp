/*
 * BasePairLib.cpp
 *
 *  Created on: 2023��8��7��
 *      Author: nuc
 */

#include "model/BasePairLib.h"

namespace NSPmodel {

BasePairLib::BasePairLib() {

	string path = NSPdataio::datapath();
	string nbCenter = path+"basePair/nb.center.dm";
	string nbCenterInfo = path+"basePair/nb.center.info";

	string nnbCenter = path+"basePair/nnb.center.dm";
	string nnbCenterInfo = path+"basePair/nnb.center.info";

	ifstream file;
	string s;
	vector<string> spt;
	int typeA, typeB, currentType, lastType, currentIndex;


	map<char, int> typeMap;
	typeMap['A'] = 0;
	typeMap['U'] = 1;
	typeMap['G'] = 2;
	typeMap['C'] = 3;


	for(int i=0;i<16;i++){
		this->nbBasePairNum[i] = 0;
		this->nnbBasePairNum[i] = 0;
		for(int j=0;j<200;j++){
			this->nbEnegy[i][j] = 0.0;
			this->nbProportion[i][j] = 0.0;
			this->nnbEnegy[i][j] = 0.0;
			this->nnbProportion[i][j] = 0.0;
		}
	}

	file.open(nbCenter.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open file: nb.center.dm" << endl;
		exit(1);
	}
	currentType = -1;
	lastType = -1;
	currentIndex = -1;
	while(getline(file, s)){
		if(s.length() < 10) continue;
		BaseDistanceMatrix dm(s);
		typeA = typeMap[s[0]];
		typeB = typeMap[s[2]];
		currentType = typeA*4+typeB;
		if(currentType != lastType)
			currentIndex = 0;
		else
			currentIndex++;
		this->nbDMClusterCenters[typeA*4+typeB][currentIndex] = dm;
		this->nbBasePairNum[currentType]++;
		lastType = currentType;
	}
	file.close();


	file.open(nbCenterInfo.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open file: nb.center.info" << endl;
		exit(1);
	}
	currentType = -1;
	lastType = -1;
	currentIndex = -1;
	while(getline(file, s)){
		if(s.length() < 10) continue;
		typeA = typeMap[s[0]];
		typeB = typeMap[s[1]];
		currentType = typeA*4+typeB;
		if(currentType != lastType)
			currentIndex = 0;
		else
			currentIndex++;
		splitString(s, " ", &spt);
		this->nbEnegy[typeA*4+typeB][currentIndex] = atof(spt[1].c_str());
		this->nbProportion[typeA*4+typeB][currentIndex] = atof(spt[2].c_str());
		lastType = currentType;
	}
	file.close();


	file.open(nnbCenter.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open file: nnb.center.dm" << endl;
		exit(1);
	}
	currentType = -1;
	lastType = -1;
	currentIndex = -1;
	while(getline(file, s)){
		if(s.length() < 10) continue;
		BaseDistanceMatrix dm(s);
		typeA = typeMap[s[0]];
		typeB = typeMap[s[2]];
		currentType = typeA*4+typeB;
		if(currentType != lastType)
			currentIndex = 0;
		else
			currentIndex++;
		this->nnbDMClusterCenters[typeA*4+typeB][currentIndex] = dm;
		this->nnbBasePairNum[currentType]++;
		lastType = currentType;
	}
	file.close();



	file.open(nnbCenterInfo.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open file: nnb.center.info" << endl;
		exit(1);
	}
	currentType = -1;
	lastType = -1;
	currentIndex = -1;
	while(getline(file, s)){
		if(s.length() < 10) continue;
		typeA = typeMap[s[0]];
		typeB = typeMap[s[1]];
		currentType = typeA*4+typeB;
		if(currentType != lastType)
			currentIndex = 0;
		else
			currentIndex++;
		splitString(s, " ", &spt);
		this->nnbEnegy[typeA*4+typeB][currentIndex] = atof(spt[1].c_str());
		this->nnbProportion[typeA*4+typeB][currentIndex] = atof(spt[2].c_str());

		lastType = currentType;
	}
	file.close();
}

int BasePairLib::getPairType(BaseDistanceMatrix dm, int typeA, int typeB, int sep){

	if(typeA > 3)
		typeA = typeA - 4;

	if(typeB > 3)
		typeB = typeB - 4;

	if(typeA < 0 || typeA > 3) {
		cout << "invalid base type: " << typeA << endl;
		return -1;
	}

	if(typeB < 0 || typeB > 3) {
		cout << "invalid base type: " << typeB << endl;
		return -1;
	}

	double minD = 999.9;
	int minIndex = -1;
	int pairType = typeA*4+typeB;
	double d;

	if(sep == 1) { //neighbor base pair
		int pairNum = nbBasePairNum[pairType];
		for(int i=0;i<pairNum;i++){
			d = nbDMClusterCenters[pairType][i].distanceTo(dm);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
	}
	else { //non-neighbor base pair
		int pairNum = nnbBasePairNum[pairType];
		for(int i=0;i<pairNum;i++){
			d = nnbDMClusterCenters[pairType][i].distanceTo(dm);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
	}

	if(minD < 1.5)
		return minIndex;
	else
		return -1;
}

double BasePairLib::distanceToClusterCenter(BaseDistanceMatrix dm, int typeA, int typeB, int sep){

	if(typeA > 3)
		typeA = typeA - 4;

	if(typeB > 3)
		typeB = typeB - 4;

	if(typeA < 0 || typeA > 3) {
		cout << "invalid base type: " << typeA << endl;
		return -1;
	}

	if(typeB < 0 || typeB > 3) {
		cout << "invalid base type: " << typeB << endl;
		return -1;
	}

	double minD = 999.9;
	int minIndex = -1;
	int pairType = typeA*4+typeB;
	double d;

	if(sep == 1) { //neighbor base pair
		int pairNum = nbBasePairNum[pairType];
		for(int i=0;i<pairNum;i++){
			d = nbDMClusterCenters[pairType][i].distanceTo(dm);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
	}
	else { //non-neighbor base pair
		int pairNum = nnbBasePairNum[pairType];
		for(int i=0;i<pairNum;i++){
			d = nnbDMClusterCenters[pairType][i].distanceTo(dm);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
	}
	return minD;
}

double BasePairLib::getEnergy(BaseDistanceMatrix dm, int typeA, int typeB, int sep){

	if(typeA > 3)
		typeA = typeA - 4;

	if(typeB > 3)
		typeB = typeB - 4;

	if(typeA < 0 || typeA > 3) {
		cout << "invalid base type: " << typeA << endl;
		return -1;
	}

	if(typeB < 0 || typeB > 3) {
		cout << "invalid base type: " << typeB << endl;
		return -1;
	}

	double minD = 999.9;
	int minIndex = -1;
	int pairType = typeA*4+typeB;
	double d;

	if(sep == 1) { //neighbor base pair
		int pairNum = nbBasePairNum[pairType];
		for(int i=0;i<pairNum;i++){
			d = nbDMClusterCenters[pairType][i].distanceTo(dm);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
	}
	else { //non-neighbor base pair
		int pairNum = nnbBasePairNum[pairType];
		for(int i=0;i<pairNum;i++){
			d = nnbDMClusterCenters[pairType][i].distanceTo(dm);
			if(d < minD){
				minD = d;
				minIndex = i;
			}
		}
	}

	double ene = 0.0;



	if(sep == 1) {
		double e = this->nbEnegy[typeA*4+typeB][minIndex];
		double p = this->nbProportion[typeA*4+typeB][minIndex];
		if(p < 0.001) p = 0.001;
		if(p > 0.3) p = 0.3;
		ene = e*0.25 - (log(p)+6.0)*0.4;
		if(ene > -0.01)
			ene = -0.01;
	}
	else if(minD < 1.5) {
		ene = (this->nnbEnegy[typeA*4+typeB][minIndex]-2.0)*0.25;
	}

	if(ene > 0)
		ene = 0;

	return ene;
}

double BasePairLib::getPairEnergy(RNABase* baseA, RNABase* baseB){
	int sep = 2;
	if(baseA->connectToNeighbor(baseB))
		sep = 1;

	double ene = 0.0;
	LocalFrame csA = baseA->getCoordSystem();
	LocalFrame csB = baseB->getCoordSystem();

	BaseDistanceMatrix dm(csA, csB);

	int pairType = getPairType(dm, baseA->baseTypeInt, baseB->baseTypeInt, sep);
	int typeA = baseA->baseTypeInt;
	int typeB = baseB->baseTypeInt;

	ene = getEnergy(dm, typeA, typeB, sep);
	if(sep == 1)
		return ene;

	int hbondNum = 0;
	PolarAtom* o2A = new PolarAtom(baseA, "O2'");
	PolarAtom* op1A = new PolarAtom(baseA, "OP1");
	PolarAtom* op2A = new PolarAtom(baseA, "OP2");

	PolarAtom* o2B = new PolarAtom(baseB, "O2'");
	PolarAtom* op1B = new PolarAtom(baseB, "OP1");
	PolarAtom* op2B = new PolarAtom(baseB, "OP2");

	vector<PolarAtom*> paListA;
	vector<PolarAtom*> paListB;

	vector<Atom*>* scListA = baseA->getSidechainAtoms();
	vector<Atom*>* scListB = baseB->getSidechainAtoms();

	for(int i=0;i<scListA->size();i++){
		Atom* a = scListA->at(i);
		if(a->getType() == "O" || a->getType() == "N")
			paListA.push_back(new PolarAtom(baseA, a->getName()));
	}

	for(int i=0;i<scListB->size();i++){
		Atom* a = scListB->at(i);
		if(a->getType() == "O" || a->getType() == "N")
			paListB.push_back(new PolarAtom(baseB, a->getName()));
	}

	double hbondEne = -7.0*0.25;

	if(o2A->hbondedTo(o2B)) ene += hbondEne;

	for(int i=0;i<paListA.size();i++){
		if(o2B->hbondedTo(paListA[i])) ene += hbondEne;
		if(op1B->hbondedTo(paListA[i])) ene += hbondEne;
		if(op2B->hbondedTo(paListA[i])) ene += hbondEne;
	}

	for(int i=0;i<paListB.size();i++){
		if(o2A->hbondedTo(paListB[i])) ene += hbondEne;
		if(op1A->hbondedTo(paListB[i])) ene += hbondEne;
		if(op2A->hbondedTo(paListB[i])) ene += hbondEne;
	}

	/*
	Atom* a1 = baseA->getAtom("O2'");
	Atom* a2 = baseA->getAtom("O4'");
	Atom* a3 = baseB->getAtom("O2'");
	Atom* a4 = baseB->getAtom("O4'");

	if(a1 != NULL && a2 != NULL && a3 != NULL && a4 != NULL){
		double d1 = a1->distance(*a4);
		double d2 = a2->distance(*a3);
		if(d1 < 4.0 && d1 > 3.0) ene += -4.5;
		if(d2 < 4.0 && d2 > 3.0) ene += -4.5;
	}
	*/

	delete o2A;
	delete op1A;
	delete op2A;
	delete o2B;
	delete op1B;
	delete op2B;

	for(int i=0;i<paListA.size();i++){
		delete paListA[i];
	}
	for(int i=0;i<paListB.size();i++){
		delete paListB[i];
	}

	return ene;
}

BasePairLib::~BasePairLib() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
