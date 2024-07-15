/*
 * AssignRNASasa.cpp
 *
 */

#include "model/AssignRNASasa.h"

namespace NSPmodel {

BaseSasaPoints::BaseSasaPoints(){
	string path = NSPdataio::datapath();
	string s;
	double x,y,z;
	this->pointNum = 200;
	vector<XYZ> aList;
	vector<XYZ> uList;
	vector<XYZ> gList;
	vector<XYZ> cList;

	string libFile = path+"sasa/A.pts";
	ifstream file;

	file.open(libFile.c_str(), ios::in);
	if(! file.is_open())
	{
		 cout << "fail to open file " << libFile << endl;
		 exit (1);
	}

	while(file >> x >> y >> z)
	{
		XYZ t(x,y,z);
		aList.push_back(t);
	}
	file.close();

	libFile = path+"sasa/U.pts";
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open())
	{
		 cout << "fail to open file " << libFile << endl;
		 exit (1);
	}

	while(file >> x >> y >> z)
	{
		XYZ t(x,y,z);
		uList.push_back(t);
	}
	file.close();

	libFile = path+"sasa/G.pts";
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open())
	{
		 cout << "fail to open file " << libFile << endl;
		 exit (1);
	}

	while(file >> x >> y >> z)
	{
		XYZ t(x,y,z);
		gList.push_back(t);
	}
	file.close();

	libFile = path+"sasa/C.pts";
	file.open(libFile.c_str(), ios::in);
	if(! file.is_open())
	{
		 cout << "fail to open file " << libFile << endl;
		 exit (1);
	}
	while(file >> x >> y >> z)
	{
		XYZ t(x,y,z);
		cList.push_back(t);
	}
	file.close();
	this->sasaPoints.push_back(aList);
	this->sasaPoints.push_back(uList);
	this->sasaPoints.push_back(gList);
	this->sasaPoints.push_back(cList);
}

BaseSasaPoints::~BaseSasaPoints(){

}

AssignRNASasa::AssignRNASasa(vector<RNABase*>& baseList, BaseSasaPoints* bs) {
	// TODO Auto-generated constructor stub
	this->bs = bs;
	this->len = baseList.size();
	this->expNum = new int[len];
	this->types = new int[len];
	this->radii = 3.2;

	RNABaseLib baseLib;
	this->localCoords.push_back(baseLib.getPolarAtomCoords(0));
	this->localCoords.push_back(baseLib.getPolarAtomCoords(1));
	this->localCoords.push_back(baseLib.getPolarAtomCoords(2));
	this->localCoords.push_back(baseLib.getPolarAtomCoords(3));

	for(int i=0;i<len;i++){
		this->expNum[i] = exposeNum(baseList[i], baseList);
		this->types[i] = baseList[i]->baseTypeInt;
	}
}

AssignRNASasa::AssignRNASasa(BaseSasaPoints* bs, vector<NuNode*>& nodes){
	this->bs = bs;
	this->len = nodes.size();
	this->expNum = new int[len];
	this->types = new int[len];
	this->radii = 3.2;

	RNABaseLib baseLib;
	this->localCoords.push_back(baseLib.getPolarAtomCoords(0));
	this->localCoords.push_back(baseLib.getPolarAtomCoords(1));
	this->localCoords.push_back(baseLib.getPolarAtomCoords(2));
	this->localCoords.push_back(baseLib.getPolarAtomCoords(3));

	for(int i=0;i<len;i++){
		this->expNum[i] = exposeNum(nodes[i], nodes);
		this->types[i] = nodes[i]->baseConf->rot->baseType;
	}
}

int AssignRNASasa::exposeNum(vector<XYZ>& terns, int type){
	bool buried[this->bs->pointNum];
	for(int i=0;i<bs->pointNum;i++){
		buried[i] = false;
	}
	int n=0;
	int pn = bs->pointNum;
	int tn = terns.size();
	double rr = radii*radii;
	for(int i=0;i<pn;i++){
		XYZ ball = bs->sasaPoints[type][i];
		int j=0;
		while(!buried[i] && j<tn){
			if(squareDistance(ball, terns[j]) < rr)
				buried[i] = true;
			j++;
		}
		if(!buried[i]) n++;
	}
	return n;
}

int AssignRNASasa::exposeNum(RNABase* base, vector<RNABase*>& baseList){
	vector<XYZ> terns;
	int type = base->baseTypeInt;
	LocalFrame cs = base->getCoordSystem();
	vector<XYZ> baseAtoms;
	for(int i=0;i<this->localCoords[type].size();i++){
		baseAtoms.push_back(cs.local2globalcrd(this->localCoords[type][i]));
	}
	int na = baseAtoms.size();
	double RR = 4*radii*radii;
	int i,j,k;
	for(i=0;i<baseList.size();i++){
		RNABase* baseB = baseList[i];
		LocalFrame cs2 = baseB->getCoordSystem();
		double od = cs.origin_.distance(cs2.origin_);
		if(od < 0.1 || od > 17.0) continue;

		vector<Atom*>* atoms = baseB->getAtomList();
		for(j=0;j<atoms->size();j++){
			if(atoms->at(j)->type == "H") continue;
			bool select = false;
			XYZ t = atoms->at(j)->coord;
			for(k=0;k<na;k++){
				if(squareDistance(t, baseAtoms[k]) < RR){
					select = true;
					break;
				}
			}
			if(select){
				terns.push_back(cs.global2localcrd(t));
			}
		}
	}

	vector<Atom*>* atoms = base->getBackboneAtoms();
	for(j=0;j<atoms->size();j++){
		XYZ t = atoms->at(j)->coord;
		terns.push_back(cs.global2localcrd(t));
	}

	//cout << base->baseSeqID << " " << terns.size() << endl;
	return exposeNum(terns, type);
}

int AssignRNASasa::exposeNum(NuNode* node, vector<NuNode*>& nodeList){
	vector<XYZ> terns;
	int type = node->baseConf->rot->baseType;
	LocalFrame cs = node->baseConf->cs1;

	vector<XYZ> baseAtoms;
	for(int i=0;i<this->localCoords[type].size();i++){
		baseAtoms.push_back(cs.local2globalcrd(this->localCoords[type][i]));
	}
	int na = baseAtoms.size();
	double RR = 4*radii*radii;
	int i,j,k;

	for(i=0;i<nodeList.size();i++){
		NuNode* node = nodeList[i];
		LocalFrame cs2 = node->baseConf->cs1;
		double od = cs.origin_.distance(cs2.origin_);
		if(od < 0.1 || od > 17.0) continue;

		for(j=0;j<node->baseConf->rot->atomNum;j++){
			terns.push_back(cs.local2globalcrd(node->baseConf->coords[j]));
		}

		for(j=0;j<node->riboseConf->rot->atomNum;j++){
			terns.push_back(cs.local2globalcrd(node->riboseConf->coords[j]));
		}

		for(j=0;j<4;j++){
			terns.push_back(cs.local2globalcrd(node->phoConf->coords[j]));
		}
	}

	return exposeNum(terns, type);
}

void AssignRNASasa::printExposeNum(){
	for(int i=0;i<this->len;i++){
		printf("%-4d %d %3d\n", i, this->types[i], this->expNum[i]);
	}
}



AssignRNASasa::~AssignRNASasa() {
	// TODO Auto-generated destructor stub
	delete [] expNum;
}

} /* namespace NSPalign */
