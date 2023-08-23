/*
 * ResPeptideRotamer.cpp
 *
 */

#include "model/ResPeptideRotamer.h"

namespace NSPmodel {
ResPeptideRotamer::ResPeptideRotamer(const string& line) {
	// TODO Auto-generated constructor stub
	ResName rn;
	vector<string> spt;
	splitString(line," ",&spt);

	this->typeA = rn.triToInt(spt[0]);
	this->typeB = rn.triToInt(spt[2]);
	this->psi = atof(spt[1].c_str());
	this->phi = atof(spt[3].c_str());
	this->tList[0] = XYZ(atof(spt[4].c_str()), atof(spt[5].c_str()), atof(spt[6].c_str()));
	this->tList[1] = XYZ(atof(spt[7].c_str()), atof(spt[8].c_str()), atof(spt[9].c_str()));
	this->tList[2] = XYZ(atof(spt[10].c_str()), atof(spt[11].c_str()), atof(spt[12].c_str()));
	this->tList[3] = XYZ(atof(spt[13].c_str()), atof(spt[14].c_str()), atof(spt[15].c_str()));
	this->tList[4] = XYZ(atof(spt[16].c_str()), atof(spt[17].c_str()), atof(spt[18].c_str()));

	LocalFrame cs3 = generateLocalFrameResidueStyle(tList[0], tList[1], tList[2]);
	LocalFrame cs1 = generateLocalFrameResidueStyle(tList[1], tList[3], tList[4]);
	this->cm31 = cs1 - cs3;
	this->cm13 = cs3 - cs1;

	this->dm = DistanceMatrixHbond(cs3, cs1);
	this->rotID = 0;
}

ResPeptideRotamer::~ResPeptideRotamer() {
	// TODO Auto-generated destructor stub
}

}
