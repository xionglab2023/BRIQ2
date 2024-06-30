
#ifndef MODEL_ASSIGNPROSSANDSASA_H_
#define MODEL_ASSIGNPROSSANDSASA_H_

#include "model/StructureModel.h"
#include "dataio/datapaths.h"
#include <fstream>

namespace NSPmodel {

using namespace NSPgeometry;
using namespace std;

class ResSasaPoints{
public:

	XYZ points[120];
	float sasaIndex[100];
	float wts[120];
	float radii;
	XYZ psdList[3];

	ResSasaPoints();
	float getSAI(vector<XYZ>& localTerns);
	float getSAI(Residue* resA, vector<Residue*>& resList);
	virtual ~ResSasaPoints();
};

class BackboneHBond{
public:
	string donerChainID;
	string acceptorChainID;
	string donerResID;
	string acceptorResID;
	int seqSeparation;
	XYZ hDoner;
	XYZ hydrogen;
	XYZ hAcceptor;
	float distance;
	float angle;

	BackboneHBond(XYZ& N, XYZ& H, XYZ& O, Residue* resA, Residue* resB);
	string toString();
	virtual ~BackboneHBond();
};


class AssignProSSAndSasa {
public:
	int resNum;
	vector<Residue*> resList;
	vector<XYZ*> NList;
	vector<XYZ*> CAList;
	vector<XYZ*> CList;
	char* ssSeq;
	map<string,int> chainIDResIDToSeqID;
	vector<BackboneHBond*> hbList;
	vector<float> saiList;
	vector<bool> backboneHBonded;

	AssignProSSAndSasa(PDB* protein);
	AssignProSSAndSasa(ProteinChain* pc);
	AssignProSSAndSasa(vector<Residue*>& resList);

private:
	void updateBBHbonds();

public:
	void updateSS();
	void updateSasa(ResSasaPoints* rsp);

	string getSS();
	float getSASA(Residue* res);
	char getSS(int seqID);
	float getSASA(int seqID);

	virtual ~AssignProSSAndSasa();
};

} /* namespace NSPmath */

#endif /* MODEL_ASSIGNPROSSANDSASA_H_ */
