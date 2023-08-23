/*
 * SidechainModelingTemplate.h
 *
 */

#ifndef PREDPROT_SIDECHAINMODELINGTEMPLATE_H_
#define PREDPROT_SIDECHAINMODELINGTEMPLATE_H_

#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamerLib.h"
#include "forcefield/ProEnergyTable.h"
#include "protein/ResNode.h"
#include "protein/EnergyCalculator.h"
#include <float.h>
#include <stdlib.h>
#include "model/StructureModel.h"

namespace NSPprotein {

using namespace NSPmodel;
using namespace std;

class SidechainModelingTemplate {

public:
	PDB* pdb;
	vector<Residue*> resList;
	vector<ResNode*> nodeList;
	vector<int> seqIDList;
	ResName* rn;
	ResBBRotamerLib* bbLib;
	ResScRotamerLib* scLib;
	AtomLib* atLib;

	vector<bool> resValidList;
	vector<int> initScIDList;
	vector<ResScRotamer*> initRotList;
	vector<vector<int>> nbList;

	SidechainModelingTemplate(const string& pdbFile);
	double mutEnergy(int pos, ResScRotamer* rot, EnergyCalculator* ec);
	void acceptMutation(int pos, ResScRotamer* rot);
	ResScRotamer* findBestRotamer(int pos, EnergyCalculator* ec);
	void runMC(EnergyCalculator* ec);
	void printRMS(ofstream& out);
	void printPDB(const string& outfile);
	void printEnergy();
	double rms1();
	double rms2();

	virtual ~SidechainModelingTemplate();
};

} /* namespace NSPmodel */

#endif /* PREDPROT_SIDECHAINMODELINGTEMPLATE_H_ */
