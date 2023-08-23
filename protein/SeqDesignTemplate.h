/*
 * SeqDesignTemplate.h
 *
 */

#ifndef PREDPROT_SEQDESIGNTEMPLATE_H_
#define PREDPROT_SEQDESIGNTEMPLATE_H_

#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamerLib.h"
#include "model/AssignProSSAndSasa.h"
#include "forcefield/ProEnergyTable.h"
#include "forcefield/ProS1S2Energy.h"
#include "protein/ResNode.h"
#include "protein/EnergyCalculator.h"
#include <float.h>
#include <stdlib.h>
#include "model/StructureModel.h"
#include "math/AAScoreArray.h"
#include "math/AAScoreMatrix.h"
#include "tools/InputParser.h"
#include "protein/ResMutator.h"

namespace NSPprotein {

using namespace NSPmodel;
using namespace NSPmath;
using namespace NSPtools;

class SeqDesignTemplate {

private:
	PDB* pdb;
	vector<Residue*> resList;
	int resNum;
	map<string, int> chainIDResIDToIndex;

	vector<ResNode*> nodeList;
	vector<int> seqIDList;

	ResName* rn;
	ResBBRotamerLib* bbLib;
	ResScRotamerLib* scLib;

	AtomLib* atLib;

	vector<bool> resValidList;
	vector<int> initScIDList;
	vector<ResScRotamer*> initRotList;
	vector<vector<int>> nbList; //neighbors for atomic energy calculation

	vector<ResInfo*> riList;
	vector<ResPairInfo*> rpList;
	vector<vector<int>> involvedRpIndex; //neighbors for S2 calculation
	vector<AAScoreArray> saList;
	vector<AAScoreMatrix> smList;

	vector<ResMutator*> resMutatorList;


	double T0;
	double T1;
	int stepNumFactor;
	ProS1S2Energy* eS1S2;
	EnergyCalculator* ec;
	string initSequence;

public:
	SeqDesignTemplate(const string& inputFile, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, AtomLib* atLib, EnergyCalculator* ec, ProS1S2Energy* eS1S2);

	void printInvolvedPairInfo();
	void loadS1S2();
	void printS2(const string& outfile);
	void loadS1S2(const string& s2File);
	void printPairInfo();
	void printBackboneInfo();
	void printResInfo();
	void loadMSAProfile(const string& file, double weight);
	void checkConf();
	void updateRotChoice(const string& resFile, double eneCutoff);
	void getPositionRanks(ofstream& out);
	double rotEnergyWithBackbone(int pos, ResScRotamer* rot);
	double totalEnergy();
	void printPolarEnergy();
	void printEnergyDetail();
	void totalEnergyDetail(double* eList);
	void printDetailScSc();
	void checkScSc();
	void printDetailMutScSc(int pos, ResScRotamer* rot);
	double mutEnergy(int pos, ResScRotamer* rot);
	void mutEnergyDetail(int pos, ResScRotamer* rot, double* eList);
	void acceptMutation(int pos, ResScRotamer* rot);
	ResScRotamer* findBestRotamer(int pos);
	void runMC();
	string getSequence();
	double identity();
	void printPDB(const string& outFile);
	virtual ~SeqDesignTemplate();
};

} /* namespace NSPmodel */

#endif /* PREDPROT_SEQDESIGNTEMPLATE_H_ */
