/*
 * SeqDesignTemplateFast.h
 *
 *  Created on: 2022Äê7ÔÂ16ÈÕ
 *      Author: pengx
 */

#ifndef PREDPROT_SEQDESIGNTEMPLATEFAST_H_
#define PREDPROT_SEQDESIGNTEMPLATEFAST_H_

#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamerLibMini.h"
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
#include "protein/EnergyArray.h"
#include "protein/EnergyMatrix.h"
#include "protein/RotSequence.h"

namespace NSPprotein {

class SeqDesignTemplateFast {
private:
	PDB* pdb;
	vector<Residue*> resList;
	int resNum;
	map<string, int> chainIDResIDToIndex;

	vector<ResNode*> nodeList;
	vector<int> seqIDList;

	ResName* rn;
	ResBBRotamerLib* bbLib;
	ResScRotamerLibMini* scLib;
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
	vector<vector<ResConformer*>> resConfList;

	vector<EnergyArray*> eaList;
	vector<EnergyMatrix*> emList;
	RotSequence* rotSeq;

	double T0;
	double T1;
	int stepNumFactor;
	ProS1S2Energy* eS1S2;
	EnergyCalculator* ec;
	string initSequence;

public:
	SeqDesignTemplateFast(const string& inputFile, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib, AtomLib* atLib, EnergyCalculator* ec, ProS1S2Energy* eS1S2);
	SeqDesignTemplateFast(const string& inputFile, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib, AtomLib* atLib, EnergyCalculator* ec, ProS1S2Energy* eS1S2, const string& s2File);

	void printInvolvedPairInfo();
	void loadS1S2();
	void loadS1S2(const string& outfile);
	void printS2(const string& outfile);

	void loadEAEM();

	void loadRnaAsLigand(RNAPDB* rpdb);

	void printPairInfo();
	void printBackboneInfo();
	void printResInfo();
	void loadMSAProfile(const string& file, double weight);
	void checkConf();
	void updateRotChoice(const string& resFile);
	void printDetailEnergy();
	void getPositionRanks(ofstream& out);
	bool rotamerValid(int pos, ResScRotamer* rot);
	double totalEnergy();
	void printNatEnergy();


	void printDetailMutScSc(int pos, ResScRotamer* rot);
	double resInvolvedEnergy(int pos, int choice);
	double mutEnergy(int pos, int choice);
	double mutEnergy(int pos, ResScRotamer* rot);
	void mutEnergyDetail(int pos, ResScRotamer* rot, double* eList);
	void acceptMutation(int pos, ResScRotamer* rot);
	ResScRotamer* findBestRotamer(int pos);
	void designMC();
	string getDesignSequence();
	void printDesignPDB(const string& outfile);

	virtual ~SeqDesignTemplateFast();
};

} /* namespace NSPmodel */

#endif /* PREDPROT_SEQDESIGNTEMPLATEFAST_H_ */
