/*
 * BuildMutationTemplate.h
 *
 */

#ifndef PREDPROT_BUILDMUTATIONTEMPLATE_H_
#define PREDPROT_BUILDMUTATIONTEMPLATE_H_

#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/AtomLib.h"
#include "model/ResScRotamerLibMini.h"
#include "forcefield/ProEnergyTable.h"
#include "forcefield/ProS1S2Energy.h"
#include "protein/ResNode.h"
#include "protein/EnergyCalculator.h"
#include <float.h>
#include "model/StructureModel.h"
#include "math/AAScoreArray.h"
#include "math/AAScoreMatrix.h"

namespace NSPprotein {

using namespace NSPmodel;
using namespace std;
using namespace NSPforcefield;

class BuildMutationTemplate {
public:

	ResNode* targetNode;
	vector<ResNode*> neighborNodes;
	int preAAType;
	vector<int> sepList;
	ResName* rn;
	AtomLib* atLib;
	ResBBRotamerLib* bbLib;
	ResScRotamerLib* scLib; //scRotLib tag: bbAll
	ResScRotamerLibMini* scLibSimp;
	ResInfo* riTarget;
	vector<int> nbAAList;
	vector<ResPairInfo*> rpList;
	vector<bool> targetFirst;
	AAScoreArray sa;
	vector<AAScoreMatrix> smList;


	BuildMutationTemplate(vector<string>& lines, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, ResScRotamerLibMini* scLibSimp, AtomLib* atLib, ResName* rn);

	void printBBInfo();
	void printNodeSep();

	void printNativeEnergy(EnergyCalculator* ec, ProS1S2Energy* eS1S2);

	void loadS1S2(ProS1S2Energy* eS1S2);
	void printS1S2(ofstream& out);
	void loadS1S2(vector<string>& lines);

	int rankS1(ProS1S2Energy* eS1S2);
	int rankS2(ProS1S2Energy* eS1S2);
	int rankS1S2(ProS1S2Energy* eS1S2);
	int rankAtomic(EnergyCalculator* ec);
	int rankAtomicABACUS(EnergyCalculator* ec);

	int nativeRankDetail(EnergyCalculator* ec, ProS1S2Energy* eS1S2);
	string nativeRank(EnergyCalculator* ec, ProS1S2Energy* eS1S2);
	int nativeRankABACUS(EnergyCalculator* ec, ProS1S2Energy* eS1S2);
	double nativeLogits(EnergyCalculator* ec, ProS1S2Energy* eS1S2);
	int singleResidueDesign(EnergyCalculator* ec, ProS1S2Energy* eS1S2);

	virtual ~BuildMutationTemplate();
};

} /* namespace NSPmodel */

#endif /* PREDPROT_BUILDMUTATIONTEMPLATE_H_ */
