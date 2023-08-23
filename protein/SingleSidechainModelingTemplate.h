/*
 * SingleSidechainModelingTemplate.h
 *
 */

#ifndef PREDPROT_SINGLESIDECHAINMODELINGTEMPLATE_H_
#define PREDPROT_SINGLESIDECHAINMODELINGTEMPLATE_H_

#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamerLib.h"
#include "forcefield/ProEnergyTable.h"
#include "protein/ResNode.h"
#include "protein/EnergyCalculator.h"
#include <float.h>
#include "model/StructureModel.h"

namespace NSPprotein {

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace std;

class SingleSidechainModelingTemplate {
public:
	ResNode* targetNode;
	vector<ResNode*> neighborNodes;
	int preAAType;
	vector<int> sepList;
	ResName* rn;
	ResBBRotamerLib* bbLib;
	ResScRotamerLib* scLib;

	SingleSidechainModelingTemplate(vector<string>& lines, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, ResName* rn);
	double bestRotamerRMSD(EnergyCalculator* ec); //search all rotamers
	double bestRotamerRMSD2(EnergyCalculator* ec); //search level1 and level2 rotamers
	void printNativeEnergy(EnergyCalculator* ec);
	void printInit(const string& output, AtomLib* atLib);
	void printPDB(const string& rawState, const string& finalState, AtomLib* atLib);
	virtual ~SingleSidechainModelingTemplate();
};

} /* namespace NSPforcefield */

#endif /* PREDPROT_SINGLESIDECHAINMODELINGTEMPLATE_H_ */
