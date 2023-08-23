/*
 * LoopModelingTemplate.h
 *
 */

#ifndef PREDPROT_LOOPMODELINGTEMPLATE_H_
#define PREDPROT_LOOPMODELINGTEMPLATE_H_

#include "model/StructureModel.h"
#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamerLib.h"
#include "forcefield/ProEnergyTable.h"
#include "protein/ResNode.h"
#include "protein/EnergyCalculator.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace std;


namespace NSPprotein {

class LoopModelingTemplate {
public:
	int loopLength;
	vector<ResNode*> targetNodes;
	vector<ResNode*> neighborNodes;
	vector<int> seqIDList;
	ResName* rn;
	ResBBRotamerLib* bbLib;
	ResScRotamerLib* scLib;

	LoopModelingTemplate(int len, vector<string>& lines, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, ResName* rn);
	void printInit(const string& output, AtomLib* atLib);

	virtual ~LoopModelingTemplate();
};

} /* namespace NSPmodel */

#endif /* PREDPROT_LOOPMODELINGTEMPLATE_H_ */
