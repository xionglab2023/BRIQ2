/*
 * ResFoldingTree.h
 *
 */

#ifndef PREDPROT_RESFOLDINGTREE_H_
#define PREDPROT_RESFOLDINGTREE_H_

#include <iostream>
#include <stdlib.h>

#include "forcefield/RnaEnergyTable.h"
#include "model/StructureModel.h"
#include "model/ResName.h"
#include "geometry/localframe.h"
#include "protein/ResNode.h"
#include "protein/ResConnection.h"

namespace NSPprotein {

using namespace std;
using namespace NSPmodel;
using namespace NSPgeometry;
using namespace NSPforcefield;

class ResFoldingTree {
public:
	int seqLen;
	int* seq;

	bool* connectToDownstream;
	bool* chainBreakPoints;

	ResNode** nodes;

	vector<int> bbFixedNodes;
	vector<int> scFixedNodes;
	vector<int> bbFlexibleNodes;
	vector<int> scFlexibleNodes;

	ResNode* pseudoNode;
	ResNode* rootNode;

	vector<ResConnection*> fixedConnectionList;
	vector<ResConnection*> flexibleConnectionList;

	int* seqSepTable;

	double* allBbScE;
	double* allBbBbE;
	double* allScScE;
	double* allRotE;
	double* allPepE;

	double* tmpBbScE;
	double* tmpBbBbE;
	double* tmpScScE;
	double* tmpRotE;
	double* tmpPepE;


	ResFoldingTree(const string& inputFile);
	virtual ~ResFoldingTree();
};

} /* namespace NSPpred */

#endif /* PREDPROT_RESFOLDINGTREE_H_ */
