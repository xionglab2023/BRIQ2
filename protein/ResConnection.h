/*
 * ResConnection.h
 *
 */

#ifndef PREDPROT_RESCONNECTION_H_
#define PREDPROT_RESCONNECTION_H_
#include "model/StructureModel.h"
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "protein/ResNode.h"
#include "model/ResBBRotamer.h"
#include "model/ResScRotamer.h"
#include "model/ResPeptideRotamer.h"

namespace NSPprotein {

using namespace NSPmodel;
using namespace NSPgeometry;

class ResConnection {

public:

	CsMove cm;
	CsMove tmp;
	ResNode* father;
	ResNode* child;
	string ctType;

	int treeSize;
	vector<int> resGroupA;
	vector<int> resGroupB;
	vector<int> chainBreaks;
	vector<ResConnection*> childConnections;
	bool isFixed;


	ResConnection(const string& ctType, ResNode* father, ResNode* child);
	void updateChildInfo(ResNode** allNodes);


	virtual ~ResConnection();
};

} /* namespace NSPpred */

#endif /* PREDPROT_RESCONNECTION_H_ */
