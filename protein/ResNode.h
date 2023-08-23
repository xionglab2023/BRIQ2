/*
 * ResNode.h
 *
 */

#ifndef PREDPROT_RESNODE_H_
#define PREDPROT_RESNODE_H_

#include "model/StructureModel.h"
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/ResBBRotamer.h"
#include "model/ResScRotamer.h"
#include "model/ResConformer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamerLib.h"
#include "model/ResScRotamerLibMini.h"
#include "model/ResPeptideRotamer.h"

namespace NSPprotein {

using namespace NSPgeometry;
using namespace NSPmodel;
using namespace std;

class ResConnection;

class ResNode {
public:

	int aaType;
	int nodeID; //from 0 to N-1, N is the number of nodes
	int seqID; //indicate sequence seperation
	float sai;
	char ss;
	int ssInt;

	LocalFrame cs1; //pre_C-N-CA, origin N
	LocalFrame cs1Tmp;
	LocalFrame cs2; //N-CA-C, origin CA
	LocalFrame cs2Tmp;
	LocalFrame cs3; //CA-C-O, origin C
	LocalFrame cs3Tmp;

	ResConformer* conf;
	ResConformer* confTmp;

	bool scRotmerNeedDelete;

	ResNode* father; //father node
	ResNode* leftChild; //C-terminal sequential node
	ResNode* rightChild; //N-terminal sequential node
	ResNode* hbondAChild; //N to O hbond child
	ResNode* hbondBChild; //O to N hbond child
	ResNode* jumpChild; //jump
	ResConnection* upConnection;


	ResNode();
	ResNode(const string& line, int nodeID, ResName* rn, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib);
	ResNode(Residue* resP, Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib);
	ResNode(Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib);
	ResNode(Residue* resP, Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib);
	ResNode(Residue* res, int seqID, int nodeID, AtomLib* atLib, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib);

	XYZ getCbCoord(){
		XYZ localCb(-0.942, 0.009, 1.208);
		return local2global(cs2, localCb);
	}

	virtual ~ResNode();
};

}



#endif /* PREDPROT_RESNODE_H_ */
