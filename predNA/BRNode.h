/*
 * BRNode.h
 *
 */

#ifndef predNA_BRNODE_H_
#define predNA_BRNODE_H_
#include "model/StructureModel.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "geometry/Angles.h"
#include "model/AtomLib.h"
#include "model/RotamerLib.h"
#include "model/BaseRotamer.h"
#include "model/RiboseRotamer.h"
#include "model/PhosphateRotamer.h"

namespace NSPpredNA {

class BRConnection;

using namespace NSPgeometry;
using namespace NSPmodel;

class BRNode {
public:

	int baseType;
	int seqID;

	BaseConformer* baseConf;
	BaseConformer* baseConfTmp;

	RiboseConformer* riboseConf;
	RiboseConformer* riboseConfTmp;

	PhosphateConformer* phoConf;
	PhosphateConformer* phoConfTmp;

	bool fixed;
	bool connectToNeighbor;
	int groupID;

	BRNode* father;

	BRNode* leftChild; //pairing base
	BRNode* midChild; //sequential base, connect by base move
	BRNode* rightChild; //jump

	BRNode* bulge13Child; //bulge13
	BRNode* bulge14Child; //bulge14
	BRNode* reverseChild; //previous base
	BRNode* revBulge13Child; //reverse bulge13
	BRNode* revBulge14Child; //reverse bulge14

	BRConnection* upConnection; //pairing move, base move, ribo move

	vector<int> chainBreaks;

	vector<int> baseGroupA;
	vector<int> baseGroupC;

	vector<int> riboseGroupA;
	vector<int> riboseGroupC;

	vector<int> phoGroupA;
	vector<int> phoGroupC;


	BRNode(int baseType, int seqID, RotamerLib* rotLib) {
		this->baseType = baseType;
		this->seqID = seqID;
		this->fixed = false;
		this->connectToNeighbor = false;
		this->groupID = -1;

		LocalFrame cs1;

		this->baseConf = new BaseConformer(rotLib->baseRotLib->baseLib[baseType], cs1);
		this->baseConfTmp = new BaseConformer(rotLib->baseRotLib->baseLib[baseType], cs1);
		this->riboseConf = new RiboseConformer();
		this->riboseConfTmp = new RiboseConformer();
		this->phoConf = new PhosphateConformer();
		this->phoConfTmp = new PhosphateConformer();

		this->father = NULL;
		this->leftChild = NULL;
		this->rightChild = NULL;
		this->midChild = NULL;
		this->reverseChild = NULL;
		this->revBulge13Child = NULL;
		this->revBulge14Child = NULL;
		this->bulge13Child = NULL;
		this->bulge14Child = NULL;
		this->upConnection = NULL;
	}

	BRNode(RNABase* base, RotamerLib* rotLib){
		this->baseType = base->baseTypeInt;
		this->seqID = base->baseSeqID;
		this->fixed = false;
		this->connectToNeighbor = false;

		this->groupID = -1;
		LocalFrame cs1 = base->getCoordSystem();

		this->baseConf = new BaseConformer(rotLib->baseRotLib->baseLib[baseType], cs1);
		this->baseConfTmp = new BaseConformer(rotLib->baseRotLib->baseLib[baseType], cs1);

		this->riboseConf = new RiboseConformer(rotLib->riboseRotLib->rotLib[baseType][0], cs1);
		this->riboseConfTmp = new RiboseConformer(rotLib->riboseRotLib->rotLib[baseType][0], cs1);

		LocalFrame cs2 = riboseConf->cs2;

		this->phoConf = new PhosphateConformer(rotLib->phoRotLib->prLib[0][0], cs2);
		this->phoConfTmp = new PhosphateConformer(rotLib->phoRotLib->prLib[0][0], cs2);


		this->father = NULL;
		this->leftChild = NULL;
		this->rightChild = NULL;
		this->midChild = NULL;
		this->reverseChild = NULL;
		this->upConnection = NULL;
		this->revBulge13Child = NULL;
		this->revBulge14Child = NULL;
		this->bulge13Child = NULL;
		this->bulge14Child = NULL;
	}



	BRNode(RNABase* base, RiboseRotamer* riboRot, PhosphateRotamer* phoRot, RotamerLib* rotLib){
		this->baseType = base->baseTypeInt;
		this->seqID = base->baseSeqID;
		this->fixed = false;
		this->connectToNeighbor = false;

		this->groupID = -1;
		LocalFrame cs1 = base->getCoordSystem();

		this->baseConf = new BaseConformer(rotLib->baseRotLib->baseLib[baseType], cs1);
		this->baseConfTmp = new BaseConformer(rotLib->baseRotLib->baseLib[baseType], cs1);

		this->riboseConf = new RiboseConformer(riboRot, cs1);
		this->riboseConfTmp = new RiboseConformer(riboRot, cs1);

		LocalFrame cs2 = riboseConf->cs2;

		this->phoConf = new PhosphateConformer(phoRot, cs2);
		this->phoConfTmp = new PhosphateConformer(phoRot, cs2);

		this->father = NULL;
		this->leftChild = NULL;
		this->rightChild = NULL;
		this->midChild = NULL;
		this->reverseChild = NULL;
		this->upConnection = NULL;
		this->revBulge13Child = NULL;
		this->revBulge14Child = NULL;
		this->bulge13Child = NULL;
		this->bulge14Child = NULL;
	}


	bool baseConsistent();
	bool riboConsistent();
	bool rotamerConsistent();

	bool consistent();
	bool phoConsistent();
	bool phoLocalConsistent();

	void updateChildInfo(BRNode** allNodes, int treeSize){

		if(this->seqID > 0 && allNodes[seqID-1]->connectToNeighbor){
			this->chainBreaks.push_back(seqID-1);
		}


		if(connectToNeighbor)
			this->chainBreaks.push_back(seqID);

		for(int i=0;i<treeSize;i++){
			if(i != this->seqID)
				this->baseGroupA.push_back(i);
			else
				this->baseGroupC.push_back(i);
		}

		int* riboseInfo = new int[treeSize];
		int* phoInfo = new int[treeSize];
		for(int i=0;i<treeSize;i++){
			riboseInfo[i] = 0;
			phoInfo[i] = 0;
		}
		riboseInfo[seqID] = 2;

		for(int i=0;i<chainBreaks.size();i++){
			phoInfo[chainBreaks[i]] = 2;
		}
		for(int i=0;i<treeSize;i++){
			if(riboseInfo[i]==0)
				this->riboseGroupA.push_back(i);
			else
				this->riboseGroupC.push_back(i);
		}
		for(int i=0;i<treeSize;i++){
			if(phoInfo[i] == 0 && allNodes[i]->connectToNeighbor)
				this->phoGroupA.push_back(i);
			if(phoInfo[i] == 2 && allNodes[i]->connectToNeighbor)
				this->phoGroupC.push_back(i);
		}
		delete riboseInfo;
		delete phoInfo;
	}

	void printPartition(){
		cout << "single node: " << this->seqID << endl;
		cout << "base groupA:";
		for(int i=0;i<baseGroupA.size();i++){
			cout << " " << baseGroupA[i];
		}
		cout << " groupC:";
		for(int i=0;i<baseGroupC.size();i++){
			cout << " " << baseGroupC[i];
		}
		cout << endl;

		cout << "ribose groupA:";
		for(int i=0;i<riboseGroupA.size();i++){
			cout << " " << riboseGroupA[i];
		}
		cout << " groupC:";
		for(int i=0;i<riboseGroupC.size();i++){
			cout << " " << riboseGroupC[i];
		}
		cout << endl;

		cout << "pho groupA:";
		for(int i=0;i<phoGroupA.size();i++){
			cout << " " << phoGroupA[i];
		}
		cout << " groupC:";
		for(int i=0;i<phoGroupC.size();i++){
			cout << " " << phoGroupC[i];
		}
		cout << endl;
	}

	BRNode& operator=(const BRNode& other);
	void copyValueFrom(const BRNode& other);

	vector<Atom*> toAtomList(AtomLib& atLib);
	vector<Atom*> toBaseAtomList(AtomLib& atLib);
	vector<Atom*> phoAtoms();
	vector<Atom*> toTmpAtomList(AtomLib& atLib);

	void checkRotamer();
	void checkTmpRotamer();

	virtual ~BRNode();
};


} /* namespace NSPpred */

#endif /* predNA_BRNODE_H_ */
