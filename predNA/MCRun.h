/*
 * MCRun.h
 *
 */

#ifndef predNA_MCRUN_H_
#define predNA_MCRUN_H_

#include "predNA/BRFoldingTree.h"

namespace NSPpredna {

using namespace NSPmodel;
using namespace std;

class MCRun {

	//BRFoldingTreeBasic* simpFt;

	//BRTreeInfoBasic* initBasic;

public:
	BRFoldingTree* ft;
	BRTreeInfo* init;
	MCRun(BRFoldingTree* ft){
		this->ft = ft;
		init = ft->getTreeInfo();
	}

	void generateDecoysRandInit(const string& output);

	void optimizeFromInit(const string& keyFile, const string& outFilePrefix, const int startID);

	void simpleMC(const string& outpdb, bool traj);

	void optimize(double t0, double kStep);

	void optimizeBackbone(const string& output);

	void debug();

	void test();

	virtual ~MCRun();
};

} /* namespace NSPpred */

#endif /* predNA_MCRUN_H_ */
