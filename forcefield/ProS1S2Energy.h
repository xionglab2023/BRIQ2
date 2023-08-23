/*
 * ProS1S2Energy.h
 *
 */

#ifndef FORCEFIELD_PROS1S2ENERGY_H_
#define FORCEFIELD_PROS1S2ENERGY_H_

#include "math/AAScoreArray.h"
#include "math/AAScoreMatrix.h"
#include "model/StructureModel.h"
#include "dataio/datapaths.h"
#include "para/DesignPara.h"

namespace NSPforcefield {

using namespace NSPmath;
using namespace NSPmodel;
using namespace std;
using namespace NSPpara;

class ProS1S2Energy {
public:
	vector<AAScoreArray> saList;
	map<string, vector<ResPairInfo*>> pairMap;
	map<string, vector<SaiPair*>> saiMap;
	map<string, int> repNums;
	vector<string> keys;
	string pathS2;
	//map<string, float> wtMap;
	float wtS1;
	float wtS2[45];
	int s2SubID;

	ProS1S2Energy(DesignPara* para, int s2ID);

private:
	void loadRepPoints(const string& key, const string& fileName);
	void loadAllRepPoints();
	void loadSaiPoints(const string& key, const string& fileName);
	void loadAllSaiPoints();

public:
	map<int, double> findSaiIndexList(ResPairInfo* rp);
	map<int, double> findRPIndexList(ResPairInfo* rp);
	AAScoreArray getS1(ResInfo* ri);
	AAScoreMatrix getS2(ResPairInfo* rp);
	AAScoreMatrix getS2NearestNb(ResPairInfo* rp);

	virtual ~ProS1S2Energy();
};

} /* namespace NSPmath */

#endif /* FORCEFIELD_PROS1S2ENERGY_H_ */
