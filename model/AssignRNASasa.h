/*
 * AssignRNASasa.h
 *
 */

#ifndef MODEL_ASSIGNRNASASA_H_
#define MODEL_ASSIGNRNASASA_H_

#include "model/RNABaseLib.h"
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "predNA/NuGraph.h"
#include "model/StructureModel.h"

namespace NSPmodel {

using namespace NSPpredNA;
using namespace NSPgeometry;
using namespace std;

class BaseSasaPoints{
public:
	int pointNum;
	vector<vector<XYZ>> sasaPoints;
	BaseSasaPoints();
	virtual ~BaseSasaPoints();
};

class AssignRNASasa {

private:
	BaseSasaPoints* bs;
	int len;
	int* expNum;
	int* types;
	double radii;
	vector<vector<XYZ>> localCoords;

public:
	AssignRNASasa(vector<RNABase*>& baseList, BaseSasaPoints* bs);
	AssignRNASasa(BaseSasaPoints* bs, vector<NuNode*>& nodeList);
	int exposeNum(vector<XYZ>& localTerns, int type);
	int exposeNum(RNABase* base, vector<RNABase*>& baseList);
	int exposeNum(NuNode* node, vector<NuNode*>& nodeList);
	int getExposeNum(int id){
		return expNum[id];
	}
	void printExposeNum();
	virtual ~AssignRNASasa();
};

} /* namespace NSPalign */

#endif /* MODEL_ASSIGNRNASASA_H_ */
