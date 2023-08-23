/*
 * AssignRNASasa.h
 *
 */

#ifndef MODEL_ASSIGNRNASASA_H_
#define MODEL_ASSIGNRNASASA_H_

#include "model/RNABaseLib.h"
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "predNA/BRNode.h"
#include "model/StructureModel.h"
namespace NSPmodel {

using namespace NSPpredna;
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
	AssignRNASasa(vector<BRNode*>& baseList, BaseSasaPoints* bs);
	int exposeNum(vector<XYZ>& localTerns, int type);
	int exposeNum(RNABase* base, vector<RNABase*>& baseList);
	int exposeNum(BRNode* node, vector<BRNode*>& nodeList);
	void printExposeNum();
	virtual ~AssignRNASasa();
};

} /* namespace NSPalign */

#endif /* MODEL_ASSIGNRNASASA_H_ */
