/*
 * LigandRotamer.h
 *
 */

#ifndef MODEL_LIGANDROTAMER_H_
#define MODEL_LIGANDROTAMER_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/LigandLib.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "model/AtomLib.h"
#include "model/StructureModel.h"
#include "tools/StringTool.h"

namespace NSPmodel {

using namespace NSPgeometry;
using namespace std;
using namespace NSPtools;

class LigandRotamer {

public:

	LigandInfo* ligInfo;
    int atomNum;
    XYZ* localTList;
 
    int polarAtomNum;
	CsMove* polarCmList;

	int rotIndex;
	double ene;

	LigandRotamer(LigandInfo* ligInfo, const string& rotLine);
	virtual ~LigandRotamer();
};

class LigandRotamerLib{

public:

	string ligandName;
	vector<LigandRotamer*> rots;

	LigandRotamerLib(LigandInfo* ligInfo);

	virtual ~LigandRotamerLib();
};


class LigandConformer{
public:

	LigandRotamer* rot;

	XYZ* coords;
	LocalFrame cs1;

	LocalFrame* csPolarList;

	LigandConformer(LigandRotamer* rot, LocalFrame& cs);
	void copyValueFrom(LigandConformer* other);
	void updateCoords(LocalFrame& cs);
	void updateRotamer(LigandRotamer* rot);
	void updateRotamerAndLocalFrame(LigandRotamer* rot, LocalFrame& cs);
	virtual ~LigandConformer();

};


}
#endif /* MODEL_LIGANDROTAMER_H_ */
