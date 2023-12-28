

#ifndef MODEL_POLARGROUPLIB_LIB_H_
#define MODEL_POLARGROUPLIB_LIB_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "geometry/xyz.h"
#include "model/StructureModel.h"
#include "tools/StringTool.h"
#include <vector>
#include <fstream>
#include "dataio/datapaths.h"

namespace NSPmodel {

using namespace NSPgeometry;
using namespace NSPtools;
using namespace NSPdataio;

class PolarGroupRotamer{
public:
    string groupName;
    int atomNum;
    vector<string> atomTypes;
    vector<XYZ> localCoords;
    
    PolarGroupRotamer(const string& groupName);
    virtual ~PolarGroupRotamer();

};

class PolarGroupConformer{
public:
    string groupName;
    PolarGroupRotamer* rot;
    vector<XYZ> coords;
    LocalFrame cs;

    PolarGroupConformer(PolarGroupRotamer* rot, LocalFrame& cs);
    virtual ~PolarGroupConformer();
};

class PolarGroupRotamerLib{
public:
    vector<string> names;
    map<string, PolarGroupRotamer*> rotMap;

    PolarGroupRotamerLib();
    virtual ~PolarGroupRotamerLib();
};

}


#endif /*MODEL_POLARGROUPLIB_LIB_H_*/