
#ifndef MODEL_LIGANDLIB_H_
#define MODEL_LIGANDLIB_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "geometry/xyz.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include "model/AtomLib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPtools;

namespace NSPmodel {

class LigandInfo{

private:
    int curPolarAtomNum;

public:

    string ligandName;
    int atomNum;
    vector<string> atomNames;
    int* uniqueIDs;
    int polarAtomNum;

	int* polarAtomIndex;
	int* supportAtomIndex1;
	int* supportAtomIndex2;

	//support type: 0,1
	//0: sup2-sup1-polar, polar atom on the terminal
	//1: sup1-polar-sup2, polar atom in the middle
    int* supportType; 

    LigandInfo(const string& ligName, int atomNum, int polarAtomNum);
    void addAtom(const string& line, AtomLib* atLib);
    void printInfo();
    virtual ~LigandInfo();
};


class LigandLib{

public:

    map<string, LigandInfo*> ligMap;
    vector<LigandInfo*> ligandList;

    LigandLib();
    void printInfo();
    virtual ~LigandLib();

};




} /* namespace NSPmodel */

#endif /* MODEL_ATOMLIB_H_ */