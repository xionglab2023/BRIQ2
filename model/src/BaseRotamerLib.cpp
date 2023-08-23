/*
 * BaseRotamerLib.cpp
 *
 */

#include "model/BaseRotamerLib.h"

namespace NSPmodel {

BaseRotamerLib::BaseRotamerLib(AtomLib* atLib){
	for(int i=0;i<8;i++){
		baseLib[i] = new BaseRotamer(i, atLib);
	}
}

void BaseRotamerLib::testRotamerLib(AtomLib* atLib){
	vector<string> names;
	names.push_back("A");
	names.push_back("U");
	names.push_back("G");
	names.push_back("C");
	names.push_back("DA");
	names.push_back("DT");
	names.push_back("DG");
	names.push_back("DC");

	LocalFrame cs0;
	for(int i=0;i<8;i++){
		BaseRotamer* rot = baseLib[i];
		int polarNum = rot->polarAtomNum;
		for(int j=0;j<polarNum;j++){
			int index = rot->polarAtomIndex[j];
			int uniqueID = rot->uniqueIDs[index];
			XYZ t = rot->coordsLocal[index];
			string name = atLib->uniqueIDToName(uniqueID);
			LocalFrame cs1 = cs0 + rot->polarCmList[j];
			double d = t.distance(cs1.origin_);
			cout << name << " " << d << endl;
		}
	}
}

BaseRotamerLib::~BaseRotamerLib(){
	for(int i=0;i<8;i++){
		delete baseLib[i];
	}
}

} /* namespace NSPforcefield */
