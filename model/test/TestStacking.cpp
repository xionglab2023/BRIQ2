

#include "model/RNABaseLib.h"
#include "geometry/RMSD.h"
#include "stdlib.h"
#include <iostream>
#include "model/StructureModel.h"


using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){
    cout << "init pdb" << endl;
	RNAPDB* pdb = new RNAPDB(string(argv[1]), "xxxx");

	vector<RNABase*> baseList = pdb->getBaseList();

    cout << "init atLib" << endl;
    AtomLib* atLib = new AtomLib();

    for(int i=0;i<baseList.size()-1;i++){
        RNABase* baseA = baseList[i];
        RNABase* baseB = baseList[i+1];


        if(baseA->isStackingTo(baseB, atLib)){
            cout << i << " " << (i+1) << " stacking" << endl; 
        }
        else {
            cout << i << " " << (i+1) << " not stacking" << endl; 
        }

    }

    delete pdb;

}
