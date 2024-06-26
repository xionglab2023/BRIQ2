
#include "model/AssignRNASS.h"
#include "model/AssignRNASasa.h"
#include "model/StructureModel.h"
#include "model/AtomLib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){
	string pdbFile = string(argv[1]);
    AtomLib* atLib = new AtomLib();

    {
        RNAPDB* pdb = new RNAPDB(pdbFile);
        vector<RNABase*> baseList = pdb->getValidBaseList(atLib);
        string s = "";
        AssignRNASS* as = new AssignRNASS(pdb, atLib);
    
        string seq = as->seq;
        string sec = as->ssSeq;
        cout << seq << endl;
        cout << sec << endl;

        BaseSasaPoints* bsp = new BaseSasaPoints();
        AssignRNASasa* sasa = new AssignRNASasa(baseList, bsp);
        sasa->printExposeNum();

    }
 




}