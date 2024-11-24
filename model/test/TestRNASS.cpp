
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

    cout << pdbFile << endl;
    
    /*
    string outFile = string(argv[2]);
    

    RNAPDB* pdb = new RNAPDB(pdbFile);

    ofstream out;
    out.open(outFile.c_str(), ios::out);
    pdb->printCIFFormat(out);
    out.close();
    */


    
    {
        RNAPDB* pdb = new RNAPDB(pdbFile);
        vector<RNABase*> baseList = pdb->getValidBaseList(atLib);
        string s = "";
        AssignRNASS* as = new AssignRNASS(pdb, atLib);
    
        string seq = as->seq;
        string sec = as->ssSeq;
        cout << seq << endl;
        cout << sec << endl;

        int breakNum = as->breakList.size();
        for(int i=0;i<breakNum;i++){
            cout << as->breakList[i] << endl;
        } 

    //    BaseSasaPoints* bsp = new BaseSasaPoints();
    //    AssignRNASasa* sasa = new AssignRNASasa(baseList, bsp);
    //    sasa->printExposeNum();

    }
    




}