#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/BaseRotamer.h"
#include "model/BasePairLib.h"
#include "predNA/NuGraph.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;

int main(int argc, char** argv){

    if(argc < 5){
        cout << "indexToPDB sp1000 $TYPEA $TYPEB $INDEX $OUTPDB" << endl;
        cout << "indexToPDB sp2000 $TYPEA $TYPEB $INDEXA $INDEXB $OUTPDB" << endl;
        exit(0);
    }

    string spType = string(argv[1]);

    OrientationIndex* oi = new OrientationIndex();
    AtomLib* atLib = new AtomLib();

    string augc = "AUGC";
    if(spType == "sp1000") {
        int typeA = atoi(argv[2]);
        int typeB = atoi(argv[3]);
        int index = atoi(argv[4]);
        string outfile = string(argv[5]);
        cout << "outFile: " << outfile << endl;
        ofstream out = ofstream(outfile.c_str(), ios::out);

        BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
        BaseRotamer* rotB = new BaseRotamer(typeB, atLib);
        
        LocalFrame csA;
        LocalFrame csB;
        CsMove cm = oi->index1000ToCsMove(index);
        csB = csA + cm;

        BaseConformer* confA = new BaseConformer(rotA, csA);
        BaseConformer* confB = new BaseConformer(rotB, csB);

        RNABase* baseA = new RNABase("1", "A", augc[typeA]);
        RNABase* baseB = new RNABase("2", "A", augc[typeB]);

    	vector<Atom*> atomListA;
        vector<string> namesA;
	    atLib->getRnaSidechainAtoms(typeA, namesA);
	    for(int i=0;i<rotA->atomNum;i++){
            baseA->addAtom(new Atom(namesA.at(i), confA->coords[i]));
	    }

    	vector<Atom*> atomListB;
        vector<string> namesB;
	    atLib->getRnaSidechainAtoms(typeB, namesB);
	    for(int i=0;i<rotB->atomNum;i++){
            baseB->addAtom(new Atom(namesB.at(i), confB->coords[i]));
	    }

        RNAChain* rc = new RNAChain("A");
        rc->addBase(baseA);
        rc->addBase(baseB);
        rc->printPDBFormat(out, 1);
        out.close();

        delete rc;
        for(int i=0;i<rotA->atomNum;i++){
            delete baseA->getAtomList()->at(i);
	    }
	    for(int i=0;i<rotB->atomNum;i++){
            delete baseB->getAtomList()->at(i);
	    }

        delete baseA;
        delete baseB;
        delete rotA;
        delete rotB;
    }
    else if(spType == "sp2000") {
        int typeA = atoi(argv[2]);
        int typeB = atoi(argv[3]);
        int indexA= atoi(argv[4]);
        int indexB = atoi(argv[5]);
        string outfile = string(argv[6]);
        ofstream out = ofstream(outfile.c_str(), ios::out);

        BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
        BaseRotamer* rotB = new BaseRotamer(typeB, atLib);
        
        LocalFrame csA;
        LocalFrame csB;
        CsMove cm = oi->index2000ToCsMove(indexA, indexB);
        csB = csA + cm;

        BaseConformer* confA = new BaseConformer(rotA, csA);
        BaseConformer* confB = new BaseConformer(rotB, csB);

        RNABase* baseA = new RNABase("1", "A", augc[typeA]);
        RNABase* baseB = new RNABase("2", "A", augc[typeB]);

    	vector<Atom*> atomListA;
        vector<string> namesA;
	    atLib->getRnaSidechainAtoms(typeA, namesA);
	    for(int i=0;i<rotA->atomNum;i++){
            baseA->addAtom(new Atom(namesA.at(i), confA->coords[i]));
	    }

    	vector<Atom*> atomListB;
        vector<string> namesB;
	    atLib->getRnaSidechainAtoms(typeB, namesB);
	    for(int i=0;i<rotB->atomNum;i++){
            baseB->addAtom(new Atom(namesB.at(i), confB->coords[i]));
	    }

        RNAChain* rc = new RNAChain("A");
        rc->addBase(baseA);
        rc->addBase(baseB);
        rc->printPDBFormat(out, 1);
        out.close();

        delete rc;
        for(int i=0;i<rotA->atomNum;i++){
            delete baseA->getAtomList()->at(i);
	    }
	    for(int i=0;i<rotB->atomNum;i++){
            delete baseB->getAtomList()->at(i);
	    }

        delete baseA;
        delete baseB;
        delete rotA;
        delete rotB;        
    }
    else {
        cout << "invalid sp type: " << spType << endl;
    }
    

}