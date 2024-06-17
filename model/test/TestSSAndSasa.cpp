
#include "model/AssignProSSAndSasa.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResBBRotamer.h"


using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){
	string pdbFile = string(argv[1]);
	string outFile = string(argv[2]);
	PDB* pdb = new PDB(pdbFile, "xxxx");
	ofstream out;
	out.open(outFile.c_str(), ios::out);
	AssignProSSAndSasa* as = new AssignProSSAndSasa(pdb);
	as->updateSS();
	ResSasaPoints* rsp = new ResSasaPoints();
	as->updateSasa(rsp);

	char xx[200];
	int n = pdb->getResList().size();

	cout << "res num: " << n << endl;

	ResBBRotamerLib* rotLib = new ResBBRotamerLib();
	AtomLib* atLib = new AtomLib();
	ResBBRotamer* rot;

	for(int i=0;i<n;i++){
		char ss = as->getSS(i);
		float sai = as->getSASA(i);
		int bbIndex = 0;
		ResBBRotamer* rot;
		if(i==0){
			rot = new ResBBRotamer(pdb->getResList()[i], atLib);
		}
		else if(pdb->getResList()[i-1]->contactToNeighbor(pdb->getResList()[i])){
			rot = new ResBBRotamer(pdb->getResList()[i-1], pdb->getResList()[i], atLib);
		}
		else {
			rot = new ResBBRotamer(pdb->getResList()[i], atLib);
		}
		bbIndex = rotLib->getRotamerIndex1W(rot);
		double dist = rot->distanceTo(rotLib->allRotLib1w[0][bbIndex]);


		sprintf(xx, "%c%-3s %c %5.3f %4d %5.3f", pdb->getResList()[i]->chainID, pdb->getResList()[i]->resID.c_str(), ss, sai, bbIndex, dist);
		out << string(xx) << endl;
		//out << rot->toString() << endl;
		//out << rotLib->allRotLib1w[0][bbIndex]->toString() << endl;
	}

	out.close();
	delete pdb;
	delete as;
	delete rsp;
}

