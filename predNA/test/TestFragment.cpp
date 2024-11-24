/*
 * TestFragment.cpp
 *
 */



#include "predNA/FragmentLibrary.h"

using namespace NSPpredNA;
using namespace NSPgeometry;


int main(int argc, char** argv){


	RotamerLib* rotLib = new RotamerLib();
	FragmentLibrary* lib = new FragmentLibrary(rotLib);


	for(int j=0;j<100;j++){
		char xx[200];
		sprintf(xx, "/export/home/s2982206/aspen/rnaModeling/fragLib/f2/loopNb/pdb/loopNb-1-%d.pdb",j);
		string outfile = string(xx);
		F2Fragment* f = lib->loopNb[1]->fragListLevel2[j];
		f->printPDB(outfile);
	}

	for(int j=0;j<100;j++){
		char xx[200];
		sprintf(xx, "/export/home/s2982206/aspen/rnaModeling/fragLib/f2/loopNb/pdb/wcPair-1-%d.pdb",j);
		string outfile = string(xx);
		F2Fragment* f = lib->wcPair[1]->fragListLevel2[j];
		f->printPDB(outfile);
	}

	for(int j=0;j<100;j++){
		char xx[200];
		sprintf(xx, "/export/home/s2982206/aspen/rnaModeling/fragLib/f2/loopNb/pdb/wcNb-1-%d.pdb",j);
		string outfile = string(xx);
		F2Fragment* f = lib->wcNb[1]->fragListLevel2[j];
		f->printPDB(outfile);
	}

	for(int j=0;j<100;j++){
		char xx[200];
		sprintf(xx, "/export/home/s2982206/aspen/rnaModeling/fragLib/f2/loopNb/pdb/nwcPair-1-%d.pdb",j);
		string outfile = string(xx);
		F2Fragment* f = lib->nwcPair[1]->fragListLevel2[j];
		f->printPDB(outfile);
	}



}
