/*
 * TestRiboseRotamer.cpp
 *
 *  Created on: 2023��11��23��
 *      Author: nuc
 */



#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include "model/RiboseRotamerLib.h"
#include "model/AtomLib.h"
#include "model/BasePairLib.h"
#include "predNA/BRNode.h"


using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;

int main(int argc, char** argv){

	RiboseRotamerLib* rotLib = new RiboseRotamerLib();

	ofstream out;
	char x[200];
	for(int type=0;type<8;type++){
		string augc = "AUGCaugc";
		out.open("/public/home/pengx/briqx/rnaRiboseRotamer/plot/"+augc.substr(type, 1)+".plot", ios::out);
		for(int i=0;i<1500;i++){
			RiboseRotamer* rot = rotLib->rotLib[type][i];
			double imp = rot->improper;
			double chi = rot->chi;
			sprintf(x, "%-8.3f %-7.3f", chi, imp);
			out << string(x) << endl;

		}
		out.close();
	}

	delete rotLib;

}
