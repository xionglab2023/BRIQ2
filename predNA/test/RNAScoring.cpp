/*
 * RNAScoring.cpp
 *
 */


#include "geometry/localframe.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "forcefield/RnaEnergyTable.h"
#include "model/StructureModel.h"
#include "predNA/BRFoldingTree.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;

void scorePDBFile(const string& pdbFile, const string& output, RnaEnergyTable& et, RiboseRotamerLib& rotLib){
/*

	AtomLib* atLib = new AtomLib();
	{
		ofstream of;
		of.open(output.c_str(), ios::out);

		char xx[200];
		RNAPDB pdb(pdbFile, "xxxx");
		vector<RNABase*> baseList = pdb.getValidBaseList(atLib);
		int seqLen = baseList.size();
		if(seqLen == 0) return;

		bool connectToDownstream[seqLen];

		for(int i=0;i<seqLen;i++){
			if(i<seqLen-1 && baseList[i]->connectToNeighbor(*baseList[i+1]))
				connectToDownstream[i] = true;
			else
			{
				//cout << "break " << i << endl;
				connectToDownstream[i] = false;
			}
		}


		vector<BRNode*> nodes;
		for(int i=0;i<seqLen;i++){
			RiboseRotamer* rot = rotLib.getNearestRotamer(baseList[i]);
			BRNode* node = new BRNode(baseList[i], rot);
			node->connectToNeighbor = connectToDownstream[i];
			node->rot = rotLib.getNearestRotamer(baseList[i]);
			nodes.push_back(node);
		}

		for(int i=0;i<seqLen-1;i++){
			BRNode* nodeA = nodes[i];
			BRNode* nodeB = nodes[i+1];
			if(!connectToDownstream[i]) continue;
			PhophateGroupLocal pho = et.pb->getPhoLocal(nodeA->cs2, nodeB->riboAtomCoords, nodeA->rot->improper, nodeB->rot->improper);
			nodeA->pho = PhophateGroup(pho, nodeA->cs2);
			nodeA->phoLocal = pho;
		}


		double bpEnergy1=0;
		double bpEnergy2=0;
		double bpEnergy3=0;

		double baseClash = 0;
		double oxyE = 0;

		double phoConnectionEnergy = 0;
		double rotEnergy = 0;

		for(int i=0;i<seqLen;i++){
			rotEnergy += nodes[i]->rot->energy;
			phoConnectionEnergy += nodes[i]->pho.e;

		}

		for(int i=0;i<seqLen;i++){
			sprintf(xx, "rot: %3d %8.3f", i,  nodes[i]->rot->energy);
			of << string(xx) << endl;
		}

		for(int i=0;i<seqLen;i++){
			sprintf(xx, "pho: %3d %8.3f", i, nodes[i]->pho.e);
			of << string(xx) << endl;
		}

		double shift = 0;

		for(int a=0;a<seqLen;a++){
			BRNode* nodeA = nodes[a];
			PhophateGroup phoA = nodeA->pho;
			for(int b=a+1;b<seqLen;b++){
				BRNode* nodeB = nodes[b];
				PhophateGroup phoB = nodeB->pho;
				int sep = 3;
				if(b == a+1 && connectToDownstream[a])
					sep = 1;
				if(b == a+2 && connectToDownstream[a] && connectToDownstream[a+1])
					sep = 2;

				int sep2 = -sep;
				if(sep2 == -3) sep2 = 3;

				int i,j;


				double bb = 0;
				double bo = 0;
				double bp = 0;

				double minDD = 990009.9;



				double e1,e2,e3,e4,e5,e6,e7,e8;
				e1 = getBaseRiboseEnergy(nodeA, nodeB, sep, &et, false);
				e2 = getBaseRiboseEnergy(nodeB, nodeA, sep2, &et, false);
				e3 = getBasePhoEnergy(nodeA, nodeB, sep, &et,false);
				e4 = getBasePhoEnergy(nodeB, nodeA, sep2, &et, false);
				e5 = getRiboseRiboseEnergy(nodeA, nodeB, sep, &et, false);
				e6 = getRibosePhoEnergy(nodeA, nodeB, sep, &et, false);
				e7 = getRibosePhoEnergy(nodeB, nodeA, sep2, &et, false);
				e8 = getPhoPhoEnergy(nodeA, nodeB, sep, &et, false);

				if(e1 < 0) bo += e1;
				if(e2 < 0) bo += e2;
				if(e3 < 0) bp += e3;
				if(e4 < 0) bp += e4;


				if(sep == 1){
					//bp energy
					bb = getBaseBaseEnergy(nodeA, nodeB, sep, &et, false);
				}
				else if(sep == 2){
					//bp energy
					bb = getBaseBaseEnergy(nodeA, nodeB, sep, &et, false);
					}
				else if(sep == 3){
					//bp energy
					bb = getBaseBaseEnergy(nodeA, nodeB, sep, &et, false);
				}


				if(abs(bb+bo+bp) > 0.001){
					sprintf(xx, "pair-%d: %3d %3d %8.3f %8.3f %8.3f", sep, a, b, bb, bo, bp);
					of << string(xx) << endl;
				}
			}
		}
		of.close();

		for(int i=0;i<seqLen;i++){
			delete nodes[i];
		}
	}
	delete atLib;

*/
}

int main(int argc, char** argv){

	//without clash energy

	/*usage:
		rnaScoring -list inputList outputList
		rnaScoring -single inputFile outputFile
	*/
	if(argc != 4){
		cout << "Usage: " << endl;
		cout << "rnaScoring -list inputList outputList" << endl;
		cout << "rnaScoring -single inputFile outputFile" << endl;
	}

	string tag = string(argv[1]);
	if(tag == "-list" || tag == "-l"){
		RnaEnergyTable et;
		RiboseRotamerLib rotLib;


		string fileList = string(argv[2]);
		string outputList = string(argv[3]);

		ifstream file;
		file.open(fileList, ios::in);
		string s;

		vector<string> inputFiles;
		vector<string> outFiles;


		while(file >> s){
			inputFiles.push_back(s);
		}
		file.close();

		file.open(outputList, ios::in);
		while(file >> s){
			outFiles.push_back(s);
		}
		file.close();

		if(inputFiles.size() != outFiles.size()){
			cout << "input output file number not equal: " << inputFiles.size() << " " << outFiles.size() << endl;
		}

		for(int k=0;k<inputFiles.size();k++)	{
			scorePDBFile(inputFiles[k], outFiles[k], et, rotLib);
		}
	}
	else if(tag == "-single" || tag == "-s"){

		string pdbFile = string(argv[2]);
		string outFile = string(argv[3]);
		RnaEnergyTable et;
		RiboseRotamerLib rotLib;
		scorePDBFile(pdbFile, outFile, et, rotLib);
	}
	else {
		cout << "Usage: " << endl;
		cout << "rnaScoring -list inputList outputList" << endl;
		cout << "rnaScoring -single inputFile outputFile" << endl;
		cout << argv[1]	 << endl;
		cout << argc << endl;
	}


}



