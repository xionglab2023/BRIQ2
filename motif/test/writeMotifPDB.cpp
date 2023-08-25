/*
 * writeMotifPDB.cpp
 *
 *  Created on: 2023年8月25日
 *      Author: nuc
 */


/*
 * writeMotifPDB.cpp
 *
 *  Created on: 2020年12月2日
 *      Author: pengx
 */



#include "motif/RNAGraph.h"
#include "model/StructureModel.h"
#include <stdlib.h>
#include <iostream>
#include <time.h>

using namespace NSPmotif;
using namespace std;
using namespace NSPmodel;

int main(int argc, char** argv){

	ifstream file;
	file.open(argv[1], ios::in);
	string line;
	vector<string> spt;
	AtomLib* atLib = new AtomLib();
	while(getline(file, line)){
		string pdbID = line.substr(0,4);
		RNAPDB* pdb = new RNAPDB("/lustre/home/pengx/lib/rna/pdb09/"+line);
		vector<RNABase*> baseList = pdb->getValidBaseList(atLib);

		ifstream f2;
		string mtfFile = "/lustre/home/pengx/motif/mtf/"+pdbID+".mtf";
		string f2Line;
		vector<string> lines;
		f2.open(mtfFile.c_str(), ios::in);
		if(f2.is_open()){
			while(getline(f2, f2Line)){
				lines.push_back(f2Line);
			}
		}
		int mtfID = 0;
		for(int i=0;i<lines.size();i+=2){
			splitString(lines[i], " ", &spt);
			string mtfSeq = lines[i+1];
			string lastChainID = "X";
			if(mtfSeq.length() != baseList.size()){
				cout << "invalid pdb " << pdbID << endl;
				continue;
			}
			double score = atof(spt[3].c_str());

			if(score < 1.0){
				RNAPDB* rp = new RNAPDB();
				RNAChain* rc = NULL;
				for(int k=0;k<mtfSeq.length();k++){
					string curChainID = mtfSeq.substr(k,1);
					if(mtfSeq[k] == '.') continue;
					if(curChainID == lastChainID){
						baseList[k]->chainID = curChainID;
						rc->addBase(baseList[k]);
						//cout << rc->getBaseList().size() << endl;
					}
					else {
						//cout << "add base " << k << endl;
						baseList[k]->chainID = curChainID;
						lastChainID = curChainID;
						rc = new RNAChain();
						rc->setChainID(lastChainID);

						rp->addChain(rc);
						rc->addBase(baseList[k]);
						//cout << rc->getBaseList().size() << endl;
					}
				}
				int baseNum = rp->getBaseList().size();
				//cout  << pdbID << " " << baseNum << endl;
				ofstream out;
				char xx[200];
				sprintf(xx, "/lustre/home/pengx/motif/mpdb/%s-%d.mtf.pdb", pdbID.c_str(), mtfID);
				mtfID++;
				string outFile = string(xx);
				out.open(outFile.c_str(), ios::out);
				rp->printPDBFormat(out);
				out.close();
			}
		}
		delete pdb;


	}
}



