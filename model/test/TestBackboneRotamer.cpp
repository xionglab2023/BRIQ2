/*
 * TestBackboneRotamer.cpp
 *
 *  Created on: 2022Äê3ÔÂ23ÈÕ
 *      Author: pengx
 */


#include "model/RNABaseLib.h"
#include "geometry/RMSD.h"
#include "stdlib.h"
#include <iostream>
#include "model/StructureModel.h"
#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResScRotamer.h"
#include "model/ResScRotamerLib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

	ResBBRotamerLib* rotLib = new ResBBRotamerLib();

	cout << "testNeighbor 1k" << endl;
	for(int i=0;i<10;i+=20){
		int id = rotLib->neighborLib1k[i][0];
		cout << "id: " << id << endl;
		ResBBRotamer* rot = rotLib->allRotLib1k[0][rotLib->neighborLib1k[i][0]];
		cout << "find rot" << rot->index1K << " " << rot->index1W << endl;
		for(int j=0;j<20;j++){
			double d= rot->distanceTo(rotLib->allRotLib1k[0][rotLib->neighborLib1k[i][j]]);
			printf("%5.3f ", d);
		}
		printf("\n");
	}

	cout << "testNeighbor 1w" << endl;
	for(int i=0;i<10;i+=20) {
		int id = rotLib->neighborLib1w[i][0];
		cout << "id: " << id << endl;
		ResBBRotamer* rot = rotLib->allRotLib1w[0][rotLib->neighborLib1w[i][0]];
		for(int j=0;j<12;j++){
			double d= rot->distanceTo(rotLib->allRotLib1w[0][rotLib->neighborLib1w[i][j]]);
			printf("%5.3f ", d);
		}
		printf("\n");
	}

	ResName rn;
	ResScRotamerLib* scLib = new ResScRotamerLib();
	vector<string> triList;
	triList.push_back("GLU");
	triList.push_back("GLN");
	triList.push_back("ARG");
	triList.push_back("LYS");
	triList.push_back("MET");

	for(int i=0;i<triList.size();i++){

		int type = rn.triToInt(triList[i]);
		int rotNum = scLib->rotClusterUnique[type][223][0];
		cout << triList[i] << " " << rotNum << endl;
		for(int j=0;j<rotNum;j++) {
			ResScRotamer* rot = scLib->rotList[type][scLib->rotClusterUnique[type][223][j+1]];
			double e = scLib->getEnergy(223, rot);
			printf("id: %3d ene: %7.3f\n", scLib->rotClusterUnique[type][223][j+1], e);
		}
	}
}



