#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairComposition.h"
#include "model/RiboseRotamerLib.h"
#include "model/AtomLib.h"
#include "predNA/BRNode.h"
#include <time.h>
#include <stdio.h>
#include "predNA/BRFoldingTree.h"
#include "model/StructureModel.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/HbondEnergy.h"
#include "model/BasePairLib.h"

using namespace NSPmodel;
using namespace NSPforcefield;

int main(int argc, char** argv){
	clock_t start = clock();

    string augc = "AUGC";

    BasePairLib* bpLib = new BasePairLib("stat");
    ifstream file;
    ofstream out;
    string line;
    int idA, idB;
    double ene;
    int idCluster;

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            int num = bpLib->nnbBasePairNum[i*4+j];
            set<int> lowEnergySet;
            for(int k=0;k<num;k++){
                if(bpLib->nnbEnergyWithOxy[i*4+j][k] < -2.0){
                    lowEnergySet.insert(k);
                }
            }
            string pairType = augc.substr(i,1) + augc.substr(j,1);     
            string fileName = "/public/home/pengx/cpp/briqx/data/pairEneCG-raw/nnb/"+pairType+".ene";
            string outfile = "/public/home/pengx/cpp/briqx/data/pairEneCG-raw/nnb/"+pairType+".lowEnergyClusters";
            file.open(fileName.c_str(), ios::in);
            out.open(outfile.c_str(), ios::out);

            while(file >> idA >> idB >> ene >> idCluster){
                if(ene > -0.5) continue;
                if(!lowEnergySet.contains(idCluster)) continue;
                out << idA << " " << idB << " " << idCluster << endl;
            }

            file.close();
            out.close();

        }
    }


}