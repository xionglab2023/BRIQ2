/*
 * xtbEnergyCalculator.cpp
 *
 *  Created on: 2023��12��19��
 *      Author: nuc
 */

#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/PolarGroupLib.h"
#include "stdlib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

    /*
     *  usage: xtbEnergy $GROUP_TYPE_A $GROUP_TYPE_B $OI_TYPE $OI_INDEX ($OI_INDEX2) 
     */

    string groupA = string(argv[1]);
    string groupB = string(argv[2]);
    string cmListFile = string(argv[3]);
    string comFile = string(argv[4]);
    PolarGroupRotamer rotA(groupA);
    PolarGroupRotamer rotB(groupB);

    ifstream f;
    f.open(cmListFile, ios::in);
    string m;
    int id = 0;
    while(getline(f, m)){
        CsMove cm(m);
        LocalFrame csA;
        LocalFrame csB = csA.add(cm);
        PolarGroupConformer conA(&rotA, csA);
        PolarGroupConformer conB(&rotB, csB);

        ofstream out;
        char xx[100];

        sprintf(xx, "%s-%d.com", comFile.substr(0, comFile.length()-4).c_str(), id);
        id++;
        string outFile = string(xx);
        out.open(outFile, ios::out);


        out << "%NProc=32" << endl;
        out << "%MEM=2GB" << endl;
        out << "#n M062X/6-31+G(d,p) SP" << endl;
        out << endl;
        out << " Title" << endl;
        out << endl;
        out << "0 1" << endl;

        for(int i=0;i<rotA.atomNum;i++) {
            sprintf(xx, "%s %9.4f %9.4f %9.4f", rotA.atomTypes[i].c_str(), conA.coords[i].x_, conA.coords[i].y_, conA.coords[i].z_);
            out << string(xx) << endl;
        }

        for(int i=0;i<rotB.atomNum;i++) {
            sprintf(xx, "%s %9.4f %9.4f %9.4f", rotB.atomTypes[i].c_str(), conB.coords[i].x_, conB.coords[i].y_, conB.coords[i].z_);
            out << string(xx) << endl;
        }
        out << endl;
        out << endl;

        out.close();
    }

}


