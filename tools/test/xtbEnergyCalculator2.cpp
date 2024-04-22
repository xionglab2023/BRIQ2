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


    string groupA = string(argv[1]);
    string groupB = string(argv[2]);
    string cmListFile = string(argv[3]);
    string xyzFile = string(argv[4]);
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

        sprintf(xx, "%s-%d.xyz", xyzFile.substr(0, xyzFile.length()-4).c_str(), id);
        id++;
        string outFile = string(xx);
        out.open(outFile, ios::out);

        int atomNum = rotA.atomNum + rotB.atomNum;

        out << atomNum << endl;
        out << endl;

        for(int i=0;i<rotA.atomNum;i++) {
            sprintf(xx, "%s %9.4f %9.4f %9.4f", conA.coords[i].x_, conA.coords[i].y_, conA.coords[i].z_, rotA.atomTypes[i].c_str());
            out << string(xx) << endl;
        }

        for(int i=0;i<rotB.atomNum;i++) {
            sprintf(xx, "%s %9.4f %9.4f %9.4f", conB.coords[i].x_, conB.coords[i].y_, conB.coords[i].z_, rotB.atomTypes[i].c_str());
            out << string(xx) << endl;
        }
        
        out << "$fix:" << endl;
        out << "    elements: C,O,N" << endl;
        out << "$end" << endl;
        out.close();
    }

}

