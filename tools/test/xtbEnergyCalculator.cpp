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

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

    /*
     *  usage: xtbEnergy $GROUP_TYPE_A $GROUP_TYPE_B $OI_TYPE $OI_INDEX ($OI_INDEX2) 
     */


    string groupA = string(argv[1]);
    string groupB = string(argv[2]);
    string xyzFile = string(argv[3]);
    string oiType = string(argv[4]);
    int index1, index2;
    CsMove cm;
    OrientationIndex x;
    if(oiType == "SP500") {
        index1 = atoi(argv[5]);
        cm = x.index500ToCsMove(index1);
    }      
    else if(oiType == "SP1000") {
        index1 = atoi(argv[5]);
        cm = x.index1000ToCsMove(index1);
    }
    else if(oiType == "SP2000") {
        index1 = atoi(argv[5]);
        index2 = atoi(argv[6]);
        cm = x.index2000ToCsMove(index1, index2);
    }

    LocalFrame csA;
    LocalFrame csB = csA.add(cm);

    PolarGroupRotamer rotA(groupA);
    PolarGroupRotamer rotB(groupB);
    
    PolarGroupConformer conA(&rotA, csA);
    PolarGroupConformer conB(&rotB, csB);

    /*
    double minD = 99.9;
    for(int i=0;i<rotA.atomNum;i++){
        for(int j=0;j<rotB.atomNum;j++){
            double d = conA.coords[i].distance(conB.coords[j]);
            if(d < 1.0) exit(0);
            if(d < minD){
                minD = d;
            }
        }
    }
    if(minD > 5.0) exit(0);
    */

    ofstream out;
    char xx[100];
    out.open(xyzFile, ios::out);

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

    out.close();


}


