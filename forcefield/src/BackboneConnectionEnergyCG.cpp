#include "forcefield/BackboneConnectionEnergyCG.h"
using namespace NSPforcefield;

BackboneConnectionEnergyCG::BackboneConnectionEnergyCG(ForceFieldPara* ffp){
    string path = NSPdataio::datapath()+"/";
    string fileName = path + "backboneCG/bbcg.ene";
    ifstream file;
    file.open(fileName.c_str(), ios::in);
    if(!file.is_open()){
        cout << "fail to open file: " << fileName << endl;
        exit(1);
    }
    int a,b,c,d;
    double e;

    while(file >> a >> b >> c >> d >> e){
        this->et[a*192000 + b*4800 + c*120 + d] = e;
    }
    file.close();
    this->para = ffp;
}

double BackboneConnectionEnergyCG::getEnergy(double d, double ang1, double ang2, double dihed){
    int iD = (int)((d-2.6)*10);
    int iAng1 = (int)((ang1-60)*0.33333);
    int iAng2 = (int)((ang2-60)*0.33333);
    int iDihed = (int)((dihed+180)*0.33333);


    double outOfBoundaryEnergy1 = 0;
    double outOfBoundaryEnergy2 = 0;
    double outOfBoundaryEnergy3 = 0;
    
    if(iD < 0) {
        outOfBoundaryEnergy1 = (2.6-d)*8.0;
        iD = 0;
    }
    else if(iD > 14){
        outOfBoundaryEnergy1 = (d-4.1)*6.0;
        iD = 14;
    }

    if(iAng1 < 0) {
        outOfBoundaryEnergy2 = (60 - ang1)*0.15;
        iAng1 = 0;
    }
    else if(iAng1 > 39){
        iAng1 = 39;
    }

    if(iAng2 < 0) {
        outOfBoundaryEnergy3 = (60 - ang2)*0.15;
        iAng2 = 0;
    }
    else if(iAng2 > 39){
        iAng2 = 39;
    }

    if(iDihed < 0) {
        iDihed = 0;
    }
    else if(iDihed > 119){
        iDihed = 120;
    }

    double e = this->et[iD*192000+iAng1*4800+iAng2*120+iDihed] + outOfBoundaryEnergy1 + outOfBoundaryEnergy2+ outOfBoundaryEnergy3;
    if(e > 0)
        return e*para->connectRescale * e;
    else    
        return this->et[iD*192000+iAng1*4800+iAng2*120+iDihed] + outOfBoundaryEnergy1 + outOfBoundaryEnergy2+ outOfBoundaryEnergy3;
}

BackboneConnectionEnergyCG::~BackboneConnectionEnergyCG(){

}