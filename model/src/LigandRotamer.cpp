#include "model/LigandRotamer.h"

namespace NSPmodel {

LigandRotamer::LigandRotamer(LigandInfo* ligInfo, const string& rotLine){
    this->ligInfo = ligInfo;
    this->atomNum = ligInfo->atomNum;
    this->polarAtomNum = ligInfo->polarAtomNum;
    this->localTList = new XYZ[ligInfo->atomNum];
    this->polarCmList = new CsMove[ligInfo->polarAtomNum];
    
    vector<string> spt;
    splitString(rotLine, " " , &spt);

    if(spt.size() != atomNum*3+2){
        cout << "rotamer line error: " << endl;
        cout << rotLine << endl;
        exit(0);
    }
    
    vector<XYZ> coordList;
    for(int i=0;i<atomNum;i++){
        coordList.push_back(XYZ(atof(spt[i*3+1].c_str()), atof(spt[i*3+2].c_str()), atof(spt[i*3+3].c_str())));
    }

    LocalFrame cs1 = LocalFrame(coordList[0], coordList[1], coordList[2]);
    for(int i=0;i<atomNum;i++){
        localTList[i] = global2local(cs1, coordList[i]);
    }

    for(int i=0;i<polarAtomNum;i++){
        int pIndex = ligInfo->polarAtomIndex[i];
        int sup1 = ligInfo->supportAtomIndex1[i];
        int sup2 = ligInfo->supportAtomIndex2[i];
        int supType = ligInfo->supportType[i];
        
        LocalFrame csP;
        if(supType == '0'){
            csP = LocalFrame(coordList[sup2], coordList[sup1], coordList[pIndex]);
        }
        else {
            csP = generateLocalFrameResidueStyle(coordList[sup1], coordList[pIndex], coordList[sup2]);
        }

        this->polarCmList[i] = csP - cs1;
    }

}

LigandRotamer::~LigandRotamer(){
    delete [] localTList;
    delete [] polarCmList;
}


LigandRotamerLib::LigandRotamerLib(LigandInfo* ligInfo){
    string path = NSPdataio::datapath();
    string ligandName = ligInfo->ligandName;
    string ligFile = path + "ligandLib/" + ligandName + ".rot";
    ifstream file;
    file.open(ligFile.c_str(), ios::in);
    if(!file.is_open()){
        cout << "fail to open file: " << ligFile << endl;
    }
    string line;
    while(getline(file, line)){
        if(line.length() < 2) continue;
        if(line[0] == '#') continue;
        LigandRotamer* rot = new LigandRotamer(ligInfo, line);
        this->rots.push_back(rot);
    }
}

LigandRotamerLib::~LigandRotamerLib(){
    for(int i=0;i<rots.size();i++){
        delete rots[i];
    }
}

LigandConformer::LigandConformer(LigandRotamer* rot, LocalFrame& cs){
    this->rot = rot;
    this->coords = new XYZ[rot->atomNum];
    this->cs1 = cs;
    this->csPolarList = new LocalFrame[rot->polarAtomNum];

}

LigandConformer::~LigandConformer(){
    delete [] coords;
    delete [] csPolarList;
}

void LigandConformer::copyValueFrom(LigandConformer* other){
    if(this->rot->atomNum != other->rot->atomNum){
        cout << "different ligand type could not copy value" << endl;
        exit(0);
    }
    else if(this->rot->polarAtomNum != other->rot->polarAtomNum){
        cout << "different ligand type could not copy value" << endl;
        exit(0);
    }

    this->rot = other->rot;
    updateCoords(other->cs1);
}

void LigandConformer::updateCoords(LocalFrame& cs){
    this->cs1 = cs;
    for(int i=0;i<rot->atomNum;i++){
        this->coords[i] = local2global(cs, this->rot->localTList[i]);
    }
    for(int i=0;i<rot->polarAtomNum;i++){
        this->csPolarList[i] = cs + rot->polarCmList[i];
    }
}

void LigandConformer::updateRotamer(LigandRotamer* rot){
    this->rot = rot;
        for(int i=0;i<rot->atomNum;i++){
        this->coords[i] = local2global(cs1, this->rot->localTList[i]);
    }
    for(int i=0;i<rot->polarAtomNum;i++){
        this->csPolarList[i] = cs1 + rot->polarCmList[i];
    }
}

void LigandConformer::updateRotamerAndLocalFrame(LigandRotamer* rot, LocalFrame& cs){
    this->cs1 = cs;
    this->rot = rot;
        for(int i=0;i<rot->atomNum;i++){
        this->coords[i] = local2global(cs1, this->rot->localTList[i]);
    }
    for(int i=0;i<rot->polarAtomNum;i++){
        this->csPolarList[i] = cs1 + rot->polarCmList[i];
    }
}


}
