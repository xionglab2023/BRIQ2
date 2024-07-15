#include "model/LigandLib.h"

namespace NSPmodel {

LigandInfo::LigandInfo(const string& ligName, int atomNum, int polarAtomNum){
    this->ligandName = ligName;
    this->atomNum = atomNum;
    this->polarAtomNum = polarAtomNum;
    this->uniqueIDs = new int[atomNum];
    this->polarAtomIndex = new int[polarAtomNum];
    this->supportAtomIndex1 = new int[polarAtomNum];
    this->supportAtomIndex2 = new int[polarAtomNum];
    this->supportType = new int[polarAtomNum];

    curPolarAtomNum = 0;
}

void LigandInfo::printInfo(){
    cout << "Ligand name: " << ligandName << endl;
    cout << "unique names: " << endl;
    for(int i=0;i<atomNum;i++){
        cout << i << " " << uniqueIDs[i] << endl;
    }

    cout << "polarAtoms: " << endl;
    for(int i=0;i<polarAtomNum;i++){
        cout << i << " " << polarAtomIndex[i] << endl;
    }
}

LigandInfo::~LigandInfo(){
    delete [] uniqueIDs;
    delete [] polarAtomIndex;
    delete [] supportAtomIndex1;
    delete [] supportAtomIndex2;
    delete [] supportType;
}

void LigandInfo::addAtom(const string& line, AtomLib* atLib){
    vector<string> spt;
    splitString(line, " ", &spt);
    int atomID = atoi(spt[0].c_str());
    string atomName = spt[1];
    this->atomNames.push_back(atomName);
    string uniqueName = "ATM-"+spt[2];
    this->uniqueIDs[atomID] = atLib->uniqueNameToID(uniqueName);

    string donorAcceptorType = spt[3];
    if(donorAcceptorType != "n"){
        //polar atom
        this->polarAtomIndex[curPolarAtomNum] = atomID;

        string supType = spt[4];
        if(supType == "t")
            this->supportType[curPolarAtomNum] = 0;
        else 
            this->supportType[curPolarAtomNum] = 1;

        int sup1 = atoi(spt[5].c_str());
        int sup2 = atoi(spt[6].c_str());

        this->supportAtomIndex1[curPolarAtomNum] = sup1;
        this->supportAtomIndex2[curPolarAtomNum] = sup2;

        curPolarAtomNum++;
    }
}

LigandLib::LigandLib(){
    string path = NSPdataio::datapath();
    string ligFile = path + "atomLib/ligand.dat";
    ifstream file;
    file.open(ligFile, ios::in);
    if(!file.is_open()){
        cout << "fail to open file: " << ligFile << endl;
        exit(0);
    }
    string line;
    vector<string> spt;
    AtomLib* atLib = new AtomLib();

    int curLigandID = -1;

    while(getline(file, line)){
        if(line.length() < 5) continue;
        if(line.substr(0,3) == "LIG"){
            splitString(line, " ", &spt);
            LigandInfo* lig = new LigandInfo(spt[1], atoi(spt[2].c_str()), atoi(spt[3].c_str()));
            this->ligandList.push_back(lig);
            curLigandID ++;
        }
        else {
            this->ligandList[curLigandID]->addAtom(line, atLib);
        }
    }

    for(int i=0;i<ligandList.size();i++){
        this->ligMap[ligandList[i]->ligandName] = ligandList[i];
    }

    delete atLib;



}

void LigandLib::printInfo(){
    for(int i=0;i<ligandList.size();i++){
        ligandList[i]->printInfo();
    }

}


LigandLib::~LigandLib(){
    int n = ligandList.size();
    for(int i=0;i<n;i++){
        delete ligandList[i];
    }
}

}