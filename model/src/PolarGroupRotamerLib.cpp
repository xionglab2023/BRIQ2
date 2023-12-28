#include <model/PolarGroupLib.h>

namespace NSPmodel{


PolarGroupRotamer::PolarGroupRotamer(const string& groupName){


    string path = NSPdataio::datapath();
    string fileName = NSPdataio::datapath() + "polarGroup/" + groupName + ".pdb";
    ifstream f;
    f.open(fileName.c_str(), ios::in);
    if(!f.is_open()) {
        cout << "fail to open polar group file: " << fileName << endl;
        exit(0);
    }


    string s;
	while(getline(f, s)){
        if(!s.starts_with("ATOM")) continue;
        Atom at(s);
        this->atomTypes.push_back(at.getType());
        this->localCoords.push_back(at.getCoord());
	}
	f.close();

    this->atomNum = this->atomTypes.size();
    this->groupName = groupName;
   
}

PolarGroupRotamer::~PolarGroupRotamer(){

}

/// @brief polar group conformer
/// @param rot 
/// @param cs 
PolarGroupConformer::PolarGroupConformer(PolarGroupRotamer* rot, LocalFrame& cs){
    this->groupName = rot->groupName;
    this->rot = rot;
    this->coords.clear();
     this->cs = cs;
    for(int i=0;i<rot->atomNum;i++){
       XYZ coord = local2global(cs, rot->localCoords[i]);
        this->coords.push_back(coord);
    }
}

PolarGroupConformer::~PolarGroupConformer(){

}

PolarGroupRotamerLib::PolarGroupRotamerLib(){
    this->names.push_back("A");
    this->names.push_back("U");
    this->names.push_back("G");
    this->names.push_back("C");
    this->names.push_back("T");

    for(int i=0;i<this->names.size();i++){
        PolarGroupRotamer* rot = new PolarGroupRotamer(this->names[i]);
        rotMap[this->names[i]] = rot;
    }
}

PolarGroupRotamerLib::~PolarGroupRotamerLib(){
    for(int i=0;i<names.size();i++){
        delete rotMap[this->names[i]];
    }
}

}