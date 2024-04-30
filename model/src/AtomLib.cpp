/*
 * AtomLib.cpp
 *
 *  Created on: 2022��8��29��
 *      Author: pengx
 */

#include "model/AtomLib.h"

namespace NSPmodel {

AtomProperty::AtomProperty(const string& line) {
        vector<string> spt;
        string s;
        splitString(line, " ", &spt);
        this->uniqueID = atoi(spt[0].c_str());
        this->atomUniqueName = spt[1]+"-"+spt[2];

        this->vdwRadius = atof(spt[3].c_str());
        this->atomType = spt[5];
        this->atomTypeInt = -1;
        this->donorIndex = -1;
        this->acceptorIndex = -1;

        s = spt[4];
        if(s == "d"){
                this->isHDonor = true;
                this->isHAcceptor = false;
        }
        else if(s == "a")
        {
                this->isHDonor = false;
                this->isHAcceptor = true;
        }
        else if(s == "b")
        {
                this->isHDonor = true;
                this->isHAcceptor = true;
        }
        else
        {
                this->isHDonor = false;
                this->isHAcceptor = false;
        }

        if(s == "i")
        	this->isIon = true;
        else
        	this->isIon = false;
}


AtomProperty::~AtomProperty() {
        // TODO Auto-generated destructor stub
}

AtomLib::AtomLib() {
	string path = NSPdataio::datapath();
	string libFile = path+"atomLib/atLib.dat";
	ifstream file;

	file.open(libFile.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open full atom library file " << libFile << endl;
		exit(1);
	}

	for(int i=0;i<20;i++) {
		vector<int> ids;
		vector<string> atomNames;
		vector<string> scAtomNames;
		this->aaUniqueIDs.push_back(ids);
		this->aaAtomNames.push_back(atomNames);
		this->aaScAtomNames.push_back(scAtomNames);
	}

	for(int i=0;i<8;i++) {
		vector<int> ids;
		vector<string> atomNames;
		vector<string> scAtomNames;
		this->baseUniqueIDs.push_back(ids);
		this->baseAtomNames.push_back(atomNames);
		this->baseScAtomNames.push_back(scAtomNames);
	}



	for(int i=0;i<8;i++){
		string name = bn.intToString(i)+"-O2'";
		donorAtomMap[name] = 0;
	}

	donorAtomMap["SER-OG"] = 1;
	donorAtomMap["THR-OG1"] = 2;
	donorAtomMap["TYR-OH"] = 3;
	donorAtomMap["HIS-ND1"] = 4;
	donorAtomMap["HIS-NE2"] = 4;

	for(int i=0;i<20;i++){
		string name = rn.intToTri(i)+"-N";
		if(name == "PRO-N") continue;
		donorAtomMap[name] = 5;
	}
	donorAtomMap["ASN-ND2"] = 6;
	donorAtomMap["GLN-NE2"] = 6;
	donorAtomMap["ARG-NH1"] = 7;
	donorAtomMap["ARG-NH2"] = 7;
	donorAtomMap["ARG-NE"] = 8;
	donorAtomMap["LYS-NZ"] = 9;
	donorAtomMap["TRP-NE1"] = 10;
	donorAtomMap["A-N6"] = 11;
	donorAtomMap["DA-N6"] = 11;
	donorAtomMap["U-N3"] = 12;
	donorAtomMap["DT-N3"] = 12;
	donorAtomMap["G-N1"] = 13;
	donorAtomMap["DG-N1"] = 13;
	donorAtomMap["G-N2"] = 14;
	donorAtomMap["DG-N2"] = 14;
	donorAtomMap["C-N4"] = 15;
	donorAtomMap["DC-N4"] = 15;

	donorAtomMap["ATM-O1"] = 0;
	donorAtomMap["ATM-N1"] = 5;
	donorAtomMap["ATM-N2"] = 7;
	donorAtomMap["ATM-N3"] = 9;

	for(int i=0;i<8;i++){
		string name = bn.intToString(i)+"-O2'";
		acceptorAtomMap[name] = 0;
	}
	acceptorAtomMap["SER-OG"] = 1;
	acceptorAtomMap["THR-OG1"] = 2;
	acceptorAtomMap["TYR-OH"] = 3;
	acceptorAtomMap["HIS-ND1"] = 4;
	acceptorAtomMap["HIS-NE2"] = 4;
	for(int i=0;i<20;i++){
		string name = rn.intToTri(i)+"-O";
		acceptorAtomMap[name] = 5;
	}
	acceptorAtomMap["ASP-OD1"] = 6;
	acceptorAtomMap["ASP-OD2"] = 6;
	acceptorAtomMap["GLU-OE1"] = 7;
	acceptorAtomMap["GLU-OE2"] = 7;
	acceptorAtomMap["ASN-OD1"] = 8;
	acceptorAtomMap["GLN-OE1"] = 8;
	acceptorAtomMap["OP1"] = 9;
	acceptorAtomMap["OP2"] = 9;
	for(int i=0;i<8;i++){
		string name1 = bn.intToString(i)+"-OP1";
		string name2 = bn.intToString(i)+"-OP2";
		acceptorAtomMap[name1] = 9;
		acceptorAtomMap[name2] = 9;
	}
	acceptorAtomMap["A-N7"] = 10;
	acceptorAtomMap["DA-N7"] = 10;
	acceptorAtomMap["A-N1"] = 11;
	acceptorAtomMap["DA-N1"] = 11;
	acceptorAtomMap["A-N3"] = 12;
	acceptorAtomMap["DA-N3"] = 12;
	acceptorAtomMap["U-O2"] = 13;
	acceptorAtomMap["U-O4"] = 13;
	acceptorAtomMap["DT-O2"] = 13;
	acceptorAtomMap["DT-O4"] = 13;
	acceptorAtomMap["G-N7"] = 14;
	acceptorAtomMap["DG-N7"] = 14;
	acceptorAtomMap["G-O6"] = 15;
	acceptorAtomMap["DG-O6"] = 15;
	acceptorAtomMap["G-N3"] = 16;
	acceptorAtomMap["DG-N3"] = 16;
	acceptorAtomMap["C-O2"] = 17;
	acceptorAtomMap["DC-O2"] = 17;
	acceptorAtomMap["C-N3"] = 18;
	acceptorAtomMap["DC-N3"] = 18;

	acceptorAtomMap["ATM-OA"] = 5;
	acceptorAtomMap["ATM-O1"] = 0;
	acceptorAtomMap["ATM-N0A"] = 10;

	map<string,int>::const_iterator it;

	string s,prefix,uniqueName,atomName;
	int uniqueID = 0;
	int aaType, baseType;
	vector<string> spt;
	vector<string> spt2;
	char type = 'a';

	while(getline(file,s)) {
		prefix = s.substr(0,4);
		if(prefix == "BASE") {
			type = 'b';
			continue;
		}
		else if(prefix == "RESI") {
			type = 'a';
			continue;
		}
		else if(prefix == "GENE"){
			type = 'c';
			continue;
		}
		else {
			AtomProperty* ap = new AtomProperty(s);

			splitString(s, " ", &spt);
			uniqueID = atoi(spt[0].c_str());
			atomName = spt[2];
			if(type == 'a'){
				aaType = rn.triToInt(spt[1]);
				this->aaUniqueIDs[aaType].push_back(uniqueID);
				this->aaAtomNames[aaType].push_back(atomName);
				if(atomName != "N" && atomName != "CA" && atomName != "C" && atomName != "O")
					this->aaScAtomNames[aaType].push_back(spt[2]);
			}
			else if(type == 'b'){
				baseType = bn.stringToInt(spt[1]);
				this->baseUniqueIDs[baseType].push_back(uniqueID);
				this->baseAtomNames[baseType].push_back(atomName);
				if(atomName.length() == 2)
					this->baseScAtomNames[baseType].push_back(atomName);
			}
			else {
				generalAtomNames.push_back(atomName);
				generalAtomTypeToInt[atomName] = uniqueID - 334;
			}

			uniqueIDToDonorID[uniqueID] = -1;
			uniqueIDToAcceptorID[uniqueID] = -1;

			uniqueName = ap->atomUniqueName;
			it = donorAtomMap.find(uniqueName);
			if(it != donorAtomMap.end()){
				ap->donorIndex = it->second;
				uniqueIDToDonorID[uniqueID] = it->second;
			}

			it = acceptorAtomMap.find(uniqueName);
			if(it != acceptorAtomMap.end()){
				ap->acceptorIndex = it->second;
				uniqueIDToAcceptorID[uniqueID] = it->second;
			}

			this->uniqueNames.push_back(uniqueName);
			this->apList.push_back(ap);
			this->uniqueNameToUniqueID[uniqueName] = uniqueID;
		}
	}


	for(int i=0;i<apList.size();i++){
		apList[i]->atomTypeInt = generalAtomTypeToInt[apList[i]->atomType];
	}
}

bool AtomLib::atomDefined(const string& uniqueName) const
{
	map<string,int>::const_iterator it;
	it = uniqueNameToUniqueID.find(uniqueName);
	if(it != uniqueNameToUniqueID.end())
		return true;
	else
		return false;
}

int AtomLib::uniqueNameToID(const string& uniqueName) const
{
	map<string,int>::const_iterator it;
	it = uniqueNameToUniqueID.find(uniqueName);
	if(it != uniqueNameToUniqueID.end())
		return it->second;
	else
	{
		cerr << "Warning: unknown uniqueName " << uniqueName << endl;
		return -1;
	}
}

string AtomLib::uniqueIDToName(int uniqueID) const
{
	if(uniqueID < 0 || uniqueID >= (int)this->apList.size())
	{
		cerr << "invalid uniqueID: " << uniqueID << endl;
		exit(0);
	}
	return apList[uniqueID]->atomUniqueName;
}

void AtomLib::getAminoAcidAtomIDs(int aaName, vector<int>& aminoAcidAtomIDs)
{
	if(aaName < 0 || aaName >= (int)this->aaUniqueIDs.size())
	{
		cerr << "invalid atom type: " << aaName << endl;
		return;
	}
	aminoAcidAtomIDs.clear();
	for(int i=0;i<this->aaUniqueIDs[aaName].size();i++)
		aminoAcidAtomIDs.push_back(this->aaUniqueIDs[aaName][i]);
}

void AtomLib::getAminoAcidAtomNames(int aaType, vector<string>& aminoAcidNames)
{
	if(aaType < 0 || aaType >= (int)this->aaUniqueIDs.size())
	{
		cerr << "invalid atom type: " << aaType << endl;
		return;
	}
	aminoAcidNames.clear();
	for(int i=0;i<this->aaAtomNames[aaType].size();i++)
		aminoAcidNames.push_back(this->aaAtomNames[aaType][i]);
}

void AtomLib::getAminoAcidSidechainAtomNames(int aaType, vector<string>& aminoAcidSidechainAtoms)
{
	if(aaType < 0 || aaType >= (int)this->aaUniqueIDs.size())
	{
		cerr << "invalid atom type: " << aaType << endl;
		return;
	}
	aminoAcidSidechainAtoms.clear();
	for(int i=0;i<this->aaScAtomNames[aaType].size();i++)
		aminoAcidSidechainAtoms.push_back(this->aaScAtomNames[aaType][i]);
}

AtomProperty* AtomLib::getAtomProperty(int uniqueID){
	if(uniqueID <0 || uniqueID > apList.size()){
		cerr << "invalid atom unique ID: " << uniqueID << endl;
		return NULL;
	}
	return this->apList[uniqueID];
}

vector<int> AtomLib::getAASidechainUniqueIDs(int type){
	vector<int> idList;
	int n = aaUniqueIDs[type].size();
	for(int i=4;i<n;i++){
		idList.push_back(aaUniqueIDs[type][i]);
	}
	return idList;
}

void AtomLib::getRnaAtomNames(int type, vector<string>& rnaAtomNames) const{
	if(type < 0 || type >= (int)this->baseUniqueIDs.size())
	{
		cerr << "invalid atom type: " << type << endl;
		return;
	}
	rnaAtomNames.clear();
	for(int i=0;i<this->baseAtomNames[type].size();i++)
		rnaAtomNames.push_back(this->baseAtomNames[type][i]);
}

void AtomLib::getRnaSidechainAtoms(int type,vector<string>& rnaScNames) const{
	if(type < 0 || type >= (int)this->baseUniqueIDs.size())
	{
		cerr << "invalid atom type: " << type << endl;
		return;
	}
	rnaScNames.clear();
	for(int i=0;i<this->baseScAtomNames[type].size();i++)
		rnaScNames.push_back(this->baseScAtomNames[type][i]);
}

AtomLib::~AtomLib() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<apList.size();i++){
		delete apList[i];
	}

/*
	for(int i=0;i<8;i++){
		delete this->baseUniqueIDs[i];
		delete this->baseAtomNames[i];
		delete this->baseScAtomNames[i];
	}

	for(int i=0;i<20;i++){
		delete this->aaUniqueIDs[i];
		delete this->aaAtomNames[i];
		delete this->aaScAtomNames[i];
	}
*/

}

} /* namespace NSPmodel */
