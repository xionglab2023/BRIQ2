/*
 * EdgeInformation.cpp
 *
 *  Created on: 2023��11��30��
 *      Author: nuc
 */

#include <predNA/EdgeInformation.h>

namespace NSPpredNA {

EdgeInformation::EdgeInformation(int sep, int typeA, int typeB, BasePairLib* pairLib){
	this->sep = sep;
	this->typeA = typeA%4;
	this->typeB = typeB%4;


	if(sep == 1)
		this->totalClusterNum = pairLib->nbBasePairNum[this->typeA*4+this->typeB];
	else if(sep == -1)
		this->totalClusterNum = pairLib->nbBasePairNum[this->typeB*4+this->typeA];
	else if(sep == 0)
		this->totalClusterNum = 1;
	else
		this->totalClusterNum = pairLib->nnbBasePairNum[this->typeA*4+this->typeB];


	this->pCluster = new double[totalClusterNum];

	for(int i=0;i<totalClusterNum;i++){
		pCluster[i] = 1.0/totalClusterNum;
	}
	this->validClusterNum = totalClusterNum;

	if(sep == 1 || sep == -1)
		this->pContact = 1.0;
	else if(sep == 0)
		this->pContact = 1.0;
	else
		this->pContact = 0.01;


	this->weight = 0.0;

	for(int i=0;i<totalClusterNum;i++){
		double e = pairLib->getEnergy(i, typeA, typeB, sep);
		this->weight += pCluster[i]*e;
	}

	this->weight = weight*pContact;
	this->fixed = false;
	this->ssSepKey = "UNK";
}


EdgeInformation::~EdgeInformation() {
	delete [] pCluster;
}

void EdgeInformation::updatePCluster(double* pList, double pContact, BasePairLib* pairLib){
	double sum = 0.0;
	for(int i=0;i<totalClusterNum;i++){
			sum += pList[i];
	}
	if(sum == 0.0) sum = 1.0;

	this->weight = 0;
	this->validClusterNum = 0;
	for(int i=0;i<totalClusterNum;i++){
		this->pCluster[i] = pList[i]/sum;
		if(pList[i] > 0) 
			this->validClusterNum++;
		double e = pairLib->getEnergy(i, typeA, typeB, sep);
		this->weight += pCluster[i]*e;
	}

	this->pContact = pContact;
	this->weight = weight*pContact;
	this->ssSepKey = "UNK";
}

void EdgeInformation::setUniqueCluster(int clusterID, BasePairLib* pairLib){

	this->ssSepKey = "UNK";
	this->validClusterNum = 1;
	if(clusterID < 0) {
		for(int i=0;i<totalClusterNum;i++){
			this->pCluster[i] = 1.0/totalClusterNum;
		}
		this->pContact = 1.0;
		this->weight = 99.9;
	}
	else if(clusterID < totalClusterNum) {
		for(int i=0;i<totalClusterNum;i++){
			this->pCluster[i] = 0.0;
		}
		this->pCluster[clusterID] = 1.0;
		this->pContact = 1.0;
		this->weight = pairLib->getEnergy(clusterID, typeA, typeB, sep);
	}
	else {
		cout << "invalid clusterID: "<< clusterID <<  " typeA: " << typeA << " typeB: " << typeB << " sep: " << sep << " totNum: " << totalClusterNum << endl;
		exit(0);
	}
}

void EdgeInformation::setToLibPCluster(const string& ssSepType, EdgeInformationLib* eiLib){
	this->ssSepKey = ssSepType;
	map<string, EdgeInformation*>::iterator it;
	it = eiLib->eiMap.find(ssSepType);
	if(it != eiLib->eiMap.end()){
		this->copyClusterFrom(it->second);
	}
	else {
		cout << "invalid ssSepType: " << ssSepType << endl;
		exit(0);
	}
}

void EdgeInformation::setClusterList(vector<int>& clusterList, vector<double>& pList, BasePairLib* pairLib){
	
	double e, pSum;
	if(clusterList.size() == 0) {
		for(int i=0;i<totalClusterNum;i++){
			this->pCluster[i] = 1.0/totalClusterNum;
		}
		this->validClusterNum = totalClusterNum;
		this->pContact = 1.0;
		this->weight = 999.9;
		if(abs(this->sep) == 1 )
			this->weight = 99.9;
		this->ssSepKey = "UNK";	
	}
	else {
		this->validClusterNum = clusterList.size();
		
		if(clusterList.size() == 1)
			this->ssSepKey = "UNK";
		else 
			this->ssSepKey = "UNK";

		for(int i=0;i<totalClusterNum;i++){
			this->pCluster[i] = 0.0;
		}

		vector<double> eneList;
		pSum = 0;
		for(int i=0;i<clusterList.size();i++){
			e = pairLib->getEnergy(clusterList[i], typeA, typeB, sep);
			eneList.push_back(e);
			pSum += exp(-e);
		}

		this->weight = 0;
		for(int i=0;i<clusterList.size();i++){
			e = eneList[i];
			this->pCluster[clusterList[i]] = pList[i];
			this->pContact = 1.0;
			this->weight += eneList[i] * pList[i];
		}
	}
}

EdgeInformationLib::EdgeInformationLib(BasePairLib* pairLib){
	string path = NSPdataio::datapath();

	ifstream file1, file2;
	string s;
	
	string fileName = path+"edgeInfoLib/list";
	file1.open(fileName.c_str(), ios::in);
	if(!file1.is_open()){
		cout << "fail to open list file : " << fileName << endl;
		exit(1);
	}


	vector<int> sepList;
	vector<int> typeAList;
	vector<int> typeBList;

	string f;
	int sep, typeA, typeB;

	while(file1 >> f >> sep >> typeA >> typeB){
		this->keyList.push_back(f);
		sepList.push_back(sep);
		typeAList.push_back(typeA);
		typeBList.push_back(typeB);
	}
	file1.close();

	for(int i=0;i<keyList.size();i++){
		fileName = path+"edgeInfoLib/" + keyList[i];
		file2.open(fileName.c_str(), ios::in);
		if(!file2.is_open()){
			cout << "fail to open cluster file : " << fileName << endl;
			exit(1);
		}
		EdgeInformation* ei = new EdgeInformation(sepList[i], typeAList[i], typeBList[i], pairLib);

		int totalClusterNum;
		if(sepList[i] == 1)
			totalClusterNum = pairLib->nbBasePairNum[typeAList[i]*4+typeBList[i]];
		else 
			totalClusterNum = pairLib->nnbBasePairNum[typeAList[i]*4+typeBList[i]];

		double pList[totalClusterNum];
		for(int j=0;j<totalClusterNum;j++){
			pList[j] = 0.0;
		}
		int clusterID;
		double p;
		double pContact = 0.0;

		while(file2 >> clusterID >> p){
			if(clusterID < 0)
				pContact = 1-p;
			else
				pList[clusterID] = p;
		}
		ei->updatePCluster(pList, pContact, pairLib);
		this->eiMap[keyList[i]] = ei;
	}

	
}

EdgeInformationLib::~EdgeInformationLib(){
	for(int i=0;i<this->keyList.size();i++){
		delete eiMap[keyList[i]];
	}
}


} /* namespace NSPpredNA */

