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

	this->weight = 999.9;
	
	if(sep == 1 || sep == -1) {
		this->pContact = 1.0;
		this->weight = 0.0;
	}

	else if(sep == 0) {
		this->pContact = 1.0;
		this->weight = 0.0;
	}
	else
		this->pContact = 0.0;


	

	for(int i=0;i<totalClusterNum;i++){
		double e = pairLib->getEnergyWithOxy(i, typeA, typeB, sep);
		this->weight += pCluster[i]*e;
	}

	this->fixed = false;
	this->ssSepKey = "UNK";

	this->pairLibType = pairLib->libType;

}


EdgeInformation::~EdgeInformation() {
	delete [] pCluster;
}

void EdgeInformation::updatePCluster(double* pList, double pContact, BasePairLib* pairLib){

	if(this->pairLibType != pairLib->libType){
		this->pairLibType = pairLib->libType;
		delete [] this->pCluster;
		if(sep == 1)
			this->totalClusterNum = pairLib->nbBasePairNum[this->typeA*4+this->typeB];
		else if(sep == -1)
			this->totalClusterNum = pairLib->nbBasePairNum[this->typeB*4+this->typeA];
		else if(sep == 0)
			this->totalClusterNum = 1;
		else
			this->totalClusterNum = pairLib->nnbBasePairNum[this->typeA*4+this->typeB];
		this->pCluster = new double[totalClusterNum];
	}

	


	double sum = 0.0;
	for(int i=0;i<totalClusterNum;i++){
			sum += pList[i];
	}
	if(sum == 0.0) {
		
		this->weight = 99.9;
		this->validClusterNum = totalClusterNum;
		for(int i=0;i<totalClusterNum;i++){
			this->pCluster[i] = 1.0/totalClusterNum;
		}
		this->pContact = 0;
		return;
	}
	
	this->weight = 0;
	this->validClusterNum = 0;
	for(int i=0;i<totalClusterNum;i++){
		this->pCluster[i] = pList[i]/sum;
		if(pList[i] > 0) 
			this->validClusterNum++;
		double e = pairLib->getEnergyWithOxy(i, typeA, typeB, sep);
		this->weight += pCluster[i]*e;
	}

	this->pContact = pContact;
}

void EdgeInformation::setUniqueCluster(int clusterID, BasePairLib* pairLib){

	if(this->pairLibType != pairLib->libType){
		this->pairLibType = pairLib->libType;
		delete [] this->pCluster;
		if(sep == 1)
			this->totalClusterNum = pairLib->nbBasePairNum[this->typeA*4+this->typeB];
		else if(sep == -1)
			this->totalClusterNum = pairLib->nbBasePairNum[this->typeB*4+this->typeA];
		else if(sep == 0)
			this->totalClusterNum = 1;
		else
			this->totalClusterNum = pairLib->nnbBasePairNum[this->typeA*4+this->typeB];
		this->pCluster = new double[totalClusterNum];
	}

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
		this->weight = pairLib->getEnergyWithOxy(clusterID, typeA, typeB, sep);
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

void EdgeInformation::setToLibReversePCluster(const string& ssSepType, EdgeInformationLib* eiLib){
	if(eiLib->reverseType.find(ssSepType) == eiLib->reverseType.end()){
		cout << "invalid ssSepType: " << this->ssSepKey << endl;
		exit(0);
	}

	this->ssSepKey = eiLib->reverseType[ssSepType];

	map<string, EdgeInformation*>::iterator it;
	it = eiLib->eiMap.find(this->ssSepKey);
	if(it != eiLib->eiMap.end()){
		this->copyClusterFrom(it->second);
	}
	else {
		cout << "invalid ssSepType: " << this->ssSepKey << endl;
		exit(0);
	}
}

void EdgeInformation::setClusterList(vector<int>& clusterList, vector<double>& pList, BasePairLib* pairLib){
	
	if(this->pairLibType != pairLib->libType){
		this->pairLibType = pairLib->libType;
		delete [] this->pCluster;
		if(sep == 1)
			this->totalClusterNum = pairLib->nbBasePairNum[this->typeA*4+this->typeB];
		else if(sep == -1)
			this->totalClusterNum = pairLib->nbBasePairNum[this->typeB*4+this->typeA];
		else if(sep == 0)
			this->totalClusterNum = 1;
		else
			this->totalClusterNum = pairLib->nnbBasePairNum[this->typeA*4+this->typeB];
		this->pCluster = new double[totalClusterNum];
	}

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
			e = pairLib->getEnergyWithOxy(clusterList[i], typeA, typeB, sep);
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

EdgeInformationLib::EdgeInformationLib(){
	string path = NSPdataio::datapath();
	BasePairLib* pairLib = new BasePairLib("stat");
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
	string revType;

	while(file1 >> f >> sep >> typeA >> typeB >> revType){
		this->keyList.push_back(f);
		sepList.push_back(sep);
		typeAList.push_back(typeA);
		typeBList.push_back(typeB);
		this->reverseType[f] = revType;
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
		else if(sepList[i] == -1)
			totalClusterNum = pairLib->nbBasePairNum[typeBList[i]*4+typeAList[i]];
		else 
			totalClusterNum = pairLib->nnbBasePairNum[typeAList[i]*4+typeBList[i]];

		
		double pList[totalClusterNum];
		for(int j=0;j<totalClusterNum;j++){
			pList[j] = 0.0;
		}
		int clusterID;
		double p;
		double pContact = 1.0;

		int n = 0;
		while(file2 >> clusterID >> p){
			n++;
			if(clusterID < 0)
				pContact = pContact-p;
			else if(p > 0)
				pList[clusterID] = p;
		}

		if(sepList[i] == 1 ){
			int contactPairNum = pairLib->nbContactBasePairNum[typeAList[i]*4+typeBList[i]];
			pContact = 0.0;
			for(int j=0;j<contactPairNum;j++){
				pContact += pList[j];
			}
		}
		else if(sepList[i] == -1){
			int contactPairNum = pairLib->nbContactBasePairNum[typeBList[i]*4+typeAList[i]];
			pContact = 0.0;
			for(int j=0;j<contactPairNum;j++){
				pContact += pList[j];
			}
		}

		ei->updatePCluster(pList, pContact, pairLib);
		this->eiMap[keyList[i]] = ei;

		file2.close();
	}

	delete pairLib;

	
}

EdgeInformationLib::~EdgeInformationLib(){
	for(int i=0;i<this->keyList.size();i++){
		delete eiMap[keyList[i]];
	}

}


} /* namespace NSPpredNA */

