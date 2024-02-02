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
	else
		this->totalClusterNum = pairLib->nnbBasePairNum[this->typeA*4+this->typeB];


	this->pCluster = new double[totalClusterNum];

	for(int i=0;i<totalClusterNum;i++){
		pCluster[i] = 1.0/totalClusterNum;
	}
	this->validClusterNum = totalClusterNum;

	if(sep == 1 || sep == -1)
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

	this->moveType = "all";
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
	this->moveType = "multiple";
}

void EdgeInformation::setUniqueCluster(int clusterID, BasePairLib* pairLib){

	this->moveType = "single";
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
		this->moveType = "all";	
	}
	else {
		this->validClusterNum = clusterList.size();
		
		if(clusterList.size() == 1)
			this->moveType = "single";
		else 
			this->moveType = "multiple";

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


} /* namespace NSPpredNA */

