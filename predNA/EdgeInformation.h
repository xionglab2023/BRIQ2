/*
 * EdgeInformation.h
 *
 *  Created on: 2023��11��30��
 *      Author: nuc
 */

#ifndef PREDNA_EDGEINFORMATION_H_
#define PREDNA_EDGEINFORMATION_H_

#include "model/BasePairLib.h"

namespace NSPpredNA {

class EdgeInformation;
class EdgeInformationLib;

using namespace NSPmodel;

class EdgeInformation{
public:
	int sep;
	int typeA;
	int typeB;

	int totalClusterNum;
	int validClusterNum;

	double* pCluster;
	double pContact;

	double weight;

	bool fixed;
	CsMove cm;

	string ssSepKey; //For example, nb-HH, rnb-HL, nnb-1-2, nnb-3-3

	EdgeInformation(int sep, int typeA, int typeB, BasePairLib* pairLib);

	void copyClusterFrom(EdgeInformation* other){

		this->totalClusterNum = other->totalClusterNum;
		this->validClusterNum = other->validClusterNum;
		for(int i=0;i<totalClusterNum;i++){
			this->pCluster[i] = other->pCluster[i];
		}
		this->pContact = other->pContact;
		this->weight = other->weight;
	}

	void updatePCluster(double* pList, double pContact, BasePairLib* pairLib);
	
	double getRandomWeight() {
		if(pContact == 1.0) 
			return weight - rand()*1.0/RAND_MAX;
		else if(sep < 2 && rand()*1.0/RAND_MAX < pContact)
			return weight - rand()*1.0/RAND_MAX;
		else if(sep == 2 && rand()*1.0/RAND_MAX < pContact*0.5)
			return weight - rand()*1.0/RAND_MAX;
		else if(sep < 2)
			return 0.0;
		else 
			return 9.9;
	}

	void setToLibPCluster(const string& ssSepType, EdgeInformationLib* eiLib);

	void setUniqueCluster(int clusterID, BasePairLib* pairLib);
	void setClusterList(vector<int>& clusterList, vector<double>& pClusters, BasePairLib* pairLib);
	void setFixed(CsMove& cm){
		this->validClusterNum = 0;
		this->fixed = true;
		this->cm = cm;
		this->weight = -999.9;
	}

	virtual ~EdgeInformation();
};

class EdgeInformationLib{
public:
	map<string, EdgeInformation*> eiMap;
	vector<string> keyList;

	EdgeInformationLib(BasePairLib* pairLib);
	virtual ~EdgeInformationLib();

};

} /* namespace NSPpredNA */

#endif /* PREDNA_EDGEINFORMATION_H_ */
