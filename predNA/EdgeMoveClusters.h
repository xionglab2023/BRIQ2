/*
 * EdgeMoveClusters.h
 *
 *  Created on: 2023��11��30��
 *      Author: nuc
 */

#ifndef PREDNA_EDGEMOVECLUSTERS_H_
#define PREDNA_EDGEMOVECLUSTERS_H_

#include "model/BasePairLib.h"

namespace NSPpredNA {

class EdgeMoveClusters;
class EdgeMoveClustersLib;

using namespace NSPmodel;

class EdgeMoveClusters{
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
	
	string pairLibType; //xtb or stat
	string ssSepKey; //For example, nb-HH, rnb-HL, nnb-1-2, nnb-3-3

	EdgeMoveClusters(int sep, int typeA, int typeB, BasePairLib* pairLib);

	void copyClusterFrom(EdgeMoveClusters* other){

		if(this->totalClusterNum != other->totalClusterNum){
			delete [] pCluster;
			this->pCluster = new double[other->totalClusterNum];
		}

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
		if(fixed)
			return -999.9;

		if(pContact == 1.0) 
			return weight - rand()*2.0/RAND_MAX;
		else if(rand()*1.0/RAND_MAX < 0.5*pContact)
			return weight - rand()*2.0/RAND_MAX;
		else if(sep < 2)
			return 0.0;
		else
			return 99.9;
	}

	void setToLibPCluster(const string& ssSepType, EdgeMoveClustersLib* eiLib);
	void setToLibReversePCluster(const string& ssSepType, EdgeMoveClustersLib* eiLib);
	void setUniqueCluster(int clusterID, BasePairLib* pairLib);
	void setClusterList(vector<int>& clusterList, vector<double>& pClusters, BasePairLib* pairLib);
	void setFixed(CsMove& cm){
		this->pContact = 0.0;
		this->validClusterNum = 0;
		this->fixed = true;
		this->cm = cm;
		this->weight = -999.9;
	}
	void print(){
		printf("edge information: sep: %d  type: %d %d\n", sep, typeA, typeB);
		for(int i=0;i<totalClusterNum;i++){
			double p = pCluster[i];
			if(p > 0) {
				printf("cluster %3d %7.5f\n", i, p);
			}
		}
	}
	virtual ~EdgeMoveClusters();
};

class EdgeMoveClustersLib{
public:
	map<string, EdgeMoveClusters*> eiMap;
	map<string, string> reverseType;

	vector<string> keyList;
	EdgeMoveClustersLib();
	virtual ~EdgeMoveClustersLib();

};

} /* namespace NSPpredNA */

#endif /* PREDNA_EDGEMOVECLUSTERS_H_ */
