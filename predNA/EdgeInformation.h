/*
 * EdgeInformation.h
 *
 *  Created on: 2023Äê11ÔÂ30ÈÕ
 *      Author: nuc
 */

#ifndef PREDNA_EDGEINFORMATION_H_
#define PREDNA_EDGEINFORMATION_H_

#include "model/BasePairLib.h"

namespace NSPpredNA {

using namespace NSPmodel;

class EdgeInformation{
public:
	int sep;
	int typeA;
	int typeB;

	int totalClusterNum;

	double* pCluster;
	double pContact;

	double weight;

	bool fixed;
	CsMove cm;

	EdgeInformation(int sep, int typeA, int typeB, BasePairLib* pairLib);
	void updatePCluster(double* pList, double pContact, BasePairLib* pairLib);
	void setUniqueCluster(int clusterID, BasePairLib* pairLib);
	void setFixed(CsMove& cm){
		this->fixed = true;
		this->cm = cm;
		this->weight = -999.9;
	}

	virtual ~EdgeInformation();
};

} /* namespace NSPpredNA */

#endif /* PREDNA_EDGEINFORMATION_H_ */
