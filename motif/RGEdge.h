/*
 * RGEdge.h
 *
 *  Created on: 2020Äê11ÔÂ24ÈÕ
 *      Author: pengx
 */

#ifndef MOTIF_RGEDGE_H_
#define MOTIF_RGEDGE_H_

#include <stdio.h>
#include <iostream>

namespace NSPmotif {

using namespace std;

class RGEdge {
public:

	//indexA < indexB
	int indexA;
	int indexB;

	int sep;

	bool isSequentialNeighbor;
	bool isWatsonCrickPair;


	double eBB;
	double eBO;
	double eBP;

	RGEdge(int i){
		this->indexA = i;
		this->indexB = i+1;
		this->isSequentialNeighbor = true;
		this->isWatsonCrickPair = false;
		this->eBB = 0.0;
		this->eBO = 0.0;
		this->eBP = 0.0;
		this->sep = 1;
	}

	RGEdge(int i, int j, double bb, double bo, double bp, int sep){
		if(i < j){
			this->indexA = i;
			this->indexB = j;
			this->eBB = bb;
			this->eBO = bo;
			this->eBP = bp;
			this->sep = sep;
		}
		else {
			this->indexA = j;
			this->indexB = i;
			this->eBB = bb;
			this->eBO = bo;
			this->eBP = bp;
			this->sep = sep;
		}
		this->isSequentialNeighbor = false;
		this->isWatsonCrickPair = false;
	}


	double getEnergy(){
		if(isSequentialNeighbor)
			return this->eBB + this->eBO + this->eBP - 0.3;
		else
			return this->eBB + this->eBO + this->eBP;
	}

	virtual ~RGEdge();
};

} /* namespace NSPmotif */

#endif /* MOTIF_RGEDGE_H_ */
