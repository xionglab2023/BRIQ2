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


	double ene;

	RGEdge(int i){
		this->indexA = i;
		this->indexB = i+1;
		this->isSequentialNeighbor = true;
		this->isWatsonCrickPair = false;
		this->ene = 0.0;
		this->sep = 1;
	}

	RGEdge(int i, int j, double ene, int sep){
		if(i < j){
			this->indexA = i;
			this->indexB = j;
			this->ene = ene;
			this->sep = sep;
		}
		else {
			this->indexA = j;
			this->indexB = i;
			this->ene = ene;
			this->sep = sep;
		}
		this->isSequentialNeighbor = false;
		this->isWatsonCrickPair = false;
		if(sep == 1)
			this->isSequentialNeighbor = true;
	}


	double getEnergy(){
		if(isSequentialNeighbor)
			return this->ene - 0.3;
		else
			return this->ene;
	}

	virtual ~RGEdge();
};

} /* namespace NSPmotif */

#endif /* MOTIF_RGEDGE_H_ */
