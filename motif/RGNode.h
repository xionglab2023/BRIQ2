/*
 * RGNode.h
 *
 *  Created on: 2020Äê11ÔÂ24ÈÕ
 *      Author: pengx
 */

#ifndef MOTIF_RGNODE_H_
#define MOTIF_RGNODE_H_

#include <iostream>
#include <vector>
#include <map>
#include <set>

class RGEdge;

namespace NSPmotif {
using namespace std;

class RGNode {
public:

	int seqID;
	int baseType;
	set<int> egSet;


	RGNode(){
		this->seqID = 0;
		this->baseType = 0;
	}

	RGNode(int seqID, char baseType){
		this->seqID = seqID;
		if(baseType == 'A' || baseType == 'a')
			this->baseType = 0;
		else if(baseType == 'U'|| baseType == 'u')
			this->baseType = 1;
		else if(baseType == 'G'|| baseType == 'g')
			this->baseType = 2;
		else if(baseType == 'C'|| baseType == 'c')
			this->baseType = 3;
		else {
			cout << "invalid base type: "<< baseType << endl;
			exit(1);
		}

	}

	virtual ~RGNode();
};

} /* namespace NSPmotif */

#endif /* MOTIF_RGNODE_H_ */
