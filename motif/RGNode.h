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

		if(baseType == 'A')
			this->baseType = 0;
		else if(baseType == 'U')
			this->baseType = 1;
		else if(baseType == 'G')
			this->baseType = 2;
		else if(baseType == 'C')
			this->baseType = 3;
		else if(baseType == 'a')
			this->baseType = 4;
		else if(baseType == 't')
			this->baseType = 5;
		else if(baseType == 'g')
			this->baseType = 6;
		else if(baseType == 'c')
			this->baseType = 7;
		else {
			this->baseType = -1;
			cout << "invalid base type: "<< baseType << endl;
			exit(1);
		}

	}

	virtual ~RGNode();
};

} /* namespace NSPmotif */

#endif /* MOTIF_RGNODE_H_ */
