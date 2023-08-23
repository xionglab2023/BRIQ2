/*
 * ResMutator.h
 *
 *  Created on: 2022Äê7ÔÂ16ÈÕ
 *      Author: pengx
 */

#ifndef PREDPROT_RESMUTATOR_H_
#define PREDPROT_RESMUTATOR_H_

#include "model/ResBBRotamer.h"
#include "model/ResBBRotamerLib.h"
#include "model/ResConformer.h"
#include "model/ResScRotamerLib.h"
#include <float.h>
#include <stdlib.h>
#include "model/StructureModel.h"
#include "model/ResName.h"

namespace NSPprotein {

using namespace NSPmodel;
using namespace NSPtools;

class ResMutator {

public:
	vector<ResScRotamer*> rotList;
	vector<int> aaList;
	vector<vector<int>> aaRotIDList;
	int aaChoiceNum;
	int choiceNum;

	ResMutator(){
		this->aaChoiceNum = 0;
		this->choiceNum = 0;
	}

	ResMutator(vector<ResScRotamer*>& scRotList){
		bool hasAA[20];
		for(int i=0;i<20;i++){
			hasAA[i] = false;
		}
		for(int i=0;i<scRotList.size();i++){
			rotList.push_back(scRotList[i]);
			hasAA[scRotList[i]->aaType] = true;
		}
		int duplicate = 1;
		for(int i=0;i<20;i++){
			aaRotIDList.push_back(vector<int>());
			if(hasAA[i])
			{
				if(i == 0 || i == 5 || i==12)
					duplicate = 1;
				else if(i == 1 || i == 15 || i == 16 || i == 17)
					duplicate = 2;
				else if(i == 2 || i == 7 || i == 9 || i == 11)
					duplicate = 3;
				else
					duplicate = 4;
				for(int j=0;j<duplicate;j++){
					aaList.push_back(i);
				}
			}
		}
		this->aaChoiceNum = aaList.size();
		this->choiceNum = scRotList.size();

		if(this->aaChoiceNum == 0){
			cout << "aa choice num is zero!" << endl;
			exit(1);
		}
		for(int i=0;i<scRotList.size();i++){
			aaRotIDList[scRotList[i]->aaType].push_back(i);
		}
	}

	ResScRotamer* getRandomRotamer(){
		int randAA = aaList[rand()%aaChoiceNum];
		int aaRotNum = aaRotIDList[randAA].size();
		return rotList[aaRotIDList[randAA][rand()%aaRotNum]];
	}

	virtual ~ResMutator();

}; /* namespace NSPredprot */

}
#endif /* PREDPROT_RESMUTATOR_H_ */
