/*
 * RotSequence.h
 *
 *  Created on: 2022Äê7ÔÂ17ÈÕ
 *      Author: pengx
 */

#ifndef PREDPROT_ROTSEQUENCE_H_
#define PREDPROT_ROTSEQUENCE_H_
#include <string>
#include <vector>
#include "stdio.h"
#include <iostream>
#include "time.h"
#include "protein/ResMutator.h"

namespace NSPprotein {

using namespace std;

class RotSequence {
public:
	int seqLen;
	int* rotNum;
	int* rotChoice;
	ResName rn;
	vector<ResMutator*> mutList;
	RotSequence(vector<ResMutator*>& resMutList);
	void copyValueFrom(RotSequence* other);
	void setRandomChoice();
	void applyMutation(int pos, int choice);
	string toAASequence();
	int getChoice(int pos);
	int getAAChoice(int pos);
	int choiceNum(int pos);
	virtual ~RotSequence();
};

} /* namespace NSPprotein */

#endif /* PREDPROT_ROTSEQUENCE_H_ */
