/*
 * RotSequence.cpp
 *
 *  Created on: 2022Äê7ÔÂ17ÈÕ
 *      Author: pengx
 */

#include "protein/RotSequence.h"

namespace NSPprotein {

RotSequence::RotSequence(vector<ResMutator*>& resMutList) {
	// TODO Auto-generated constructor stub
	this->seqLen = resMutList.size();
	if(seqLen == 0) {
		cout << "empty resMutList" << endl;
		exit(1);
	}

	this->rotChoice = new int[seqLen];
	this->rotNum = new int[seqLen];
	for(int i=0;i<seqLen;i++){
		mutList.push_back(resMutList[i]);
		this->rotNum[i] = mutList[i]->rotList.size();
		this->rotChoice[i] = 0;
	}

}

void RotSequence::copyValueFrom(RotSequence* other){
	if(this->seqLen != other->seqLen){
		cerr << "rotamer sequence length not equal" << endl;
		exit(1);
	}

	this->mutList.clear();
	for(int i=0;i<seqLen;i++){
		this->rotChoice[i] = other->rotChoice[i];
		this->mutList.push_back(other->mutList[i]);
	}
}

void RotSequence::setRandomChoice(){
	for(int i=0;i<this->seqLen;i++){
		int randNum = rand()%rotNum[i];
		this->rotChoice[i] = randNum;
	}
}

void RotSequence::applyMutation(int pos, int choice){
	this->rotChoice[pos] = choice;
}

string RotSequence::toAASequence(){

	char seq[this->seqLen+1];
	for(int i=0;i<seqLen;i++){
		seq[i] = rn.intToSin(this->mutList[i]->rotList[rotChoice[i]]->aaType);
	}
	seq[this->seqLen] = '\0';
	string s = seq;
	return s;
}

int RotSequence::getChoice(int pos){
	return rotChoice[pos];
}

int RotSequence::getAAChoice(int pos){
	return mutList[pos]->rotList[rotChoice[pos]]->aaType;
}

int RotSequence::choiceNum(int pos){
	return rotNum[pos];
}

RotSequence::~RotSequence() {
	// TODO Auto-generated destructor stub
	delete [] rotNum;
	delete [] rotChoice;
}

} /* namespace NSPprotein */
