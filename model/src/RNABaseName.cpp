/*
 * RNABaseName.cpp
 *
 */

#include "model/RNABaseName.h"

namespace NSPmodel {

RNABaseName::RNABaseName() {
	// TODO Auto-generated constructor stub
	this->rnaSeq = "AUGCatgcN";
	this->charToIntMap['A'] = 0;
	this->charToIntMap['U'] = 1;
	this->charToIntMap['G'] = 2;
	this->charToIntMap['C'] = 3;
	this->charToIntMap['a'] = 4;
	this->charToIntMap['t'] = 5;
	this->charToIntMap['g'] = 6;
	this->charToIntMap['c'] = 7;
	this->charToIntMap['N'] = 8;
	this->charToIntMap['X'] = 8;

	this->stringToIntMap["A"] = 0;
	this->stringToIntMap["U"] = 1;
	this->stringToIntMap["G"] = 2;
	this->stringToIntMap["C"] = 3;
	this->stringToIntMap["DA"] = 4;
	this->stringToIntMap["DT"] = 5;
	this->stringToIntMap["DG"] = 6;
	this->stringToIntMap["DC"] = 7;

	baseTypes.push_back("A");
	baseTypes.push_back("U");
	baseTypes.push_back("G");
	baseTypes.push_back("C");
	baseTypes.push_back("DA");
	baseTypes.push_back("DT");
	baseTypes.push_back("DG");
	baseTypes.push_back("DC");
	baseTypes.push_back("N");
}

char RNABaseName::intToChar(int i) {
	if(i>=0 && i <8)
		return rnaSeq[i];
	return 'N';
}

int RNABaseName::charToInt(char sin) {
	map<char,int>::const_iterator it = charToIntMap.find(sin);
	if(it != charToIntMap.end())
		return it->second;
	return 8;
}

string RNABaseName::intToString(int i) {
	if(i >=0 && i < 8)
		return baseTypes[i];
	return "N";
}

int RNABaseName::stringToInt(const string& name){
	map<string, int>::const_iterator it = stringToIntMap.find(name);
	if(it != stringToIntMap.end())
		return it->second;
	return 8;
}

char RNABaseName::stringToChar(const string& name){
	map<string, int>::const_iterator it = stringToIntMap.find(name);
	if(it != stringToIntMap.end())
		return rnaSeq[it->second];
	return 'N';
}

string RNABaseName::charToString(char sin){
	map<char,int>::const_iterator it = charToIntMap.find(sin);
	if(it != charToIntMap.end())
		return baseTypes[it->second];
	return "N";
}

RNABaseName::~RNABaseName() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
