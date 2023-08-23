/*
 * ResName.h
 *
 */

#ifndef DESIGNSEQ_RESNAME_H_
#define DESIGNSEQ_RESNAME_H_

#include <string>
#include <vector>
#include <map>

namespace NSPmodel {

using namespace std;
class ResName {
private:
	string aaSeq;
	vector<float> vols;
	vector<string> triNameList;
	map<string,int> triToIntMap;
	map<char,int> sinToIntMap;
	map<string, int> rnaNameToInt;
	vector<string> augc;

public:
	ResName();
	char intToSin(int i) const;
	int sinToInt(char sin) const;
	char triToSin(const string& tri) const;
	string sinToTri(char sin) const;
	int triToInt(const string& tri) const;
	string intToTri(int i) const;
	float getVol(int aa) const;
	bool isStandardAminoAcid(const string& tri);
	bool isAminoAcid(const string& tri);
	int chiNum(const string& tri) const;
	bool isRNABase(const string& baseName) const;
	string toStandardBase(const string& baseName) const;
	virtual ~ResName();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_RESNAME_H_ */
