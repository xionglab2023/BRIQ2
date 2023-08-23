/*
 * RNABaseName.h
 *
 */

#ifndef MODEL_RNABASENAME_H_
#define MODEL_RNABASENAME_H_

#include <string>
#include <vector>
#include <map>

namespace NSPmodel {

using namespace std;
class RNABaseName {
private:
	string rnaSeq;
	vector<string> baseTypes;
	map<char,int> charToIntMap;
	map<string, int> stringToIntMap;
public:

	RNABaseName();
	char intToChar(int i);
	int charToInt(char sin);
	string intToString(int i);
	int stringToInt(const string& name);
	char stringToChar(const string& name);
	string charToString(char sin);
	virtual ~RNABaseName();
};

} /* namespace NSPmodel */

#endif /* MODEL_RNABASENAME_H_ */
