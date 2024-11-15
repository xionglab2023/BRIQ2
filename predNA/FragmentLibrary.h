/*
 * FragmentLibrary.h
 *
 */

#ifndef predNA_FRAGMENTLIBRARY_H_
#define predNA_FRAGMENTLIBRARY_H_

#include "geometry/CsMove.h"
#include "geometry/TransMatrix.h"
#include "model/RiboseRotamer.h"
#include "model/RotamerLib.h"
#include "model/RNABaseLib.h"
#include "model/PhosphateRotamer.h"
#include "model/BaseDistanceMatrix.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include <vector>
#include <fstream>
#include <time.h>



namespace NSPpredNA {

using namespace std;
using namespace NSPgeometry;
using namespace NSPdataio;
using namespace NSPmodel;


class RiboConnectLib {
public:
	vector<CsMove> cmList;
	int num;
	RiboConnectLib();
	CsMove getRandomMove();
	virtual ~RiboConnectLib();
};


class F2Fragment {
public:
	CsMove cm;
	bool hasRibose;
	bool isReverse;
	bool baseFlip;


	int typeA;
	int typeB;
	PhosphateRotamer* pl;

	int level1ID;
	int level2ID;
	string mvType;


	F2Fragment(int typeA, int typeB);
	F2Fragment(const string& line, int typeA, int typeB, RotamerLib* rotLib);
	double distanceTo(F2Fragment* other);
	void printPDB(const string& outfile);
	virtual ~F2Fragment();
};

class F3Fragment {
public:
	CsMove cmAB;
	CsMove cmBC;
	RiboseRotamer* rotA;
	RiboseRotamer* rotB;
	RiboseRotamer* rotC;
	int typeA;
	int typeB;
	int typeC;
	PhosphateRotamer* plA;
	PhosphateRotamer* plB;
	F3Fragment(const string& line, int typeA, int typeB, int typeC, RotamerLib* rotLib);
	void printPDB(const string& outfile);
	virtual ~F3Fragment();
};


class F2FragmentLib {
public:
	int rotNum;

	vector<F2Fragment*> fragListLevel0;
	vector<F2Fragment*> fragListLevel1;
	vector<F2Fragment*> fragListLevel2;

	bool hasRibose;
	bool isReverse;

	F2FragmentLib(const string& tag, RotamerLib* rotLib);
	F2FragmentLib(const string& tag, int typeA, int typeB, RotamerLib* rotLib);


	F2Fragment* getRandomFrag(){
		return fragListLevel2[rand()%rotNum*rotNum + rand()%rotNum];
	}

	F2Fragment* getRandomFragLevel2(int lv1Index){
		return fragListLevel2[lv1Index*rotNum + rand()%rotNum];
	}

	F2Fragment* getFrag(int idA, int idB){
		return fragListLevel2[idA*rotNum+idB];
	}

	virtual ~F2FragmentLib();
};

class F3FragmentLib {
public:
	int rotNum;
	vector<F3Fragment*> f3List;
	F3FragmentLib(const string& tag, int typeA, int typeB, int typeC, RotamerLib* rotLib);
	F3Fragment* getRandomFrag(){
		return f3List[rand()%rotNum];
	}
	virtual ~F3FragmentLib();
};

class FragmentLibrary {
public:

	vector<F2FragmentLib*> wcNb;
	vector<F2FragmentLib*> nwcNb;
	vector<F2FragmentLib*> loopNb;
	vector<F2FragmentLib*> revNb;
	vector<F2FragmentLib*> revWcNb;
	vector<F2FragmentLib*> wcPair;
	vector<F2FragmentLib*> nwcPair;
	vector<F2FragmentLib*> bulge13;
	//vector<F2FragmentLib*> bulge14;
	//vector<F2FragmentLib*> revBulge13;
	//vector<F2FragmentLib*> revBulge14;

	F2FragmentLib* agLib;
	F2FragmentLib* gaLib;

	vector<F3FragmentLib*> acF3Lib;
	vector<F3FragmentLib*> allF3Lib;
	vector<F2FragmentLib*> jumpLib;

	RiboConnectLib* ribLib;

	FragmentLibrary(RotamerLib* rotLib);
	virtual ~FragmentLibrary();
};

} /* namespace NSPforcefield */

#endif /* predNA_FRAGMENTLIBRARY_H_ */
