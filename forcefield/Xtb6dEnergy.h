/*
 * Xtb6dEnergy.h
 */

#ifndef FORCEFIELD_XTB6DENERGY_H_
#define FORCEFIELD_XTB6DENERGY_H_

#include "dataio/datapaths.h"
#include "forcefield/ForceFieldPara.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include <time.h>
#include <math.h>
#include <map>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;

class CsMoveTo6DKeySP500 {

private:

	map<int, int> sphereKeyMap;
	map<int, int>::iterator it;

	map<int, int> keyIdToListID; //key id range: 0~400*400*400

	int spIndex1[753848];
	int spIndex2[753848];
	int spIndex3[753848];

	double spWt1[753848];
	double spWt2[753848];
	double spWt3[753848];
	int spherePointNum;

	int biIndex[800000]; //interpolation: 500*400*4
	double biWt[800000]; //interpolation: 500*400*4

public:

	CsMoveTo6DKeySP500();
	void getIndexAndWeight(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[5], double outWeights[4]);
	void getIndexAndWeight6DInterpolation(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[13], double outWeights[13]);
	int toIndex(const LocalFrame& csA, const LocalFrame& csB, double len);

	virtual ~CsMoveTo6DKeySP500();

};

class Xtb6dEnergy {
public:

	map<int, double> nnbKeysEnergy[16]; // distacne 50 bins, dihedral 40 bins, sphere 500*500

	CsMoveTo6DKeySP500 cm2Key;
	map<int,double>::iterator it;


	Xtb6dEnergy(ForceFieldPara* para);

	double getEnergy(const LocalFrame csA, const LocalFrame csB, int typeA, int typeB, double minDistance){
		//6D interpolation

		if(minDistance < 1.5) return 0.0;
		if(minDistance >= 5.0) return 0.0;

		double len = csA.origin_.distance(csB.origin_);
		if(len >= 15.0)
			return 0.0;

		int index[13];
		double weights[13];
		double ene[36];

		cm2Key.getIndexAndWeight6DInterpolation(csA, csB, len, index, weights);


		for(int i=0;i<4;i++){ //len-ang index
			for(int j=0;j<9;j++){ //spA-spB index
                int pairIndex = typeA*4+typeB;
                int fullIndex = index[i]*250000+index[4+j];

				it = nnbKeysEnergy[pairIndex].find(fullIndex);
				if(it != nnbKeysEnergy[pairIndex].end()){
					ene[i*9+j] = it->second;
				}
				else
					ene[i*9+j] = 0.0;
			}
		}

		double e = 0;
		for(int i=0;i<4;i++){
			for(int j=0;j<9;j++){
				e += weights[i] * weights[4+j] * ene[i*9+j];
				//printf("%-2d %6.3f %5.4f\n", i, ene[i*9+j], weights[i]*weights[4+j]);
			}
		}
		return e;
	}

	virtual ~Xtb6dEnergy();
};

}

#endif /* FORCEFIELD_XTB6DENERGY_H_*/