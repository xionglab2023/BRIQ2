/*
 * BasePair6DEnergyTable.h
 *
 */

#ifndef FORCEFIELD_BASEPAIR6DENERGYTABLE_H_
#define FORCEFIELD_BASEPAIR6DENERGYTABLE_H_


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

class CsMoveTo6DKey {

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

	int biIndex[900000]; //interpolation: 500*450*4
	double biWt[900000]; //interpolation: 500*450*4

public:

	CsMoveTo6DKey();
	void getIndexAndWeight(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[5], double outWeights[4]);
	void getIndexAndWeight6DInterpolation(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[13], double outWeights[13]);
	pair<int,int> toIndexPair(const LocalFrame& csA, const LocalFrame& csB, double len);

	/*
	string toKey(const LocalFrame& csA, const LocalFrame& csB, double len);
	string toKey2(const LocalFrame& csA, const LocalFrame& csB, double len);
	int toIndex(const LocalFrame& csA, const LocalFrame& csB, double len);
	int toIndex2(const LocalFrame& csA, const LocalFrame& csB, double len);
	*/
	virtual ~CsMoveTo6DKey();
};

class BasePair6DEnergyTable {
public:

	map<int, double> nbKeysEnergy[36000]; //16*2250, distacne 50 bins, dihedral 45 bins, sphere 2000*2000
	map<int, double> nnbKeysEnergy[36000];

	CsMoveTo6DKey cm2Key;
	map<int,double>::iterator it;


	double wtNb, wtNnb;

	BasePair6DEnergyTable(ForceFieldPara* para);

	double getEnergyBiInterpolation(const LocalFrame csA, const LocalFrame csB, int typeA, int typeB, int sep, double minDistance){
		//bilinear interpolation

		if(minDistance < 1.5) return 0.0;
		if(minDistance >= 5.0) return 0.0;

		double len = csA.origin_.distance(csB.origin_);
		if(len >= 15.0)
			return 0.0;


		int index[5];
		double weights[4];
		double ene[4];

		cm2Key.getIndexAndWeight(csA, csB, len, index, weights);

		for(int i=0;i<4;i++){
			int d2Index = (typeA*4+typeB)*2250+index[i];
			if(d2Index <0 || d2Index >= 36000){
				cout << "invalid d2Index: " << d2Index << endl;
				exit(0);
			}
			if(sep == 1){
				it = nbKeysEnergy[d2Index].find(index[4]);
				if(it != nbKeysEnergy[d2Index].end()){
					ene[i] = wtNb*it->second;
				}
				else
					ene[i] = 0.0;
			}
			else {
				it = nnbKeysEnergy[d2Index].find(index[4]);
				if(it != nnbKeysEnergy[d2Index].end()){
					ene[i] = wtNnb*it->second;
				}
				else
					ene[i] = 0.0;
			}

		}

		double e = weights[0]*ene[0]+weights[1]*ene[1]+weights[2]*ene[2]+weights[3]*ene[3];
		return e;
	}


	double getEnergy(const LocalFrame csA, const LocalFrame csB, int typeA, int typeB, int sep, double minDistance){
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


		for(int i=0;i<4;i++){
			for(int j=0;j<9;j++){
				int d2Index = (typeA*4+typeB)*2250+index[i];
				if(d2Index <0 || d2Index >= 36000){
					cout << "invalid d2Index: " << d2Index << " " << i << " " << index[i] << endl;
					exit(0);
				}
				if(sep == 1){
					it = nbKeysEnergy[d2Index].find(index[4+j]);
					if(it != nbKeysEnergy[d2Index].end()){
						ene[i*9+j] = wtNb*it->second;
					}
					else
						ene[i*9+j] = 0.0;
				}

				else {
					it = nnbKeysEnergy[d2Index].find(index[4+j]);
					if(it != nnbKeysEnergy[d2Index].end()){
						ene[i*9+j] = wtNnb*it->second;
					}
					else
						ene[i*9+j] = 0.0;
				}
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

	double getEnergyNearestNeighbor(const LocalFrame csA, const LocalFrame csB, int typeA, int typeB, int sep, double minDistance){
		if(minDistance >= 5.0) return 0;
		if(minDistance < 1.5) return 0;

		double len = csA.origin_.distance(csB.origin_);
		if(len >= 15.0)
			return 0;

	    pair<int,int> p = cm2Key.toIndexPair(csA, csB, len);
		int pairType = typeA*4+typeB;
		int mapIndex = pairType*2250 + p.first;

		if(sep == 1){
			it = nbKeysEnergy[mapIndex].find(p.second);
			if(it != nbKeysEnergy[mapIndex].end()){
				//printf("indexA: %8d indexB: %8d\n", mapIndex, p.second);
				return wtNb*it->second;
			}
			else
				return 0.0;
		}
		else {
			it = nnbKeysEnergy[mapIndex].find(p.second);
			if(it != nnbKeysEnergy[mapIndex].end()){
				//printf("indexA: %8d indexB: %8d\n", mapIndex, p.second);
				return wtNnb*it->second;
			}
			else
				return 0;
		}

	}
	virtual ~BasePair6DEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_BASEPAIR6DENERGYTABLE_H_ */
