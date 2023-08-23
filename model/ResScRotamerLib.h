/*
 * ResScRotamerLib.h
 *
 */

#ifndef MODEL_RESSCROTAMERLIB_H_
#define MODEL_RESSCROTAMERLIB_H_

#include "geometry/CsMove.h"
#include "geometry/TransMatrix.h"
#include "model/ResScRotamer.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include <vector>
#include <fstream>
#include <time.h>


namespace NSPmodel {


class ResScRotamerLib {
public:

	AtomLib* atLib;
	vector<ResScRotamer*> rotList[20];
	vector<int> rotNum;

	//backbone rotamer index: amino acid specific
	double eALA[1000][10]; //0
	double eCYS[1000][50]; //1
	double eASP[1000][400]; //2
	double eGLU[1000][1000]; //3
	double ePHE[1000][1000]; //4
	double eHIS[1000][1000]; //6
	double eILE[1000][400]; //7
	double eLYS[1000][1000]; //8
	double eLEU[1000][400]; //9
	double eMET[1000][1000]; //10
	double eASN[1000][600]; //11
	double ePRO[1000][50]; //12
	double eGLN[1000][1000]; //13
	double eARG[1000][1000]; //14
	double eSER[1000][50]; //15
	double eTHR[1000][50]; //16
	double eVAL[1000][50]; //17
	double eTRP[1000][2000]; //18
	double eTYR[1000][1200]; //19

	/*
	 * clusterNum:
	 * ASP 50
	 * ASN 50
	 * ILE 50
	 * LEU 50
	 * HIS 100
	 * PHE 100
	 * TYR 100
	 * TRP 100
	 * GLU 1000
	 * GLN 1000
	 * LYS 1000
	 * MET 1000
	 * ARG 1000
	 */

	/*
	 * rot library clustering
	 * rotCluster[x][y][0] is the member number in each cluster
	 */
	int rotCluster[20][1000][210];  //coverage*4
	int rotClusterUnique[20][1000][60]; //coverage*1


	/*
	 * rotamer index to cluster-x1 index
	 */
	int GLUIndex[10000];
	int GLNIndex[10000];
	int LYSIndex[10000];
	int METIndex[10000];
	int ARGIndex[20000];

	/*
	 * rotamer index to cluster-x1 center
	 */
	int GLUCenter[10000];
	int GLNCenter[10000];
	int LYSCenter[10000];
	int METCenter[10000];
	int ARGCenter[20000];

	/*
	 * backbone independent rotamer energy minus center energy
	 */
	double GLUELocal[10000];
	double GLNELocal[10000];
	double LYSELocal[10000];
	double METELocal[10000];
	double ARGELocal[20000];

	ResScRotamerLib();
	ResScRotamerLib(const string& tag);

	int getRotamerID(ResScRotamer* rot);
	int getAminoAcidRotamerClusterNum(int aaType);

	void checkResScRotamerPolarGroups();


	double getEnergy(int bbIndex, ResScRotamer* rot){

		if(bbIndex < 0 || bbIndex > 999) {
			cout << "invalid bbIndex" << endl;
			exit(0);
		}
		int type = rot->aaType;
		if(type == 0)
			return eALA[bbIndex][rot->rotID];
		else if(type == 1)
			return eCYS[bbIndex][rot->rotID];
		else if(type == 2)
			return eASP[bbIndex][rot->rotID];
		else if(type == 3)
			return eGLU[bbIndex][GLUIndex[rot->rotID]] + GLUELocal[rot->rotID];
		else if(type == 4)
			return ePHE[bbIndex][rot->rotID];
		else if(type == 6)
			return eHIS[bbIndex][rot->rotID];
		else if(type == 7)
			return eILE[bbIndex][rot->rotID];
		else if(type == 8)
			return eLYS[bbIndex][LYSIndex[rot->rotID]] + LYSELocal[rot->rotID];
		else if(type == 9)
			return eLEU[bbIndex][rot->rotID];
		else if(type == 10)
			return eMET[bbIndex][METIndex[rot->rotID]] + METELocal[rot->rotID];
		else if(type == 11)
			return eASN[bbIndex][rot->rotID];
		else if(type == 12)
			return ePRO[bbIndex][rot->rotID];
		else if(type == 13)
			return eGLN[bbIndex][GLNIndex[rot->rotID]] + GLNELocal[rot->rotID];
		else if(type == 14)
			return eARG[bbIndex][ARGIndex[rot->rotID]] + ARGELocal[rot->rotID];
		else if(type == 15)
			return eSER[bbIndex][rot->rotID];
		else if(type == 16)
			return eTHR[bbIndex][rot->rotID];
		else if(type == 17)
			return eVAL[bbIndex][rot->rotID];
		else if(type == 18)
			return eTRP[bbIndex][rot->rotID];
		else if(type == 19)
			return eTYR[bbIndex][rot->rotID];
		return 0;
	}

	void checkRotamerUniqueID() {
		for(int i=0;i<20;i++){
			int rotNum = this->rotList[i].size();
			for(int j=0;j<rotNum;j++){
				ResScRotamer* rot = rotList[i][j];
				int atomNum = rot->atomNum;
				for(int k=0;k<atomNum;k++){
					int uniqueID = rot->uniqueIDs[k];
					if(uniqueID > 166 || uniqueID < 0){
						cout << "invalid uniqueID " << uniqueID << endl;
						cout << rot->aaType << endl;
						cout << rot->rotID << endl;
						exit(0);
					}
				}
			}
		}
	}

	virtual ~ResScRotamerLib();
};

}
#endif /* MODEL_RESSCROTAMERLIB_H_ */
