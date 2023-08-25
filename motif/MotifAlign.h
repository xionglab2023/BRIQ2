/*
 * MotifAlign.h
 *
 *  Created on: 2020Äê12ÔÂ3ÈÕ
 *      Author: pengx
 */

#ifndef MOTIF_MOTIFALIGN_H_
#define MOTIF_MOTIFALIGN_H_

#include "model/StructureModel.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "geometry/RMSD.h"

namespace NSPmotif {

using namespace std;
using namespace NSPmodel;
using namespace NSPgeometry;

typedef struct idDM{
	int index;
	double dd;
	inline bool operator < (const idDM &x) const {
		return dd < x.dd ;
    }
}idDM;


class MotifAlign {
public:

	int nA;
	int nB;
	idDM* dm;
	int* align;
	double* alignDD;

	double d0;
	double xLN;


	double bestRM;
	double bestRMS;
	int* bestAlign;

	vector<LocalFrame> csAList;
	vector<LocalFrame> csBList;
	vector<LocalFrame> csBListTrans;

	XYZ* initListA;
	XYZ* initListB;
	XYZ* transedB;
	bool* alignedA;
	bool* alignedB;


	MotifAlign(RNAPDB* pdbA, RNAPDB* pdbB){
		AtomLib* atLib = new AtomLib();
		vector<RNABase*> baseListA = pdbA->getValidBaseList(atLib);
		vector<RNABase*> baseListB = pdbB->getValidBaseList(atLib);

		if(baseListA.size() < baseListB.size()){
			baseListA = pdbB->getValidBaseList(atLib);
			baseListB = pdbA->getValidBaseList(atLib);
		}

		this->nA = baseListA.size();
		this->nB = baseListB.size();

		double LN = 0.5*(nA+nB);
		this->xLN = 1.0/LN;

		this->d0 = (nA+nB)/22.3 - 0.28;
		if(d0 < 0.05)
			d0 = 0.05;

		/*
		cout << "nA: " << nA << endl;
		cout << "nB: " << nB << endl;
		cout << "LN: " << LN << endl;
		cout << "d0: " << d0 << endl;
		*/

		this->dm = new idDM[nA*nB];
		this->align = new int[nA];
		this->alignDD = new double[nA];
		this->initListA = new XYZ[nA*3];
		this->initListB = new XYZ[nB*3];
		this->transedB = new XYZ[nB*3];
		this->alignedA = new bool[nA];
		this->alignedB = new bool[nB];
		this->bestAlign = new int[nA];

		for(int i=0;i<nA*nB;i++){
			this->dm[i]= {i, 0.0};
		}

		XYZ xa(0 , 0 , 0);
		XYZ xb(3.464 , 2 , 0);
		XYZ xc(3.464 , -2 , 0);

		this->bestRM = 0.0;
		this->bestRMS = 99.9;

		for(int i=0;i<baseListA.size();i++){
			LocalFrame csA = baseListA[i]->getCoordSystem();
			initListA[i*3] = local2global(csA, xa);
			initListA[i*3+1] = local2global(csA, xb);
			initListA[i*3+2] = local2global(csA, xc);
			csAList.push_back(csA);
		}

		for(int i=0;i<baseListB.size();i++){
			LocalFrame csB = baseListB[i]->getCoordSystem();
			initListB[i*3] = local2global(csB, xa);
			initListB[i*3+1] = local2global(csB, xb);
			initListB[i*3+2] = local2global(csB, xc);
			csBList.push_back(csB);
			csBListTrans.push_back(csB);
		}
	}

	void alignAtPosition(int idA, int idB){
		LocalFrame csA = csAList[idA];
		LocalFrame csB = csBList[idB];

		for(int i=0;i<nB*3;i++){
			XYZ localT = global2local(csB, this->initListB[i]);
			this->transedB[i] = local2global(csA, localT);
		}
	}

	void updateCsMove(CsMove& cm){
		LocalFrame cs;
		cs = cs + cm;
		XYZ xa(0 , 0 , 0);
		XYZ xb(1.732 , 1 , 0);
		XYZ xc(1.732 , -1 , 0);
		for(int i=0;i<nB;i++){
			LocalFrame csTrans= csBList[i] + cm;
			this->transedB[i*3] = local2global(csTrans, xa);
			this->transedB[i*3+1] = local2global(csTrans, xb);
			this->transedB[i*3+2] = local2global(csTrans, xc);
		}
	}

	double calculateRMScore(){
		double dd;
		for(int i=0;i<nA;i++){
			for(int j=0;j<nB;j++){
				dd = squareDistance(initListA[i*3], transedB[j*3]) + squareDistance(initListA[i*3+1], transedB[j*3+1]) + squareDistance(initListA[i*3+2], transedB[j*3+2]);
				this->dm[i*nB+j] = {i*nB+j, dd};
			}
		}

		/*
		cout << "coordA" << endl;
		for(int i=0;i<nA*3;i++){
			cout << this->initListA[i].toString() << endl;
		}

		cout << "coordB" << endl;
		for(int i=0;i<nB*3;i++){
			cout << this->initListB[i].toString() << endl;
		}

		cout << "transB" << endl;
		for(int i=0;i<nB*3;i++){
			cout << this->transedB[i].toString() << endl;
		}
		*/

		vector<int> listA;
		vector<int> listB;
		vector<double> listC;

		for(int i=0;i<nA*nB;i++){
			listA.push_back(i/nB);
			listB.push_back(i%nB);
			listC.push_back(dm[i].dd);
		}


		sort(dm, dm+nA*nB);

		vector<int> listD;
		vector<int> listE;
		vector<double> listF;

		for(int i=0;i<nA*nB;i++){

			listD.push_back(dm[i].index/nB);
			listE.push_back(dm[i].index%nB);
			listF.push_back(dm[i].dd);
		}

		for(int i=0;i<nA*nB;i++){
			int a = i/nB;
			int b = i%nB;
			//printf("%-2d %-2d %6.4f %-2d %-2d %6.4f\n", listA[i], listB[i], listC[i], listD[i], listE[i], listF[i] );
		}

		for(int i=0;i<nA;i++){
			align[i] = -1;
			alignDD[i] = 99999.9;
			alignedA[i] = false;
		}

		for(int i=0;i<nB;i++){
			alignedB[i] = false;
		}

		int N = nA*nB;
		int idA, idB;
		int alignedBase = 0;

		for(int i=0;i<N;i++){
			idB = this->dm[i].index%nB;
			idA = this->dm[i].index/nB;
			dd = this->dm[i].dd;
		//	cout << "idA: " << idA << " idB: " << idB << " " << dd << endl;
			if(dd > 30.0) break;
			if(!alignedA[idA] && !alignedB[idB]){
				//cout << "assign align: " << idA << " " << idB << endl;
				align[idA] = idB;
				alignDD[idA] =dd;
				alignedA[idA] = true;
				alignedB[idB] = true;
				alignedBase++;
			}
		}

	//	cout << "aligned base number: " << alignedBase << endl;

		double RG = 0;
		double xD0D0 = 1.0/d0/d0;

		double rmsd = 0;

		for(int i=0;i<nA;i++){
			if(alignedA[i]){
				rmsd += alignDD[i];
				RG += xLN*(1.0/(1+alignDD[i]*0.3333333*xD0D0));
			}
		}

		cout << "RG: " << RG << endl;

		rmsd = sqrt(rmsd/alignedBase);

		if(RG > bestRM){
			bestRM = RG;
			bestRMS = rmsd;
			for(int i=0;i<nA;i++){

				bestAlign[i] = align[i];
			}
		}

	//	cout << "RG score: " << RG << endl;
		return RG;
	}

	void getTransformFromAlign(){
		vector<XYZ> points1;
		vector<XYZ> points2;

		for(int i=0;i<nA;i++){
			if(align[i] >= 0) {
				int j = align[i];
				cout << i << " " << j << endl;
				points1.push_back(initListA[i*3]);
				points1.push_back(initListA[i*3+1]);
				points1.push_back(initListA[i*3+2]);

				points2.push_back(initListB[j*3]);
				points2.push_back(initListB[j*3+1]);
				points2.push_back(initListB[j*3+2]);

			}
		}

		XYZ Acog = getCOG(points1);
		XYZ Bcog = getCOG(points2);

		int len = points1.size();
		if(len == 0) return;

		vector<XYZ> listA;
		vector<XYZ> listB;
		for(int i=0;i<len;i++){
			XYZ a = points1[i] - Acog;
			XYZ b = points2[i] - Bcog;
			listA.push_back(a);
			listB.push_back(b);
		}

		TransForm tf = buildRotation(listA, listB);

		for(int i=0;i<nB*3;i++){
			XYZ localT = this->initListB[i] - Bcog;
			this->transedB[i] = tf.transform(localT) + Acog;
		}
	}

	void generateAlign(){
		for(int i=0;i<nA;i++){
			//if(i > 1) continue;
			LocalFrame csA = this->csAList[i];
			for(int j=0;j<nB;j++){
				//if(j > 1) continue;
				LocalFrame csB = this->csBList[j];

				//cout << "csMove: " << endl;
			//	cm1.print();
				//cout << "update move" << endl;

				alignAtPosition(i, j);
				calculateRMScore();

				getTransformFromAlign();
				calculateRMScore();

				getTransformFromAlign();
				calculateRMScore();

				/*
				CsMove cm3 = getCMFromAlign();
				cout << "csMove: " << endl;
				cm3.print();
				double rm3 = updateCsMove(cm3);

				cout << "round3" << endl;
				printResult();
				CsMove cm4 = getCMFromAlign();
				cout << "csMove: " << endl;
				cm4.print();
				double rm4= updateCsMove(cm4);

				cout << "round4" << endl;
				printResult();
				*/

			}
		}
	}

	double getRMScore() {
		return this->bestRM;
	}

	void printResult(){

		int n = 0;
		for(int i=0;i<this->nA;i++){
			if(bestAlign[i] >= 0){
				n++;
			}
		}
		printf("RMscore= %5.3f nA= %2d nB= %2d d0= %5.3f\n", this->bestRM, nA, nB, d0);
	}

	virtual ~MotifAlign();
};

} /* namespace NSPalign */

#endif /* MOTIF_MOTIFALIGN_H_ */
