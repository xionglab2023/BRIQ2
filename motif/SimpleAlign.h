/*
 * SimpleAlign.h
 *
 *  Created on: 2023年8月25日
 *      Author: nuc
 */

/*
 * SimpleAlign.h
 *
 *  Created on: 2021年1月4日
 *      Author: pengx
 */

#ifndef MOTIF_SIMPLEALIGN_H_
#define MOTIF_SIMPLEALIGN_H_

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

class SimpleAlign {
public:

	int nA;
	int nB;

	int* align;
	double* alignDD;

	double d0;
	double xLN;


	vector<LocalFrame> csAList;
	vector<LocalFrame> csBList;
	vector<LocalFrame> csBListTrans;

	XYZ* initListA;
	XYZ* initListB;
	XYZ* transedB;



	SimpleAlign(RNAPDB* pdbA, RNAPDB* pdbB){
		AtomLib* atLib = new AtomLib();
		vector<RNABase*> baseListA = pdbA->getValidBaseList(atLib);
		vector<RNABase*> baseListB = pdbB->getValidBaseList(atLib);

		if(baseListA.size() < baseListB.size()){
			baseListA = pdbB->getValidBaseList(atLib);
			baseListB = pdbA->getValidBaseList(atLib);
		}

		this->nA = baseListA.size();
		this->nB = baseListB.size();

		if(nA != nB){
			cout << "sequence length not equal " << nA << " " << nB << endl;
			exit(0);
		}

		double LN = 0.5*(nA+nB);
		this->xLN = 1.0/LN;

		this->d0 = (nA+nB)/22.3 - 0.28;
		if(d0 < 0.05)
			d0 = 0.05;

		this->d0 = d0;
		/*
		cout << "nA: " << nA << endl;
		cout << "nB: " << nB << endl;
		cout << "LN: " << LN << endl;
		cout << "d0: " << d0 << endl;
		*/


		this->align = new int[nA];
		this->alignDD = new double[nA];
		this->initListA = new XYZ[nA*3];
		this->initListB = new XYZ[nB*3];
		this->transedB = new XYZ[nB*3];


		XYZ xa(0 , 0 , 0);
		XYZ xb(3.464 , 2 , 0);
		XYZ xc(3.464 , -2 , 0);

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

	void getTransformFromAlign(){
		vector<XYZ> points1;
		vector<XYZ> points2;

		for(int i=0;i<nA;i++){

			points1.push_back(initListA[i*3]);
			points1.push_back(initListA[i*3+1]);
			points1.push_back(initListA[i*3+2]);

			points2.push_back(initListB[i*3]);
			points2.push_back(initListB[i*3+1]);
			points2.push_back(initListB[i*3+2]);

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

	double getRMScore(){

		getTransformFromAlign();
		double gdt = 0;
		double xD0D0 = 1.0/d0/d0;
		double rmsd = 0;
		double dd;
		int n1 = 0;
		int n2 = 0;
		int n3 = 0;
		int n4 = 0;
		for(int i=0;i<nA;i++){
			dd = squareDistance(initListA[i*3], transedB[i*3]) + squareDistance(initListA[i*3+1], transedB[i*3+1]) + squareDistance(initListA[i*3+2], transedB[i*3+2]);
			dd = sqrt(dd*0.333333333);
			if(dd < 0.5)
				n1 ++;
			if(dd < 1.0)
				n2 ++;
			if(dd < 2.0)
				n3++;
			if(dd < 4.0)
				n4++;
		}

		gdt = (n1+n2+n3+n4)*0.25/nA;
		return gdt;
	}

	void printResult(){

		double RG = 0;
		double xD0D0 = 1.0/d0/d0;
		double rmsd = 0;
		double dd;
		for(int i=0;i<nA;i++){
			dd = squareDistance(initListA[i*3], transedB[i*3]) + squareDistance(initListA[i*3+1], transedB[i*3+1]) + squareDistance(initListA[i*3+2], transedB[i*3+2]);
			RG += xLN*(1.0/(1+dd*0.3333333*xD0D0));
			rmsd += dd;
		}

		rmsd = sqrt(rmsd/nA/3.0);

		printf("GDT= %5.3f RMSD= %6.3f\n", getRMScore(), rmsd);
	}

	virtual ~SimpleAlign();
};

} /* namespace NSPalign */

#endif /* MOTIF_SIMPLEALIGN_H_ */
