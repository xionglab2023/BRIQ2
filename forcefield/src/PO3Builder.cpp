/*
 * PO3Builder.cpp
 *
 */

#include <forcefield/PO3Builder.h>

namespace NSPforcefield {

PO3Builder::PO3Builder(ForceFieldPara* para) {
	this->rotLib = new PhosphateRotamerLib();
	XYZ c2,c3,o3,p,op1,op2,o5,c5,c4,o4;

	double len1 = 1.605;
	double len2 = 1.592;
	double ang1 = 120.1;
	double ang2 = 103.5;

	this->para = para;

	for(int i=0;i<360;i++){
		for(int j=0;j<360;j++){
			LocalFrame cs0;
			LocalFrame cs1 = cs0.csNext(len1, ang1, i+0.5);
			LocalFrame cs2 = cs1.csNext(len2, ang2, j+0.5);
			XYZ p = cs1.origin_;
			XYZ o5 = cs2.origin_;
			XYZ op1 = local2global(cs2, XYZ(-2.062, 0.570, 1.278));
			XYZ op2 = local2global(cs2, XYZ(-2.054, 0.589, -1.275));

			this->pList.push_back(p);
			this->o5List.push_back(o5);
			this->op1List.push_back(op1);
			this->op2List.push_back(op2);
		}
	}


	string path = NSPdataio::datapath()+"/";
	ifstream file;
	string s;
	string fileName = path+"rnaDihedEnergy/impD1D2.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	vector<string> spt;
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eImpD1D2.push_back(energyRescale(atof(spt[i].c_str())));
		}
	}
	file.close();

	fileName = path+"rnaDihedEnergy/impD4D5.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eImpD4D5.push_back(energyRescale(atof(spt[i].c_str())));
		}
	}
	file.close();

	fileName = path+"rnaDihedEnergy/D2D4D3.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eD2D4D3.push_back(energyRescale(atof(spt[i].c_str())));
		}
	}
	file.close();


	c2 = XYZ(-2.008, -1.402, 0);
	c3 = XYZ(-1.418, 0, 0);
	o3 = XYZ(0, 0, 0);

	fileName = path+"d1d2Lib/libA/rot-level1.txt";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	getline(file,s);
	double x,y;
	while(getline(file,s)){
		NSPtools::splitString(s," ",&spt);
		x = atof(spt[0].c_str());
		y = atof(spt[1].c_str());
		this->d1d2Lib1A.push_back(XYZ(x,y,0));
	}
	file.close();

	for(int i=0;i<d1d2Lib1A.size();i++){
		char ss[20];
		sprintf(ss, "%d", i);
		fileName = path + "d1d2Lib/libA/rot-level2-"+string(ss)+".txt";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()){
			cout << "fail to open file: " << fileName << endl;
			exit(1);
		}
		getline(file,s);
		NSPtools::splitString(s," ",&spt);
		double error = atof(spt[2].c_str());
		this->lib2ErrorA.push_back(error);
		double x,y;
		vector<XYZ> subList;
		while(getline(file,s)){
			NSPtools::splitString(s," ",&spt);
			x = atof(spt[0].c_str());
			y = atof(spt[1].c_str());
			subList.push_back(XYZ(x,y,0));
		}
		file.close();
		this->d1d2Lib2A.push_back(subList);
	}

	fileName = path+"d1d2Lib/libB/rot-level1.txt";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	getline(file,s);
	while(getline(file,s)){
		NSPtools::splitString(s," ",&spt);
		x = atof(spt[0].c_str());
		y = atof(spt[1].c_str());
		this->d1d2Lib1B.push_back(XYZ(x,y,0));
	}
	file.close();

	for(int i=0;i<d1d2Lib1B.size();i++){
		char ss[20];
		sprintf(ss, "%d", i);
		fileName = path + "d1d2Lib/libB/rot-level2-"+string(ss)+".txt";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()){
			cout << "fail to open file: " << fileName << endl;
			exit(1);
		}
		getline(file,s);
		NSPtools::splitString(s," ",&spt);
		double error = atof(spt[2].c_str());
		this->lib2ErrorB.push_back(error);
		vector<XYZ> subList;
		while(getline(file,s)){
			NSPtools::splitString(s," ",&spt);
			x = atof(spt[0].c_str());
			y = atof(spt[1].c_str());
			subList.push_back(XYZ(x,y,0));
		}
		file.close();
		this->d1d2Lib2B.push_back(subList);
	}

	fileName = path+"dihed.wt.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	double e1,e2,e3,e4,e5;

	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		e1 = atof(spt[1].c_str());
		e2 = atof(spt[2].c_str());
		e3 = atof(spt[3].c_str());
		e4 = atof(spt[4].c_str());
		e5 = atof(spt[5].c_str());

		eDihed1.push_back(e1);
		eDihed2.push_back(e2);
		eDihed3.push_back(e3);
		eDihed4.push_back(e4);
		eDihed5.push_back(e5);
	}
	eDihed1.push_back(e1);
	eDihed2.push_back(e2);
	eDihed3.push_back(e3);
	eDihed4.push_back(e4);
	eDihed5.push_back(e5);

	file.close();
}

void PO3Builder::buildPhosphate(RiboseConformer* riboConfA, RiboseConformer* riboConfB, PhosphateConformer* outPhoConf){

	XYZ o3,p,op1,op2,o5,c5,c4,o4, o2;
	LocalFrame cs2A = riboConfA->cs2;

	o3 = global2local(cs2A, riboConfA->coords[5]);
	o2 = global2local(cs2A, riboConfA->coords[7]);

	//distance between atom O3' and C5'
	double d0 = cs2A.origin_.distance(riboConfB->coords[6]);

	LocalFrame cs0;
	c5 = global2local(cs2A, riboConfB->coords[6]);
	c4 = global2local(cs2A, riboConfB->coords[3]);
	o4 = global2local(cs2A, riboConfB->coords[4]);

	LocalFrame cs1;
	LocalFrame cs2;


	LocalFrame bestCs1, bestCs2;
	int bestDihed1, bestDihed2;
	int localBestD1, localBestD2;


	double minE = 999999999.9;
	double e, u;

	int dihed1, dihed2;
	int idA, idB;

	double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5;

	int impIndexA = (int)((riboConfA->rot->improper+60)*0.166666667);
	if(impIndexA < 0) impIndexA = 0;
	if(impIndexA > 19) impIndexA = 19;

	char impType = 'A';
	if(riboConfA->rot->improper > 0) impType = 'B';

	int impIndexB = (int)((riboConfB->rot->improper+60)*0.166666667);
	if(impIndexB < 0) impIndexB = 0;
	if(impIndexB > 19) impIndexB = 19;

	int indexDD;

	int d1d2LibSize = this->d1d2Lib1A.size();
	

	int bestIndex1=0;

	int regionIndexA;
	int regionIndexB;
	int regionIndexC;

	if(d0 > 4.5){
		//If the distance between atom O3' and C5' is larger than 4.5 angstrom, we use default dihedral angles to build PO3
		dihed1 = 100;
		dihed2 = 290;
		indexDD = (dihed1 * 360) + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];
		op1 = this->op1List[indexDD];
		op2 = this->op2List[indexDD];

		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;

		if(xdihed3 < 0 || xdihed3 >= 360 || isnan(xdihed3)) {
			cout << "xdihed3: " << xdihed3 << endl;
			cout << "o3: " << o3.toString() << endl;
			cout << "p: " << p.toString() << endl;
			cout << "o5: " << o5.toString() << endl;
			cout << "c5: " << c5.toString() << endl;
			exit(0);
		}

		if(xdihed4 < 0 || xdihed4 >= 360 || isnan(xdihed4)) {
			cout << "xdihed4: " << xdihed4 << endl;
			cout << "p: " << p.toString() << endl;
			cout << "o5: " << o5.toString() << endl;
			cout << "c5: " << c5.toString() << endl;
			cout << "c4: " << c4.toString() << endl;
			exit(0);
		}

		if(xdihed5 < 0 || xdihed5 >= 360 || isnan(xdihed5)) {
			cout << "xdihed5: " << xdihed5 << endl;
		
			cout << "o5: " << o5.toString() << endl;
			cout << "c5: " << c5.toString() << endl;
			cout << "c4: " << c4.toString() << endl;
			cout << "o4: " << o4.toString() << endl;
			exit(0);
		}

		e = 0;
		u = (xd3-len3)*para->rnaKBond;

		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;

		

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (op1.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		u = (op2.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;
		


		if(impIndexA < 10 && dihed2 > 250)
			regionIndexA = 0;
		else if(impIndexA < 10 && dihed1 < 210)
			regionIndexA = 1;
		else if(impIndexA < 10)
			regionIndexA = 2;
		else if(dihed2 > 250)
			regionIndexA = 3;
		else if(dihed1 < 210)
			regionIndexA = 4;
		else
			regionIndexA = 5;


		if(impIndexB < 10 && xdihed5 > 240)
			regionIndexB = 0;
		else if(impIndexB < 10 && xdihed5 > 120)
			regionIndexB = 1;
		else if(impIndexB < 10)
			regionIndexB = 2;
		else if(xdihed5 > 240)
			regionIndexB = 3;
		else if(xdihed5 > 120)
			regionIndexB = 4;
		else
			regionIndexB = 5;

		if(dihed2 > 250 && xdihed3 > 240)
			regionIndexC = 0;
		else if(dihed2 > 250 && xdihed3 > 120)
			regionIndexC = 1;
		else if(dihed2 > 250)
			regionIndexC = 2;
		else if(xdihed3 > 240)
			regionIndexC = 3;
		else if(xdihed3 > 120)
			regionIndexC = 4;
		else
			regionIndexC = 5;

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];


		if(e < minE){
			minE = e;
		}

		outPhoConf->updateLocalFrameAndRotamer(cs2A, rotLib->prLib[dihed1][dihed2], minE*para->wtPho*para->connectRescale);
		return;
	}

	for(int i=0;i<d1d2LibSize;i++){
		if(impType == 'A') {
			dihed1 = (int)d1d2Lib1A[i].x_;
			dihed2 = (int)d1d2Lib1A[i].y_;
		}
		else {
			dihed1 = (int)d1d2Lib1B[i].x_;
			dihed2 = (int)d1d2Lib1B[i].y_;
		}
		indexDD = (dihed1 * 360) + dihed2;
		//cout << "indexDD: " << indexDD << endl;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];
		op1 = this->op1List[indexDD];
		op2 = this->op2List[indexDD];

		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;

		if(xdihed3 < 0 || xdihed3 >= 360 || isnan(xdihed3)) {
			cout << "xdihed3: " << xdihed3 << endl;
			cout << "o3: " << o3.toString() << endl;
			cout << "p: " << p.toString() << endl;
			cout << "o5: " << o5.toString() << endl;
			cout << "c5: " << c5.toString() << endl;
			exit(0);
		}

		if(xdihed4 < 0 || xdihed4 >= 360 || isnan(xdihed4)) {
			cout << "xdihed4: " << xdihed4 << endl;
			cout << "p: " << p.toString() << endl;
			cout << "o5: " << o5.toString() << endl;
			cout << "c5: " << c5.toString() << endl;
			cout << "c4: " << c4.toString() << endl;
			exit(0);
		}

		if(xdihed5 < 0 || xdihed5 >= 360 || isnan(xdihed5)) {
			cout << "xdihed5: " << xdihed5 << endl;
		
			cout << "o5: " << o5.toString() << endl;
			cout << "c5: " << c5.toString() << endl;
			cout << "c4: " << c4.toString() << endl;
			cout << "o4: " << o4.toString() << endl;
			exit(0);
		}

		e = 0;
		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;

		

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (op1.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		u = (op2.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		if(impIndexA < 10 && dihed2 > 250)
			regionIndexA = 0;
		else if(impIndexA < 10 && dihed1 < 210)
			regionIndexA = 1;
		else if(impIndexA < 10)
			regionIndexA = 2;
		else if(dihed2 > 250)
			regionIndexA = 3;
		else if(dihed1 < 210)
			regionIndexA = 4;
		else
			regionIndexA = 5;


		if(impIndexB < 10 && xdihed5 > 240)
			regionIndexB = 0;
		else if(impIndexB < 10 && xdihed5 > 120)
			regionIndexB = 1;
		else if(impIndexB < 10)
			regionIndexB = 2;
		else if(xdihed5 > 240)
			regionIndexB = 3;
		else if(xdihed5 > 120)
			regionIndexB = 4;
		else
			regionIndexB = 5;

		if(dihed2 > 250 && xdihed3 > 240)
			regionIndexC = 0;
		else if(dihed2 > 250 && xdihed3 > 120)
			regionIndexC = 1;
		else if(dihed2 > 250)
			regionIndexC = 2;
		else if(xdihed3 > 240)
			regionIndexC = 3;
		else if(xdihed3 > 120)
			regionIndexC = 4;
		else
			regionIndexC = 5;

		
		/*
			cout << "impIndexA: " << impIndexA << endl;
			cout << "impIndexB: " << impIndexB << endl;
			cout << "dihed1: " << dihed1 << endl;
			cout << "dihed2: " << dihed2 << endl;
			cout << "dihed3: " << xdihed3 << endl;
			cout << "dihed4: " << xdihed4 << endl;
			cout << "dihed5: " << xdihed5 << endl;
		*/

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			bestIndex1 = i;
			minE = e;
		}

	}


	if(d0 > 3.9){
		//If the distance between atom O3' and C5' is larger than 3.8 angstrom, we use default dihedral angles to build PO3
		outPhoConf->updateLocalFrameAndRotamer(cs2A, rotLib->prLib[bestDihed1][bestDihed2], minE*para->wtPho*para->connectRescale);
		return;
	}

	for(int i=0;i<d1d2LibSize;i++){
		if(impType == 'A') {
			dihed1 = d1d2Lib2A[bestIndex1][i].x_;
			dihed2 = d1d2Lib2A[bestIndex1][i].y_;
		}
		else {
			dihed1 = d1d2Lib2B[bestIndex1][i].x_;
			dihed2 = d1d2Lib2B[bestIndex1][i].y_;
		}
		indexDD = dihed1 * 360 + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];
		op1 = this->op1List[indexDD];
		op2 = this->op2List[indexDD];
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;


		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;
		

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (op1.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		u = (op2.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		if(impIndexA < 10 && dihed2 > 250)
			regionIndexA = 0;
		else if(impIndexA < 10 && dihed1 < 210)
			regionIndexA = 1;
		else if(impIndexA < 10)
			regionIndexA = 2;
		else if(dihed2 > 250)
			regionIndexA = 3;
		else if(dihed1 < 210)
			regionIndexA = 4;
		else
			regionIndexA = 5;


		if(impIndexB < 10 && xdihed5 > 240)
			regionIndexB = 0;
		else if(impIndexB < 10 && xdihed5 > 120)
			regionIndexB = 1;
		else if(impIndexB < 10)
			regionIndexB = 2;
		else if(xdihed5 > 240)
			regionIndexB = 3;
		else if(xdihed5 > 120)
			regionIndexB = 4;
		else
			regionIndexB = 5;

		if(dihed2 > 250 && xdihed3 > 240)
			regionIndexC = 0;
		else if(dihed2 > 250 && xdihed3 > 120)
			regionIndexC = 1;
		else if(dihed2 > 250)
			regionIndexC = 2;
		else if(xdihed3 > 240)
			regionIndexC = 3;
		else if(xdihed3 > 120)
			regionIndexC = 4;
		else
			regionIndexC = 5;

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];


		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			minE = e;
		}

	}


	double libErr;
	if(impType == 'A')
		libErr = lib2ErrorA[bestDihed1*d1d2LibSize+bestDihed2];
	else
		libErr = lib2ErrorB[bestDihed1*d1d2LibSize+bestDihed2];

	if(libErr < 2.5) {
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<5;i++){
			dihed1 = localBestD1 + (i-2);
			if(dihed1 >= 360 || dihed1 < 0) continue;
			for(int j=0;j<5;j++){
				dihed2 = localBestD2 + (j-2);
				if(dihed2 >= 360 || dihed2 < 0) continue;
				indexDD = dihed1 * 360 + dihed2;
				p = this->pList[indexDD];
				o5 = this->o5List[indexDD];
				op1 = this->op1List[indexDD];
				op2 = this->op2List[indexDD];
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;

				
			if(xdihed3 < 0 || xdihed3 >= 360 || isnan(xdihed3)) {
				cout << "xdihed3: " << xdihed3 << endl;
				cout << "o3: " << o3.toString() << endl;
				cout << "p: " << p.toString() << endl;
				cout << "o5: " << o5.toString() << endl;
				cout << "c5: " << c5.toString() << endl;
				exit(0);
			}

			if(xdihed4 < 0 || xdihed4 >= 360 || isnan(xdihed4)) {
				cout << "xdihed4: " << xdihed4 << endl;
				cout << "p: " << p.toString() << endl;
				cout << "o5: " << o5.toString() << endl;
				cout << "c5: " << c5.toString() << endl;
				cout << "c4: " << c4.toString() << endl;
				exit(0);
			}

			if(xdihed5 < 0 || xdihed5 >= 360 || isnan(xdihed5)) {
				cout << "xdihed5: " << xdihed5 << endl;
		
				cout << "o5: " << o5.toString() << endl;
				cout << "c5: " << c5.toString() << endl;
				cout << "c4: " << c4.toString() << endl;
				cout << "o4: " << o4.toString() << endl;
				exit(0);
			}


				e = 0;
				u = (xd3-len3)*para->rnaKBond;
				if(u < -1)
					e += (-2*u-1);
				else if(u > 1)
					e += 2*u-1;
				else if(u <= 1)
					e += u;
				else 
					e += -u;

				u = (xang3-ang3)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				u = (op1.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				u = (op2.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				if(impIndexA < 10 && dihed2 > 250)
					regionIndexA = 0;
				else if(impIndexA < 10 && dihed1 < 210)
					regionIndexA = 1;
				else if(impIndexA < 10)
					regionIndexA = 2;
				else if(dihed2 > 250)
					regionIndexA = 3;
				else if(dihed1 < 210)
					regionIndexA = 4;
				else
					regionIndexA = 5;


				if(impIndexB < 10 && xdihed5 > 240)
					regionIndexB = 0;
				else if(impIndexB < 10 && xdihed5 > 120)
					regionIndexB = 1;
				else if(impIndexB < 10)
					regionIndexB = 2;
				else if(xdihed5 > 240)
					regionIndexB = 3;
				else if(xdihed5 > 120)
					regionIndexB = 4;
				else
					regionIndexB = 5;

				if(dihed2 > 250 && xdihed3 > 240)
					regionIndexC = 0;
				else if(dihed2 > 250 && xdihed3 > 120)
					regionIndexC = 1;
				else if(dihed2 > 250)
					regionIndexC = 2;
				else if(xdihed3 > 240)
					regionIndexC = 3;
				else if(xdihed3 > 120)
					regionIndexC = 4;
				else
					regionIndexC = 5;

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}

			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		outPhoConf->updateLocalFrameAndRotamer(cs2A, rotLib->prLib[bestDihed1][bestDihed2], minE*para->wtPho*para->connectRescale);

		return;

	}
	else {
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<9;i++){
			dihed1 = localBestD1 + (i-4)*2;
			if(dihed1 >= 360 || dihed1 < 0) continue;
			for(int j=0;j<9;j++){
				dihed2 = localBestD2 + (j-4)*2;
				if(dihed2 >= 360 || dihed2 < 0) continue;
				indexDD = dihed1 * 360 + dihed2;
				p = this->pList[indexDD];
				o5 = this->o5List[indexDD];
				op1 = this->op1List[indexDD];
				op2 = this->op2List[indexDD];
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;

				if(xdihed3 < 0 || xdihed3 >= 360 || isnan(xdihed3)) {
					cout << "xdihed3: " << xdihed3 << endl;
					cout << "o3: " << o3.toString() << endl;
					cout << "p: " << p.toString() << endl;
					cout << "o5: " << o5.toString() << endl;
					cout << "c5: " << c5.toString() << endl;
					exit(0);
				}

				if(xdihed4 < 0 || xdihed4 >= 360 || isnan(xdihed4)) {
					cout << "xdihed4: " << xdihed4 << endl;
					cout << "p: " << p.toString() << endl;
					cout << "o5: " << o5.toString() << endl;
					cout << "c5: " << c5.toString() << endl;
					cout << "c4: " << c4.toString() << endl;
					exit(0);
				}

				if(xdihed5 < 0 || xdihed5 >= 360 || isnan(xdihed5)) {
					cout << "xdihed5: " << xdihed5 << endl;
		
					cout << "o5: " << o5.toString() << endl;
					cout << "c5: " << c5.toString() << endl;
					cout << "c4: " << c4.toString() << endl;
					cout << "o4: " << o4.toString() << endl;
					exit(0);
				}

				e = 0;
				u = (xd3-len3)*para->rnaKBond;
				if(u < -1)
					e += (-2*u-1);
				else if(u > 1)
					e += 2*u-1;
				else if(u <= 1)
					e += u;
				else 
					e += -u;
				u = (xang3-ang3)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				u = (op1.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				u = (op2.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				if(impIndexA < 10 && dihed2 > 250)
					regionIndexA = 0;
				else if(impIndexA < 10 && dihed1 < 210)
					regionIndexA = 1;
				else if(impIndexA < 10)
					regionIndexA = 2;
				else if(dihed2 > 250)
					regionIndexA = 3;
				else if(dihed1 < 210)
					regionIndexA = 4;
				else
					regionIndexA = 5;


				if(impIndexB < 10 && xdihed5 > 240)
					regionIndexB = 0;
				else if(impIndexB < 10 && xdihed5 > 120)
					regionIndexB = 1;
				else if(impIndexB < 10)
					regionIndexB = 2;
				else if(xdihed5 > 240)
					regionIndexB = 3;
				else if(xdihed5 > 120)
					regionIndexB = 4;
				else
					regionIndexB = 5;

				if(dihed2 > 250 && xdihed3 > 240)
					regionIndexC = 0;
				else if(dihed2 > 250 && xdihed3 > 120)
					regionIndexC = 1;
				else if(dihed2 > 250)
					regionIndexC = 2;
				else if(xdihed3 > 240)
					regionIndexC = 3;
				else if(xdihed3 > 120)
					regionIndexC = 4;
				else
					regionIndexC = 5;

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		outPhoConf->updateLocalFrameAndRotamer(cs2A, rotLib->prLib[bestDihed1][bestDihed2], minE*para->wtPho*para->connectRescale);
		return;
	}

}

double PO3Builder::getEnergy(RiboseConformer* riboConfA, RiboseConformer* riboConfB){
	XYZ o3,p,op1,op2,o5,c5,c4,o4, o2;

	LocalFrame cs2A = riboConfA->cs2;

	o3 = global2local(cs2A, riboConfA->coords[5]);
	o2 = global2local(cs2A, riboConfA->coords[7]);

	//distance between atom O3' and C5'
	double d0 = cs2A.origin_.distance(riboConfB->coords[6]);

	LocalFrame cs0;
	c5 = global2local(cs2A, riboConfB->coords[6]);
	c4 = global2local(cs2A, riboConfB->coords[3]);
	o4 = global2local(cs2A, riboConfB->coords[4]);
	LocalFrame cs1;
	LocalFrame cs2;


	LocalFrame bestCs1, bestCs2;
	int bestDihed1, bestDihed2;
	int localBestD1, localBestD2;


	double minE = 999999999.9;
	double e, u;

	int dihed1, dihed2;
	int idA, idB;


	double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5;



	int impIndexA = (int)((riboConfA->rot->improper+60)*0.166666667);
	if(impIndexA < 0) impIndexA = 0;
	if(impIndexA > 19) impIndexA = 19;

	char impType = 'A';
	if(riboConfA->rot->improper > 0) impType = 'B';

	int impIndexB = (int)((riboConfB->rot->improper+60)*0.166666667);
	if(impIndexB < 0) impIndexB = 0;
	if(impIndexB > 19) impIndexB = 19;

	int regionIndexA, regionIndexB, regionIndexC;


	int indexDD;

	int d1d2LibSize = this->d1d2Lib1A.size();
	int bestIndex1=0;

	if(d0 > 4.5){
		//If the distance between atom O3' and C5' is larger than 4.5 angstrom, we use default dihedral angles to build PO3
		dihed1 = 100;
		dihed2 = 290;
		indexDD = (dihed1 * 360) + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];
		op1 = this->op1List[indexDD];
		op2 = this->op2List[indexDD];

		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;
		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (op1.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		u = (op2.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		if(impIndexA < 10 && dihed2 > 250)
			regionIndexA = 0;
		else if(impIndexA < 10 && dihed1 < 210)
			regionIndexA = 1;
		else if(impIndexA < 10)
			regionIndexA = 2;
		else if(dihed2 > 250)
			regionIndexA = 3;
		else if(dihed1 < 210)
			regionIndexA = 4;
		else
			regionIndexA = 5;


		if(impIndexB < 10 && xdihed5 > 240)
			regionIndexB = 0;
		else if(impIndexB < 10 && xdihed5 > 120)
			regionIndexB = 1;
		else if(impIndexB < 10)
			regionIndexB = 2;
		else if(xdihed5 > 240)
			regionIndexB = 3;
		else if(xdihed5 > 120)
			regionIndexB = 4;
		else
			regionIndexB = 5;

		if(dihed2 > 250 && xdihed3 > 240)
			regionIndexC = 0;
		else if(dihed2 > 250 && xdihed3 > 120)
			regionIndexC = 1;
		else if(dihed2 > 250)
			regionIndexC = 2;
		else if(xdihed3 > 240)
			regionIndexC = 3;
		else if(xdihed3 > 120)
			regionIndexC = 4;
		else
			regionIndexC = 5;

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

		return minE*para->wtPho*para->connectRescale;
	}



	for(int i=0;i<d1d2LibSize;i++){
		if(impType == 'A') {
			dihed1 = (int)d1d2Lib1A[i].x_;
			dihed2 = (int)d1d2Lib1A[i].y_;
		}
		else {
			dihed1 = (int)d1d2Lib1B[i].x_;
			dihed2 = (int)d1d2Lib1B[i].y_;
		}
		indexDD = (dihed1 * 360) + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];
		op1 = this->op1List[indexDD];
		op2 = this->op2List[indexDD];

		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;
		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (op1.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		u = (op2.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		if(impIndexA < 10 && dihed2 > 250)
			regionIndexA = 0;
		else if(impIndexA < 10 && dihed1 < 210)
			regionIndexA = 1;
		else if(impIndexA < 10)
			regionIndexA = 2;
		else if(dihed2 > 250)
			regionIndexA = 3;
		else if(dihed1 < 210)
			regionIndexA = 4;
		else
			regionIndexA = 5;


		if(impIndexB < 10 && xdihed5 > 240)
			regionIndexB = 0;
		else if(impIndexB < 10 && xdihed5 > 120)
			regionIndexB = 1;
		else if(impIndexB < 10)
			regionIndexB = 2;
		else if(xdihed5 > 240)
			regionIndexB = 3;
		else if(xdihed5 > 120)
			regionIndexB = 4;
		else
			regionIndexB = 5;

		if(dihed2 > 250 && xdihed3 > 240)
			regionIndexC = 0;
		else if(dihed2 > 250 && xdihed3 > 120)
			regionIndexC = 1;
		else if(dihed2 > 250)
			regionIndexC = 2;
		else if(xdihed3 > 240)
			regionIndexC = 3;
		else if(xdihed3 > 120)
			regionIndexC = 4;
		else
			regionIndexC = 5;

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			bestIndex1 = i;
			minE = e;
		}
	}

	if(d0 > 3.9){
		//If the distance between atom O3' and C5' is larger than 3.8 angstrom, we use default dihedral angles to build PO3
		return minE*para->wtPho*para->connectRescale;
	}

	for(int i=0;i<d1d2LibSize;i++){
		if(impType == 'A') {
			dihed1 = d1d2Lib2A[bestIndex1][i].x_;
			dihed2 = d1d2Lib2A[bestIndex1][i].y_;
		}
		else {
			dihed1 = d1d2Lib2B[bestIndex1][i].x_;
			dihed2 = d1d2Lib2B[bestIndex1][i].y_;
		}
		indexDD = dihed1 * 360 + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];
		op1 = this->op1List[indexDD];
		op2 = this->op2List[indexDD];
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;
		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;
		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (op1.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		u = (op2.distance(o2) - 2.7);
		if(u < 0)
			e += u*u*100;

		if(impIndexA < 10 && dihed2 > 250)
			regionIndexA = 0;
		else if(impIndexA < 10 && dihed1 < 210)
			regionIndexA = 1;
		else if(impIndexA < 10)
			regionIndexA = 2;
		else if(dihed2 > 250)
			regionIndexA = 3;
		else if(dihed1 < 210)
			regionIndexA = 4;
		else
			regionIndexA = 5;


		if(impIndexB < 10 && xdihed5 > 240)
			regionIndexB = 0;
		else if(impIndexB < 10 && xdihed5 > 120)
			regionIndexB = 1;
		else if(impIndexB < 10)
			regionIndexB = 2;
		else if(xdihed5 > 240)
			regionIndexB = 3;
		else if(xdihed5 > 120)
			regionIndexB = 4;
		else
			regionIndexB = 5;

		if(dihed2 > 250 && xdihed3 > 240)
			regionIndexC = 0;
		else if(dihed2 > 250 && xdihed3 > 120)
			regionIndexC = 1;
		else if(dihed2 > 250)
			regionIndexC = 2;
		else if(xdihed3 > 240)
			regionIndexC = 3;
		else if(xdihed3 > 120)
			regionIndexC = 4;
		else
			regionIndexC = 5;

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			minE = e;
		}
	}

	double libErr;
	if(impType == 'A')
		libErr = lib2ErrorA[bestDihed1*d1d2LibSize+bestDihed2];
	else
		libErr = lib2ErrorB[bestDihed1*d1d2LibSize+bestDihed2];

	if(libErr < 2.5) {
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<5;i++){
			dihed1 = localBestD1 + (i-2);
			if(dihed1 >= 360 || dihed1 < 0) continue;
			for(int j=0;j<5;j++){
				dihed2 = localBestD2 + (j-2);
				if(dihed2 >= 360 || dihed2 < 0) continue;
				indexDD = dihed1 * 360 + dihed2;
				p = this->pList[indexDD];
				o5 = this->o5List[indexDD];
				op1 = this->op1List[indexDD];
				op2 = this->op2List[indexDD];
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;
				e = 0;
				u = (xd3-len3)*para->rnaKBond;
				if(u < -1)
					e += (-2*u-1);
				else if(u > 1)
					e += 2*u-1;
				else if(u <= 1)
					e += u;
				else 
					e += -u;

				u = (xang3-ang3)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				u = (op1.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				u = (op2.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				if(impIndexA < 10 && dihed2 > 250)
					regionIndexA = 0;
				else if(impIndexA < 10 && dihed1 < 210)
					regionIndexA = 1;
				else if(impIndexA < 10)
					regionIndexA = 2;
				else if(dihed2 > 250)
					regionIndexA = 3;
				else if(dihed1 < 210)
					regionIndexA = 4;
				else
					regionIndexA = 5;


				if(impIndexB < 10 && xdihed5 > 240)
					regionIndexB = 0;
				else if(impIndexB < 10 && xdihed5 > 120)
					regionIndexB = 1;
				else if(impIndexB < 10)
					regionIndexB = 2;
				else if(xdihed5 > 240)
					regionIndexB = 3;
				else if(xdihed5 > 120)
					regionIndexB = 4;
				else
					regionIndexB = 5;

				if(dihed2 > 250 && xdihed3 > 240)
					regionIndexC = 0;
				else if(dihed2 > 250 && xdihed3 > 120)
					regionIndexC = 1;
				else if(dihed2 > 250)
					regionIndexC = 2;
				else if(xdihed3 > 240)
					regionIndexC = 3;
				else if(xdihed3 > 120)
					regionIndexC = 4;
				else
					regionIndexC = 5;

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		return minE*para->wtPho*para->connectRescale;

	}
	else {
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<9;i++){
			dihed1 = localBestD1 + (i-4)*2;
			if(dihed1 >= 360 || dihed1 < 0) continue;
			for(int j=0;j<9;j++){
				dihed2 = localBestD2 + (j-4)*2;
				if(dihed2 >= 360 || dihed2 < 0) continue;
				indexDD = dihed1 * 360 + dihed2;
				p = this->pList[indexDD];
				o5 = this->o5List[indexDD];
				op1 = this->op1List[indexDD];
				op2 = this->op2List[indexDD];
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;
				e = 0;
				u = (xd3-len3)*para->rnaKBond;
				if(u < -1)
					e += (-2*u-1);
				else if(u > 1)
					e += 2*u-1;
				else if(u <= 1)
					e += u;
				else 
					e += -u;

				u = (xang3-ang3)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->rnaKAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				u = (op1.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				u = (op2.distance(o2) - 2.7);
				if(u < 0)
					e += u*u*100;

				if(impIndexA < 10 && dihed2 > 250)
					regionIndexA = 0;
				else if(impIndexA < 10 && dihed1 < 210)
					regionIndexA = 1;
				else if(impIndexA < 10)
					regionIndexA = 2;
				else if(dihed2 > 250)
					regionIndexA = 3;
				else if(dihed1 < 210)
					regionIndexA = 4;
				else
					regionIndexA = 5;


				if(impIndexB < 10 && xdihed5 > 240)
					regionIndexB = 0;
				else if(impIndexB < 10 && xdihed5 > 120)
					regionIndexB = 1;
				else if(impIndexB < 10)
					regionIndexB = 2;
				else if(xdihed5 > 240)
					regionIndexB = 3;
				else if(xdihed5 > 120)
					regionIndexB = 4;
				else
					regionIndexB = 5;

				if(dihed2 > 250 && xdihed3 > 240)
					regionIndexC = 0;
				else if(dihed2 > 250 && xdihed3 > 120)
					regionIndexC = 1;
				else if(dihed2 > 250)
					regionIndexC = 2;
				else if(xdihed3 > 240)
					regionIndexC = 3;
				else if(xdihed3 > 120)
					regionIndexC = 4;
				else
					regionIndexC = 5;

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		return minE*para->wtPho*para->connectRescale;
	}
}

double PO3Builder::getEnergyFast(RiboseConformer* riboConfA, RiboseConformer* riboConfB) {
	XYZ o3,p,op1,op2,o5,c5,c4,o4;

	LocalFrame cs2A = riboConfA->cs2;
	o3 = global2local(cs2A, riboConfA->coords[5]);

	//distance between atom O3' and C5'
	double d0 = cs2A.origin_.distance(riboConfB->coords[6]);

	LocalFrame cs0;
	c5 = global2local(cs2A, riboConfB->coords[6]);
	c4 = global2local(cs2A, riboConfB->coords[3]);
	o4 = global2local(cs2A, riboConfB->coords[4]);
	LocalFrame cs1;
	LocalFrame cs2;


	LocalFrame bestCs1, bestCs2;
	int bestDihed1, bestDihed2;
	int localBestD1, localBestD2;


	double minE = 999999999.9;
	double e, u;

	int dihed1, dihed2;
	int idA, idB;

	double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5;

	int impIndexA = (int)((riboConfA->rot->improper+60)*0.166666667);
	if(impIndexA < 0) impIndexA = 0;
	if(impIndexA > 19) impIndexA = 19;

	char impType = 'A';
	if(riboConfA->rot->improper > 0) impType = 'B';

	int impIndexB = (int)((riboConfB->rot->improper+60)*0.166666667);
	if(impIndexB < 0) impIndexB = 0;
	if(impIndexB > 19) impIndexB = 19;

	int indexDD;

	int d1d2LibSize = this->d1d2Lib1A.size();
	int bestIndex1=0;


	int regionIndexA;
	int regionIndexB;
	int regionIndexC;


	if(d0 > 4.5){
		//If the distance between atom O3' and C5' is larger than 4.5 angstrom, we use default dihedral angles to build PO3
		dihed1 = 100;
		dihed2 = 290;
		indexDD = (dihed1 * 360) + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];

		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;
		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)];

		return e;
	}

	for(int i=0;i<d1d2LibSize;i++){
		if(impType == 'A') {
			dihed1 = (int)d1d2Lib1A[i].x_;
			dihed2 = (int)d1d2Lib1A[i].y_;
		}
		else {
			dihed1 = (int)d1d2Lib1B[i].x_;
			dihed2 = (int)d1d2Lib1B[i].y_;
		}
		indexDD = (dihed1 * 360) + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];

		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;

		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)];

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			bestIndex1 = i;
			minE = e;
		}	
	}


	if(d0 > 3.9){
		//If the distance between atom O3' and C5' is larger than 3.8 angstrom, we use default dihedral angles to build PO3
		return minE*para->wtPho*para->connectRescale;
	}

	for(int i=0;i<d1d2LibSize;i++){
		if(impType == 'A') {
			dihed1 = d1d2Lib2A[bestIndex1][i].x_;
			dihed2 = d1d2Lib2A[bestIndex1][i].y_;
		}
		else {
			dihed1 = d1d2Lib2B[bestIndex1][i].x_;
			dihed2 = d1d2Lib2B[bestIndex1][i].y_;
		}
		indexDD = dihed1 * 360 + dihed2;
		p = this->pList[indexDD];
		o5 = this->o5List[indexDD];
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;
		u = (xd3-len3)*para->rnaKBond;
		if(u < -1)
			e += (-2*u-1);
		else if(u > 1)
			e += 2*u-1;
		else if(u <= 1)
			e += u;
		else 
			e += -u;

		u = (xang3-ang3)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->rnaKAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		if(impIndexA < 10 && dihed2 > 250)
			regionIndexA = 0;
		else if(impIndexA < 10 && dihed1 < 210)
			regionIndexA = 1;
		else if(impIndexA < 10)
			regionIndexA = 2;
		else if(dihed2 > 250)
			regionIndexA = 3;
		else if(dihed1 < 210)
			regionIndexA = 4;
		else
			regionIndexA = 5;


		if(impIndexB < 10 && xdihed5 > 240)
			regionIndexB = 0;
		else if(impIndexB < 10 && xdihed5 > 120)
			regionIndexB = 1;
		else if(impIndexB < 10)
			regionIndexB = 2;
		else if(xdihed5 > 240)
			regionIndexB = 3;
		else if(xdihed5 > 120)
			regionIndexB = 4;
		else
			regionIndexB = 5;

		if(dihed2 > 250 && xdihed3 > 240)
			regionIndexC = 0;
		else if(dihed2 > 250 && xdihed3 > 120)
			regionIndexC = 1;
		else if(dihed2 > 250)
			regionIndexC = 2;
		else if(xdihed3 > 240)
			regionIndexC = 3;
		else if(xdihed3 > 120)
			regionIndexC = 4;
		else
			regionIndexC = 5;

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] + para->rnaDihedImpD1D2Shift[regionIndexA];
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] + para->rnaDihedImpD4D5Shift[regionIndexB];
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] + para->rnaDihedD2D3D4Shift[regionIndexC];

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			minE = e;
		}
	}
	return minE*para->wtPho*para->connectRescale;

}

PO3Builder::~PO3Builder() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
