/*
 * PO3Builder.cpp
 *
 */

#include <forcefield/PO3Builder.h>

namespace NSPforcefield {

PO3Builder::PO3Builder(XPara* para) {
	this->rotLib = new PhosphateRotamerLib();
	XYZ c2,c3,o3,p,op1,op2,o5,c5,c4,o4;

	double len1 = 1.605;
	double len2 = 1.592;
	double ang1 = 120.1;
	double ang2 = 103.5;

	this->para = para;

	string path = NSPdataio::datapath()+"/";
	ifstream file;
	string s;
	string fileName = path+"dihedEnergy/"+para->dihedEneType + "/impD1D2.ene";
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

	fileName = path+"dihedEnergy/"+para->dihedEneType + "/impD4D5.ene";
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

	fileName = path+"dihedEnergy/"+para->dihedEneType + "/D2D4D3.ene";
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
	XYZ c2,c3,o3,p,op1,op2,o5,c5,c4,o4;

	LocalFrame cs2A = riboConfA->cs2;

	double d0 = cs2A.origin_.distance(riboConfB->coords[7]);

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
		p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
		o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];

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
		u = (xd3-len3)*para->kBond;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang3-ang3)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			bestIndex1 = i;
			minE = e;
		}
	}

	if(d0 > 4.2){
		outPhoConf->updateLocalFrameAndRotamer(cs2A, rotLib->prLib[bestDihed1][bestDihed2], minE);
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
		p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
		o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];
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
		u = (xd3-len3)*para->kBond;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang3-ang3)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

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
				p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
				o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];
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
				u = (xd3-len3)*para->kBond;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang3-ang3)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		outPhoConf->updateLocalFrameAndRotamer(cs2A, rotLib->prLib[bestDihed1][bestDihed2], minE);
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
				p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
				o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];
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
				u = (xd3-len3)*para->kBond;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang3-ang3)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		outPhoConf->updateLocalFrameAndRotamer(cs2A, rotLib->prLib[bestDihed1][bestDihed2], minE);
	}
}

double PO3Builder::getEnergy(RiboseConformer* riboConfA, RiboseConformer* riboConfB){
	XYZ c2,c3,o3,p,op1,op2,o5,c5,c4,o4;

	LocalFrame cs2A = riboConfA->cs2;

	double d0 = cs2A.origin_.distance(riboConfB->coords[7]);

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
		p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
		o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];

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
		u = (xd3-len3)*para->kBond;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang3-ang3)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			bestIndex1 = i;
			minE = e;
		}
	}

	if(d0 > 4.2){
		return minE;
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
		p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
		o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];
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
		u = (xd3-len3)*para->kBond;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang3-ang3)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

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
				p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
				o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];
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
				u = (xd3-len3)*para->kBond;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang3-ang3)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		return minE;
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
				p = this->rotLib->prLib[dihed1][dihed2]->localCoords[0];
				o5 = this->rotLib->prLib[dihed1][dihed2]->localCoords[1];
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
				u = (xd3-len3)*para->kBond;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang3-ang3)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
		indexDD = bestDihed1 * 360 + bestDihed2;
		return minE;
	}
}

PO3Builder::~PO3Builder() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
