/*
 * ExtendHashkey.cpp
 *
 */


#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include "tools/StringTool.h"
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "geometry/TransMatrix.h"
#include "geometry/Angles.h"
#include "model/BaseDistanceMatrix.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

void extend(const string inputFile, const string outputFile, double dLen, double dAng) {
	clock_t start = clock();

	char tb1[80][60];
	char tb2[80][60];

	int lineIndex = 0;
	string s;
	ifstream f1;
	f1.open("/export/home/s2982206/aspen/rnaModeling/ribose/stat/tb1.txt", ios::in);
	while(getline(f1,s)){
		for(int i=0;i<60;i+=2){
			tb1[lineIndex][i] = s[i];
		}
		lineIndex++;
	}
	f1.close();

	ifstream f2;
	lineIndex=0;
	f2.open("/export/home/s2982206/aspen/rnaModeling/ribose/stat/tb2.txt", ios::in);
	while(getline(f2,s)){
		for(int i=0;i<60;i+=2){
			tb2[lineIndex][i] = s[i];
		}
		lineIndex++;
	}
	f2.close();


	ifstream f;
	f.open(inputFile.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open file: " << inputFile << endl;
		exit(1);
	}
	ofstream of;
	of.open(outputFile, ios::out);
	if(!of.is_open()){
		cout << "fail to open file: " << outputFile << endl;
		exit(1);
	}

	//int size = 6859;

	//char allKey[size][9];


	vector<string> spt;
	int i,j,k;
	double len1, len2, len3, ang1, ang2, ang3;
	double sinx, cosx, d;

	XYZ a1 = XYZ(0 ,  0 , 0);  //O3'
	XYZ b1 = XYZ(-2.008 ,  -1.402  , 0); //C2'
	XYZ c1; //O2'
	XYZ d1 = XYZ(-1.418, 0, 0); //C3'

	XYZ localA2 = XYZ(0, 0, 0); //C5'
	XYZ localB2 = XYZ(-1.508, 0, 0); //C4'
	XYZ localC2 = XYZ(-1.998, -1.365, 0); //O4'

	XYZ a2, b2, c2;

	int x1,x2,x3,x4,x5,x6, x;

	double dm1[9];
	double dm2[9];

	int diffCount[10];
	for(i=0;i<10;i++)
		diffCount[i] = 0;

	string preKey = "";
	string curKey = "";

	int index;
	int n;
	char ss1[10];
	char ss[10];
	ss1[9] = '\0';
	ss[9] = '\0';

	int fivePow[9] = {390625, 78125, 15625, 3125, 625, 125, 25, 5, 1};
	int neibhorKeys[1953125];

	while(getline(f,s)){
		map<string,string> keySet;
		for(i=0;i<1953125;i++)
			neibhorKeys[i] = 0;
		CsMove mv(s);
		vector<string> spt;
		string dim = " ";
		NSPtools::splitString(s, dim, &spt);
		c1 = XYZ(atof(spt[12].c_str()),atof(spt[13].c_str()) , atof(spt[14].c_str()));

		a2 = localA2 * mv.tm + mv.oriMove;
		b2 = localB2 * mv.tm + mv.oriMove;
		c2 = localC2 * mv.tm + mv.oriMove;

		dm1[0] = sqrt((a1.x_-a2.x_)*(a1.x_-a2.x_) + (a1.y_-a2.y_)*(a1.y_-a2.y_) + (a1.z_-a2.z_)*(a1.z_-a2.z_));
		dm1[1] = sqrt((a1.x_-b2.x_)*(a1.x_-b2.x_) + (a1.y_-b2.y_)*(a1.y_-b2.y_) + (a1.z_-b2.z_)*(a1.z_-b2.z_));
		dm1[2] = sqrt((a1.x_-c2.x_)*(a1.x_-c2.x_) + (a1.y_-c2.y_)*(a1.y_-c2.y_) + (a1.z_-c2.z_)*(a1.z_-c2.z_));
		dm1[3] = sqrt((b1.x_-a2.x_)*(b1.x_-a2.x_) + (b1.y_-a2.y_)*(b1.y_-a2.y_) + (b1.z_-a2.z_)*(b1.z_-a2.z_));
		dm1[4] = sqrt((b1.x_-b2.x_)*(b1.x_-b2.x_) + (b1.y_-b2.y_)*(b1.y_-b2.y_) + (b1.z_-b2.z_)*(b1.z_-b2.z_));
		dm1[5] = sqrt((b1.x_-c2.x_)*(b1.x_-c2.x_) + (b1.y_-c2.y_)*(b1.y_-c2.y_) + (b1.z_-c2.z_)*(b1.z_-c2.z_));
		dm1[6] = sqrt((c1.x_-a2.x_)*(c1.x_-a2.x_) + (c1.y_-a2.y_)*(c1.y_-a2.y_) + (c1.z_-a2.z_)*(c1.z_-a2.z_));
		dm1[7] = sqrt((c1.x_-b2.x_)*(c1.x_-b2.x_) + (c1.y_-b2.y_)*(c1.y_-b2.y_) + (c1.z_-b2.z_)*(c1.z_-b2.z_));
		dm1[8] = sqrt((c1.x_-c2.x_)*(c1.x_-c2.x_) + (c1.y_-c2.y_)*(c1.y_-c2.y_) + (c1.z_-c2.z_)*(c1.z_-c2.z_));

		BaseDistanceMatrix dmA(dm1);
		string keyA = dmA.riboConnectHashKey();
		//cout << keyA << endl;
		//cout << a1.toString() << " " << b1.toString() << " " << c1.toString() << endl;
		//cout << a2.toString() << " " << b2.toString() << " " << c2.toString() << endl;

		//cout << dmA.toString() << endl;
		for(int i=0;i<9;i++){
			ss1[i] = keyA[i];
		}

		for(x1=0;x1<41;x1++) {
			ang1 = dAng*(x1-20);
			sinx = sin(ang1*0.0174);
			cosx = cos(ang1*0.0174);
			TransMatrix tm1;
			tm1.mtx[1][1] = cosx;
			tm1.mtx[1][2] = sinx;
			tm1.mtx[2][1] = -sinx;
			tm1.mtx[2][2] = cosx;
			for(x2=0;x2<41;x2++) {
				ang2 = dAng*(x2-20);
				sinx = sin(ang2*0.0174);
				cosx = cos(ang2*0.0174);
				TransMatrix tm2;
				tm2.mtx[0][0] = cosx;
				tm2.mtx[0][2] = -sinx;
				tm2.mtx[2][0] = sinx;
				tm2.mtx[2][2] = cosx;
				for(x3=0;x3<40;x3++) {
					ang3 = dAng*(x3-20);
					sinx = sin(ang3*0.0174);
					cosx = cos(ang3*0.0174);
					TransMatrix tm3;
					tm3.mtx[0][0] = cosx;
					tm3.mtx[0][1] = sinx;
					tm3.mtx[1][0] = -sinx;
					tm3.mtx[1][1] = cosx;


					for(x4=0;x4<15;x4++){
						len1 = dLen*(x4-7);
						for(x5=0;x5<15;x5++){
							len2 = dLen*(x5-7);
							for(x6=0;x6<15;x6++) {
								len3 = dLen*(x6-7);

								XYZ t(len1, len2, len3);
								if(t.length() > 0.5) continue;
								CsMove mv2(t, tm1*tm2*tm3);
								mv2 = mv+mv2;
								a2 = localA2 * mv2.tm + mv2.oriMove;
								b2 = localB2 * mv2.tm + mv2.oriMove;
								c2 = localC2 * mv2.tm + mv2.oriMove;
								dm2[0] = sqrt((a1.x_-a2.x_)*(a1.x_-a2.x_) + (a1.y_-a2.y_)*(a1.y_-a2.y_) + (a1.z_-a2.z_)*(a1.z_-a2.z_));
								dm2[1] = sqrt((a1.x_-b2.x_)*(a1.x_-b2.x_) + (a1.y_-b2.y_)*(a1.y_-b2.y_) + (a1.z_-b2.z_)*(a1.z_-b2.z_));
								dm2[2] = sqrt((a1.x_-c2.x_)*(a1.x_-c2.x_) + (a1.y_-c2.y_)*(a1.y_-c2.y_) + (a1.z_-c2.z_)*(a1.z_-c2.z_));
								dm2[3] = sqrt((b1.x_-a2.x_)*(b1.x_-a2.x_) + (b1.y_-a2.y_)*(b1.y_-a2.y_) + (b1.z_-a2.z_)*(b1.z_-a2.z_));
								dm2[4] = sqrt((b1.x_-b2.x_)*(b1.x_-b2.x_) + (b1.y_-b2.y_)*(b1.y_-b2.y_) + (b1.z_-b2.z_)*(b1.z_-b2.z_));
								dm2[5] = sqrt((b1.x_-c2.x_)*(b1.x_-c2.x_) + (b1.y_-c2.y_)*(b1.y_-c2.y_) + (b1.z_-c2.z_)*(b1.z_-c2.z_));
								dm2[6] = sqrt((c1.x_-a2.x_)*(c1.x_-a2.x_) + (c1.y_-a2.y_)*(c1.y_-a2.y_) + (c1.z_-a2.z_)*(c1.z_-a2.z_));
								dm2[7] = sqrt((c1.x_-b2.x_)*(c1.x_-b2.x_) + (c1.y_-b2.y_)*(c1.y_-b2.y_) + (c1.z_-b2.z_)*(c1.z_-b2.z_));
								dm2[8] = sqrt((c1.x_-c2.x_)*(c1.x_-c2.x_) + (c1.y_-c2.y_)*(c1.y_-c2.y_) + (c1.z_-c2.z_)*(c1.z_-c2.z_));

								if(dm2[0] < 2.487 || dm2[0] > 3.999) continue;
								if(dm2[1] < 2.705 || dm2[1] > 5.157) continue;
								if(dm2[2] < 2.447 || dm2[2] > 6.254) continue;
								if(dm2[3] < 2.701 || dm2[3] > 6.372) continue;
								if(dm2[4] < 3.201 || dm2[4] > 7.403) continue;
								if(dm2[5] < 2.617 || dm2[5] > 8.414) continue;
								if(dm2[6] < 2.324 || dm2[6] > 6.723) continue;
								if(dm2[7] < 2.768 || dm2[7] > 7.767) continue;
								if(dm2[8] < 2.114 || dm2[8] > 8.867) continue;


								BaseDistanceMatrix dmB(dm2);
								string keyB = dmB.riboConnectHashKey();
								index = 0;
								for(int i=0;i<9;i++){
									ss[i] = keyB[i];
									k = ss[i] - ss1[i];
									if(index>=0){
										if(k< -2 || k > 2)
											index = -1;
										else
											index += (k+2)*fivePow[i];
									}
									else if(index == -1) {
										if(k< -3 || k > 3)
											index = -2;
									}
								}


								if(index >= 0 && neibhorKeys[index] == 0){
									double ang1 = angleX(d1, a1, a2);
									double ang2 = angleX(a1, a2, b2);
									if(ang1 > 60 && ang2 > 60) {
										int idA = (int)((dm2[0]-2.4)*50);
										int idB = (int)((ang1-60)*0.5);
										int idC = (int)((ang2-60)*0.5);
										if(tb1[idA][idB] == '1' && tb2[idA][idB] == '1'){
											keySet[keyB] = mv2.toString();
											neibhorKeys[index] = 1;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		clock_t end = clock();

		map<string,string>::iterator it;

		for(it = keySet.begin();it!=keySet.end();++it) {
			of << it->first << " " << it->second << endl;
		}
	}
	of.close();

}

int main(int argc, char** argv){
	extend(string(argv[1]), string(argv[2]), 0.045, 1.2);
}
