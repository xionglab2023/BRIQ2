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

using namespace std;
using namespace NSPgeometry;

void extend(const string inputFile, const string outputFile, double dLen, double dAng) {
	clock_t start = clock();
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

	string s;
	vector<string> spt;
	int i,j,k;
	double len1, len2, len3, ang1, ang2, ang3;
	double sinx, cosx, d;

	XYZ a1 = XYZ(5.722 ,  4.303 , -1.758);
	XYZ b1 = XYZ(-2.635 ,  2.489  , 0.351);
	XYZ c1 = XYZ(6.605 , -1.044  , 3.374);
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
		set<string> keySet;
		for(i=0;i<1953125;i++)
			neibhorKeys[i] = 0;
		CsMove mv(s);
		a2 = a1 * mv.tm + mv.oriMove;
		b2 = b1 * mv.tm + mv.oriMove;
		c2 = c1 * mv.tm + mv.oriMove;

		dm1[0] = (a1.x_-a2.x_)*(a1.x_-a2.x_) + (a1.y_-a2.y_)*(a1.y_-a2.y_) + (a1.z_-a2.z_)*(a1.z_-a2.z_);
		dm1[1] = (a1.x_-b2.x_)*(a1.x_-b2.x_) + (a1.y_-b2.y_)*(a1.y_-b2.y_) + (a1.z_-b2.z_)*(a1.z_-b2.z_);
		dm1[2] = ((a1.x_-c2.x_)*(a1.x_-c2.x_) + (a1.y_-c2.y_)*(a1.y_-c2.y_) + (a1.z_-c2.z_)*(a1.z_-c2.z_));
		dm1[3] = ((b1.x_-a2.x_)*(b1.x_-a2.x_) + (b1.y_-a2.y_)*(b1.y_-a2.y_) + (b1.z_-a2.z_)*(b1.z_-a2.z_));
		dm1[4] = ((b1.x_-b2.x_)*(b1.x_-b2.x_) + (b1.y_-b2.y_)*(b1.y_-b2.y_) + (b1.z_-b2.z_)*(b1.z_-b2.z_));
		dm1[5] = ((b1.x_-c2.x_)*(b1.x_-c2.x_) + (b1.y_-c2.y_)*(b1.y_-c2.y_) + (b1.z_-c2.z_)*(b1.z_-c2.z_));
		dm1[6] = ((c1.x_-a2.x_)*(c1.x_-a2.x_) + (c1.y_-a2.y_)*(c1.y_-a2.y_) + (c1.z_-a2.z_)*(c1.z_-a2.z_));
		dm1[7] = ((c1.x_-b2.x_)*(c1.x_-b2.x_) + (c1.y_-b2.y_)*(c1.y_-b2.y_) + (c1.z_-b2.z_)*(c1.z_-b2.z_));
		dm1[8] = ((c1.x_-c2.x_)*(c1.x_-c2.x_) + (c1.y_-c2.y_)*(c1.y_-c2.y_) + (c1.z_-c2.z_)*(c1.z_-c2.z_));

		{
			// get init key
			for(i=0;i<9;i++){
				d = dm1[i];
				if(d<368.64) {
					if(d<92.16) {
						if(d < 23.04) {
							if(d < 5.76) {
								if(d < 1.44)
									ss1[i] = 'A';
								else
									ss1[i] = 'B';
							}
							else {
								if(d < 12.96)
									ss1[i] = 'C';
								else
									ss1[i] = 'D';
							}
						}
						else {
							if(d < 51.84) {
								if(d < 36)
									ss1[i] = 'E';
								else
									ss1[i] = 'F';
							}
							else {
								if(d < 70.56)
									ss1[i] = 'G';
								else
									ss1[i] = 'H';
							}
						}
					}
					else{
						if(d < 207.36) {
							if(d < 144) {
								if(d < 116.64)
									ss1[i] = 'I';
								else
									ss1[i] = 'J';
							}
							else {
								if(d < 174.24)
									ss1[i] = 'K';
								else
									ss1[i] = 'L';
							}
						}
						else {
							if(d < 282.24) {
								if(d < 243.36)
									ss1[i] = 'M';
								else
									ss1[i] = 'N';
							}
							else {
								if(d < 324)
									ss1[i] = 'O';
								else
									ss1[i] = 'P';
							}
						}
					}
				}
				else {
					if(d<829.44) {
						if(d < 576) {
							if(d < 466.56) {
								if(d < 416.16)
									ss1[i] = 'Q';
								else
									ss1[i] = 'R';
							}
							else {
								if(d < 519.84)
									ss1[i] = 'S';
								else
									ss1[i] = 'T';
							}
						}
						else {
							if(d < 696.96) {
								if(d < 635.04)
									ss1[i] = 'U';
								else
									ss1[i] = 'V';
							}
							else {
								if(d < 761.76)
									ss1[i] = 'W';
								else
									ss1[i] = 'X';
							}
						}
					}
					else if(d < 900)
						ss1[i] = 'Y';
					else
						ss1[i] = 'Z';
				}
			}
		}

		for(x1=0;x1<29;x1++) {
			ang1 = dAng*(x1-14);
			sinx = sin(ang1*0.0174);
			cosx = cos(ang1*0.0174);
			TransMatrix tm1;
			tm1.mtx[1][1] = cosx;
			tm1.mtx[1][2] = sinx;
			tm1.mtx[2][1] = -sinx;
			tm1.mtx[2][2] = cosx;
			for(x2=0;x2<29;x2++) {
				ang2 = dAng*(x2-14);
				sinx = sin(ang2*0.0174);
				cosx = cos(ang2*0.0174);
				TransMatrix tm2;
				tm2.mtx[0][0] = cosx;
				tm2.mtx[0][2] = -sinx;
				tm2.mtx[2][0] = sinx;
				tm2.mtx[2][2] = cosx;
				for(x3=0;x3<29;x3++) {
					ang3 = dAng*(x3-14);
					sinx = sin(ang3*0.0174);
					cosx = cos(ang3*0.0174);
					TransMatrix tm3;
					tm3.mtx[0][0] = cosx;
					tm3.mtx[0][1] = sinx;
					tm3.mtx[1][0] = -sinx;
					tm3.mtx[1][1] = cosx;


					for(x4=0;x4<29;x4++){
						len1 = dLen*(x4-14);
						for(x5=0;x5<29;x5++){
							len2 = dLen*(x5-14);
							for(x6=0;x6<29;x6++) {
								len3 = dLen*(x6-14);
								if(len1*len1+len2*len2+len3*len3 > 8.0) continue;
								//n++;
								//cout << n << endl;
								XYZ t(len1, len2, len3);
								CsMove mv2(t, tm1*tm2*tm3);
								mv2 = mv+mv2;
								a2 = a1 * mv2.tm + mv2.oriMove;
								b2 = b1 * mv2.tm + mv2.oriMove;
								c2 = c1 * mv2.tm + mv2.oriMove;
								dm2[0] = ((a1.x_-a2.x_)*(a1.x_-a2.x_) + (a1.y_-a2.y_)*(a1.y_-a2.y_) + (a1.z_-a2.z_)*(a1.z_-a2.z_));
								dm2[1] = ((a1.x_-b2.x_)*(a1.x_-b2.x_) + (a1.y_-b2.y_)*(a1.y_-b2.y_) + (a1.z_-b2.z_)*(a1.z_-b2.z_));
								dm2[2] = ((a1.x_-c2.x_)*(a1.x_-c2.x_) + (a1.y_-c2.y_)*(a1.y_-c2.y_) + (a1.z_-c2.z_)*(a1.z_-c2.z_));
								dm2[3] = ((b1.x_-a2.x_)*(b1.x_-a2.x_) + (b1.y_-a2.y_)*(b1.y_-a2.y_) + (b1.z_-a2.z_)*(b1.z_-a2.z_));
								dm2[4] = ((b1.x_-b2.x_)*(b1.x_-b2.x_) + (b1.y_-b2.y_)*(b1.y_-b2.y_) + (b1.z_-b2.z_)*(b1.z_-b2.z_));
								dm2[5] = ((b1.x_-c2.x_)*(b1.x_-c2.x_) + (b1.y_-c2.y_)*(b1.y_-c2.y_) + (b1.z_-c2.z_)*(b1.z_-c2.z_));
								dm2[6] = ((c1.x_-a2.x_)*(c1.x_-a2.x_) + (c1.y_-a2.y_)*(c1.y_-a2.y_) + (c1.z_-a2.z_)*(c1.z_-a2.z_));
								dm2[7] = ((c1.x_-b2.x_)*(c1.x_-b2.x_) + (c1.y_-b2.y_)*(c1.y_-b2.y_) + (c1.z_-b2.z_)*(c1.z_-b2.z_));
								dm2[8] = ((c1.x_-c2.x_)*(c1.x_-c2.x_) + (c1.y_-c2.y_)*(c1.y_-c2.y_) + (c1.z_-c2.z_)*(c1.z_-c2.z_));

								index = 0;
								for(i=0;i<9;i++){
									d = dm2[i];
									if(d<368.64) {
										if(d<92.16) {
											if(d < 23.04) {
												if(d < 5.76) {
													if(d < 1.44)
														ss[i] = 'A';
													else
														ss[i] = 'B';
												}
												else {
													if(d < 12.96)
														ss[i] = 'C';
													else
														ss[i] = 'D';
												}
											}
											else {
												if(d < 51.84) {
													if(d < 36)
														ss[i] = 'E';
													else
														ss[i] = 'F';
												}
												else {
													if(d < 70.56)
														ss[i] = 'G';
													else
														ss[i] = 'H';
												}
											}
										}
										else{
											if(d < 207.36) {
												if(d < 144) {
													if(d < 116.64)
														ss[i] = 'I';
													else
														ss[i] = 'J';
												}
												else {
													if(d < 174.24)
														ss[i] = 'K';
													else
														ss[i] = 'L';
												}
											}
											else {
												if(d < 282.24) {
													if(d < 243.36)
														ss[i] = 'M';
													else
														ss[i] = 'N';
												}
												else {
													if(d < 324)
														ss[i] = 'O';
													else
														ss[i] = 'P';
												}
											}
										}
									}
									else {
										if(d<829.44) {
											if(d < 576) {
												if(d < 466.56) {
													if(d < 416.16)
														ss[i] = 'Q';
													else
														ss[i] = 'R';
												}
												else {
													if(d < 519.84)
														ss[i] = 'S';
													else
														ss[i] = 'T';
												}
											}
											else {
												if(d < 696.96) {
													if(d < 635.04)
														ss[i] = 'U';
													else
														ss[i] = 'V';
												}
												else {
													if(d < 761.76)
														ss[i] = 'W';
													else
														ss[i] = 'X';
												}
											}
										}
										else if(d < 900)
											ss[i] = 'Y';
										else
											ss[i] = 'Z';
									}

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
									keySet.insert(string(ss));
									neibhorKeys[index] = 1;
								}

							}
						}
					}
				}
			}
		}

		clock_t end = clock();
		for(string key : keySet) {
			of << key << endl;
		}
	}
	of.close();

}

int main(int argc, char** argv){

	extend(string(argv[1]), string(argv[2]), 0.18, 2.2);


}
