/*
 * KeyToHbondCsMove.cpp
 *
 */


/*
 * KeyToCsMove.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: s2982206
 */



#include "model/BaseRotamer.h"
#include "model/BaseRotamerLib.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "model/StructureModel.h"

using namespace NSPmodel;
using namespace std;

LocalFrame getCsA(XYZ t, double dihed, double dist){
	double x = t.x_;
	double y = t.y_;
	double z;
	z = sqrt(1-x*x-y*y);
	if(t.z_ < 0) z = -z;
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);

	double sinD1, cosD1, sinD2, cosD2;

	sinD1 = z/sx/sy;
	cosD1 = -x*y/sx/sy;


	double ang1 = atan2(sinD1, cosD1);

	double ang0 = dihed*0.008726646;
	double tm[3][3];
	tm[0][0] = x;
	tm[0][1] = y;
	tm[0][2] = z;
	tm[1][0] = sx*sin(ang0);
	tm[1][1] = sy*sin(ang0-ang1);
	tm[2][0] = sx*cos(ang0);
	tm[2][1] = sy*cos(ang0-ang1);

	XYZ t1(x, tm[1][0], tm[2][0]);
	XYZ t2(y, tm[1][1], tm[2][1]);
	XYZ t3 = t1^t2;
	tm[1][2] = t3.y_;
	tm[2][2] = t3.z_;

	XYZ ori(-dist*0.5, 0, 0);
	LocalFrame cs(ori, tm);
	return cs;
}

LocalFrame getCsB(XYZ t, double dihed, double dist){
	double x = t.x_;
	double y = t.y_;
	double z;
	z = sqrt(1-x*x-y*y);
	if(t.z_ < 0) z = -z;
	double sx = sqrt(1-x*x);
	double sy = sqrt(1-y*y);
	double sz = sqrt(1-z*z);

	double sinD1, cosD1, sinD2, cosD2;

	sinD1 = z/sx/sy;
	cosD1 = -x*y/sx/sy;


	double ang1 = atan2(sinD1, cosD1);

	double ang0 = -dihed*0.008726646;
	double tm[3][3];
	tm[0][0] = x;
	tm[0][1] = y;
	tm[0][2] = z;
	tm[1][0] = sx*sin(ang0);
	tm[1][1] = sy*sin(ang0-ang1);
	tm[2][0] = sx*cos(ang0);
	tm[2][1] = sy*cos(ang0-ang1);

	XYZ t1(x, tm[1][0], tm[2][0]);
	XYZ t2(y, tm[1][1], tm[2][1]);
	XYZ t3 = t1^t2;
	tm[1][2] = t3.y_;
	tm[2][2] = t3.z_;

	XYZ ori(dist*0.5, 0, 0);
	LocalFrame cs(ori, tm);
	return cs;
}


void getPolarMove1(int id){

	vector<XYZ> localTList;
	ifstream input;
	input.open("/user/xiongpeng/cpp/ProteinModeling/data/sphere/sphere500",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;

	char xx[20];
	sprintf(xx, "%d", id);

	out.open("/user/xiongpeng/proteinModeling/6dKeys/keyMove1/km-"+string(xx), ios::out);

	XYZ atomCoordLocal[4];
	atomCoordLocal[0] = XYZ(0.000 ,  0.000 ,  0.000); //CG
	atomCoordLocal[1] = XYZ(0.670 , -1.150 ,  0.000); //OD1
	atomCoordLocal[2] = XYZ(0.670 ,  1.150 ,  0.000); //OD2
	atomCoordLocal[3] = XYZ(-1.331 , 0.000 ,  0.000); //CB


	char ss[7];
	double d, ang;
	ss[6] = '\0';
	for(int i=0;i<500;i++) {
		if(i%50 != id) continue;
		ss[0] = i/30 + '!';
		ss[1] = i%30 + '!';
		XYZ t1 = localTList[i];
		for(int j=0;j<500;j++){
			ss[2] = j/30 + '!';
			ss[3] = j%30 + '!';
			XYZ t2 = localTList[j];
			for(int k=0;k<80;k++){
				ss[4] = k + '!';
				d = k*0.2;
				for(int l=0;l<60;l++){
					ss[5] = l + '!';
					ang = l*6.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);


					double minD = 9.9;
					for(int a=0;a<4;a++){
						for(int b=0;b<4;b++){
							double xd = local2global(csA, atomCoordLocal[a]).distance(local2global(csB, atomCoordLocal[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}

					if(minD < 4.5)
					{
						string key = string(ss);
						CsMove mv = csB - csA;
						out << key << " " << mv.toString() << endl;
					}
				}
			}
		}
	}
	out.close();
}

void getPolarMove2(int id){

	vector<XYZ> localTList;
	ifstream input;
	input.open("/user/xiongpeng/cpp/ProteinModeling/data/sphere/sphere500",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;

	char xx[20];
	sprintf(xx, "%d", id);

	out.open("/user/xiongpeng/proteinModeling/6dKeys/keyMove2/km-"+string(xx), ios::out);



	//TYR-TYR
	XYZ atomCoordLocal[7];
	atomCoordLocal[0] = XYZ(0.000 ,  0.000 ,  0.000);
	atomCoordLocal[1] = XYZ(2.109 , -1.205 ,  0.000);
	atomCoordLocal[2] = XYZ(2.104 ,  1.202 ,  0.000);
	atomCoordLocal[3] = XYZ(2.794 ,  0.000 ,  0.000);
	atomCoordLocal[4] = XYZ(4.168 ,  0.006 ,  0.000);
	atomCoordLocal[5] = XYZ(0.719 , -1.198 ,  0.000);
	atomCoordLocal[6] = XYZ(0.714 ,  1.196 ,  0.000);



	char ss[7];
	double d, ang;
	ss[6] = '\0';
	for(int i=0;i<500;i++) {
		if(i%50 != id) continue;
		ss[0] = i/30 + '!';
		ss[1] = i%30 + '!';
		XYZ t1 = localTList[i];
		for(int j=0;j<500;j++){
			ss[2] = j/30 + '!';
			ss[3] = j%30 + '!';
			XYZ t2 = localTList[j];
			for(int k=0;k<80;k++){
				ss[4] = k + '!';
				d = k*0.2;
				for(int l=0;l<60;l++){
					ss[5] = l + '!';
					ang = l*6.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);


					double minD = 9.9;
					for(int a=0;a<7;a++){
						for(int b=0;b<7;b++){
							double xd = local2global(csA, atomCoordLocal[a]).distance(local2global(csB, atomCoordLocal[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}

					if(minD < 4.5)
					{
						string key = string(ss);
						CsMove mv = csB - csA;
						out << key << " " << mv.toString() << endl;
					}
				}
			}
		}
	}
	out.close();
}

void getPolarMove3(int id){

	vector<XYZ> localTList;
	ifstream input;
	input.open("/user/xiongpeng/cpp/ProteinModeling/data/sphere/sphere500",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;

	char xx[20];
	sprintf(xx, "%d", id);

	out.open("/user/xiongpeng/proteinModeling/6dKeys/keyMove3/km-"+string(xx), ios::out);



	//TYR-TYR
	XYZ atomCoordLocal[9];
	atomCoordLocal[0] = XYZ(0.000 ,  0.000 ,  0.000);
	atomCoordLocal[1] = XYZ(1.940 , -1.153 ,  0.000);
	atomCoordLocal[2] = XYZ(3.608 ,  2.145 ,  0.000);
	atomCoordLocal[3] = XYZ(0.566 , -1.251 ,  0.000);
	atomCoordLocal[4] = XYZ(1.088 ,  0.930 ,  0.000);
	atomCoordLocal[5] = XYZ(2.288 ,  0.174 ,  0.000);
	atomCoordLocal[6] = XYZ(1.168 ,  2.329 ,  0.000);
	atomCoordLocal[7] = XYZ(3.555 ,  0.774 ,  0.000);
	atomCoordLocal[8] = XYZ(2.435 ,  2.925 ,  0.000);

	char ss[7];
	double d, ang;
	ss[6] = '\0';
	for(int i=0;i<500;i++) {
		if(i%50 != id) continue;
		ss[0] = i/30 + '!';
		ss[1] = i%30 + '!';
		XYZ t1 = localTList[i];
		for(int j=0;j<500;j++){
			ss[2] = j/30 + '!';
			ss[3] = j%30 + '!';
			XYZ t2 = localTList[j];
			for(int k=0;k<80;k++){
				ss[4] = k + '!';
				d = k*0.2;
				for(int l=0;l<60;l++){
					ss[5] = l + '!';
					ang = l*6.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);


					double minD = 9.9;
					for(int a=0;a<9;a++){
						for(int b=0;b<9;b++){
							double xd = local2global(csA, atomCoordLocal[a]).distance(local2global(csB, atomCoordLocal[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}

					if(minD < 4.5)
					{
						string key = string(ss);
						CsMove mv = csB - csA;
						out << key << " " << mv.toString() << endl;
					}
				}
			}
		}
	}
	out.close();
}


void getPolarMove4(int id){

	vector<XYZ> localTList;
	ifstream input;
	input.open("/user/xiongpeng/cpp/ProteinModeling/data/sphere/sphere500",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;

	char xx[20];
	sprintf(xx, "%d", id);

	out.open("/user/xiongpeng/proteinModeling/6dKeys/keyMove4/km-"+string(xx), ios::out);



	XYZ atomCoordLocal1[4];
	atomCoordLocal1[0] = XYZ(0.000 ,  0.000 ,  0.000); //CG
	atomCoordLocal1[1] = XYZ(0.670 , -1.150 ,  0.000); //OD1
	atomCoordLocal1[2] = XYZ(0.670 ,  1.150 ,  0.000); //OD2
	atomCoordLocal1[3] = XYZ(-1.331 , 0.000 ,  0.000); //CB

	XYZ atomCoordLocal2[7];
	atomCoordLocal2[0] = XYZ(0.000 ,  0.000 ,  0.000);
	atomCoordLocal2[1] = XYZ(2.109 , -1.205 ,  0.000);
	atomCoordLocal2[2] = XYZ(2.104 ,  1.202 ,  0.000);
	atomCoordLocal2[3] = XYZ(2.794 ,  0.000 ,  0.000);
	atomCoordLocal2[4] = XYZ(4.168 ,  0.006 ,  0.000);
	atomCoordLocal2[5] = XYZ(0.719 , -1.198 ,  0.000);
	atomCoordLocal2[6] = XYZ(0.714 ,  1.196 ,  0.000);



	char ss[7];
	double d, ang;
	ss[6] = '\0';
	for(int i=0;i<500;i++) {
		if(i%50 != id) continue;
		ss[0] = i/30 + '!';
		ss[1] = i%30 + '!';
		XYZ t1 = localTList[i];
		for(int j=0;j<500;j++){
			ss[2] = j/30 + '!';
			ss[3] = j%30 + '!';
			XYZ t2 = localTList[j];
			for(int k=0;k<80;k++){
				ss[4] = k + '!';
				d = k*0.2;
				for(int l=0;l<60;l++){
					ss[5] = l + '!';
					ang = l*6.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);


					double minD = 9.9;
					for(int a=0;a<4;a++){
						for(int b=0;b<7;b++){
							double xd = local2global(csA, atomCoordLocal1[a]).distance(local2global(csB, atomCoordLocal2[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}

					if(minD < 4.5)
					{
						string key = string(ss);
						CsMove mv = csB - csA;
						out << key << " " << mv.toString() << endl;
					}
				}
			}
		}
	}
	out.close();
}

void getPolarMove5(int id){

	vector<XYZ> localTList;
	ifstream input;
	input.open("/user/xiongpeng/cpp/ProteinModeling/data/sphere/sphere500",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;

	char xx[20];
	sprintf(xx, "%d", id);

	out.open("/user/xiongpeng/proteinModeling/6dKeys/keyMove5/km-"+string(xx), ios::out);



	XYZ atomCoordLocal1[4];
	atomCoordLocal1[0] = XYZ(0.000 ,  0.000 ,  0.000); //CG
	atomCoordLocal1[1] = XYZ(0.670 , -1.150 ,  0.000); //OD1
	atomCoordLocal1[2] = XYZ(0.670 ,  1.150 ,  0.000); //OD2
	atomCoordLocal1[3] = XYZ(-1.331 , 0.000 ,  0.000); //CB

	XYZ atomCoordLocal2[9];
	atomCoordLocal2[0] = XYZ(0.000 ,  0.000 ,  0.000);
	atomCoordLocal2[1] = XYZ(1.940 , -1.153 ,  0.000);
	atomCoordLocal2[2] = XYZ(3.608 ,  2.145 ,  0.000);
	atomCoordLocal2[3] = XYZ(0.566 , -1.251 ,  0.000);
	atomCoordLocal2[4] = XYZ(1.088 ,  0.930 ,  0.000);
	atomCoordLocal2[5] = XYZ(2.288 ,  0.174 ,  0.000);
	atomCoordLocal2[6] = XYZ(1.168 ,  2.329 ,  0.000);
	atomCoordLocal2[7] = XYZ(3.555 ,  0.774 ,  0.000);
	atomCoordLocal2[8] = XYZ(2.435 ,  2.925 ,  0.000);



	char ss[7];
	double d, ang;
	ss[6] = '\0';
	for(int i=0;i<500;i++) {
		if(i%50 != id) continue;
		ss[0] = i/30 + '!';
		ss[1] = i%30 + '!';
		XYZ t1 = localTList[i];
		for(int j=0;j<500;j++){
			ss[2] = j/30 + '!';
			ss[3] = j%30 + '!';
			XYZ t2 = localTList[j];
			for(int k=0;k<80;k++){
				ss[4] = k + '!';
				d = k*0.2;
				for(int l=0;l<60;l++){
					ss[5] = l + '!';
					ang = l*6.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);


					double minD = 9.9;
					for(int a=0;a<4;a++){
						for(int b=0;b<7;b++){
							double xd = local2global(csA, atomCoordLocal1[a]).distance(local2global(csB, atomCoordLocal2[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}

					if(minD < 4.5)
					{
						string key = string(ss);
						CsMove mv = csB - csA;
						out << key << " " << mv.toString() << endl;
					}
				}
			}
		}
	}
	out.close();
}

void getPolarMove6(int id){

	vector<XYZ> localTList;
	ifstream input;
	input.open("/user/xiongpeng/cpp/ProteinModeling/data/sphere/sphere500",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;

	char xx[20];
	sprintf(xx, "%d", id);

	out.open("/user/xiongpeng/proteinModeling/6dKeys/keyMove6/km-"+string(xx), ios::out);



	XYZ atomCoordLocal1[7];
	atomCoordLocal1[0] = XYZ(0.000 ,  0.000 ,  0.000);
	atomCoordLocal1[1] = XYZ(2.109 , -1.205 ,  0.000);
	atomCoordLocal1[2] = XYZ(2.104 ,  1.202 ,  0.000);
	atomCoordLocal1[3] = XYZ(2.794 ,  0.000 ,  0.000);
	atomCoordLocal1[4] = XYZ(4.168 ,  0.006 ,  0.000);
	atomCoordLocal1[5] = XYZ(0.719 , -1.198 ,  0.000);
	atomCoordLocal1[6] = XYZ(0.714 ,  1.196 ,  0.000);

	XYZ atomCoordLocal2[9];
	atomCoordLocal2[0] = XYZ(0.000 ,  0.000 ,  0.000);
	atomCoordLocal2[1] = XYZ(1.940 , -1.153 ,  0.000);
	atomCoordLocal2[2] = XYZ(3.608 ,  2.145 ,  0.000);
	atomCoordLocal2[3] = XYZ(0.566 , -1.251 ,  0.000);
	atomCoordLocal2[4] = XYZ(1.088 ,  0.930 ,  0.000);
	atomCoordLocal2[5] = XYZ(2.288 ,  0.174 ,  0.000);
	atomCoordLocal2[6] = XYZ(1.168 ,  2.329 ,  0.000);
	atomCoordLocal2[7] = XYZ(3.555 ,  0.774 ,  0.000);
	atomCoordLocal2[8] = XYZ(2.435 ,  2.925 ,  0.000);



	char ss[7];
	double d, ang;
	ss[6] = '\0';
	for(int i=0;i<500;i++) {
		if(i%50 != id) continue;
		ss[0] = i/30 + '!';
		ss[1] = i%30 + '!';
		XYZ t1 = localTList[i];
		for(int j=0;j<500;j++){
			ss[2] = j/30 + '!';
			ss[3] = j%30 + '!';
			XYZ t2 = localTList[j];
			for(int k=0;k<80;k++){
				ss[4] = k + '!';
				d = k*0.2;
				for(int l=0;l<60;l++){
					ss[5] = l + '!';
					ang = l*6.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);


					double minD = 9.9;
					for(int a=0;a<4;a++){
						for(int b=0;b<7;b++){
							double xd = local2global(csA, atomCoordLocal1[a]).distance(local2global(csB, atomCoordLocal2[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}

					if(minD < 4.5)
					{
						string key = string(ss);
						CsMove mv = csB - csA;
						out << key << " " << mv.toString() << endl;
					}
				}
			}
		}
	}
	out.close();
}

int main(int argc, char** argv){

	int id = atoi(argv[2]);
	string type = string(argv[1]);

	cout << "type: " << type << endl;
	cout << "id: " << id << endl;

	if(type == "km1")
		getPolarMove1(id);
	else if(type == "km2")
		getPolarMove2(id);
	else if(type == "km3")
		getPolarMove3(id);
	else if(type == "km4")
		getPolarMove4(id);
	else if(type == "km5")
		getPolarMove5(id);
	else if(type == "km6")
		getPolarMove6(id);

}



