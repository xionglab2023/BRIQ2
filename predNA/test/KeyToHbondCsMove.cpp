/*
 * KeyToHbondCsMove.cpp
 *
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


int main(int argc, char** argv){

	int id = atoi(argv[1]);

	vector<XYZ> localTList;
	ifstream input;
	input.open("/lustre/home/pengx/cpp/RNAModeling/data/sphere500",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;
	out.open("/lustre/home/pengx/rnaModeling/6dKeys/hbondKeyMove/km-"+string(argv[1]), ios::out);

	XYZ atomCoordLocal[4];

	atomCoordLocal[0] = XYZ(0.992 , -1.260 , -0.000); //O3'
	atomCoordLocal[1] = XYZ(0.985 ,  1.251 , -0.000); //O5'
	atomCoordLocal[2] = XYZ(-0.738 , -0.012 ,  1.267); //OP1
	atomCoordLocal[3] = XYZ(-0.750 ,  0.008 , -1.260); //OP2


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
			for(int k=0;k<50;k++){
				ss[4] = k + '!';
				d = k*0.2;
				for(int l=0;l<45;l++){
					ss[5] = l + '!';
					ang = l*8.0;
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


					if(d < 2) continue;

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


