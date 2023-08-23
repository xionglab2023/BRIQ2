/*
 * KeyToCsMove.cpp
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
#include "predNA/MCRun.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
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

	int bpType = atoi(argv[1]);
	int id = atoi(argv[2]);


	vector<XYZ> localTList;
	ifstream input;
	input.open("/public/home/pengx/cpp/briqx/data/sphere/sphere2000",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;
	out.open("/public/home/pengx/rnaModeling/6dKeys/keyMoves/pair" + string(argv[1]) +"/km-"+string(argv[2]), ios::out);

	AtomLib* atLib = new AtomLib();
	BaseRotamerLib* baseLib = new BaseRotamerLib(atLib);
	int typeA = bpType/4;
	int typeB = bpType%4;
	BaseRotamer* rotA = baseLib->baseLib[typeA];
	BaseRotamer* rotB = baseLib->baseLib[typeB];
	int nA = rotA->atomNum;
	int nB = rotB->atomNum;

	double d, ang;
	int i,j,k,l,a,b;
	double minD,  xd;

	int indexA, indexB;

	char xx[200];

	for(i=0;i<2000;i++) {
		if(i != id) continue;

		XYZ t1 = localTList[i];
		for(j=0;j<2000;j++){
			XYZ t2 = localTList[j];
			for(k=0;k<50;k++){
				d = k*0.3;
				for(l=0;l<45;l++){
					ang = l*8.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);

					minD = 9.9;
					for(a=0;a<nA;a++){
						for(b=0;b<nB;b++){
							xd = local2global(csA, rotA->coordsLocal[a]).distance(local2global(csB, rotB->coordsLocal[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}
					for(a=0;a<nA;a++){
						xd = local2global(csA,  rotA->coordsLocal[a]).distance(csB.origin_);
						if(xd < minD){
							minD = xd;
						}
					}

					for(b=0;b<nB;b++){
						xd = local2global(csB, rotB->coordsLocal[b]).distance(csA.origin_);
						if(xd < minD){
							minD = xd;
						}
					}

					if(minD < 5 && minD > 2.3)
					{
						indexA = k*45+l;
						indexB = i*2000+j;
						CsMove mv = csB - csA;
						sprintf(xx, "%-4d %-7d %s", indexA, indexB, mv.toString().c_str());
						out << string(xx) << endl;
					}
				}
			}
		}
	}
	out.close();
}

