/*
 * KeyToCsMove2.cpp
 *
 */

#include "model/BaseRotamer.h"
#include "model/BaseRotamerLib.h"
#include "forcefield/PO3Builder.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

#include "model/StructureModel.h"
#include "predNA/MCRun.h"
#include "forcefield/XPara.h"

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

double calPhoEnergy(double len, double xang3, double xang4, int dihed1, int dihed2, int dihed3, int dihed4, int dihed5, PO3Builder& pb){
	double e = 0;
	double u;

	double len1 = 1.605;
	double len2 = 1.592;
	double len3 = 1.422;
	double ang1 = 120.1;
	double ang2 = 103.5;
	double ang3 = 120.7;
	double ang4 = 111.1;

	u = (len-len3)*3.5;
	e += u*u;

	u = (xang3-ang3)*0.07;
	e += u*u;

	u = (xang4-ang4)*0.07;
	e += u*u;

	e += pb.eDihed1[dihed1];
	e += pb.eDihed2[dihed2];
	e += pb.eDihed3[dihed3];
	e += pb.eDihed4[dihed4];
	e += pb.eDihed5[dihed5];
	return e;
}


int main(int argc, char** argv){

	/*
	int id = atoi(argv[1]);

	vector<XYZ> localTList;
	ifstream input;
	input.open("/export/home/s2982206/aspen/cppWorkplace/RNAModeling/data/sphere300",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	cout << localTList.size() << endl;
	input.close();

	ofstream out;
	out.open("/export/home/s2982206/aspen/rnaModeling/6dKeys/riboKeyMove/km-"+string(argv[1]), ios::out);

	XPara* para = new XPara();
	RnaAtomicEnergyTable* et = new RnaAtomicEnergyTable(para);
	RiboseOxygenEnergyTable* rET = new RiboseOxygenEnergyTable();
	PO3Builder* pb = new PO3Builder(para);

	XYZ t1(-2.004, -1.403, 0);
	XYZ t2(-1.417, 0, 0);
	XYZ t3(0,0,0);
	XYZ t4(-1.997, -1.365, 0);
	XYZ t5(-1.507, 0, 0);
	XYZ t6(0, 0, 0);

	char ss[7];
	double d, ang;
	ss[6] = '\0';
	for(int i=0;i<300;i++) {
		if(i%50 != id) continue;
		ss[0] = i/30 + '!';
		ss[1] = i%30 + '!';
		XYZ t1 = localTList[i];
		for(int j=0;j<300;j++){
			ss[2] = j/30 + '!';
			ss[3] = j%30 + '!';
			XYZ t2 = localTList[j];
			for(int k=0;k<10;k++){
				ss[4] = k + '!';
				d = k*0.2+2.5;
				for(int l=0;l<36;l++){
					ss[5] = l + '!';
					ang = l*10.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);
					string key = string(ss);
					CsMove mv = csB - csA;

					XYZ c2 = local2global(csA, t1);
					XYZ c3 = local2global(csA, t2);
					XYZ o3 = csA.origin_;
					XYZ o4 = local2global(csB, t4);
					XYZ c4 = local2global(csB, t5);
					XYZ c5 = csB.origin_;
					XYZ p, op1, op2, o5;
					int d1d2LibSize = pb->d1d2Lib1A.size();

					int d1Index,d2Index, d1d2Index;
					double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5;
					double u, e;

					double minE = 99999999.9;
					int bestDihed1, bestDihed2, bestIndex1, bestIndex2, localBestD1, localBestD2;
					int dihed1, dihed2;

					for(int i=0;i<d1d2LibSize;i++){
						d1Index = (int)pb->d1d2Lib1A[i].x_;
						d2Index = (int)pb->d1d2Lib1A[i].y_;
						d1d2Index = d1Index*360+d2Index;
						p = local2global(csA, pb->pList[d1d2Index]);
						op1 = local2global(csA,pb->op1List[d1d2Index]);
						op2 = local2global(csA,pb->op2List[d1d2Index]);
						o5 = local2global(csA,pb->o5List[d1d2Index]);
						xd3 = o5.distance(c5);
						xang3 = angleX(p, o5, c5);
						xang4 = angleX(o5, c5, c4);
						xdihed3 = dihedral(o3, p, o5, c5);
						xdihed4 = dihedral(p, o5, c5, c4);
						xdihed5 = dihedral(o5, c5, c4, o4);
						if(xdihed3 < 0) xdihed3 += 360;
						if(xdihed4 < 0) xdihed4 += 360;
						if(xdihed5 < 0) xdihed5 += 360;
						e = calPhoEnergy(xd3, xang3, xang4, d1Index, d2Index, (int)xdihed3, (int)xdihed4, (int)xdihed5, *pb);
						if(e < minE){
							bestDihed1 = d1Index;
							bestDihed2 = d2Index;
							bestIndex1 = i;
							minE = e;
						}
					}

					for(int i=0;i<d1d2LibSize;i++){
						d1Index = (int)pb->d1d2Lib2A[bestIndex1][i].x_;
						d2Index = (int)pb->d1d2Lib2A[bestIndex1][i].y_;
						d1d2Index = d1Index*360+d2Index;
						p = local2global(csA, pb->pList[d1d2Index]);
						op1 = local2global(csA,pb->op1List[d1d2Index]);
						op2 = local2global(csA,pb->op2List[d1d2Index]);
						o5 = local2global(csA,pb->o5List[d1d2Index]);
						xd3 = o5.distance(c5);
						xang3 = angleX(p, o5, c5);
						xang4 = angleX(o5, c5, c4);
						xdihed3 = dihedral(o3, p, o5, c5);
						xdihed4 = dihedral(p, o5, c5, c4);
						xdihed5 = dihedral(o5, c5, c4, o4);
						if(xdihed3 < 0) xdihed3 += 360;
						if(xdihed4 < 0) xdihed4 += 360;
						if(xdihed5 < 0) xdihed5 += 360;
						e = calPhoEnergy(xd3, xang3, xang4, d1Index, d2Index, (int)xdihed3, (int)xdihed4, (int)xdihed5, *pb);
						if(e < minE){
							bestDihed1 = d1Index;
							bestDihed2 = d2Index;
							bestIndex2 = i;
							minE = e;
						}
					}

					{
						localBestD1 = bestDihed1;
						localBestD2 = bestDihed2;
						for(int i=0;i<31;i++){
							dihed1 = localBestD1 + (i-15);
							if(dihed1 < 0) dihed1 = dihed1 +360;
							if(dihed1 > 359) dihed1 = dihed1 - 360;

							for(int j=0;j<31;j++){
								dihed2 = localBestD2 + (j-15);
								if(dihed2 < 0) dihed2 = dihed2 +360;
								if(dihed2 > 359) dihed2 = dihed2 - 360;
								d1d2Index = dihed1*360+dihed2;
								p = local2global(csA, pb->pList[d1d2Index]);
								op1 = local2global(csA,pb->op1List[d1d2Index]);
								op2 = local2global(csA,pb->op2List[d1d2Index]);
								o5 = local2global(csA,pb->o5List[d1d2Index]);
								xd3 = o5.distance(c5);
								xang3 = angleX(p, o5, c5);
								xang4 = angleX(o5, c5, c4);
								xdihed3 = dihedral(o3, p, o5, c5);
								xdihed4 = dihedral(p, o5, c5, c4);
								xdihed5 = dihedral(o5, c5, c4, o4);
								if(xdihed3 < 0) xdihed3 += 360;
								if(xdihed4 < 0) xdihed4 += 360;
								if(xdihed5 < 0) xdihed5 += 360;
								e = calPhoEnergy(xd3, xang3, xang4, d1Index, d2Index, (int)xdihed3, (int)xdihed4, (int)xdihed5, *pb);
								if(e < minE){
									bestDihed1 = dihed1;
									bestDihed2 = dihed2;
									minE = e;
								}
							}
						}
					}

					char xx[200];
					sprintf(xx, "%s %s %3d %3d %7.3f",key.c_str(), mv.toString().c_str(), bestDihed1, bestDihed2, minE);
					out << string(xx) << endl;
				}
			}
		}
	}
	out.close();
	*/
}


