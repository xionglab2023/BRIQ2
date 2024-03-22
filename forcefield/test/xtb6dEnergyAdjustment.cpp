#include "forcefield/BasePair6DEnergyTable.h"
#include <sys/time.h>
#include "model/BaseRotamer.h"
#include "model/BaseRotamerLib.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include "model/BasePairLib.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "geometry/xyz.h"
#include "model/StructureModel.h"
#include "predNA/MCRun.h"


using namespace NSPforcefield;
using namespace NSPgeometry;
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

    int typeA = atoi(argv[1]);
    int typeB = atoi(argv[2]);
    int index = atoi(argv[3]);


	cout << "start" << endl;

    char xx[200];
	string augc = "AUGC";

	sprintf(xx, "%c%c", augc[typeA], augc[typeB]);
	string pairType= string(xx);

	sprintf(xx, "%c%c-%d", augc[typeA], augc[typeB], index);
	string outName = string(xx);

    vector<XYZ> originShiftList;
    vector<double> energyShiftList;

    ifstream in;
	cout << "read energy align file" << endl;
    string clusterShiftFile = "/public/home/pengx/briqx/xtb/sp2000/cluster/align/"+pairType+".align";
    in.open(clusterShiftFile.c_str(), ios::in);
    string line;
	double sx, sy, sz, se;
    while(in >> sx >> sy >> sz >> se){
		XYZ os(sx, sy, sz);
		originShiftList.push_back(os);
		energyShiftList.push_back(se);
    }
	in.close();

    int i,j,k,l,a,b,nA, indexA, indexB, nB;
    double d, ang, minD, xd, ene;

	cout << "open output file " << endl;
    ofstream out;
    string outfile = "/public/home/pengx/briqx/xtb/sp2000/ene/adjusted/"+pairType+"/"+outName;
    out.open(outfile.c_str(), ios::out);


	cout << "init para" << endl;
    ForceFieldPara* para = new ForceFieldPara();
	cout << "para xtb" << endl;
    para->bwTag = "xtb";
	para->wtBp1 = 1.0;
	para->wtBp2 = 1.0;

	cout << "init et" << endl;
    BasePair6DEnergyTable* et = new BasePair6DEnergyTable(para, true, 1);
	et->load(para);

	AtomLib* atLib = new AtomLib();
	BaseRotamerLib* baseLib = new BaseRotamerLib(atLib);
	BaseRotamer* rotA = baseLib->baseLib[typeA];
	BaseRotamer* rotB = baseLib->baseLib[typeB];
    nA = rotA->atomNum;
    nB = rotB->atomNum;

    BasePairLib* bpLib = new BasePairLib();

	vector<XYZ> localTList;
	ifstream input;
	input.open("/public/home/pengx/cpp/briqx/data/sphere/sphere2000",ios::in);
	double x,y,z;
	while(input >> x >> y >> z){
		localTList.push_back(XYZ(x,y,z));
	}
	input.close();


	

	
	cout << "start energy calculation: " << endl;
	for(i=0;i<2000;i++) {
		if(i != index) continue;
		
		XYZ t1 = localTList[i];
		for(j=0;j<2000;j++){
		
			XYZ t2 = localTList[j];
			for(k=0;k<50;k++){
				
				d = k*0.3;
				for(l=0;l<45;l++){
				
					ang = l*8.0;
					LocalFrame csA = getCsA(t1, ang, d);
					LocalFrame csB = getCsB(t2, ang, d);

                    BaseDistanceMatrix dm(csA, csB);
                    int clusterID = bpLib->getPairType(dm, typeA, typeB, 2);
                    XYZ originShift = originShiftList[clusterID];
                    double eneShift = energyShiftList[clusterID];

                    csB.origin_ = csB.origin_ - originShift;

					minD = 9.9;
					for(a=0;a<nA;a++){
						for(b=0;b<nB;b++){
							xd = local2global(csA, rotA->coordsLocal[a]).distance(local2global(csB, rotB->coordsLocal[b]));
							if(xd < minD){
								minD = xd;
							}
						}
					}

                    ene = et->getEnergy(csA, csB, typeA, typeB, 2, minD);

					if(indexA == 1620 && indexB == 44943) 
						cout << ene << endl;

					if(ene >= -0.1) continue;

                    ene = 0.5*(ene - eneShift);

                    if(ene < -0.001) {

                        indexA = k*45+l;
					    indexB = i*2000+j;
					    CsMove mv = csB - csA;
					    sprintf(xx, "%-4d %-7d %7.4f %d", indexA, indexB, ene, clusterID);
					    out << string(xx) << endl;
                    }
				}
			}
		}
	}
	
	out.close();

}