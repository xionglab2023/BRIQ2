

#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/PolarGroupLib.h"
#include "model/BaseRotamerLib.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

    /*
     *  usage: xtbEnergy $GROUP_TYPE_A $GROUP_TYPE_B $OI_TYPE $OI_INDEX ($OI_INDEX2) 
     */

    string augc = "AUGC";

    int typeA = atoi(argv[1]);
    int typeB = atoi(argv[2]);
    
    AtomLib* atLib = new AtomLib();
    BaseRotamerLib baseLib(atLib);
    OrientationIndex x;
    double minD, xd;
    int i,j,k,l,a,b;
    int index, id2;
    CsMove cm;
    LocalFrame csA, csB;
    int nA, nB;
    {
        BaseRotamer* rotA = baseLib.baseLib[typeA];
        nA = rotA->atomNum;
        {
            BaseRotamer* rotB = baseLib.baseLib[typeB];
            nB = rotB->atomNum;
            for(i=1;i<50;i++){
                for(j=0;j<40;j++){
                    for(k=0;k<500;k++){
                        for(l=0;l<500;l++){
                            index = i*10000000+j*250000+k*500+l;
                     
                            cm = x.index500ToCsMove(index);
                            id2 = x.moveToIndex500(cm);
                            if(index != id2) {
                                cout << "error   " << endl;
                            }

                            csA = LocalFrame();
					        csB = csA.add(cm);

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

                            if(minD < 2.0 || minD > 5.0) 
                                continue;

                            cout << augc[typeA] << " " << augc[typeB] << " " << index << endl;
                        }
                    }
                }
            }
        }

    }

}


