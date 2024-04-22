
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <random>
#include "model/StructureModel.h"
#include "predNA/NuGraph.h"
#include <string>
#include <vector>

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;
using namespace std;


	double clash(double d0, double d, double lamda, double shift){
		double u, e, slop, e1;

		if(d < d0+shift - 0.4){
			u = -0.4*lamda;
			slop = 4*u*u*u*lamda;
			e1 = u*u*u*u;
			e = e1 + slop*(d-d0-shift+0.4);
		}
		else if(d < d0+shift) {
			u = (d-d0-shift)*lamda;
			e = u*u*u*u;
		}
		else
			e = 0;
		return e;
	}

int main(int argc, char** argv){	


	double e;
	double d0 = 2.9;
	double wd = -2.0;
    ForceFieldPara* para = new ForceFieldPara();
    para->hbLamda1 = 0.1;
    para->hbLamda2 = 0.25;

    cout << para->hbLamda1 << endl;
    cout << para->hbLamda2 << endl;

    

    for(double len=2.0;len<5.0;len+=0.02) {
		double e = clash(d0, len, 2.5, 0.0);

	/*
	    if(len < d0)
	    {
		    double u = (len-d0)/d0/para->hbLamda1;
		    e = u*u+wd;
	    }
	    else
	    {
		    double u = (len-d0)/d0/para->hbLamda2;
		    if(u>4)
			    e = 0.0;
		    else
			    e = wd*exp(-u*u);
	    }
	*/
        printf("%4.2f %8.4f\n", len, e);
    }
	


}