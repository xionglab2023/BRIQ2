/*
 * CalculateConnectEnergy.cpp
 *
 *  Created on: 2024,1,9
 *      Author: pengx
 */

#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/PolarGroupLib.h"
#include "forcefield/HbondEnergy.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/BackboneConnectionEnergyCG.h"
#include "forcefield/BasePair6DEnergyTable.h"
#include "model/BaseRotamer.h"

#include "model/RiboseRotamerLib.h"
#include "forcefield/PO3Builder.h"
#include "predNA/NuEnergyCalculator.h"
#include "predNA/NuGraph.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredNA;

double calculateEnergy(BaseConformer* baseConfA, BaseConformer* baseConfB, RiboseConformer* riboseConfA, RiboseConformer* riboseConfB, AtomicClashEnergy* acET, PO3Builder* pbET, HbondEnergy* hbET){

    double e = 0.0;
	int i,j;
	double dd;
	LocalFrame csA, csB;
	int O2UniqueID = 177; //uniqueID: A-O2'

    /*
     * baseA-RiboseB 
     */

	if(squareDistance(baseConfA->coords[0], riboseConfB->coords[2]) < 144.0){

			for(i=0;i<baseConfA->rot->polarAtomNum;i++){
				e += hbET->getEnergy(baseConfA->rot->polarAtomUniqueID[i], baseConfA->csPolar[i], O2UniqueID, riboseConfB->o2Polar);
			}

			for(i=0;i<baseConfA->rot->atomNum;i++){
				for(j=0;j<riboseConfB->rot->atomNum;j++){
					dd = squareDistance(baseConfA->coords[i], riboseConfB->coords[j]);
					if(dd < 16)
						e += acET->getBaseRiboseEnergy(baseConfA->rot->baseType, i, j, dd, 1);
				}
			}
	}    

    

    /*
     * baseB-riboseA
     */
	if(squareDistance(baseConfB->coords[0], riboseConfA->coords[2]) < 144.0){
			
			for(i=0;i<baseConfB->rot->polarAtomNum;i++){
				e += hbET->getEnergy(baseConfB->rot->polarAtomUniqueID[i], baseConfB->csPolar[i], O2UniqueID, riboseConfA->o2Polar);
			}

			for(i=0;i<baseConfB->rot->atomNum;i++){
				for(j=0;j<riboseConfA->rot->atomNum;j++){
					dd = squareDistance(baseConfB->coords[i], riboseConfA->coords[j]);
					if(dd < 16)
						e += acET->getBaseRiboseEnergy(baseConfB->rot->baseType, i, j, dd, -1);
				}
			}
	}  


    /*
     * ribose-ribose
     */

    int nA, nB;
	if(squareDistance(riboseConfA->coords[2], riboseConfB->coords[2]) < 81.0) {

		double dO4O2C2AB = riboseConfA->coords[4].distance((riboseConfB->coords[7] + riboseConfB->coords[1])*0.5);
		double dO4O2C2BA = riboseConfB->coords[4].distance((riboseConfA->coords[7] + riboseConfA->coords[1])*0.5);

		if(dO4O2C2AB < 5.0)
			e += hbET->getO4O2C2Energy(dO4O2C2AB, 1);
		if(dO4O2C2BA < 5.0)
			e += hbET->getO4O2C2Energy(dO4O2C2BA, 1);

		nA = riboseConfA->rot->atomNum;
		nB = riboseConfB->rot->atomNum;
		for(i=0;i<nA;i++){
			for(j=0;j<nB;j++){
				dd = squareDistance(riboseConfA->coords[i], riboseConfB->coords[j]);
				if(dd < 16)
					e += acET->getRiboseRiboseEnergy(i,j,dd, 1);
			}
		}
	}  

    if(e > 5.0) return 99.9;

    double phoEne = pbET->getEnergyFast(riboseConfA, riboseConfB);

    if(phoEne > 10.0) return 99.9;

    return e + phoEne;  
}

void scan(int typeA, int typeB, int index1, const string& outfile){



    int index2;

    char xx[100];
    OrientationIndex x;
    ofstream out;
    out.open(outfile, ios::out);
    int i,j,k,l;
    CsMove cm;
    LocalFrame csA;
    LocalFrame csB;

    cout << "para" << endl;
    ForceFieldPara* para = new ForceFieldPara();//
    cout << "ribo rot lib  " << endl;
    RiboseRotamerLib* riboRotLib = new RiboseRotamerLib(para);//
    cout << "bbcg " << endl;
    BackboneConnectionEnergyCG* bbcg = new BackboneConnectionEnergyCG(para);//
    cout << "acET" << endl;
    AtomicClashEnergy* acET = new AtomicClashEnergy(para);//
    cout << "bpET" << endl;
    NeighborPairEnergyTable* bpET = new NeighborPairEnergyTable(para, typeA, typeB);//
    cout << "atLib " << endl;
   
    PO3Builder* pbET = new PO3Builder(para);//
    
    HbondEnergy* hbET = new HbondEnergy(para);//

    AtomLib* atLib = new AtomLib();//

    BaseRotamer* rotA = new BaseRotamer(typeA, atLib);//
    BaseRotamer* rotB = new BaseRotamer(typeB, atLib);//
    BaseConformer* confA;
    BaseConformer* confB;

    RiboseConformer* riboConfA;
    RiboseConformer* riboConfB;


    double dd;
    double e, totalE;
    int hbNum;

    double minE;
    int minIndexA, minIndexB;
    double eBBCG;

    int rotIndexA, rotIndexB;
    XYZ c1A, o3A, c5B, c1B;
    double d, ang1, ang2, dihed;

    double bpEnergy = 0.0;
	double clashEnergy = 0.0;
    double baseBaseEnergy = 0.0;

	double minDD;
	int nA, nB;

    int spA;
    int spB;

    for(spA=0;spA<2000;spA++){
        cout << spA << endl;
        for(spB=0;spB<2000;spB++){
            index2 = spA*2000+spB;
            cm = x.index2000ToCsMove(index1, index2);
            csB = csA.add(cm);

            if(csA.origin_.distance(csB.origin_) > 10.5) continue;

            confA = new BaseConformer(rotA, csA); //
            confB = new BaseConformer(rotB, csB); //

            clashEnergy = 0.0;
            bpEnergy = 0.0;

	        if(squareDistance(confA->coords[0], confB->coords[0]) < 225.0) {
		        minDD = 999999.9;
		        nA = confA->rot->atomNum;
		        nB = confB->rot->atomNum;
		        for(i=0;i<nA;i++){
			        for(j=0;j<nB;j++){
				        dd = squareDistance(confA->coords[i], confB->coords[j]);
				        if(dd < minDD){
					        minDD = dd;
				        }
				        if(dd < 16.0) {
					        clashEnergy += acET->getBaseBaseEnergy(confA->rot->baseType, i, confB->rot->baseType, j, dd, 1);
				        }
			        }
		        }

		        if(minDD < 20.25){
			        bpEnergy = bpET->getEnergy(confA->cs1, confB->cs1, 1, sqrt(minDD));
		        }
	        }

            baseBaseEnergy = clashEnergy + bpEnergy;

            if(baseBaseEnergy > 4.0) {
                delete confA;
                delete confB;
                continue;
            }

            minE = 999.99;

            for(rotIndexA=0;rotIndexA<20;rotIndexA++){
                for(rotIndexB=0;rotIndexB<20;rotIndexB++){
                    riboConfA = new RiboseConformer(riboRotLib->rotLib20[typeA][rotIndexA], csA); //
                    riboConfB = new RiboseConformer(riboRotLib->rotLib20[typeB][rotIndexB], csB); //

                    e = calculateEnergy(confA, confB, riboConfA, riboConfB, acET, pbET, hbET) + (riboConfA->rot->energy + riboConfB->rot->energy)*0.5;
                    if(e < minE) {
                        minE = e;
                        minIndexA = rotIndexA;
                        minIndexB = rotIndexB;
                    }

                    delete riboConfA;
                    delete riboConfB;
                }
            }

            riboConfB = new RiboseConformer(riboRotLib->rotLib20[typeB][minIndexB], csB); //
            minE = 999.9;
            minIndexA = 0;
            for(rotIndexA=0;rotIndexA<100;rotIndexA++) {
                riboConfA = new RiboseConformer(riboRotLib->rotLib100[typeA][rotIndexA], csA);
                e = calculateEnergy(confA, confB, riboConfA, riboConfB, acET, pbET, hbET)+ (riboConfA->rot->energy + riboConfB->rot->energy)*0.5;
                if(e < minE) {
                    minE = e;
                    minIndexA = rotIndexA;
                }
                delete riboConfA; //
            }
            delete riboConfB;

            riboConfA = new RiboseConformer(riboRotLib->rotLib100[typeA][minIndexA], csA); //
            minE = 999.9;
            minIndexB = 0;
            for(rotIndexB=0;rotIndexB<100;rotIndexB++) {
                riboConfB = new RiboseConformer(riboRotLib->rotLib100[typeB][rotIndexB], csB); //
                e = calculateEnergy(confA, confB, riboConfA, riboConfB, acET, pbET, hbET)+ (riboConfA->rot->energy + riboConfB->rot->energy)*0.5;
                if(e < minE) {
                    minE = e;
                    minIndexB = rotIndexB;
                }
                delete riboConfB;
            }
            delete riboConfA;

            delete confA;
            delete confB;

            if(minE + baseBaseEnergy < 15.0){
                sprintf(xx, "%d %d %8.4f", index1, index2, minE + baseBaseEnergy);
                out << string(xx) << endl;
            }
        }
    }

    out.close();

    delete atLib;
    delete rotA;
    delete rotB;
    delete para;
    delete acET;
    delete bpET;
    delete riboRotLib;
    delete bbcg;
    delete pbET;
    delete hbET;  

}

void single(int typeA, int typeB, int index1, int index2, const string& outfile){
    char xx[100];
    OrientationIndex x;

    int i,j,k,l;
    CsMove cm;
    LocalFrame csA;
    LocalFrame csB;

    cout << "para" << endl;
    ForceFieldPara* para = new ForceFieldPara();
    cout << "ribo rot lib  " << endl;
    RiboseRotamerLib* riboRotLib = new RiboseRotamerLib(para);
    cout << "bbcg " << endl;
    BackboneConnectionEnergyCG* bbcg = new BackboneConnectionEnergyCG(para);
    cout << "acET" << endl;
    AtomicClashEnergy* acET = new AtomicClashEnergy(para);
    cout << "bpET" << endl;
    NeighborPairEnergyTable* bpET = new NeighborPairEnergyTable(para, typeA, typeB);
    cout << "atLib " << endl;
   
    PO3Builder* pbET = new PO3Builder(para);
    
    HbondEnergy* hbET = new HbondEnergy(para);

    AtomLib* atLib = new AtomLib();

    BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
    BaseRotamer* rotB = new BaseRotamer(typeB, atLib);
    BaseConformer* confA;
    BaseConformer* confB;

    RiboseConformer* riboConfA;
    RiboseConformer* riboConfB;


    double dd;
    double e, totalE;
    int hbNum;

    double minE;
    int minIndexA, minIndexB;
    double eBBCG;

    int rotIndexA, rotIndexB;
    XYZ c1A, o3A, c5B, c1B;
    double d, ang1, ang2, dihed;

    double bpEnergy = 0.0;
	double clashEnergy = 0.0;
    double baseBaseEnergy = 0.0;

	double minDD;
	int nA, nB;

    int spA;
    int spB;

   {
        {
            cm = x.index2000ToCsMove(index1, index2);
            csB = csA.add(cm);

           
            confA = new BaseConformer(rotA, csA);
            confB = new BaseConformer(rotB, csB);

            clashEnergy = 0.0;
            bpEnergy = 0.0;

	        if(squareDistance(confA->coords[0], confB->coords[0]) < 225.0) {
		        minDD = 999999.9;
		        nA = confA->rot->atomNum;
		        nB = confB->rot->atomNum;
		        for(i=0;i<nA;i++){
			        for(j=0;j<nB;j++){
				        dd = squareDistance(confA->coords[i], confB->coords[j]);
				        if(dd < minDD){
					        minDD = dd;
				        }
				        if(dd < 16.0) {
					        clashEnergy += acET->getBaseBaseEnergy(confA->rot->baseType, i, confB->rot->baseType, j, dd, 1);
				        }
			        }
		        }

		        if(minDD < 20.25){
			        bpEnergy = bpET->getEnergy(confA->cs1, confB->cs1, 1, sqrt(minDD));
		        }
	        }

            baseBaseEnergy = clashEnergy + bpEnergy;


            minE = 99999.99;

            for(rotIndexA=0;rotIndexA<20;rotIndexA++){
                for(rotIndexB=0;rotIndexB<20;rotIndexB++){
                    riboConfA = new RiboseConformer(riboRotLib->rotLib20[typeA][rotIndexA], csA);
                    riboConfB = new RiboseConformer(riboRotLib->rotLib20[typeB][rotIndexB], csB);

                    e = calculateEnergy(confA, confB, riboConfA, riboConfB, acET, pbET, hbET) + (riboConfA->rot->energy + riboConfB->rot->energy)*0.5;
                    if(e < minE) {
                        minE = e;
                        minIndexA = rotIndexA;
                        minIndexB = rotIndexB;
                    }

                    delete riboConfA;
                    delete riboConfB;
                }
            }

            cout << "rotA: " << minIndexA << " rotB: " << minIndexB << " " << minE << endl;

            riboConfB = new RiboseConformer(riboRotLib->rotLib20[typeB][minIndexB], csB);
            minE = 99999.9;
            minIndexA = 0;
            for(rotIndexA=0;rotIndexA<100;rotIndexA++) {
                riboConfA = new RiboseConformer(riboRotLib->rotLib100[typeA][rotIndexA], csA);
                e = calculateEnergy(confA, confB, riboConfA, riboConfB, acET, pbET, hbET)+ (riboConfA->rot->energy + riboConfB->rot->energy)*0.5;
                if(e < minE) {
                    minE = e;
                    minIndexA = rotIndexA;
                }
                delete riboConfA;
            }
            delete riboConfB;

            cout << "rotA: " << minIndexA << " rotB: " << minIndexB << " " << minE << endl;

            riboConfA = new RiboseConformer(riboRotLib->rotLib100[typeA][minIndexA], csA);
            minE = 99.9;
            minIndexB = 0;
            for(rotIndexB=0;rotIndexB<100;rotIndexB++) {
                riboConfB = new RiboseConformer(riboRotLib->rotLib100[typeB][rotIndexB], csB);
                e = calculateEnergy(confA, confB, riboConfA, riboConfB, acET, pbET, hbET)+ (riboConfA->rot->energy + riboConfB->rot->energy)*0.5;
                if(e < minE) {
                    minE = e;
                    minIndexB = rotIndexB;
                }
                delete riboConfB;
            }
            delete riboConfA;

            cout << "rotA: " << minIndexA << " rotB: " << minIndexB << " " << minE << endl;

            riboConfA = new RiboseConformer(riboRotLib->rotLib100[typeA][minIndexA], csA);
            riboConfB = new RiboseConformer(riboRotLib->rotLib100[typeB][minIndexB], csB);

            cout << "minE " << minE << endl;

            cout << "nodeA" << endl;
            NuNode* nodeA = new NuNode(0, typeA, csA, confA->rot, riboConfA->rot, atLib);
            cout << "nodeB" << endl;
            NuNode* nodeB = new NuNode(1, typeB, csB, confB->rot, riboConfB->rot, atLib); 

            cout << "pho" << endl;
            pbET->buildPhosphate(nodeA->riboseConf, nodeB->riboseConf, nodeA->phoConf);
            cout << "phoEnergy: " << nodeA->phoConf->ene << endl;

            double phoEne1 = pbET->getEnergy(nodeA->riboseConf, nodeB->riboseConf);
             cout << "phoEnergy1: " << phoEne1 << endl;
            double phoEne2 = pbET->getEnergyFast(nodeA->riboseConf, nodeB->riboseConf);
            cout << "phoEnergyFast: " << phoEne2 << endl;

            nodeA->phoConfTmp->copyValueFrom(nodeA->phoConf);

            NuNode* nodes[2];
            nodes[0] = nodeA;
            nodes[1] = nodeB;

            bool connectToDownstream[2];
            connectToDownstream[0] = true;
            connectToDownstream[1] = false;
            int seq[2];
            seq[0] = typeA;
            seq[1] = typeB;

            cout << "graph info" << endl;
            graphInfo* gi = new graphInfo(2, seq, connectToDownstream, nodes, 0.0, atLib);

            cout << "print " << endl;
            gi->printPDB(outfile);



            delete confA;
            delete confB;

            


        }
    }



    delete atLib;
    delete rotA;
    delete rotB;
    delete para;
    delete acET;
    delete bpET;
    delete riboRotLib;
    delete bbcg;
    delete pbET;
    delete hbET;
}

int main(int argc, char** argv){
    string jobType = string(argv[1]);

    if(jobType == "scan") {
        int typeA = atoi(argv[2]);
        int typeB = atoi(argv[3]);
        int index1 = atoi(argv[4]);
        string outfile = string(argv[5]);
        scan(typeA, typeB, index1, outfile);
    }
    else if(jobType == "single") {
        int typeA = atoi(argv[2]);
        int typeB = atoi(argv[3]);
        int index1 = atoi(argv[4]);
        int index2 = atoi(argv[5]);
        string outfile = string(argv[6]);
        single(typeA, typeB, index1, index2, outfile);
    }


}


