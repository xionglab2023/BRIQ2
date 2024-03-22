#include <vector>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "geometry/OrientationIndex.h"
#include "model/PolarGroupLib.h"
#include "forcefield/HbondEnergy.h"
#include "forcefield/AtomicClashEnergy.h"
#include "forcefield/Xtb6dEnergy.h"
#include "forcefield/BasePair6DEnergyTable.h"
#include "model/BaseRotamer.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPforcefield;


int main(int argc, char** argv){

    
    string nbnnb = "nnb";

    cout << "start" << endl;
    vector<string> scanList;
    scanList.push_back("scanD");
    scanList.push_back("scanDihed");
    scanList.push_back("scanCsAX");
    scanList.push_back("scanCsAY");
    scanList.push_back("scanCsAZ");
    scanList.push_back("scanCsBX");
    scanList.push_back("scanCsBY");
    scanList.push_back("scanCsBZ");

    char xx[200];
    
    cout << "load atLib" << endl;
    AtomLib* atLib = new AtomLib();
    string augc = "AUGC";

    cout << "load et1" << endl;
    ForceFieldPara* para1 = new ForceFieldPara();
    para1->bwTag = "xtb";
    BasePair6DEnergyTable* et1 = new BasePair6DEnergyTable(para1,true, 1);
    et1->load(para1);

    cout << "load et2" << endl;
    ForceFieldPara* para2 = new ForceFieldPara();
    para2->bwTag = "adj";
    BasePair6DEnergyTable* et2 = new BasePair6DEnergyTable(para2, true, 1);
    et2->load(para2);

    cout << "loat et5" << endl;
    cout << "para" << endl;
    ForceFieldPara* para5 = new ForceFieldPara();
    para5->bwTag = "default";
    para5->wtBp1 = 1.0;
    para5->wtBp2 = 1.0;
    para5->wtHb = 1.0;

    cout << "et" << endl;
    BasePair6DEnergyTable* et5 = new BasePair6DEnergyTable(para5, true, 1);
    et5->load(para5);

    {
        
                vector<double> eneList1; //xtb
                vector<double> eneList2; //adj
               // vector<double> eneList3; //xtb Energy table
 
              //  vector<double> eneList40; //hbond energy kHbOri = 0.2


                vector<double> eneList5; //6dEnergy default
        int typeA = 0;
        int typeB = 1;
        BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
        BaseRotamer* rotB = new BaseRotamer(typeB, atLib);
        string cmFile = "/public/home/pengx/briqx/xtb/sp2000/cluster/AA.cm";
        vector<CsMove> cmList;



        ifstream file;
        file.open(cmFile.c_str(), ios::in);
        string s;
        while(getline(file, s)){
            CsMove cm(s);
            cmList.push_back(cm);
        }
        file.close();



        for(int i=0;i<cmList.size();i++){
            LocalFrame csA;
            LocalFrame csB;
            csB = csA + cmList[i];
            BaseConformer* confA = new BaseConformer(rotA, csA);
            BaseConformer* confB = new BaseConformer(rotB, csB);
            double minD = 99.9;
            for(int j=0;j<rotA->atomNum;j++){
                for(int k=0;k<rotB->atomNum;k++){
                    double d = confA->coords[j].distance(confB->coords[k]);
                    if(d < minD){
                        minD = d;
                    }
                }
            }
        
           

            eneList1.push_back(et1->getEnergy(csA, csB, typeA, typeB, 2, minD));
            eneList2.push_back(et2->getEnergy(csA, csB, typeA, typeB, 2, minD));
            eneList5.push_back(et5->getEnergy(csA, csB, typeA, typeB, 2, minD));
                        
            delete confA;
            delete confB;
        }

        ofstream out;

        string outfile = "/public/home/pengx/briqx/xtb/sp2000/cluster/compareResult/AA.xtb";
        out.open(outfile.c_str(), ios::out);
        for(int i=0;i<eneList1.size();i++){
            out << eneList1[i]*0.25 << endl;
        }
        out.close();

        outfile = "/public/home/pengx/briqx/xtb/sp2000/cluster/compareResult/AA.adj";
        out.open(outfile.c_str(), ios::out);
        for(int i=0;i<eneList2.size();i++){
            out << eneList2[i] << endl;
        }
        out.close();

        outfile = "/public/home/pengx/briqx/xtb/sp2000/cluster/compareResult/AA.default";
        out.open(outfile.c_str(), ios::out);
        for(int i=0;i<eneList5.size();i++){
             out << eneList5[i] << endl;
        }
        out.close();
    }


    
    for(int typeA=0;typeA<4;typeA++){
        for(int typeB=0;typeB<4;typeB++){
            BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
            BaseRotamer* rotB = new BaseRotamer(typeB, atLib);
            for(int clusterID=0;clusterID<4;clusterID++){
                string outlines[168];

                vector<double> eneList1; //xtb
                vector<double> eneList2; //adj
               // vector<double> eneList3; //xtb Energy table
 
              //  vector<double> eneList40; //hbond energy kHbOri = 0.2


                vector<double> eneList5; //6dEnergy default
               
                for(int scan=0;scan<scanList.size();scan++){
                    
                    sprintf(xx, "%d", clusterID);
                    string cmFile = "/public/home/pengx/briqx/basePairScan/"+nbnnb+"/"+augc.substr(typeA,1) + augc.substr(typeB,1) + "-" + string(xx) + "-" + scanList[scan] + ".cm";
                    cout << cmFile << endl;

                    vector<CsMove> cmList;



                    ifstream file;
                    file.open(cmFile.c_str(), ios::in);
                    string s;
                    while(getline(file, s)){
                        CsMove cm(s);
                        cmList.push_back(cm);
                    }
                    file.close();



                    for(int i=0;i<cmList.size();i++){
                        LocalFrame csA;
                        LocalFrame csB;
                        csB = csA + cmList[i];
                        BaseConformer* confA = new BaseConformer(rotA, csA);
                        BaseConformer* confB = new BaseConformer(rotB, csB);
                        double minD = 99.9;
                        for(int j=0;j<rotA->atomNum;j++){
                            for(int k=0;k<rotB->atomNum;k++){
                                double d = confA->coords[j].distance(confB->coords[k]);
                                if(d < minD){
                                    minD = d;
                                }
                            }
                        }
        
                        /*
                        double hbEne0 = 0.0;
                        for(int j=0;j<rotA->polarAtomNum;j++){    
                            for(int k=0;k<rotB->polarAtomNum;k++){
                                hbEne0 += hbET0->getEnergy(rotA->polarAtomUniqueID[j], confA->csPolar[j], rotB->polarAtomUniqueID[k], confB->csPolar[k]);
                            }
                        }

                        double hbEne1 = 0.0;
                        for(int j=0;j<rotA->polarAtomNum;j++){    
                            for(int k=0;k<rotB->polarAtomNum;k++){
                                hbEne1 += hbET1->getEnergy(rotA->polarAtomUniqueID[j], confA->csPolar[j], rotB->polarAtomUniqueID[k], confB->csPolar[k]);
                            }
                        }

                        double hbEne2 = 0.0;
                        for(int j=0;j<rotA->polarAtomNum;j++){    
                            for(int k=0;k<rotB->polarAtomNum;k++){
                                hbEne2 += hbET2->getEnergy(rotA->polarAtomUniqueID[j], confA->csPolar[j], rotB->polarAtomUniqueID[k], confB->csPolar[k]);
                            }
                        }

                        double hbEne3 = 0.0;
                        for(int j=0;j<rotA->polarAtomNum;j++){    
                            for(int k=0;k<rotB->polarAtomNum;k++){
                                hbEne3 += hbET3->getEnergy(rotA->polarAtomUniqueID[j], confA->csPolar[j], rotB->polarAtomUniqueID[k], confB->csPolar[k]);
                            }
                        }                        

                        double hbEne4 = 0.0;
                        for(int j=0;j<rotA->polarAtomNum;j++){    
                            for(int k=0;k<rotB->polarAtomNum;k++){
                                hbEne4 += hbET4->getEnergy(rotA->polarAtomUniqueID[j], confA->csPolar[j], rotB->polarAtomUniqueID[k], confB->csPolar[k]);
                            }
                        }

                        double hbEne5 = 0.0;
                        for(int j=0;j<rotA->polarAtomNum;j++){    
                            for(int k=0;k<rotB->polarAtomNum;k++){
                                hbEne5 += hbET5->getEnergy(rotA->polarAtomUniqueID[j], confA->csPolar[j], rotB->polarAtomUniqueID[k], confB->csPolar[k]);
                            }
                        }                        
                        eneList3.push_back(xtbEt->getEnergy(csA, csB, typeA, typeB, minD) * 0.25);
                        
                        eneList40.push_back(hbEne0);
                        eneList41.push_back(hbEne1);
                        eneList42.push_back(hbEne2);
                        eneList43.push_back(hbEne3);
                        eneList44.push_back(hbEne4);
                        eneList45.push_back(hbEne5);
                        */

                        eneList1.push_back(et1->getEnergy(csA, csB, typeA, typeB, 2, minD));
                        eneList2.push_back(et2->getEnergy(csA, csB, typeA, typeB, 2, minD));
                        eneList5.push_back(et5->getEnergy(csA, csB, typeA, typeB, 2, minD));
                        //eneList8.push_back(et8->getEnergy(csA, csB, typeA, typeB, 2, minD));

                        delete confA;
                        delete confB;
                    }
                }

                ofstream out;


                sprintf(xx, "%c%c-%d.ene", augc[typeA], augc[typeB], clusterID);


                string outfile = "/public/home/pengx/briqx/basePairScan/xtb/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList1.size();i++){
                    out << eneList1[i]*0.25 << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/adj/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList2.size();i++){
                    out << eneList2[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/default/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList5.size();i++){
                    out << eneList5[i] << endl;
                }
                out.close();

                /*
                outfile = "/public/home/pengx/briqx/basePairScan/hbond2/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList42.size();i++){
                    out << eneList42[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/hbond3/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList43.size();i++){
                    out << eneList43[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/hbond4/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList44.size();i++){
                    out << eneList44[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/hbond5/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList45.size();i++){
                    out << eneList45[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/merge6D/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList5.size();i++){
                    out << eneList5[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/bw1/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList6.size();i++){
                    out << eneList6[i] << endl;
                }
                out.close();


                outfile = "/public/home/pengx/briqx/basePairScan/bw2/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList7.size();i++){
                    out << eneList7[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/bw3/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList8.size();i++){
                    out << eneList8[i] << endl;
                }
                out.close();
                */

            }
            delete rotA;
            delete rotB;
        }
    }

    delete atLib;
    delete para5;

    delete et5;



}