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

    ForceFieldPara* para6 = new ForceFieldPara();
    para6->bwTag = "bw1";
    para6->wtBp1 = 1.0;
    para6->wtBp2 = 1.0;
    para6->wtHb = 1.0;

    cout << "load et6" << endl;
    BasePair6DEnergyTable* et6 = new BasePair6DEnergyTable(para6, true, 1);
    et6->load(para6);

    ForceFieldPara* para7 = new ForceFieldPara();

    cout << "load et7" << endl;
    para7->bwTag = "bw2";
    para7->wtBp1 = 1.0;
    para7->wtBp2 = 1.0;
    para7->wtHb = 1.0;
    BasePair6DEnergyTable* et7 = new BasePair6DEnergyTable(para7, true, 1);
    et7->load(para7);
    ForceFieldPara* para8 = new ForceFieldPara();

    cout << "load et8" << endl;
    para8->bwTag = "bw3";
    para8->wtBp1 = 1.0;
    para8->wtBp2 = 1.0;
    para8->wtHb = 1.0;
    BasePair6DEnergyTable* et8 = new BasePair6DEnergyTable(para8, true, 1);
    et8->load(para8);

    cout << "load hbet" << endl;
    ForceFieldPara* paraHb0 = new ForceFieldPara();
    paraHb0->kHbOri = 0.2;
    HbondEnergy* hbET0 = new HbondEnergy(paraHb0);

    ForceFieldPara* paraHb1 = new ForceFieldPara();
    paraHb1->kHbOri = 0.4;
    HbondEnergy* hbET1 = new HbondEnergy(paraHb1);

    ForceFieldPara* paraHb2 = new ForceFieldPara();
    paraHb2->kHbOri = 0.6;
    HbondEnergy* hbET2 = new HbondEnergy(paraHb2);

    ForceFieldPara* paraHb3 = new ForceFieldPara();
    paraHb3->kHbOri = 0.8;
    HbondEnergy* hbET3 = new HbondEnergy(paraHb3);

    ForceFieldPara* paraHb4 = new ForceFieldPara();
    paraHb4->kHbOri = 1.0;
    HbondEnergy* hbET4 = new HbondEnergy(paraHb4);

    ForceFieldPara* paraHb5 = new ForceFieldPara();
    paraHb5->kHbOri = 1.2;
    HbondEnergy* hbET5 = new HbondEnergy(paraHb5);

    cout << "load xtbet " << endl;
    Xtb6dEnergy* xtbEt = new Xtb6dEnergy(para5);

    for(int typeA=0;typeA<4;typeA++){
        for(int typeB=0;typeB<4;typeB++){
            BaseRotamer* rotA = new BaseRotamer(typeA, atLib);
            BaseRotamer* rotB = new BaseRotamer(typeB, atLib);
            for(int clusterID=0;clusterID<4;clusterID++){
                string outlines[168];


                vector<double> eneList3; //xtb Energy table
 
                vector<double> eneList40; //hbond energy kHbOri = 0.2
                vector<double> eneList41; //hbond energy kHbOri = 0.4
                vector<double> eneList42; //hbond energy kHbOri = 0.6
                vector<double> eneList43; //hbond energy kHbOri = 0.8
                vector<double> eneList44; //hbond energy kHbOri = 1.0
                vector<double> eneList45; //hbond energy kHbOri = 1.2
                vector<double> eneList46; //hbond energy kHbOri = 1.4

                vector<double> eneList5; //6dEnergy default
                vector<double> eneList6; //bw1
                vector<double> eneList7; //bw2
                vector<double> eneList8; //bw3
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

                        eneList5.push_back(et5->getEnergy(csA, csB, typeA, typeB, 2, minD));
                        eneList6.push_back(et6->getEnergy(csA, csB, typeA, typeB, 2, minD));
                        eneList7.push_back(et7->getEnergy(csA, csB, typeA, typeB, 2, minD));
                        eneList8.push_back(et8->getEnergy(csA, csB, typeA, typeB, 2, minD));

                        delete confA;
                        delete confB;
                    }
                }

                ofstream out;
                sprintf(xx, "%c%c-%d.ene", augc[typeA], augc[typeB], clusterID);
                string outfile = "/public/home/pengx/briqx/basePairScan/xtbEt/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList3.size();i++){
                    out << eneList3[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/hbond0/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList40.size();i++){
                    out << eneList40[i] << endl;
                }
                out.close();

                outfile = "/public/home/pengx/briqx/basePairScan/hbond1/nnb/" + string(xx);
                out.open(outfile.c_str(), ios::out);
                for(int i=0;i<eneList41.size();i++){
                    out << eneList41[i] << endl;
                }
                out.close();

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

            }
            delete rotA;
            delete rotB;
        }
    }

    delete atLib;
    delete para5;
    delete para6;
    delete para7;
    delete para8;
    delete et5;
    delete et6;
    delete et7;
    delete et8;
    delete hbET0;
    delete hbET1;
    delete hbET2;
    delete hbET3;
    delete hbET4;
    delete hbET5;
    delete xtbEt;


}