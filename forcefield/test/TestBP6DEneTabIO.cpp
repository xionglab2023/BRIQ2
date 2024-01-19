/**
 * @file TestBP6DEneTabIO.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief test IO performance of BasePair6DEnergyTable with txt and binary data sources.
 * @version 0.1
 * @date 2024-01-11
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#include "forcefield/BasePair6DEnergyTable.h"
#include <sys/time.h>

using namespace NSPforcefield;
using namespace std;

int main() {
    struct timeval start, end;

    cout << "Initializing  ForceFieldPara from txt" << endl;
    gettimeofday(&start, NULL);
    auto* para = new ForceFieldPara();
    gettimeofday(&end, NULL);
    double timeParaTxt = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeParaTxt << " seconds" << endl;

    cout << "Initializing BasePair6DEnergyTable from txt" << endl;
    gettimeofday(&start, NULL);
    auto* set2 = new BasePair6DEnergyTable(para, false);
    gettimeofday(&end, NULL);
    double timeTxt = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeTxt << " seconds" << endl;
    cout << "cm2Key[2233]: " << endl;
    set2->cm2Key.printElem(2333);
    cout << "BasePair6DEnergyTable nnbKeysEnergy AG 763 74 ene: " << set2->nnbKeysEnergy[2*2250+763][74] << endl;
    set2->dump();
    delete set2;

    cout << "Initializing BasePair6DEnergyTable from binaryTable" << endl;
    gettimeofday(&start, NULL);
    set2 = new BasePair6DEnergyTable(para, true,2);
    gettimeofday(&end, NULL);
    double timeBinaryTable = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeBinaryTable << " seconds" << endl;
    cout << "cm2Key[2233]: " << endl;
    set2->cm2Key.printElem(2333);
    cout << "BasePair6DEnergyTable nnbKeysEnergy AG 763 74 ene: " << set2->nnbKeysEnergy[2*2250+763][74] << endl;
    delete set2;
    
    cout << "Initializing BasePair6DEnergyTable from binaryCache" << endl;
    gettimeofday(&start, NULL);
    set2 = new BasePair6DEnergyTable(para, true,1);
    set2->load();
    gettimeofday(&end, NULL);
    double timeBinaryCache = end.tv_sec - start.tv_sec + (double)(end.tv_usec-start.tv_usec) /1e6;
    cout << "time consumed: " << timeBinaryCache << " seconds" << endl;
    cout << "cm2Key[2233]: " << endl;
    set2->cm2Key.printElem(2333);
    cout << "BasePair6DEnergyTable nnbKeysEnergy AG 763 74 ene: " << set2->nnbKeysEnergy[2*2250+763][74] << endl;
    delete set2;
    delete para;
}