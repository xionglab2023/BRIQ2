/**
 * @file txt2binaryMove.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Convert NuMoveSet library (pairMove2) to binary file.
 * @version 0.1
 * @date 2024-01-07
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#include "dataio/datapaths.h"
#include "dataio/binaryTable.h"
#include "dataio/dirOperations.h"
#include "tools/CmdArgs.h"
#include <dirent.h>
#include <cstring>
#include <time.h>
// #include <direct.h>  // For Windows only
#include <sys/stat.h>

using namespace NSPdataio;
using namespace NSPtools;
using namespace std;


int main(int argc, char** argv) {
    CmdArgs cmdArgs{argc, argv};
    bool bwb = false;  // flag of write binary
    if(cmdArgs.specifiedOption("-wb")) {
        bwb = true;
    }

    vector<int> ct = {BINARY_COL_TYPE_INT, BINARY_COL_TYPE_INT, BINARY_COL_TYPE_DOUBLE, BINARY_COL_TYPE_INT};
	string augc = "AUGC";
    string path = datapath() + "pairEneCG";
    string outPath = datapath() + "../binaryData/pairEneCG";
    if(!makeDirs(outPath)) {
        throw("[Error]Unable to create " + outPath);
    }
    string sepType[2] = {"nb", "nnb"}; 
    clock_t start,end;
    for(string st:sepType) {
        double timeTxt, timeBinary;
        string curPath = path + "/" + st;

        cout << "reading " << curPath  << " as TXT..."<<endl;
        start = clock();
        auto* bb = new BinaryBook(1, 0, ct);
        for(int i=0; i<4; i++) {
            for(int j=0; j<4; j++) {
                string pairType = augc.substr(i,1) + augc.substr(j,1);
                string tN = pairType + ".ene";
                vector<variant<vector<int>,vector<double>>*> data_vec;
                int lct = ct.size();
                for(int k=0; k<lct; k++) {
                    switch(ct[k]) {
                        case BINARY_COL_TYPE_INT :
                            data_vec.emplace_back(new variant<vector<int>,vector<double>>{vector<int>{}});
                            break;
                        case BINARY_COL_TYPE_DOUBLE :
                            data_vec.emplace_back(new variant<vector<int>,vector<double>>{vector<double>{}});
                            break;
                        default:
                            throw("[Error] Unrecognized column type");
                    }
                }
                ifstream ins;
                ins.open(curPath + "/" + tN, ios::in);
                int indexA, indexB, clusterID;
                double ene;
                while(ins >> indexA >> indexB >> ene >> clusterID) {
                    get_if<vector<int> >(data_vec[0])->emplace_back(indexA);
                    get_if<vector<int> >(data_vec[1])->emplace_back(indexB);
                    get_if<vector<double> >(data_vec[2])->emplace_back(ene);
                    get_if<vector<int> >(data_vec[3])->emplace_back(clusterID);
                }
                ins.close();
                bb->addTable(tN, move(data_vec));
                for(int k=0; k<lct; k++) {
                    delete data_vec[k];
                }
            }
        }
        end = clock();
        timeTxt = (double)(end-start) / CLOCKS_PER_SEC;
        cout << "time consumed: " << timeTxt << " seconds" << endl;
        cout << "AG.ene Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at("AG.ene")->cols[0]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("AG.ene")->cols[1]))[333]
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at("AG.ene")->cols[2]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("AG.ene")->cols[3]))[333] << endl;
        string outFile = outPath + "/" + st;
        if(bwb) {
            ofstream outs;
            outs.open(outFile, ios::out | ios::binary);
            if(!outs.is_open()) {
                throw("Unable to open " + outFile);
            }
            bb->write(outs);
            outs.close();
        }
        delete bb;
        cout << "reading " << outFile  << " as binary..."<<endl;
        start = clock();
        ifstream ins2;
        ins2.open(outFile,ios::in | ios::binary);
        if(!ins2.is_open()) {
            throw("Unable to open " + outFile);
        }
        bb = new BinaryBook;
        bb->read(ins2);
        ins2.close();
        end = clock();
        timeBinary = (double)(end-start) / CLOCKS_PER_SEC;
        cout << "time consumed: " << timeBinary << " seconds" << endl;
        cout << "AG.ene Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at("AG.ene")->cols[0]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("AG.ene")->cols[1]))[333]
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at("AG.ene")->cols[2]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("AG.ene")->cols[3]))[333] << endl;
        cout << "Reading TXT is " << timeTxt/timeBinary << "x slower than reading binary." << endl;
    }
}