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
    string path = datapath() + "sphere";
    string outPath = datapath() + "../binaryData/sphere";
    if(!makeDirs(outPath)) {
        throw("[Error]Unable to create " + outPath);
    }
    vector<string> spherePointNum{"100", "200", "300", "500", "600", "800", "1000", "2000", "5000"};
    clock_t start, end;
    double timeTxt, timeBinary;
    vector<int> ct = {BINARY_COL_TYPE_INT, BINARY_COL_TYPE_INT};  // spherexx.index
    int lspn = spherePointNum.size();
    std::cout << "reading " << path  << " .index as TXT..."<<endl;
    for(int i=0; i<lspn; i++) {
        string curPath = path + "/sphere" + spherePointNum[i] + ".index";
        auto* bb = new BinaryBook(1, 0, ct);
        string tN = "sphere" + spherePointNum[i] + ".index";
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
        start = clock();
        ifstream ins;
        ins.open(curPath, ios::in);
        int keyIndex, spIndex;
        while(ins >> keyIndex >> spIndex) {
            get_if<vector<int> >(data_vec[0])->emplace_back(keyIndex);
            get_if<vector<int> >(data_vec[1])->emplace_back(spIndex);
        }
        ins.close();
        bb->addTable(tN, move(data_vec));
        for(int k=0; k<lct; k++) {
            delete data_vec[k];
        }
        if(!spherePointNum[i].compare("1000")) {
            cout << "sphere1000.index Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at("sphere1000.index")->cols[0]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("sphere1000.index")->cols[1]))[333] << endl;
        }
        end = clock();
        timeTxt = (double)(end-start) / CLOCKS_PER_SEC;
        cout << "time consumed for sphere" << spherePointNum[i] <<" index: " << timeTxt << " seconds" << endl;
        string outFile = outPath + "/" + tN;
        if(bwb) {
            ofstream outs;
            outs.open(outFile, ios::out | ios::binary);
            bb->write(outs);
            outs.close();
        }
        delete bb;
        cout << "reading " << outFile  << " as binary..."<<endl;
        start = clock();
        ifstream ins2;
        ins2.open(outFile,ios::in | ios::binary);
        bb = new BinaryBook;
        bb->read(ins2);
        ins2.close();
        end = clock();
        timeBinary = (double)(end-start) / CLOCKS_PER_SEC;
        cout << "time consumed: " << timeBinary << " seconds" << endl;
        if(!spherePointNum[i].compare("1000")) {
            cout << "sphere1000.index Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at("sphere1000.index")->cols[0]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("sphere1000.index")->cols[1]))[333] << endl;
        }
        cout << tN << ": Reading TXT is " << timeTxt/timeBinary << "x slower than reading binary." << endl;
    }

    ct = {BINARY_COL_TYPE_INT, BINARY_COL_TYPE_INT, BINARY_COL_TYPE_INT, BINARY_COL_TYPE_INT, 
          BINARY_COL_TYPE_DOUBLE, BINARY_COL_TYPE_DOUBLE, BINARY_COL_TYPE_DOUBLE};
    cout << "reading " << path  << " .index_nearest3Pts as TXT..."<<endl;
    for(int i=0; i<lspn; i++) {
        string curPath = path + "/sphere" + spherePointNum[i] + ".index_nearest3Pts";
        auto* bb = new BinaryBook(1, 0, ct);
        string tN = "sphere" + spherePointNum[i] + ".index_nearest3Pts";
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
        start = clock();
        ifstream ins;
        ins.open(curPath, ios::in);
        int keyIndex, sp1, sp2, sp3;
        double wt1, wt2, wt3;
        while(ins >> keyIndex >> sp1 >> sp2 >> sp3 >> wt1 >> wt2 >> wt3) {
            get_if<vector<int> >(data_vec[0])->emplace_back(keyIndex);
            get_if<vector<int> >(data_vec[1])->emplace_back(sp1);
            get_if<vector<int> >(data_vec[2])->emplace_back(sp2);
            get_if<vector<int> >(data_vec[3])->emplace_back(sp3);
            get_if<vector<double> >(data_vec[4])->emplace_back(wt1);
            get_if<vector<double> >(data_vec[5])->emplace_back(wt2);
            get_if<vector<double> >(data_vec[6])->emplace_back(wt3);
        }
        ins.close();
        bb->addTable(tN, move(data_vec));
        for(int k=0; k<lct; k++) {
            delete data_vec[k];
        }
        if(!spherePointNum[i].compare("1000")) {
            cout << "sphere1000.index_nearest3Pts Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[0]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[1]))[333] 
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[2]))[333] 
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[3]))[333] 
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at(tN)->cols[4]))[333] 
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at(tN)->cols[5]))[333] 
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at(tN)->cols[6]))[333] 
             << endl;
        }
        end = clock();
        timeTxt = (double)(end-start) / CLOCKS_PER_SEC;
        cout << "time consumed for sphere" << spherePointNum[i] <<" index_nearest3Pts: " << timeTxt << " seconds" << endl;
        string outFile = outPath + "/" + tN;
        if(bwb) {
            ofstream outs;
            outs.open(outFile, ios::out | ios::binary);
            bb->write(outs);
            outs.close();
        }
        delete bb;
        cout << "reading " << outFile  << " as binary..."<<endl;
        start = clock();
        ifstream ins2;
        ins2.open(outFile,ios::in | ios::binary);
        bb = new BinaryBook;
        bb->read(ins2);
        ins2.close();
        end = clock();
        timeBinary = (double)(end-start) / CLOCKS_PER_SEC;
        cout << "time consumed: " << timeBinary << " seconds" << endl;
        if(!spherePointNum[i].compare("1000")) {
            cout << "sphere1000.index_nearest3Pts Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[0]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[1]))[333] 
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[2]))[333] 
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at(tN)->cols[3]))[333] 
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at(tN)->cols[4]))[333] 
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at(tN)->cols[5]))[333] 
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at(tN)->cols[6]))[333] 
             << endl;
        }
        cout << tN << ": Reading TXT is " << timeTxt/timeBinary << "x slower than reading binary." << endl;
    }

    ct = {BINARY_COL_TYPE_INT, BINARY_COL_TYPE_INT, BINARY_COL_TYPE_DOUBLE, BINARY_COL_TYPE_INT};
	string augc = "AUGC";
    path = datapath() + "pairEne";
    outPath = datapath() + "../binaryData/pairEne";
    if(!makeDirs(outPath)) {
        throw("[Error]Unable to create " + outPath);
    }
    string sepType[2] = {"nb", "nnb"}; 
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