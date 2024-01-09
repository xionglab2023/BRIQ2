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
#include "tools/CmdArgs.h"
#include <dirent.h>
#include <cstring>
#include <time.h>
// #include <direct.h>  // For Windows only
#include <sys/stat.h>

using namespace NSPdataio;
using namespace NSPtools;
using namespace std;

/**
 * @brief find all files in directory and fill fileNames (without path) to fileNames
 * 
 * @param path 
 * @param fileNames 
 */
void findFilesInDir(string& path, vector<string>& fileNames) {
    DIR* pDIR;
    struct dirent* pDirent;
    if(!(pDIR = opendir(path.c_str()))) {
        throw("Fail to open " + path);
    }
    while(pDirent = readdir(pDIR)) {
        if(strcmp(pDirent->d_name, ".") && strcmp(pDirent->d_name, "..")) {
            fileNames.emplace_back(pDirent->d_name);
        }
    }
    closedir(pDIR);
}

/**
 * @brief Iteratively create directories
 * 
 * @param dir full directory path
 * @return true Directories successfully created. 
 * @return false Creation of any of the directories failed.
 */
bool makeDirs(string dir) {
    string curPath, substr;
    substr.reserve(dir.size());
    for(auto it = dir.begin(); it != dir.end(); ++it) {
        substr.push_back(*it);
        if (*it == '/' || it == dir.end() -1) {
            curPath.append(substr);
            if(!opendir(curPath.c_str())) {  // curPath not exist
                if(mkdir(curPath.c_str(), S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH)) {  // create curPath fail
                    return false;
                }
            }
            substr.clear();
        }
    }
    return true;
}


int main(int argc, char** argv) {
    CmdArgs cmdArgs{argc, argv};
    bool bwb = false;  // flag of write binary
    if(cmdArgs.specifiedOption("-wb")) {
        bwb = true;
    }
    string path = datapath() + "pairMove2";
    string outPath = datapath() + "../binaryData/pairMove2";
    if(!makeDirs(outPath)) {
        throw("[Error]Unable to create " + outPath);
    }
    vector<string> pairTypes{"nb","nnb"};
    int npt = pairTypes.size();
    clock_t start, end;
    vector<int> ct = {BINARY_COL_TYPE_INT, BINARY_COL_TYPE_DOUBLE, BINARY_COL_TYPE_INT};
    for(int i=0; i<npt; i++) {
        vector<string> fileNames;
        double timeTxt, timeBinary;
        string curPath = path + "/" + pairTypes[i];
        cout << "reading " << curPath  << " as TXT..."<<endl;
        findFilesInDir(curPath, fileNames);
        auto* bb = new BinaryBook(1, 0, ct);
        start = clock();
        int lf = fileNames.size();
        for(int j=0; j<lf; j++) {
            string tN = fileNames[j];
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
            int moveID, binID;
            double p;
            while(ins >> moveID >> p >> binID) {
                get_if<vector<int> >(data_vec[0])->emplace_back(moveID);
                get_if<vector<double> >(data_vec[1])->emplace_back(p);
                get_if<vector<int> >(data_vec[2])->emplace_back(binID);
            }
            ins.close();
            bb->addTable(tN, move(data_vec));
            for(int k=0; k<lct; k++) {
                delete data_vec[k];
            }
        }
        end = clock();
        timeTxt = (double)(end-start) / CLOCKS_PER_SEC;
        cout << "time consumed: " << timeTxt << " seconds" << endl;
        cout << "AG34.move Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at("AG34.move")->cols[0]))[333]
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at("AG34.move")->cols[1]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("AG34.move")->cols[2]))[333] << endl;
        string outFile = outPath + "/" + pairTypes[i];
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
        cout << "AG34.move Line 334: " << get<BinaryColumn<int>>(*(bb->tables_map.at("AG34.move")->cols[0]))[333]
             << "," << get<BinaryColumn<double>>(*(bb->tables_map.at("AG34.move")->cols[1]))[333]
             << "," << get<BinaryColumn<int>>(*(bb->tables_map.at("AG34.move")->cols[2]))[333] << endl;
        cout << "Reading TXT is " << timeTxt/timeBinary << "x slower than reading binary." << endl;
    }

}