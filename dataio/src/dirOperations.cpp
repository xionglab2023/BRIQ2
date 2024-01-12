/**
 * @file dirOperations.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Implementation of dirOperations.h
 * @version 0.1
 * @date 2024-01-11
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#include "dataio/dirOperations.h"
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace NSPdataio {

    using namespace std;

    bool findFilesInDir(string& path, vector<string>& fileNames) {
        DIR* pDIR;
        struct dirent* pDirent;
        if(!(pDIR = opendir(path.c_str()))) {
            cerr << "Fail to open " << path <<endl;
            return false;
        }
        while(pDirent = readdir(pDIR)) {
            if(strcmp(pDirent->d_name, ".") && strcmp(pDirent->d_name, "..")) {
                fileNames.emplace_back(pDirent->d_name);
            }
        }
        closedir(pDIR);
        return true;
    }

    
    bool makeDirs(string& dir) {
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
}