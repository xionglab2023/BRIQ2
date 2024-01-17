/**
 * @file dirOperation.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Directory operation functions
 * @version 0.1
 * @date 2024-01-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef DATAIO_DIROPERATION_H_
#define DATAIO_DIROPERATION_H_

#include <dirent.h>
#include <vector>
#include <string>

namespace NSPdataio {
    using namespace std;
    
    /**
     * @brief find all files in directory and fill fileNames (without path) to fileNames.
     * The search is not iterative.
     * 
     * @param[in] path: The directory to search
     * @param[out] fileNames: The vector to store the filenames found 
     * @return true: Search completed successfully
     * @return false: Can not open directory
     */
    bool findFilesInDir(string& path, vector<string>& fileNames);

    /**
     * @brief Iteratively create directories
     * 
     * @param[in] dir: The target directory in full path, the intermediate directories will
     * be created automatically.
     * @return true: Successfully create all directories. 
     * @return false: Creation failure on any of the directories
     */
    bool makeDirs(string& dir);
}
#endif /* DATAIO_DIROPERATION_H_ */