/**
 * @file test.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Test C++ usage
 * @version udef
 * @date: 2023/08/09
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

#include <vector>
#include <iostream>
#include <map>
#include <set>

using namespace std;

int test1() {
vector<int> v1 = {1, 2, 3};
vector<int> v2 = {4, 5, 6};
vector<vector<int> > v3 = {{11,12,13},
                           {21,22,23},
                           {31,32,33}};
vector<vector<int> > v4;

v2=v1;
v4=v3;
cout << "test complete" << endl;
return EXIT_SUCCESS;
}

int main() {
    map<int, set<int> > testMap;
    testMap[1].emplace(10);
    testMap[1].emplace(20);
    for (auto iter = testMap[1].begin(); iter != testMap[1].end(); ++iter) {
        cout << *iter << endl;
    }
    return EXIT_SUCCESS;
}