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
#include <array>

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
cout << "v2 values after =" << endl;
for(auto& iter: v2) {
    cout << iter <<endl;
}
v2[0] = 10;
cout << "v2[0] = " << v2[0] << ", v1[0] = " << v1[0] <<endl;

auto v5 = move(v1); 
cout << "v5 values after move =" << endl;
for(auto& iter: v5) {
    cout << iter <<endl;
}
cout << "v1 values after move =" << endl;
for(auto& iter: v1) {
    cout << iter <<endl;
}
cout << "test complete" << endl;
return EXIT_SUCCESS;
}

int test2() {
    map<int, set<int> > testMap;
    testMap[1].emplace(10);
    testMap[1].emplace(20);
    for (auto iter = testMap[1].begin(); iter != testMap[1].end(); ++iter) {
        cout << *iter << endl;
    }
    return EXIT_SUCCESS;
}

int test3() {
    int a = 1;
    int b = 2;
    array<int*,2> ar1{&a, &b}, ar2{&a, &b};
    map<array<int*,2>, double> map1;
    map1.emplace(ar1,0.314);
    cout << map1.at(ar2) << endl;
    return EXIT_SUCCESS;

}

// int main() {
//     test3();
//     return EXIT_SUCCESS;
// }

int main(int argc, char** argv) {
    if(argc == 1) {
        cout<< "No arg given" << endl;
        return EXIT_SUCCESS;
    }
}