/*
 * test.cpp
 */

#include <iostream>
#include <vector>

#include "math/KDTree.h"

    using point_t = std::vector<double>;
    using pointVec = std::vector<double*>;
    using namespace std;

// point_t pt(2);

int main() {
    pointVec points;

    double pts[] = {2.0, 3.0, 5.0, 4.0, 9.0, 6.0, 4.0, 7.0, 8.0, 1.0, 7.0, 2.0};

    for(int i=0;i<6;i++) {
    	double *pt = new double[2];
    	pt[0] = pts[i*2];
    	pt[1] = pts[i*2+1];
    	points.push_back(pt);
    }

    {
        NSPmath::KDTree tree(points,2);
        tree.printTree();
        std::cout << "nearest test\n";


        double* pt = new double[2];
        pt[0] = 7.5;
        pt[1] = 2.5;

        auto res = tree.nearest_(pt);
        for(int i=0;i<res->dim;i++){
        	double s = res->x[i];
        	std::cout << s << std::endl;
        }
    }

    return 0;
}
