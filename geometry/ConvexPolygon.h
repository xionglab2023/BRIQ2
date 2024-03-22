#ifndef CONVEXPOLYGON_H_
#define CONVEXPOLYGON_H_

#include "stdio.h"
#include "geometry/Line2D.h"
#include "geometry/xy.h"

namespace NSPgeometry {

using namespace std;

class ConvexPolygon{
public:
    int vertexNum;
    vector<XY> vertexes;
    vector<Line2D> edges;

    ConvexPolygon() {
        this->vertexes;
    }

    ConvexPolygon& operator=(const ConvexPolygon& other){
        vertexNum = other.vertexNum;
        for(int i=0;i<other.vertexNum;i++){
            this->vertexes.push_back(other.vertexes[i]);
            this->edges.push_back(other.edges[i]);
        }
        return *this;
    }

    ConvexPolygon(vector<XY>& pointList);

    int insideTheConvex(XY& point);
    double getArea();
    bool overlap(ConvexPolygon& other);
    void print(){
        cout << "Convex Polygon: " << endl;
        cout << "edge num: " << vertexNum << endl;
        for(int i=0;i<vertexNum;i++){
            printf("%7.3f %7.3f\n", vertexes[i].x_, vertexes[i].y_);
        }
    }
};

}


#endif