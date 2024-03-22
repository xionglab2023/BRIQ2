#ifndef XY_H_
#define XY_H_

#include <string>
#include <vector>
#include <cmath>
#include "stdio.h"
#include <iostream>

namespace NSPgeometry{


struct XY{
    double x_;
    double y_;

    XY(double ix=0.0, double iy=0.0) : 
        x_(ix), y_(iy) {
            ;
    }

    XY& operator=(const XY& other){
        x_ = other.x_;
        y_ = other.y_;
        return *this;
    }

    double distance(XY& other){
        double dx = x_ - other.x_;
        double dy = y_ - other.y_;
        return sqrt(dx*dx+dy*dy);
    }

};

}

#endif