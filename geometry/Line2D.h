

#ifndef GEOMETRY_LINE2D_H_
#define GEOMETRY_ANGLES_H_

#include "stdio.h"
#include "geometry/xy.h"

namespace NSPgeometry{

using namespace std;

class Line2D{
public:
    XY from;
    XY to;

    Line2D(XY from, XY to){
        this->from = from;
        this->to = to;
    }

    Line2D& operator=(const Line2D& other){
        from = other.from;
        to = other.to;
        return *this;
    }

    bool containPoint(XY p) {
		double a = this->from.y_ - this->to.y_;
		double b = this->to.x_ - this->from.x_;
		double c = this->to.x_*this->from.y_ - this->from.x_*this->to.y_;
		if(p.x_ > from.x_ && p.x_ > to.x_)
			return false;
		if(p.x_ < from.x_ && p.x_ < to.x_)
			return false;
		if(p.y_ > from.y_ && p.y_ > to.y_)
			return false;
		if(p.y_ < from.y_ && p.y_ < to.y_)
			return false;
		if(abs(a*p.x_+b*p.y_-c) < 0.0000001)
			return true;
		return false;       
    }

    int jointPointNum(Line2D& other){
		if(containPoint(other.from))
			return 1;
		else if(containPoint(other.to))
			return 1;
		else if(other.containPoint(from))
			return 1;
		else if(other.containPoint(to))
			return 1;
		
		double a = this->from.y_ - this->to.y_;
		double b = this->to.x_ - this->from.x_;
		double c = this->to.x_*this->from.y_ - this->from.x_*this->to.y_;
		
		double d = other.from.y_ - other.to.y_;
		double e = other.to.x_ - other.from.x_;
		double f = other.to.x_*other.from.y_ - other.from.x_*other.to.y_;
		
		if(abs(b*d - a*e) < 0.0000001)
			return 0;
		
		double y = (c*d - a*f)/(b*d - a*e);
		double x = (b*f - c*e)/(b*d - a*e);
		
		if(x > from.x_ && x > to.x_)
			return 0;
		if(x < from.x_ && x < to.x_)
			return 0;
		if(y  > from.y_ && y > to.y_)
			return 0;
		if(y < from.y_ && y < to.y_)
			return 0;
		if(x > other.from.x_ && x > other.to.x_)
			return 0;
		if(x < other.from.x_ && x < other.to.x_)
			return 0;
		if(y  > other.from.y_ && y > other.to.y_)
			return 0;
		if(y < other.from.y_ && y < other.to.y_)
			return 0;
		return 1;     
    }

    virtual ~Line2D();


};


}

#endif