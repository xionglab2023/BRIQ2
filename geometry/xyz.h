/*
 * xyz.h
 */

#ifndef XYZ_H_
#define XYZ_H_
#include <string>
#include <vector>
#include <cmath>
#include "stdio.h"
#include <iostream>

namespace NSPgeometry {

struct XYZ {
	double x_;
	double y_;
	double z_;


	XYZ(double ix = 0.0, double iy = 0.0, double iz = 0.0) :
			x_(ix), y_(iy), z_(iz) {
			;
	}

	XYZ (const std::vector<double> & p) : x_(p[0]), y_(p[1]),z_(p[2]){;}

	double squarednorm() const {
		return x_ * x_ + y_ * y_ + z_ * z_;
	}

	void copyValue(const XYZ& other){
		this->x_= other.x_;
		this->y_=other.y_;
		this->z_=other.z_;
	}

	double length() const{
		return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
	}

	double squaredDistance(const XYZ& other){
		return (x_-other.x_)*(x_-other.x_)+(y_-other.y_)*(y_-other.y_)+(z_-other.z_)*(z_-other.z_);
	}

	double distance(const XYZ& other) const{
		return sqrt((x_-other.x_)*(x_-other.x_)+(y_-other.y_)*(y_-other.y_)+(z_-other.z_)*(z_-other.z_));
	}


	XYZ operator-() const {
		return XYZ(-x_,-y_,-z_);
	}

	XYZ& operator=(const XYZ& other);
	double & operator[] (int i) { if(i==0) return x_; else if(i==1) return y_; else if(i==2) return z_;
	else {std::cout <<"XYZ index out of range" <<std::endl;exit(1);}}
	const double & operator[] (int i) const { if(i==0) return x_; else if(i==1) return y_; else if(i==2) return z_;
	else {std::cout <<"XYZ index out of range" <<std::endl;exit(1);}}

	std::string toString() const {
		char s[30];
		sprintf(s, "%8.3f%8.3f%8.3f", x_, y_, z_);
		return (std::string(s));
	}

	~XYZ() {
	}
};

inline double len(const XYZ& p){
	return sqrt(p.x_*p.x_+ p.y_*p.y_+ p.z_*p.z_);
}

inline XYZ& XYZ::operator=(const XYZ& other){
	x_ = other.x_;
	y_ = other.y_;
	z_ = other.z_;
	return *this;
}

inline XYZ operator -(const XYZ &p1, const XYZ &p2) {
	return XYZ(p1.x_ - p2.x_, p1.y_ - p2.y_, p1.z_ - p2.z_);
}

inline XYZ operator +(const XYZ &p1, const XYZ &p2) {
	return XYZ(p1.x_ + p2.x_, p1.y_ + p2.y_, p1.z_ + p2.z_);
}

inline XYZ operator *(const XYZ &p, double d) {
	return XYZ(p.x_ *d, p.y_*d,p.z_*d);
}
inline XYZ operator *(double d, const XYZ &p) {
	return XYZ(p.x_ *d, p.y_*d,p.z_*d);
}
inline XYZ operator /(const XYZ &p, double d) {
	return p * (1.0 / d);
}

inline double operator *(const XYZ &p1, const XYZ &p2){
	return p1.x_ * p2.x_ + p1.y_ * p2.y_ + p1.z_ * p2.z_;
}

inline XYZ operator ^(const XYZ& p1, const XYZ& p2){
	return XYZ(p1.y_ * p2.z_ - p1.z_ * p2.y_,p1.z_ * p2.x_ - p1.x_ * p2.z_,p1.x_ * p2.y_ - p1.y_ * p2.x_);
}

inline XYZ operator ~(const XYZ& p){

     double len = sqrt(p.x_*p.x_+p.y_*p.y_+p.z_*p.z_);
     if(len == 0){
    	 std::cout << "zero vector can't make unit" << std::endl;
    	 exit(1);
     }
     return p*(1.0/len);
}

inline double squareDistance(const XYZ &p1, const XYZ &p2) {
	return (p1.x_ - p2.x_)*(p1.x_ - p2.x_)+(p1.y_ - p2.y_)*(p1.y_ - p2.y_)+(p1.z_ - p2.z_)*(p1.z_ - p2.z_);
}

inline bool isNeighbor(const XYZ &p1, const XYZ &p2, double distanceCutoff) {
	return (p1.x_ - p2.x_)*(p1.x_ - p2.x_)+(p1.y_ - p2.y_)*(p1.y_ - p2.y_)+(p1.z_ - p2.z_)*(p1.z_ - p2.z_) < distanceCutoff*distanceCutoff;
}

inline double dot(const XYZ &p1, const XYZ &p2) {
	return p1.x_ * p2.x_ + p1.y_ * p2.y_ + p1.z_ * p2.z_;
}

inline XYZ cross(const XYZ &p1, const XYZ &p2) {
	XYZ res;
	res.x_ = p1.y_ * p2.z_ - p1.z_ * p2.y_;
	res.y_ = p1.z_ * p2.x_ - p1.x_ * p2.z_;
	res.z_ = p1.x_ * p2.y_ - p1.y_ * p2.x_;
	return res;
}

} //namespace geometry

#endif /* XYZ_H_ */
