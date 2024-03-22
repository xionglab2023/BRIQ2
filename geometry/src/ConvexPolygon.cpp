#include "geometry/ConvexPolygon.h"

namespace NSPgeometry{

    ConvexPolygon::ConvexPolygon(vector<XY>& pointList){
        for(int i=0;i<pointList.size();i++){
            this->vertexes.push_back(pointList[i]);
        }
        this->vertexNum = pointList.size();

		if(vertexNum == 0) 
			return;
        for(int i=0;i<pointList.size()-1;i++){
            this->edges.push_back(Line2D(vertexes[i], vertexes[i+1]));
        }
        this->edges.push_back(Line2D(vertexes[vertexNum-1], vertexes[0]));
    }

    int ConvexPolygon::insideTheConvex(XY& point){
		for(int i=0;i<vertexNum;i++) {
			if(edges[i].containPoint(point))
				return 0;
		}
		
		double maxX = 0.0;
		double maxY = 0.0;
		
		for(int i=0;i<vertexNum;i++) {
			if(vertexes[i].x_ > maxX) maxX = vertexes[i].x_;
			if(vertexes[i].y_ > maxY) maxY = vertexes[i].y_;
		}
		
		maxX += 1;
		maxY += 1;
		
		Line2D line = Line2D(point, XY(maxX,maxY));
		int jointNum = 0;
		for(int i=0;i<vertexNum;i++) {
			jointNum += this->edges[i].jointPointNum(line);
		}
		
		if(jointNum == 1)
			return 1;
		else
			return -1;
    }

    double ConvexPolygon::getArea(){
		double area = 0.0;
		for(int i=0;i<vertexNum;i++) {
			XY a = edges[i].from;
			XY b = edges[i].to;
			area += (b.x_ - a.x_) * (b.y_ + a.y_) * 0.5;
		}
		return abs(area);	
    }

    bool ConvexPolygon::overlap(ConvexPolygon& other){
        if(this->vertexNum == 0 || other.vertexNum == 0) return false;

		int aInsideB = 0;
		int bInsideA = 0;
		
		for(int i=0;i<this->vertexNum;i++) {
			if(other.insideTheConvex(this->vertexes[i]) > 0)
				aInsideB ++;
		}
		for(int i=0;i<other.vertexNum;i++) {
			if(insideTheConvex(other.vertexes[i]) > 0)
				bInsideA ++;
		}
		
		if(aInsideB > 0 || bInsideA > 0)
			return true;
		
		int jointNum = 0;
		
		for(int i=0;i<this->vertexNum;i++) {
			Line2D edgeA = this->edges[i];
			for(int j=0;j<other.vertexNum;j++) {
				Line2D edgeB = other.edges[j];
				jointNum += edgeA.jointPointNum(edgeB);
			
			}
		}
		if(jointNum > 0 )
			return true;
		else 
			return false;        
    }
}