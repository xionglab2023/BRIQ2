/*
 * PhosphateRotamerLib.cpp
 *
 *  Created on: 2023��8��14��
 *      Author: nuc
 */

#include <model/PhosphateRotamerLib.h>

namespace NSPmodel {

PhosphateRotamerLib::PhosphateRotamerLib() {
	// TODO Auto-generated constructor stub
	for(int i=0;i<360;i++){
		for(int j=0;j<360;j++){
			this->prLib[i][j] = new PhosphateRotamer(i+0.5, j+0.5);
		}
	}
}

PhosphateRotamerLib::~PhosphateRotamerLib() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<360;i++){
		for(int j=0;j<360;j++){
			delete prLib[i][j];
		}
	}

}

} /* namespace NSPmodel */
