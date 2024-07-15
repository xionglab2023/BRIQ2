/*
 * TestBaseLib.cpp
 *
 */


#include "model/LigandLib.h"
#include "geometry/RMSD.h"
#include "stdlib.h"
#include <iostream>

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;

int main(int argc, char** argv){

	LigandLib* lib = new LigandLib();
	lib->printInfo();
	

}



