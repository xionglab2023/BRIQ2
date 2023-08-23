/*
 * BaseRotamer.cpp
 *
 */

#include <model/BaseRotamer.h>

namespace NSPmodel {

	BaseRotamer::BaseRotamer(int baseType, AtomLib* atLib){
		vector<string> baseTypes;

		baseTypes.push_back("A");
		baseTypes.push_back("U");
		baseTypes.push_back("G");
		baseTypes.push_back("C");
		baseTypes.push_back("DA");
		baseTypes.push_back("DT");
		baseTypes.push_back("DG");
		baseTypes.push_back("DC");

		vector<string>* scNames = atLib->baseScAtomNames[baseType];
		this->baseType = baseType;

		for(int i=0;i<scNames->size();i++){
			string name = baseTypes[baseType]+"-"+scNames->at(i);
			uniqueIDs[i] = atLib->uniqueNameToID(name);
		}

		LocalFrame cs0;

		if(baseType == 0 || baseType == 4){
			atomNum = 10; //N9, C8, N7, C5, C6, N6, N1, C2, N3, C4
			coordsLocal[0] = XYZ(1.468,   0.000,   0.000); //N9
			coordsLocal[1] = XYZ(2.306,  -1.084,   0.000); //C8
			coordsLocal[2] = XYZ(3.577,  -0.770,   0.000); //N7
			coordsLocal[3] = XYZ(3.576,   0.615,   0.000); //C5
			coordsLocal[4] = XYZ(4.614,   1.556,   0.000); //C6
			coordsLocal[5] = XYZ(5.904,   1.230,   0.000); //N6
			coordsLocal[6] = XYZ(4.276,   2.861,   0.000); //N1
			coordsLocal[7] = XYZ(2.980,   3.184,   0.000); //C2
			coordsLocal[8] = XYZ(1.914,   2.391,   0.000); //N3
			coordsLocal[9] = XYZ(2.285,   1.103,   0.000); //C4

			this->polarAtomNum = 4;
			polarAtomIndex[0] = 2; //N7
			polarAtomIndex[1] = 5; //N6
			polarAtomIndex[2] = 6; //N1
			polarAtomIndex[3] = 8; //N3
			polarAtomUniqueID[0] = uniqueIDs[polarAtomIndex[0]];
			polarAtomUniqueID[1] = uniqueIDs[polarAtomIndex[1]];
			polarAtomUniqueID[2] = uniqueIDs[polarAtomIndex[2]];
			polarAtomUniqueID[3] = uniqueIDs[polarAtomIndex[3]];

			LocalFrame cs1 = generateLocalFrameResidueStyle(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
			LocalFrame cs2 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
			LocalFrame cs3 = generateLocalFrameResidueStyle(coordsLocal[4], coordsLocal[6], coordsLocal[7]);
			LocalFrame cs4 = generateLocalFrameResidueStyle(coordsLocal[7], coordsLocal[8], coordsLocal[9]);
			this->polarCmList[0] = cs1 - cs0;
			this->polarCmList[1] = cs2 - cs0;
			this->polarCmList[2] = cs3 - cs0;
			this->polarCmList[3] = cs4 - cs0;

		}
		else if(baseType == 1){
			atomNum = 8; //N1, C2, O2, N3, C4, O4, C5, C6
			coordsLocal[0] = XYZ(1.478,   0.000,   0.000); //N1
			coordsLocal[1] = XYZ(2.122,   1.221,   0.000); //C2
			coordsLocal[2] = XYZ(1.528,   2.282,   0.000); //O2
			coordsLocal[3] = XYZ(3.491,   1.159,   0.000); //N3
			coordsLocal[4] = XYZ(4.265,   0.020,   0.000); //C4
			coordsLocal[5] = XYZ(5.490,   0.123,   0.000); //O4
			coordsLocal[6] = XYZ(3.526,  -1.204,   0.000); //C5
			coordsLocal[7] = XYZ(2.191,  -1.173,   0.000); //C6

			this->polarAtomNum = 3;
			polarAtomIndex[0] = 2;
			polarAtomIndex[1] = 3;
			polarAtomIndex[2] = 5;
			polarAtomUniqueID[0] = uniqueIDs[polarAtomIndex[0]];
			polarAtomUniqueID[1] = uniqueIDs[polarAtomIndex[1]];
			polarAtomUniqueID[2] = uniqueIDs[polarAtomIndex[2]];
			LocalFrame cs1 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[2]);
			LocalFrame cs2 = generateLocalFrameResidueStyle(coordsLocal[1], coordsLocal[3], coordsLocal[4]);
			LocalFrame cs3 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
			this->polarCmList[0] = cs1 - cs0;
			this->polarCmList[1] = cs2 - cs0;
			this->polarCmList[2] = cs3 - cs0;
		}
		else if(baseType == 5){
			atomNum = 9; //N1, C2, O2, N3, C4, O4, C5, C6
			coordsLocal[0] = XYZ(1.478,   0.000,   0.000); //N1
			coordsLocal[1] = XYZ(2.122,   1.221,   0.000); //C2
			coordsLocal[2] = XYZ(1.528,   2.282,   0.000); //O2
			coordsLocal[3] = XYZ(3.491,   1.159,   0.000); //N3
			coordsLocal[4] = XYZ(4.265,   0.020,   0.000); //C4
			coordsLocal[5] = XYZ(5.490,   0.123,   0.000); //O4
			coordsLocal[6] = XYZ(3.526,  -1.204,   0.000); //C5
			coordsLocal[7] = XYZ(2.191,  -1.173,   0.000); //C6
			coordsLocal[8] = XYZ(4.273,  -2.521,   0.000); //C7
			this->polarAtomNum = 3;
			polarAtomIndex[0] = 2;
			polarAtomIndex[1] = 3;
			polarAtomIndex[2] = 5;
			polarAtomUniqueID[0] = uniqueIDs[polarAtomIndex[0]];
			polarAtomUniqueID[1] = uniqueIDs[polarAtomIndex[1]];
			polarAtomUniqueID[2] = uniqueIDs[polarAtomIndex[2]];
			LocalFrame cs1 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[2]);
			LocalFrame cs2 = generateLocalFrameResidueStyle(coordsLocal[1], coordsLocal[3], coordsLocal[4]);
			LocalFrame cs3 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
			this->polarCmList[0] = cs1 - cs0;
			this->polarCmList[1] = cs2 - cs0;
			this->polarCmList[2] = cs3 - cs0;
		}
		else if(baseType == 2 || baseType == 6){
			atomNum = 11;
			coordsLocal[0] = XYZ(1.468,   0.000,   0.000); //N9
			coordsLocal[1] = XYZ(2.295,  -1.094,   0.000); //C8
			coordsLocal[2] = XYZ(3.560,  -0.779,   0.000); //N7
			coordsLocal[3] = XYZ(3.570,   0.606,   0.000); //C5
			coordsLocal[4] = XYZ(4.655,   1.514,   0.000); //C6
			coordsLocal[5] = XYZ(5.864,   1.265,   0.000); //O6
			coordsLocal[6] = XYZ(4.221,   2.832,   0.000); //N1
			coordsLocal[7] = XYZ(2.909,   3.225,   0.000); //C2
			coordsLocal[8] = XYZ(2.690,   4.543,   0.000); //N2
			coordsLocal[9] = XYZ(1.886,   2.389,   0.000); //N3
			coordsLocal[10] = XYZ(2.287,   1.103,   0.000); //C4
			this->polarAtomNum = 5;
			polarAtomIndex[0] = 2;
			polarAtomIndex[1] = 5;
			polarAtomIndex[2] = 6;
			polarAtomIndex[3] = 8;
			polarAtomIndex[4] = 9;
			polarAtomUniqueID[0] = uniqueIDs[polarAtomIndex[0]];
			polarAtomUniqueID[1] = uniqueIDs[polarAtomIndex[1]];
			polarAtomUniqueID[2] = uniqueIDs[polarAtomIndex[2]];
			polarAtomUniqueID[3] = uniqueIDs[polarAtomIndex[3]];
			polarAtomUniqueID[4] = uniqueIDs[polarAtomIndex[4]];
			LocalFrame cs1 = generateLocalFrameResidueStyle(coordsLocal[1], coordsLocal[2], coordsLocal[3]);
			LocalFrame cs2 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
			LocalFrame cs3 = generateLocalFrameResidueStyle(coordsLocal[4], coordsLocal[6], coordsLocal[7]);
			LocalFrame cs4 = LocalFrame(coordsLocal[6], coordsLocal[7], coordsLocal[8]);
			LocalFrame cs5 = generateLocalFrameResidueStyle(coordsLocal[7], coordsLocal[9], coordsLocal[10]);
			this->polarCmList[0] = cs1 - cs0;
			this->polarCmList[1] = cs2 - cs0;
			this->polarCmList[2] = cs3 - cs0;
			this->polarCmList[3] = cs4 - cs0;
			this->polarCmList[4] = cs5 - cs0;
		}
		else if(baseType == 3 || baseType == 7){
			atomNum = 8;
			coordsLocal[0] = XYZ(1.478,   0.000,   0.000); //N1
			coordsLocal[1] = XYZ(2.151,   1.224,   0.000); //C2
			coordsLocal[2] = XYZ(1.490,   2.271,   0.000); //O2
			coordsLocal[3] = XYZ(3.503,   1.239,   0.000); //N3
			coordsLocal[4] = XYZ(4.178,   0.091,   0.000); //C4
			coordsLocal[5] = XYZ(5.508,   0.150,   0.000); //N4
			coordsLocal[6] = XYZ(3.519,  -1.170,   0.000); //C5
			coordsLocal[7] = XYZ(2.181,  -1.170,   0.000); //C6
			this->polarAtomNum = 3;
			polarAtomIndex[0] = 2;
			polarAtomIndex[1] = 3;
			polarAtomIndex[2] = 5;
			polarAtomUniqueID[0] = uniqueIDs[polarAtomIndex[0]];
			polarAtomUniqueID[1] = uniqueIDs[polarAtomIndex[1]];
			polarAtomUniqueID[2] = uniqueIDs[polarAtomIndex[2]];
			LocalFrame cs1 = LocalFrame(coordsLocal[0], coordsLocal[1], coordsLocal[2]);
			LocalFrame cs2 = generateLocalFrameResidueStyle(coordsLocal[1], coordsLocal[3], coordsLocal[4]);
			LocalFrame cs3 = LocalFrame(coordsLocal[3], coordsLocal[4], coordsLocal[5]);
			this->polarCmList[0] = cs1 - cs0;
			this->polarCmList[1] = cs2 - cs0;
			this->polarCmList[2] = cs3 - cs0;
		}
		else {
			cout << "invalid base type: " << baseType << endl;
			exit(0);
		}

		if(atomNum != scNames->size()){
			cout << "atom number error: " << atomNum << " != "	 << scNames->size() << endl;
		}
	}

	BaseRotamer::~BaseRotamer(){

	}


	BaseConformer::BaseConformer(BaseRotamer* rot, LocalFrame& cs){
		this->rot = rot;
		this->cs1 = cs;

		for(int i=0;i<rot->atomNum;i++){
			coords[i] = local2global(cs1, rot->coordsLocal[i]);
		}
		for(int i=0;i<rot->polarAtomNum;i++){
			csPolar[i] = cs1 + rot->polarCmList[i];
		}
	}

	void BaseConformer::copyValueFrom(BaseConformer* other){
		this->rot = other->rot;
		this->cs1 = other->cs1;

		for(int i=0;i<11;i++){
			this->coords[i] = other->coords[i];
		}
		for(int i=0;i<5;i++){
			this->csPolar[i] = other->csPolar[i];
		}
	}

	void BaseConformer::updateCoords(LocalFrame& cs){
		this->cs1 = cs;
		for(int i=0;i<rot->atomNum;i++){
			coords[i] = local2global(cs, rot->coordsLocal[i]);
		}
		for(int i=0;i<rot->polarAtomNum;i++){
			csPolar[i] = cs1 + rot->polarCmList[i];
		}
	}

	BaseConformer::~BaseConformer(){

	}

} /* namespace NSPforcefield */
