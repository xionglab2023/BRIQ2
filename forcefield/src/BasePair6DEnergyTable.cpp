/*
 * BasePair6DEnergyTable.cpp
 *
 */

#include <forcefield/BasePair6DEnergyTable.h>
#include "dataio/dirOperations.h"
#include <string.h>

namespace NSPforcefield {

using namespace NSPdataio;

CsMoveTo6DKey::CsMoveTo6DKey(bool withBinary, int binaryMode) {
	// TODO Auto-generated constructor stub

	#ifdef TIMING
	clock_t start, end;
	start = clock();
	double indexTime, indexN3PTime;
	#endif 

	this->spherePointNum = 2000;
	ifstream file;
	char xx[20];
	sprintf(xx, "%d", spherePointNum);
	if(!withBinary) {
		string fileName = NSPdataio::datapath() + "sphere/sphere" + string(xx) + ".index";
		file.open(fileName, ios::in);
		if(!file.is_open()) {
			cout << "can't open file: " << fileName << endl;
		}
		int keyIndex;
		int spIndex;
		while(file >> keyIndex >> spIndex){
			this->sphereKeyMap[keyIndex] = spIndex;
		}
		file.close();

		#ifdef TIMING
		end = clock();
		indexTime = (double) (end-start)/CLOCKS_PER_SEC;
		cout << "sphere index reading time: " << indexTime << " seconds" << endl;
		#endif

		fileName = NSPdataio::datapath() + "sphere/sphere" + string(xx) + ".index_nearest3Pts";
		file.open(fileName, ios::in);
		if(!file.is_open()) {
			cout << "can't open file: " << fileName << endl;
		}

		int sp1, sp2, sp3;
		double wt1, wt2, wt3;
		int lineId = 0;
		while(file >> keyIndex >> sp1 >> sp2 >> sp3 >> wt1 >> wt2 >> wt3){
			this->keyIdToListID[keyIndex] = lineId;
			this->spIndex1[lineId] = sp1;

			this->spIndex2[lineId] = sp2;
			this->spIndex3[lineId] = sp3;
			this->spWt1[lineId] = wt1;
			this->spWt2[lineId] = wt2;
			this->spWt3[lineId] = wt3;
			lineId ++;
		} 
		file.close();

		#ifdef TIMING
		end = clock();
		indexN3PTime = (double) (end-start)/CLOCKS_PER_SEC - indexTime;
		cout << "sphere index_nearest3Pts reading time: " << indexN3PTime << " seconds" << endl;
		#endif
	} else if(binaryMode==2) {
		string fileName = NSPdataio::datapath() + "../binaryData/sphere/sphere" + string(xx) + ".index";
		file.open(fileName, ios::in|ios::binary);
		if(!file.is_open()) {
			throw("can't open file: " + fileName);
		}
		auto* bb = new BinaryBook;
 		bb->read(file);
		file.close();
		auto* btab = bb->tables_vec[0];
		auto& keyIndexCol = get<BinaryColumn<int>>(*(btab->cols[0]));
		auto& spIndexCol = get<BinaryColumn<int>>(*(btab->cols[1]));
		for(int i=0;i<btab->nRow; i++) {
			this->sphereKeyMap[keyIndexCol[i]] = spIndexCol[i];
		}
		delete bb;

		fileName = NSPdataio::datapath() + "../binaryData/sphere/sphere" + string(xx) + ".index_nearest3Pts";
		file.open(fileName, ios::in|ios::binary);
		if(!file.is_open()) {
			throw("can't open file: " + fileName);
		}
		bb = new BinaryBook;
 		bb->read(file);
		file.close();
		btab = bb->tables_vec[0];
		auto& keyIndexCol2 = get<BinaryColumn<int>>(*(btab->cols[0]));
		for(int i=0;i<btab->nRow; i++) {
			this->keyIdToListID[keyIndexCol2[i]] = i;
		}
		memcpy(this->spIndex1, &(get<BinaryColumn<int>>(*(btab->cols[1])).getVal()[0]), sizeof(int)*btab->nRow);
		memcpy(this->spIndex2, &(get<BinaryColumn<int>>(*(btab->cols[2])).getVal()[0]), sizeof(int)*btab->nRow);
		memcpy(this->spIndex3, &(get<BinaryColumn<int>>(*(btab->cols[3])).getVal()[0]), sizeof(int)*btab->nRow);
		memcpy(this->spWt1, &(get<BinaryColumn<double>>(*(btab->cols[4])).getVal()[0]), sizeof(double)*btab->nRow);
		memcpy(this->spWt2, &(get<BinaryColumn<double>>(*(btab->cols[5])).getVal()[0]), sizeof(double)*btab->nRow);
		memcpy(this->spWt3, &(get<BinaryColumn<double>>(*(btab->cols[6])).getVal()[0]), sizeof(double)*btab->nRow);
		delete bb;
	} else if(binaryMode==1) {
		return;
	}


	int indexA, indexB, indexC, indexD;
	double dA, dB, dC, dD;
	double wA, wB, wC, wD;
	for(int i=0;i<500;i++){
		double d = i*0.03+0.015;
		if(i >= 490){
			indexA = 49;
			indexB = 0;
			dA = (d-indexA*0.3)*3.3333333333;
			dB = 1-dA;
		}
		else {
			indexA = i/10;
			indexB = indexA+1;
			dA = (d-indexA*0.3)*3.3333333333;
			dB = 1-dA;
		}
		for(int j=0;j<450;j++){
			double ang = j*0.8+0.04;
			if(j >= 440){
				indexC = 44;
				indexD = 0;
				dC = (ang - indexC*8)*0.125;
				dD = 1 - dC;
			}
			else {
				indexC = j/10;
				indexD = indexC+1;
				dC = (ang - indexC*8)*0.125;
				dD = 1 - dC;
			}

			//indexA, indexB: index for distance
			//indexC, indexD: index for angle

			this->biIndex[i*1800 + j*4] = indexA*45+indexC;
			this->biIndex[i*1800 + j*4 + 1] = indexA*45+indexD;
			this->biIndex[i*1800 + j*4 + 2] = indexB*45+indexC;
			this->biIndex[i*1800 + j*4 + 3] = indexB*45+indexD;

			this->biWt[i*1800 + j*4] = dB * dD;
			this->biWt[i*1800 + j*4 + 1] = dB * dC;
			this->biWt[i*1800 + j*4 + 2] = dA * dD;
			this->biWt[i*1800 + j*4 + 3] = dA * dC;
		}
	}
	
	#ifdef TIMING
	end = clock();
	double othTime = (double) (end-start)/CLOCKS_PER_SEC - indexTime - indexN3PTime;
	cout << "Other initialization time: " << othTime << " seconds" << endl;
	#endif
}

int CsMoveTo6DKey::dump(ostream& outs) {
	char tint[4] = {'\0'};
	int skmSize = sphereKeyMap.size();
	int ktlSize = keyIdToListID.size();
	outs.write(reinterpret_cast<char*>(&skmSize), sizeof(int));
	outs.write(reinterpret_cast<char*>(&ktlSize), sizeof(int));
	outs.write(reinterpret_cast<char*>(&spherePointNum), sizeof(int));
	outs.write(tint, sizeof(int));
	{
		int amap[skmSize*2];
		int ii = 0;
		for(auto& it:sphereKeyMap){
			amap[ii] = it.first;
			amap[ii+1] = it.second;
			ii+=2;
		}
		outs.write(reinterpret_cast<char*>(amap), sizeof(int)*2*skmSize);
	}
	{
		int amap[ktlSize*2];
		int ii = 0;
		for(auto& it:keyIdToListID){
			amap[ii] = it.first;
			amap[ii+1] = it.second;
			ii+=2;
		}
		outs.write(reinterpret_cast<char*>(amap), sizeof(int)*2*ktlSize);
	}
	outs.write(reinterpret_cast<char*>(spIndex1), 753848*sizeof(int));
	outs.write(reinterpret_cast<char*>(spIndex2), 753848*sizeof(int));
	outs.write(reinterpret_cast<char*>(spIndex3), 753848*sizeof(int));
	outs.write(reinterpret_cast<char*>(spWt1), 753848*sizeof(double));
	outs.write(reinterpret_cast<char*>(spWt2), 753848*sizeof(double));
	outs.write(reinterpret_cast<char*>(spWt3), 753848*sizeof(double));
	outs.write(reinterpret_cast<char*>(biIndex), 900000*sizeof(int));
	outs.write(reinterpret_cast<char*>(biWt), 900000*sizeof(double));

	return EXIT_SUCCESS;
}

int CsMoveTo6DKey::load(istream& ins) {
	int sizes[4]; // skmSize, ktlSize, spherePointNum and filled zeros
	ins.read(reinterpret_cast<char*>(sizes), 4*sizeof(int));
	spherePointNum = sizes[2];
	{
		int memSize = sizes[0]*2;
		int amap[memSize];
		ins.read(reinterpret_cast<char*>(&amap), sizeof(int)*memSize);
		for(int i=0; i<memSize; i=i+2) {
			sphereKeyMap.emplace(amap[i], amap[i+1]);
		}
	}
	{
		int memSize = sizes[1]*2;
		int amap[memSize];
		ins.read(reinterpret_cast<char*>(&amap), sizeof(int)*memSize);
		for(int i=0; i<memSize; i=i+2) {
			keyIdToListID.emplace(amap[i], amap[i+1]);
		}
	}
	ins.read(reinterpret_cast<char*>(spIndex1), sizeof(int)*753848);
	ins.read(reinterpret_cast<char*>(spIndex2), sizeof(int)*753848);
	ins.read(reinterpret_cast<char*>(spIndex3), sizeof(int)*753848);
	ins.read(reinterpret_cast<char*>(spWt1), sizeof(double)*753848);
	ins.read(reinterpret_cast<char*>(spWt2), sizeof(double)*753848);
	ins.read(reinterpret_cast<char*>(spWt3), sizeof(double)*753848);
	ins.read(reinterpret_cast<char*>(biIndex), 900000*sizeof(int));
	ins.read(reinterpret_cast<char*>(biWt), 900000*sizeof(double));

	return EXIT_SUCCESS;
}

void CsMoveTo6DKey::getIndexAndWeight(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[5], double outWeights[4]){

	if(len == 0) {
		int index = 0;

		outIndex[0] = this->biIndex[index];
		outIndex[1] = this->biIndex[index+1];
		outIndex[2] = this->biIndex[index+2];
		outIndex[3] = this->biIndex[index+3];

		outIndex[4] = 0;
		outWeights[0] = this->biWt[index];
		outWeights[1] = this->biWt[index+1];
		outWeights[2] = this->biWt[index+2];
		outWeights[3] = this->biWt[index+3];
		return;
	}

	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	int idX1 = (int)((tBInA.x_ + 1.0)*200);
	if(idX1 == 400) idX1 = 399;
	int idY1 = (int)((tBInA.y_ + 1.0)*200);
	if(idY1 == 400) idY1 = 399;
	int idZ1 = (int)((tBInA.z_ + 1.0)*200);
	if(idZ1 == 400) idZ1 = 399;

	int idX2 = (int)((tAInB.x_ + 1.0)*200);
	if(idX2 == 400) idX2 = 399;
	int idY2 = (int)((tAInB.y_ + 1.0)*200);
	if(idY2 == 400) idY2 = 399;
	int idZ2 = (int)((tAInB.z_ + 1.0)*200);
	if(idZ2 == 400) idZ2 = 399;

	it = sphereKeyMap.find(idX1*160000+idY1*400+idZ1);
	if(it == sphereKeyMap.end()){
		cout << "len: " << len << endl;
		cout << "B in A: " << tBInA.toString() << endl;
		cout << "invalid sphere index: " << idX1 << " " << idY1 << " " << idZ1 << endl;
		exit(1);
	}
	int spIndexA = it->second;

	it = sphereKeyMap.find(idX2*160000+idY2*400+idZ2);
	if(it == sphereKeyMap.end()){
		cout << "len: " << len << endl;
		cout << "A in B: " << tAInB.toString() << endl;
		cout << "invalid sphere index: " << idX2 << " " << idY2 << " " << idZ2 << endl;
		exit(1);
	}
	int spIndexB = it->second;

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;

	int indexLen = (int)(len*33.3333);
	int indexAng = (int)(ang*1.25);
	if(indexAng == 450) indexAng = 0;

	if(indexLen <0 || indexLen >= 500) {
		cout << "invalid indexLen: " << indexLen << endl;
		exit(0);
	}

	if(indexAng < 0 || indexAng >= 450) {
		cout << "invalid index ang: " << indexAng << endl;
		exit(0);
	}

	int index = indexLen*1800 + indexAng*4;

	outIndex[0] = this->biIndex[index];
	outIndex[1] = this->biIndex[index+1];
	outIndex[2] = this->biIndex[index+2];
	outIndex[3] = this->biIndex[index+3];

	outIndex[4] = spIndexA*this->spherePointNum+spIndexB;
	outWeights[0] = this->biWt[index];
	outWeights[1] = this->biWt[index+1];
	outWeights[2] = this->biWt[index+2];
	outWeights[3] = this->biWt[index+3];

	/*
	printf("len interpolation: %6.3f %7.4f %6.3f %6.3f\n", len, indexLen*0.03+0.015, (outIndex[0]/45)*0.3, (outIndex[2]/45)*0.3);
	printf("ang interpolation: %6.2f %7.3f %6.2f %6.2f\n", ang, indexAng*0.8+0.4, (outIndex[0]%45)*8.0, (outIndex[1]%45)*8.0);
	printf("wts: %6.4f %6.4f %6.4f %6.4f\n", outWeights[0], outWeights[1], outWeights[2], outWeights[3]);
	*/
}

void CsMoveTo6DKey::getIndexAndWeight6DInterpolation(const LocalFrame& csA, const LocalFrame& csB, double len, int outIndex[13], double outWeights[13]){


	if(len == 0){
		int index = 0;
		int idA = 0;
		int idB = 0;

		outIndex[0] = this->biIndex[index];
		outIndex[1] = this->biIndex[index+1];
		outIndex[2] = this->biIndex[index+2];
		outIndex[3] = this->biIndex[index+3];

		outIndex[4] = spIndex1[idA]*this->spherePointNum+spIndex1[idB];
		outIndex[5] = spIndex1[idA]*this->spherePointNum+spIndex2[idB];
		outIndex[6] = spIndex1[idA]*this->spherePointNum+spIndex3[idB];
		outIndex[7] = spIndex2[idA]*this->spherePointNum+spIndex1[idB];
		outIndex[8] = spIndex2[idA]*this->spherePointNum+spIndex2[idB];
		outIndex[9] = spIndex2[idA]*this->spherePointNum+spIndex3[idB];
		outIndex[10] = spIndex3[idA]*this->spherePointNum+spIndex1[idB];
		outIndex[11] = spIndex3[idA]*this->spherePointNum+spIndex2[idB];
		outIndex[12] = spIndex3[idA]*this->spherePointNum+spIndex3[idB];


		outWeights[0] = this->biWt[index];
		outWeights[1] = this->biWt[index+1];
		outWeights[2] = this->biWt[index+2];
		outWeights[3] = this->biWt[index+3];
		outWeights[4] = spWt1[idA]*spWt1[idB];
		outWeights[5] = spWt1[idA]*spWt2[idB];
		outWeights[6] = spWt1[idA]*spWt3[idB];
		outWeights[7] = spWt2[idA]*spWt1[idB];
		outWeights[8] = spWt2[idA]*spWt2[idB];
		outWeights[9] = spWt2[idA]*spWt3[idB];
		outWeights[10] = spWt3[idA]*spWt1[idB];
		outWeights[11] = spWt3[idA]*spWt2[idB];
		outWeights[12] = spWt3[idA]*spWt3[idB];

		return;
	}

	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	int idX1 = (int)((tBInA.x_ + 1.0)*200);
	if(idX1 == 400) idX1 = 399;
	int idY1 = (int)((tBInA.y_ + 1.0)*200);
	if(idY1 == 400) idY1 = 399;
	int idZ1 = (int)((tBInA.z_ + 1.0)*200);
	if(idZ1 == 400) idZ1 = 399;


	int idX2 = (int)((tAInB.x_ + 1.0)*200);
	if(idX2 == 400) idX2 = 399;
	int idY2 = (int)((tAInB.y_ + 1.0)*200);
	if(idY2 == 400) idY2 = 399;
	int idZ2 = (int)((tAInB.z_ + 1.0)*200);
	if(idZ2 == 400) idZ2 = 399;

	it = keyIdToListID.find(idX1*160000+idY1*400+idZ1);
	if(it == keyIdToListID.end()){
		cout << "len: " << len << endl;
		cout << "B in A: " << tBInA.toString() << endl;
		cout << "invalid sphere index: " << idX1 << " " << idY1 << " " << idZ1 << endl;
		exit(1);
	}
	int idA = it->second;

	it = keyIdToListID.find(idX2*160000+idY2*400+idZ2);
	if(it == keyIdToListID.end()){
		cout << "len: " << len << endl;
		cout << "A in B: " << tAInB.toString() << endl;
		cout << "invalid sphere index: " << idX2 << " " << idY2 << " " << idZ2 << endl;
		exit(1);
	}
	int idB = it->second;

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;

	int indexLen = (int)(len*33.3333);
	int indexAng = (int)(ang*1.25);
	if(indexAng == 450) indexAng = 0;

	if(indexLen <0 || indexLen >= 500) {
		cout << "invalid indexLen: " << indexLen << endl;
		exit(0);
	}

	if(indexAng < 0 || indexAng >= 450) {
		cout << "invalid index ang: " << indexAng << endl;
		exit(0);
	}

	int index = indexLen*1800 + indexAng*4;

	outIndex[0] = this->biIndex[index];
	outIndex[1] = this->biIndex[index+1];
	outIndex[2] = this->biIndex[index+2];
	outIndex[3] = this->biIndex[index+3];


	outIndex[4] = spIndex1[idA]*this->spherePointNum+spIndex1[idB];
	outIndex[5] = spIndex1[idA]*this->spherePointNum+spIndex2[idB];
	outIndex[6] = spIndex1[idA]*this->spherePointNum+spIndex3[idB];
	outIndex[7] = spIndex2[idA]*this->spherePointNum+spIndex1[idB];
	outIndex[8] = spIndex2[idA]*this->spherePointNum+spIndex2[idB];
	outIndex[9] = spIndex2[idA]*this->spherePointNum+spIndex3[idB];
	outIndex[10] = spIndex3[idA]*this->spherePointNum+spIndex1[idB];
	outIndex[11] = spIndex3[idA]*this->spherePointNum+spIndex2[idB];
	outIndex[12] = spIndex3[idA]*this->spherePointNum+spIndex3[idB];


	outWeights[0] = this->biWt[index];
	outWeights[1] = this->biWt[index+1];
	outWeights[2] = this->biWt[index+2];
	outWeights[3] = this->biWt[index+3];
	outWeights[4] = spWt1[idA]*spWt1[idB];
	outWeights[5] = spWt1[idA]*spWt2[idB];
	outWeights[6] = spWt1[idA]*spWt3[idB];
	outWeights[7] = spWt2[idA]*spWt1[idB];
	outWeights[8] = spWt2[idA]*spWt2[idB];
	outWeights[9] = spWt2[idA]*spWt3[idB];
	outWeights[10] = spWt3[idA]*spWt1[idB];
	outWeights[11] = spWt3[idA]*spWt2[idB];
	outWeights[12] = spWt3[idA]*spWt3[idB];

	/*
	printf("len interpolation: %6.3f %7.4f %6.3f %6.3f\n", len, indexLen*0.03+0.015, (outIndex[0]/45)*0.3, (outIndex[2]/45)*0.3);
	printf("ang interpolation: %6.2f %7.3f %6.2f %6.2f\n", ang, indexAng*0.8+0.4, (outIndex[0]%45)*8.0, (outIndex[1]%45)*8.0);
	printf("wts: %6.4f %6.4f %6.4f %6.4f\n", outWeights[0], outWeights[1], outWeights[2], outWeights[3]);
	*/
}

pair<int,int> CsMoveTo6DKey::toIndexPair(const LocalFrame& csA, const LocalFrame& csB, double len){
	if(len == 0){
		pair<int,int> p(0, 0);
		return p;
	}

	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	char ss[7];
	ss[6] = '\0';

	int idX1 = (int)((tBInA.x_ + 1.0)*200);
	if(idX1 == 400) idX1 = 399;
	int idY1 = (int)((tBInA.y_ + 1.0)*200);
	if(idY1 == 400) idY1 = 399;
	int idZ1 = (int)((tBInA.z_ + 1.0)*200);
	if(idZ1 == 400) idZ1 = 399;

	int idX2 = (int)((tAInB.x_ + 1.0)*200);
	if(idX2 == 400) idX2 = 399;
	int idY2 = (int)((tAInB.y_ + 1.0)*200);
	if(idY2 == 400) idY2 = 399;
	int idZ2 = (int)((tAInB.z_ + 1.0)*200);
	if(idZ2 == 400) idZ2 = 399;

	int spIndexA = sphereKeyMap[idX1*160000+idY1*400+idZ1];
	int spIndexB = sphereKeyMap[idX2*160000+idY2*400+idZ2];

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;


	int indexAng = (int)(ang*0.125 + 0.5);
	if(indexAng == 45) indexAng = 0;

	int indexLen = (int)(len*3.33333333333+0.5);
	if(indexLen == 50) indexLen = 49;
	pair<int,int> p(indexLen*45+indexAng, spIndexA*this->spherePointNum+spIndexB);
	return p;
}

/*
string CsMoveTo6DKey::toKey(const LocalFrame& csA, const LocalFrame& csB, double len){
	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	char ss[7];
	ss[6] = '\0';

	int idX1 = (int)((tBInA.x_ + 1.0)*100);
	if(idX1 == 200) idX1 = 199;
	int idY1 = (int)((tBInA.y_ + 1.0)*100);
	if(idY1 == 200) idY1 = 199;
	int idZ1 = (int)((tBInA.z_ + 1.0)*100);
	if(idZ1 == 200) idZ1 = 199;

	int idX2 = (int)((tAInB.x_ + 1.0)*100);
	if(idX2 == 200) idX2 = 199;
	int idY2 = (int)((tAInB.y_ + 1.0)*100);
	if(idY2 == 200) idY2 = 199;
	int idZ2 = (int)((tAInB.z_ + 1.0)*100);
	if(idZ2 == 200) idZ2 = 199;

	int spIndexA = sphereKeyMap[idX1*40000+idY1*200+idZ1];
	int spIndexB = sphereKeyMap[idX2*40000+idY2*200+idZ2];

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;


	int indexAng = (int)(ang*0.125 + 0.5);
	if(indexAng == 45) indexAng = 0;

	ss[0] = spIndexA/30 + '!';
	ss[1] = spIndexA%30 + '!';
	ss[2] = spIndexB/30 + '!';
	ss[3] = spIndexB%30 + '!';
	ss[4] = (int)(len*3.33333333333+0.5) + '!';
	ss[5] = indexAng + '!';
	return string(ss);
}

string CsMoveTo6DKey::toKey2(const LocalFrame& csA, const LocalFrame& csB, double len){
	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	char ss[7];
	ss[6] = '\0';

	int idX1 = (int)((tBInA.x_ + 1.0)*100);
	if(idX1 == 200) idX1 = 199;
	int idY1 = (int)((tBInA.y_ + 1.0)*100);
	if(idY1 == 200) idY1 = 199;
	int idZ1 = (int)((tBInA.z_ + 1.0)*100);
	if(idZ1 == 200) idZ1 = 199;

	int idX2 = (int)((tAInB.x_ + 1.0)*100);
	if(idX2 == 200) idX2 = 199;
	int idY2 = (int)((tAInB.y_ + 1.0)*100);
	if(idY2 == 200) idY2 = 199;
	int idZ2 = (int)((tAInB.z_ + 1.0)*100);
	if(idZ2 == 200) idZ2 = 199;

	int spIndexA = sphereKeyMap2[idX1*40000+idY1*200+idZ1];
	int spIndexB = sphereKeyMap2[idX2*40000+idY2*200+idZ2];

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;


	int indexAng = (int)(ang*0.041666667 + 0.5);
	if(indexAng == 15) indexAng = 0;

	ss[0] = spIndexA/30 + '!';
	ss[1] = spIndexA%30 + '!';
	ss[2] = spIndexB/30 + '!';
	ss[3] = spIndexB%30 + '!';
	ss[4] = (int)(len+0.5) + '!';
	ss[5] = indexAng + '!';
	return string(ss);
}

int CsMoveTo6DKey::toIndex(const LocalFrame& csA, const LocalFrame& csB, double len){
	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	int idX1 = (int)((tBInA.x_ + 1.0)*100);
	if(idX1 == 200) idX1 = 199;
	int idY1 = (int)((tBInA.y_ + 1.0)*100);
	if(idY1 == 200) idY1 = 199;
	int idZ1 = (int)((tBInA.z_ + 1.0)*100);
	if(idZ1 == 200) idZ1 = 199;

	int idX2 = (int)((tAInB.x_ + 1.0)*100);
	if(idX2 == 200) idX2 = 199;
	int idY2 = (int)((tAInB.y_ + 1.0)*100);
	if(idY2 == 200) idY2 = 199;
	int idZ2 = (int)((tAInB.z_ + 1.0)*100);
	if(idZ2 == 200) idZ2 = 199;

	int spIndexA = sphereKeyMap[idX1*40000+idY1*200+idZ1];
	int spIndexB = sphereKeyMap[idX2*40000+idY2*200+idZ2];

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;

	int indexLen = (int)(len*3.333333333+0.5);
	int indexAng = (int)(ang*0.125 + 0.5);
	if(indexAng == 45) indexAng = 0;

	return spIndexA*1125000+spIndexB*2250+indexLen*45+indexAng;
}

int CsMoveTo6DKey::toIndex2(const LocalFrame& csA, const LocalFrame& csB, double len){
	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	int idX1 = (int)((tBInA.x_ + 1.0)*100);
	if(idX1 == 200) idX1 = 199;
	int idY1 = (int)((tBInA.y_ + 1.0)*100);
	if(idY1 == 200) idY1 = 199;
	int idZ1 = (int)((tBInA.z_ + 1.0)*100);
	if(idZ1 == 200) idZ1 = 199;

	int idX2 = (int)((tAInB.x_ + 1.0)*100);
	if(idX2 == 200) idX2 = 199;
	int idY2 = (int)((tAInB.y_ + 1.0)*100);
	if(idY2 == 200) idY2 = 199;
	int idZ2 = (int)((tAInB.z_ + 1.0)*100);
	if(idZ2 == 200) idZ2 = 199;

	int spIndexA = sphereKeyMap2[idX1*40000+idY1*200+idZ1];
	int spIndexB = sphereKeyMap2[idX2*40000+idY2*200+idZ2];

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;

	int indexLen = (int)(len+0.5);
	int indexAng = (int)(ang*0.0416666667 + 0.5);
	if(indexAng == 15) indexAng = 0;

	return spIndexA*22500+spIndexB*225+indexLen*15+indexAng;
}
*/

CsMoveTo6DKey::~CsMoveTo6DKey() {
	// TODO Auto-generated destructor stub
}

BasePair6DEnergyTable::BasePair6DEnergyTable(ForceFieldPara* para, bool withBinary, int binaryMode):
cm2Key(withBinary, binaryMode) 
{
	this->wtNb = para->wtBp1;
	this->wtNnb = para->wtBp2;

	string path = NSPdataio::datapath();

	int indexA, indexB, clusterID;
	double ene;
	ifstream file;
	string augc = "AUGC";


	if(withBinary && binaryMode==2) {
		string fileName = path + "../binaryData/pairEne/nb";
		if(para->bwTag != "default") {
			fileName = path + "../binaryData/pairEne/nb-" + para->bwTag;
		}
		ifstream ins;
        ins.open(fileName,ios::in | ios::binary);
        if(!ins.is_open()) {
            throw("Unable to open " + fileName);
        }
        auto* bb = new BinaryBook;
        bb->read(ins);
        ins.close();
		for(int i=0; i<4; i++) {
			for(int j=0;j<4;j++) {
				string tN = string(augc.substr(i,1)) + augc.substr(j,1) + ".ene";
				BinaryTable* btab = bb->tables_map.at(tN);
				auto& indexACol = get<BinaryColumn<int>>(*(btab->cols[0]));
				auto& indexBCol = get<BinaryColumn<int>>(*(btab->cols[1]));
				auto& eneCol = get<BinaryColumn<double>>(*(btab->cols[2]));
				// nbKeysEnergy 在元素map中随机赋值，且各map size 未知，无法并行
				for(int k=0;k<btab->nRow;k++) {
					// this->nbKeysEnergy[(i*4+j)*2250+indexACol[k]][indexBCol[k]] = eneCol[k];
					this->nbKeysEnergy[(i*4+j)*2250+indexACol[k]].emplace(indexBCol[k],eneCol[k]);
				}
			}
		}
		delete bb;

		fileName = path + "../binaryData/pairEne/nnb";
		if(para->bwTag != "default") {
			fileName = path + "../binaryData/pairEne/nnb-" + para->bwTag;
		}
        ins.open(fileName,ios::in | ios::binary);
        if(!ins.is_open()) {
            throw("Unable to open " + fileName);
        }
        bb = new BinaryBook;
        bb->read(ins);
        ins.close();
		for(int i=0; i<4; i++) {
			for(int j=i;j<4;j++) {
				string tN = string(augc.substr(i,1)) + augc.substr(j,1) + ".ene";
				BinaryTable* btab = bb->tables_map.at(tN);
				auto& indexACol = get<BinaryColumn<int>>(*(btab->cols[0]));
				auto& indexBCol = get<BinaryColumn<int>>(*(btab->cols[1]));
				auto& eneCol = get<BinaryColumn<double>>(*(btab->cols[2]));
				for(int k=0;k<btab->nRow;k++) {
					// this->nnbKeysEnergy[(i*4+j)*2250+indexACol[k]][indexBCol[k]] = eneCol[k];
					this->nnbKeysEnergy[(i*4+j)*2250+indexACol[k]].emplace(indexBCol[k], eneCol[k]);
				}
			}
		}
		delete bb;		
	} else if(withBinary && binaryMode == 1) {
		return;
	} else {
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				string pairType = augc.substr(i,1) + augc.substr(j,1);
				string fileName = path + "pairEne/nb/"+pairType+".ene";
				if(para->bwTag != "default") {
					fileName = path + "pairEne/nb/"+pairType+".ene-" + para->bwTag;
				}
				file.open(fileName.c_str());
				if(!file.is_open()) {
					cout << "can't open file " << fileName << endl;
				}
				while(file >> indexA >> indexB >> ene >> clusterID){
					this->nbKeysEnergy[(i*4+j)*2250+indexA][indexB] = ene;
				}
				file.close();
			}
		}

		for(int i=0;i<4;i++) {
			for(int j=i;j<4;j++) {
				string pairType = augc.substr(i,1) + augc.substr(j,1);
				string fileName = path + "pairEne/nnb/"+pairType+".ene";
				if(para->bwTag != "default") {
					fileName = path + "pairEne/nnb/"+pairType+".ene-" + para->bwTag;
				}
				file.open(fileName.c_str());
				if(!file.is_open()) {
					cout << "can't open file " << fileName << endl;
				}
				while(file >> indexA >> indexB >> ene >> clusterID){
					this->nnbKeysEnergy[(i*4+j)*2250+indexA][indexB] = ene;
				}
				file.close();
			}
		}
	}
}

int BasePair6DEnergyTable::dump(ForceFieldPara* para) {
	// serialized dump method
	string outpath = datapath() + "../binaryCache";
	string fileName = "BasePair6DEnergyTable";
	if(para->bwTag != "default") {
		fileName = "BasePair6DEnergyTable-"+para->bwTag;
	}
	char z0[4] = {'\0'};
	if(! makeDirs(outpath)) {
		throw("[Error]Unable to create " + outpath);
	}
	ofstream outs;
	outs.open(outpath + "/" + fileName, ios::out|ios::binary);
	if(!outs.is_open()) {
		throw("[Error]Fail to open " + outpath + "/" + fileName);
	}
	cm2Key.dump(outs);
	outs.write(reinterpret_cast<char*>(&wtNb), sizeof(double));
	outs.write(reinterpret_cast<char*>(&wtNnb), sizeof(double));
	int nbMapSizes[36000], nnbMapSizes[36000];
	for(int i=0; i<36000; i++) {
		nbMapSizes[i] = nbKeysEnergy[i].size();
		nnbMapSizes[i] = nnbKeysEnergy[i].size();
	}
	outs.write(reinterpret_cast<char*>(nbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nbMapSizes[i]];
		double vals[nbMapSizes[i]];
		int j = 0;
		for(auto& it1:nbKeysEnergy[i]) {
			keys[j] = it1.first;
			vals[j] = it1.second;
			j++;
		}
		outs.write(reinterpret_cast<char*>(keys), sizeof(int)*nbMapSizes[i]);
		if(nbMapSizes[i]%2) {
			outs.write(reinterpret_cast<char*>(z0), sizeof(int));  // align Bytes
		}
		outs.write(reinterpret_cast<char*>(vals), sizeof(double)*nbMapSizes[i]);
	}
	outs.write(reinterpret_cast<char*>(nnbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nnbMapSizes[i]];
		double vals[nnbMapSizes[i]];
		int j = 0;
		for(auto& it1:nnbKeysEnergy[i]) {
			keys[j] = it1.first;
			vals[j] = it1.second;
			j++;
		}
		outs.write(reinterpret_cast<char*>(keys), sizeof(int)*nnbMapSizes[i]);
		if(nnbMapSizes[i]%2) {
			outs.write(reinterpret_cast<char*>(z0), sizeof(int));  // align Bytes
		}
		outs.write(reinterpret_cast<char*>(vals), sizeof(double)*nnbMapSizes[i]);
	}

	return EXIT_SUCCESS;
}

int BasePair6DEnergyTable::load(ForceFieldPara* para) {
	ifstream ins;
	string fileName = datapath() + "../binaryCache/BasePair6DEnergyTable";
	if(para->bwTag != "default") {
		
		fileName = datapath() + "../binaryCache/BasePair6DEnergyTable-"+para->bwTag;
		cout << "load: " << fileName << endl;
	}
	ins.open(fileName, ios::in|ios::binary);
	if(!ins.is_open()) {
		throw("[Error]Fail to open " + fileName);
	}
	int tint[4];
	cm2Key.load(ins);
	ins.read(reinterpret_cast<char*>(&wtNb), sizeof(double));
	ins.read(reinterpret_cast<char*>(&wtNnb), sizeof(double));
	int nbMapSizes[36000], nnbMapSizes[36000];
	ins.read(reinterpret_cast<char*>(nbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nbMapSizes[i]];
		double vals[nbMapSizes[i]];
		if(nbMapSizes[i]%2) {
			int memSize = nbMapSizes[i]+1;
			int tint2[memSize];
			ins.read(reinterpret_cast<char*>(tint2), sizeof(int)*memSize);
			memcpy(keys, tint2, nbMapSizes[i]*sizeof(int));
		} else {
			ins.read(reinterpret_cast<char*>(keys), sizeof(int)*nbMapSizes[i]);
		}
		ins.read(reinterpret_cast<char*>(vals), sizeof(double)*nbMapSizes[i]);
		for(int j=0; j<nbMapSizes[i]; j++) {
			nbKeysEnergy[i].emplace(keys[j], vals[j]);
		}
	}
	ins.read(reinterpret_cast<char*>(nnbMapSizes),sizeof(int)*36000);
	for(int i=0; i<36000; i++) {
		int keys[nnbMapSizes[i]];
		double vals[nnbMapSizes[i]];
		if(nnbMapSizes[i]%2) {
			int memSize = nnbMapSizes[i]+1;
			int tint2[memSize];
			ins.read(reinterpret_cast<char*>(tint2), sizeof(int)*memSize);
			memcpy(keys, tint2, nnbMapSizes[i]*sizeof(int));
		} else {
			ins.read(reinterpret_cast<char*>(keys), sizeof(int)*nnbMapSizes[i]);
		}
		ins.read(reinterpret_cast<char*>(vals), sizeof(double)*nnbMapSizes[i]);
		for(int j=0; j<nnbMapSizes[i]; j++) {
			nnbKeysEnergy[i].emplace(keys[j], vals[j]);
		}
	}
	wtNb = para->wtBp1;
	wtNnb = para->wtBp2;
	return EXIT_SUCCESS;
}


BasePair6DEnergyTable::~BasePair6DEnergyTable() {
	// TODO Auto-generated destructor stub
}


NeighborPairEnergyTable::NeighborPairEnergyTable(ForceFieldPara* para, int typeA, int typeB) {

	this->wtNb = para->wtBp1;
	this->wtNnb = para->wtBp2;
	this->typeA = typeA;
	this->typeB = typeB;

	string path = NSPdataio::datapath();

	int indexA, indexB, clusterID;
	double ene;
	ifstream file;
	string augc = "AUGC";

	int i = typeA;
	int j = typeB;
	string pairType = augc.substr(i,1) + augc.substr(j,1);
	file.open(path + "pairEne/nb/"+pairType+".ene");

	if(!file.is_open()) {
		cout << "can't open file " << path + "pairEne/nb/" +pairType+".ene" << endl;
	}
	else {
		cout << "load " << path + "pairEne/nb/" +pairType+".ene" << endl;
	}
			
	while(file >> indexA >> indexB >> ene >> clusterID){
		this->nbKeysEnergy[indexA][indexB] = ene;
	}
	file.close();
}

MergeThreeNnbEnergyTable::MergeThreeNnbEnergyTable(int typeA, int typeB) {

	string path = NSPdataio::datapath();

	int indexA, indexB, clusterID;
	double ene;
	ifstream file;
	string augc = "AUGC";

	string pairType = augc.substr(typeA,1) + augc.substr(typeB,1);
	string s;
	int hbNum;

	cout << "start read file" << endl;
	int n1=0;
	int n2=0;
	int n3=0;

	file.open(path + "pairEne/nnb/"+pairType+".ene");
	if(!file.is_open()) {
		cout << "can't open file " << path + "pairEne/nnb/" +pairType+".ene" << endl;
	}
	while(file >> indexA >> indexB >> ene >> clusterID){
		this->nnbKeysEnergy1[indexA][indexB] = ene;
		this->nnbKeysCluster1[indexA][indexB] = clusterID;
		n1++;
	}
	file.close();


	file.open("/public/home/pengx/briqx/xtb/sp2000/stackingEne/"+ pairType + ".ene", ios::in);
	if(!file.is_open()) {
		cout << "can't open file " <<  "/public/home/pengx/briqx/xtb/sp2000/stackingEne/"+ pairType + ".ene" << endl;
	}
	while(file >> indexA >> indexB >> ene) {
		this->nnbKeysEnergy2[indexA][indexB] = (ene+4.0)*0.3;
		n2++;
	}
	file.close();


	file.open("/public/home/pengx/briqx/basePair/hbPair/ene/"+ pairType + ".ene", ios::in);
	if(!file.is_open()) {
		cout << "can't open file " <<  "/public/home/pengx/briqx/basePair/hbPair/ene/"+ pairType + ".ene" << endl;
	}
	while(file >> indexA >> indexB >> ene >> hbNum) {
		this->nnbKeysEnergy3[indexA][indexB] = ene;
		n3++;
	}
	file.close();
	
	cout << n1 << endl;
	cout << n2 << endl;
	cout << n3 << endl;

	file.open("/public/home/pengx/briqx/basePair/mergeBasePairEnergy/fixedNnb", ios::in);
	if(!file.is_open()){
		cout << "can't open /public/home/pengx/briqx/basePair/mergeBasePairEnergyfixedNnb" << endl;
	}
	while(file >> s >> clusterID) {
		if(pairType == s) {
			fixedSet.insert(clusterID);
		}
	}
	file.close();
}

void MergeThreeNnbEnergyTable::mergeEnergy(const string& outfile){
	 
	double e1, e2, e3, em;
	int clusterID;
	int spB;
	char xx[200];
	ofstream out;
	out.open(outfile.c_str(), ios::out);
	if(!out.is_open()){
		cout << "fail to open outfile: " << outfile << endl;
	}

	int pointNum = 0;

	 for(int i=0;i<2250;i++){

		set<int> keySet;
		for(it = nnbKeysEnergy1[i].begin();it != nnbKeysEnergy1[i].end();++it){
			keySet.insert(it->first);
		}
		for(it = nnbKeysEnergy2[i].begin();it != nnbKeysEnergy2[i].end();++it){
			keySet.insert(it->first);
		}
		for(it = nnbKeysEnergy3[i].begin();it != nnbKeysEnergy3[i].end();++it){
			keySet.insert(it->first);
		}

		pointNum += keySet.size();

		for(it2 = keySet.begin();it2 != keySet.end();++it2){
			e1 = 0.0;
			e2 = 0.0;
			e3 = 0.0;
			clusterID = -1;
			spB = *it2;

			it = nnbKeysEnergy1[i].find(spB);
			if(it!=nnbKeysEnergy1[i].end()){
				e1 = it->second;
				it3 = nnbKeysCluster1[i].find(spB);
				clusterID = it3->second;
			}

			it = nnbKeysEnergy2[i].find(spB);
			if(it!=nnbKeysEnergy2[i].end()){
				e2 = it->second;
			}

			it = nnbKeysEnergy3[i].find(spB);
			if(it!=nnbKeysEnergy3[i].end()){
				e3 = it->second;
			}

			if(e1 > 0) e1 = 0;
			if(e2 > 0) e2 = 0;
			if(e3 > 0) e3 = 0;

			if(clusterID >=0 && fixedSet.find(clusterID) != fixedSet.end()){
				em = e1;
			}
			else if(e1 < e2 && e1 < e3)
				em = e1;
			else if(e2 < e1 && e2 < e3)
				em = e2;
			else 
				em = e3;

			em = em + 0.1;
			if(em >= 0.0) continue;

			sprintf(xx, "%d %d %8.4f", i, spB, em);
			out << string(xx) << endl;			
		}
	 }

	 cout << "total point num: " << pointNum << endl;
	 out.close();
}

MergeThreeNnbEnergyTable::~MergeThreeNnbEnergyTable(){

}

} /* namespace NSPforcefield */
