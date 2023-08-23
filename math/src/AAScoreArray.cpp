/*
 * AAScoreArray.cpp
 *
 */

#include "math/AAScoreArray.h"

namespace NSPmath {

AAScoreArray::AAScoreArray() {
	for(int i=0;i<20;i++){
		sa[i] = 0.0;
	}
	saBg[0] = 2.463;
	saBg[1] = 4.423;
	saBg[2] = 2.829;
	saBg[3] = 2.696;
	saBg[4] = 3.209;
	saBg[5] = 2.609;
	saBg[6] = 3.742;
	saBg[7] = 2.851;
	saBg[8] = 2.895;
	saBg[9] = 2.356;
	saBg[10] = 3.990;
	saBg[11] = 3.168;
	saBg[12] = 3.081;
	saBg[13] = 3.300;
	saBg[14] = 2.955;
	saBg[15] = 2.823;
	saBg[16] = 2.904;
	saBg[17] = 2.631;
	saBg[18] = 4.276;
	saBg[19] = 3.352;
}

AAScoreArray::AAScoreArray(const string& line){
	vector<string> spt;
	NSPtools::splitString(line, " ", &spt);
	if(spt.size() != 20){
		cout << "invalid data line: " << line << endl;
		exit(1);
	}
	for(int i=0;i<20;i++){
		this->sa[i] = atof(spt[i].c_str());
	}
	saBg[0] = 2.463;
	saBg[1] = 4.423;
	saBg[2] = 2.829;
	saBg[3] = 2.696;
	saBg[4] = 3.209;
	saBg[5] = 2.609;
	saBg[6] = 3.742;
	saBg[7] = 2.851;
	saBg[8] = 2.895;
	saBg[9] = 2.356;
	saBg[10] = 3.990;
	saBg[11] = 3.168;
	saBg[12] = 3.081;
	saBg[13] = 3.300;
	saBg[14] = 2.955;
	saBg[15] = 2.823;
	saBg[16] = 2.904;
	saBg[17] = 2.631;
	saBg[18] = 4.276;
	saBg[19] = 3.352;
}

AAScoreArray::AAScoreArray(AAProbabilityArray* pa){
	for(int i=0;i<20;i++){
		double p = pa->pa[i];
		if(p <= 0) p = 0.00000001;
		sa[i] = -log(p);
		if(sa[i] > 19.9)
			sa[i] = 19.9;
	}
	saBg[0] = 2.463;
	saBg[1] = 4.423;
	saBg[2] = 2.829;
	saBg[3] = 2.696;
	saBg[4] = 3.209;
	saBg[5] = 2.609;
	saBg[6] = 3.742;
	saBg[7] = 2.851;
	saBg[8] = 2.895;
	saBg[9] = 2.356;
	saBg[10] = 3.990;
	saBg[11] = 3.168;
	saBg[12] = 3.081;
	saBg[13] = 3.300;
	saBg[14] = 2.955;
	saBg[15] = 2.823;
	saBg[16] = 2.904;
	saBg[17] = 2.631;
	saBg[18] = 4.276;
	saBg[19] = 3.352;
}

void AAScoreArray::multipy(float wt){
	for(int i=0;i<20;i++){
		this->sa[i] = this->sa[i]*wt;
	}
}

AAScoreArray::~AAScoreArray() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmath */
