
#include "model/AssignProSSAndSasa.h"

namespace NSPmodel {

ResSasaPoints::ResSasaPoints(){

	this->radii = 3.5;
	psdList[0] = XYZ(-1.040, 0.006, 1.345); /* CB */
	psdList[1] = XYZ(-2.332, -0.663, 0.981); /* CG */
	psdList[2] = XYZ(-3.372, -0.657, 2.325); /* CD */
	points[  0] = XYZ(-2.6710,-0.6280, 5.7540);
	points[  1] = XYZ(-1.2450,-3.2320, 3.3700);
	points[  2] = XYZ(-3.3530, 2.5480, 3.7310);
	points[  3] = XYZ(-4.7830, 0.0620, 5.4460);
	points[  4] = XYZ(-4.6190, 2.5790, 1.8550);
	points[  5] = XYZ(-3.7880,-0.1390, 5.7610);
	points[  6] = XYZ(-5.8400,-2.2830, 4.1990);
	points[  7] = XYZ(-6.6910,-0.5180, 3.4270);
	points[  8] = XYZ(-2.4760,-1.6860, 5.5480);
	points[  9] = XYZ(-6.3290,-1.4460, 0.6270);
	points[ 10] = XYZ(-3.6090,-1.1990, 5.7750);
	points[ 11] = XYZ(-5.6290, 1.9950, 1.9740);
	points[ 12] = XYZ(-5.4320,-3.4710, 2.6240);
	points[ 13] = XYZ(-2.2400,-3.7520, 3.5040);
	points[ 14] = XYZ(-6.8040, 0.0290, 2.3440);
	points[ 15] = XYZ(-3.3670,-3.9570, 3.4910);
	points[ 16] = XYZ(-3.3990,-2.1810, 5.4760);
	points[ 17] = XYZ(-5.4040,-1.6740, 4.9870);
	points[ 18] = XYZ(-1.6050,-3.7420, 2.4780);
	points[ 19] = XYZ(-2.6490,-4.0750, 2.5430);
	points[ 20] = XYZ(-6.2230,-0.0030, 4.2470);
	points[ 21] = XYZ(-5.9660, 1.5970, 2.9880);
	points[ 22] = XYZ(-5.2380,-3.2420, 3.7690);
	points[ 23] = XYZ(-5.2530,-0.4700,-0.9300);
	points[ 24] = XYZ(-4.2700,-3.3670,-0.1000);
	points[ 25] = XYZ(-3.5080, 2.8030, 1.8130);
	points[ 26] = XYZ(-3.8100,-2.9670, 4.9180);
	points[ 27] = XYZ(-6.3970, 1.0950, 2.1530);
	points[ 28] = XYZ(-5.6150,-2.0260, 0.0120);
	points[ 29] = XYZ(-1.8740, 0.0610, 5.4060);
	points[ 30] = XYZ(-5.3060, 1.8480, 3.8190);
	points[ 31] = XYZ(-2.9230, 0.4370, 5.6190);
	points[ 32] = XYZ(-6.1130,-2.8130, 2.0330);
	points[ 33] = XYZ(-1.6450,-2.1700, 4.9660);
	points[ 34] = XYZ(-6.5050,-1.6370, 3.5390);
	points[ 35] = XYZ(-6.5630,-2.0660, 2.6080);
	points[ 36] = XYZ(-2.9370, 2.9100, 0.8800);
	points[ 37] = XYZ(-4.4930,-3.6930, 0.9930);
	points[ 38] = XYZ(-4.6650,-1.0010, 5.5590);
	points[ 39] = XYZ(-2.2550, 2.4150, 3.5760);
	points[ 40] = XYZ(-6.0140,-2.5020, 0.9590);
	points[ 41] = XYZ(-4.1760, 0.9060, 5.3520);
	points[ 42] = XYZ(-2.2180, 1.1290, 5.1050);
	points[ 43] = XYZ(-3.2870,-4.0720, 1.5610);
	points[ 44] = XYZ(-2.1030,-4.1080, 1.5560);
	points[ 45] = XYZ(-2.9780,-3.4790, 4.3580);
	points[ 46] = XYZ(-3.9100, 1.9530, 4.5940);
	points[ 47] = XYZ(-6.5990,-1.7970, 1.5950);
	points[ 48] = XYZ(-2.7700, 2.0200, 4.4980);
	points[ 49] = XYZ(-6.3450,-0.3440, 0.5050);
	points[ 50] = XYZ(-4.1920,-3.5410, 4.1310);
	points[ 51] = XYZ(-1.6430,-1.0530, 5.3420);
	points[ 52] = XYZ(-6.1750,-1.1230, 4.3680);
	points[ 53] = XYZ(-5.9970, 0.6180, 0.3930);
	points[ 54] = XYZ(-5.0550, 2.0990, 0.9760);
	points[ 55] = XYZ(-6.7610,-0.7040, 1.4510);
	points[ 56] = XYZ(-0.8650,-1.4760, 4.6250);
	points[ 57] = XYZ(-4.0230, 2.7490, 2.7970);
	points[ 58] = XYZ(-5.9380, 1.0350, 4.0000);
	points[ 59] = XYZ(-4.4390,-2.0460, 5.3550);
	points[ 60] = XYZ(-6.0310,-2.7460, 3.2270);
	points[ 61] = XYZ(-5.3610, 0.4090,-0.4000);
	points[ 62] = XYZ(-5.8020,-0.9000,-0.1800);
	points[ 63] = XYZ(-5.9500, 1.3940, 1.1440);
	points[ 64] = XYZ(-1.8620,-3.1070, 4.3160);
	points[ 65] = XYZ(-6.8460,-1.0500, 2.4990);
	points[ 66] = XYZ(-4.6040,-3.9170, 2.0060);
	points[ 67] = XYZ(-4.8550, 1.4500, 4.6940);
	points[ 68] = XYZ(-5.5820,-0.5450, 5.0370);
	points[ 69] = XYZ(-3.2620, 1.3520, 5.1890);
	points[ 70] = XYZ(-2.9080, 2.7870, 2.7440);
	points[ 71] = XYZ(-5.2530, 1.3820, 0.1910);
	points[ 72] = XYZ(-5.4620, 0.5850, 4.8430);
	points[ 73] = XYZ(-6.5550, 0.4050, 1.3280);
	points[ 74] = XYZ(-0.5900,-3.3240, 2.4420);
	points[ 75] = XYZ(-5.1700,-2.9760, 0.4170);
	points[ 76] = XYZ(-5.4170,-3.3610, 1.4540);
	points[ 77] = XYZ(-2.5980,-2.7310, 5.0360);
	points[ 78] = XYZ(-3.7210,-4.1350, 2.5070);
	points[ 79] = XYZ(-4.9140,-2.7440, 4.6740);
	points[ 80] = XYZ(-5.0850, 2.3500, 2.8450);
	points[ 81] = XYZ(-6.5190, 0.5890, 3.2160);
	points[ 82] = XYZ(-4.5370,-3.8510, 3.1570);
	points[ 83] = XYZ(-0.9730,-0.3590, 4.8560);
	points[ 84] = XYZ(-1.2810, 0.7430, 4.7580);
	points[ 85] = XYZ(-0.9640,-2.4660, 4.1080);
	points[ 86] = XYZ(-4.4110, 2.3650, 3.7520);
	points[ 87] = XYZ(-3.9900, 2.4770, 0.8940);
	points[ 88] = XYZ(-4.8100,-2.5580,-0.6000);
	points[ 89] = XYZ(-4.7200, 0.3280,-1.3700);
	points[ 90] = XYZ(-4.3740,-2.0310,-1.5100);
	points[ 91] = XYZ(-3.0870,-2.8230,-1.6600);
	points[ 92] = XYZ(-4.3650, 2.0060,-0.0100);
	points[ 93] = XYZ(-1.4730,-4.0380, 0.6320);
	points[ 94] = XYZ(-3.9620, 0.0250,-2.0300);
	points[ 95] = XYZ(-3.5820,-3.9070, 0.5760);
	points[ 96] = XYZ(-2.5950,-3.4790,-1.0800);
	points[ 97] = XYZ(-4.6710, 1.2560,-0.7700);
	points[ 98] = XYZ(-4.6290,-0.8860,-1.6500);
	points[ 99] = XYZ(-0.9800,-3.8380, 1.5630);
	points[100] = XYZ(-3.4950,-1.9930,-2.0400);
	points[101] = XYZ(-2.5680,-4.1320, 0.5780);
	points[102] = XYZ(-2.0450,-3.9070,-0.3000);
	points[103] = XYZ(-3.9030,-3.0130,-1.0800);
	points[104] = XYZ(-3.7520, 1.9060,-0.9200);
	points[105] = XYZ(-5.1500,-1.5590,-0.8900);
	points[106] = XYZ(-3.7970,-1.0310,-2.1700);
	points[107] = XYZ(-3.3860, 2.5100,-0.0500);
	points[108] = XYZ(-3.9480, 1.0610,-1.6000);
	points[109] = XYZ(-3.2600,-3.7430,-0.3900);
	points[110] = XYZ(-0.7130, 1.8570, 4.2970);
	points[111] = XYZ(-2.4390, 3.1750, 1.8450);
	points[112] = XYZ(-0.2180,-2.7780, 3.3000);
	points[113] = XYZ(-0.0440,-0.0260, 4.7000);
	points[114] = XYZ(-1.7430, 1.7200, 4.3140);
	points[115] = XYZ(-1.8510, 3.0750, 2.8200);
	points[116] = XYZ(-0.3110, 0.9590, 4.6330);
	points[117] = XYZ( 0.0160,-1.9920, 4.0180);
	points[118] = XYZ(-1.2390, 2.6000, 3.6870);
	points[119] = XYZ( 0.0720,-1.0370, 4.4960);
	sasaIndex[ 0] = 0.994;
	sasaIndex[ 1] = 0.982;
	sasaIndex[ 2] = 0.972;
	sasaIndex[ 3] = 0.964;
	sasaIndex[ 4] = 0.957;
	sasaIndex[ 5] = 0.949;
	sasaIndex[ 6] = 0.942;
	sasaIndex[ 7] = 0.936;
	sasaIndex[ 8] = 0.929;
	sasaIndex[ 9] = 0.923;
	sasaIndex[10] = 0.918;
	sasaIndex[11] = 0.912;
	sasaIndex[12] = 0.907;
	sasaIndex[13] = 0.902;
	sasaIndex[14] = 0.896;
	sasaIndex[15] = 0.891;
	sasaIndex[16] = 0.885;
	sasaIndex[17] = 0.879;
	sasaIndex[18] = 0.871;
	sasaIndex[19] = 0.863;
	sasaIndex[20] = 0.851;
	sasaIndex[21] = 0.839;
	sasaIndex[22] = 0.829;
	sasaIndex[23] = 0.821;
	sasaIndex[24] = 0.814;
	sasaIndex[25] = 0.807;
	sasaIndex[26] = 0.801;
	sasaIndex[27] = 0.794;
	sasaIndex[28] = 0.789;
	sasaIndex[29] = 0.783;
	sasaIndex[30] = 0.777;
	sasaIndex[31] = 0.771;
	sasaIndex[32] = 0.765;
	sasaIndex[33] = 0.760;
	sasaIndex[34] = 0.754;
	sasaIndex[35] = 0.749;
	sasaIndex[36] = 0.744;
	sasaIndex[37] = 0.738;
	sasaIndex[38] = 0.733;
	sasaIndex[39] = 0.728;
	sasaIndex[40] = 0.722;
	sasaIndex[41] = 0.717;
	sasaIndex[42] = 0.712;
	sasaIndex[43] = 0.707;
	sasaIndex[44] = 0.702;
	sasaIndex[45] = 0.698;
	sasaIndex[46] = 0.693;
	sasaIndex[47] = 0.688;
	sasaIndex[48] = 0.683;
	sasaIndex[49] = 0.678;
	sasaIndex[50] = 0.673;
	sasaIndex[51] = 0.668;
	sasaIndex[52] = 0.663;
	sasaIndex[53] = 0.658;
	sasaIndex[54] = 0.654;
	sasaIndex[55] = 0.649;
	sasaIndex[56] = 0.644;
	sasaIndex[57] = 0.639;
	sasaIndex[58] = 0.634;
	sasaIndex[59] = 0.629;
	sasaIndex[60] = 0.624;
	sasaIndex[61] = 0.619;
	sasaIndex[62] = 0.615;
	sasaIndex[63] = 0.610;
	sasaIndex[64] = 0.604;
	sasaIndex[65] = 0.599;
	sasaIndex[66] = 0.594;
	sasaIndex[67] = 0.589;
	sasaIndex[68] = 0.584;
	sasaIndex[69] = 0.579;
	sasaIndex[70] = 0.573;
	sasaIndex[71] = 0.568;
	sasaIndex[72] = 0.562;
	sasaIndex[73] = 0.556;
	sasaIndex[74] = 0.551;
	sasaIndex[75] = 0.545;
	sasaIndex[76] = 0.539;
	sasaIndex[77] = 0.532;
	sasaIndex[78] = 0.526;
	sasaIndex[79] = 0.519;
	sasaIndex[80] = 0.512;
	sasaIndex[81] = 0.505;
	sasaIndex[82] = 0.498;
	sasaIndex[83] = 0.490;
	sasaIndex[84] = 0.482;
	sasaIndex[85] = 0.474;
	sasaIndex[86] = 0.465;
	sasaIndex[87] = 0.455;
	sasaIndex[88] = 0.446;
	sasaIndex[89] = 0.435;
	sasaIndex[90] = 0.423;
	sasaIndex[91] = 0.410;
	sasaIndex[92] = 0.396;
	sasaIndex[93] = 0.380;
	sasaIndex[94] = 0.361;
	sasaIndex[95] = 0.339;
	sasaIndex[96] = 0.313;
	sasaIndex[97] = 0.278;
	sasaIndex[98] = 0.238;
	sasaIndex[99] = 0.109;
	wts[  0] = 0.002140;
	wts[  1] = 0.000319;
	wts[  2] = 0.014558;
	wts[  3] = 0.016743;
	wts[  4] = 0.006264;
	wts[  5] = 0.018719;
	wts[  6] = 0.023217;
	wts[  7] = 0.017526;
	wts[  8] = 0.004407;
	wts[  9] = 0.014522;
	wts[ 10] = 0.001675;
	wts[ 11] = 0.018883;
	wts[ 12] = 0.003624;
	wts[ 13] = 0.012337;
	wts[ 14] = 0.005408;
	wts[ 15] = 0.017818;
	wts[ 16] = 0.002513;
	wts[ 17] = 0.023317;
	wts[ 18] = 0.007366;
	wts[ 19] = 0.011690;
	wts[ 20] = 0.017763;
	wts[ 21] = 0.016434;
	wts[ 22] = 0.016261;
	wts[ 23] = 0.013921;
	wts[ 24] = 0.008467;
	wts[ 25] = 0.006601;
	wts[ 26] = 0.019311;
	wts[ 27] = 0.017690;
	wts[ 28] = 0.010734;
	wts[ 29] = 0.000164;
	wts[ 30] = 0.022680;
	wts[ 31] = 0.004325;
	wts[ 32] = 0.011190;
	wts[ 33] = 0.020003;
	wts[ 34] = 0.006492;
	wts[ 35] = 0.000546;
	wts[ 36] = 0.006847;
	wts[ 37] = 0.006892;
	wts[ 38] = 0.005490;
	wts[ 39] = 0.001429;
	wts[ 40] = 0.009615;
	wts[ 41] = 0.006910;
	wts[ 42] = 0.003760;
	wts[ 43] = 0.007256;
	wts[ 44] = 0.007766;
	wts[ 45] = 0.002832;
	wts[ 46] = 0.003023;
	wts[ 47] = 0.004662;
	wts[ 48] = 0.002495;
	wts[ 49] = 0.006783;
	wts[ 50] = 0.002249;
	wts[ 51] = 0.023390;
	wts[ 52] = 0.023062;
	wts[ 53] = 0.004771;
	wts[ 54] = 0.018218;
	wts[ 55] = 0.007639;
	wts[ 56] = 0.000127;
	wts[ 57] = 0.011144;
	wts[ 58] = 0.023681;
	wts[ 59] = 0.023253;
	wts[ 60] = 0.008558;
	wts[ 61] = 0.006619;
	wts[ 62] = 0.006920;
	wts[ 63] = 0.003150;
	wts[ 64] = 0.024537;
	wts[ 65] = 0.006501;
	wts[ 66] = 0.005709;
	wts[ 67] = 0.004452;
	wts[ 68] = 0.007429;
	wts[ 69] = 0.019575;
	wts[ 70] = 0.007083;
	wts[ 71] = 0.005445;
	wts[ 72] = 0.007976;
	wts[ 73] = 0.005517;
	wts[ 74] = 0.000346;
	wts[ 75] = 0.007994;
	wts[ 76] = 0.010598;
	wts[ 77] = 0.001785;
	wts[ 78] = 0.006055;
	wts[ 79] = 0.005945;
	wts[ 80] = 0.023454;
	wts[ 81] = 0.007821;
	wts[ 82] = 0.003341;
	wts[ 83] = 0.000737;
	wts[ 84] = 0.000765;
	wts[ 85] = 0.000319;
	wts[ 86] = 0.019611;
	wts[ 87] = 0.007238;
	wts[ 88] = 0.011581;
	wts[ 89] = 0.003997;
	wts[ 90] = 0.003187;
	wts[ 91] = 0.000856;
	wts[ 92] = 0.006747;
	wts[ 93] = 0.005554;
	wts[ 94] = 0.000656;
	wts[ 95] = 0.007666;
	wts[ 96] = 0.000419;
	wts[ 97] = 0.004807;
	wts[ 98] = 0.001002;
	wts[ 99] = 0.000692;
	wts[100] = 0.001329;
	wts[101] = 0.008194;
	wts[102] = 0.004134;
	wts[103] = 0.003815;
	wts[104] = 0.018073;
	wts[105] = 0.013001;
	wts[106] = 0.000328;
	wts[107] = 0.012091;
	wts[108] = 0.001211;
	wts[109] = 0.001165;
	wts[110] = 0.000373;
	wts[111] = 0.006865;
	wts[112] = 0.000519;
	wts[113] = 0.017763;
	wts[114] = 0.001511;
	wts[115] = 0.001074;
	wts[116] = 0.002486;
	wts[117] = 0.001420;
	wts[118] = 0.000920;
	wts[119] = 0.011162;
}

float ResSasaPoints::getSAI(vector<XYZ>& localTerns){
	int n = localTerns.size();
	bool bury[120];
	float score = 0;
	for(int i=0;i<120;i++){
		bury[i] = false;
	}
	for(int i=0;i<120;i++){
		XYZ ball = points[i];
		int j = 0;
		while(!bury[i] && j<n){
			if(localTerns[j].distance(ball) < this->radii)
				bury[i] = true;
			j++;
		}
		if(bury[i])
			score += wts[i];
	}
	score = score * 100;
	return sasaIndex[(int)score];
}

float ResSasaPoints::getSAI(Residue* resA, vector<Residue*>& resList){
	vector<XYZ> ternList;
	LocalFrame cs = resA->getCoordSystem();
	XYZ CB = cs.local2globalcrd(psdList[0]);
	XYZ CG = cs.local2globalcrd(psdList[1]);
	XYZ CD = cs.local2globalcrd(psdList[2]);

	double cutoff = 7.0;
	double d0, d1, d2, d3;
	int N = resList.size();
	int n;
	Residue* res;
	Atom* a;
	for(int i=0;i<N;i++){
		res = resList[i];
		double d0 = res->getCbCoord().distance(resA->getCbCoord());
		if(d0 < 0.0001 || d0 > 18.0) continue;
		vector<Atom*>* bbAtomList = res->getBackboneAtoms();
		for(int j=0;j<bbAtomList->size();j++){
			a = bbAtomList->at(j);
			XYZ t = a->getCoord();
			d1 = t.distance(CB);
			d2 = t.distance(CG);
			d3 = t.distance(CD);
			if(d1 > cutoff && d2 > cutoff && d3 > cutoff)
				continue;
			ternList.push_back(cs.global2localcrd(t));
		}
		if(res->intName == 5) n = 0;
		else if(res->intName == 0 || res->intName == 12) n = 1;
		else if(res->intName == 1 || res->intName == 15 || res->intName == 16 || res->intName == 17) n = 2;
		else n = 3;

		LocalFrame cs2 = res->getCoordSystem();
		for(int k=0;k<n;k++){
			XYZ t = psdList[k];
			XYZ tGlobal = cs2.local2globalcrd(t);
			d1 = tGlobal.distance(CB);
			d2 = tGlobal.distance(CG);
			d3 = tGlobal.distance(CD);
			if(d1 > cutoff && d2 > cutoff && d3 > cutoff) continue;
			ternList.push_back(cs.global2localcrd(tGlobal));
		}
	}

	vector<Atom*>* bbAtomList = resA->getBackboneAtoms();

	for(int j=0;j<bbAtomList->size();j++){
		a = bbAtomList->at(j);
		XYZ t = a->getCoord();
		d1 = t.distance(CB);
		d2 = t.distance(CG);
		d3 = t.distance(CD);
		if(d1 > cutoff && d2 > cutoff && d3 > cutoff) continue;
		ternList.push_back(cs.global2localcrd(t));
	}
	return getSAI(ternList);
}


ResSasaPoints::~ResSasaPoints(){

}

BackboneHBond::BackboneHBond(XYZ& N, XYZ& H, XYZ& O, Residue* resA, Residue* resB) {
	this->donerChainID = resA->chainID;
	this->acceptorChainID = resB->chainID;
	this->donerResID = resA->resID;
	this->acceptorResID = resB->resID;
	if(donerChainID != acceptorChainID)
		this->seqSeparation = 999;
	else
		this->seqSeparation = abs(resA->resSeqID - resB->resSeqID);

	this->hDoner = N;
	this->hydrogen = H;
	this->hAcceptor = O;
	this->distance = hydrogen.distance(hAcceptor);
	this->angle = angleX(N,H,O);
}

string BackboneHBond::toString()
{
	char s[100];
	sprintf(s,"%-4d %c%-4s %c%-4s %5.3f %6.2f", this->seqSeparation, donerChainID, donerResID.c_str(), acceptorChainID, acceptorResID.c_str(), distance, angle);
	return s;
}

BackboneHBond::~BackboneHBond(){

}


AssignProSSAndSasa::AssignProSSAndSasa(PDB* protein) {
	vector<Residue*> tmpResList = protein->getResList();
	int index = 0;
	char s[10];
	//cout << "tmpList " << tmpResList.size() << endl;
	for(unsigned int i=0;i<tmpResList.size();i++)
	{
		Residue* res = tmpResList[i];
		if(res->hasThreeCoreAtoms() && res->hasAtom("O"))
		{
			resList.push_back(res);
			NList.push_back(&res->getAtom("N")->getCoord());
			CAList.push_back(&res->getAtom("CA")->getCoord());
			CList.push_back(&res->getAtom("C")->getCoord());
			sprintf(s,"%c%s",res->chainID,res->resID.c_str());
			chainIDResIDToSeqID[string(s)] = index;
			index++;
		}
	}

	this->resNum = resList.size();
	this->ssSeq = new char[resNum+1];
	for(int i=0;i<resNum;i++)
	{
		ssSeq[i] = 'C';
	}
	ssSeq[resNum] = '\0';
}

AssignProSSAndSasa::AssignProSSAndSasa(ProteinChain* pc){
	vector<Residue*> tmpResList = pc->getResList();
	int index = 0;
	char s[10];
	//cout << "tmpList " << tmpResList.size() << endl;
	for(unsigned int i=0;i<tmpResList.size();i++)
	{
		Residue* res = tmpResList[i];
		if(res->hasThreeCoreAtoms() && res->hasAtom("O"))
		{
			resList.push_back(res);
			NList.push_back(&res->getAtom("N")->getCoord());
			CAList.push_back(&res->getAtom("CA")->getCoord());
			CList.push_back(&res->getAtom("C")->getCoord());
			sprintf(s,"%c%s",res->chainID,res->resID.c_str());
			chainIDResIDToSeqID[string(s)] = index;
			index++;
		}
	}

	this->resNum = resList.size();
	this->ssSeq = new char[resNum+1];
	for(int i=0;i<resNum;i++)
	{
		ssSeq[i] = 'C';
	}
	ssSeq[resNum] = '\0';
}

AssignProSSAndSasa::AssignProSSAndSasa(vector<Residue*>& resList){

	int index = 0;
	char s[10];
	//cout << "tmpList " << tmpResList.size() << endl;
	for(unsigned int i=0;i<resList.size();i++)
	{
		Residue* res = resList[i];
		if(res->hasThreeCoreAtoms() && res->hasAtom("O"))
		{
			this->resList.push_back(res);
			NList.push_back(&res->getAtom("N")->getCoord());
			CAList.push_back(&res->getAtom("CA")->getCoord());
			CList.push_back(&res->getAtom("C")->getCoord());
			sprintf(s,"%c%s",res->chainID,res->resID.c_str());
			chainIDResIDToSeqID[string(s)] = index;
			index++;
		}
	}

	this->resNum = resList.size();
	this->ssSeq = new char[resNum+1];
	for(int i=0;i<resNum;i++)
	{
		ssSeq[i] = 'C';
	}
	ssSeq[resNum] = '\0';
}

void AssignProSSAndSasa::updateBBHbonds(){
	if(this->hbList.size() != 0)
		return;

	this->backboneHBonded.clear();
	for(int i=0;i<resNum;i++){
		this->backboneHBonded.push_back(false);
	}

	vector<XYZ*> HNList;
	vector<XYZ*> OList;

	for(int i=0;i<resNum;i++){
		Residue* res = resList.at(i);
		Atom* O = res->getAtom("O");
		if(O != NULL)
			OList.push_back(&(O->getCoord()));
		else
			OList.push_back(NULL);

		if(i==0 || res->triName == "PRO")
			HNList.push_back(NULL);
		else{
			Atom* aC1 = resList.at(i-1)->getAtom("C");
			Atom* aC2 = resList.at(i)->getAtom("CA");
			if(aC1->getCoord().distance(aC2->getCoord()) > 3.0){
				HNList.push_back(NULL);
			}
			else{
				XYZ NC1 = ~(aC1->getCoord() - *NList[i]);
				XYZ NC2 = ~(aC2->getCoord() - *NList[i]);
				XYZ NH = ~(NC1 + NC2);
				XYZ* H = new XYZ(*NList[i] - NH); //need delete
				HNList.push_back(H);
			}
		}
	}

	int i,j;
	XYZ *NI, *NJ, *OI, *OJ, *HI, *HJ;
	float d1,d2;
	for(i=1;i<resNum;i++){
		NI = NList.at(i);
		OI = OList.at(i);
		HI = HNList.at(i);

		for(j=i+2;j<resNum;j++){
			NJ = NList.at(j);
			OJ = OList.at(j);
			HJ = HNList.at(j);

			if(HI != NULL && OJ != NULL){
				d1 = HI->distance(*OJ);
				if(d1 < 2.5 && d1 > 1.6){
					BackboneHBond* hb = new BackboneHBond(*NI,*HI, *OJ, resList.at(i), resList.at(j));
					if(hb->angle > 100){
						hbList.push_back(hb);
						this->backboneHBonded[i] = true;
						this->backboneHBonded[j] = true;
					}
					else
						delete hb;
				}
			}
			if(HJ != NULL && OI != NULL){
				d2 = HJ->distance(*OI);
				if(d2 < 2.5 && d2 > 1.6){
					BackboneHBond* hb = new BackboneHBond(*NJ, *HJ, *OI, resList.at(j), resList.at(i));
					if(hb->angle > 100){
						hbList.push_back(hb);
						this->backboneHBonded[i] = true;
						this->backboneHBonded[j] = true;
					}
					else
						delete hb;
				}
			}
		}
	}
	XYZ* t;
	for(i=0;i<HNList.size();i++)
	{
		t = HNList[i];
		if(t != NULL)
			delete t;
	}
}

void AssignProSSAndSasa::updateSS(){
	updateBBHbonds();

	int NHbList[resNum];
	int OHbList[resNum];
	for(int i=0;i<resNum;i++){
		NHbList[i] = 0;
		OHbList[i] = 0;
	}

	int resASeqID, resBSeqID;
	char s1[10], s2[10];

	for(unsigned int i=0;i<hbList.size();i++){
		BackboneHBond *hb = hbList.at(i);
		sprintf(s1,"%c%s",hb->donerChainID, hb->donerResID.c_str());
		sprintf(s2, "%c%s", hb->acceptorChainID, hb->acceptorResID.c_str());
		resASeqID = chainIDResIDToSeqID[s1];
		resBSeqID = chainIDResIDToSeqID[s2];
		if(NHbList[resASeqID] != 0){
			int oldResBSeqID = resASeqID + NHbList[resASeqID];
			if(abs(resASeqID - resBSeqID) > abs(resASeqID - oldResBSeqID))
				continue;
		}

		if(OHbList[resBSeqID] != 0){
			int oldSeqASeqID =	resBSeqID + OHbList[resBSeqID];
			if(abs(resASeqID - resBSeqID) > abs(oldSeqASeqID - resBSeqID))
				continue;
		}

		NHbList[resASeqID] = resBSeqID - resASeqID;
		OHbList[resBSeqID] = resASeqID - resBSeqID;
	}

	/*
	 * assign Helix
	 */
	for(int i=1;i<resNum-1;i++)
	{
		if(OHbList[i+1] == 4 && OHbList[i] == 4)
			this->ssSeq[i] = 'H';
		else if(NHbList[i] == -4 && NHbList[i+1] != -4)
			this->ssSeq[i] = 'C';
		else if(this->ssSeq[i-1] == 'H')
			this->ssSeq[i] = 'H';
		else if(OHbList[i-1] == 3 && OHbList[i] == 3)
			this->ssSeq[i] = 'G';
		else if(NHbList[i] == -3 && NHbList[i+1] != -3)
			this->ssSeq[i] = 'C';
		else if(this->ssSeq[i-1] == 'G')
			this->ssSeq[i] = 'G';
	}

	/*
	 * assign beta
	 */
	for(int i=1;i<resNum-1;i++)
	{
		if(NHbList[i] == OHbList[i] && abs(NHbList[i]) > 3)
			this->ssSeq[i] = 'E';
		if(NHbList[i]+2 == OHbList[i] && abs(NHbList[i]) > 3)
			this->ssSeq[i] = 'E';
	}

	for(int i=1;i<resNum-1;i++)
	{
		if(this->ssSeq[i] == 'C' && (this->ssSeq[i-1] == 'E' || this->ssSeq[i+1] == 'E'))
			this->ssSeq[i] = 'e';
	}

	for(int i=1;i<resNum-1;i++)
	{
		if(this->ssSeq[i] == 'e')
			this->ssSeq[i] = 'E';
	}
}

string AssignProSSAndSasa::getSS(){
	return string(this->ssSeq);
}

float AssignProSSAndSasa::getSASA(Residue* res){
	char s[10];
	sprintf(s,"%c%s",res->chainID,res->resID.c_str());
	map<string,int>::const_iterator it;
	it = chainIDResIDToSeqID.find(string(s));
	if(it != chainIDResIDToSeqID.end())
	{
		int id = it->second;
		return this->saiList.at(id);
	}
	else
	{
		cerr << "can't find sai of Residue: " << s << endl;
		return 0.5;
	}
}

char AssignProSSAndSasa::getSS(int seqID){
	if(seqID < 0 || seqID >= resNum){
		cout << "out of index: " << seqID << endl;
		exit(1);
	}
	char c = ssSeq[seqID];
	if(c == 'H' || c == 'G')
		return 'H';
	else if(c == 'C')
		return 'C';
	else
		return 'E';
}

float AssignProSSAndSasa::getSASA(int seqID){
	if(seqID < 0 || seqID >= resNum){
		cout << "out of index: " << seqID << endl;
		exit(1);
	}
	return saiList.at(seqID);
}

void AssignProSSAndSasa::updateSasa(ResSasaPoints* rsp){
	for(int i=0;i<this->resNum;i++){
		float sai = rsp->getSAI(resList[i], resList);
		this->saiList.push_back(sai);
	}
}

AssignProSSAndSasa::~AssignProSSAndSasa() {

	delete ssSeq;
	for(int i=0;i<hbList.size();i++){
		delete hbList[i];
	}
}

} /* namespace NSPmath */
