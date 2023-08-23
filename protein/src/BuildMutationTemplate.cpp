/*
 * BuildMutationTemplate.cpp
 *
 */

#include "protein/BuildMutationTemplate.h"

namespace NSPprotein {

BuildMutationTemplate::BuildMutationTemplate(vector<string>& lines, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, ResScRotamerLibMini* scLibSimp, AtomLib* atLib, ResName* rn) {
	// TODO Auto-generated constructor stub
	this->rn = rn;
	this->bbLib = bbLib;
	this->scLib = scLib;
	this->atLib = atLib;
	this->scLibSimp = scLibSimp;
	this->preAAType = 0;

	this->targetNode = new ResNode(lines[0], 0, rn, bbLib, scLib);
	vector<string> spt;

	splitString(lines[0], " ", &spt);
	int bbIndex1w = atoi(spt[4].c_str());
	int bbIndex1k = bbLib->index1wTo1k[bbIndex1w];

	for(int i=1;i<lines.size();i++){
		splitString(lines[i], " ", &spt);
		int sep = atoi(spt[1].c_str());
		this->sepList.push_back(sep);

		this->neighborNodes.push_back(new ResNode(lines[i],i,  rn, bbLib, scLib));
	}

	for(int i=0;i<neighborNodes.size();i++){
		if(sepList[i] == -1){
			preAAType = neighborNodes[i]->aaType;
		}
	}

	this->riTarget = new ResInfo(targetNode->aaType, targetNode->ss, targetNode->sai, bbIndex1k);

	XYZ targetCB = targetNode->getCbCoord();
	ResNode* nd;
	for(int i=0;i<neighborNodes.size();i++){
		nd = neighborNodes[i];
		XYZ nbCB = nd->getCbCoord();
		if(targetCB.distance(nbCB) < 8){
			nbAAList.push_back(nd->aaType);
			int index = nd->seqID;
			int sep = abs(index);
			if(sep > 5) sep = 5;
			if(index < 0) {
				CsMove cm = targetNode->cs2 - nd->cs2;
				rpList.push_back(new ResPairInfo(index, 0, nd->aaType, targetNode->aaType, nd->ss, targetNode->ss, nd->sai, targetNode->sai, cm, sep));
				targetFirst.push_back(false);
			}
			else {
				CsMove cm = nd->cs2 - targetNode->cs2;
				rpList.push_back(new ResPairInfo(0, index, targetNode->aaType, nd->aaType, targetNode->ss, nd->ss, targetNode->sai, nd->sai, cm, sep));
				targetFirst.push_back(true);
			}
		}
	}
}



void BuildMutationTemplate::loadS1S2(ProS1S2Energy* eS1S2){
	this->sa = eS1S2->getS1(riTarget);

	this->smList.clear();
	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}
}

void BuildMutationTemplate::printS1S2(ofstream& out){
	this->sa.print(out);
	int n = nbAAList.size();
	out << n << endl;
	for(int i=0;i<n;i++){
		this->smList[i].print(out);
	}
}

void BuildMutationTemplate::loadS1S2(vector<string>& lines){
	this->sa = AAScoreArray(lines[0]);
	this->smList.clear();
	int n = atoi(lines[1].c_str());
	double* matrix = new double[400];
	vector<string> spt;
	int begin;
	for(int i=0;i<n;i++){
		begin = i*21 + 2;
		for(int j=0;j<20;j++){
			splitString(lines[j+begin], " ", &spt);
			for(int k=0;k<20;k++){
				matrix[j*20+k] = atof(spt[k].c_str());
			}
		}
		AAScoreMatrix sm(matrix, 1.0);
		this->smList.push_back(sm);
	}
	delete [] matrix;
}

int BuildMutationTemplate::rankS1(ProS1S2Energy* eS1S2){
	vector<double> aaEnergyList;
	AAScoreArray sa = eS1S2->getS1(riTarget);
	ResName rn;
	for(int aa=0;aa<20;aa++){
		double s1 = sa.sa[aa];
		aaEnergyList.push_back(s1);
	}

//	cout << "s1 nat type: " << rn.intToTri(this->targetNode->aaType) << endl;
//	cout << riTarget->ss << " " << riTarget->sai << " " << riTarget->bbIndex << " " << riTarget->aaType << endl;
	int rank = 1;
	for(int i=0;i<20;i++){
//		printf("%3s %7.3f\n", rn.intToTri(i).c_str(), aaEnergyList[i]);
		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}
	return rank;
}

int BuildMutationTemplate::rankS2(ProS1S2Energy* eS1S2) {
	vector<double> aaEnergyList;

	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;
	string nat = rn.intToTri(this->targetNode->aaType);

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();
		double s2 = 0;

		double pairEne = 0;
		for(int i=0;i<n;i++){
			ResPairInfo* rp = this->rpList[i];
			//cout << rp->toString() << endl;
			if(targetFirst[i]){
				pairEne = smList[i].sm[aa][nbAAList[i]];
				s2 += smList[i].sm[aa][nbAAList[i]];
			}
			else{
				pairEne = smList[i].sm[nbAAList[i]][aa];
				s2 += smList[i].sm[nbAAList[i]][aa];
			}
			//printf("nat %s pair: %2d %s %s %6.3f\n", nat.c_str(), this->rpList[i]->sep, rn.intToTri(aa).c_str(), rn.intToTri(nbAAList[i]).c_str(), pairEne);
		}

		aaEnergyList.push_back(s2);
	}

	//cout << "s2 nat type: " << rn.intToTri(this->targetNode->aaType) << endl;
	int rank = 1;
	for(int i=0;i<20;i++){
		//printf("%3s %7.3f\n", rn.intToTri(i).c_str(), aaEnergyList[i]);
		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}


	return rank;
}



int BuildMutationTemplate::rankS1S2(ProS1S2Energy* eS1S2){
	vector<double> aaEnergyList;

	AAScoreArray sa = eS1S2->getS1(riTarget);
	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;

	vector<double> s1List;
	vector<double> s2List;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double s1 = sa.sa[aa];
		double s2 = 0;

		for(int i=0;i<n;i++){
			if(targetFirst[i]){
				s2 += smList[i].sm[aa][nbAAList[i]];
			}
			else{
				s2 += smList[i].sm[nbAAList[i]][aa];
			}
		}

		s1List.push_back(s1);
		s2List.push_back(s2);
		aaEnergyList.push_back(s1+s2);
	}

//	cout << "s1s2 nat type: " << rn.intToTri(this->targetNode->aaType) << endl;
	int rank = 1;
	for(int i=0;i<20;i++){
//		printf("%3s %7.3f %7.3f %7.3f\n", rn.intToTri(i).c_str(), s1List[i], s2List[i], aaEnergyList[i]);
		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}
//	cout << "rank: " << rank << endl;
	return rank;
}

int BuildMutationTemplate::rankAtomic(EnergyCalculator* ec){
	vector<double> aaEnergyList;
	ResName rn;

	//cout << "nat " << rn.intToTri(this->targetNode->aaType) << endl;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double minEne = DBL_MAX;
		int minID = 0;
		int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aa);
		if(clusterNum == 0) {
			for(int i=0;i<rotNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[i]);
				//double e = scLib->getEnergy(targetNode->conBB->rot->index1K, targetNode->conScTmp->rot) * ec->pa->wtRot;
				double e = 0;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					//double e1 = ec->getEnergyScScDesign(targetNode->conScTmp, neighborNodes[j]->conSc, meanSai, neighborNodes[j]->seqID) ;
					//double e2 = ec->getEnergyBbScDesign(neighborNodes[j]->conBB, targetNode->conScTmp, meanSai, -neighborNodes[j]->seqID);
					//printf("sai: %5.3f rot: %3d neighbor: %2d e1: %7.3f e2: %7.3f\n", meanSai, i, j, e1, e2);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					minID = i;
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		else {
			//search level1
			int level1MinID = 0;
			for(int i=0;i<clusterNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][i][1]]);
				//double e = scLib->getEnergy(targetNode->conBB->rot->index1K, targetNode->conScTmp->rot) * ec->pa->wtRot;
				double e = 0;

				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					//double e1 = ec->getEnergyScScDesign(targetNode->conScTmp, neighborNodes[j]->conSc, meanSai, neighborNodes[j]->seqID);
					//double e2 = ec->getEnergyBbScDesign(neighborNodes[j]->conBB, targetNode->conScTmp, meanSai, -neighborNodes[j]->seqID);
					//printf("sai: %5.3f rot: %3d neighbor: %2d e1: %7.3f e2: %7.3f\n", meanSai, i, j, e1, e2);

					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					level1MinID = i;
					minID = scLib->rotCluster[aa][i][1];
				}
			}
			//search level2
			int clusterMemberNum = scLib->rotCluster[aa][level1MinID][0];
			for(int i=0;i<clusterMemberNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][level1MinID][i+1]]);
				//double e = scLib->getEnergy(targetNode->conBB->rot->index1K, targetNode->conScTmp->rot) * ec->pa->wtRot;
				double e = 0;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = scLib->rotCluster[aa][level1MinID][i+1];
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		aaEnergyList.push_back(minEne);
		//printf("%s %7.3f\n", rn.intToTri(aa).c_str(), minEne);
	}

	int rank = 1;
	string natAA = rn.intToTri(this->targetNode->aaType);
	string minEneAA = "";
	double minEne = DBL_MAX;

	for(int i=0;i<20;i++){
		if(aaEnergyList[i] < minEne){
			minEne = aaEnergyList[i];
			minEneAA = rn.intToTri(i);
		}

		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}

	//char x[200];
	//sprintf(x, "%-2d %s %s %5.3f", rank, natAA.c_str(), minEneAA.c_str(), targetNode->sai);

	return rank;
}

int BuildMutationTemplate::rankAtomicABACUS(EnergyCalculator* ec){
	vector<double> aaEnergyList;
	ResName rn;

	//cout << "nat " << rn.intToTri(this->targetNode->aaType) << endl;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double minEne = DBL_MAX;
		int minID = 0;
		int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aa);
		if(clusterNum == 0) {
			for(int i=0;i<rotNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[i]);
				//double e = scLib->getEnergy(targetNode->conBB->rot->index1K, targetNode->conScTmp->rot) * ec->pa->wtRot;
				double e = 0;
				double de = 0;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					//double e1 = ec->getEnergyScScDesign(targetNode->conScTmp, neighborNodes[j]->conSc, meanSai, neighborNodes[j]->seqID) ;
					//double e2 = ec->getEnergyBbScDesign(neighborNodes[j]->conBB, targetNode->conScTmp, meanSai, -neighborNodes[j]->seqID);
					//printf("sai: %5.3f rot: %3d neighbor: %2d e1: %7.3f e2: %7.3f\n", meanSai, i, j, e1, e2);
					e += ec->getEnergyABACUS(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);

				}

				if(e < minEne){
					minEne = e;
					minID = i;
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		else {
			//search level1
			int level1MinID = 0;
			for(int i=0;i<clusterNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][i][1]]);
				//double e = scLib->getEnergy(targetNode->conBB->rot->index1K, targetNode->conScTmp->rot) * ec->pa->wtRot;
				double e = 0;

				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					//double e1 = ec->getEnergyScScDesign(targetNode->conScTmp, neighborNodes[j]->conSc, meanSai, neighborNodes[j]->seqID);
					//double e2 = ec->getEnergyBbScDesign(neighborNodes[j]->conBB, targetNode->conScTmp, meanSai, -neighborNodes[j]->seqID);
					//printf("sai: %5.3f rot: %3d neighbor: %2d e1: %7.3f e2: %7.3f\n", meanSai, i, j, e1, e2);

					e += ec->getEnergyABACUS(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					level1MinID = i;
					minID = scLib->rotCluster[aa][i][1];
				}
			}
			//search level2
			int clusterMemberNum = scLib->rotCluster[aa][level1MinID][0];
			for(int i=0;i<clusterMemberNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][level1MinID][i+1]]);
				//double e = scLib->getEnergy(targetNode->conBB->rot->index1K, targetNode->conScTmp->rot) * ec->pa->wtRot;
				double e = 0;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyABACUS(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = scLib->rotCluster[aa][level1MinID][i+1];
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		aaEnergyList.push_back(minEne);
		//printf("%s %7.3f\n", rn.intToTri(aa).c_str(), minEne);
	}

	int rank = 1;
	string natAA = rn.intToTri(this->targetNode->aaType);
	string minEneAA = "";
	double minEne = DBL_MAX;

	for(int i=0;i<20;i++){
		if(aaEnergyList[i] < minEne){
			minEne = aaEnergyList[i];
			minEneAA = rn.intToTri(i);
		}

		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}

	//char x[200];
	//sprintf(x, "%-2d %s %s %5.3f", rank, natAA.c_str(), minEneAA.c_str(), targetNode->sai);

	return rank;
}

int BuildMutationTemplate::nativeRankDetail(EnergyCalculator* ec, ProS1S2Energy* eS1S2){
	vector<double> aaEnergyList;

	AAScoreArray sa = eS1S2->getS1(riTarget);
	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;

	//printf("neighbor %2d\n", neighborNodes.size());

	//cout << "nat " << rn.intToTri(this->targetNode->aaType) << endl;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double s1 = sa.sa[aa];
		double s2 = 0;

		for(int i=0;i<n;i++){
			if(targetFirst[i]){
				s2 += smList[i].sm[aa][nbAAList[i]];
			}
			else{
				s2 += smList[i].sm[nbAAList[i]][aa];
			}
		}

		double minEne = DBL_MAX;
		int minID = 0;
		int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aa);
		if(clusterNum == 0) {
			for(int i=0;i<rotNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[i]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = i;
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		else {
			//search level1
			int level1MinID = 0;
			for(int i=0;i<clusterNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][i][1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					level1MinID = i;
					minID = scLib->rotCluster[aa][i][1];
				}
			}
			//search level2
			int clusterMemberNum = scLib->rotCluster[aa][level1MinID][0];
			for(int i=0;i<clusterMemberNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][level1MinID][i+1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = scLib->rotCluster[aa][level1MinID][i+1];
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		aaEnergyList.push_back(s1 + s2 + minEne);

	}

	int rank = 1;
	for(int i=0;i<20;i++){
		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}
	return rank;
}


void BuildMutationTemplate::printBBInfo(){
	XYZ t1 = this->targetNode->conf->bbConf->cs1.origin_;
	XYZ t2 = this->targetNode->conf->bbConf->cs2.origin_;
	XYZ t3 = this->targetNode->conf->bbConf->cs3.origin_;
	XYZ t4 = this->targetNode->conf->bbConf->cs4.origin_;

	cout << "N:  " << t1.toString() << endl;
	cout << "CA: " << t2.toString() << endl;
	cout << "C:  " << t3.toString() << endl;
	cout << "O:  " << t4.toString() << endl;
}

void BuildMutationTemplate::printNodeSep() {

	XYZ targetCB = targetNode->getCbCoord();
	ResNode* nd;
	for(int i=0;i<neighborNodes.size();i++){
		nd = neighborNodes[i];
		XYZ nbCB = nd->getCbCoord();
		if(targetCB.distance(nbCB) < 8){
			nbAAList.push_back(nd->aaType);
			int index = nd->seqID;
			int sep = abs(index);
			if(sep > 5) sep = 5;
			if(index < 0) {
				CsMove cm = targetNode->cs2 - nd->cs2;
				int s1 = sep;
				int s2 = this->sepList[i];
				printf("neighbor %d %d %d\n", i, s1, s2);
			}
			else {
				CsMove cm = nd->cs2 - targetNode->cs2;
				int s1 = sep;
				int s2 = this->sepList[i];
				printf("neighbor %d %d %d\n", i, s1, s2);
			}
		}
	}
}

void BuildMutationTemplate::printNativeEnergy(EnergyCalculator* ec, ProS1S2Energy* eS1S2) {
	AAScoreArray sa = eS1S2->getS1(riTarget);
	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;

	//cout << "nat " << rn.intToTri(this->targetNode->aaType) << endl;

	{
		int aa = riTarget->aaType;
		ResScRotamer* rot = this->targetNode->conf->scConf->rot;
		int neighborNum = neighborNodes.size();

		double s1 = sa.sa[aa];
		printf("S1: %6.3f\n", s1);
		double s2 = 0;

		for(int i=0;i<n;i++){
			if(targetFirst[i]){
				s2 = smList[i].sm[aa][nbAAList[i]];
			}
			else{
				s2 = smList[i].sm[nbAAList[i]][aa];
			}
			printf("S2: %2d %2d %6.3f\n", rpList[i]->indexA, rpList[i]->indexB, s2);
		}

		double minEne = DBL_MAX;
		int minID = 0;
		int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aa);

		targetNode->conf->updateScRotamer(rot);
		double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
		printf("ROT: %6.3f\n", e);
		for(int j=0;j<neighborNum;j++){
			float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
			e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
			printf("sai: %5.3f neighbor: %2d ene: %7.3f\n", meanSai,  j, e);
		}
	}
}


string BuildMutationTemplate::nativeRank(EnergyCalculator* ec, ProS1S2Energy* eS1S2){
	vector<double> aaEnergyList;
	AAScoreArray sa = eS1S2->getS1(riTarget);
	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;

	//cout << "nat " << rn.intToTri(this->targetNode->aaType) << endl;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLibSimp->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double s1 = sa.sa[aa];

		double s2 = 0;

		for(int i=0;i<n;i++){
			if(targetFirst[i]){
				s2 += smList[i].sm[aa][nbAAList[i]];
			}
			else{
				s2 += smList[i].sm[nbAAList[i]][aa];
			}
		}

		double minEne = DBL_MAX;
		int minID = 0;

		for(int i=0;i<rotNum;i++){
			targetNode->confTmp->updateScRotamer(rotList[i]);
			double e = ec->dp->ref[20*targetNode->ssInt + aa];
			//printf("ref: %s %6.3f\n", rn.intToTri(aa).c_str(), e);
			e += scLibSimp->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
			for(int j=0;j<neighborNum;j++){
				float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
				e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
			}
			if(e < minEne){
				minEne = e;
				minID = i;
			}
		}
		targetNode->confTmp->updateScRotamer(rotList[minID]);

		aaEnergyList.push_back(s1 + s2 + minEne);
	}

	int rank = 1;
	string natAA = rn.intToTri(this->targetNode->aaType);
	string minEneAA = "";
	double minEne = DBL_MAX;

	for(int i=0;i<20;i++){
		if(aaEnergyList[i] < minEne){
			minEne = aaEnergyList[i];
			minEneAA = rn.intToTri(i);
		}

		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}

	char xx[20];
	sprintf(xx, "%s %s %2d", natAA.c_str(), minEneAA.c_str(), rank);
	return string(xx);
}


int BuildMutationTemplate::nativeRankABACUS(EnergyCalculator* ec, ProS1S2Energy* eS1S2){
	vector<double> aaEnergyList;
	AAScoreArray sa = eS1S2->getS1(riTarget);
	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;

	//cout << "nat " << rn.intToTri(this->targetNode->aaType) << endl;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double s1 = sa.sa[aa];
		double s2 = 0;

		for(int i=0;i<n;i++){
			if(targetFirst[i]){
				s2 += smList[i].sm[aa][nbAAList[i]];
			}
			else{
				s2 += smList[i].sm[nbAAList[i]][aa];
			}
		}

		double minEne = DBL_MAX;
		int minID = 0;
		int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aa);
		if(clusterNum == 0) {
			for(int i=0;i<rotNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[i]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyABACUS(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = i;
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		else {
			//search level1
			int level1MinID = 0;
			for(int i=0;i<clusterNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][i][1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;

				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyABACUS(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					level1MinID = i;
					minID = scLib->rotCluster[aa][i][1];
				}
			}
			//search level2
			int clusterMemberNum = scLib->rotCluster[aa][level1MinID][0];
			for(int i=0;i<clusterMemberNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][level1MinID][i+1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyABACUS(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = scLib->rotCluster[aa][level1MinID][i+1];
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		aaEnergyList.push_back(s1 + s2 + minEne);
	}

	int rank = 1;
	string natAA = rn.intToTri(this->targetNode->aaType);
	string minEneAA = "";
	double minEne = DBL_MAX;

	for(int i=0;i<20;i++){
		if(aaEnergyList[i] < minEne){
			minEne = aaEnergyList[i];
			minEneAA = rn.intToTri(i);
		}

		if(aaEnergyList[i] < aaEnergyList[this->targetNode->aaType]){
			rank ++;
		}
	}
	return rank;
}

double BuildMutationTemplate::nativeLogits(EnergyCalculator* ec, ProS1S2Energy* eS1S2){
	vector<double> aaEnergyList;
	AAScoreArray sa = eS1S2->getS1(riTarget);
	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;

	//cout << "nat " << rn.intToTri(this->targetNode->aaType) << endl;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double s1 = sa.sa[aa];

		double s2 = 0;

		for(int i=0;i<n;i++){
			if(targetFirst[i]){
				s2 += smList[i].sm[aa][nbAAList[i]];
			}
			else{
				s2 += smList[i].sm[nbAAList[i]][aa];
			}
		}

		double minEne = DBL_MAX;
		int minID = 0;
		int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aa);
		if(clusterNum == 0) {
			for(int i=0;i<rotNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[i]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = i;
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		else {
			//search level1
			int level1MinID = 0;
			for(int i=0;i<clusterNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][i][1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;

				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					level1MinID = i;
					minID = scLib->rotCluster[aa][i][1];
				}
			}
			//search level2
			int clusterMemberNum = scLib->rotCluster[aa][level1MinID][0];
			for(int i=0;i<clusterMemberNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][level1MinID][i+1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}
				if(e < minEne){
					minEne = e;
					minID = scLib->rotCluster[aa][level1MinID][i+1];
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		aaEnergyList.push_back(s1 + s2 + minEne);
	}

	double par = 0;
	for(int i=0;i<20;i++){
		par += exp(-1.0*aaEnergyList[i]);
	}
	double p = exp(-1.0*aaEnergyList[targetNode->aaType])/par;
	return -log(p);
}

int BuildMutationTemplate::singleResidueDesign(EnergyCalculator* ec, ProS1S2Energy* eS1S2){

	vector<double> aaEnergyList;
	vector<ResScRotamer*> aaBestRotList;

	AAScoreArray sa = eS1S2->getS1(riTarget);
	vector<AAScoreMatrix> smList;

	int n = nbAAList.size();
	for(int i=0;i<n;i++){
		smList.push_back(eS1S2->getS2(rpList[i]));
	}

	ResName rn;

	for(int aa=0;aa<20;aa++){
		vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
		int rotNum = rotList.size();
		int neighborNum = neighborNodes.size();

		double s1 = sa.sa[aa];

		double s2 = 0;

		for(int i=0;i<n;i++){
			if(targetFirst[i]){
				s2 += smList[i].sm[aa][nbAAList[i]];
			}
			else{
				s2 += smList[i].sm[nbAAList[i]][aa];
			}
		}

		double minEne = DBL_MAX;
		int minID = 0;
		int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aa);
		ResScRotamer* bestRot;

		if(clusterNum == 0) {
			for(int i=0;i<rotNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[i]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					minID = i;
					bestRot = rotList[i];
				}

			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		else {
			//search level1
			int level1MinID = 0;
			for(int i=0;i<clusterNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][i][1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					level1MinID = i;
					minID = scLib->rotCluster[aa][i][1];
					bestRot = rotList[scLib->rotCluster[aa][i][1]];
				}
			}
			//search level2
			int clusterMemberNum = scLib->rotCluster[aa][level1MinID][0];
			for(int i=0;i<clusterMemberNum;i++){
				targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[aa][level1MinID][i+1]]);
				double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->dp->wtRot;
				for(int j=0;j<neighborNum;j++){
					float meanSai = 0.5*(targetNode->sai+neighborNodes[j]->sai);
					e += ec->getEnergyDesign(targetNode->confTmp, neighborNodes[j]->conf, meanSai, neighborNodes[j]->seqID);
				}

				if(e < minEne){
					minEne = e;
					minID = scLib->rotCluster[aa][level1MinID][i+1];
					bestRot = rotList[scLib->rotCluster[aa][level1MinID][i+1]];
				}
			}
			targetNode->confTmp->updateScRotamer(rotList[minID]);
		}
		aaEnergyList.push_back(s1 + s2 + minEne);
		aaBestRotList.push_back(bestRot);
	}

	double totMin = 999999.9;
	int minEAA = 0;
	for(int i=0;i<20;i++){
		if(aaEnergyList[i] < totMin){
			totMin = aaEnergyList[i];
			minEAA = i;
		}
	}
	return minEAA;
}

BuildMutationTemplate::~BuildMutationTemplate() {

	delete targetNode;
	delete riTarget;
	for(int i=0;i<neighborNodes.size();i++){
		delete neighborNodes[i];
	}
	for(int i=0;i<rpList.size();i++){
		delete rpList[i];
	}
}

} /* namespace NSPmodel */
