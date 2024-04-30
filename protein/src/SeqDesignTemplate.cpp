/*
 * SeqDesignTemplate.cpp
 *
 */

#include "protein/SeqDesignTemplate.h"

namespace NSPprotein {

SeqDesignTemplate::SeqDesignTemplate(const string& inputFile, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, AtomLib* atLib, EnergyCalculator* ec, ProS1S2Energy* eS1S2) {

	cout << "inputfile: " << inputFile << endl;
	InputParser* input = new InputParser(inputFile);
	string pdbFile = input->getValue("PDB");
	string resFile = input->getValue("RESFILE");
	cout << pdbFile << endl;
	cout << resFile << endl;

	this->bbLib = bbLib;
	this->scLib = scLib;
	this->atLib = atLib;
	this->ec = ec;
	this->eS1S2 = eS1S2;

	this->T0 = 30.0;
	this->T1 = 0.1;
	this->stepNumFactor = 1000;

	if(input->specifiedOption("T0"))
		this->T0 = atof(input->getValue("T0").c_str());
	if(input->specifiedOption("T1"))
		this->T1 = atof(input->getValue("T1").c_str());
	if(input->specifiedOption("stepPerRes"))
		this->stepNumFactor = atof(input->getValue("stepPerRes").c_str());

	this->rn = new ResName();


	this->pdb = new PDB(pdbFile, "xxxx");

	vector<Residue*> allRes = pdb->getResList();
	int n = allRes.size();
	for(int i=0;i<n;i++){
		Residue* res = allRes[i];
		if(res->hasThreeCoreAtoms() && res->getAtom("O") != NULL)
		{
			resList.push_back(res);
			if(res->sidechainComplete(atLib))
				resValidList.push_back(true);
			else
				resValidList.push_back(false);
		}
		else {
			cout << "RES " << res->chainID << " " << res->resID << " backbone incomplete, ignored" << endl;
		}
	}

	char xxx[resList.size()+1];
	for(int i=0;i<resList.size();i++){
		xxx[i] = rn->intToSin(resList[i]->intName);
	}
	xxx[resList.size()] = '\0';
	this->initSequence = string(xxx);

	int id = 0;
	seqIDList.push_back(0);
	for(int i=1;i<resList.size();i++) {
		if(resList[i-1]->contactToNeighbor(resList[i])){
			id++;
			seqIDList.push_back(id);
		}
		else{
			id+=5;
			seqIDList.push_back(id);
		}
	}

	nodeList.push_back(new ResNode(resList[0], seqIDList[0], 0,  this->atLib, this->bbLib, scLib));
	for(int i=1;i<resList.size();i++) {
		if(seqIDList[i] == seqIDList[i-1]+1){
			nodeList.push_back(new ResNode(resList[i-1], resList[i], seqIDList[i],i,  this->atLib, this->bbLib, scLib));
		}
		else {
			nodeList.push_back(new ResNode(resList[i], seqIDList[i], i, this->atLib, this->bbLib, scLib));
		}
	}


	for(int i=0;i<nodeList.size();i++){
		initScIDList.push_back(nodeList[i]->conf->scConf->rot->rotID);
		initRotList.push_back(nodeList[i]->conf->scConf->rot);
	}

	vector<XYZ> cbList;
	for(int i=0;i<resList.size();i++){
		cbList.push_back(resList[i]->getCbCoord());
	}

	for(int i=0;i<resList.size();i++){
		vector<int> cList;
		for(int j=0;j<resList.size();j++){
			if(i==j) continue;
			if(cbList[i].distance(cbList[j]) < 16)
				cList.push_back(j);
		}
		nbList.push_back(cList);
	}

	this->resNum = resList.size();
	char xx[20];
	for(int i=0;i<resNum;i++){
		sprintf(xx,"%c%s",resList[i]->chainID,resList[i]->resID.c_str());
		string key = string(xx);
		chainIDResIDToIndex[key] = i;
	}

	ResSasaPoints* rsp = new ResSasaPoints();
	AssignProSSAndSasa* ssa = new AssignProSSAndSasa(pdb);
	if(ssa->resNum != resList.size()) {
		cout << "pdb residue num is not equal to secondary structure residue num: " << resList.size() << " != " << ssa->resNum << endl;
	}

	ssa->updateSS();
	ssa->updateSasa(rsp);
	Residue* res;
	ResBBRotamer* bbRot;
	for(int i=0;i<resList.size();i++){
		res = resList[i];
		bbRot = nodeList[i]->conf->bbConf->rot;
		ResInfo* ri = new ResInfo(res->intName, ssa->getSS(i), ssa->getSASA(i), this->bbLib->getRotamerIndex1K(bbRot));
		this->nodeList[i]->ss = ri->ss;
		this->nodeList[i]->sai = ri->sai;
		this->riList.push_back(ri);
	}

	for(int i=0;i<resList.size();i++){
		involvedRpIndex.push_back(vector<int>());
	}

	int rpIndex = 0;
	for(int i=0;i<resList.size();i++){
		LocalFrame csA = resList[i]->coordSys;
		XYZ cbA = resList[i]->getCbCoord();
		for(int j=i+1;j<resList.size();j++){
			LocalFrame csB = resList[j]->coordSys;
			XYZ cbB = resList[j]->getCbCoord();
			CsMove cm = csB - csA;
			if(cbA.distance(cbB) > 8.0) continue;
			ResPairInfo* rp = new ResPairInfo(i, j, resList[i]->intName, resList[j]->intName, riList[i]->ss, riList[j]->ss, riList[i]->sai, riList[j]->sai, cm, seqIDList[j] - seqIDList[i]);
			involvedRpIndex[i].push_back(rpIndex);
			involvedRpIndex[j].push_back(rpIndex);
			rpIndex++;
			this->rpList.push_back(rp);
		}
	}

	//cout << "load S1S2" << endl;
	loadS1S2();

	double profWeight = 1.0;
	if(input->specifiedOption("PRIFILE_WEIGHT"))
		profWeight = atof(input->getValue("PRIFILE_WEIGHT").c_str());
	if(input->specifiedOption("PROFILE"))
		loadMSAProfile(input->getValue("PROFILE"), profWeight);


	double rotChoiceEnergyCutoff = 35.0;
	if(input->specifiedOption("ROT_ENERGY_CUTOFF"))
		rotChoiceEnergyCutoff = atof(input->getValue("ROT_ENERGY_CUTOFF").c_str());
	updateRotChoice(resFile, rotChoiceEnergyCutoff);


	delete rsp;
	delete ssa;
	delete input;
}

void SeqDesignTemplate::printInvolvedPairInfo(){
	for(int i=0;i<involvedRpIndex.size();i++){
		cout << i << " ";
		for(int j=0;j<involvedRpIndex[i].size();j++){
			cout << involvedRpIndex[i][j] << " ";
		}
		cout << endl;
	}
}

void SeqDesignTemplate::loadS1S2(){
	saList.clear();
	smList.clear();
	for(int i=0;i<riList.size();i++){
		this->saList.push_back(eS1S2->getS1(riList[i]));
	}
	cout << "pairNum: " << rpList.size() << endl;
	for(int i=0;i<rpList.size();i++){
		cout << i << endl;
		this->smList.push_back(eS1S2->getS2(rpList[i]));
	}
}

void SeqDesignTemplate::printS2(const string& outfile){

	ofstream out;
	out.open(outfile.c_str(), ios::out);

	char xx[200];
	for(int pid=0;pid<rpList.size();pid++){
		int indexA = rpList[pid]->indexA;
		int indexB = rpList[pid]->indexB;
		int aaA = nodeList[indexA]->aaType;
		int aaB = nodeList[indexB]->aaType;
		double s2 = smList[pid].sm[aaA][aaB];
		sprintf(xx, "pair %3d %s %3d %s %6.3f\n", indexA, rn->intToTri(aaA).c_str(), indexB, rn->intToTri(aaB).c_str(), s2);
		out << string(xx);
		for(int i=0;i<20;i++){
			for(int j=0;j<20;j++){
				sprintf(xx, "%6.3f ", smList[pid].sm[i][j]);
				out << string(xx);
			}
			out << endl;
		}
	}
	out.close();
}

void SeqDesignTemplate::loadS1S2(const string& s2File){
	saList.clear();
	smList.clear();

	//cout << "get S1" << endl;
	for(int i=0;i<riList.size();i++){
		this->saList.push_back(eS1S2->getS1(riList[i]));
	}

	ifstream f;
	f.open(s2File.c_str(), ios::in);
	string s;
	vector<string> spt;
	double sm[400];

	//cout << "get S2" << endl;
	for(int i=0;i<rpList.size();i++){

		getline(f, s);
		//cout << s << endl;
		for(int j=0;j<20;j++){
			getline(f, s);
			splitString(s, " ", &spt);
			for(int k=0;k<20;k++){
				sm[j*20+k] = atof(spt[k].c_str());
			}
		}

		AAScoreMatrix matrix(sm, -1);
		this->smList.push_back(matrix);
	}
	//cout << rpList.size() << endl;
	//cout << smList.size() << endl;
	f.close();
}

void SeqDesignTemplate::printBackboneInfo(){

	for(int i=0;i<5;i++) {
		cout << "res " << i	 <<endl;
		Residue* res = resList[i];
		XYZ n = res->getAtom("N")->getCoord();
		XYZ ca = res->getAtom("CA")->getCoord();
		XYZ c = res->getAtom("C")->getCoord();
		XYZ o = res->getAtom("O")->getCoord();

		cout << "N  " << n.toString() << endl;
		cout << "CA " << ca.toString() << endl;
		cout << "C  " << c.toString() << endl;
		cout << "O  " << o.toString() << endl;

		LocalFrame cs1 = nodeList[i]->conf->bbConf->cs1;
		LocalFrame cs2 = nodeList[i]->conf->bbConf->cs2;
		LocalFrame cs3 = nodeList[i]->conf->bbConf->cs3;
		LocalFrame cs4 = nodeList[i]->conf->bbConf->cs4;

		cout << "cs1 " << cs1.origin_.toString() << endl;
		cout << "cs2 " << cs2.origin_.toString() << endl;
		cout << "cs3 " << cs3.origin_.toString() << endl;
		cout << "cs4 " << cs4.origin_.toString() << endl;
	}
}

void SeqDesignTemplate::printResInfo(){
	for(int i=0;i<resNum;i++){
		char ss = riList[i]->ss;
		double sai = riList[i]->sai;
		printf("%c %5.3f\n", ss, sai);
	}
}

void SeqDesignTemplate::printPairInfo(){
	for(int pid=0;pid<rpList.size();pid++){
		int indexA = rpList[pid]->indexA;
		int indexB = rpList[pid]->indexB;
		ResPairInfo* rp = this->rpList[pid];
		cout << rp->toString() << endl;
	}
}

void SeqDesignTemplate::loadMSAProfile(const string& file, double weight){
	ifstream in;
	in.open(file.c_str(), ios::in);
	string s;
	vector<string> spt;
	int index = 0;

	cout << "load MSA prof: " << file << endl;
	cout << "prof weight: " << weight << endl;
	float profEne[20];

	float eProf;
	while(getline(in, s)){
		splitString(s, " ", &spt);
		for(int i=0;i<20;i++){
			 eProf = atof(spt[i].c_str());
			 this->saList[index].sa[i] += eProf;
		}
		index++;
	}
	if(index != this->resList.size()){
		cout << "MSA profile length not equal to residue number" << endl;
		exit(1);
	}
	in.close();
}

void SeqDesignTemplate::updateRotChoice(const string& resFile,  double eneCutoff){
	ifstream rf;
	rf.open(resFile.c_str(),ios::in);
	if(!rf.is_open()){
		cout << "fail to open file " << resFile << endl;
		exit(1);
	}
	string s;

	cout << "selected amino acid choice from resFile" << endl;

	vector<string> posAAChoice;
	for(int i=0;i<resNum;i++){
		posAAChoice.push_back("ADEFGHIKLMNPQRSTVWY");
	}

	string allButCys = "ADEFGHIKLMNPQRSTVWY";
	string all = "ACDEFGHIKLMNPQRSTVWY";
	string defaultValue = "nat";
	set<int> fixed; //rotamer fixed
	vector<string> spt;


	while(getline(rf,s)){
		splitString(s," ",&spt);
		if(spt.size() == 2 && spt[0] == "default") {
			if(spt[1] == "nat") {
				for(int i=0;i<resNum;i++){
					int aaType = this->resList[i]->intName;
					char xx[2];
					xx[0] = rn->intToSin(aaType);
					xx[1] = '\0';
					posAAChoice[i] = string(xx);
				}
			}
			else if(spt[1] == "all") {
				for(int i=0;i<resNum;i++){
					posAAChoice[i] = all;
				}
			}
			else if(spt[1] == "allButCys") {
				for(int i=0;i<resNum;i++){
					posAAChoice[i] = allButCys;
				}
			}
			else {
				for(int i=0;i<resNum;i++){
					posAAChoice[i] = spt[1];
				}
			}
			continue;
		}
		if(spt.size() != 3){
			cout << "invalid res choice: " << s << endl;
			exit(0);
		}

		string key = spt[0] + spt[1];
		if(chainIDResIDToIndex.find(key) == chainIDResIDToIndex.end()){
			cout << "invalid resID: " << s << endl;
			exit(0);
		}

		int index = chainIDResIDToIndex[key];
		if(spt[2] == "all")
			posAAChoice[index] = all;
		else if(spt[2] == "allButCys")
			posAAChoice[index] = allButCys;
		else if(spt[2] == "nat"){
			int aaType = this->resList[index]->intName;
			char xx[2];
			xx[0] = rn->intToSin(aaType);
			xx[1] = '\0';
			posAAChoice[index] = string(xx);
		}
		else if(spt[2] == "rot") {
			int aaType = this->resList[index]->intName;
			char xx[2];
			xx[0] = rn->intToSin(aaType);
			xx[1] = '\0';
			posAAChoice[index] = string(xx);
			fixed.insert(index);
		}
		else
			posAAChoice[index] = spt[2];
	}

	rf.close();

	int aaType;
	double e;
	for(int pos=0;pos<resNum;pos++){
		string aaChoice = posAAChoice[pos];
		if(aaChoice.size() == 0) {
			cout << "invalid aa choice: " << aaChoice << endl;
			exit(0);
		}
		vector<ResScRotamer*> rotList;
		double cutoff = eneCutoff;
		vector<ResScRotamer*> selectRotList;

		if(fixed.find(pos) != fixed.end()){
			selectRotList.push_back(this->initRotList[pos]);
			this->resMutatorList.push_back(new ResMutator(selectRotList));
			cout << "pos: " << pos << " choice: " << aaChoice  << " rotNum: " << selectRotList.size() << endl;
			continue;
		}

		while(selectRotList.size() == 0){
			selectRotList.clear();
			for(int i=0;i<aaChoice.length();i++){
				aaType = rn->sinToInt(aaChoice[i]);
				if(aaType < 0 || aaType > 19){
					cout << "invalid aa choice: " << aaChoice << endl;
					exit(0);
				}
				vector<ResScRotamer*> rotList = this->scLib->rotList[aaType];
				int rotNum = rotList.size();
				for(int i=0;i<rotNum;i++){
					e = rotEnergyWithBackbone(pos, rotList[i]);
					if(e < eneCutoff){
						selectRotList.push_back(rotList[i]);
					}
				}

				/*
				int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(aaType);
				if(clusterNum == 0){
					for(int i=0;i<rotNum;i++){
						e = rotEnergyWithBackbone(pos, rotList[i]);
						if(e < eneCutoff){
							selectRotList.push_back(rotList[i]);
						}
					}
				}
				else {
					for(int i=0;i<clusterNum;i++){
						e = rotEnergyWithBackbone(pos, rotList[scLib->rotClusterUnique[aaType][i][1]]);
						if(e < cutoff) {
							int clusterMemberNum = scLib->rotClusterUnique[aaType][i][0];
							for(int j=0;j<clusterMemberNum;j++){
								selectRotList.push_back(rotList[scLib->rotClusterUnique[aaType][i][j+1]]);
							}
						}
					}
				}
				*/
			}
			cutoff += 5.0;
		}
		this->resMutatorList.push_back(new ResMutator(selectRotList));
		cout << "pos: " << pos << " choice: " << aaChoice  << " rotNum: " << selectRotList.size() << endl;
	}
}


void SeqDesignTemplate::getPositionRanks(ofstream& out){
	ResName rn;
	for(int pos=0;pos<resNum;pos++){
		string tri = rn.intToTri(this->nodeList[pos]->aaType);

		vector<double> aaEnergyList;

		vector<double> s1EnergyList;
		vector<double> s2EnergyList;
		vector<double> s1s2EnergyList;

		vector<double> atomicEnergyList;
		vector<double> vdwEnergyList;

		AAScoreArray sa = this->saList[pos];

		for(int aa=0;aa<20;aa++){

			vector<ResScRotamer*> rotList = this->scLib->rotList[aa];
			int rotNum = rotList.size();

			double s1 = sa.sa[aa];

			double ds1 = sa.sa[aa] - sa.sa[nodeList[pos]->aaType];
			double s2 = 0;
			double ds2 = 0;

			int neighborNum = involvedRpIndex[pos].size();
			ResPairInfo* rp;
			int pairIndex, aaA0, aaA1, aaB0, aaB1;


			for(int i=0;i<neighborNum;i++){
				pairIndex = involvedRpIndex[pos][i];
				rp = rpList[pairIndex];

				//cout << rp->toString() << endl;
				if(rp->indexA == pos){
					aaA0 = aa;
					aaA1 = nodeList[pos]->aaType;
					aaB0 = nodeList[rp->indexB]->conf->aaType;
					s2 += smList[pairIndex].sm[aaA0][aaB0];
					ds2 +=  smList[pairIndex].sm[aaA0][aaB0] - smList[pairIndex].sm[aaA1][aaB0];
					//printf("nat %s Pair: %2d-%2d AA: %s-%s S2: %6.3f\n", tri.c_str(), pos, rp->indexB, rn.intToTri(aaA0).c_str(), rn.intToTri(aaB0).c_str(), smList[pairIndex].sm[aaA0][aaB0]);
				}
				else if(rp->indexB == pos){
					aaA0 = nodeList[rp->indexA]->conf->aaType;
					aaB0 = aa;
					aaB1 = nodeList[pos]->aaType;
					s2 += smList[pairIndex].sm[aaA0][aaB0];
					ds2 +=  smList[pairIndex].sm[aaA0][aaB0] - smList[pairIndex].sm[aaA0][aaB1];
					//printf("nat %s Pair: %2d-%2d AA: %s-%s S2: %6.3f\n", tri.c_str(), rp->indexA, pos, rn.intToTri(aaA0).c_str(), rn.intToTri(aaB0).c_str(), smList[pairIndex].sm[aaA0][aaB0]);
				}
				else {
					cout << "rp: indexA: " << rp->indexA << " indexB: " << rp->indexB << " " << pos << endl;
					cout << "rp error" << endl;
				}
			}
			s1EnergyList.push_back(s1);
			s2EnergyList.push_back(s2);
			s1s2EnergyList.push_back(s1+s2);

			double minEne = 9999.9;
			int minIndex = 0;

			double minEVdw = 999.9;

			for(int k=0;k<rotList.size();k++){

				ResScRotamer* rot = rotList[k];
				nodeList[pos]->confTmp->updateScRotamer(rot);
				double e0 = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->confTmp->scConf->rot) * ec->dp->wtRot;
				double vdwE = 0;
				int neighborNum = nbList[pos].size();
				float meanSai;
				for(int j=0;j<neighborNum;j++){

					meanSai = 0.5*(riList[pos]->sai + riList[nbList[pos][j]]->sai);
					float vdw = ec->getEnergyDesign(nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->conf, meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);

					e0 += vdw;
					vdwE += vdw;

					//if(pos == 3){
					//	ec->printEnergyDesign(pos, nbList[pos][j], nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->conf, meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
					//}
				}

				if(e0 < minEne) {
					minEne = e0;
					minIndex = k;
				}

				if(vdwE < minEVdw){
					minEVdw = vdwE;
				}

			}

			mutEnergy(pos, rotList[minIndex]);
			aaEnergyList.push_back(s1+s2 + minEne);
			atomicEnergyList.push_back(minEne);
			vdwEnergyList.push_back(minEVdw);

		}

		int rank1 = 1;
		int rank2 = 1;
		int rank3 = 1;
		int rank4 = 1;
		int rankV = 1;
		int rank = 1;
		int minEType=-1;
		double minE = 999.9;
		for(int i=0;i<20;i++){

			if(aaEnergyList[i] < aaEnergyList[this->nodeList[pos]->aaType]){
				rank ++;
			}

			if(aaEnergyList[i] < minE){
				minE = aaEnergyList[i];
				minEType = i;
			}


			if(s1EnergyList[i] < s1EnergyList[this->nodeList[pos]->aaType]){
				rank1++;
			}

			if(s2EnergyList[i] < s2EnergyList[this->nodeList[pos]->aaType]){
				rank2++;
			}

			if(s1s2EnergyList[i] < s1s2EnergyList[this->nodeList[pos]->aaType]){
				rank3++;
			}

			if(atomicEnergyList[i] < atomicEnergyList[this->nodeList[pos]->aaType]){
				rank4++;
			}

			if(vdwEnergyList[i] < vdwEnergyList[this->nodeList[pos]->aaType]){
				rankV++;
			}

		}

		out << rn.intToTri(this->nodeList[pos]->aaType) << " " << rn.intToTri(minEType) << " " << rank << " " << rankV << endl;
		//out << rn.intToTri(this->nodeList[pos]->aaType) << " " << rn.intToTri(minEType) << " " << rank1 << " " << rank2 << " " << rank3 << " " <<rankV << " " << rank << endl;;

		//printf("%c %s %s %s %2d %2d %2d %2d %2d\n", this->resList[pos]->chainID, this->resList[pos]->resID.c_str(), rn.intToTri(this->nodeList[pos]->aaType).c_str(), rn.intToTri(minEType).c_str(), rank1, rank2, rank3, rank4, rank);

		//printf("THREE RANKS: %-3d %s %2d %2d %2d\n", pos, rn.intToTri(this->nodeList[pos]->aaType).c_str(), rank1, rank2, rank3);
	}
}

double SeqDesignTemplate::rotEnergyWithBackbone(int pos, ResScRotamer* rot){
	nodeList[pos]->confTmp->updateScRotamer(rot);
	double e = saList[pos].sa[rot->aaType];
	e += scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->confTmp->scConf->rot) * ec->dp->wtRot;
	int neighborNum = nbList[pos].size();
	for(int j=0;j<neighborNum;j++){
		float meanSai = 0.5*(nodeList[nbList[pos][j]]->sai+nodeList[pos]->sai);
		e += ec->getEnergyBbScDesign(nodeList[nbList[pos][j]]->conf, nodeList[pos]->confTmp, meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
	}
	nodeList[pos]->conf->updateScRotamer(nodeList[pos]->conf->scConf->rot);
	return e;
}

double SeqDesignTemplate::totalEnergy(){

	double e = 0;
	double meanSai;
	// rotamer energy
	for(int pos=0;pos<nodeList.size();pos++){
		 e += scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->conf->scConf->rot) * ec->dp->wtRot;
	}

	for(int pos=0;pos<nodeList.size();pos++){
		int neighborNum = nbList[pos].size();
		for(int j=0;j<neighborNum;j++){
			meanSai = 0.5*(nodeList[pos]->sai + nodeList[nbList[pos][j]]->sai);
			e += ec->getEnergyDesign(nodeList[nbList[pos][j]]->conf, nodeList[pos]->conf,meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
		}
	}

	double s1s2=0;
	//s1 energy
	for(int pos=0;pos<nodeList.size();pos++){
		s1s2 += saList[pos].sa[nodeList[pos]->conf->aaType];
	}

	//s2 energy
	ResPairInfo* rp;
	int pairIndex, aaA0, aaA1, aaB0, aaB1;
	for(int pos=0;pos<nodeList.size();pos++){
		int neighborNum = involvedRpIndex[pos].size();
		for(int i=0;i<neighborNum;i++){
			pairIndex = involvedRpIndex[pos][i];
			rp = rpList[pairIndex];
			if(rp->indexA == pos){
				aaA0 = nodeList[pos]->conf->aaType;
				aaB0 = nodeList[rp->indexB]->conf->aaType;
				s1s2 += smList[pairIndex].sm[aaA0][aaB0];
			}
		}
	}

	return e+s1s2;
}

void SeqDesignTemplate::printPolarEnergy() {
	// atomic energy
	float meanSai;
	for(int pos=0;pos<nodeList.size();pos++){
		int neighborNum = nbList[pos].size();
		for(int j=0;j<neighborNum;j++){
			meanSai = 0.5*(nodeList[pos]->sai + nodeList[nbList[pos][j]]->sai);
			if(pos < nbList[pos][j])
			{
				ResConformer* confA = nodeList[pos]->conf;
				ResConformer* confB = nodeList[nbList[pos][j]]->conf;

			}
		}
	}
}

void SeqDesignTemplate::printEnergyDetail() {
	double e = 0;
	float s1, erot;

	for(int pos=0;pos<nodeList.size();pos++){
		string typeA = rn->intToTri(nodeList[pos]->conf->aaType);
		string typeB = rn->intToTri(nodeList[pos]->conf->aaType);
		printf("id: %2d typeBB: %s typeSc: %s\n", pos, typeA.c_str(), typeB.c_str());
	}

	//s1 energy
	for(int pos=0;pos<nodeList.size();pos++){
		s1 = saList[pos].sa[nodeList[pos]->conf->aaType];
		erot = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->conf->scConf->rot) * ec->dp->wtRot;
		printf("pos: %2d s1: %7.3f erot: %7.3f\n", pos, s1, erot);
	}

	//s2 energy
	ResPairInfo* rp;
	int pairIndex, aaA0, aaA1, aaB0, aaB1;

	for(int i=0;i<rpList.size();i++){
		rp = rpList[i];
		int aaA = nodeList[rp->indexA]->conf->aaType;
		int aaB = nodeList[rp->indexB]->conf->aaType;
		double s2 = smList[i].sm[aaA][aaB];
		printf("s2: %2d %2d %7.3f\n", rp->indexA, rp->indexB, s2);
	}

	// atomic energy
	float meanSai;
	for(int pos=0;pos<nodeList.size();pos++){
		int neighborNum = nbList[pos].size();
		for(int j=0;j<neighborNum;j++){
			meanSai = 0.5*(nodeList[pos]->sai + nodeList[nbList[pos][j]]->sai);
			if(pos < nbList[pos][j])
				ec->printEnergyDesign(pos, nbList[pos][j], nodeList[pos]->conf, nodeList[nbList[pos][j]]->conf, meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
		}
	}
}


double SeqDesignTemplate::mutEnergy(int pos, ResScRotamer* rot){

	nodeList[pos]->confTmp->updateScRotamer(rot);

	double e0 = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->conf->scConf->rot) * ec->dp->wtRot;
	int neighborNum = nbList[pos].size();
	float meanSai;
	for(int j=0;j<neighborNum;j++){
		meanSai = 0.5*(nodeList[pos]->sai + nodeList[nbList[pos][j]]->sai);
		if(pos < nbList[pos][j])
			e0 += ec->getEnergyDesign(nodeList[pos]->conf, nodeList[nbList[pos][j]]->conf, meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
		else
			e0 += ec->getEnergyDesign(nodeList[nbList[pos][j]]->conf, nodeList[pos]->conf, meanSai, nodeList[pos]->seqID - nodeList[nbList[pos][j]]->seqID);
	}

	double e1 = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->conf->scConf->rot) * ec->dp->wtRot;
	for(int j=0;j<neighborNum;j++){
		meanSai = 0.5*(nodeList[pos]->sai + nodeList[nbList[pos][j]]->sai);
		if(pos < nbList[pos][j])
			e1 += ec->getEnergyDesign(nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->conf, meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
		else
			e1 += ec->getEnergyDesign(nodeList[nbList[pos][j]]->conf, nodeList[pos]->confTmp, meanSai, nodeList[pos]->seqID - nodeList[nbList[pos][j]]->seqID);
	}

	int aa0 = nodeList[pos]->conf->aaType;
	int aa1 = rot->aaType;

	if(aa0 == aa1)
		return e1 - e0;

	double s0 = saList[pos].sa[aa0];
	double s1 = saList[pos].sa[aa1];

	neighborNum = involvedRpIndex[pos].size();
	ResPairInfo* rp;
	int pairIndex, aaA0, aaA1, aaB0, aaB1;

	for(int i=0;i<neighborNum;i++){
		pairIndex = involvedRpIndex[pos][i];
		rp = rpList[pairIndex];
		if(rp->indexA == pos){
			aaA0 = nodeList[pos]->conf->aaType;
			aaA1 = nodeList[pos]->confTmp->aaType;
			aaB0 = nodeList[rp->indexB]->conf->aaType;
			s0 += smList[pairIndex].sm[aaA0][aaB0];
			s1 += smList[pairIndex].sm[aaA1][aaB0];
		}
		else {
			aaA0 = nodeList[rp->indexA]->conf->aaType;
			aaB0 = nodeList[pos]->conf->aaType;
			aaB1 = nodeList[pos]->confTmp->aaType;
			s0 += smList[pairIndex].sm[aaA0][aaB0];
			s1 += smList[pairIndex].sm[aaA0][aaB1];
		}
	}

	//printf("mutEnergy Pos: %2d S1S2: %5.3f Atomic: %6.3f\n", pos, 2.2*(s1-s0), e1-e0);
	return e1 - e0 + (s1-s0);
}


void SeqDesignTemplate::runMC(){
	int step = resList.size()*this->stepNumFactor;
	cout << "step number: " << step << endl;
	int randPos;
	double deltaE = 0;

	cout << "aaa" << endl;
	cout << "resNum: " << resNum << endl;
	cout << "resMutator: " << this->resMutatorList.size() << endl;


	vector<int> mutPosList;
	for(int i=0;i<resNum;i++){
		cout << "res " <<  i << " rotNum: ";
		int rotNum = this->resMutatorList[i]->rotList.size();
		cout << rotNum << endl;
		if(rotNum > 1){
			mutPosList.push_back(i);
		}
	}

	int mutPosNum = mutPosList.size();

	cout << "mutPosNum: " << mutPosNum << endl;

	double currentE = totalEnergy();

	cout << "init energy: " << currentE << endl;

	for(double T=T0;T>T1;T=T*0.9){
		int ac = 0;
		for(int i=0;i<step;i++){
			randPos = mutPosList[rand()%mutPosNum];
			ResScRotamer* randRot = this->resMutatorList[randPos]->getRandomRotamer();
			deltaE = mutEnergy(randPos, randRot);
			if(deltaE < 0 || rand()*exp(deltaE/T) < RAND_MAX){
				currentE += deltaE;
				nodeList[randPos]->conf->copyValueFrom(nodeList[randPos]->confTmp);

				ac++;
			}
			else {
				nodeList[randPos]->confTmp->copyValueFrom(nodeList[randPos]->conf);
			}
		}
		double totalE = totalEnergy();
		printf("T=%7.4f Accept: %7d CurE: %10.3f TotE: %10.3f idt: %5.3f\n", T, ac, currentE, totalE, identity());
	}
}

double SeqDesignTemplate::identity(){
	int n = initSequence.length();
	string curSeq = getSequence();
	int idt = 0;
	for(int i=0;i<n;i++){
		if(curSeq[i] == initSequence[i])
			idt++;
	}
	return idt*1.0/n;
}

string SeqDesignTemplate::getSequence(){
	char seq[nodeList.size()+1];
	for(int i=0;i<nodeList.size();i++){
		seq[i] = rn->intToSin(nodeList[i]->conf->aaType);
	}
	seq[nodeList.size()] = '\0';
	return string(seq);
}

void SeqDesignTemplate::printPDB(const string& outFile){
	vector<string> bbNames;
	bbNames.push_back("N");
	bbNames.push_back("CA");
	bbNames.push_back("C");
	bbNames.push_back("O");
	ofstream out;
	out.open(outFile.c_str(), ios::out);
	if(!out.is_open()){
		cout << "fail to open file: " << outFile << endl;
		exit(1);
	}
	int atomID = 1;
	char s[100];
	for(int i=0;i<this->nodeList.size();i++){
		char chainID = this->resList[i]->chainID;
		string resID = this->resList[i]->resID;
		char c = resID[resID.length()-1];
		vector<string> scNames;
		this->atLib->getAminoAcidSidechainAtomNames(nodeList[i]->conf->aaType, scNames);

		if(c >='0' && c <= '9'){
			for(int j=0;j<4;j++){
				string name = bbNames[j];
				string type = name.substr(0,1);
				string tri = rn->intToTri(nodeList[i]->conf->aaType);
				XYZ coord = nodeList[i]->conf->bbConf->coords[j+1];
				sprintf(s,"ATOM%7d  %-4s%3s %c%4s    %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,name.c_str(),tri.c_str(),chainID,resID.c_str(),coord[0],coord[1],coord[2],type.c_str());
				out << s << endl;
				atomID++;
			}
			for(int j=0;j<scNames.size();j++){
				string name = scNames.at(j);
				string type = name.substr(0,1);
				string tri = rn->intToTri(nodeList[i]->conf->aaType);
				XYZ coord = nodeList[i]->conf->scConf->coords[j];
				sprintf(s,"ATOM%7d  %-4s%3s %c%4s    %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,name.c_str(),tri.c_str(),chainID,resID.c_str(),coord[0],coord[1],coord[2],type.c_str());
				out << s << endl;
				atomID++;
			}
		}
		else {
			for(int j=0;j<4;j++){
				string name = bbNames[j];
				string type = name.substr(0,1);
				string tri = rn->intToTri(nodeList[i]->conf->aaType);
				XYZ coord = nodeList[i]->conf->bbConf->coords[j+1];
	            sprintf(s,"ATOM%7d  %-4s%3s %c%5s   %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,name.c_str(),tri.c_str(),chainID,resID.c_str(),coord[0],coord[1],coord[2],type.c_str());
				out << s << endl;
				atomID++;
			}
			for(int j=0;j<scNames.size();j++){
				string name = scNames.at(j);
				string type = name.substr(0,1);
				string tri = rn->intToTri(nodeList[i]->conf->aaType);
				XYZ coord = nodeList[i]->conf->scConf->coords[j];
				sprintf(s,"ATOM%7d  %-4s%3s %c%5s   %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,name.c_str(),tri.c_str(),chainID,resID.c_str(),coord[0],coord[1],coord[2],type.c_str());
				out << s << endl;
				atomID++;
			}
		}
	}
	out.close();
}

SeqDesignTemplate::~SeqDesignTemplate() {
	// TODO Auto-generated destructor stub
	delete this->rn;
	delete this->pdb;
	for(int i=0;i<nodeList.size();i++){
		delete nodeList[i];
	}
	for(int i=0;i<riList.size();i++) {
		delete riList[i];
	}
	for(int i=0;i<rpList.size();i++) {
		delete rpList[i];
	}
}

} /* namespace NSPmodel */
