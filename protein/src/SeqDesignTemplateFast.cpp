/*
 * SeqDesignTemplateFast.cpp
 *
 *  Created on: 2022��7��16��
 *      Author: pengx
 */

#include "protein/SeqDesignTemplateFast.h"

namespace NSPprotein {

SeqDesignTemplateFast::SeqDesignTemplateFast(const string& inputFile, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib, AtomLib* atLib, EnergyCalculator* ec, ProS1S2Energy* eS1S2) {
	// TODO Auto-generated constructor stub
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
			if(cbList[i].distance(cbList[j]) < 12)
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
			if(cbA.distance(cbB) >= 12.0) continue;
			ResPairInfo* rp = new ResPairInfo(i, j, resList[i]->intName, resList[j]->intName, riList[i]->ss, riList[j]->ss, riList[i]->sai, riList[j]->sai, cm, seqIDList[j] - seqIDList[i]);
			involvedRpIndex[i].push_back(rpIndex);
			involvedRpIndex[j].push_back(rpIndex);
			rpIndex++;
			this->rpList.push_back(rp);
		}
	}

	cout << "load S1S2" << endl;
	loadS1S2();

	double profWeight = 1.0;
	if(input->specifiedOption("PRIFILE_WEIGHT"))
		profWeight = atof(input->getValue("PRIFILE_WEIGHT").c_str());
	if(input->specifiedOption("PROFILE"))
		loadMSAProfile(input->getValue("PROFILE"), profWeight);

	cout << "update rotamer choice" << endl;


	updateRotChoice(resFile);

	cout << "init rotSequence" << endl;
	this->rotSeq = new RotSequence(resMutatorList);

	delete rsp;
	delete ssa;
	delete input;
}

SeqDesignTemplateFast::SeqDesignTemplateFast(const string& inputFile, ResBBRotamerLib* bbLib, ResScRotamerLibMini* scLib, AtomLib* atLib, EnergyCalculator* ec, ProS1S2Energy* eS1S2, const string& s2File) {
	// TODO Auto-generated constructor stub
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
			if(cbList[i].distance(cbList[j]) < 12)
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
			if(cbA.distance(cbB) >= 12.0) continue;
			ResPairInfo* rp = new ResPairInfo(i, j, resList[i]->intName, resList[j]->intName, riList[i]->ss, riList[j]->ss, riList[i]->sai, riList[j]->sai, cm, seqIDList[j] - seqIDList[i]);
			involvedRpIndex[i].push_back(rpIndex);
			involvedRpIndex[j].push_back(rpIndex);
			rpIndex++;
			this->rpList.push_back(rp);
		}
	}

	cout << "load S1S2" << endl;
	loadS1S2(s2File);

	double profWeight = 1.0;
	if(input->specifiedOption("PRIFILE_WEIGHT"))
		profWeight = atof(input->getValue("PRIFILE_WEIGHT").c_str());
	if(input->specifiedOption("PROFILE"))
		loadMSAProfile(input->getValue("PROFILE"), profWeight);

	cout << "update rotamer choice" << endl;
	updateRotChoice(resFile);

	cout << "init rotSequence" << endl;
	this->rotSeq = new RotSequence(resMutatorList);

	delete rsp;
	delete ssa;
	delete input;
}

void SeqDesignTemplateFast::printInvolvedPairInfo(){
	for(int i=0;i<involvedRpIndex.size();i++){
		cout << i << " ";
		for(int j=0;j<involvedRpIndex[i].size();j++){
			cout << involvedRpIndex[i][j] << " ";
		}
		cout << endl;
	}
}

void SeqDesignTemplateFast::loadS1S2(){
	saList.clear();
	smList.clear();
	for(int i=0;i<riList.size();i++){
		this->saList.push_back(eS1S2->getS1(riList[i]));
	}

	for(int i=0;i<rpList.size();i++){
		cout << i << endl;
		if(rpList[i]->cbDistance > 8.0)
			this->smList.push_back(AAScoreMatrix());
		else
			this->smList.push_back(eS1S2->getS2(rpList[i]));
	}

	cout << "pairNum: " << rpList.size() << endl;
	cout << "smNum: " << smList.size() << endl;
}

void SeqDesignTemplateFast::loadS1S2(const string& s2File){
	saList.clear();
	smList.clear();

	for(int i=0;i<riList.size();i++){
		this->saList.push_back(eS1S2->getS1(riList[i]));
	}

	ifstream f;
	f.open(s2File.c_str(), ios::in);
	string s;
	vector<string> spt;
	double sm[400];
	for(int i=0;i<rpList.size();i++){

		getline(f, s);
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
	f.close();
}

void SeqDesignTemplateFast::printS2(const string& outfile){

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

void SeqDesignTemplateFast::printBackboneInfo(){

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

void SeqDesignTemplateFast::printResInfo(){
	for(int i=0;i<resNum;i++){
		char ss = riList[i]->ss;
		double sai = riList[i]->sai;
		printf("%c %5.3f\n", ss, sai);
	}
}

void SeqDesignTemplateFast::loadEAEM(){
	//load energy array
	float s1, s2, eRot, eRef, eVdw;
	ResConformer* confA;
	ResConformer* confB;
	for(int pos=0;pos<resNum;pos++){
		int cn = resMutatorList[pos]->choiceNum;
		EnergyArray* ea = new EnergyArray(resMutatorList[pos]->choiceNum);
		for(int j=0;j<cn;j++){
			s1 = saList[pos].sa[resMutatorList[pos]->rotList[j]->aaType];
			eRot = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, resMutatorList[pos]->rotList[j]) * ec->dp->wtRot;
			eRef = ec->dp->ref[riList[pos]->ssType*20 + resMutatorList[pos]->rotList[j]->aaType];
			ea->setEnergy(j, s1+eRot+eRef);
		}
		this->eaList.push_back(ea);
	}

	for(int i=0;i<rpList.size();i++){
		int posA = rpList[i]->indexA;
		int posB = rpList[i]->indexB;
		int cnA = resMutatorList[posA]->choiceNum;
		int cnB = resMutatorList[posB]->choiceNum;
		double meanSai = (rpList[i]->saiA + rpList[i]->saiB)*0.5;

		EnergyMatrix* em = new EnergyMatrix(cnA, cnB, posA, posB);
		for(int j=0;j<cnA;j++){
			confA = resConfList[posA][j];
			for(int k=0;k<cnB;k++){
				confB = resConfList[posB][k];
				s2 = smList[i].sm[confA->aaType][confB->aaType];
				eVdw = ec->getEnergyDesign(confA, confB, meanSai, nodeList[posA]->seqID - nodeList[posB]->seqID);
				em->em[j*em->choiceBNum + k] = s2 + eVdw;
			}
		}
		this->emList.push_back(em);
	}
}

void SeqDesignTemplateFast::loadRnaAsLigand(RNAPDB* rpdb){

	RNABaseName rn;
	AtomLib* atLib = new AtomLib();
	vector<XYZ> atomCoords;
	vector<AtomProperty*> apList;
	char ss[20];
	vector<RNABase*> baseList = rpdb->getBaseList();
	for(int i=0;i<baseList.size();i++){
		RNABase* b = baseList[i];
		for(int j=0;j<b->getAtomList()->size();j++){
			Atom* a = b->getAtomList()->at(j);
			atomCoords.push_back(a->coord);
			string name = a->getName();
			sprintf(ss, "%s-%s", rn.intToString(b->baseTypeInt).c_str(), name.c_str());
			string uniqueName = string(ss);
			int uniqueID = atLib->uniqueNameToUniqueID[uniqueName];
			apList.push_back(atLib->getAtomProperty(uniqueID));
		}
	}
	int nRnaAtoms = atomCoords.size();

	float s1, s2, eRot, eRef, eVdw;
	ResConformer* confA;
	float d, d0;

	for(int pos=0;pos<resNum;pos++){
		int cn = resMutatorList[pos]->choiceNum;

		for(int j=0;j<cn;j++){
			eVdw = 0;
			confA = resConfList[pos][j];
			//sc - rna
			for(int k=0;k<confA->scConf->atomNum;k++){
				XYZ t = confA->scConf->coords[k];
				int uniqueID = confA->scConf->rot->uniqueIDs[k];
				AtomProperty* ap = this->atLib->getAtomProperty(uniqueID);
				for(int m=0;m<nRnaAtoms;m++){
					d = t.distance(atomCoords[m]);
					if(d > 8.0) continue;
					d0 = apList[m]->vdwRadius + ap->vdwRadius;
					eVdw += this->ec->et->etAt->vdwEnergy(d, d0, 0.1, -1.0, 2.8);
				}
			}

			this->eaList[pos]->setEnergy(j, this->eaList[pos]->getEnergy(j)+eVdw);
		}
	}
}

void SeqDesignTemplateFast::printPairInfo(){
	for(int pid=0;pid<rpList.size();pid++){
		int indexA = rpList[pid]->indexA;
		int indexB = rpList[pid]->indexB;
		ResPairInfo* rp = this->rpList[pid];
		cout << rp->toString() << endl;
	}
}

void SeqDesignTemplateFast::loadMSAProfile(const string& file, double weight){
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

void SeqDesignTemplateFast::updateRotChoice(const string& resFile){
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
		vector<ResScRotamer*> selectRotList;

		if(fixed.find(pos) != fixed.end()){
			selectRotList.push_back(this->initRotList[pos]);
			this->resMutatorList.push_back(new ResMutator(selectRotList));
			cout << "pos: " << pos << " choice: " << aaChoice  << " rotNum: " << selectRotList.size() << endl;
			continue;
		}

		for(int i=0;i<aaChoice.length();i++){
			aaType = rn->sinToInt(aaChoice[i]);
			if(aaType < 0 || aaType > 19){
				cout << "invalid aa choice: " << aaChoice << endl;
				exit(0);
			}
			vector<ResScRotamer*> rotList = this->scLib->rotList[aaType];
			int rotNum = rotList.size();
			for(int i=0;i<rotNum;i++){
				if(rotamerValid(pos, rotList[i]))
					selectRotList.push_back(rotList[i]);
			}
		}

		if(selectRotList.size() == 0) {
			for(int i=0;i<aaChoice.length();i++){
				aaType = rn->sinToInt(aaChoice[i]);
				if(aaType < 0 || aaType > 19){
					cout << "invalid aa choice: " << aaChoice << endl;
					exit(0);
				}
				vector<ResScRotamer*> rotList = this->scLib->rotList[aaType];
				int rotNum = rotList.size();
				for(int i=0;i<rotNum;i++){
					selectRotList.push_back(rotList[i]);
				}
			}
		}

		this->resMutatorList.push_back(new ResMutator(selectRotList));
		cout << "pos: " << pos << " choice: " << aaChoice  << " rotNum: " << selectRotList.size() << endl;
	}

	for(int pos=0;pos<resNum;pos++){
		vector<ResConformer*> confList;
		int choiceNum = resMutatorList[pos]->choiceNum;
		for(int i=0;i<choiceNum;i++){
			ResScRotamer* rot = resMutatorList[pos]->rotList[i];
			ResConformer* conf = new ResConformer(rot->aaType);
			conf->init(rot, nodeList[pos]->conf->bbConf->rot, nodeList[pos]->cs2);
			confList.push_back(conf);
		}
		resConfList.push_back(confList);
	}
}

void SeqDesignTemplateFast::printDetailEnergy(){
	ResConformer* confA;
	ResConformer* confB;
	double s2, eVdw;

	for(int i=0;i<rpList.size();i++){
		int posA = rpList[i]->indexA;
		int posB = rpList[i]->indexB;
		confA = nodeList[posA]->conf;
		confB = nodeList[posB]->conf;
		double meanSai = (rpList[i]->saiA + rpList[i]->saiB)*0.5;
		s2 = smList[i].sm[confA->aaType][confB->aaType];
		eVdw = ec->getEnergyDesign(confA, confB, meanSai, nodeList[posA]->seqID - nodeList[posB]->seqID);
		printf("confA %3d confB %3d s2: %6.3f vdw: %7.3f\n", posA, posB, s2, eVdw);
		ec->printEnergyDesign(posA, posB, confA, confB, meanSai, nodeList[posA]->seqID - nodeList[posB]->seqID);
	}
}


void SeqDesignTemplateFast::getPositionRanks(ofstream& out){
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
		for(int i=0;i<20;i++){
			string tri = rn.intToTri(i);
			double e = aaEnergyList[i];
			out << tri << " " << e << endl;
		}
		//out << rn.intToTri(this->nodeList[pos]->aaType) << " " << rn.intToTri(minEType) << " " << rank1 << " " << rank2 << " " << rank3 << " " <<rankV << " " << rank << endl;;

		//printf("%c %s %s %s %2d %2d %2d %2d %2d\n", this->resList[pos]->chainID, this->resList[pos]->resID.c_str(), rn.intToTri(this->nodeList[pos]->aaType).c_str(), rn.intToTri(minEType).c_str(), rank1, rank2, rank3, rank4, rank);

		//printf("THREE RANKS: %-3d %s %2d %2d %2d\n", pos, rn.intToTri(this->nodeList[pos]->aaType).c_str(), rank1, rank2, rank3);
	}
}

bool SeqDesignTemplateFast::rotamerValid(int pos, ResScRotamer* rot){
	bool valid = true;

	nodeList[pos]->confTmp->updateScRotamer(rot);
	double e = saList[pos].sa[rot->aaType] - saList[pos].saBg[rot->aaType];

	if(e > 1.0 && rot->aaType != 5)
		valid = false;


	nodeList[pos]->confTmp->updateScRotamer(rot);
	double eRot = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->confTmp->scConf->rot);

	if(eRot > 4.0)
		valid = false;

	double clash = 0;
	int neighborNum = nbList[pos].size();
	for(int j=0;j<neighborNum;j++){
		float meanSai = 0.5*(nodeList[nbList[pos][j]]->sai+nodeList[pos]->sai);
		clash += ec->getEnergyBbScDesign(nodeList[nbList[pos][j]]->conf, nodeList[pos]->confTmp, meanSai, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
	}
	nodeList[pos]->confTmp->updateScRotamer(nodeList[pos]->conf->scConf->rot);
	if(clash > 5.0)
		valid = false;

	return valid;
}

double SeqDesignTemplateFast::totalEnergy(){

	double e = 0;
	for(int pos=0;pos<resNum;pos++){
		e += eaList[pos]->ea[this->rotSeq->rotChoice[pos]];
	}
	for(int i=0;i<emList.size();i++){
		int choiceA = rotSeq->rotChoice[emList[i]->posA];
		int choiceB = rotSeq->rotChoice[emList[i]->posB];
		e += emList[i]->em[choiceA * emList[i]->choiceBNum + choiceB];
	}
	return e;
}

void SeqDesignTemplateFast::printNatEnergy(){

	cout << "res num: " << this->resNum << endl;
	cout << "ea num: " << this->eaList.size() << endl;
	int* natChoice = new int[this->resNum];
	for(int i=0;i<resNum;i++){
		int rotNum = this->resMutatorList[i]->rotList.size();
		ResScRotamer* scRot = this->initRotList[i];
		double minD = 99.9;
		int minID = 0;
		for(int j=0;j<rotNum;j++){
			ResScRotamer* rot = this->resMutatorList[i]->rotList[j];
			double d = scRot->distanceTo(rot);
			if(d < minD){
				minD = d;
				minID = j;
			}
		}
		cout << "res: " << i << " " << minD << endl;
		natChoice[i] = minID;
	}

	for(int i=0;i<resNum;i++){
		double e1 = eaList[i]->getEnergy(natChoice[i]);
		printf("%-2d %8.3f\n", i, e1);
	}

	cout << "rp num: " << rpList.size() << endl;
	cout << "em num: " << emList.size() << endl;

	for(int i=0;i<rpList.size();i++){
		ResPairInfo* rp = rpList[i];
		int posA = rp->indexA;
		int posB = rp->indexB;
		double e2 = emList[i]->getEnergy(natChoice[posA], natChoice[posB]);
		printf("%-2d %-2d %8.3f\n", posA, posB, e2);
	}
}

double SeqDesignTemplateFast::resInvolvedEnergy(int pos, int choice){
	double e1 = eaList[pos]->ea[choice];
	double e2 = 0;
	EnergyMatrix* em;
	int posA, posB;
	int pairNum = involvedRpIndex[pos].size();
	for(int i=0;i<pairNum;i++){
		int pairID = involvedRpIndex[pos][i];
		em = emList[pairID];
		posA = this->rpList[pairID]->indexA;
		posB = this->rpList[pairID]->indexB;
		if(pos == posA){
			e2 += em->em[choice * em->choiceBNum + this->rotSeq->rotChoice[posB]];
		}
		else if(pos == posB){
			e2 += em->em[this->rotSeq->rotChoice[posA] * em->choiceBNum + choice];
		}
		else {
			cerr << "posA: " << posA << " posB: " << posB << " mutPos: " << pos << " resInvolved Energy error" << endl;
			exit(1);
		}
	}
	return e1+e2;
}

double SeqDesignTemplateFast::mutEnergy(int pos, int choice){
	int oldChoice = rotSeq->rotChoice[pos];
	if(choice == oldChoice) return 0.0;
	return resInvolvedEnergy(pos, choice) - resInvolvedEnergy(pos, oldChoice);
}

double SeqDesignTemplateFast::mutEnergy(int pos, ResScRotamer* rot){

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
	return e1 - e0 + (s1-s0);
}


void SeqDesignTemplateFast::designMC(){

	if(this->eaList.size() == 0 && this->emList.size()==0)
		loadEAEM();

	int step = resList.size()*3000;
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
	rotSeq->setRandomChoice();
	double currentE = totalEnergy();
	cout << "init energy: " << currentE << endl;

	int choiceNum, randChoice;
	for(double T=T0;T>T1;T=T*0.9){
		int ac = 0;
		for(int i=0;i<step;i++){
			randPos = mutPosList[rand()%mutPosNum];
			choiceNum = resMutatorList[randPos]->choiceNum;
			randChoice = rand()%choiceNum;
			deltaE = mutEnergy(randPos, randChoice);
			if(deltaE < 0 || rand()*exp(deltaE/T) < RAND_MAX){
				currentE += deltaE;
				rotSeq->rotChoice[randPos] = randChoice;
				ac++;
			}
		}
		double totalE = totalEnergy();
		printf("T=%7.4f Accept: %7d CurE: %10.3f TotE: %10.3f\n", T, ac, currentE, totalE);
	}
}


string SeqDesignTemplateFast::getDesignSequence() {
	return this->rotSeq->toAASequence();
}

void SeqDesignTemplateFast::printDesignPDB(const string& outFile){
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

		int choice = this->rotSeq->getChoice(i);
		ResScRotamer* scRot = this->resMutatorList[i]->rotList[choice];

		nodeList[i]->conf->updateScRotamer(scRot);
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

SeqDesignTemplateFast::~SeqDesignTemplateFast() {
	// TODO Auto-generated destructor stub
	cout << "delete rn" << endl;
	delete this->rn;
	cout << "delete pdb" << endl;
	delete this->pdb;

	cout << "delete node" << endl;
	for(int i=0;i<nodeList.size();i++){
		delete nodeList[i];
	}

	cout << "delete ri" << endl;
	for(int i=0;i<riList.size();i++) {
		delete riList[i];
	}

	cout << "delete rp" << endl;
	for(int i=0;i<rpList.size();i++) {
		delete rpList[i];
	}

	cout << "delete resConf" << endl;
	for(int i=0;i<resNum;i++){
		for(int j=0;j<resConfList[i].size();j++){
			delete resConfList[i][j];
		}
	}

	cout << "delete rotSeq" << endl;
	delete rotSeq;

	for(int i=0;i<eaList.size();i++){
		delete eaList[i];
	}

	for(int i=0;i<emList.size();i++){
		delete emList[i];
	}
}

} /* namespace NSPmodel */
