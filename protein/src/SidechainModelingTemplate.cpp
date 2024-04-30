/*
 * SidechainModelingTemplate.cpp
 *
 */

#include "protein/SidechainModelingTemplate.h"

namespace NSPprotein {

SidechainModelingTemplate::SidechainModelingTemplate(const string& pdbFile) {
	// TODO Auto-generated constructor stub
	this->rn = new ResName();
	this->bbLib = new ResBBRotamerLib();
	this->scLib = new ResScRotamerLib();
	this->atLib = new AtomLib();

	this->pdb = new PDB(pdbFile, "xxxx");
	vector<Residue*> pdbResList = pdb->getResList();
	for(int i=0;i<pdbResList.size();i++) {
		Residue* res = pdbResList[i];
		if(res->hasThreeCoreAtoms() && res->getAtom("O") != NULL)
		{
			resList.push_back(res);
			if(res->sidechainComplete(atLib))
				resValidList.push_back(true);
			else
				resValidList.push_back(false);
		}
	}

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

	nodeList.push_back(new ResNode(resList[0], seqIDList[0], 0, atLib, bbLib, scLib));
	for(int i=1;i<resList.size();i++) {
		if(seqIDList[i] == seqIDList[i-1]+1){
			nodeList.push_back(new ResNode(resList[i-1], resList[i], seqIDList[i], i,  atLib, bbLib, scLib));
		}
		else {
			nodeList.push_back(new ResNode(resList[i], seqIDList[i], i,  atLib, bbLib, scLib));
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
	srand(time(0));
}


double SidechainModelingTemplate::mutEnergy(int pos, ResScRotamer* rot, EnergyCalculator* ec){

	nodeList[pos]->confTmp->updateScRotamer(rot);

	double e0 = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->conf->scConf->rot) * ec->dp->wtRot;
	int neighborNum = nbList[pos].size();

	for(int j=0;j<neighborNum;j++){
		e0 += ec->getEnergyScSc(nodeList[pos]->conf, nodeList[nbList[pos][j]]->conf, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
		e0 += ec->getEnergyBbSc(nodeList[nbList[pos][j]]->conf, nodeList[pos]->conf, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
	}

	double e1 = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->confTmp->scConf->rot) * ec->pa->wtRot;
	for(int j=0;j<neighborNum;j++){
		e1 += ec->getEnergyScSc(nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->conf, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
		e1 += ec->getEnergyBbSc(nodeList[nbList[pos][j]]->conf, nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
	}
	return e1 - e0;
}

void SidechainModelingTemplate::acceptMutation(int pos, ResScRotamer* rot) {
	nodeList[pos]->conf->updateScRotamer(rot);
	nodeList[pos]->confTmp->updateScRotamer(rot);
}

ResScRotamer* SidechainModelingTemplate::findBestRotamer(int pos, EnergyCalculator* ec){
	vector<ResScRotamer*> rotList = this->scLib->rotList[this->nodeList[pos]->aaType];
	int rotNum = rotList.size();
	int neighborNum = nbList[pos].size();
	double minEne = DBL_MAX;
	int minID = 0;

	int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(this->nodeList[pos]->aaType);

	if(clusterNum == 0) {
		for(int i=0;i<rotNum;i++){
			nodeList[pos]->confTmp->updateScRotamer(rotList[i]);
			double e = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->confTmp->scConf->rot) * ec->pa->wtRot;
			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->conf, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
				e += ec->getEnergyBbSc(nodeList[nbList[pos][j]]->conf, nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
			}

			if(e < minEne){
				minEne = e;
				minID = i;
			}
		}
		nodeList[pos]->confTmp->updateScRotamer(rotList[minID]);
	}
	else {
		//search level1
		int level1MinID = 0;
		for(int i=0;i<clusterNum;i++){
			nodeList[pos]->confTmp->updateScRotamer(rotList[scLib->rotCluster[nodeList[pos]->aaType][i][1]]);
			double e = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->confTmp->scConf->rot) * ec->pa->wtRot;
			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->conf, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
				e += ec->getEnergyBbSc(nodeList[nbList[pos][j]]->conf, nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
			}

			if(e < minEne){
				minEne = e;
				level1MinID = i;
				minID = scLib->rotCluster[nodeList[pos]->aaType][i][1];
			}
		}
		//search level2
		int clusterMemberNum = scLib->rotCluster[nodeList[pos]->aaType][level1MinID][0];
		for(int i=0;i<clusterMemberNum;i++){
			nodeList[pos]->confTmp->updateScRotamer(rotList[scLib->rotCluster[nodeList[pos]->aaType][level1MinID][i+1]]);
			double e = scLib->getEnergy(nodeList[pos]->conf->bbConf->rot->index1K, nodeList[pos]->confTmp->scConf->rot) * ec->pa->wtRot;
			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->conf, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
				e += ec->getEnergyBbSc(nodeList[nbList[pos][j]]->conf, nodeList[pos]->confTmp, nodeList[nbList[pos][j]]->seqID - nodeList[pos]->seqID);
			}
			if(e < minEne){
				minEne = e;
				minID = scLib->rotCluster[nodeList[pos]->aaType][level1MinID][i+1];
			}
		}
	}
	return rotList[minID];
}

void SidechainModelingTemplate::runMC(EnergyCalculator* ec) {

	vector<int> randPosIndex;
	for(int i=0;i<nodeList.size();i++){
		int aa = nodeList[i]->aaType;
		int rot = this->scLib->rotList[aa].size();
		if(rot > 4000) rot = 4000;
		rot = (int)sqrt(rot);
		for(int j=0;j<rot;j++){
			randPosIndex.push_back(i);
		}
	}
	int posIndexSize = randPosIndex.size();

	/*
	 * random init
	 */
	for(int i=0;i<nodeList.size();i++){
		int randRotChoice = rand()%scLib->rotNum[nodeList[i]->aaType];
		ResScRotamer* rot = scLib->rotList[nodeList[i]->aaType][randRotChoice];
		nodeList[i]->conf->updateScRotamer(rot);
		nodeList[i]->confTmp->updateScRotamer(rot);
	}

	double T0 = 40.0;
	double T1 = 0.04;
	int step = resList.size() * 400;
	int randPos;
	int currentChoice;
	int randRotChoice;
	double deltaE = 0;

	for(double T = T0;T>T1;T=T*0.9){

		int ac = 0;
		for(int i=0;i<step;i++){
			randPos = randPosIndex[rand()%posIndexSize];
			currentChoice = nodeList[randPos]->conf->scConf->rot->rotID;
			randRotChoice = rand()%scLib->rotNum[nodeList[randPos]->aaType];

			if(randRotChoice == currentChoice) continue;
			deltaE = mutEnergy(randPos, scLib->rotList[nodeList[randPos]->aaType][randRotChoice], ec);

			if(deltaE < 0 || rand()*exp(deltaE/T) < RAND_MAX)
			{
				nodeList[randPos]->conf->copyValueFrom(nodeList[randPos]->confTmp);
				ac++;
			}
			else
			{
				nodeList[randPos]->confTmp->copyValueFrom(nodeList[randPos]->conf);
			}
		}
		printf("T=%7.4f Accept: %7d / %7d rms: %6.4f %6.4f\n", T, ac,step, rms1(), rms2());
	}

	for(int i=0;i<nodeList.size();i++){
		ResScRotamer* bestRot = findBestRotamer(i, ec);
		acceptMutation(i, bestRot);
	}

	printf("rms: %6.4f %6.4f\n", rms1(), rms2());

	for(int i=0;i<nodeList.size();i++){
		ResScRotamer* bestRot = findBestRotamer(i, ec);
		acceptMutation(i, bestRot);
	}
	printf("rms: %6.4f %6.4f\n", rms1(), rms2());
}

void SidechainModelingTemplate::printPDB(const string& outfile){
	for(int i=0;i<resList.size();i++){
		Residue* res = resList[i];
		vector<string> scNames;
		atLib->getAminoAcidSidechainAtomNames(res->intName, scNames);
		ResConformer* conf = nodeList[i]->conf;
		int n = scNames.size();
		for(int j=0;j<n;j++){
			string name = scNames.at(j);
			if(res->getAtom(name) != NULL){
				res->getAtom(name)->setCoord(conf->scConf->coords[j]);
			}
			else {
				res->addAtom(new Atom(name, conf->scConf->coords[j]));
			}
		}
	}

	ProteinChain* pc = new ProteinChain();
	for(int i=0;i<resList.size();i++){
		pc->addResidue(resList[i]);
	}
	ofstream out;
	out.open(outfile, ios::out);
	pc->printPDBFormat(out, 1);
	out.close();

}

double SidechainModelingTemplate::rms1(){
	double rr = 0.0;
	int n = 0;
	ResScRotamer* nat;
	ResScRotamer* pred;
	double d;

	for(int i=0;i<resList.size();i++){
		if(resValidList[i]) {
			int aa = nodeList[i]->aaType;
			nat = initRotList[i];
			pred = nodeList[i]->conf->scConf->rot;
			d = nat->distanceTo(pred);
			n++;
			rr += d*d;
		}
	}
	return sqrt(rr/n);
}

double SidechainModelingTemplate::rms2(){
	double rr = 0.0;
	int n = 0;

	for(int i=0;i<resList.size();i++){
		vector<string> names;
		atLib->getAminoAcidSidechainAtomNames(resList[i]->intName, names);
		int m = names.size();
		double dd = 0;

		for(int j=0;j<m;j++){
			Atom* a = resList[i]->getAtom(names.at(j));
			XYZ coordPred = nodeList[i]->conf->scConf->coords[j];
			if(a != NULL){
				dd += a->coord.squaredDistance(nodeList[i]->conf->scConf->coords[j]);
			}
		}
		if(m > 0) {
			double d = sqrt(dd/m);
			n++;
			rr += d*d;
		}
	}
	return sqrt(rr/n);
}

void SidechainModelingTemplate::printRMS(ofstream& out){

	char xx[200];
	for(int i=0;i<resList.size();i++){
		vector<string> names;
		atLib->getAminoAcidSidechainAtomNames(resList[i]->intName, names);
		int m = names.size();
		double dd = 0;

		for(int j=0;j<m;j++){
			Atom* a = resList[i]->getAtom(names.at(j));
			XYZ coordPred = nodeList[i]->conf->scConf->coords[j];
			if(a != NULL){
				dd += a->coord.squaredDistance(nodeList[i]->conf->scConf->coords[j]);
			}
		}
		if(m > 0) {
			double rms = sqrt(dd/m);
			sprintf(xx, "%s %6.4f", resList[i]->triName.c_str(), rms);
			out << string(xx) << endl;
		}
	}
}


void SidechainModelingTemplate::printEnergy() {

}

SidechainModelingTemplate::~SidechainModelingTemplate() {
	// TODO Auto-generated destructor stub
	delete rn;
	delete bbLib;
	delete scLib;
	delete atLib;
	delete pdb;
	for(int i=0;i<nodeList.size();i++){
		delete nodeList[i];
	}
}

} /* namespace NSPmodel */
