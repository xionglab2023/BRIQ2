/*
 * SingleSidechainModelingTemplate.cpp
 *
 */

#include "protein/SingleSidechainModelingTemplate.h"

namespace NSPprotein {

SingleSidechainModelingTemplate::SingleSidechainModelingTemplate(vector<string>& lines, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, ResName* rn) {
	// TODO Auto-generated constructor stub
	this->rn = rn;
	this->bbLib = bbLib;
	this->scLib = scLib;
	this->preAAType = 0;

	this->targetNode = new ResNode(lines[0], 0,  rn, bbLib, scLib);
	vector<string> spt;

	for(int i=1;i<lines.size();i++){
		splitString(lines[i], " ", &spt);
		int sep = atoi(spt[1].c_str());
		this->sepList.push_back(sep);

		this->neighborNodes.push_back(new ResNode(lines[i], i, rn, bbLib, scLib));
	}

	for(int i=0;i<neighborNodes.size();i++){
		if(sepList[i] == -1){
			preAAType = neighborNodes[i]->aaType;
		}
	}
}

double SingleSidechainModelingTemplate::bestRotamerRMSD(EnergyCalculator* ec){

	vector<ResScRotamer*> rotList = this->scLib->rotList[this->targetNode->aaType];
	int rotNum = rotList.size();
	int neighborNum = neighborNodes.size();

	double minEne = DBL_MAX;
	int minID = 0;

	for(int i=0;i<rotNum;i++){
		targetNode->confTmp->updateScRotamer(rotList[i]);
		double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;
		for(int j=0;j<neighborNum;j++){
			e += ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
			e += ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);
		}

		if(e < minEne){
			minEne = e;
			minID = i;
		}
	}

	targetNode->confTmp->updateScRotamer(rotList[minID]);

	/*
	 * print final energy
	 */

	/*
	{
		double e = scLib->getEnergy(targetNode->conBB->rot->index1K, targetNode->conScTmp->rot) * ec->para3->wtRot;
		printf("eRot: %7.3f\n", e);
		for(int j=0;j<neighborNum;j++){
			double eScSc = ec->getEnergyScSc(targetNode->conScTmp, neighborNodes[j]->conSc, neighborNodes[j]->seqID);
			double ePolar = ec->getEnergyScScPolar(targetNode->conScTmp, neighborNodes[j]->conSc, neighborNodes[j]->seqID);
			double eBbSc = ec->getEnergyBbSc(neighborNodes[j]->conBB, targetNode->conScTmp, -neighborNodes[j]->seqID);
			printf("%-3d %7.3f %7.3f %7.3f\n", j, eScSc, ePolar, eBbSc);
		}
	}
	*/

	double rms = targetNode->conf->scConf->rot->distanceTo(rotList[minID]);
	return rms;
}

double SingleSidechainModelingTemplate::bestRotamerRMSD2(EnergyCalculator* ec){

	//search in two level rotamer libraries
	vector<ResScRotamer*> rotList = this->scLib->rotList[this->targetNode->aaType];
	int rotNum = rotList.size();
	int neighborNum = neighborNodes.size();
	double minEne = DBL_MAX;
	int minID = 0;

	int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(this->targetNode->aaType);
	if(clusterNum == 0) {
		for(int i=0;i<rotNum;i++){
			targetNode->confTmp->updateScRotamer(rotList[i]);
			double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;

			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
				e += ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);
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
			targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[targetNode->aaType][i][1]]);
			double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;
			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
				e += ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);
			}

			if(e < minEne){
				minEne = e;
				level1MinID = i;
				minID = scLib->rotCluster[targetNode->aaType][i][1];
			}
		}
		//search level2
		int clusterMemberNum = scLib->rotCluster[targetNode->aaType][level1MinID][0];
		for(int i=0;i<clusterMemberNum;i++){
			targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[targetNode->aaType][level1MinID][i+1]]);
			double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;
			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
				e += ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);
			}
			if(e < minEne){
				minEne = e;
				minID = scLib->rotCluster[targetNode->aaType][level1MinID][i+1];
			}
		}
		targetNode->confTmp->updateScRotamer(rotList[minID]);
	}

	double rms = targetNode->conf->scConf->rot->distanceTo(rotList[minID]);
	return rms;
}

void SingleSidechainModelingTemplate::printNativeEnergy(EnergyCalculator* ec){

	vector<ResScRotamer*> rotList = this->scLib->rotList[this->targetNode->aaType];
	int rotNum = rotList.size();
	int neighborNum = neighborNodes.size();
	double minEne = DBL_MAX;
	int minID = 0;

	int clusterNum = this->scLib->getAminoAcidRotamerClusterNum(this->targetNode->aaType);
	if(clusterNum == 0) {
		for(int i=0;i<rotNum;i++){
			targetNode->confTmp->updateScRotamer(rotList[i]);
			double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;
			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
				e += ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);
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
			targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[targetNode->aaType][i][1]]);
			double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;

			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
				e += ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);
			}

			if(e < minEne){
				minEne = e;
				level1MinID = i;
				minID = scLib->rotCluster[targetNode->aaType][i][1];
			}
		}
		//search level2
		int clusterMemberNum = scLib->rotCluster[targetNode->aaType][level1MinID][0];
		for(int i=0;i<clusterMemberNum;i++){
			targetNode->confTmp->updateScRotamer(rotList[scLib->rotCluster[targetNode->aaType][level1MinID][i+1]]);
			double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;
			for(int j=0;j<neighborNum;j++){
				e += ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
				e += ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);
			}

			if(e < minEne){
				minEne = e;
				minID = scLib->rotCluster[targetNode->aaType][level1MinID][i+1];
			}
		}
		targetNode->confTmp->updateScRotamer(rotList[minID]);
	}

	double e = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->conf->scConf->rot) * ec->pa->wtRot;
	double E = scLib->getEnergy(targetNode->conf->bbConf->rot->index1K, targetNode->confTmp->scConf->rot) * ec->pa->wtRot;
	cout << endl;
	printf("target node: %s", rn->intToTri(targetNode->aaType).c_str());
	printf("rotamer energy: %7.3f %7.3f\n",  e, E);
	for(int j=0;j<neighborNum;j++){
		double e1 = ec->getEnergyScSc(targetNode->conf, neighborNodes[j]->conf, neighborNodes[j]->seqID);
		double e2 = ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->conf, -neighborNodes[j]->seqID);

		double E1 = ec->getEnergyScSc(targetNode->confTmp, neighborNodes[j]->conf, neighborNodes[j]->seqID);
		double E2 = ec->getEnergyBbSc(neighborNodes[j]->conf, targetNode->confTmp, -neighborNodes[j]->seqID);

		printf("nat  neighbor %2d %s eScSc: %7.3f eBbSc: %7.3f\n", neighborNodes[j]->seqID, rn->intToTri(neighborNodes[j]->aaType).c_str(), e1, e2);
		printf("pred neighbor %2d %s eScSc: %7.3f eBbSc: %7.3f\n", neighborNodes[j]->seqID, rn->intToTri(neighborNodes[j]->aaType).c_str(), E1, E2);


	}
	cout << endl;
}

void SingleSidechainModelingTemplate::printInit(const string& output, AtomLib* atLib){
	ResName rn;

	ProteinChain* pc1 = new ProteinChain('A');

	Residue* resTarget1 = new Residue("0", 'A', rn.intToTri(targetNode->aaType));
	resTarget1->addAtom(new Atom("N", targetNode->conf->bbConf->coords[1]));
	resTarget1->addAtom(new Atom("CA", targetNode->conf->bbConf->coords[2]));
	resTarget1->addAtom(new Atom("C", targetNode->conf->bbConf->coords[3]));
	resTarget1->addAtom(new Atom("O", targetNode->conf->bbConf->coords[4]));
	vector<string>* scNames = atLib->getAminoAcidSidechainAtomNames(targetNode->aaType);

	for(int i=0;i<targetNode->conf->scConf->atomNum;i++){
		resTarget1->addAtom(new Atom(scNames->at(i), targetNode->conf->bbConf->coords[i]));
	}
	pc1->addResidue(resTarget1);



	vector<Residue*> nbResList;
	char xx[200];
	for(int i=0;i<neighborNodes.size();i++){
		sprintf(xx, "%d", (i+1));
		string resid = string(xx);
		Residue* res = new Residue(resid, 'A', rn.intToTri(neighborNodes[i]->aaType));
		res->addAtom(new Atom("N", neighborNodes[i]->conf->bbConf->coords[1]));
		res->addAtom(new Atom("CA", neighborNodes[i]->conf->bbConf->coords[2]));
		res->addAtom(new Atom("C", neighborNodes[i]->conf->bbConf->coords[3]));
		res->addAtom(new Atom("O", neighborNodes[i]->conf->bbConf->coords[4]));

		vector<string>* scList = atLib->getAminoAcidSidechainAtomNames(neighborNodes[i]->aaType);
		for(int j=0;j<neighborNodes[i]->conf->scConf->atomNum;j++){
			res->addAtom(new Atom(scList->at(j), neighborNodes[i]->conf->scConf->coords[j]));
		}
		nbResList.push_back(res);
		pc1->addResidue(res);
	}


	ofstream file;
	file.open(output.c_str(), ios::out);
	pc1->printPDBFormat(file, 1);
	file.close();

	Residue* p;

	p = resTarget1;
	for(int i=0;i<p->getAtomList()->size();i++){
		delete p->getAtomList()->at(i);
	}
	delete p;

	for(int k=0;k<nbResList.size();k++) {
		p = nbResList[k];
		for(int i=0;i<p->getAtomList()->size();i++){
			delete p->getAtomList()->at(i);
		}
		delete p;
	}

	delete pc1;
}

void SingleSidechainModelingTemplate::printPDB(const string& rawState, const string& finalState, AtomLib* atLib){
	ResName rn;

	ProteinChain* pc1 = new ProteinChain('A');
	ProteinChain* pc2 = new ProteinChain('A');

	Residue* resTarget1 = new Residue("0", 'A', rn.intToTri(targetNode->aaType));
	resTarget1->addAtom(new Atom("N", targetNode->conf->bbConf->coords[1]));
	resTarget1->addAtom(new Atom("CA", targetNode->conf->bbConf->coords[2]));
	resTarget1->addAtom(new Atom("C", targetNode->conf->bbConf->coords[3]));
	resTarget1->addAtom(new Atom("O", targetNode->conf->bbConf->coords[4]));
	vector<string>* scNames = atLib->getAminoAcidSidechainAtomNames(targetNode->aaType);

	for(int i=0;i<targetNode->conf->scConf->atomNum;i++){
		resTarget1->addAtom(new Atom(scNames->at(i), targetNode->conf->scConf->coords[i]));
	}
	pc1->addResidue(resTarget1);

	Residue* resTarget2 = new Residue("0", 'A', rn.intToTri(targetNode->aaType));
	resTarget2->addAtom(new Atom("N", targetNode->conf->bbConf->coords[1]));
	resTarget2->addAtom(new Atom("CA", targetNode->conf->bbConf->coords[2]));
	resTarget2->addAtom(new Atom("C", targetNode->conf->bbConf->coords[3]));
	resTarget2->addAtom(new Atom("O", targetNode->conf->bbConf->coords[4]));

	for(int i=0;i<targetNode->conf->scConf->atomNum;i++){
		resTarget2->addAtom(new Atom(scNames->at(i), targetNode->confTmp->scConf->coords[i]));
	}
	pc2->addResidue(resTarget2);


	vector<Residue*> nbResList;
	char xx[200];
	for(int i=0;i<neighborNodes.size();i++){
		sprintf(xx, "%d", (i+1));
		string resid = string(xx);
		Residue* res = new Residue(resid, 'A', rn.intToTri(neighborNodes[i]->aaType));
		res->addAtom(new Atom("N", neighborNodes[i]->conf->bbConf->coords[1]));
		res->addAtom(new Atom("CA", neighborNodes[i]->conf->bbConf->coords[2]));
		res->addAtom(new Atom("C", neighborNodes[i]->conf->bbConf->coords[3]));
		res->addAtom(new Atom("O", neighborNodes[i]->conf->bbConf->coords[4]));

		vector<string>* scList = atLib->getAminoAcidSidechainAtomNames(neighborNodes[i]->aaType);
		for(int j=0;j<neighborNodes[i]->conf->scConf->atomNum;j++){
			res->addAtom(new Atom(scList->at(j), neighborNodes[i]->conf->scConf->coords[j]));
		}
		nbResList.push_back(res);
		pc1->addResidue(res);
		pc2->addResidue(res);
	}


	ofstream file;
	file.open(rawState.c_str(), ios::out);
	pc1->printPDBFormat(file, 1);
	file.close();

	file.open(finalState.c_str(), ios::out);
	pc2->printPDBFormat(file, 1);
	file.close();

	Residue* p;

	p = resTarget1;
	for(int i=0;i<p->getAtomList()->size();i++){
		delete p->getAtomList()->at(i);
	}
	delete p;
	p = resTarget2;
	for(int i=0;i<p->getAtomList()->size();i++){
		delete p->getAtomList()->at(i);
	}
	delete p;
	for(int k=0;k<nbResList.size();k++) {
		p = nbResList[k];
		for(int i=0;i<p->getAtomList()->size();i++){
			delete p->getAtomList()->at(i);
		}
		delete p;
	}

	delete pc1;
	delete pc2;
}

SingleSidechainModelingTemplate::~SingleSidechainModelingTemplate() {

	delete targetNode;
	for(int i=0;i<neighborNodes.size();i++){
		delete neighborNodes[i];
	}
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
