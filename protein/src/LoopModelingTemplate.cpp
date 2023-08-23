/*
 * LoopModelingTemplate.cpp
 *
 */

#include "protein/LoopModelingTemplate.h"

namespace NSPprotein {

LoopModelingTemplate::LoopModelingTemplate(int len, vector<string>& lines, ResBBRotamerLib* bbLib, ResScRotamerLib* scLib, ResName* rn) {
	// TODO Auto-generated constructor stub
	this->rn = rn;
	this->rn = rn;
	this->bbLib = bbLib;
	this->scLib = scLib;
	this->loopLength = len;

	cout << "loop len: " << len << endl;
	for(int i=0;i<loopLength;i++){
		cout << "node: " << i << endl;
		this->targetNodes.push_back(new ResNode(lines[i], i, rn, bbLib, scLib));
	}


	vector<string> spt;
	for(int i=loopLength;i<lines.size();i++) {
		cout << "nb: " << i << endl;
		splitString(lines[i], " ", &spt);
		int sep = atoi(spt[1].c_str());
		this->seqIDList.push_back(sep);
		this->neighborNodes.push_back(new ResNode(lines[i],i, rn, bbLib, scLib));
	}
}

void LoopModelingTemplate::printInit(const string& output, AtomLib* atLib){
	ResName rn;

	ProteinChain* pc1 = new ProteinChain('A');
	vector<Residue*> targetResList;
	char xx[200];

	for(int k=0;k<loopLength;k++){
		sprintf(xx, "%d", (k+1));
		Residue* resTarget1 = new Residue(string(xx), 'A', rn.intToTri(targetNodes[k]->aaType));

		resTarget1->addAtom(new Atom("N", targetNodes[k]->conf->bbConf->coords[1]));
		resTarget1->addAtom(new Atom("CA", targetNodes[k]->conf->bbConf->coords[2]));
		resTarget1->addAtom(new Atom("C", targetNodes[k]->conf->bbConf->coords[3]));
		resTarget1->addAtom(new Atom("O", targetNodes[k]->conf->bbConf->coords[4]));
		vector<string>* scNames = atLib->getAminoAcidSidechainAtomNames(targetNodes[k]->aaType);

		for(int i=0;i<targetNodes[k]->conf->scConf->atomNum;i++){
			resTarget1->addAtom(new Atom(scNames->at(i), targetNodes[k]->conf->scConf->coords[i]));
		}
		pc1->addResidue(resTarget1);
		targetResList.push_back(resTarget1);
	}

	ProteinChain* pc2 = new ProteinChain('B');
	vector<Residue*> nbResList;

	for(int i=0;i<neighborNodes.size();i++){
		sprintf(xx, "%d", (i+1));
		string resid = string(xx);
		Residue* res = new Residue(resid, 'B', rn.intToTri(neighborNodes[i]->aaType));

		res->addAtom(new Atom("N", neighborNodes[i]->conf->bbConf->coords[1]));
		res->addAtom(new Atom("CA", neighborNodes[i]->conf->bbConf->coords[2]));
		res->addAtom(new Atom("C", neighborNodes[i]->conf->bbConf->coords[3]));
		res->addAtom(new Atom("O", neighborNodes[i]->conf->bbConf->coords[4]));

		vector<string>* scList = atLib->getAminoAcidSidechainAtomNames(neighborNodes[i]->aaType);
		for(int j=0;j<neighborNodes[i]->conf->scConf->atomNum;j++){
			res->addAtom(new Atom(scList->at(j), neighborNodes[i]->conf->scConf->coords[j]));
		}
		nbResList.push_back(res);
		pc2->addResidue(res);
	}


	ofstream file;
	file.open(output.c_str(), ios::out);
	pc1->printPDBFormat(file, 1);
	pc2->printPDBFormat(file, 1000);
	file.close();

	Residue* p;

	for(int k=0;k<loopLength;k++){
		p = targetResList[k];
		for(int i=0;i<p->getAtomList()->size();i++){
			delete p->getAtomList()->at(i);
		}
		delete p;
	}

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

LoopModelingTemplate::~LoopModelingTemplate() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
