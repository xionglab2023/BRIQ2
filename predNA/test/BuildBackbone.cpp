/*
 * BuidBackbone.cpp
 *
 */


#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "model/StructureModel.h"
#include "predNA/BackboneModelingTemplate.h"

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpredna;
using namespace std;

int main(int argc, char** argv){

	string pdbFile = string(argv[1]);
	//string paraFile = string(argv[2]);
	string paraTag  = string(argv[2]);
	string output = string(argv[3]);
	//string paraFile = string(argv[3]);




	RNAPDB* pdb = new RNAPDB(pdbFile);
	int len = pdb->getBaseList().size();

	if(paraTag == "init") {
		ForceFieldPara* para = new ForceFieldPara();
		BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
		BRTreeInfo* bi = bm->toTreeInfo();
		bi->printPDB(output+".pdb");
	}
	else if(paraTag == "sd") {
		ForceFieldPara* para = new ForceFieldPara();
		BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
		double rms = bm->runMC();
		BRTreeInfo* bi = bm->toTreeInfo();
		bi->printPDB(output+".pdb");
	}
	if(paraTag == "clashNb") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=1.0;p<5.0;p=p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->clashLamdaNb = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "clashNnb") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=1.0;p<5.0;p=p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->clashLamdaNnb = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "ang") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=0.01;p<0.4;p+=0.01){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaKAng = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "bond") {
		ofstream out;
		out.open(output, ios::out);	
		for(double p=1.0;p<20;p+=1.0){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaKBond = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "dihed") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.1;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaDihedImpD1D2Shift[0] = p;
			para->rnaDihedImpD1D2Shift[3] = p;
			para->rnaDihedImpD4D5Shift[0] = p;
			para->rnaDihedImpD4D5Shift[3] = p;
			para->rnaDihedD2D3D4Shift[0] = p;
			para->rnaDihedD2D3D4Shift[3] = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "dihed1") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.1;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaDihedImpD1D2Shift[0] = p;
			para->rnaDihedImpD1D2Shift[3] = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "dihed2") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.1;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaDihedImpD4D5Shift[0] = p;
			para->rnaDihedImpD4D5Shift[3] = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "dihed3") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.1;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaDihedD2D3D4Shift[0] = p;
			para->rnaDihedD2D3D4Shift[3] = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "rot") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=0.1;p<2.01;p+=0.1) {
			ForceFieldPara* para = new ForceFieldPara();
			//para->wtRibose = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "rot1") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.01;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaRiboseRotamerShift[0] += p;
			para->rnaRiboseRotamerShift[1] += p;

			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "rot2") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.01;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaRiboseRotamerShift[4] += p;
			para->rnaRiboseRotamerShift[5] += p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "rot3") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.01;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaRiboseRotamerShift[8] += p;
			para->rnaRiboseRotamerShift[9] += p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "rot4") {
		ofstream out;
		out.open(output, ios::out);
		for(double p=-2.0;p<2.01;p+=0.2){
			ForceFieldPara* para = new ForceFieldPara();
			para->rnaRiboseRotamerShift[12] += p;
			para->rnaRiboseRotamerShift[13] += p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}
	else if(paraTag == "pho"){
		ofstream out;
		out.open(output, ios::out);
		for(double p=0.1;p<1.51;p+=0.1) {
			ForceFieldPara* para = new ForceFieldPara();
			para->wtPho = p;
			BackboneModelingTemplate* bm = new BackboneModelingTemplate(pdbFile, para);
			double rms = bm->runMC();
			out << p << " " << rms << " " << len << endl;
			delete bm;
		}
		out.close();
	}




}

