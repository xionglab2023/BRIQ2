/*
 * RNAGraph.cpp
 *
 *  Created on: 2020Äê11ÔÂ24ÈÕ
 *      Author: pengx
 */

#include "motif/RNAGraph.h"

namespace NSPmotif {

SubRNAGraph::~SubRNAGraph(){
	delete [] this->connectToDownstream;
	delete [] this->fragIndex;
	delete [] this->inGraph;
	delete [] this->choosen;
}

HelixGraph::~HelixGraph(){
	//cout << "delete helix graph" << endl;
	delete [] interactionsOutsideHelix;
}


RNAGraph::RNAGraph(const string& ssFile, const string& eneFile) {
	// TODO Auto-generated constructor stub
	ifstream file;
	file.open(ssFile.c_str(), ios::in);
	if(!file.is_open()){
		cout << "can't open file "<< ssFile << endl;
		exit(0);
	}
	string line;
	getline(file, line);
	if(line.length() < 5){
		this->seqLen = 0;
		this->seq = "";
		this->ss = "";
		return;
	}

	if(line.substr(0,3) != "seq"){
		cout << "invalid file: " << endl;
		cout << line << endl;
		exit(0);
	}

	this->seq = line.substr(4, line.length()-4);

	getline(file, line);
	if(line.length() < 5){
		this->seqLen = 0;
		this->seq = "";
		this->ss = "";
		return;
	}
	if(line.substr(0,3) != "sec"){
		cout << "invalid file: " << endl;
		cout << line << endl;
		exit(0);
	}

	this->ss = line.substr(4, line.length()-4);
	string ssString = this->ss;

	vector<string> spt;
	vector<int> chainBreaks;
	getline(file,line);
	if(getline(file, line)){
		if(line.length() > 7 && line.substr(0,5) == "break"){
			splitString(line, " ", &spt);
			for(int i=1;i<spt.size();i++){
				chainBreaks.push_back(atoi(spt[i].c_str()));
			}
		}
	}

	file.close();

	this->seqLen = seq.length();
	if(seq.length() != ss.length()){
		cout << "sequence length not equal to ss length" << endl;
		exit(0);
	}
	this->allNodes = new RGNode*[this->seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->pairingIndex = new int[seqLen];
	this->subGraphIndex = new int[seqLen];
	this->contactEnergyOfNonNeighbor = new double[seqLen];
	this->allEnergy = new double[seqLen*seqLen];
	this->inHelix = new bool[seqLen];
	this->choosen = new bool[seqLen];
	this->initialGraphIndex = new int[seqLen];
	this->initialClusterNumber = 0;

	for(int i=0;i<seqLen;i++){
		this->connectToDownstream[i] = true;
		this->subGraphIndex[i] = 0;
		this->initialGraphIndex[i] = -1;
		this->contactEnergyOfNonNeighbor[i] = 0.0;
		this->inHelix[i] = false;
		this->choosen[i] = false;
	}

	int n1 = seqLen*seqLen;
	for(int i=0;i<n1;i++){
		this->allEnergy[i] = 0.0;
	}

	this->connectToDownstream[seqLen-1] = false;
	for(int i=0;i<chainBreaks.size();i++){
		this->connectToDownstream[chainBreaks[i]] = false;
	}

	for(int i=0;i<seqLen-1;i++){
		if(connectToDownstream[i]){
			this->allEnergy[i*seqLen+i+1] = -0.0;
			this->allEnergy[(i+1)*seqLen+i] = -0.0;
		}
	}

	for(int i=0;i<seqLen;i++){
		pairingIndex[i] = -1;
	}
	for(int i=0;i<seqLen;i++) {
		char c = ss[i];
		if(c == ')') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '(') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cerr << "invalid ssSeq: " << ss << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			pairingIndex[i] = preIndex;
			pairingIndex[preIndex] = i;
		}
		if(c == ']') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '[') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cerr << "invalid ssSeq: " << ss << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			pairingIndex[i] = preIndex;
			pairingIndex[preIndex] = i;
		}
		if(c == '}') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '{') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cerr << "invalid ssSeq: " << ss << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			pairingIndex[i] = preIndex;
			pairingIndex[preIndex] = i;
		}
		if(c == '>') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '<') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cerr << "invalid ssSeq: " << ss << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			pairingIndex[i] = preIndex;
			pairingIndex[preIndex] = i;
		}
	}

	for(int i=0;i<seqLen;i++){
		if(ss[i] != '.')
		{
			cerr << "invalid ssSeq: " << ss << endl;
			exit(1);
		}
	}


	this->ss = ssString;

	for(int i=0;i<seqLen;i++){
		this->allNodes[i] = new RGNode(i, seq[i]);
	}

	vector<string> rotLines;
	vector<string> phoLines;
	vector<string> pairLines;
	file.open(eneFile.c_str(), ios::in);
	if(!file.is_open()){
		cout << "can't open file: " << eneFile << endl;
		exit(0);
	}
	while(getline(file, line)){
		if(line.length() < 5) continue;
		if(line.substr(0, 3) == "rot")
			rotLines.push_back(line);
		else if(line.substr(0,3) == "pho")
			phoLines.push_back(line);
		else if(line.substr(0,4) == "pair")
			pairLines.push_back(line);
	}

	file.close();
	if(rotLines.size() != seqLen){
		cout << "eneFile and ssFile not consistent" << endl;
		exit(0);
	}

	int idA, idB;
	double bb, bo, bp, e;
	for(int i=0;i<pairLines.size();i++){
		splitString(pairLines[i], " ", &spt);
		if(spt.size() < 5) continue;
		idA = atoi(spt[1].c_str());
		idB = atoi(spt[2].c_str());
		bb = atof(spt[3].c_str());
		bo = atof(spt[4].c_str());
		bp = atof(spt[5].c_str());
		e = bb + bo + bp;



		if(idA < 0 || idA >= seqLen){
			cout << "invalid line: " << pairLines[i] << endl;
			exit(0);
		}
		if(idB < 0 || idB >= seqLen){
			cout << "invalid line: " << pairLines[i] << endl;
			exit(0);
		}

		/*
		if(idB == idA+1 && connectToDownstream[idA] && pairingIndex[idA] >=0 && pairingIndex[idB] >=0 && pairingIndex[idA] == pairingIndex[idB]+1 && connectToDownstream[pairingIndex[idB]]){
			//helix stacking pair
			e = e*0.5;
		}
		*/

		allEnergy[idA*seqLen+idB] += e;
		allEnergy[idB*seqLen+idA] += e;
		if(bb + bo + bp< 0){
			int sep = idB - idA;
			if(idB == idA+1 && !connectToDownstream[idA])
				sep = 5;
			if(idB == idA+2 && (!connectToDownstream[idA] || !connectToDownstream[idA+1]))
				sep = 5;

			if(sep > 1){
				contactEnergyOfNonNeighbor[idA] += 0.5*e;
				contactEnergyOfNonNeighbor[idB] += 0.5*e;
			}

			if(e< -0.01) {
				RGEdge* eg = new RGEdge(idA, idB, bb, bo,bp, sep);
				if(idB == idA+1 && connectToDownstream[idA]) eg->isSequentialNeighbor = true;
				if(pairingIndex[idA] == idB) eg->isWatsonCrickPair = true;
				this->egList.push_back(eg);
			}
		}
	}
}


void RNAGraph::splitHelix(int i){

	if(i < 0 || i >= seqLen-1) return;
	if(pairingIndex[i] < 0 || pairingIndex[i+1] < 0) return;
	if(!connectToDownstream[i] || !connectToDownstream[pairingIndex[i+1]]) return;

	removeEdge(i, i+1);
	removeEdge(pairingIndex[i+1], pairingIndex[i]);
	removeEdge(i, pairingIndex[i+1]);
	removeEdge(i+1, pairingIndex[i]);

}

void RNAGraph::removeEdgeOnNode(int i){
	set<int>::iterator it;
	for(it = allNodes[i]->egSet.begin();it != allNodes[i]->egSet.end();++it){
		int egIndex = *it;
		RGEdge* eg = egList[egIndex];
		allNodes[eg->indexA]->egSet.erase(egIndex);
		allNodes[eg->indexB]->egSet.erase(egIndex);
	}
}

void RNAGraph::addEdge(){
	for(int i=0;i<egList.size();i++){
		RGEdge* eg = egList[i];
		allNodes[eg->indexA]->egSet.insert(i);
		allNodes[eg->indexB]->egSet.insert(i);
		baseIndexToEgIndex[eg->indexA*seqLen+eg->indexB] = i;
		baseIndexToEgIndex[eg->indexB*seqLen+eg->indexA] = i;
	}

	map<int,int>::iterator it;
	for(int i=0;i<seqLen-1;i++){
		if(!connectToDownstream[i]) continue;
		it = this->baseIndexToEgIndex.find(i*seqLen+i+1);
		if(it == baseIndexToEgIndex.end()){
			RGEdge* eg = new RGEdge(i);
			int egIndex = egList.size();
			baseIndexToEgIndex[i*seqLen+i+1] = egIndex;
			baseIndexToEgIndex[(i+1)*seqLen+i] = egIndex;

			egList.push_back(eg);
			allNodes[i]->egSet.insert(egIndex);
			allNodes[i+1]->egSet.insert(egIndex);
		}
	}
}



void RNAGraph::removeEdge(int i, int j){
	map<int,int>::iterator it;
	it = baseIndexToEgIndex.find(i*seqLen+j);
	if(it == baseIndexToEgIndex.end()){
		//cout << "there is no edge from "<<i << " to " << j << endl;
		return;
	}

	allNodes[i]->egSet.erase(it->second);
	allNodes[j]->egSet.erase(it->second);
}

void RNAGraph::findHelix(){
	int helixStart = 0;
	int lastPairingID = -1;
	int helixLen = 0;
	bool inHelix = false;
	for(int i=0;i<this->seqLen;i++){
		if(!inHelix && pairingIndex[i] > i){
			//helix start
			inHelix = true;
			helixStart = i;
			lastPairingID = pairingIndex[i];
			helixLen = 1;
		}
		else if(inHelix && pairingIndex[i] > i && pairingIndex[i] == lastPairingID-1){
			//helix extend
			lastPairingID = pairingIndex[i];
			helixLen++;
		}
		else if(inHelix && pairingIndex[i] > i && pairingIndex[i] != lastPairingID-1){
			//helix turn, start new helix
			if(helixLen > 1){
				HelixGraph* hg = new HelixGraph(helixStart, helixLen);
				this->helixList.push_back(hg);
			}
			helixLen = 1;
			helixStart = i;
			lastPairingID = pairingIndex[i];
		}
		else if(inHelix && pairingIndex[i] < i ){
			//helix end
			if(helixLen > 1){
				HelixGraph* hg = new HelixGraph(helixStart, helixLen);
				this->helixList.push_back(hg);
			}
			helixLen = 0;
			inHelix = false;
		}
	}

	for(int i=0;i<helixList.size();i++){
		HelixGraph* hg = helixList[i];
		if(hg->length < 3) continue;
		int leftBegin = hg->beginIndex;
		int leftEnd = leftBegin+hg->length-1;
		int rightBegin = this->pairingIndex[leftEnd];
		int rightEnd = this->pairingIndex[leftBegin];

		for(int i=leftBegin;i<=leftEnd;i++){
			this->inHelix[i] = true;
		}
		for(int i=rightBegin;i<=rightEnd;i++){
			this->inHelix[i] = true;
		}

	}

	set<int>::iterator it;
	RGEdge* eg;
	int idA, idB;
	int egIndex;
	for(int i=0;i<helixList.size();i++){

		HelixGraph* hg = helixList[i];
		int leftBegin = hg->beginIndex;
		int leftEnd = leftBegin+hg->length-1;
		int rightBegin = this->pairingIndex[leftEnd];
		int rightEnd = this->pairingIndex[leftBegin];

		//cout << "helix: "<<i << " " << leftBegin << " " << leftEnd << " " << rightBegin << " " << rightEnd << endl;
		bool contact[this->seqLen];
		for(int k=0;k<seqLen;k++)
		{
			contact[k] = false;
		}

		for(int k=0;k<hg->length;k++){
			hg->interactionsOutsideHelix[k] = 0.0;
		}

		for(int j=leftBegin;j<=leftEnd;j++){
			double e = 0;
			int pairIndexJ = this->pairingIndex[j];

			for(it=this->allNodes[j]->egSet.begin();it!=this->allNodes[j]->egSet.end();++it){
				egIndex = *it;
				eg = egList[egIndex];
				idA = eg->indexA;
				idB = eg->indexB;

				if(idA == j && (getSep(idB, j) < 3 || getSep(idB, pairIndexJ) < 3 || (idB >= leftBegin && idB <= leftEnd) || (idB >= rightBegin && idB <= rightEnd))) continue;
				if(idB == j && (getSep(idA, j) < 3 || getSep(idA, pairIndexJ) < 3 || (idA >= leftBegin && idA <= leftEnd) || (idA >= rightBegin && idA <= rightEnd))) continue;


				e += eg->eBB + eg->eBO;
				if(eg->eBB + eg->eBO < -0.5){
					if(idA == j)
						contact[idB] = true;
					else
						contact[idA] = true;
				}
			}
			hg->interactionsOutsideHelix[j-leftBegin] += e;

			e = 0;
			for(it=this->allNodes[pairingIndex[j]]->egSet.begin();it!=this->allNodes[pairingIndex[j]]->egSet.end();++it){
				egIndex = *it;
				eg = egList[egIndex];
				idA = eg->indexA;
				idB = eg->indexB;

				if(idA == pairIndexJ && (getSep(idB, pairIndexJ) < 3 || getSep(idB, j) < 3 || (idB >= leftBegin && idB <= leftEnd) || (idB >= rightBegin && idB <= rightEnd))) continue;
				if(idB == pairIndexJ && (getSep(idA, pairIndexJ) < 3 || getSep(idA, j) < 3 || (idA >= leftBegin && idA <= leftEnd) || (idA >= rightBegin && idA <= rightEnd))) continue;

				e += eg->eBB + eg->eBO;
				if(eg->eBB + eg->eBO < -0.5){
					if(idA == pairingIndex[j])
						contact[idB] = true;
					else
						contact[idA] = true;
				}
			}
			hg->interactionsOutsideHelix[j-leftBegin] += e;

		}
		for(int i=0;i<seqLen;i++){
			if(contact[i]){
				hg->contactBaseList.push_back(i);
			}
		}
	}
}

void RNAGraph::initialAssign(int nodeID, int graphID){
	if(initialGraphIndex[nodeID] >= 0)
		return;
	this->initialGraphIndex[nodeID] = graphID;
	set<int>::iterator it;
	int idA, idB;
	for(it = allNodes[nodeID]->egSet.begin();it!= allNodes[nodeID]->egSet.end();++it){
		int egIndex = *it;
		if(egList[egIndex]->getEnergy() < -0.2)
		{
			idA = egList[egIndex]->indexA;
			idB = egList[egIndex]->indexB;
			if(idA == nodeID){
				initialAssign(idB, graphID);
			}
			else{
				initialAssign(idA, graphID);
			}
		}

	}
}

void RNAGraph::assign(int nodeID, int graphID, double T){
	if(subGraphIndex[nodeID] >= 0 || T <= 0)
		return;
	this->subGraphIndex[nodeID] = graphID;
	set<int>::iterator it;
	int idA, idB;
	for(it = allNodes[nodeID]->egSet.begin();it!= allNodes[nodeID]->egSet.end();++it){
		int egIndex = *it;
		double p = 1/(1+exp((egList[egIndex]->getEnergy()+4*T)/T));
		if((double)rand()/RAND_MAX < p){
			idA = egList[egIndex]->indexA;
			idB = egList[egIndex]->indexB;
			if(idA == nodeID){
				assign(idB, graphID, T);
			}
			else{
				assign(idA, graphID, T);
			}
		}

	}
}

void RNAGraph::updateSubGraph(double T, int initialClusterID){

	for(int i=0;i<this->subGraphs.size();i++){
		delete this->subGraphs[i];
	}
	subGraphs.clear();
	for(int i=0;i<seqLen;i++){
		this->subGraphIndex[i] = -1;
	}

	int currentIndex = 0;
	for(int i=0;i<seqLen;i++){
		if(choosen[i]) continue;
		if(initialGraphIndex[i] != initialClusterID) continue;
		if(subGraphIndex[i] < 0){
			assign(i, currentIndex, T);
			currentIndex++;
		}
	}

	int totalSubgraphNum = currentIndex;

	double scoreCutoff = 1.0;

	for(int i=0;i<totalSubgraphNum;i++){
		vector<int> idList;
		for(int j=0;j<seqLen;j++){
			if(subGraphIndex[j] == i){
				idList.push_back(j);
			}
		}
		if(idList.size() > 3){
			SubRNAGraph* subRG1 = new SubRNAGraph(idList, seqLen, this->connectToDownstream, this->choosen);

			vector<SubRNAGraph*> varList = subRG1->generateAllVariant();

			double minScore = getGraphScore(subRG1);
			int selectIndex = -1;

			for(int k=0;k<varList.size();k++){
				double score = getGraphScore(varList[k]);
				if(score < minScore){
					minScore = score;
					selectIndex = k;
				}
			}

			if(minScore > scoreCutoff){
				for(int k=0;k<varList.size();k++){
					delete varList[k];
				}
				delete subRG1;
			}
			else {
				if(selectIndex < 0){
					for(int k=0;k<varList.size();k++){
						delete varList[k];
					}
					this->subGraphs.push_back(subRG1);
				}
				else {
					delete subRG1;
					for(int k=0;k<varList.size();k++){
						if(k != selectIndex)
							delete varList[k];
					}
					this->subGraphs.push_back(varList[selectIndex]);
				}
			}
		}
	}

}

double RNAGraph::getEdgeEnergy(int i, int j){
	return this->allEnergy[i*seqLen+j];
}

double RNAGraph::getGraphScore(SubRNAGraph* sg){
	double totalE = 0;
	double wcE = 0;
	double nbE = 0;
	double otherE = 0;

	for(int i=0;i<seqLen;i++){
		if(!sg->inGraph[i]) continue;
		for(int j=i+1;j<seqLen;j++){
			if(!sg->inGraph[j]) continue;
			if(pairingIndex[i] == j && inHelix[i])
				wcE += this->allEnergy[i*seqLen + j];
			else if(j==i+1 && connectToDownstream[i])
				nbE += this->allEnergy[i*seqLen + j];
			else
				otherE += this->allEnergy[i*seqLen + j];
		}
	}

	totalE += wcE*0.5 + nbE*0.7 + otherE;

	int B = 0;
	for(int i=0;i<seqLen;i++){
		if(sg->inGraph[i]) B++;
	}

	double F = 0.0;
	if(sg->inGraph[0])
		F = 1.0;

	for(int i=0;i<seqLen-1;i++){
		if(!sg->inGraph[i] && sg->inGraph[i+1])
			F += 1.0;
	}

	int helixBaseNum = 0;
	for(int i=0;i<seqLen;i++){
		if(sg->inGraph[i] && pairingIndex[i] >=0 && sg->inGraph[pairingIndex[i]]){
			helixBaseNum++;
		}
	}

	if(helixBaseNum == B)
		return 99.9;

	if(B < 4)
		return 9.9;

	double ff = 2.0*F*F;

	/*
	if(F > 4){
		ff = 10.0*F - 20;
	}
	*/

	if(B > 25)
		ff += 2*(B - 25);

	return (totalE+ff)/B;
}

void RNAGraph::printConnect(){
	set<int>::iterator it;
	for(int i=0;i<seqLen;i++){
		for(it=allNodes[i]->egSet.begin();it!= allNodes[i]->egSet.end();++it){
			int egIndex = *it;
			RGEdge* eg = egList[egIndex];
			if(eg->indexA == i){
				double e = eg->eBB + eg->eBO;
				if(e < -0.5){
					printf("edge: %3d %3d %8.3f\n", i, eg->indexB, e);
				}
			}
		}
	}
}

void RNAGraph::processLoop(){
	int loopStart = 0;
	int loopLength = 0;
	bool inLoop = false;
	double eneCutoff = -0.2;
	int loopLengthCutoff = 3;
	for(int i=0;i<seqLen;i++){

		if(!inLoop && contactEnergyOfNonNeighbor[i] > eneCutoff){ //loop start
			loopLength = 1;
			loopStart = i;
			inLoop = true;
		}
		else if(inLoop && contactEnergyOfNonNeighbor[i] > eneCutoff){
			// loop continue
			loopLength++;
		}
		else if(inLoop && contactEnergyOfNonNeighbor[i] <= eneCutoff) {
			// loop end
			if(loopLength >= loopLengthCutoff){
				for(int k=loopStart;k<loopStart+loopLength;k++){
					removeEdgeOnNode(k);
					choosen[k] = true;
				}
			}

			inLoop = false;
			loopLength = 0;
		}
	}

	if(inLoop && loopLength >= loopLengthCutoff){
		for(int k=loopStart;k<loopStart+loopLength;k++){
			removeEdgeOnNode(k);
			choosen[k] = true;
		}
	}
}

void RNAGraph::processHelix(){

	for(int i=0;i<this->helixList.size();i++){
		HelixGraph* hg = helixList[i];

		for(int j=hg->beginIndex;j<hg->beginIndex+hg->length-1;j++){
			if(hg->interactionsOutsideHelix[j] > -0.1 || hg->interactionsOutsideHelix[j+1] > -0.1){
				//cout << "split at" << j << endl;
				splitHelix(j);
			}
		}

		/*
		for(int j=hg->beginIndex;j<hg->beginIndex+hg->length-1;j++){
			if(hg->interactionsOutsideHelix[j] > -0.3 && hg->interactionsOutsideHelix[j+1] > -0.3){
				splitHelix(j);
			}
		}
		*/

		if(hg->isIndependantHelix()){
			int leftBegin = hg->simpleHelixBegin;
			int leftEnd = hg->simpleHelixBegin+ hg->simpleHelixLen-1;
			int rightBegin = pairingIndex[leftEnd];
			int rightEnd = pairingIndex[leftBegin];

			for(int j=leftBegin;j<=leftEnd;j++){
				//cout << "remove node: " << j << endl;
				removeEdgeOnNode(j);
				this->choosen[j] = true;
			}

			for(int j=rightBegin;j<=rightEnd;j++){
				//cout << "remove node: " << j << endl;
				removeEdgeOnNode(j);
				this->choosen[j] = true;
			}
		}
	}
}

void RNAGraph::generateInitialClusters(){

	for(int i=0;i<seqLen;i++){
		this->initialGraphIndex[i] = -1;
		if(choosen[i])
			this->initialGraphIndex[i] = 9999999;
	}

	int currentIndex = 0;
	for(int i=0;i<seqLen;i++){
		if(choosen[i]) continue;
		if(initialGraphIndex[i] >= 0) continue;
		initialAssign(i, currentIndex);
		currentIndex++;
	}

	this->initialClusterNumber = currentIndex;
}

void RNAGraph::printInitialCluster(int id){
	for(int i=0;i<seqLen;i++){
		if(initialGraphIndex[i] == id){
			cout << seq[i];
		}
	}
	cout << endl;
	for(int i=0;i<seqLen;i++){
		if(initialGraphIndex[i] == id){
			cout << ss[i];
		}
	}
	cout << endl;
}

SubRNAGraph* RNAGraph::findBestSubGraph(int initialClusterID){
	double minE = 999.9;
	SubRNAGraph* sg = NULL;

	double bestT = 0;

	for(double T=0.01; T < 1.0; T += 0.002){
		//cout << "T " << T << endl;
		updateSubGraph(T, initialClusterID);
		for(int i=0;i<this->subGraphs.size();i++){
			double e = getGraphScore(subGraphs[i]);
			if(e < minE){
				minE = e;
				if(sg!= NULL) delete sg;
				sg = subGraphs[i]->copy();
				bestT = T;
			}
		}
	}
	//cout << "bestT: " << bestT  << " " << initialClusterID << endl;
	return sg;
}


vector<SubRNAGraph*> RNAGraph::findAllGraphs(){

	vector<SubRNAGraph*> sgList;

	for(int i=0;i<this->initialClusterNumber;i++){
		//printInitialCluster(i);
		SubRNAGraph* sg = findBestSubGraph(i);
		while(sg != NULL){
			sgList.push_back(sg);
			deleteSubGraph(sg);
			sg = findBestSubGraph(i);
		}
	}
	return sgList;

}

void RNAGraph::deleteSubGraph(SubRNAGraph* sg){
	for(int i=0;i<seqLen;i++){
		if(sg->inGraph[i]){
			choosen[i] = true;
			removeEdgeOnNode(i);
		}
	}
}

RNAGraph::~RNAGraph() {
	//cout << "delete rna graph" << endl;
	for(int i=0;i<this->seqLen;i++){
		delete allNodes[i];
	}
	delete [] this->allNodes;
	delete [] connectToDownstream;
	delete [] this->pairingIndex;
	delete [] this->subGraphIndex;
	delete [] this->initialGraphIndex;
	delete [] allEnergy;
	delete [] inHelix;
	delete [] contactEnergyOfNonNeighbor;
	delete [] choosen;

	for(int i=0;i<this->egList.size();i++){
		delete this->egList[i];
	}
	for(int i=0;i<this->helixList.size();i++){
		delete this->helixList[i];
	}
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmotif */
