/*
 * NuGraph.cpp
 *
 *  Created on: 2023Äê11ÔÂ15ÈÕ
 *      Author: nuc
 */

#include <predNA/NuGraph.h>

namespace NSPpredNA {

NuNode::NuNode(int id, int baseType, LocalFrame& cs1, RiboseRotamer* riboRot, AtomLib* atLib){


	this->seqID = id;
	this->baseType = baseType;
	this->connectToNeighbor = false;

	//cout << "init base rot" << endl;
	BaseRotamer* baseRot = new BaseRotamer(baseType, atLib);
	BaseRotamer* baseRotTmp = new BaseRotamer(baseType, atLib);

	//cout << "init base conf" << endl;
	this->baseConf = new BaseConformer(baseRot, cs1);
	this->baseConfTmp = new BaseConformer(baseRotTmp, cs1);

	//cout << "init ribose conf" << endl;
	this->riboseConf = new RiboseConformer(riboRot, cs1);
	this->riboseConfTmp = new RiboseConformer(riboRot, cs1);

	//cout << "init pho conf" << endl;
	this->phoConf = new PhosphateConformer();
	this->phoConfTmp = new PhosphateConformer();

	//cout << "init base rot CG" << endl;
	BaseRotamerCG* baseRotCG = new BaseRotamerCG(baseType, atLib);
	BaseRotamerCG* baseRotCGTmp = new BaseRotamerCG(baseType, atLib);

	//cout << "init base conf CG" << endl;
	this->baseConfCG = new BaseConformerCG(baseRotCG, cs1);
	this->baseConfCGTmp = new BaseConformerCG(baseRotCGTmp, cs1);

	//cout << "init ribose rot CG" << endl;
	RiboseRotamerCG* riboRotCG = new RiboseRotamerCG(riboRot);
	RiboseRotamerCG* riboRotCGTmp = new RiboseRotamerCG(riboRot);

	//cout << "init ribose conf CG" << endl;
	this->riboseConfCG = new RiboseConformerCG(riboRotCG, cs1);
	this->riboseConfCGTmp = new RiboseConformerCG(riboRotCGTmp, cs1);

	this->graph = NULL;
	this->ene = 0;
	this->eneTmp = 0;
	this->eneCG = 0;
	this->eneCGTmp = 0;
}

NuNode::~NuNode(){

	delete this->baseConf->rot;
	delete this->baseConfTmp->rot;

	delete this->baseConf;
	delete this->baseConfTmp;
	delete this->riboseConf;
	delete this->riboseConfTmp;
	delete this->phoConf;
	delete this->phoConfTmp;

	delete this->baseConfCG->rot;
	delete this->baseConfCGTmp->rot;
	delete this->riboseConfCG->rot;
	delete this->riboseConfCGTmp->rot;
	delete this->baseConfCG;
	delete this->baseConfCGTmp;
	delete this->riboseConfCG;
	delete this->riboseConfCGTmp;
}

void NuNode::updateNodeInformation(NuTree* tree){

	this->graph = tree->graph;

	int i,j;
	this->neighborList.clear();
	this->connectionBreakPoints.clear();


	for(i=0;i<graph->len;i++){
		if(i == seqID) continue;
		if(tree->adjMtx[seqID*tree->graph->len+i])
			this->neighborList.push_back(i);
	}

	if(i>0 && tree->graph->connectToDownstream[i-1])
		this->connectionBreakPoints.push_back(i-1);

	if(tree->graph->connectToDownstream[i])
		this->connectionBreakPoints.push_back(i);

	/*
	 * update energy
	 */
}

void NuNode::updateRiboseRotamer(RiboseRotamer* rot){
	this->riboseConfTmp->updateRotamer(rot);
	/*
	 * update energy tmp
	 */
}

void NuNode::acceptRotMutation(){
	this->riboseConf->copyValueFrom(riboseConfTmp);
	this->phoConf->copyValueFrom(phoConfTmp);
	this->ene = this->eneTmp;
}

void NuNode::clearRotMutation() {
	this->riboseConfTmp->copyValueFrom(riboseConf);
	this->phoConfTmp->copyValueFrom(phoConf);
	this->eneTmp = this->ene;
}

void NuNode::updateCoordinate(LocalFrame& cs){
	this->baseConfTmp->updateCoords(cs);
	this->riboseConfTmp->updateLocalFrame(cs);
	this->phoConfTmp->updateLocalFrame(riboseConfTmp->cs2);
}

void NuNode::acceptCoordMove(){
	this->baseConf->copyValueFrom(baseConfTmp);
	this->riboseConf->copyValueFrom(riboseConfTmp);
	this->phoConf->copyValueFrom(phoConfTmp);
}

void NuNode::clearCoordMove() {
	this->baseConfTmp->copyValueFrom(baseConf);
	this->riboseConfTmp->copyValueFrom(riboseConf);
	this->phoConfTmp->copyValueFrom(phoConf);
}

double NuNode::mutEnergy(){
	return this->eneTmp - this->ene;
}

void NuNode::updateRiboseRotamerCG(RiboseRotamerCG* rot){
	this->riboseConfCGTmp->updateRotamer(rot);
	/*
	 * update energy CG tmp
	 */
}

void NuNode::acceptRotMutationCG(){
	this->riboseConfCG->copyValueFrom(this->riboseConfCGTmp);
	this->eneCG = this->eneCGTmp;
}

void NuNode::clearRotMutationCG(){
	this->riboseConfCGTmp->copyValueFrom(this->riboseConfCG);
	this->eneCGTmp = this->eneCG;
}

void NuNode::updateCoordinateCG(LocalFrame& cs){
	this->baseConfCGTmp->updateCoords(cs);
	this->riboseConfCGTmp->updateLocalFrame(cs);
}

void NuNode::acceptCoordMoveCG(){
	this->baseConfCG->copyValueFrom(baseConfCGTmp);
	this->riboseConfCG->copyValueFrom(riboseConfCGTmp);
}

void NuNode::clearCoordMoveCG(){
	this->baseConfCGTmp->copyValueFrom(baseConfCG);
	this->riboseConfCGTmp->copyValueFrom(riboseConfCG);
}

double NuNode::mutEnergyCG(){
	return this->eneCGTmp - this->eneCG;
}



NuEdge::NuEdge(NuNode* nodeA, NuNode* nodeB, NuGraph* graph){
	this->graph = graph;
	this->indexA = nodeA->seqID;
	this->indexB = nodeB->seqID;
	this->nodeA = nodeA;
	this->nodeB = nodeB;
	this->cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->cmTmp = this->cm;
	this->sep = graph->sepTable[indexA*graph->len+indexB];

	int typeA = nodeA->baseType%4;
	int typeB = nodeB->baseType%4;

	this->ei = new EdgeInformation(sep, typeA, typeB, graph->pairLib);
	this->moveSet = new MixedNuPairCluster(sep, typeA*4+typeB, graph->moveLib);
	this->moveSet->updateEdgeInformation(this->ei);
	this->weight = this->ei->weight;

	this->weightRand = weight;

	this->pairLib = graph->pairLib;
}

NuEdge::NuEdge(NuNode* nodeA, NuNode* nodeB, int sep, BasePairLib* pairLib, NuPairMoveSetLibrary* moveLib){
	this->graph = NULL;
	this->indexA = nodeA->seqID;
	this->indexB = nodeB->seqID;
	this->nodeA = nodeA;
	this->nodeB = nodeB;
	this->cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->cmTmp = this->cm;
	this->sep = sep;
	int typeA = nodeA->baseType%4;
	int typeB = nodeB->baseType%4;

	this->ei = new EdgeInformation(sep, typeA, typeB, pairLib);
	this->moveSet = new MixedNuPairCluster(sep, typeA*4+typeB, moveLib);
	this->moveSet->updateEdgeInformation(this->ei);
	this->weight = this->ei->weight;
	this->weightRand = weight;

	this->pairLib = pairLib;
}


NuEdge::~NuEdge(){
	delete this->ei;
	delete this->moveSet;
}

void NuEdge::initNativeMoveSet(){

	BaseDistanceMatrix dm(nodeA->baseConf->cs1, nodeB->baseConf->cs1);
	int clusterID = pairLib->getPairType(dm, nodeA->baseType, nodeB->baseType, sep);
	this->ei->setUniqueCluster(clusterID, pairLib);
	this->moveSet->updateEdgeInformation(ei);
	this->weight = this->ei->weight;
	this->weightRand = this->weight;
}

void NuEdge::fixNaiveMove(){
	CsMove cm = nodeB->baseConf->cs1 - nodeA->baseConf->cs1;
	this->ei->setFixed(cm);
	this->moveSet->fixNativeMove(cm);
	this->weight = this->ei->weight;
	this->weightRand = this->weight;
}

void NuEdge::updateSubTrees(NuTree* tree){
	//BFS algorithm

	int i,j,k;
	queue<NuNode*> qA;
	queue<NuNode*> qB;
	NuNode* vn;
	NuNode* vw;

	bool visited[graph->len];
	for(i=0;i<graph->len;i++){
		visited[i] = false;
	}

	qA.push(nodeA);
	qB.push(nodeB);

	nodeListA.push_back(nodeA);
	nodeListB.push_back(nodeB);

	visited[nodeA->seqID] = true;
	visited[nodeB->seqID] = true;

	while(!qA.empty()){
		vn = qA.front();
		qA.pop();
		for(i=0;i<vn->neighborList.size();i++){
			if(visited[vn->neighborList[i]]) continue;
			vw = graph->allNodes[vn->neighborList[i]];
			qA.push(vw);
			nodeListA.push_back(vw);
			edgeListA.push_back(graph->allEdges[vn->seqID*graph->len+vw->seqID]);
			visited[vw->seqID] = true;
		}
	}

	while(!qB.empty()){
		vn = qB.front();
		qB.pop();
		for(i=0;i<vn->neighborList.size();i++){
			if(visited[vn->neighborList[i]]) continue;
			vw = graph->allNodes[vn->neighborList[i]];
			qB.push(vw);
			nodeListB.push_back(vw);
			edgeListB.push_back(graph->allEdges[vn->seqID*graph->len+vw->seqID]);
			visited[vw->seqID] = true;
		}
	}
}

void NuEdge::updateCsMove(CsMove& cm){
	this->cmTmp = cm;
	LocalFrame cs1 = nodeA->baseConfTmp->cs1 + cm;
	nodeB->updateCoordinate(cs1);

	for(int i=0;i<edgeListB.size();i++){
		cs1 = edgeListB[i]->nodeA->baseConfTmp->cs1 + edgeListB[i]->cmTmp;
		edgeListB[i]->nodeB->updateCoordinate(cs1);
	}
}

double NuEdge::mutEnergy(){
	return 0.0;
}

void NuEdge::acceptMutation(){
	this->cm = this->cmTmp;
	for(int i=0;i<nodeListB.size();i++){
		nodeListB[i]->acceptCoordMove();
	}
}

void NuEdge::clearMutation(){
	this->cmTmp = this->cm;
	for(int i=0;i<nodeListB.size();i++){
		nodeListB[i]->clearCoordMove();
	}
}

void NuEdge::updateCsMoveCG(CsMove& cm) {
	this->cmTmp = cm;
	LocalFrame cs1 = nodeA->baseConfCGTmp->cs1 + cm;
	nodeB->updateCoordinateCG(cs1);
	for(int i=0;i<edgeListB.size();i++){
		cs1 = edgeListB[i]->nodeA->baseConfCGTmp->cs1 + edgeListB[i]->cmTmp;
		edgeListB[i]->nodeB->updateCoordinateCG(cs1);
	}
}

double NuEdge::mutEnergyCG(){
	return 0.0;
}

void NuEdge::acceptMutationCG(){
	this->cm = this->cmTmp;
	for(int i=0;i<nodeListB.size();i++){
		nodeListB[i]->acceptCoordMoveCG();
	}
}

void NuEdge::clearMutationCG(){
	this->cmTmp = this->cm;
	for(int i=0;i<nodeListB.size();i++){
		nodeListB[i]->clearCoordMoveCG();
	}
}

NuTree::NuTree(NuGraph* graph){
	this->graph = graph;
	this->adjMtx = new bool[graph->len*graph->len];

}

void NuTree::printEdges(){
	for(int i=0;i<geList.size();i++){
		NuEdge* e = geList[i];
		printf("%-3d %-3d %7.3f\n", e->indexA, e->indexB, e->weight);
	}
}

NuTree::~NuTree(){
	delete [] adjMtx;
}

NuGraph::NuGraph(const string& inputFile){

	this->pairLib = new BasePairLib();
	this->rotLib = new RotamerLib();
	this->atLib = new AtomLib();
	this->moveLib = NULL;
	//this->moveLib = new NuPairMoveSetLibrary();

	init(inputFile);
}



NuGraph::~NuGraph() {
	delete this->pairLib;
	delete this->rotLib;
	delete this->atLib;
	delete this->moveLib;

	delete [] seq;
	delete [] wcPairPosID;
	delete [] connectToDownstream;
	delete [] sepTable;
	for(int i=0;i<len;i++){
		delete allNodes[i];
	}
	for(int i=0;i<len*len;i++){
		delete allEdges[i];
	}
	delete [] allNodes;
	delete [] allEdges;
	delete [] connectToDownstream;
	for(int i=0;i<len;i++){
		delete initRiboseRotList[i];
	}

}

void NuGraph::init(const string& inputFile){
	int i,j,k;

	NSPtools::InputParser input(inputFile);

	input.printOptions();

	/*
	 * task:
	 * 		predict: ab initial prediction
	 * 		refinement: fixed cluster type refinement
	 */

	string task = input.getValue("task");
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string cst = input.getValue("cst");
	string chainBreak = input.getValue("break");
	vector<string> templates = input.getMultiValues("template");
	vector<string> templatesAlignA = input.getMultiValues("alignNat");
	vector<string> templatesAlignB = input.getMultiValues("alignTmp");
	vector<string> templatesType = input.getMultiValues("tempType");

	len = baseSeq.length();

	cout << "sequence length: " <<  len << endl;

	this->seq = new int[len];
	this->wcPairPosID = new int[len];
	this->connectToDownstream = new bool[len];
	this->sepTable = new int[len*len];
	this->allNodes = new NuNode*[len];
	this->allEdges = new NuEdge*[len*len];

	cout << "read chain break" << endl;

	/*
	 * read chain break points
	 */
	for(i=0;i<len;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[len-1] = false;

	vector<string> spt;
	splitString(chainBreak, " ", &spt);
	for(i=0;i<spt.size();i++){
		int pos = atoi(spt[i].c_str());
		if(pos >= baseSeq.length() || pos < 0){
			cout << "invalid chain break position: " << pos << endl;
			exit(0);
		}
		connectToDownstream[atoi(spt[i].c_str())] = false;
	}

	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			int ij = i*len+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else sepTable[ij] = 2;
		}
	}

	cout << "init nodes" << endl;

	/*
	 * init NuNodes
	 */
	RNAPDB pdb(pdbFile, "pdbid");
	vector<RNABase*> baseList = pdb.getBaseList();

	if(baseList.size() != baseSeq.length()) {
		cout << "pdb size: " << baseList.size() << " not equal to sequence length: " << baseSeq.length() << endl;
		exit(0);
	}
	else if(baseSeq.length() != baseSec.length()) {
		cout << "seq length not equal to sec length" << endl;
		exit(0);
	}
	else if(cst.length() != baseSeq.length()){
		cout << "invalid cst " << endl;
		cout << cst << endl;
		exit(1);
	}

	for(i=0;i<len;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		this->wcPairPosID[i] = -1;

		RiboseRotamer* rot = new RiboseRotamer();

		if(!baseList[i]->backboneComplete()){
			rot->copyValueFrom(rotLib->riboseRotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt));
		}
		else {
			rot->copyValueFrom(rotLib->riboseRotLib->getNearestRotamer(baseList[i]));
		}

		this->initRiboseRotList.push_back(rot);
	}


	for(i=0;i<len;i++){
		LocalFrame cs1 = baseList[i]->getCoordSystem();
		this->allNodes[i] = new NuNode(i, baseList[i]->baseTypeInt, cs1, initRiboseRotList[i], atLib);
		this->allNodes[i]->graph = this;
	}

	cout << "parse secondary structure" << endl;
	/*
	 * parse secondary structure information
	 */

	char ss[len];
	for(int i=0;i<len;i++){
		ss[i] = baseSec[i];
	}

	map<char,char> brackets;
	brackets[')'] = '(';
	brackets[']'] = '[';
	brackets['}'] = '{';
	brackets['>'] = '<';
	brackets['a'] = 'A';
	brackets['b'] = 'B';
	brackets['c'] = 'C';
	brackets['d'] = 'D';
	brackets['e'] = 'E';
	brackets['f'] = 'F';
	brackets['g'] = 'G';
	brackets['h'] = 'H';
	brackets['i'] = 'I';
	brackets['j'] = 'J';
	map<char,char>::iterator it;

	for(i=0;i<len;i++) {
		char c = ss[i];
		it = brackets.find(c);
		if(it == brackets.end()) continue;
		int preIndex = -1;
		for(j=i-1;j>=0;j--) {
			if(ss[j] == it->second) {
				preIndex = j;
				break;
			}
		}
		if(preIndex < 0) {
			cout << "invalid ssSeq: " << baseSec << endl;
			exit(1);
		}
		ss[i] = '.';
		ss[preIndex] = '.';
		wcPairPosID[i] = preIndex;
		wcPairPosID[preIndex] = i;
	}

	cout << "init edges" << endl;
	/*
	 * init NuEdges
	 */
	for(int i=0;i<len;i++){
		for(int j=0;j<len;j++){
			this->allEdges[i*len+j] = new NuEdge(allNodes[i], allNodes[j], this);
			this->allEdges[i*len+j]->graph = this;
			this->allEdges[i*len+j]->weight = 0.0;
			if(task == "refinement"){
				this->allEdges[i*len+j]->initNativeMoveSet();
			}
		}
	}

	cout << "init constraints" << endl;
	/*
	 * init constaints
	 */

	string cstString = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	for(i=0;i<len;i++){
		char c = cst[i];
		for(j=0;j<cstString.length();j++){
			if(c == cstString[j]) {
				for(k=i+1;k<len;k++){
					char d = cst[k];
					if(c == d){
						this->allEdges[i*len+k]->weight = -999.9;
						this->allEdges[k*len+i]->weight = -999.9;
					}
				}
			}
		}
	}

	for(i=0;i<len;i++){
		for(j=i+1;j<len;j++){
			this->geList.push_back(allEdges[i*len+j]);
		}
	}

	cout << "finish init" << endl;
}

bool MST_cmp_weight(NuEdge* e1, NuEdge* e2){
	return (e1->weightRand < e2->weightRand);
}

int MST_find(int* parent, int f){
	while(parent[f] >= 0){
		f = parent[f];
	}
	return f;
}

void NuGraph::initRandWeight(){
	for(int i=0;i<len;i++){
		for(int j=0;j<len;j++){
			allEdges[i*len+j]->weightRand = allEdges[i*len+j]->weight - rand()*1.0/RAND_MAX;
		}
	}
}

void NuGraph::MST_kruskal(NuTree* outputTree){
	int i,j,n,m,a,b;
	int parent[len];
	for(i=0;i<len;i++){
		parent[i] = -1;
	}

	int ne = geList.size();

	NuEdge* sortedGeList[ne];

	for(i=0;i<ne;i++){
		sortedGeList[i] = geList[i];
	}

	sort(sortedGeList, sortedGeList+ne, MST_cmp_weight);

	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			outputTree->adjMtx[i*len+j] = false;
		}
	}
	outputTree->geList.clear();

	for(i=0;i<ne;i++){


		n = MST_find(parent, sortedGeList[i]->indexA);
		m = MST_find(parent, sortedGeList[i]->indexB);

		if(n != m) {
			a = sortedGeList[i]->indexA;
			b = sortedGeList[i]->indexB;
			parent[n] = m;
			outputTree->adjMtx[a*len+b] = true;
			outputTree->adjMtx[a*len+b] = true;
			if(a<b)
				outputTree->geList.push_back(allEdges[a*len+b]);
			else
				outputTree->geList.push_back(allEdges[b*len+a]);
		}
	}
}

void NuGraph::printAllEdge(){
	for(int i=0;i<len;i++){
		for(int j=i+1;j<len;j++){
			NuEdge* e = allEdges[i*len+j];
			printf("%-3d %-3d %7.3f\n", i, j, e->weight);
		}
	}
}


} /* namespace NSPpredNA */
