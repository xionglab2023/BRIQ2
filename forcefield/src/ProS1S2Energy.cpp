/*
 * ProS1S2Energy.cpp
 *
 */

#include "forcefield/ProS1S2Energy.h"

namespace NSPforcefield {

ProS1S2Energy::ProS1S2Energy(DesignPara* para, int s2ID) {

	this->s2SubID = s2ID;
	string path = NSPdataio::datapath();
	string s1File = path+"S1/all.sa";
	this->pathS2 = path + "S2/";
	ifstream f;
	f.open(s1File.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open file " << s1File << endl;
		exit(1);
	}

	string line;
	while((getline(f,line))){
		saList.push_back(AAScoreArray(line));
	}
	f.close();

	vector<string> spt;
	string repNumFile = path + "S2/repNum";
	f.open(repNumFile.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open file " << repNumFile << endl;
		exit(1);
	}
	while((getline(f,line))){
		splitString(line, " ", &spt);
		repNums[spt[0]] = atoi(spt[1].c_str());
		keys.push_back(spt[0]);
		pairMap[spt[0]] = vector<ResPairInfo*>();
		saiMap[spt[0]] = vector<SaiPair*>();
	}
	f.close();

	this->wtS1 = para->wS1;

	for(int i=0;i<45;i++){
		this->wtS2[i] = para->wS245[i];
		printf("S2 sep weight: %5.3f\n", this->wtS2[i]);
	}

	for(int i=0;i<saList.size();i++){
		saList[i].multipy(wtS1);
	}

	loadAllRepPoints();
	loadAllSaiPoints();
}

void ProS1S2Energy::loadRepPoints(const string& key, const string& fileName){

	ifstream f;
	f.open(fileName.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open " << fileName << endl;
		exit(1);
	}
	string line;
	getline(f, line);
	while((getline(f, line))){
		ResPairInfo* rp = new ResPairInfo(line);
		this->pairMap[key].push_back(rp);
	}
	f.close();
}

void ProS1S2Energy::loadAllRepPoints(){
	string listFile = pathS2 + "pairRep/list";
	ifstream f;
	f.open(listFile.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open " << listFile << endl;
		exit(1);
	}
	string line, fileName, key;
	while(getline(f, line)){
		fileName = pathS2 + "pairRep/" + line;
		key = line.substr(0,3);
		loadRepPoints(key, fileName);
	}
	f.close();
}

void ProS1S2Energy::loadSaiPoints(const string& key, const string& fileName){
	ifstream f;
	f.open(fileName.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open " << fileName << endl;
		exit(1);
	}
	string line;
	float saiA, saiB;
	vector<string> spt;
	while(getline(f, line)){
		splitString(line, " ", &spt);
		saiA = atof(spt[0].c_str());
		saiB = atof(spt[1].c_str());
		SaiPair* sp = new SaiPair(saiA, saiB);
		this->saiMap[key].push_back(sp);
	}
	//printf("load sasa point: %s %2d\n", key.c_str(), this->saiMap[key].size());
	f.close();
}

void ProS1S2Energy::loadAllSaiPoints(){
	string listFile = pathS2 + "saiRep/list";
	ifstream f;
	f.open(listFile.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open " << listFile << endl;
		exit(1);
	}
	string line, fileName, key;
	while(getline(f, line)){
		fileName = pathS2 + "saiRep/" + line;
		key = line.substr(0,3);
		loadSaiPoints(key, fileName);
	}
	f.close();
}

map<int, double> ProS1S2Energy::findSaiIndexList(ResPairInfo* rp) {
	SaiPair sp(rp->saiA, rp->saiB);
	vector<SaiPair*> spList = saiMap[rp->key];

	double minD = 999999.9;
	int saiIndex = 0;
	int i,j, k;
	double d;
	int n = 5;
	int indexList[n];
	double distList[n];
	for(int i=0;i<n;i++) {
		indexList[i] = -1;
		distList[i] = 9999999.9;
	}

	//cout << rp->key << endl;
	//printf("sp size: %d\n", spList.size());

	for(int i=0;i<spList.size();i++){
		d = spList[i]->distanceSquare(&sp);
		if(d < 0.00001) d = 0.00001;
		if(d < distList[n-1]){
			distList[n-1] = d;
			indexList[n-1] = i;
		}
		else
			continue;

		/*
		 * bubble sort
		 */
        for(int j=n-2;j>=0;j--) {
            if(distList[j+1] < distList[j]) {
                d = distList[j];
                distList[j] = distList[j+1];
                distList[j+1] = d;
                k = indexList[j];
                indexList[j] = indexList[j+1];
                indexList[j+1] = k;
            }
            else
                break;
        }
	}
	map<int,double> saiMap;
	for(int i=0;i<n;i++){
		saiMap[indexList[i]] = distList[i];
	}
	return saiMap;
}

map<int, double> ProS1S2Energy::findRPIndexList(ResPairInfo* rp) {
	vector<ResPairInfo*> rpList = pairMap[rp->key];
	double minD = 999999.9;
	int rpIndex = 0;
	int i,j, k;
	double d;
	int n = 5;
	int indexList[n];
	double distList[n];
	for(int i=0;i<n;i++) {
		indexList[i] = -1;
		distList[i] = 9999999.9;
	}

	for(int i=0;i<rpList.size();i++){
		d = rpList[i]->distance(rp);
		if(d < 0.00001) d = 0.00001;
		if(d < distList[n-1]){
			distList[n-1] = d;
			indexList[n-1] = i;
		}
		else
			continue;

		/*
		 * bubble sort
		 */
        for(int j=n-2;j>=0;j--) {
            if(distList[j+1] < distList[j]) {
                d = distList[j];
                distList[j] = distList[j+1];
                distList[j+1] = d;
                k = indexList[j];
                indexList[j] = indexList[j+1];
                indexList[j+1] = k;
            }
            else
                break;
        }
	}
	map<int,double> rpoMap;
	for(int i=0;i<n;i++){
		rpoMap[indexList[i]] = distList[i];
	}
	return rpoMap;
}

AAScoreArray ProS1S2Energy::getS1(ResInfo* ri){
	return this->saList[ri->repIndex];
}

AAScoreMatrix ProS1S2Energy::getS2(ResPairInfo* rp){
	map<int, double> saiMap = findSaiIndexList(rp);
	map<int, double> rpoMap = findRPIndexList(rp);
	map<int,double>::iterator it1,it2;

	double sigmaRPO = 1.5;
	double sigmaSAI = 1.0;
	if(rp->sep == 5)
		sigmaRPO = 2.0;

	double wtSum = 0;
	double wt;
	for(it1=saiMap.begin();it1!=saiMap.end();++it1){
		for(it2=rpoMap.begin();it2!=rpoMap.end();++it2){
			wt = pow(it2->second,-sigmaRPO) * pow(it1->second,-sigmaSAI);
			wtSum += wt;
		}
	}

	double* matrix = new double[400];
	for(int i=0;i<400;i++){
		matrix[i] = 0.0;
	}

	int repNum, binNum, rpIndex, fileIndex, saiIndex;
	char xx[200];
	string line;
	vector<string> spt;
	int i,j;
	ifstream f;

	for(it1=saiMap.begin();it1!=saiMap.end();++it1){
		for(it2=rpoMap.begin();it2!=rpoMap.end();++it2){
			wt = pow(it2->second,-sigmaRPO) * pow(it1->second,-sigmaSAI) /wtSum * wtS2[rp->keyID];
			//printf("sai %2d w1: %5.3f rpo %2d w2: %5.3f weight: %6.4f\n", it1->first, it1->second, it2->first, it2->second, wt);
			saiIndex = it1->first;
			rpIndex = it2->first;
			repNum = this->repNums[rp->key];
			binNum = repNum/10;
			fileIndex = rpIndex/binNum;

			sprintf(xx, "%ssub%d/%s-%d-%d", this->pathS2.c_str(), this->s2SubID, rp->key.c_str(), saiIndex, fileIndex);
			string file = string(xx);
			int start = rpIndex%binNum*22+2;


			f.open(file.c_str(), ios::in);
			if(!f.is_open()){
				cout << "fail to open file: " << file << endl;
				exit(1);
			}

			//ignore lines
			for(i=0;i<start;i++){
				getline(f,line);
			}

			//read matrix
			for(i=0;i<20;i++){
				getline(f,line);
				splitString(line, " ", &spt);
				for(j=0;j<20;j++){
					matrix[i*20+j] += atof(spt[j].c_str())*wt;
				}
			}
			f.close();
		}
	}
	AAScoreMatrix sm(matrix, 1.0);
	delete matrix;
	return sm;
}

AAScoreMatrix ProS1S2Energy::getS2NearestNb(ResPairInfo* rp){
	vector<ResPairInfo*> rpList = pairMap[rp->key];
	double minD = 99.9;
	int rpIndex = 0;
	int i,j;
	double d;
	for(int i=0;i<rpList.size();i++){
		d = rpList[i]->distance(rp);
		if(d < minD){
			minD = d;
			rpIndex = i;
		}
	}

	SaiPair sp(rp->saiA, rp->saiB);
	vector<SaiPair*> spList = saiMap[rp->key];
	minD = 99.9;
	int saiIndex;
	for(i=0;i<25;i++){
		d = spList[i]->distanceSquare(&sp);
		if(d < minD){
			minD = d;
			saiIndex = i;
		}
	}

	int repNum = this->repNums[rp->key];
	int binNum = repNum/10;
	int fileIndex = rpIndex/binNum;


	char xx[200];
	sprintf(xx, "%ssub%d/%s-%d-%d", this->pathS2.c_str(), this->s2SubID, rp->key.c_str(), saiIndex, fileIndex);
	string file = string(xx);
	int start = rpIndex%binNum*22+2;

	ifstream f;
	f.open(file.c_str(), ios::in);
	if(!f.is_open()){
		cout << "fail to open file: " << file << endl;
		exit(1);
	}

	string line;
	for(i=0;i<start;i++){
		getline(f,line);
	}

	float wt = wtS2[rp->keyID];

	vector<string> spt;
	double* matrix = new double[400];
	for(i=0;i<20;i++){
		getline(f,line);
		splitString(line, " ", &spt);
		for(j=0;j<20;j++){
			matrix[i*20+j] = atof(spt[j].c_str())*wt;
		}
	}
	f.close();
	AAScoreMatrix sm(matrix, 1.0);
	delete matrix;
	return sm;
}


ProS1S2Energy::~ProS1S2Energy() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<keys.size();i++){
		vector<ResPairInfo*> rps = pairMap[keys[i]];
		vector<SaiPair*> sps = saiMap[keys[i]];
		for(int j=0;j<rps.size();j++){
			delete rps[j];
		}
		for(int j=0;j<sps.size();j++){
			delete sps[j];
		}
	}


}

} /* namespace NSPmath */
