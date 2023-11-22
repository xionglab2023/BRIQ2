/*
 * AtomLib.h
 *
 *  Created on: 2022Äê8ÔÂ29ÈÕ
 *      Author: pengx
 */

#ifndef MODEL_ATOMLIB_H_
#define MODEL_ATOMLIB_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "model/ResName.h"
#include "model/RNABaseName.h"
#include "geometry/xyz.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"

using namespace std;
using namespace NSPgeometry;
using namespace NSPtools;

namespace NSPmodel {

class AtomProperty {
public:
        string atomUniqueName;
        int uniqueID;

        string atomType;
        int atomTypeInt;

        float vdwRadius;

        int donorIndex;
        int acceptorIndex;

        bool isHDonor;
        bool isHAcceptor;
        bool isIon;

        AtomProperty(const string& line);
        void setAtomTypeInt(int type){
        	this->atomTypeInt = type;
        }

        void setDonorType(int type){
        	this->donorIndex = type;
        }

        void setAcceptorType(int type){
        	this->acceptorIndex = type;
        }

        virtual ~AtomProperty();
};


class AtomLib {
public:

	vector<AtomProperty*> apList;
	vector<string> uniqueNames;
	map<string,int> uniqueNameToUniqueID;

	vector<vector<int>*> aaUniqueIDs;
	vector<vector<string>*> aaAtomNames;
	vector<vector<string>*> aaScAtomNames;

	vector<vector<int>*> baseUniqueIDs;
	vector<vector<string>*> baseAtomNames;
	vector<vector<string>*> baseScAtomNames;

	vector<string> generalAtomNames;
	map<string, int> generalAtomTypeToInt; //index start from 0

	map<string, int> donorAtomMap;
	map<string, int> acceptorAtomMap;

	map<int, int> uniqueIDToDonorID;
	map<int, int> uniqueIDToAcceptorID;

	ResName rn;
	RNABaseName bn;

	/*
	 * donor atom list
	 * 0: O2'
	 * 1: SER-OG
	 * 2: THR-OG1
	 * 3: TYR-OH
	 * 4: HIS-ND1 HIS-NE2
	 * 5: N
	 * 6: ASN-ND2 GLN-NE2
	 * 7: ARG-NH1 ARG-NH2
	 * 8: ARG-NE
	 * 9: LYS-NZ
	 * 10: TRP-NE1
	 * 11: A-N6 DA-N6
	 * 12: U-N3 DT-N3
	 * 13: G-N1 DG-N1
	 * 14: G-N2 DG-N2
	 * 15: C-N4 DC-N4
	 */

	/*
	 * acceptor atom list
	 * 0: O2'
	 * 1: SER-OG
	 * 2: THR-OG1
	 * 3: TYR-OH
	 * 4: HIS-ND1 HIS-NE2
	 * 5: O
	 * 6: ASP-OD1 ASP-OD2
	 * 7: GLU-OE1 GLU-OE2
	 * 8: ASN-OD1 GLN-OE1
	 * 9: OP1 OP2
	 * 10: A-N7 DA-N7
	 * 11: A-N1 DA-N1
	 * 12: A-N3 DA-N3
	 * 13: U-O2 U-O4 DT-O2 DT-O4
	 * 14: G-N7 DG-N7
	 * 15: G-O6 DG-O6
	 * 16: G-N3 DG-N3
	 * 17: C-O2 DC-O2
	 * 18: C-N3 DC-N3
	 */

	//in total 304 hbond types (16*19)
	// protein-protein: 80
	// protein-NA: 158
	// NA-NA: 66

	/*
	 * polar groups
	 * 0: BBN
	 * 1: BBO
	 * 2: ASP
	 * 3: GLU
	 * 4: HIS
	 * 5: LYS
	 * 6: ASN
	 * 7: GLN
	 * 8: ARG
	 * 9: SER
	 * 10: THR
	 * 11: TRP
	 * 12: TYR
	 * 13: A DA
	 * 14: U
	 * 15: G DG
	 * 16: C DC
	 * 17: DT
	 * 18: OP
	 * 19: O2'
	 */

	AtomLib();

	bool atomDefined(const string& uniqueName) const;
	int uniqueNameToID(const string& uniqueName) const;
	string uniqueIDToName(int uniqueID) const;

	vector<int>* getAminoAcidAtomIDs(int intName);
	vector<string>* getAminoAcidAtomNames(int intName);
	vector<string>* getAminoAcidSidechainAtomNames(int intName);
	AtomProperty* getAtomProperty(int uniqueID);
	vector<int> getAASidechainUniqueIDs(int aaType);


	vector<string>* getRnaAtomNames(int type) const;
	vector<string>* getRnaSidechainAtoms(int type) const;


	virtual ~AtomLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_ATOMLIB_H_ */
