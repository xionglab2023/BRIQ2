/*
 * FragmentAssembling.h
 *
 *  Created on: 2022Äê5ÔÂ11ÈÕ
 *      Author: pengx
 */

#ifndef predNA_FRAGMENTASSEMBLING_H_
#define predNA_FRAGMENTASSEMBLING_H_

#include "model/RNABaseName.h"
#include "model/RiboseRotamerLib.h"
#include "forcefield/PO3Builder.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include <iostream>
#include <set>

#include "forcefield/RnaEnergyTable.h"
#include "predNA/BaseMoveLibrary.h"
#include "predNA/BRConnection.h"
#include "predNA/BRNode.h"
#include "predNA/FragmentLibrary.h"
#include "predNA/MoveMutator.h"


namespace NSPpredna {

using namespace std;
using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPgeometry;


class FragmentAssembling {
public:
	FragmentAssembling();
	virtual ~FragmentAssembling();
};

} /* namespace NSPpredna */

#endif /* predNA_FRAGMENTASSEMBLING_H_ */
