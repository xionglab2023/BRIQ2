/*
 * SubstitutionMatrix.h
 *
 */

#ifndef MATH_SUBSTITUTIONMATRIX_H_
#define MATH_SUBSTITUTIONMATRIX_H_

#include "tools/StringTool.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include <stdlib.h>
#include <vector>
#include <fstream>


namespace NSPmath {

using namespace std;

class SubstitutionMatrix {
public:
	double subM[20][20];
	SubstitutionMatrix();
	virtual ~SubstitutionMatrix();
};

} /* namespace NSPmath */

#endif /* MATH_SUBSTITUTIONMATRIX_H_ */
