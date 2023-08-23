/*
 * Stat.h
 *
 */

#ifndef MATH_STAT_H_
#define MATH_STAT_H_

#include <vector>
#include <math.h>
#include <algorithm>

namespace NSPmath {

using namespace std;

double sd(vector<double>& dList);
double mid(vector<double>& dList);
double mean(vector<double>& dList);
double sum(vector<double>& dList);
double top(vector<double>& dList, double p);
float top(vector<float>& dList, double p);

} /* namespace NSPmodel */

#endif /* MATH_STAT_H_ */
