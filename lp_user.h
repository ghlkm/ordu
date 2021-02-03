
#ifndef IPREF_LP_USER_H
#define IPREF_LP_USER_H
#include <vector>
#include <lp_types.h>

using namespace std;

void lpModel(lprec* lp, int& dimen);
bool isFeasible(vector<vector<double>> &r1,  vector<vector<double>> &r2, vector<vector<double>> &r1_r2);
#endif //IPREF_LP_USER_H
