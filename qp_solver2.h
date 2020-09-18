//
// Created by 12859 on 2020/9/17.
//

#ifndef IPREF_QP_SOLVER2_H
#define IPREF_QP_SOLVER2_H
#include <vector>
#include "QuadProg++.hh"
using namespace std;

double qp_solver2(const vector<float>& w, const vector<vector<double>>& H);
double qp_solver2(const vector<float>& w, const vector<vector<double>>& H1, const vector<vector<double>>& H2);

#endif //IPREF_QP_SOLVER2_H
