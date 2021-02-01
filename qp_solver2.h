// an interface for quadprog

#ifndef IPREF_QP_SOLVER2_H
#define IPREF_QP_SOLVER2_H
#include <vector>
#include "quadProg_lib/QuadProg++.hh"
using namespace std;
vector<float> find_point_in_region(const vector<float>& w, const vector<vector<double>>& H);
vector<double> find_point_in_region(const vector<float>& w, const vector<vector<double>>& H1, const vector<vector<double>>& H2);
double qp_solver2(const vector<float>& w, const vector<vector<double>>& H);
double qp_solver2(const vector<float>& w, const vector<vector<double>>& H1, const vector<vector<double>>& H2);
double qp_solver2(const vector<float>& w, const vector<vector<double>>& H1, int opt, const vector<int>& cmp, float **PG);
double qp_solver2(const vector<float>& w, const vector<vector<double>>& H1, int opt, const vector<int>& cmp, float **PG, vector<double> &retv);
#endif //IPREF_QP_SOLVER2_H
