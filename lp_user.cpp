//
// Created by 12859 on 2020/9/9.
//

#include "lp_user.h"
#include <lp_lib.h>
#include <vector>
using namespace std;

void lpModel(lprec* lp, int& dimen)
{
    if (lp == nullptr)
    {
        fprintf(stderr, "Unable to create new LP model\n");
        exit(0);
    }
    // constraints in each dimension
    for (int di = 0; di < dimen; di++) {
        set_bounds(lp, di+1, 0.0, 1.0);
    }
}

bool isFeasible(vector<vector<double>> &r1,  vector<vector<double>> &r2){
    int dim;
    if(r1.empty() || r2.empty()){
        return true;
    }else{
        dim=r1[0].size();
    }
    lprec *lp = make_lp(0, dim);
    set_verbose(lp, IMPORTANT);
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
    set_add_rowmode(lp, TRUE);
    lpModel(lp, dim);
    vector<double> tmp(dim+1);
    for (vector<double> &r1i:r1) {
        for (int i = 0; i <dim ; ++i) {
            tmp[i+1]=r1i[i];
        }
        add_constraint(lp, tmp.data(), LE, 0.0);
    }
    for (vector<double> &r2i:r2) {
        for (int i = 0; i <dim ; ++i) {
            tmp[i+1]=r2i[i];
        }
        add_constraint(lp, tmp.data(), LE, 0.0);
    }

    set_add_rowmode(lp, FALSE);
    set_timeout(lp, 1);
    int ret = solve(lp);
    delete_lp(lp);
    return ret<=1;
}
