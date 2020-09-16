//
// Created by 12859 on 2020/9/9.
//

#include "lp_user.h"
#include <lp_lib.h>
#include <vector>
#include <cassert>
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

bool isFeasible(vector<vector<double>> &r1,  vector<vector<double>> &r2, vector<vector<double>> &r1_r2){
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
    vector<double *> store;
    for (vector<double> &r1i:r1) {
        double *tmp=new double[dim+1];
        for (int i = 0; i <dim ; ++i) {
            tmp[i+1]=r1i[i];
        }
        store.push_back(tmp);
        add_constraint(lp, tmp, LE, 0.0);
    }
    for (vector<double> &r2i:r2) {
        double *tmp=new double[dim+1];
        for (int i = 0; i <dim ; ++i) {
            tmp[i+1]=r2i[i];
        }
        store.push_back(tmp);
        add_constraint(lp, tmp, LE, 0.0);
    }
//    double *tmp=new double[dim+1];
//    for (int i = 0; i <dim ; ++i) {
//        tmp[i+1]=1.0;
//    }
//    store.push_back(tmp);
//    add_constraint(lp, tmp, LE, 1.0);
    set_add_rowmode(lp, FALSE);
//    set_presolve(lp, PRESOLVE_ROWS, 0);
//    set_presolve(lp,
//            PRESOLVE_LINDEP |
//            PRESOLVE_IMPLIEDFREE|
//            PRESOLVE_ROWDOMINATE|
//            PRESOLVE_BOUNDS|
//            PRESOLVE_COLS|
//            PRESOLVE_REDUCEGCD, get_presolveloops(lp));
    set_timeout(lp, 1);
    int ret = solve(lp);
    if(ret<=1){
        int ccnt=get_Nrows(lp);
        assert(ccnt>0);
        vector<int> simplify(ccnt);
        for (int i = 1; i <=ccnt ; ++i) {
            simplify[i-1]=get_orig_index(lp, i)-1;
        }
//        vector<int> simplify(r1.size()+r2.size());
//        for (int i = 1; i <=r1.size()+r2.size(); ++i) {
//            simplify[i-1]=i-1;
//        }
        for (int idx: simplify) {
//            if(idx>=r1.size()+r2.size()){
//                continue;
//            }
            if(idx<r1.size()){
                r1_r2.push_back(r1[idx]);
            }else{
                r1_r2.push_back(r2[idx-r1.size()]);
            }
        }
    }
    delete_lp(lp);
    for(double *tmp: store){
        delete [] (tmp);
    }
    return ret<=1;
}
