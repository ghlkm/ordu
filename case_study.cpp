
#include "case_study.h"
#include "skyline.h"
#include "math_lib.h"
double jaccard(set<int> &a, set<int>&b){
    // jaccard similarity
    set<int> all(a.begin(), a.end());
    for(int bi:b){
        all.insert(bi);
    }
    double cnt=0;
    for(int bi:b){
        if(a.find(bi)!=a.end()){
            cnt++;
        }
    }
    return cnt/all.size();
}

double jaccard(vector<int>::iterator ab, vector<int>::iterator ae, vector<int>::iterator bb, vector<int>::iterator be){
    // \para ab iterable a begin
    // \para ae iterable a end
    // \para bb iterable b begin
    // \para be iterable b end

    // jaccard similarity
    set<int> a(ab, ae);
    set<int> b(bb, be);
    return jaccard(a, b);
}


double precision(vector<int>::iterator gtb, vector<int>::iterator gte, vector<int>::iterator rb, vector<int>::iterator re){
    set<int> gt(gtb, gte);
    set<int> r(rb, re);
    return precision(gt, r);
}

double precision(set<int> &groundTruth, set<int>&result){
    // TP, true positive
    // FP, false positive
    // precision, TP/(TP+FP)
    int TP=0;
    for(int r:result){
        if(groundTruth.find(r)!=groundTruth.end()){
            ++TP;
        }
    }
    return (double)TP/result.size();
}

double recall(vector<int>::iterator gtb, vector<int>::iterator gte, vector<int>::iterator rb, vector<int>::iterator re){
    set<int> gt(gtb, gte);
    set<int> r(rb, re);
    return recall(gt, r);
}

double recall(set<int> &groundTruth, set<int>&result){
    // TP, true positive
    // FN, false negative
    // recall, TP/(TP+FN)
    int TP=0;
    for(int r:result){
        if(groundTruth.find(r)!=groundTruth.end()){
            ++TP;
        }
    }
    return (double)TP/groundTruth.size();
}


int total_do_cnt(map<int, unordered_set<int>> &do_id, vector<int> &test_this){
    // \para input do_id
    //     do_id[1]={2, 3, 4} means option 1 dominates option 2, 3, and 4
    unordered_set<int> s;
    for (int i:test_this) {
        for (int j: do_id.find(i)->second) {
            s.insert(j);
        }
    }
    return s.size();
}

vector<int> OSS_skyline(int objCnt, int r, Rtree*rtree, float **PointSet, int dim){
    vector<long> one_skyband;
    kskyband(dim, *rtree, one_skyband, PointSet, 1);

    // build dominate relations
    map<int, unordered_set<int>> do_id;
    for (int id:one_skyband) {
        do_id[id]=unordered_set<int>();
    }
    for (int id: one_skyband) {
        for (int i = 1; i <= objCnt; ++i) {
            if(i==id){
                continue;
            }
            if(v1_dominate_v2(PointSet[id], PointSet[i], dim)){
                do_id[id].insert(i);
            }
        }
    }

    // select r from n do combination and check total dominate cnt (TDC)
    // return the the combination with maximal TDC
    // This will generate correct result for OSS-skyline, but a brute force solution
    // This practical fast since 1-skyline is small when d is small
    int n=one_skyband.size();
    std::vector<bool> v(n);
    std::fill(v.begin(), v.begin() + r, true);

    int max_tdo=0;
    vector<int> best_cb;
    do {
        vector<int> one_cb;
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                one_cb.push_back(one_skyband[i]);
            }
        }
        int cur=total_do_cnt(do_id, one_cb);
        if(cur>max_tdo){
            max_tdo=cur;
            best_cb=one_cb;
        }
    } while (std::prev_permutation(v.begin(), v.end())); // forward next permutation
    return best_cb;
}