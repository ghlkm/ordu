
#include "case_study.h"

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
//    int FP=result.size()-TP;
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
    // precision, TP/(TP+FN)
    int TP=0;
    for(int r:result){
        if(groundTruth.find(r)!=groundTruth.end()){
            ++TP;
        }
    }
    return (double)TP/groundTruth.size();
}
