
#include "case_study.h"

double dist(set<int> &a, set<int>&b){
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

double dist(vector<int>::iterator ab, vector<int>::iterator ae, vector<int>::iterator bb, vector<int>::iterator be){
    set<int> a(ab, ae);
    set<int> b(bb, be);
    return dist(a, b);
}

