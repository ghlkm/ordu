
#ifndef IPREF_CASE_STUDY_H
#define IPREF_CASE_STUDY_H
#include <vector>
#include <set>
using namespace std;

double jaccard(vector<int>::iterator ab, vector<int>::iterator ae, vector<int>::iterator bb, vector<int>::iterator be);

double jaccard(set<int> &a, set<int>&b);

double precision(vector<int>::iterator gtb, vector<int>::iterator gte, vector<int>::iterator rb, vector<int>::iterator re);

double precision(set<int> &groundTruth, set<int>&result);

double recall(vector<int>::iterator gtb, vector<int>::iterator gte, vector<int>::iterator rb, vector<int>::iterator re);

double recall(set<int> &groundTruth, set<int>&result);

#endif //IPREF_CASE_STUDY_H
