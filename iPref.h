#ifndef _IPREF_H_
#define _IPREF_H_

#include "header.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "virtualRNode.h"
#include "filemem.h"
#include "global.h"
#include "skyline.h"

extern unordered_map<long int, RtreeNode*> ramTree;

//compute Top-K result set
vector<long int> computeTopK(const int dim, float* PG[], vector<long int> skyband, vector<float>& weight, int k);

// Is pj dominate pi in traditional sense
bool IsPjdominatePi(const int dimen, float* PG[], long int pi, long int pj);

// test whether pj is incomparable with pi, and pj w > pi w
bool incomparableset(float* PG[], long int pi, long int pj, vector<float>& weight);

// compute the hyperplane of incomparable record pair: pi and pj
vector<float> computePairHP(const int dimen, float* PG[], long int pi, long int pj);

// compute the distance from w to hyperplane HS
float computeDis(vector<float> tmpHS, vector<float> userpref);

// sort pair
bool sortbysec(const pair<long int, float> &a, const pair<long int, float> &b);

// optimized algorithm
float computeRho(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[]);

// unknown x baseline
float computeRho_unknownX_basic(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[]);

// unknown x efficient
float computeRho_unknownX_efficient(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[]);

// compute the radius in phase 2, optimized algorithm
float computeradius(const int k, const int dim, long int pi, vector<float>& userpref, vector<long int>& incompSet, float* PG[]);

#endif