#ifndef _IPREF_H_
#define _IPREF_H_

#include "header.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "virtualRNode.h"
#include "filemem.h"
#include "global.h"

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

#endif