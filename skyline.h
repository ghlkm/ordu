#ifndef _SKYLINE_H_
#define _SKYLINE_H_

#include "rtree_lib/header.h"
#include "rtree_lib/rtree.h"
#include "rtree_lib/rnode.h"
#include "rtree_lib/rentry.h"
#include "rtree_lib/virtualRNode.h"
#include "rtree_lib/filemem.h"
#include "rtree_lib/global.h"

float minDist(float p1[], float p2[], int dimen);

void aggregateRecords(Rtree& a_rtree);

int countRecords(Rtree& a_rtree, int pageID);

bool dominatedByK(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k);

bool dominatedByK_noSL(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k);

void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k);

void rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);

#endif
