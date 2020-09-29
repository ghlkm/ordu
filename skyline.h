#ifndef _SKYLINE_H_
#define _SKYLINE_H_

#include "header.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "virtualRNode.h"
#include "filemem.h"
#include "global.h"

#define MAXPAGEID 49999999

//define pageid and virtual rtree node pair
typedef multimap<long int, VirtualRNode*>::value_type PintVRN;

//define float and int vaule pair
typedef multimap<float, long int>::value_type PfltINT;

typedef map<int, int>::value_type PintInt;

//define cell range for two dim case
typedef multimap<float, float>::value_type PfltFLT;


float minDist(float p1[], float p2[], int dimen);

void aggregateRecords(Rtree& a_rtree);

int countRecords(Rtree& a_rtree, int pageID);

void checkaggreRtree(Rtree& a_rtree);

bool IsDominatedBy(const int dimen, const float pt[], vector<long> a_skylines, float* PG[]);

bool isDominateBy(const int dimen, const float pt[], Point& focal);

bool isDominateByFocal(const int dimen, const float pt[], Point& focal);

void GetSkylines(const int dimen, Rtree& a_rtree, std::multimap<long int, VirtualRNode*>& NonResultEntry, std::vector<long int>& PrunedNodes, set<long int>& a_skylines, float* PG[]);

int incomparableWindowsHD(Rtree& a_rtree, float* PG[], Point& a_pt, std::multimap<long int, VirtualRNode*>& NonResultEntry, vector<long int>& a_resultID);

int countDominator(Rtree& a_rtree, float* PG[], Point& a_pt, unordered_set<long int>& ignoreRecord);

void updateSkylines(const int dimen, Rtree& a_rtree, unordered_set<long int>& a_skylines, float* PG[], unordered_set<long int>& ignoreRecords);

int countkDominator(const int dimen, const float pt[], vector<long> kskyband, float* PG[]);

void Getkskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, Point& a_pt, float* PG[], const int k);

void updateSkylines_ram(const int dimen, Rtree& a_rtree, set<long int>& a_skylines, float* PG[], multimap<float, int>& heap, multimap<long int, RtreeNodeEntry*>& RecordEntry);

void computeHP(const int dimen, float* PG[], Point& a_pt, unordered_set<long int>& initSkyline, vector<long int>& newAddSL);

void computeHP(const int dimen, float* pG[], Point& a_pt, const int objcnt, vector<long int>& newAddSL);

void computeHPwithoutReduce(const int dimen, float* PG[], Point& a_pt, unordered_set<long int>& initSkyline, vector<long int>& newAddSL);

void computeHPforkskyband(const int dimen, float* PG[], Point& a_pt, vector<long int>& kskyband, vector<long int>& newAddSL);

int countDominee(float* PG[], Point& pt, int rid, const int& objCnt);

bool dominatedByK(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k);

bool dominatedByK_noSL(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k);

void computeHP(const int dimen, float* PG[], Point& a_pt, const int objcnt, vector<long int>& newAddSL, const int NoofConstraints);

// added for utk

void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k);

void onionlayer(vector<long int>& skyband, float* PG[], const int k, vector<long int >& klayers, int& dim);

//void onionlayer(Rtree& a_rtree, float* PG[], const int k, vector<long int >& klayers, int& dim);

void initHS(const int dimen, vector<float>& regions);

void rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);

void initRecordEntry(Rtree& rt, multimap<long int, RtreeNodeEntry*>& RecordEntry);

#endif
