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

template<typename T>
bool dominatedByK(const int dimen, const float pt[], vector<T> &kskyband, float* PG[], int k)
{
    // see if pt is dominated by k options of kskyband
    if (kskyband.empty())
        return false;

    int count = 0;
    for (auto pid : kskyband)
    {
        bool dominated = true;
        for (int i = 0; i < dimen; i++)
        {
            if ((PG[pid][i] + PG[pid][i+dimen])/2 < pt[i])
            {
                dominated = false;
                break;
            }
        }
        if (dominated) {
            count++;
            if(count>=k){
                return true;
            }
        }
    }
    return false;
}

template<typename T>
bool dominatedByK(const int dimen, const float pt[], vector<T> &kskyband, float* PG[], int k, vector<size_t> &ret)
{
    // see if pt is dominated by k options of kskyband
    if (kskyband.empty()) return false;
    int count = 0;
    for (auto pid : kskyband){
        bool dominated = true;
        for (int i = 0; i < dimen; i++){
            if ((PG[pid][i] + PG[pid][i+dimen])/2 < pt[i]){
                dominated = false;
                break;
            }
        }
        if (dominated) {
            ret.push_back(pid);
            count++;
            if(count>=k){
                ret.clear();
                return true;
            }
        }
    }
    return false;
}

bool dominatedByK_noSL(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k);

void kskyband(const int dimen, Rtree& a_rtree, vector<int>& kskyband, float* PG[], const int k);

void rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);

#endif
