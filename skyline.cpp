#include "skyline.h"

extern unordered_map<long int, RtreeNode*> ramTree;

float minDist(float p1[], float p2[], int dimen)
{
	float mindist = 0;

	for (int i = 0; i < dimen; i++)
	{
		float dist = p1[i] - p2[i];
		mindist += (dist * dist);
	}
	return (float)sqrt(mindist);
}

bool dominatedByK(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k)
{
    // see if pt is dominated by k options of kskyband
    if (kskyband.empty())
        return false;

    int count = 0;
    for (vector<long>::iterator iter = kskyband.begin(); iter != kskyband.end(); iter++)
    {
        long pid = *iter;
        bool dominated = true;
        for (int i = 0; i < dimen; i++)
        {
            if (PG[pid][i] + SIDELEN < pt[i])
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

bool dominatedByK_noSL(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k)
{
    if (kskyband.empty())
        return false;

    int count = 0;
    for (vector<long>::iterator iter = kskyband.begin(); iter != kskyband.end(); iter++)
    {
        long pid = *iter;
        bool dominated = true;
        for (int i = 0; i < dimen; i++)
        {
            if (PG[pid][i]< pt[i])
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

void aggregateRecords(Rtree& a_rtree)
{
	queue<long int> H;
	RtreeNode* node;
	H.push(a_rtree.m_memory.m_rootPageID);
	long int pageID;

	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = a_rtree.m_memory.loadPage(pageID);

		if (node->isLeaf() == false)
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				node->m_entry[i]->num_records = countRecords(a_rtree, node->m_entry[i]->m_id);
				H.push(node->m_entry[i]->m_id);
			}
		}
		a_rtree.m_memory.writePage(pageID, node);
	}
}

int countRecords(Rtree& a_rtree, int pageID)
{
	int sumRecords = 0;
	RtreeNode* node = a_rtree.m_memory.loadPage(pageID);
	if (node->isLeaf())
	{
		sumRecords = node->m_usedspace;
	}
	else
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			sumRecords += countRecords(a_rtree, node->m_entry[i]->m_id);
		}
	}
	delete node;
	return sumRecords;
}

void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k)
{
	RtreeNode* node;
	multimap<float, int> heap;
	multimap<float, int>::iterator heapIter;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;

	heap.emplace(INFINITY, a_rtree.m_memory.m_rootPageID);

	while (!heap.empty())
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			for (int d = 0; d < dimen; d++)
				pt[d] = (PG[pageID - MAXPAGEID][d] + PG[pageID - MAXPAGEID][d + dimen])/2;
            if (!dominatedByK(dimen, pt, kskyband, PG, k))
			{
				kskyband.push_back(pageID - MAXPAGEID);
			}
		}
		else
		{
			//node = a_rtree.m_memory.loadPage(pageID);
			node = ramTree[pageID];
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
					{
						pt[d] = node->m_entry[i]->m_hc.getLower()[d] + SIDELEN;
					}
					if (!dominatedByK(dimen, pt, kskyband, PG, k))
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.emplace(mindist, node->m_entry[i]->m_id + MAXPAGEID);
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
                    if (!dominatedByK(dimen, pt, kskyband, PG, k))
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.emplace(mindist, node->m_entry[i]->m_id);
					}
				}
			}
		}
	}
}

void rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree)
{
	ramTree.clear();
	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = rt.m_memory.loadPage(pageID);
		ramTree[pageID] = node;
		if (node->isLeaf() == false)
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				H.push(node->m_entry[i]->m_id);
			}
		}
	}
}
