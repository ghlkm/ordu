#include "skyline.h"

extern vector<vector<float>> HalfSpaces;
extern unordered_map<long int, long int> RecordIDtoHalfPlaneID;
extern double totalIO;
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

bool IsDominatedBy(const int dimen, const float pt[], vector<long> a_skylines, float* PG[])
{
	vector<long>::iterator iter;
	if (a_skylines.size() == 0)
		return false;

	for (iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		long pid = *iter;
		float s[MAXDIMEN];
		bool dominated = true;
		for (int i = 0; i < dimen; i++)
		{
			if (PG[pid][i] + SIDELEN < pt[i])
			{
				dominated = false;
				break;
			}
		}
		if (dominated)
			return dominated;
	}
	return false;
}

bool isDominateByFocal(const int dimen, const float pt[], Point& focal)
{
	bool dominated = true;
	for (int i = 0; i < dimen; i++)
	{
		if (focal.m_coor[i] + SIDELEN < pt[i])
		{
			dominated = false;
			break;
		}
	}
	return dominated;
}

int countkDominator(const int dimen, const float pt[], vector<long> kskyband, float* PG[])
{
	vector<long>::iterator iter;
	if (kskyband.size() == 0)
		return false;

	int count = 0;
	for (iter = kskyband.begin(); iter != kskyband.end(); iter++)
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
		if (dominated)
			count++;
	}
	return count;
}

bool dominatedByK(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k)
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

void checkaggreRtree(Rtree& a_rtree)
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
				cout << node->m_entry[i]->num_records << endl;
				H.push(node->m_entry[i]->m_id);
			}
		}
		delete node;
	}
}

int countDominee(float* PG[], Point& pt, int rid, const int& objCnt)
{
	float rd[MAXDIMEN];
	int dimen = pt.m_dimen;
	int ret = 0;
	for (int i = 0; i < objCnt; i++)
	{
		if (i == rid - 1)
		{
			continue;
		}
		for (int d = 0; d < dimen; d++)
		{
			rd[d] = (PG[i + 1][d] + PG[i + 1][d + dimen]) / 2;
		}

		if (dominateRecords(pt.m_coor, rd, dimen))
		{
			ret++;
		}
	}
	return ret;
}

int countDominator(Rtree& a_rtree, float* PG[], Point& a_pt, unordered_set<long int>& ignoreRecord)
{
	queue<long int> H;
	float pt[MAXDIMEN];
	float rd[MAXDIMEN];
	int dimen = a_rtree.m_dimen;
	int NoOfDominators = 0;
	int NoofDominees = 0;
	RtreeNode* node;
	RtreeNodeEntry* e0;
	long int NegPageid, pageid;

	float cl[MAXDIMEN], cu[MAXDIMEN];
	for (int d = 0; d < dimen; d++)
	{
		cl[d] = 0;
		cu[d] = a_pt[d];
	}
	Hypercube Dominee_hc(dimen, cl, cu);
	for (int d = 0; d < dimen; d++)
	{
		cl[d] = a_pt[d];
		cu[d] = 1;
	}
	Hypercube Dominator_hc(dimen, cl, cu);

	H.push(a_rtree.m_memory.m_rootPageID);
	while (!H.empty())
	{
		pageid = H.front();
		H.pop();

		node = a_rtree.m_memory.loadPage(pageid);
		if (node->isLeaf())
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				for (int d = 0; d < dimen; d++)
					rd[d] = (PG[node->m_entry[i]->m_id][d] + PG[node->m_entry[i]->m_id][dimen + d]) / 2.0;
				Point tmpPt(dimen, rd);
				if (tmpPt == a_pt)
					continue;
				else if (Dominator_hc.enclose(tmpPt) == true)   //current point lies in the dominator window
				{
					NoOfDominators++;
					ignoreRecord.insert(node->m_entry[i]->m_id);
				}
				else if (Dominee_hc.enclose(tmpPt))
				{
					NoofDominees++;
				}
			}
		}
		else
		{
			e0 = node->genNodeEntry();
			if (Dominator_hc.enclose(e0->m_hc) == true)
			{
				for (int i = 0; i < node->m_usedspace; i++)
					NoOfDominators += node->m_entry[i]->num_records;
			}
			else if (Dominee_hc.enclose(e0->m_hc) == false)
			{
				for (int i = 0; i < node->m_usedspace; i++)
					H.push(node->m_entry[i]->m_id);
			}
			else if (Dominee_hc.enclose(e0->m_hc) == true)
			{
				for (int i = 0; i < node->m_usedspace; i++)
					NoofDominees += node->m_entry[i]->num_records;
			}
			delete e0;
		}
		delete node;
	}
	return NoOfDominators;
}

void updateSkylines(const int dimen, Rtree& a_rtree, unordered_set<long int>& a_skylines, float* PG[], unordered_set<long int>& ignoreRecords)
{
	vector<long int> newskyline;
	RtreeNode* node;
	multimap<float, int> heap;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;
	multimap<float, int>::iterator heapIter;

	
	for (auto iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
		newskyline.push_back(*iter);

	heap.insert(PfltINT(INFINITY, a_rtree.m_memory.m_rootPageID));
	while (!heap.empty())
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			for (int d = 0; d < dimen; d++)
				pt[d] = (PG[pageID - MAXPAGEID][d] + PG[pageID - MAXPAGEID][d + dimen]) / 2;
			if (IsDominatedBy(dimen, pt, newskyline, PG) == false)
			{
				newskyline.push_back(pageID - MAXPAGEID);
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					if (ignoreRecords.find(node->m_entry[i]->m_id) == ignoreRecords.end())
					{
						for (int d = 0; d < dimen; d++)
							pt[d] = node->m_entry[i]->m_hc.getLower()[d] + SIDELEN;
						if (IsDominatedBy(dimen, pt, newskyline, PG) == false)
						{
							mindist = minDist(pt, ORIGNIN, dimen);
							heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
						}
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
					if (IsDominatedBy(dimen, pt, newskyline, PG) == false)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
					}
				}
			}
		}
	}

	for (int i = 0; i < newskyline.size(); i++)
	{
		a_skylines.insert(newskyline[i]);
	}
}

void Getkskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, Point& a_pt, float* PG[], const int k)
{
	RtreeNode* node;
	multimap<long int, RtreeNodeEntry*> RecordEntry;
	multimap<long int, RtreeNodeEntry*>::iterator recordIter;
	multimap<float, int> heap;
	int NegPageid;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;
	multimap<float, int>::iterator heapIter;

	node = a_rtree.m_memory.loadPage(a_rtree.m_memory.m_rootPageID);
	totalIO++;
	if (node->isLeaf())
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
			
			NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
			RtreeNodeEntry* Nentry = node->m_entry[i]->clone();
			RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
		}
	}
	else
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
		}
	}
	delete node;

	while (heap.size() != 0)
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			recordIter = RecordEntry.find(pageID);
			if (recordIter != RecordEntry.end())
			{
				for (int d = 0; d < dimen; d++)
				{
					pt[d] = recordIter->second->m_hc.getLower()[d] + SIDELEN;
				}
				if (countkDominator(dimen, pt, kskyband, PG) <= k && isDominateByFocal(dimen, pt, a_pt) == false)
				{
					kskyband.push_back(pageID - MAXPAGEID);
				}
			}
			else
			{
				cout << "an error incured in Getkskyband" << endl;
				exit(0);
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
					{
						pt[d] = node->m_entry[i]->m_hc.getLower()[d] + SIDELEN;
					}
					if (countkDominator(dimen, pt, kskyband, PG) <= k) 
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
						
						NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
						RtreeNodeEntry* Nentry = node->m_entry[i]->clone();
						RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
					}
				}
			}
			delete node;
		}
	}
}

void updateSkylines_ram(const int dimen, Rtree& a_rtree, set<long int>& a_skylines, float* PG[], multimap<float, int>& heap, multimap<long int, RtreeNodeEntry*>& RecordEntry)
{
	vector<long int> newskyline;
	RtreeNode* node;
	RtreeNodeEntry* e0;
	multimap<long int, RtreeNodeEntry*>::iterator recordIter;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;
	multimap<float, int>::iterator heapIter;

	for (set<long int>::iterator iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		newskyline.push_back(*iter);
	}

	for (multimap<long int, RtreeNodeEntry*>::iterator iter = RecordEntry.begin(); iter!= RecordEntry.end(); iter++)
	{
		for (int j = 0; j < dimen; j++)
			pt[j] = iter->second->m_hc.getLower()[j] + SIDELEN;
		mindist = minDist(pt, ORIGNIN, dimen);
		heap.insert(PfltINT(mindist, iter->first));
	}

	while (heap.size() != 0)
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			recordIter = RecordEntry.find(pageID);
			if (recordIter != RecordEntry.end())
			{
				for (int d = 0; d < dimen; d++)
					pt[d] = recordIter->second->m_hc.getLower()[d] + SIDELEN;
				if (IsDominatedBy(dimen, pt, newskyline, PG) == false)
				{
					newskyline.push_back(pageID - MAXPAGEID);
				}
			}
		}
		else
		{
			cout << "please check where are the bugs" << endl;
		}
	}

	a_skylines.clear();
	for (int i = 0; i < newskyline.size(); i++)
	{
		recordIter = RecordEntry.find(newskyline[i] + MAXPAGEID);
		if (recordIter != RecordEntry.end())
		{
			RecordEntry.erase(recordIter);
		}
		a_skylines.insert(newskyline[i]);
	}

	//delete node;
	//delete e0;
}

void initRecordEntry(Rtree& a_rtree, multimap<long int, RtreeNodeEntry*>& RecordEntry)
{
	queue<long int> H;
	int dimen = a_rtree.m_dimen;
	int NoOfDominators = 0;
	int NoofDominees = 0;
	RtreeNode* node;
	RtreeNodeEntry* e0;
	long int NegPageid, pageid;

	H.push(a_rtree.m_memory.m_rootPageID);
	while (!H.empty())
	{
		pageid = H.front();
		H.pop();

		node = a_rtree.m_memory.loadPage(pageid);
		if (node->isLeaf())
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
				RtreeNodeEntry* Nentry = node->m_entry[i]->clone();
				RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
			}
		}
		else
		{
			for (int i = 0; i < node->m_usedspace; i++)
				H.push(node->m_entry[i]->m_id);
		}
	}
}

void initHS(const int dimen, vector<float>& regions)
{
	for (int i = 0; i < regions.size(); i++)
	{
		vector<float> tmpHS;
		int loc = i / 2;
		int flag = i % 2;

		for (int di = 0; di < dimen - 1; di++)
		{
			if (di == loc)
			{
				if (flag == 1)
				{
					tmpHS.push_back(1.0);
				}
				else
				{
					tmpHS.push_back(-1.0);
				}
				
			}
			else
			{
				tmpHS.push_back(0.0);
			}
		}
		if (flag == 1)
		{
			tmpHS.push_back(regions[i]);
		}
		else
		{
			tmpHS.push_back(-regions[i]);
		}
		tmpHS.push_back(-(1 + i)); // init R extent as -1
		HalfSpaces.push_back(tmpHS); // store to Halfspaces
		RecordIDtoHalfPlaneID.insert(PintInt(-(1 + i), HalfSpaces.size()));
	}
}

void computeHP(const int dimen, float* PG[], Point& a_pt, unordered_set<long int>& initSkyline, vector<long int>& newAddSL)
{
	newAddSL.clear();
	for (auto sIter = initSkyline.begin(); sIter != initSkyline.end(); sIter++)
	{
		long int id = (*sIter);
		if (RecordIDtoHalfPlaneID.find(id) == RecordIDtoHalfPlaneID.end()) // if this record has hyperplane in halfspaces, then we skip it.
		{
			float rd[MAXDIMEN];
			int tmpDim = 0;
			for (int d = 0; d < dimen; d++)
			{
				rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;
				if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
			}
			if (tmpDim == dimen)
			{
				continue;
			}

			// Reduce d dimension record to d-1 halfspace;
			vector<float> tmpHS;
			float rd_d = rd[dimen - 1];
			float p_d = a_pt[dimen - 1];
			for (int d = 0; d < dimen - 1; d++)
			{
				tmpHS.push_back((rd[d] - rd_d) - (a_pt[d] - p_d));
			}
			tmpHS.push_back(p_d - rd_d);
			tmpHS.push_back(id);            //store the ID of incomparable record in HalfSpace
			HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p	
			RecordIDtoHalfPlaneID.insert(PintInt(id, HalfSpaces.size()));
			newAddSL.push_back(id);
		}
	}
}

void computeHPwithoutReduce(const int dimen, float* PG[], Point& a_pt, unordered_set<long int>& initSkyline, vector<long int>& newAddSL)
{
	newAddSL.clear();
	for (auto sIter = initSkyline.begin(); sIter != initSkyline.end(); sIter++)
	{
		long int id = (*sIter);
		if (RecordIDtoHalfPlaneID.find(id) == RecordIDtoHalfPlaneID.end()) // if this record has hyperplane in halfspaces, then we skip it.
		{
			float rd[MAXDIMEN];
			int tmpDim = 0;
			for (int d = 0; d < dimen; d++)
			{
				rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;
				if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
			}
			if (tmpDim == dimen)
			{
				continue;
			}

			vector<float> tmpHS;
			for (int d = 0; d < dimen; d++)
			{
				tmpHS.push_back(rd[d] -a_pt[d]);
			}
			tmpHS.push_back(0);
			tmpHS.push_back(id);            //store the ID of incomparable record in HalfSpace
			HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p	
			RecordIDtoHalfPlaneID.insert(PintInt(id, HalfSpaces.size()));
			newAddSL.push_back(id);
		}
	}
}

int incomparableWindowsHD(Rtree& a_rtree, float* PG[], Point& a_pt, std::multimap<long int, VirtualRNode*>& NonResultEntry, vector<long int>& a_resultID)
{
	multimap<long int, VirtualRNode*> DataEntryQ;  //queue storing rtree node which is simply a data entry
	multimap<long int, VirtualRNode*>::iterator IntVRN_Iter;
	vector<long int> H;
	bool isAnObject;
	float pt[MAXDIMEN];
	int dimen = a_rtree.m_dimen;
	long int NoOfDominators = 0;
	long int NoofDominees = 0;

	float cl[MAXDIMEN], cu[MAXDIMEN];
	for (int d = 0; d < dimen; d++)
	{
		cl[d] = 0;
		cu[d] = a_pt[d];
	}
	Hypercube Dominee_hc(dimen, cl, cu);
	for (int d = 0; d < dimen; d++)
	{
		cl[d] = a_pt[d];
		cu[d] = 1;
	}
	Hypercube Dominator_hc(dimen, cl, cu);

	H.push_back(a_rtree.m_memory.m_rootPageID);
	while (!H.empty())
	{
		long int pageid = H.back();
		H.pop_back();

		isAnObject = false;

		VirtualRNode* VirNode = new VirtualRNode;
		RtreeNodeEntry* e0;// create node entry e0 for node n, so that its MBR can be obtained   
		if (pageid >= MAXPAGEID) //current element in H is a data entry node (its pageid equals to m_id+MAXPAGED)
		{
			isAnObject = true;

			IntVRN_Iter = DataEntryQ.find(pageid);
			if (IntVRN_Iter == DataEntryQ.end())
			{
				cout << "Error! there is no node " << pageid << " in DataEntryQ!" << endl;
			}
			else
			{
				VirNode->copyData(IntVRN_Iter->second);
				if (VirNode->m_usedspace > 1)
				{
					cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
					return -1;
				}
				else
					e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node           	

				delete IntVRN_Iter->second; //free the virtualRnode; Note, this operation is VERY important
				DataEntryQ.erase(IntVRN_Iter);
			}
			pageid = pageid - MAXPAGEID;
		}
		else
		{
			RtreeNode* node = a_rtree.m_memory.loadPage(pageid);
			VirNode->copyData(*node);
			VirNode->copyEntries(*node, node->m_usedspace);
			e0 = node->genNodeEntry(); //compute the enclosing MBR for current index node        
			delete node;
		}

		if (VirNode->isLeaf())
		{
			if (VirNode->m_usedspace > 1)//current node is a leaf node, so all its data entries must be put into priority queue H
			{
				bool intersected = false;
				bool inDomineeWindow = false;
				bool inDominatorWindow = false;
				if (Dominator_hc.enclose(e0->m_hc) == true)
				{
					inDominatorWindow = true;
					NoOfDominators = NoOfDominators + VirNode->m_usedspace;
				}
				else if (Dominator_hc.isIntersected(Dominator_hc, e0->m_hc) == true)
				{
					intersected = true;
				}
				if (Dominee_hc.enclose(e0->m_hc) == true)
				{
					inDomineeWindow = true;
					NoofDominees = NoofDominees + VirNode->m_usedspace;
				}
				else if (Dominee_hc.isIntersected(Dominee_hc, e0->m_hc) == true)
				{
					intersected = true;
				}

				if (intersected)
				{
					for (int i = 0; i < VirNode->m_usedspace; i++)
					{
						long int NegPageid = VirNode->m_entry[i]->m_id + MAXPAGEID;
						VirtualRNode* node = new VirtualRNode; //for each data entries of n, insert it into data-entry queue
						RtreeNodeEntry* Nentry = VirNode->m_entry[i]->clone();
						node->insertEntry(Nentry);
						DataEntryQ.insert(PintVRN(NegPageid, node));
						delete Nentry;
						H.push_back(NegPageid);
					}
				}
				else if (!inDomineeWindow && !inDominatorWindow)
				{
					a_resultID.push_back(pageid);  //keep the id of index node
				}
			}
			else    //current node is a data entry node
			{
				long int id = VirNode->m_entry[0]->m_id;
				if (id != pageid) cout << "error!! page ids are not equal!" << endl;

				float rd[MAXDIMEN];
				for (int d = 0; d < dimen; d++)
					rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;

				Point tmpPt(dimen, rd);

				if (tmpPt == a_pt)
				{
					continue;
				}

				if (Dominator_hc.enclose(tmpPt) == true)   //current point lies in the dominator window
				{
					NoOfDominators++;
				}
				else
				{
					if (Dominee_hc.enclose(tmpPt) == false)  //current point lies in some incomparable window
					{
						VirtualRNode* vNode = new VirtualRNode;
						vNode->copyData(VirNode);
						NonResultEntry.insert(PintVRN(id + MAXPAGEID, vNode));
					}
					else
					{
						NoofDominees++;
					}
				}
			}
		}
		else//current entry is an index, so insert its entries into the heap              
		{
			bool intersected = false;
			bool inDomineeWindow = false;
			bool inDominatorWindow = false;

			if (Dominee_hc.enclose(e0->m_hc) == true)
				inDomineeWindow = true;
			else if (Dominee_hc.isIntersected(Dominee_hc, e0->m_hc) == true)
				intersected = true;
			if (Dominator_hc.enclose(e0->m_hc) == true)
				inDominatorWindow = true;
			else if (Dominator_hc.isIntersected(Dominator_hc, e0->m_hc) == true)
				intersected = true;

			if (intersected)
			{
				for (int i = 0; i < VirNode->m_usedspace; i++)
				{
					if (Dominee_hc.enclose(VirNode->m_entry[i]->m_hc) == false)
					{
						H.push_back(VirNode->m_entry[i]->m_id);
					}
				}
			}
			else if (!inDomineeWindow && !inDominatorWindow)
			{
				a_resultID.push_back(pageid);  //keep the id of index node              
			}
		}//end check whether it's leaf or index node
		delete VirNode;
		delete e0;
	}
	return NoOfDominators;
}

void GetSkylines(const int dimen, Rtree& a_rtree, std::multimap<long int, VirtualRNode*>& NonResultEntry, std::vector<long int>& PrunedNodes, set<long int>& a_skylines, float* PG[])
{
	multimap<float, long int> H0;
	vector<long int> skylines;

	float pt[MAXDIMEN];
	bool isAnObject;
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	RtreeNodeEntry* e0;

	set<long int>::iterator s0Iter, s1Iter;
	multimap<long int, VirtualRNode*>::iterator nretIter;
	vector<long int>::iterator pIter;

	for (pIter = PrunedNodes.begin(); pIter != PrunedNodes.end(); pIter++)
	{
		RtreeNode* rNode = a_rtree.m_memory.loadPage(*pIter);
		e0 = rNode->genNodeEntry();
		for (int i = 0; i < dimen; i++)
		{
			pt[i] = e0->m_hc.getUpper()[i];
		}
		mindist = minDist(pt, ORIGNIN, dimen);
		H0.insert(PfltINT(mindist, *pIter));
		delete rNode;
		delete e0;
	}
	PrunedNodes.clear();

	vector<long int> nIDToDelete;
	for (nretIter = NonResultEntry.begin(); nretIter != NonResultEntry.end(); nretIter++)
	{
		if (a_skylines.size() > 0)
		{
			s0Iter = a_skylines.find(nretIter->first - MAXPAGEID);
			if (s0Iter != a_skylines.end())
			{
				nIDToDelete.push_back(nretIter->first);
				continue;
			}
		}

		for (int i = 0; i < dimen; i++)
		{
			pt[i] = (nretIter->second)->m_entry[0]->m_hc.getUpper()[i];
		}
		mindist = minDist(pt, ORIGNIN, dimen);
		H0.insert(PfltINT(mindist, nretIter->first));
	}

	for (pIter = nIDToDelete.begin(); pIter != nIDToDelete.end(); pIter++)
	{
		nretIter = NonResultEntry.find(*pIter);
		delete nretIter->second;
		NonResultEntry.erase(nretIter);
	}

	multimap<float, long int>::iterator fIter;
	while (H0.size() != 0)
	{
		isAnObject = false;

		fIter = H0.begin();
		long int pageid = fIter->second;
		float distTemp = fIter->first;
		H0.erase(fIter);

		VirtualRNode* VirNode = new VirtualRNode;
		RtreeNodeEntry* e0;            // create node entry e0 for node n, so that its MBR can be obtained    	
		if (pageid >= MAXPAGEID)     //current element in H is a data entry node (for distinction,its pageid is equal to m_id+MAXPAGED)
		{
			isAnObject = true;
			nretIter = NonResultEntry.find(pageid);
			if (nretIter == NonResultEntry.end())
			{
				cout << "Error! there is no node " << pageid << " in NonResultEntry!" << endl;
			}
			else
			{
				VirNode->copyData(nretIter->second);
				if (VirNode->m_usedspace > 1)
				{
					cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
				}
				else
					e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node           	
			}
			pageid = pageid - MAXPAGEID;
		}
		else
		{
			RtreeNode* node = a_rtree.m_memory.loadPage(pageid);
			VirNode->copyData(*node);
			VirNode->copyEntries(*node, node->m_usedspace);
			e0 = node->genNodeEntry();     //compute the enclosing MBR for current index node        
			delete node;
		}

		if (isAnObject)
		{
			for (int i = 0; i < dimen; i++)
				pt[i] = VirNode->m_entry[0]->m_hc.getLower()[i] + SIDELEN;
		}
		else
		{
			for (int i = 0; i < dimen; i++)
				pt[i] = e0->m_hc.getUpper()[i];
		}
		bool dominated = IsDominatedBy(dimen, pt, skylines, PG);

		if (!dominated)
		{
			if (VirNode->isLeaf())
			{
				if (VirNode->m_usedspace>1)
				{
					for (int i = 0; i < VirNode->m_usedspace; i++)
					{
						for (int j = 0; j < dimen; j++)
						{
							pt[j] = VirNode->m_entry[i]->m_hc.getUpper()[j];
						}

						bool dominated = IsDominatedBy(dimen, pt, skylines, PG);

						long int NewPageId = VirNode->m_entry[i]->m_id + MAXPAGEID;
						VirtualRNode* node = new VirtualRNode;
						RtreeNodeEntry* newEntry = VirNode->m_entry[i]->clone();
						node->insertEntry(newEntry);
						NonResultEntry.insert(PintVRN(NewPageId, node));
						delete newEntry;

						if (dominated)
							continue;

						mindist = minDist(pt, ORIGNIN, dimen);
						H0.insert(PfltINT(mindist, NewPageId));
					}
				}
				else
				{
					skylines.push_back(pageid);
				}
			}
			else
			{
				for (int i = 0; i < VirNode->m_usedspace; i++)
				{
					for (int j = 0; j < dimen; j++)
					{
						pt[j] = VirNode->m_entry[i]->m_hc.getUpper()[j];
					}

					bool dominated = IsDominatedBy(dimen, pt, skylines, PG);
					if (!dominated)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						H0.insert(PfltINT(mindist, VirNode->m_entry[0]->m_id));
					}
					else
					{
						PrunedNodes.push_back(VirNode->m_entry[i]->m_id);
					}
				}
			}
		}
		else
		{
			if (VirNode->isLeaf())
			{
				if (VirNode->m_usedspace > 1)
				{
					PrunedNodes.push_back(pageid);
				}
			}
			else
			{
				PrunedNodes.push_back(pageid);
			}
		}
		delete VirNode;
		delete e0;
	}
	a_skylines.clear();
	for (vector<long>::iterator iter = skylines.begin(); iter != skylines.end(); iter++)
	{
		a_skylines.insert(*iter);
	}
}

void computeHP(const int dimen, float* PG[], Point& a_pt, const int objcnt, vector<long int>& newAddSL)
{
	newAddSL.clear();
	int docount = 0;
	int decount = 0;

	for (int id = 1; id <= objcnt; id++)
	{
		docount = 0;
		decount = 0;
		float rd[MAXDIMEN];
		int tmpDim = 0;
		for (int d = 0; d < dimen; d++)
		{
			rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;
			if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
		}
		if (tmpDim == dimen)
		{
			continue;
		}


		for (int di = 0; di < dimen; di++)
		{
			if (rd[di] < a_pt.m_coor[di])
			{
				docount++;
			}
			else if (rd[di] > a_pt.m_coor[di])
			{
				decount++;
			}
		}

		if (docount != dimen&& decount != dimen)
		{
			// Reduce d dimension record to d-1 halfspace;
			vector<float> tmpHS;
			float rd_d = rd[dimen - 1];
			float p_d = a_pt[dimen - 1];
			for (int d = 0; d < dimen - 1; d++)
			{
				tmpHS.push_back((rd[d] - rd_d) - (a_pt[d] - p_d));
			}
			tmpHS.push_back(p_d - rd_d);
			tmpHS.push_back(id);            //store the ID of incomparable record in HalfSpace
			HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p	
			RecordIDtoHalfPlaneID.insert(PintInt(id, HalfSpaces.size()));
			newAddSL.push_back(id);
		}
	}
}

void computeHPforkskyband(const int dimen, float* PG[], Point& a_pt, vector<long int>& kskyband, vector<long int>& newAddSL)
{
	newAddSL.clear();
	int docount;
	int decount;

	for (int id = 0; id < kskyband.size(); id++)
	{
		docount = 0;
		decount = 0;

		float rd[MAXDIMEN];
		int tmpDim = 0;
		for (int d = 0; d < dimen; d++)
		{
			rd[d] = (PG[kskyband[id]][d] + PG[kskyband[id]][dimen + d]) / 2.0;
			if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
		}
		if (tmpDim == dimen)
		{
			continue;
		}


		for (int di = 0; di < dimen; di++)
		{
			if (rd[di] < a_pt.m_coor[di])
			{
				docount++;
			}
			else if (rd[di] > a_pt.m_coor[di])
			{
				decount++;
			}
		}

		if (docount != dimen && decount != dimen)
		{
			// Reduce d dimension record to d-1 halfspace;
			vector<float> tmpHS;
			float rd_d = rd[dimen - 1];
			float p_d = a_pt[dimen - 1];
			for (int d = 0; d < dimen - 1; d++)
			{
				tmpHS.push_back((rd[d] - rd_d) - (a_pt[d] - p_d));
			}
			tmpHS.push_back(p_d - rd_d);
			tmpHS.push_back(kskyband[id]);            //store the ID of incomparable record in HalfSpace
			HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p	
			RecordIDtoHalfPlaneID.insert(PintInt(kskyband[id], HalfSpaces.size()));
			newAddSL.push_back(kskyband[id]);
		}
	}
}

void computeHP(const int dimen, float* PG[], Point& a_pt, const int objcnt, vector<long int>& newAddSL, const int NoofConstraints)
{
	newAddSL.clear();
	int docount = 0;
	int decount = 0;

	for (int id = 1; id <= objcnt && newAddSL.size() <= NoofConstraints; id++)
	{
		docount = 0;
		decount = 0;
		float rd[MAXDIMEN];
		int tmpDim = 0;
		for (int d = 0; d < dimen; d++)
		{
			rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;
			if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
		}
		if (tmpDim == dimen)
		{
			continue;
		}


		for (int di = 0; di < dimen; di++)
		{
			if (rd[di] < a_pt.m_coor[di])
			{
				docount++;
			}
			else if (rd[di] > a_pt.m_coor[di])
			{
				decount++;
			}
		}

		if (docount != dimen&& decount != dimen)
		{
			// Reduce d dimension record to d-1 halfspace;
			vector<float> tmpHS;
			float rd_d = rd[dimen - 1];
			float p_d = a_pt[dimen - 1];
			for (int d = 0; d < dimen - 1; d++)
			{
				tmpHS.push_back((rd[d] - rd_d) - (a_pt[d] - p_d));
			}
			tmpHS.push_back(p_d - rd_d);
			tmpHS.push_back(id);            //store the ID of incomparable record in HalfSpace
			HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p	
			RecordIDtoHalfPlaneID.insert(PintInt(id, HalfSpaces.size()));
			newAddSL.push_back(id);
		}
	}
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

void onionlayer(vector<long int>& skyband, float* PG[], const int k, vector<long int >& klayers, int& dim)
{
	char* tmpfile = "tmp_qhull.txt";
	
	unordered_map<int, long int> records;
	unordered_map<long int, int> value_records;

	vector<long int> dataset = skyband;

	for (int ki = 0; ki < k; ki++)
	{
		if (dataset.size() == 0)
		{
			break;
		}
		if (dataset.size() == 1)
		{
			klayers.push_back(dataset[0]);
			break;
		}
		records.clear();
		value_records.clear();
		for (int i = 0; i < dataset.size(); i++)
		{
			records[i] = dataset[i];
			value_records[dataset[i]] = i;
		}
		FILE* fp;
		if ((fp = fopen(tmpfile, "w")) != NULL)
		{
			fprintf(fp, "%d\n", dim);
			fprintf(fp, "%d\n", records.size());

			for (unordered_map<int, long int>::iterator iter = records.begin(); iter != records.end(); iter++)
			{
				for (int di = 0; di < dim; di++)
				{
					fprintf(fp, "%.4f ", (PG[iter->second][di] + PG[iter->second][di + dim]) / 2);
				}
				fprintf(fp, "\n");
			}

			fclose(fp);
		}

		char command[2048] = "D:\\SIGMOD2018\\utk\\utk\\qhull\\qconvex.exe";
		strcat(command, " < ");
		strcat(command, tmpfile);
		strcat(command, " Fx > tmp.ret");
		system(command);

		vector<int> layerrecords;
		fstream fpdata;
		fpdata.open("tmp.ret", ios::in);
		
		int count;
		fpdata >> count;
		int id;
		for (int i = 0; i < count; i++)
		{
			fpdata >> id;
			layerrecords.push_back(records[id]);
		}
		fpdata.close();

		for (int i = 0; i < layerrecords.size(); i++)
		{
			bool isDominated = false;
			for (int j = 0; j < layerrecords.size(); j++)
			{
				if (i != j)
				{
					int dimcount = 0;
					for (int di = 0; di < dim; di++)
					{
						if (PG[layerrecords[i]][di] < PG[layerrecords[j]][di])
						{
							dimcount++;
						}
					}
					if (dimcount == dim)
					{
						isDominated = true;
						break;
					}
				}
			}
			if (isDominated == false)
			{
				klayers.push_back(layerrecords[i]);
				records.erase(value_records[layerrecords[i]]);
			}
		}

		dataset.clear();
		for (auto iter = records.begin(); iter != records.end(); iter++)
		{
			dataset.push_back(iter->second);
		}
	}
}

// comment it as it is an expensive but correct method
/*
void onionlayer(vector<long int>& skyband, float* PG[], const int k, vector<long int >& klayers, int& dim)
{
	unordered_set<long int> skylinerecords;
	for (int i = 0; i < skyband.size(); i++)
	{
		skylinerecords.insert(skyband[i]);
	}

	vector<long int> skylineID;
	unordered_map<int, long int> records;
	multimap<float, int> heap;

	
	for (int i = 0; i < k; i++)
	{
		if (skylinerecords.size() == 0)
		{
			break;
		}

		string tmpfile = "tmp_qhull_";
		string output = "out_";

		tmpfile += to_string(i);
		tmpfile += ".txt";

		skylineID.clear();
		records.clear();
		int recordID = 0;
		for (auto iter = skylinerecords.begin(); iter != skylinerecords.end(); iter++)
		{
			skylineID.push_back(*iter);
			records[++recordID] = *iter;
		}
		
		FILE* fp;
		if ((fp = fopen(tmpfile.c_str(), "w")) != NULL)
		{
			fprintf(fp, "%d\n", dim);
			fprintf(fp, "%d\n", skylineID.size()+1);

			for (int di = 0; di < dim; di++)
			{
				fprintf(fp, "0.0 ");
			}
			fprintf(fp, "\n");

			for (auto iter = skylineID.begin(); iter != skylineID.end(); iter++)
			{
				for (int di = 0; di < dim; di++)
				{
					fprintf(fp, "%.4f ", (PG[*iter][di] + PG[*iter][di + dim]) / 2);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);
		}

		char command[2048] = "D:\\SIGMOD2018\\utk\\utk\\qhull\\qconvex.exe";
		strcat(command, " < ");
		strcat(command, tmpfile.c_str());
		strcat(command, " Fx >");
		output = output + to_string(i) + ".ret";
		strcat(command, output.c_str());
		system(command);

		vector<int> layerrecords;
		fstream fpdata;
		fpdata.open(output.c_str(), ios::in);

		int count;
		fpdata >> count;
		int id;
		for (int i = 0; i < count; i++)
		{
			fpdata >> id;
			if (id != 0)
			{
				klayers.push_back(records[id]);
				skylinerecords.erase(records[id]);
			}
		}
		fpdata.close();
	}
	
}
*/