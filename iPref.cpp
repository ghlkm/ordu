#include "iPref.h"

vector<long int> computeTopK(const int dim, float* PG[], vector<long int> skyband, vector<float>& weight, int k)
{
	multimap<float, long int> heap;

	for (int i = 0; i < skyband.size(); i++)
	{
		float score = 0;
		float l1_sumWeight = 0;
		for (int d = 0; d < dim-1; d++)
		{
			score += (PG[skyband[i]][d] + PG[skyband[i]][d + dim]) / 2*weight[d];
            l1_sumWeight += weight[d];
		}
        score += (PG[skyband[i]][dim - 1] + PG[skyband[i]][dim - 1 + dim]) / 2 * (1 - l1_sumWeight);

		if (heap.size() < k)
		{
			heap.insert(multimap<float, long int>::value_type(score, skyband[i]));
			//could use emplace， e.g. “heap.emplace(score, skyband[i])”
		}
		else if (heap.size() == k && heap.begin()->first < score)
		{
			heap.erase(heap.begin());
			heap.emplace(score, skyband[i]);
		}
	}

	vector<long int> topkRet;
	for (auto heapIter = heap.rbegin(); heapIter != heap.rend(); ++heapIter)
	{
		topkRet.push_back(heapIter->second);
	}
	return topkRet;
}

bool IsPjdominatePi(const int dimen, float* PG[], long int pi, long int pj)
{
	float pid;
	float pjd;
	int count = 0;
	for (int d = 0; d < dimen; d++)
	{
		pid = (PG[pi][d] + PG[pi][dimen + d]) / 2.0;
		pjd = (PG[pj][d] + PG[pj][dimen + d]) / 2.0;
		if (pjd >= pid)
			count++;
	}


	if (count == dimen)
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool incomparableset(float* PG[], long int pi, long int pj, vector<float>& weight)
{
	int dimen = weight.size() + 1;
	float l1_wd = 0;
	float spi = 0, spj = 0;
	int cpos = 0;
	int cneg = 0;

	float piv[MAXDIMEN];
	float pjv[MAXDIMEN];
	for (int d = 0; d < dimen; d++)
	{
		piv[d] = (PG[pi][d] + PG[pi][dimen + d]) / 2.0;
		pjv[d] = (PG[pj][d] + PG[pj][dimen + d]) / 2.0;
	}

	for (int i = 0; i < dimen - 1; i++)
	{
		spi += piv[i] * weight[i];
		spj += pjv[i] * weight[i];
        l1_wd += weight[i];

		if (piv[i] <= pjv[i])
		{
			cneg++;
		}
		else
		{
			cpos++;
		}
	}

	spi += (1 - l1_wd) * piv[dimen - 1];
	spj += (1 - l1_wd) * pjv[dimen - 1];

	if (piv[dimen - 1] <= pjv[dimen - 1])
	{
		cneg++;
	}
	else
	{
		cpos++;
	}

	if (spj >= spi && cneg != 0 && cpos != 0)
		return true;
	else
		return false;
}


vector<float> computePairHP(const int dimen, float* PG[], long int pi, long int pj)
{
	vector<float> retHS;
	float piv[MAXDIMEN];
	float pjv[MAXDIMEN];
	for (int d = 0; d < dimen; d++)
	{
		piv[d] = (PG[pi][d] + PG[pi][dimen + d]) / 2.0;
		pjv[d] = (PG[pj][d] + PG[pj][dimen + d]) / 2.0;
	}

	float pj_d = pjv[dimen - 1];
	float pi_d = piv[dimen - 1];
	for (int d = 0; d < dimen - 1; d++)
	{
		retHS.push_back((pjv[d] - pj_d) - (piv[d] - pi_d));
	}
	retHS.push_back(pi_d - pj_d);
	return retHS;
}


float computeDis(vector<float> tmpHS, vector<float> userpref)
{
	//TODO replace all these function to find correct radius of pi
	float normVec = 0;
	float dotRet = 0;
	for (int i = 0; i < tmpHS.size() - 1; i++)
	{
		dotRet += userpref[i] * tmpHS[i];
		normVec += tmpHS[i] * tmpHS[i];
	}
//  wrong
//	return abs(dotRet - tmpHS[tmpHS.size() - 1]) / normVec;
//  right for min distance of point to hyperplane
    return abs(dotRet - tmpHS[tmpHS.size() - 1]) / sqrt(normVec);
}

bool sortbysec(const pair<long int, float> &a, const pair<long int, float> &b)
{
	return (a.second < b.second);
}

float computeradius(const int k, const int dim, long int pi, vector<float>& userpref, vector<long int>& incompSet, float* PG[])
{
	vector<float> radiusSKI;
	vector<float> tmpHS;
	float tmpDis;

	for (int ri = 0; ri < incompSet.size(); ri++)
	{
		if (IsPjdominatePi(dim, PG, pi, incompSet[ri]))  // for these dominators, we directly input as FLTMAX
		{
			radiusSKI.push_back(FLT_MAX);
		}
		else
		{
			tmpHS.clear();
			tmpHS = computePairHP(dim, PG, pi, incompSet[ri]);
			tmpDis = computeDis(tmpHS, userpref);
			radiusSKI.push_back(tmpDis);
		}
	}
	sort(radiusSKI.begin(), radiusSKI.end());
	return radiusSKI[radiusSKI.size() - k];
}

int countDominator(Rtree& a_rtree, float* PG[], Point& a_pt)
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

		//node = a_rtree.m_memory.loadPage(pageid);
		node = ramTree[pageid];
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
			delete e0;
		}
		delete node;
	}
	return NoOfDominators;
}

int countNodeDominator(const int dim, float* pt, vector<long int>& incompSet, float* PG[]) // I am not sure it is correct if we only consider incompSet.
{
	int cnt = 0;
	for (int i = 0; i < incompSet.size(); i++)
	{
		bool isDominator = true;
		for (int di = 0; di < dim; di++)
		{
			if ((PG[incompSet[i]][di] + PG[incompSet[i]][dim + di]) / 2.0 < pt[di])
			{
				isDominator = false;
				break;
			}
		}
		if (isDominator == true)
		{
			cnt++;
		}
	}
	return cnt;
}

int countNodeRDominator(const int dim, float* pt, const float radius, vector<float>& userpref, vector<long int>& incompSet, float* PG[]) // I am not sure it is correct if we only consider incompSet.
{
	//TODO
	// we have troubles here, as we do not know how to compute the MBR of userpref with radius.
	//int cnt = 0;
	//float tmp = 0;
	//float ptScore = 0, tmpScore = 0;
	//bool isDominator = true;

	//for (int i = 0; i < incompSet.size(); i++)
	//{
	//	isDominator = true;
	//	ptScore = 0; // the score of pt
	//	tmpScore = 0; // the score of incomparable records
	//	for (int di = 0; di < dim; di++)
	//	{
	//		tmp = PG[incompSet[i]][di] + PG[incompSet[i]][dim + di] / 2.0;
	//		pt
	//	}
	//	if (isDominator == true)
	//	{
	//		cnt++;
	//	}
	//}
	//return cnt;
	return 0;
}


float computeRho(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[])
{
	float raduis = FLT_MAX;
	vector<pair<long int, float>> interval; // top-k result of w: T
	vector<long int> incompSet;
	pair<float, int> candidateOpt;

	RtreeNode* node;
	priority_queue<pair<float, int>> heap;
	priority_queue<pair<float, int>> candidateRet;

	float pt[MAXDIMEN];
	float maxscore;
	int pageID;
	float tmpScore; // compute the score of option or mbr e w.r.t. userpref
	float tmpDis; // compute the inflection radius of point pt
	float tmpRadius; // the inflection radius of point Pi.

	heap.push(make_pair(INFINITY, a_rtree.m_memory.m_rootPageID));

	while (!heap.empty())
	{
		tmpScore = heap.top().first;
		pageID = heap.top().second;
		heap.pop();

		if (pageID > MAXPAGEID)  // option processing
		{
			if (interval.size() < k)  // Phase I
			{
				interval.emplace_back(pageID - MAXPAGEID, 0);
				incompSet.push_back(pageID - MAXPAGEID);
			}
			else
			{
				if (candidateRet.size() < X - k)  // Phase II
				{
					tmpRadius = computeradius(k, dimen, pageID - MAXPAGEID, userpref, incompSet, PG);
					candidateRet.push(make_pair(tmpRadius, pageID-MAXPAGEID));
					incompSet.push_back(pageID - MAXPAGEID);
				}
				else   // Phase III
				{
					tmpRadius = computeradius(k, dimen, pageID - MAXPAGEID, userpref, incompSet, PG);
					if (tmpRadius < candidateRet.top().first)
					{
						candidateRet.push(make_pair(tmpRadius, pageID - MAXPAGEID));
						candidateOpt = candidateRet.top();
						candidateRet.pop();
						raduis = (raduis < candidateOpt.first) ? candidateOpt.first : raduis;
					}
				}
			}
		}  
		else // internal and leaf nodes processing
		{
			node = ramTree[pageID];
			if (node->isLeaf())
			{	
				for (int i = 0; i < node->m_usedspace; i++)
				{
					tmpScore = 0;
					for (int j = 0; j < dimen; j++)
					{
						pt[j] = node->m_entry[i]->m_hc.getCenter()[j];
						tmpScore += pt[j] * userpref[j];
					}
					//if (countNodeDominator(dimen, pt, incompSet, PG) < k)
					Point a_pt(dimen, pt);
					if (raduis == FLT_MAX)
					{
						if (countDominator(a_rtree, PG, a_pt) < k)
						{
							heap.push(make_pair(tmpScore, node->m_entry[i]->m_id + MAXPAGEID));
						}
					}
					else
					{
						//TODO
						//if (countNodeRDominator()) 
						//{
						//
						//}
						
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					tmpScore = 0;
					for (int j = 0; j < dimen; j++)
					{
						pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
						tmpScore += pt[j] * userpref[j];
					}
					//if (countNodeDominator(dimen, pt, incompSet, PG) < k)
					Point a_pt(dimen, pt);
					if (countDominator(a_rtree, PG, a_pt) < k)
					{
						heap.push(make_pair(tmpScore, node->m_entry[i]->m_id));
					}
					
				}
			}
		}
	}
	return raduis;
}

















































// comment Keming's code, please do not insert function into main.cpp
// We use only left main function in it as it is more clear for others to read and compile the main structure of the code.
// PLEASE do not use "auto"

//template<typename V>
//inline V proj(const V& u, const V &v)
//{
//	/*
//	* <u, v> / <u, u> * u
//	*/
//	float dot = 0;
//	assert(u.size() == v.size());
//	for (int i = 0; i < u.size(); ++i) {
//		dot += u[i] * v[i];
//	}
//	float u_mod2 = 0;
//	for (int i = 0; i < u.size(); ++i) {
//		u_mod2 += u[i] * u[i];
//	}
//	float length = dot / u_mod2;
//	V ret(u);
//	for (int i = 0; i < u.size(); ++i){
//		ret[i] *= length;
//	}
//	return ret;
//}
//
//template<typename V>
//float domin_r_ij2(const V &w, const V &h_ij) {
//	//TODO this fun is in utk_math_lib, merge utk_math_lib later
//	//TODO when vertical point is out of domain, use QP solver
//	/*
//	*     V ones(h_ij.size(), 1);
//	*     V ones_n=ones-proj(h_ij, ones);
//	*     V d=proj(h_ij, WA)+proj(ones_n, WA);
//	*     return sqrt(d*d);
//	*/
//	assert(h_ij.size() >= 2); //dim should be no less than 2
//	V ones(h_ij.size(), 1.0);
//	V ones_proj_on_hij = proj(h_ij, ones); //proj(h_ij, ones)
//	for (int i = 0; i < ones.size(); ++i) {
//		ones[i] -= ones_proj_on_hij[i];   // V ones_n=ones-proj(h_ij, ones);
//	}
//	V point_of_hij(h_ij.size(), 0.0);// point A
//	float sum = -h_ij[0] + h_ij[1];
//	if (sum != 0){
//		point_of_hij[0] = h_ij[1] / sum;
//		point_of_hij[1] = -h_ij[0] / sum; // normalize A such that \Sigma A[i]=1
//	}
//	else{
//		point_of_hij[0] = 0.5;
//		point_of_hij[1] = 0.5;
//	}
//	// point A init finish
//	for (int k = 0; k < w.size(); ++k) {
//		point_of_hij[k] -= w[k]; //A-->WA
//	}
//	V A_proj_on_hij = proj(h_ij, point_of_hij); //proj(h_ij, WA)
//	V A_proj_on_onesn = proj(ones, point_of_hij); //proj(ones_n, WA)
//	for (int j = 0; j < A_proj_on_hij.size(); ++j) {
//		A_proj_on_hij[j] += A_proj_on_onesn[j]; //d=proj(h_ij, WA)+proj(ones_n, WA)
//	}
//	float WA_proj_on_hij_ones_mod2 = 0;
//	for (int j = 0; j < A_proj_on_hij.size(); ++j) {
//		WA_proj_on_hij_ones_mod2 += A_proj_on_hij[j] * A_proj_on_hij[j]; // d*d
//	}
//	return sqrt(WA_proj_on_hij_ones_mod2); //return sqrt(d*d)
//}