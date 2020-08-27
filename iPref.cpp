#include "iPref.h"
#include "qp_solver.h"
#include "utk_math_lib.h"
vector<long int> computeTopK(const int dim, float* PG[], vector<long int> skyband, vector<float>& weight, int k)
{
	multimap<float, long int> heap;

	for (int i = 0; i < skyband.size(); i++)
	{
		float score = 0;
		for (int d = 0; d < dim; d++)
		{
			score += (PG[skyband[i]][d] + PG[skyband[i]][d + dim]) / 2*weight[d];
		}

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
	int dimen = weight.size();
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

	for (int i = 0; i < dimen; i++)
	{
		spi += piv[i] * weight[i];
		spj += pjv[i] * weight[i];

		if (piv[i] <= pjv[i])
		{
			cneg++;
		}
		else
		{
			cpos++;
		}
	}

	if (spj >= spi && cneg != 0 && cpos != 0)
		return true;
	else
		return false;
}


vector<float> computePairHP(const int dimen, float* PG[], long int pi, long int pj)
{
	vector<float> retHS(dimen);
    for (int i = 0; i < dimen; ++i) {
        retHS[i]=PG[pi][i]-PG[pj][i];
    }
//	float piv[MAXDIMEN];
//	float pjv[MAXDIMEN];
//	for (int d = 0; d < dimen; d++)
//	{
//		piv[d] = (PG[pi][d] + PG[pi][dimen + d]) / 2.0;
//		pjv[d] = (PG[pj][d] + PG[pj][dimen + d]) / 2.0;
//	}
//
//	float pj_d = pjv[dimen - 1];
//	float pi_d = piv[dimen - 1];
//	for (int d = 0; d < dimen - 1; d++)
//	{
//		retHS.push_back((pjv[d] - pj_d) - (piv[d] - pi_d));
//	}
//	retHS.push_back(pi_d - pj_d);
	return retHS;
}

extern qp_solver qp;

float computeDis(vector<float> tmpHS, vector<float> userpref)
{
    float ret=qp.update_w_h_solve(userpref, tmpHS);
    return ret;
}

bool sortbysec(const pair<long int, float> &a, const pair<long int, float> &b)
{
	return (a.second < b.second);
}

float computeradius(const int k, const int dim, long int pi, vector<float>& userpref, vector<long int>& incompSet, float* PG[])
{
	vector<float> radiusSKI;
	vector<float> tmpHS(dim);
	float tmpDis;

	for (int ri = 0; ri < incompSet.size(); ri++)
	{
		if (IsPjdominatePi(dim, PG, pi, incompSet[ri]))  // for these dominators, we directly input as FLTMAX
		{
			radiusSKI.push_back(INFINITY);
		}
		else
		{
//			tmpHS.clear();
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

bool v2_r_dominate_v1_float(float* &v1, float* &v2, const vector<float> &w, const vector<vector<c_float>> &r_domain_vec, const float &rho) {
    for (const vector<c_float> &v:r_domain_vec) {
//        vector<float> tmp_w = w + rho * v;
//        if (v1 * tmp_w < v2 * tmp_w) {
//            return false;
//        }
        vector<float> tmp_w(w.size());
        for (int i = 0; i < tmp_w.size(); ++i) {
            tmp_w[i]=w[i]+rho*v[i];
        }
        float s1=0, s2=0;
        for (int i = 0; i < tmp_w.size(); ++i) {
            s1+=tmp_w[i]*v1[i];
        }
        for (int i = 0; i < tmp_w.size(); ++i) {
            s2+=tmp_w[i]*v2[i];
        }
        if(s1>s2){
            return false;
        }
    }
    return true;
}

int countNodeRDominator(const int dim, float* pt, const float radius, vector<float>& userpref, vector<long int>& incompSet, float* PG[]) // I am not sure it is correct if we only consider incompSet.
{
    int r_dominate_cnt=0;
    for (const long int&v:incompSet) {
        if(v2_r_dominate_v1_float(pt, PG[v], userpref, g_r_domain_vec, radius)){
            ++r_dominate_cnt;
//            if(r_dominate_cnt>=k){ // TODO use it to prune, pass k to this function
//                break;
//            }
        }
    }
    return r_dominate_cnt;
}


float computeRho(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[])
{
	float raduis = INFINITY;
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

	heap.emplace(INFINITY, a_rtree.m_memory.m_rootPageID);

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
					candidateRet.emplace(tmpRadius, pageID-MAXPAGEID);
					incompSet.push_back(pageID - MAXPAGEID);
				}
				else   // Phase III
				{
					tmpRadius = computeradius(k, dimen, pageID - MAXPAGEID, userpref, incompSet, PG);
					if (tmpRadius < candidateRet.top().first)
					{
						candidateRet.emplace(tmpRadius, pageID - MAXPAGEID);
						candidateRet.pop();
                        candidateOpt = candidateRet.top();
                        raduis = candidateOpt.first;
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
					if (raduis == INFINITY)
					{
//						if (countDominator(a_rtree, PG, a_pt) < k)
//						{
							heap.emplace(tmpScore, node->m_entry[i]->m_id + MAXPAGEID);
//						}
					}
					else
					{
                        if (countNodeRDominator(dimen, pt, raduis, userpref, incompSet, PG)<k)
						{
                            heap.emplace(tmpScore, node->m_entry[i]->m_id + MAXPAGEID);
                        }
						
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
                    if (raduis == INFINITY)
                    {
//                        if (countDominator(a_rtree, PG, a_pt) < k)
//                        {
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
//                        }
                    }
                    else
                    {
                        if (countNodeRDominator(dimen, pt, raduis, userpref, incompSet, PG)<k)
                        {
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
                        }
                    }
					
				}
			}
		}
	}
	return raduis;
}
