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
			heap.emplace(score, skyband[i]);
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
//			delete e0;
		}
//		delete node;
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

int countNodeRDominator(const int dim, float* pt, const float radius,
        vector<float>& userpref, vector<long int>& incompSet, float* PG[])
{
    int r_dominate_cnt=0;
    for (const long int&v:incompSet) {
        if(v2_r_dominate_v1_float(pt, PG[v], userpref, g_r_domain_vec, radius)){
            ++r_dominate_cnt;
        }
    }
    return r_dominate_cnt;
}

bool r_dominatedByK(const int dim, float* pt, const float radius,
                        vector<float>& userpref, vector<long int>& incompSet, float* PG[], int k)
{
    int r_dominate_cnt=0;
    for (const long int&v:incompSet) {
        if(v2_r_dominate_v1_float(pt, PG[v], userpref, g_r_domain_vec, radius)){
            ++r_dominate_cnt;
            if(r_dominate_cnt>=k){
                break;
            }
        }
    }
    return r_dominate_cnt>=k;
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
			}
			else
			{
				if (candidateRet.size() < X - k)  // Phase II
				{
					tmpRadius = computeradius(k, dimen, pageID - MAXPAGEID, userpref, incompSet, PG);
					candidateRet.emplace(tmpRadius, pageID-MAXPAGEID);
					if(candidateRet.size()==X-k)
					{
                        raduis = candidateRet.top().first;
                    }
				}
				else if(X<=k)
				{
                    assert(X==k);// if fails there is a problem in data
				    raduis=0;
				    break;
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
            incompSet.push_back(pageID - MAXPAGEID);
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
					Point a_pt(dimen, pt);
					if (raduis == INFINITY)
					{
                        if (!dominatedByK(dimen, pt, incompSet, PG, k))
						{
							heap.emplace(tmpScore, node->m_entry[i]->m_id + MAXPAGEID);
						}
					}
					else
					{
                        if (!r_dominatedByK(dimen, pt, raduis, userpref, incompSet, PG, k))
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
					Point a_pt(dimen, pt);
                    if (raduis == INFINITY)
                    {
                        if (!dominatedByK(dimen, pt, incompSet, PG, k))
                        {
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
                        }
                    }
                    else
                    {
                        if (!r_dominatedByK(dimen, pt, raduis, userpref, incompSet, PG, k))
                        {
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
                        }
                    }
					
				}
			}
		}
	}
	if(interval.size()+candidateRet.size()<X){
	    return INFINITY;
	}else{
        return raduis;
    }
}


class unknown_X_node{
private:
    int dim;
    multiset<float> topK_dominate_radius;// size with k
public:
    int tau;// used in unknown efficient
    const float *data;
    int page_id;
    bool fetched;
    unknown_X_node(int d,  const float *dt, int pid){
        dim=d;
        data=dt;
        page_id=pid;
        fetched=false;
        tau=0;// used in unknown efficient
    }
    float radius_LB(){
//        assert(!topK_dominate_radius.empty());
        return *topK_dominate_radius.begin();
    }

    void init_radius(vector<long int> &fetched_options, float **PointSet, vector<float> &w, int k, float rho=INFINITY){
        for (int i = 0; i < k; ++i) {
            update_radius(PointSet[fetched_options[i]], w);
        }
        for (int j = k; j < fetched_options.size(); ++j) {
            if(this->radius_LB()>=rho){
                break;
            }
            update_radius_erase(PointSet[fetched_options[j]], w);
        }
    }

    inline bool update(int need_to_update) const{// used in unknown efficient
        return this->tau>=need_to_update;
    }

    void update_radius(vector<long int>::iterator begin, vector<long int>::iterator end, float **PointSet, vector<float> &w, float rho=INFINITY){
        for (; begin!=end; ++begin) {
            update_radius_erase(PointSet[*begin], w);
            if(this->radius_LB()>=rho){ // lazy update // used in unknown efficient
                return;
            }
        }
    }

    void update_radius(const float* other_option, vector<float> &w){
        ++tau;// used in unknown efficient
        if(v1_dominate_v2(other_option, data, dim)){
            topK_dominate_radius.insert(INFINITY);
        }else{
            vector<float> tmpHS(data, data+dim);
            for (int i = 0; i < dim; ++i) {
                tmpHS[i]-=other_option[i];
            }
            float tmpDis = computeDis(tmpHS, w);
            topK_dominate_radius.insert(tmpDis);
        }
    }

    void update_radius_erase(const float* other_option, vector<float> &w){
        update_radius(other_option, w);
        topK_dominate_radius.erase(topK_dominate_radius.begin());
    }
};

float computeRho_unknownX_basic(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[])
{
    // phase (1) if get_next_time<=k, from BBS fetch topK
    // phase (2) for each element in BBS, calculate their inflection radius and store them into unordered_set S
    // phase (3) while the the top of S is not an option:
    //          (a) pop and push node A out and into BBS heap
    //          (b) if node A is an option:
    //                  not update S
    //                  marked A as fetched
    //              else
    //                  update S to remove node
    //          (c) update all nodes in S such that not fetched by BBS yet (use a flag FETCH to record), update LB
    vector<pair<int, float>> interval; // return
    multimap<float, unknown_X_node*, greater<float>> heap; //BBS heap, <w\cdot node, node>, if node is an option then node=id+MAXPAGEID
    unordered_set<unknown_X_node*> S;
    vector<long int> incompSet;
    float pt[MAXDIMEN];
    vector<float> ones(dimen, 1);
    vector<float> zeros(dimen, 0);
    unknown_X_node *zeros_node=new unknown_X_node(dimen, zeros.data(), -1);
    unknown_X_node *rt=new unknown_X_node(dimen, ones.data(), a_rtree.m_memory.m_rootPageID);
    heap.emplace(INFINITY, rt);
    while (!heap.empty() && interval.size()<X)
    {
        unknown_X_node *popped_node=heap.begin()->second;
        popped_node->fetched=true;
        int pageID = popped_node->page_id;
        heap.erase(heap.begin());
        if (pageID > MAXPAGEID)  // option processing
        {
            if (interval.size() < k)  // Phase (1)
            {
                interval.emplace_back(pageID - MAXPAGEID, 0);
                if(interval.size()==k) // Phase (2)
                {
                    //init S
                    for (pair<int, float> &option:interval) {
                        incompSet.push_back(option.first);
                    }
                    for(pair<const float, unknown_X_node *> &ele:heap)
                    {
                        unknown_X_node *node=ele.second;
                        node->init_radius(incompSet, PG, userpref, k);
                        S.insert(node);
                    }
                }
            }
            else
            {
                if (interval.size() < X )  // should get_next X times
                {
                    //update all element in S
                    assert(!S.empty());
                    for (unknown_X_node* node:S)
                    {
                        if (!node->fetched)
                        {
                            node->update_radius_erase(popped_node->data, userpref);
                        }
                    }
                }
                else if(X<=k)
                {
                    assert(X==k);// if fails there is a problem in data
                    break;
                }
                else   // interval.size() == X, should begin to return
                {
                    break;
                }
                incompSet.push_back(pageID - MAXPAGEID);
            }
        }
        else // internal and leaf nodes processing
        {
            S.erase(popped_node);
            delete(popped_node);
            RtreeNode* node = ramTree[pageID];
            if (node->isLeaf())
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    float tmpScore = 0;
                    for (int j = 0; j < dimen; j++)
                    {
                        pt[j] = node->m_entry[i]->m_hc.getCenter()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    unknown_X_node *tmp_node=new unknown_X_node(dimen, PG[node->m_entry[i]->m_id], node->m_entry[i]->m_id + MAXPAGEID);
                    heap.emplace(tmpScore, tmp_node);
                    if(interval.size()>=k) {
                        tmp_node->init_radius(incompSet, PG, userpref, k);
                        S.insert(tmp_node);
                    }
                }
            }
            else
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    float tmpScore = 0;
                    for (int j = 0; j < dimen; j++)
                    {
                        pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    const float *ptr=node->m_entry[i]->m_hc.getUpper().m_coor;
                    unknown_X_node *tmp_node=new unknown_X_node(dimen, ptr, node->m_entry[i]->m_id);
                    heap.emplace(tmpScore, tmp_node);
                    if(interval.size()>=k) {
                        tmp_node->init_radius(incompSet, PG, userpref, k);
                        S.insert(tmp_node);
                    }
                }
            }

        }
        if(interval.size()>=k) {
            pair<float, unknown_X_node *> LB(INFINITY, zeros_node);
            for (unknown_X_node *node:S) {
                if (node->radius_LB() < LB.first) {
                    LB.first = node->radius_LB();
                    LB.second = node;
                }
            }
            while (LB.second->page_id > MAXPAGEID && LB.second->fetched) {
                // if LB is an option and is updated, add it to result list "interval"
                interval.emplace_back(LB.second->page_id - MAXPAGEID, LB.second->radius_LB());
                S.erase(LB.second);
                delete(LB.second);
                if (interval.size() == X) {
                    break;
                }
                LB.first = INFINITY;
                LB.second = zeros_node;
                for (unknown_X_node *node:S) {
                    if (node->radius_LB() < LB.first) {
                        LB.first = node->radius_LB();
                        LB.second = node;
                    }
                }
            }
        }
    }

    for(unknown_X_node *node:S){
        delete (node);
    }
    delete(zeros_node);
    if(X<=k){
        return 0;
    }else if(interval.size()<X){
        return INFINITY;
    }else{
        return interval.back().second;
    }
}

float computeRho_unknownX_efficient(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[])
{
    // phase (1) if get_next_time<=k, from BBS fetch topK
    // phase (2) for each element in BBS, calculate their inflection radius and store them into set S,  min_heap Q


    typedef float RADIUS;  // is the inflection radius
    typedef int PAGE_ID;  // page id in rtree
    typedef int PG_IDX;  // idx in PG
    typedef long PG_IDX_L;  // idx in PG
    typedef float DOT;   // is the result of dot product with w
    typedef unknown_X_node* NODE_PTR;

    vector<pair<PG_IDX, RADIUS>> interval; // return
    multimap<DOT, NODE_PTR, greater<float>> heap; //BBS max_heap, <w\cdot node, node>, if node is an option then node=id+MAXPAGEID
    unordered_set<NODE_PTR> S; // reflect which nodes and options in BBS
    multimap<RADIUS, NODE_PTR, less<float>> Q; // min_heap based on inflection radius, lazy update for S
    multimap<RADIUS, NODE_PTR, less<float>> C; // min_heap based on inflection radius, candidate list
    vector<PG_IDX_L> incompSet;
    float pt[MAXDIMEN];
    float radius = INFINITY;
    vector<float> ones(dimen, 1);
    NODE_PTR rt=new unknown_X_node(dimen, ones.data(), a_rtree.m_memory.m_rootPageID);
    heap.emplace(INFINITY, rt);

    while (!heap.empty() && interval.size()<X)
    {
        NODE_PTR popped_node=heap.begin()->second;
        PAGE_ID pageID = popped_node->page_id;
        heap.erase(heap.begin());
        if(interval.size()>=k){
            S.erase(popped_node);
        }

        if (pageID > MAXPAGEID)  // option processing
        {
            if (interval.size() < k)  // Phase (1)
            {
                interval.emplace_back(pageID - MAXPAGEID, 0);
                if(interval.size()==k) // Phase (2)
                {
                    //init S, init Q
                    for (pair<PG_IDX, RADIUS> &option:interval) {
                        incompSet.push_back(option.first);
                    }
                    for(pair<const DOT, NODE_PTR> &ele:heap)
                    {
                        NODE_PTR node=ele.second;
                        node->init_radius(incompSet, PG, userpref, k);
                        S.insert(node);
                        Q.emplace(node->radius_LB(), node);
                    }
                }
            }
            else
            {
                if (interval.size() < X )  // should get_next X times
                {
                    popped_node->update_radius(incompSet.begin()+popped_node->tau, incompSet.end(), PG, userpref);
                    C.emplace(popped_node->radius_LB(), popped_node);
                    incompSet.push_back(pageID - MAXPAGEID);
                }
                else if(X<=k)
                {
                    assert(X==k);// if fails there is a problem in data
                    radius=0;
                    break;
                }
                else   // interval.size() == X, should begin to return
                {
                    radius=interval.back().second;
                    break;
                }
            }
        }
        else // internal and leaf nodes processing
        {
            float possible_next_radius=INFINITY;
            if(!C.empty()){
                possible_next_radius=C.begin()->first;
            }
            RtreeNode* node = ramTree[pageID];
            if (node->isLeaf())
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    DOT tmpScore = 0;
                    for (int j = 0; j < dimen; j++)
                    {
                        pt[j] = node->m_entry[i]->m_hc.getCenter()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    unknown_X_node *tmp_node=new unknown_X_node(dimen, PG[node->m_entry[i]->m_id], node->m_entry[i]->m_id + MAXPAGEID);
                    heap.emplace(tmpScore, tmp_node);
                    if(interval.size()>=k) {
                        S.insert(tmp_node);
                        tmp_node->init_radius(incompSet, PG, userpref, k, possible_next_radius);
                        Q.emplace(tmp_node->radius_LB(), tmp_node);
                    }
                }
            }
            else
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    DOT tmpScore = 0;
                    for (int j = 0; j < dimen; j++)
                    {
                        pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    const float *ptr=node->m_entry[i]->m_hc.getUpper().m_coor;
                    unknown_X_node *tmp_node=new unknown_X_node(dimen, ptr, node->m_entry[i]->m_id);
                    heap.emplace(tmpScore, tmp_node);
                    if(interval.size()>=k) {
                        S.insert(tmp_node);
                        tmp_node->init_radius(incompSet, PG, userpref, k, possible_next_radius);
                        Q.emplace(tmp_node->radius_LB(), tmp_node);
                    }
                }
            }

        }
        if(interval.size()>=k && !C.empty()) {
            _Rb_tree_iterator<pair<const RADIUS, NODE_PTR>> possible_next = C.begin();

            // a lazy update with S
            while (!Q.empty() && S.find(Q.begin()->second) == S.end()) {
                if (Q.begin()->second->page_id <= MAXPAGEID) { // directly delete an inter node
                    delete (Q.begin()->second);
                }
                Q.erase(Q.begin());
            }

            // make sure the top of Q is updated or its inflection radius is larger than possible_next.inflection radius
            while (!Q.empty() &&
                   !(Q.begin()->second->update(incompSet.size()) || Q.begin()->first >= possible_next->first)) {
                NODE_PTR qnode = Q.begin()->second;

                // a lazy update with tau
                qnode->update_radius(incompSet.begin() + qnode->tau, incompSet.end(), PG, userpref,
                                     possible_next->first);
//                qnode->update_radius(incompSet.begin() + qnode->tau, incompSet.end(), PG, userpref,
//                                    INFINITY);

                Q.erase(Q.begin());
                Q.emplace(qnode->radius_LB(), qnode);
                while (!Q.empty() && S.find(Q.begin()->second) == S.end()) { // a lazy update with S
                    if (Q.begin()->second->page_id <= MAXPAGEID) { // directly delete an inter node
                        delete (Q.begin()->second);
                    }
                    Q.erase(Q.begin());
                }
            }


            // if Q empty, add element in C into interval one by one and return this function
            if (Q.empty()) {
                // begin getting ready to return
                while (interval.size() < X && !C.empty()) {
                    interval.emplace_back(C.begin()->second->page_id - MAXPAGEID, C.begin()->first);
                    delete (C.begin()->second);
                    C.erase(C.begin());
                }
                if (interval.size() >= X) {
                    radius = interval.back().second;
                }
                break;
            }
            // if Q not empty, then check whether possible_next.if_r is lower than Q.top.if_r
            // if possible_next.if_r is lower than Q.top.if_r:
            //     add possible_next into result list "interval"
            // else:
            //     continue fetch nodes or options with BBS
            if (possible_next->first <= Q.begin()->first) {
                interval.emplace_back(C.begin()->second->page_id - MAXPAGEID, C.begin()->first);
                delete (C.begin()->second);
                C.erase(C.begin());
                // if still can fetch from candidate list C
                while (!C.empty() && interval.size() < X) {
                    possible_next = C.begin();
                    if (possible_next->first <= Q.begin()->first) {
                        interval.emplace_back(C.begin()->second->page_id - MAXPAGEID, C.begin()->first);
                        delete (C.begin()->second);
                        C.erase(C.begin());
                    } else {
                        break;
                    }
                }
            }
        }
    }

    for(NODE_PTR node:S){
        delete (node);
    }
    for(pair<const RADIUS , NODE_PTR> node_iter:C){
        delete (node_iter.second);
    }
    if(interval.size()==X && X>=k){
        radius=interval.back().second;
    }
    return radius;
}


