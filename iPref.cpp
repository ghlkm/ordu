#include "iPref.h"
#include "qp_solver.h"
#include "utk_math_lib.h"
#include "qhull_user.h"
#include "lp_user.h"


inline double max(double a, float b){
    return a>b?a:b;
}

inline double max(float b, double a){
    return a>b?a:b;
}

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

float computeradius(const int k, const int dim, long int pi, vector<float>& userpref, vector<long int>& incompSet, float* PG[], float rho)
{
	multiset<float> radiusSKI;
	vector<float> tmpHS(dim);
	float tmpDis;
	for (int ri = 0; ri < incompSet.size(); ri++)
	{
		if (IsPjdominatePi(dim, PG, pi, incompSet[ri]))  // for these dominators, we directly input as FLTMAX
		{
			radiusSKI.insert(INFINITY);
		}
		else
		{
			tmpHS = computePairHP(dim, PG, pi, incompSet[ri]);
			tmpDis = computeDis(tmpHS, userpref);
			radiusSKI.insert(tmpDis);
		}
		if(radiusSKI.size()>k){
		    radiusSKI.erase(radiusSKI.begin());
		    if(*radiusSKI.begin()>=rho){
		        break;
		    }
		}
	}
	return *radiusSKI.begin();
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

float computeRho(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[],
        vector<pair<long int, float>> &interval, float radius)
{
//	float raduis = INFINITY;
//	vector<pair<long int, float>> interval; // top-k result of w: T
	vector<long int> incompSet;
	pair<float, int> candidateOpt;

	RtreeNode* node;
	multimap<float, int, greater<float>> heap;
    multimap<float, int, greater<float>> candidateRet;

	float pt[MAXDIMEN];
	float maxscore;
	int pageID;
	float tmpScore; // compute the score of option or mbr e w.r.t. userpref
	float tmpDis; // compute the inflection radius of point pt
	float tmpRadius; // the inflection radius of point Pi.

	heap.emplace(INFINITY, a_rtree.m_memory.m_rootPageID);

	// begin: use for DEBUG
//    bool flag=false;
    // end: use for DEBUG

    while (!heap.empty())
	{
		tmpScore = heap.begin()->first;
		pageID = heap.begin()->second;
		heap.erase(heap.begin());

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
					if(tmpRadius!=INFINITY){
                        candidateRet.emplace(tmpRadius, pageID-MAXPAGEID);
                        if(candidateRet.size()==X-k)
                        {
                            radius = candidateRet.begin()->first;
                        }
					}
					// begin: use for DEBUG
//					else{
//					    if(!flag){
//					        cout<<interval.back().second<<endl;
//					        flag=true;
//					    }
//					}
                    // end: use for DEBUG
                }
				else if(X<=k)
				{
                    assert(X==k);// if fails there is a problem in data
                    radius=0;
				    break;
				}
				else   // Phase III
				{
					tmpRadius = computeradius(k, dimen, pageID - MAXPAGEID, userpref, incompSet, PG, candidateRet.begin()->first);
					if (tmpRadius < candidateRet.begin()->first)
					{
						candidateRet.emplace(tmpRadius, pageID - MAXPAGEID);
						candidateRet.erase(candidateRet.begin());
                        candidateOpt = *candidateRet.begin();
                        radius = candidateOpt.first;
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
					if (radius == INFINITY)
					{
                        if (!dominatedByK(dimen, pt, incompSet, PG, k))
						{
							heap.emplace(tmpScore, node->m_entry[i]->m_id + MAXPAGEID);
						}
					}
					else
					{
                        if (!r_dominatedByK(dimen, pt, radius, userpref, incompSet, PG, k))
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
                    if (radius == INFINITY)
                    {
                        if (!dominatedByK(dimen, pt, incompSet, PG, k))
                        {
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
                        }
                    }
                    else
                    {
                        if (!r_dominatedByK(dimen, pt, radius, userpref, incompSet, PG, k))
                        {
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
                        }
                    }
					
				}
			}
		}
	}

    for (auto riter=candidateRet.rbegin();riter!=candidateRet.rend();++riter) {
        interval.emplace_back(riter->second, riter->first);
    }
//    for (pair<long int, float> &p:interval) {
//        cout<<p.second<<endl;
//    }
	if(interval.size()<X){
	    return INFINITY;
	}else{
        return radius;
    }
}



unknown_X_node::unknown_X_node(int d,  const float *dt, int pid){
        dim=d;
        data=dt;
        page_id=pid;
        fetched=false;
        tau=0;// used in unknown efficient
}

float unknown_X_node::radius_LB(){
//        assert(!topK_dominate_radius.empty());
    return *topK_dominate_radius.begin();
}

void unknown_X_node::init_radius(vector<long int> &fetched_options, float **PointSet, vector<float> &w, int k, float rho){
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

inline bool unknown_X_node::update(int need_to_update) const{// used in unknown efficient
    return this->tau>=need_to_update;
}

//inline bool update() const{
//    return
//}


void unknown_X_node::update_radius(vector<long int>::iterator begin, vector<long int>::iterator end, float **PointSet, vector<float> &w, float rho){
    for (; begin!=end; ++begin) {
        update_radius_erase(PointSet[*begin], w);
        if(this->radius_LB()>=rho){ // lazy update // used in unknown efficient
            return;
        }
    }
}

void unknown_X_node::update_radius(const float* other_option, vector<float> &w){
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

void unknown_X_node::update_radius_erase(const float* other_option, vector<float> &w){
    update_radius(other_option, w);
    topK_dominate_radius.erase(topK_dominate_radius.begin());
}

float computeRho_unknownX_basic(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float** PG)
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



unknown_x_efficient::unknown_x_efficient(const int dim, int K, vector<float> &userPref, Rtree &aRtree, float **pg
) : a_rtree(aRtree) {
    this->dimen=dim;
    this->k=K;
    this->userpref=userPref;
    this->PG=pg;
    vector<float> ones(dimen, 1);
    NODE_PTR rt=new unknown_X_node(dimen, ones.data(), a_rtree.m_memory.m_rootPageID);
    heap.emplace(INFINITY, rt);
}

unknown_x_efficient::~unknown_x_efficient(){
    for(NODE_PTR node:S){
        delete (node);
    }
    for(pair<const RADIUS , NODE_PTR> node_iter:C){
        delete (node_iter.second);
    }
}

pair<int, float> unknown_x_efficient::get_next() {
    while (!heap.empty() || !C.empty()) {
        if (interval.size() >= k && !C.empty()) {
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
            // if Q not empty, then check whether possible_next.if_r is lower than Q.top.if_r
            // if possible_next.if_r is lower than Q.top.if_r:
            //     add possible_next into result list "interval"
            // else:
            //     continue fetch nodes or options with BBS
            if (Q.empty() || possible_next->first <= Q.begin()->first) {
                // begin getting ready to return
                interval.emplace_back(C.begin()->second->page_id - MAXPAGEID, C.begin()->first);
                delete (C.begin()->second);
                C.erase(C.begin());
                return interval.back();
            }
        }
        NODE_PTR popped_node = heap.begin()->second;
        PAGE_ID pageID = popped_node->page_id;
        heap.erase(heap.begin());
        if (interval.size() >= k) {
            S.erase(popped_node);
        }
        if (pageID > MAXPAGEID)  // option processing
        {
            if (interval.size() < k)  // Phase (1)
            {
                interval.emplace_back(pageID - MAXPAGEID, 0);
                if (interval.size() == k) // Phase (2)
                {
                    //init S, init Q
                    for (pair<PG_IDX, RADIUS> &option:interval) {
                        incompSet.push_back(option.first);
                    }
                    for (pair<const DOT, NODE_PTR> &ele:heap) {
                        NODE_PTR node = ele.second;
                        node->init_radius(incompSet, PG, userpref, k);
                        S.insert(node);
                        Q.emplace(node->radius_LB(), node);
                    }
                }
                return interval.back();
            } else {
                // todo UPDATE WHEN needed
                popped_node->update_radius(incompSet.begin() + popped_node->tau, incompSet.end(), PG, userpref);
                C.emplace(popped_node->radius_LB(), popped_node);
                incompSet.push_back(pageID - MAXPAGEID);
            }
        } else // internal and leaf nodes processing
        {
            float possible_next_radius = INFINITY;
            if (!C.empty()) {
                possible_next_radius = C.begin()->first;
            }
            RtreeNode *node = ramTree[pageID];
            if (node->isLeaf()) {
                for (int i = 0; i < node->m_usedspace; i++) {
                    DOT tmpScore = 0;
                    for (int j = 0; j < dimen; j++) {
                        pt[j] = node->m_entry[i]->m_hc.getCenter()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    unknown_X_node *tmp_node = new unknown_X_node(dimen, PG[node->m_entry[i]->m_id],
                                                                  node->m_entry[i]->m_id + MAXPAGEID);
                    heap.emplace(tmpScore, tmp_node);
                    if (interval.size() >= k) {
                        S.insert(tmp_node);
                        tmp_node->init_radius(incompSet, PG, userpref, k, possible_next_radius);
                        Q.emplace(tmp_node->radius_LB(), tmp_node);
                    }
                }
            } else {
                for (int i = 0; i < node->m_usedspace; i++) {
                    DOT tmpScore = 0;
                    for (int j = 0; j < dimen; j++) {
                        pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    const float *ptr = node->m_entry[i]->m_hc.getUpper().m_coor;
                    unknown_X_node *tmp_node = new unknown_X_node(dimen, ptr, node->m_entry[i]->m_id);
                    heap.emplace(tmpScore, tmp_node);
                    if (interval.size() >= k) {
                        S.insert(tmp_node);
                        tmp_node->init_radius(incompSet, PG, userpref, k, possible_next_radius);
                        Q.emplace(tmp_node->radius_LB(), tmp_node);
                    }
                }
            }

        }

    }
    return {-1, INFINITY};
}

vector<int> build_qhull(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes){
    // \tpara ITERATABLE set<int>, unorder_set<int>, vector<int> and other iteratable STL<INT> CLASS
    int dim=square_vertexes[0].size();
    for (int opt_idx:opt_idxes) {
        for (int i = 0; i < dim; ++i) {
            square_vertexes[i][i] = max(square_vertexes[i][i],PG[opt_idx][i]);
        }
    }
    string s = to_string(dim - 1) + " " + to_string(opt_idxes.size() + square_vertexes.size()) + " ";
    for(int opt_idx:opt_idxes){
        for (int i = 0; i <dim ; ++i) {
            s += to_string(PG[opt_idx][i]) + " ";
        }
    }
    for (vector<double> & square_vertex : square_vertexes){
        for (float j : square_vertex){
            s += to_string(j) + " ";
        }
    }
    istringstream is(s);
    RboxPoints rbox;
    rbox.appendPoints(is);
    vector<int> CH;
    Qhull q(rbox, "");
    for(auto qiter=q.beginVertex();qiter!=q.endVertex();qiter=qiter.next()){
        if(qiter.id()<opt_idxes.size()){
            CH.push_back(opt_idxes[qiter.id()]);
        }
    }
    return CH;
}

vector<pair<int, vector<vector<double>>>> top_region(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes){
    //
    int dim=square_vertexes[0].size();
    for (int opt_idx:opt_idxes) {
        for (int i = 0; i < dim; ++i) {
            square_vertexes[i][i] = max(square_vertexes[i][i], PG[opt_idx][i]);
        }
    }
    string s = to_string(dim - 1) + " " + to_string(opt_idxes.size() + square_vertexes.size()) + " ";
    for(int opt_idx:opt_idxes){
        for (int i = 0; i <dim ; ++i) {
            s += to_string(PG[opt_idx][i]) + " ";
        }
    }
    for (vector<double> & square_vertex : square_vertexes){
        for (float j : square_vertex){
            s += to_string(j) + " ";
        }
    }
    istringstream is(s);
    RboxPoints rbox;
    rbox.appendPoints(is);
    vector<pair<int, vector<vector<double>>>> CH;
    Qhull q(rbox, "");
    for(auto qiter=q.beginVertex();qiter!=q.endVertex();qiter=qiter.next()){
        if(qiter.id()>=opt_idxes.size()){
            continue;
        }
        auto hyperplanes = qiter.neighborFacets().toStdVector();
        vector<vector<double>> norms;
        for(auto & hyperplane : hyperplanes){
            double *h=hyperplane.hyperplane().coordinates();
            norms.emplace_back(h, h+dim);
        }
        CH.emplace_back(opt_idxes[qiter.id()], norms);
    }
    return CH;
}


void build_qhull(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes, Qhull *q_ptr){
    // \tpara ITERATABLE set<int>, unorder_set<int>, vector<int> and other iteratable STL<INT> CLASS
    int dim=square_vertexes[0].size();
    for (int opt_idx:opt_idxes) {
        for (int i = 0; i < dim; ++i) {
            square_vertexes[i][i] = max(square_vertexes[i][i], PG[opt_idx][i]);
        }
    }
    string s = to_string(dim - 1) + " " + to_string(opt_idxes.size() + square_vertexes.size()) + " ";
    for(int opt_idx:opt_idxes){
        for (int i = 0; i <dim ; ++i) {
            s += to_string(PG[opt_idx][i]) + " ";
        }
    }
    for (vector<double> & square_vertex : square_vertexes){
        for (float j : square_vertex){
            s += to_string(j) + " ";
        }
    }
    istringstream is(s);
    RboxPoints rbox;
    rbox.appendPoints(is);
    unordered_map<int, vector<int>> CH;
    q_ptr->runQhull(rbox, "");
//    return q_ptr;
}

inline void update_square_vertexes(vector<vector<double>> &square_vertexes, float *new_ele, int dim){
    for (int i=0;i<dim;++i){
        square_vertexes[i][i]=max(square_vertexes[i][i], new_ele[i]);
    }
}

class ch{
    unordered_set<int> rest;
    unordered_map<int, QhullVertex*> pdtid_v;
    unordered_map<int, int> pdtid_layer;
    unordered_map<QhullVertex*, int> v_pdtid;
    unordered_map<int, vector<int>> A_p;
    vector<vector<int>> chs;
    vector<Qhull *> ch_q;
    vector<int> EMPTY;
    float** pointSet;
    int d;
    public:
    ch(vector<int> &idxes, float** &point_set, int dim){
        this->rest=unordered_set<int>(idxes.begin(), idxes.end());
        this->pointSet=point_set;
        this->d=dim;
    }
    const vector<int>& get_next_layer(){
        vector<vector<double>> square_vertexes(d+1, vector<double>(d));
        vector<int> rest_v(rest.begin(), rest.end());
        Qhull *ptr=new Qhull();
        build_qhull(rest_v, pointSet, square_vertexes, ptr);
        vector<int> ch;
        for(auto viter=ptr->beginVertex();viter!=ptr->endVertex();viter=viter.next()){
            ch.push_back(rest_v[viter.id()]);
        }
        for(auto viter=ptr->beginVertex();viter!=ptr->endVertex();viter=viter.next()){
            pdtid_v[rest_v[viter.id()]]=&viter;
            v_pdtid[&viter]=rest_v[viter.id()];
            pdtid_layer[rest_v[viter.id()]]=chs.size()+1;
        }
        chs.push_back(ch);
        ch_q.push_back(ptr);
        for (int idx:ch) {
            rest.erase(idx);
        }
        return chs.back();
    }

    int get_option_layer(int option){
        auto iter=pdtid_layer.find(option);
        if(iter!=pdtid_layer.end()){
            return iter->second;
        }else{
            return -1; // not in current i layers
        }
    }

    vector<int> get_neighbor_vertex(int option){
        auto lazy_get=A_p.find(option);
        if(lazy_get!=A_p.end()){
            return lazy_get->second;
        }
        auto iter=pdtid_v.find(option);
        if(iter==pdtid_v.end()){
            vector<int> ret;
            return ret;
        }
        QhullVertex *v=iter->second;
        QhullFacetSet fs=v->neighborFacets();
        set<int> before_ret;
        for (auto fiter=fs.begin();fiter!=fs.end();++fiter) {
            auto f=*fiter;
            auto vs=f.vertices();
            for (auto viter=vs.begin();viter!=vs.end();++viter) {
                auto nv=*viter;
                auto pdt_id=v_pdtid.find(&nv);
                if(pdt_id!=v_pdtid.end()){
                    before_ret.insert(pdt_id->second);
                }
            }
        }
        before_ret.erase(option);
        vector<int> ret(before_ret.begin(), before_ret.end());
        A_p[option]=ret;
        return ret;
    }

    const vector<int>& get_layer(int which_layer){
        // layer is starting from 1
        while(chs.size()<which_layer && !rest.empty()){
            this->get_next_layer();
        }
        if(chs.size()<which_layer && rest.empty()){
            return EMPTY;
        }
        return this->chs[which_layer-1];
    }

    ~ch(){
        for(Qhull *q:ch_q){
            delete (q);
        }
    }
};

class region{
public:
    vector<int> topk;
    double radius;
    vector<vector<double>> cone;
public:
    region(vector<int> &topK, vector<vector<double>> &Cone){
        topk=topK;
        cone=Cone;
    }
    void setRadius(float new_r){
        this->radius=new_r;
    }
    void setRadius(double new_r){
        this->radius=new_r;
    }
};


bool region_overlap(vector<vector<double>> &r1, vector<vector<double>> &r2){
    return isFeasible(r1, r2);
}

float dist_region_w(vector<vector<double>> &region, vector<float> &w){
    qp_solver q(w, region);
    return q.solve();
}

void topRegions(vector<vector<double>> &parent_region, const vector<int> &CH_upd, ch &ch_obj,
                multimap<double, region*> &id_radius, int deepest_layer, float **PG, int dim,
                const int k, vector<int> &topi, vector<float> &w, vector<int> &neighbors){
    //dfs
//    vector<pair<int, vector<vector<float>>>> top_region(vector<int> &opt_idxes, float **PG, vector<vector<float>> &square_vertexes){
    if(topi.size()==k){
        region *r=new region(topi, parent_region);
        id_radius.emplace(dist_region_w(parent_region, w), r);
    }
    int d=dim;
    vector<vector<double>> square_vertexes(d+1, vector<double>(d));
    vector<pair<int, vector<vector<double>>>> tops_region=top_region(CH_upd, PG, square_vertexes);
    for(pair<int, vector<vector<double>>> &opt_r: tops_region){
        if(region_overlap(parent_region, opt_r.second)){
            // update topi
            vector<int> topi1=topi;
            topi1.push_back(opt_r.first);
            // append neighbor, should be a set in nature
            vector<int> new_nei=neighbors;
            vector<int> aps=ch_obj.get_neighbor_vertex(opt_r.first);
            for(int ap:aps){
                new_nei.push_back(ap);
            }
            // sub-region
            vector<vector<double>> new_region=parent_region;
            for(vector<double> &opt:opt_r.second){
                new_region.push_back(opt);
            }
            // judge whether need next CH layers
            int deeper=max(deepest_layer, ch_obj.get_option_layer(opt_r.first));
            vector<int> update_s(new_nei.begin(), new_nei.end());
            vector<int> CH_mp1=ch_obj.get_layer(deeper+1);
            for(int next_layer_opt:CH_mp1){
                update_s.push_back(next_layer_opt);
            }
            topRegions(new_region, CH_mp1, ch_obj,id_radius, deeper, PG, dim,k, topi1, w, new_nei);
        }
    }

}


vector<vector<double>> points_to_halfspace(vector<vector<double>> &points){
    // be careful such that the norms are pointing out the convex cone
    // which means the convex cone is represented as
    // n_1 \cdot w <= 0
    // n_2 \cdot w <= 0
    // ...
    // n_a \cdot w <= 0
    int dim=points[0].size();
    vector<vector<double>> square_vertexes;
    square_vertexes.emplace_back(dim, 0.0);
    square_vertexes.emplace_back(dim, 1.0);

    string s = to_string(dim - 1) + " " + to_string(points.size() + square_vertexes.size()) + " ";
    for (vector<double> & square_vertex : square_vertexes){
        for (float j : square_vertex){
            s += to_string(j) + " ";
        }
    }
    for(vector<double> &point: points){
        for (float i:point) {
            s += to_string(i) + " ";
        }
    }
    istringstream is(s);
    RboxPoints rbox;
    rbox.appendPoints(is);
    Qhull q(rbox, "");
    auto oviter= q.endVertex();
    for(auto viter=q.beginVertex();viter!=q.endVertex();viter=viter.next()){
        if(viter.id()==0){
            oviter=viter;
        }
    }
    assert(oviter!=q.endVertex());
    auto hyperplanes = oviter.neighborFacets().toStdVector();
    vector<vector<double>> norms;
    for(auto & hyperplane : hyperplanes){
        double *h=hyperplane.hyperplane().coordinates();
        norms.emplace_back(h, h+dim);
    }
    return norms;
}

void utk_basic(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
               vector<pair<int, double>> &utk_option_ret,
               vector<pair<vector<int>, vector<vector<double>>>> &utk_cones_ret){

    // 2 return value
    // 1. array of <option, topk_radius>
    // 2. array of <topk, region>
    // "apply k=1"
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, PointSet);
    pair<int, float> next={-1, INFINITY};
    //fetch top X options
    vector<int> CH_1_X_opt;
    while(CH_1_X_opt.size() < X){  // fetch top X
        next=get_next_obj.get_next();
        CH_1_X_opt.push_back(next.first);
    }
    // qhull class in lib qhull
    const int square_vertex_cnt=dim+1;
    vector<vector<double>> square_vertexes(square_vertex_cnt, vector<double>(dim));

    // init qhull with top X options
    CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);

    // a 3-d example of square_vertex, square_vertex_cnt=4
    // point 0: (max(points[:, 0]), 0, 0)
    // point 1: (0, max(points[:, 1]), 0)
    // point 2: (0, 0, max(points[:, 1]))
    // point 3: \vec{0}

    // rho_star is computed as when CH_1.size()=X
    while(CH_1_X_opt.size() < X){ // while(CH_1.size<X)
        while(CH_1_X_opt.size() < X){
            next=get_next_obj.get_next();
            update_square_vertexes(square_vertexes, PointSet[next.first], dim);
            CH_1_X_opt.push_back(next.first);
        }
        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
    }
    // for now, qhull_obj contains points 0~3 and convex hull vertexes
    float rho_star=next.second;

    // use known X version code to fetch rskyband options,
    // bear in mind such that we init \rho as \rho_star and X as INFINITY
    vector<pair<long int, float>> interval;
    computeRho(dim, k, INFINITY, w, *rtree, PointSet, interval, rho_star);
    vector<int> rskyband_CS;
    for (pair<long int, float> &p:interval) {
        rskyband_CS.push_back(p.first);
    }
    ch ch_obj(rskyband_CS, PointSet, dim);
    vector<int> top1_idxes=ch_obj.get_layer(1);
    vector<vector<double>> tmp;
    double rho_star_d=rho_star; // change from float to double
    for (vector<c_float> &e:g_r_domain_vec) {
        tmp.push_back(rho_star_d*e+w);
    }
    vector<vector<double>> begin_region=points_to_halfspace(tmp);
    multimap<double, region*> id_radius; // <radius, region>
//    class region{
//    public:
//        vector<int> topk;
//        float radius;
//        vector<vector<float>> cone;
//    };
    vector<int> init_topi;
    vector<int> init_neighbors;
    topRegions(begin_region, CH_1_X_opt, ch_obj, id_radius, 1,  PointSet, dim,
                    k, init_topi, w, init_neighbors);
    // until X different options
//    vector<pair<int, float>> utk_option_ret;
//    vector<pair<vector<int>, vector<vector<float>>>> utk_cones_ret;
    assert(!id_radius.empty());
    _Rb_tree_iterator<pair<const double, region *>> iter=id_radius.begin();
    unordered_set<int> options;
    while(options.size()<X){
        for (int option_idx: iter->second->topk) {
            if(options.find(option_idx)==options.end()){ // new option
               options.insert(option_idx);
               utk_option_ret.emplace_back(option_idx, iter->first);
            }
        }
        utk_cones_ret.emplace_back(iter->second->topk, iter->second->cone);
        ++iter;
    }

}
