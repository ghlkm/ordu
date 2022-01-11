#include "iPref.h"
//#include "qp_solver.h"
#include "math_lib.h"
#include "qhull_user.h"
#include "lp_user.h"
#include <chrono>
#include "qp_solver2.h"

extern int objCnt;

extern qp_solver qp;

inline double max(double a, float b){
    return a>b?a:b;
}

inline double max(float b, double a){
    return a>b?a:b;
}

inline long min(int a, long b){
    return a<b?a:b;
}

inline long min(long a, int b){
    return a<b?a:b;
}


vector<float> computePairHP(const int dimen, float* PG[], long int pi, long int pj)
{
    vector<float> retHS(dimen);
    for (int i = 0; i < dimen; ++i) {
        retHS[i]=PG[pi][i]-PG[pj][i];
    }
    return retHS;
}

float computeDis(const vector<float> &tmpHS, const vector<float> &userpref)
{
    float ret=qp.update_w_h_solve(userpref, tmpHS);
    return ret;
}

float computeradius(const int k, const int dim, long int pi, vector<float>& userpref, vector<long int>& incompSet, float* PG[], float rho){
    multiset<float> radiusSKI;
    float tmpDis;
    for (int ri = 0; ri < incompSet.size(); ri++){
        if (v1_dominate_v2(PG[incompSet[ri]],PG[pi], dim)){// for these dominators, we directly input as FLTMAX
            radiusSKI.insert(INFINITY);
        }
        else{
            vector<float> tmpHS = computePairHP(dim, PG, pi, incompSet[ri]);
            tmpDis = computeDis(tmpHS, userpref);
            radiusSKI.insert(tmpDis);
        }
        if(radiusSKI.size()>k){
            radiusSKI.erase(radiusSKI.begin());
            if(*radiusSKI.begin()>rho){
                // the default rho is INF
                // if current radius is larger than rho, than prune it
                break;
            }
        }
    }
    if(radiusSKI.size()>=k){
        return *radiusSKI.begin();
    }else{
        return 0;
    }
}

bool v2_r_dominate_v1_float(const float* v1, const float* v2, const vector<float> &w, const vector<vector<c_float>> &r_domain_vec, const float &rho) {
    for (const vector<c_float> &v:r_domain_vec) {
//        vector<float> tmp_w = w + rho * v;
//        if (v1 * tmp_w < v2 * tmp_w) {
//            return false;
//        }
//        the below code is same as above
        double atc_rho=rho;
        for (int i = 0; i < v.size(); ++i) {
            if(v[i]<0){
                atc_rho=min(atc_rho, -w[i]/v[i]); // in case of w[i] + \rho * v[i] <0 or >1
            }
        }
        vector<float> tmp_w(w.size());
        for (int i = 0; i < tmp_w.size(); ++i) {
            tmp_w[i]=w[i]+atc_rho*v[i];
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
    vector<long int> incompSet;
    pair<float, int> candidateOpt;

    RtreeNode* node;
    multimap<float, int, greater<>> heap;
    multimap<float, int, greater<>> candidateRet;

    float pt[dimen];
    int pageID;
    float tmpScore; // compute the score of option or mbr e w.r.t. userpref
    float tmpRadius; // the inflection radius of point Pi.

    heap.emplace(INFINITY, a_rtree.m_memory.m_rootPageID);

    while (!heap.empty()){
        tmpScore = heap.begin()->first;
        pageID = heap.begin()->second;
        heap.erase(heap.begin());
        if (pageID > MAXPAGEID){  // option processing
            if (interval.size() < k){  // Phase I
                interval.emplace_back(pageID - MAXPAGEID, 0);
                incompSet.push_back(pageID - MAXPAGEID);
            }
            else{
                if (candidateRet.size() < X - k){  // Phase II
                    tmpRadius = computeradius(k, dimen, pageID - MAXPAGEID, userpref, incompSet, PG, radius);
                    if(tmpRadius!=INFINITY){
                        candidateRet.emplace(tmpRadius, pageID-MAXPAGEID);
                        if(candidateRet.size()==X-k){
                            radius = min(radius, candidateRet.begin()->first);
                        }
                    }
                }
                else if(X<=k){
                    assert(X==k);// if fails there is a problem in data
                    radius=0;
                    break;
                }
                else{   // Phase III
                    tmpRadius = computeradius(k, dimen, pageID - MAXPAGEID, userpref, incompSet, PG, candidateRet.begin()->first);
                    if (tmpRadius < candidateRet.begin()->first){
                        candidateRet.emplace(tmpRadius, pageID - MAXPAGEID);
                        candidateRet.erase(candidateRet.begin());
                        candidateOpt = *candidateRet.begin();
                        radius = min(candidateOpt.first, radius);
                    }
                }
                if(tmpRadius!=INFINITY){
                    incompSet.push_back(pageID - MAXPAGEID);
                }
            }
//			incompSet.push_back(pageID - MAXPAGEID);

        }
        else{ // internal and leaf nodes processing
            node = ramTree[pageID];
            if (node->isLeaf()){
                for (int i = 0; i < node->m_usedspace; i++){
                    tmpScore = 0;
                    for (int j = 0; j < dimen; j++){
                        pt[j] = node->m_entry[i]->m_hc.getCenter()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    if (radius == INFINITY){
                        if (!dominatedByK(dimen, pt, incompSet, PG, k)){
                            heap.emplace(tmpScore, node->m_entry[i]->m_id + MAXPAGEID);
                        }
                    }
                    else{
                        if (!r_dominatedByK(dimen, pt, radius, userpref, incompSet, PG, k)){
                            heap.emplace(tmpScore, node->m_entry[i]->m_id + MAXPAGEID);
                        }
                    }
                }
            }
            else{
                for (int i = 0; i < node->m_usedspace; i++){
                    tmpScore = 0;
                    for (int j = 0; j < dimen; j++){
                        pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
                        tmpScore += pt[j] * userpref[j];
                    }
                    if (radius == INFINITY){
                        if (!dominatedByK(dimen, pt, incompSet, PG, k)){
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
                        }
                    }
                    else{
                        if (!r_dominatedByK(dimen, pt, radius, userpref, incompSet, PG, k)){
                            heap.emplace(tmpScore, node->m_entry[i]->m_id);
                        }
                    }
                }
            }
        }
    }
    for (auto riter=candidateRet.rbegin();riter!=candidateRet.rend();++riter) {
        if(riter->first>radius){
            break;
        }
        interval.emplace_back(riter->second, riter->first);
    }
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

unknown_x_baseline::unknown_x_baseline(const int dim, const int K, vector<float>& userPref, Rtree& a_rtree, float** pg){
    dimen=dim;
    ones=vector<float>(dimen, 1);
    zeros=vector<float>(dimen, 0);
    zeros_node=new unknown_X_node(dimen, zeros.data(), -1);
    rt=new unknown_X_node(dimen, ones.data(), a_rtree.m_memory.m_rootPageID);
    heap.emplace(INFINITY, rt);
    to_dl.push_back(zeros_node);
    to_dl.push_back(rt);
    k=K;
    userpref=userPref;
    PG=pg;
    next=0;
}

pair<int, float> unknown_x_baseline::get_next(){
    next++;
    if(interval.size()>=next){
        return interval[next-1];
    }
    while (!heap.empty())
    {
        if(interval.size()>=next){
            return interval[next-1];
        }
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
                return interval.back();
            }
            else
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
                incompSet.push_back(pageID - MAXPAGEID);
            }
        }
        else // internal and leaf nodes processing
        {
            S.erase(popped_node);
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
                    to_dl.push_back(tmp_node);
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
                    to_dl.push_back(tmp_node);
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
}

unknown_x_baseline::~unknown_x_baseline(){
    for(unknown_X_node *node:to_dl){
        delete (node);
    }
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
    vector<unknown_X_node*> to_dl; // todo change into intel pointer
    to_dl.push_back(zeros_node);
    to_dl.push_back(rt);
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
//            delete(popped_node);
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
                    to_dl.push_back(tmp_node);
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
                    to_dl.push_back(tmp_node);
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
//                delete(LB.second);
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

    for(unknown_X_node *node:to_dl){
        delete (node);
    }
//    delete(zeros_node);
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
    multimap<DOT, NODE_PTR, greater<DOT>> heap; //BBS max_heap, <w\cdot node, node>, if node is an option then node=id+MAXPAGEID
    unordered_set<NODE_PTR> S; // reflect which nodes and options in BBS
    multimap<RADIUS, NODE_PTR, less<RADIUS>> Q; // min_heap based on inflection radius, lazy update for S
    multimap<RADIUS, NODE_PTR, less<RADIUS>> C; // min_heap based on inflection radius, candidate list
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

inline bool unknown_X_node::update() const{
    return this->tau>=this->needed_to_update;
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
//                assert(C.begin()->second->page_id>MAXPAGEID);
                delete (C.begin()->second);
                C.erase(C.begin());
                while (C.size() > 1 && !C.begin()->second->update()) {
                    unknown_X_node *ni = C.begin()->second;
                    float lb = (++C.begin())->first;
                    ni->update_radius(incompSet.begin() + ni->tau, incompSet.begin() + ni->needed_to_update, PG, userpref, lb);
                    C.erase(C.begin());
                    C.emplace(ni->radius_LB(), ni);
                }
                if (!C.empty() &&!C.begin()->second->update()) {
                    unknown_X_node *ni = C.begin()->second;
                    ni->update_radius(incompSet.begin() + ni->tau, incompSet.begin() + ni->needed_to_update, PG, userpref);
                    C.erase(C.begin());
                    C.emplace(ni->radius_LB(), ni);
                }
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
                delete (popped_node);
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
                // UPDATE WHEN needed
                if(!C.empty()) {
                    popped_node->update_radius(incompSet.begin() + popped_node->tau, incompSet.end(), PG, userpref,
                                               C.begin()->first);
                }else{
                    popped_node->update_radius(incompSet.begin() + popped_node->tau, incompSet.end(), PG, userpref);
                }
                popped_node->needed_to_update=incompSet.size();
                C.emplace(popped_node->radius_LB(), popped_node);
                while (C.size() > 1 && !C.begin()->second->update()) {
                    unknown_X_node *ni = C.begin()->second;
                    float lb = (++C.begin())->first;
                    ni->update_radius(incompSet.begin() + ni->tau, incompSet.begin() + ni->needed_to_update, PG, userpref, lb);
                    C.erase(C.begin());
                    C.emplace(ni->radius_LB(), ni);
                }
                if (!C.empty() &&!C.begin()->second->update()) {
                    unknown_X_node *ni = C.begin()->second;
                    ni->update_radius(incompSet.begin() + ni->tau, incompSet.begin() + ni->needed_to_update, PG, userpref);
                    C.erase(C.begin());
                    C.emplace(ni->radius_LB(), ni);
                }
                if(popped_node->radius_LB()!=INFINITY)
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

vector<int> non_dominated(const vector<int> &opt_idxes, float **PG, int dim){
    vector<int> ret;
    for (int i: opt_idxes) {
        bool flag=true;
        for (int j: opt_idxes) {
            if(i==j){
                continue;
            }
            if(v1_dominate_v2(PG[j], PG[i], dim)){
                flag=false;
                break;
            }
        }
        if(flag){
            ret.push_back(i);
        }
    }
    return ret;
}

vector<int> build_qhull(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes){
    int dim=square_vertexes[0].size();
    auto begin = chrono::steady_clock::now();
    square_vertexes.clear();
    square_vertexes.emplace_back(dim);
    vector<int> nd_opt_idxes=non_dominated(opt_idxes, PG, dim); // must non-dominated each other
    string s = to_string(dim) + " " + to_string(nd_opt_idxes.size() + square_vertexes.size()) + " ";
    for(int opt_idx:nd_opt_idxes){
        for (int i = 0; i <dim ; ++i) {
            if(PG[opt_idx][i]+SIDELEN < 1e-6){
                s += to_string(SIDELEN)+ " ";// in case of precision problem
            }else{
                s += to_string(PG[opt_idx][i]+SIDELEN) + " ";
            }
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
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    cout<<"input time:"<<elapsed_seconds.count()<<endl;
    // 性能瓶颈
    Qhull q(rbox, "QJ");
    // 性能瓶颈
    auto now2 = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds2= now2-begin;
    cout<<"qhull time:"<<elapsed_seconds2.count()<<endl;
    qhull_user qu;
    return qu.get_CH_pointID(q, nd_opt_idxes);
}

vector<int> build_qhull(const vector<int> &opt_idxes, vector<vector<float>> &PG, vector<vector<double>> &square_vertexes){
    // \tpara ITERATABLE set<int>, unorder_set<int>, vector<int> and other iteratable STL<INT> CLASS
    int dim=square_vertexes[0].size();
    square_vertexes.clear();
    square_vertexes.emplace_back(dim);
    for (int opt_idx:opt_idxes) {
        for (int i = 0; i < dim; ++i) {
            vector<double> tmp(PG[opt_idx].begin(), PG[opt_idx].end());
            for(auto &t:tmp){
                t+=SIDELEN;
                t=max(t, 0);
            }
            tmp[i]=0;
            square_vertexes.push_back(tmp);
        }
    }
    string s = to_string(dim) + " " + to_string(opt_idxes.size() + square_vertexes.size()) + " ";
    for(int opt_idx:opt_idxes){
//        assert(opt_idx<=objCnt);
        for (int i = 0; i <dim ; ++i) {
            if(PG[opt_idx][i]+SIDELEN < 1e-6){
                s += to_string(SIDELEN)+ " ";// in case of precision problem
            }else{
                s += to_string(PG[opt_idx][i]+SIDELEN) + " ";
            }
        }
    }
    for (vector<double> & square_vertex : square_vertexes){
//        assert(square_vertex.size()==dim);
        for (float j : square_vertex){
            s += to_string(j) + " ";
        }
    }
    istringstream is(s);
    RboxPoints rbox;
    rbox.appendPoints(is);
    Qhull q(rbox, "QJ");
    qhull_user qu;
    return qu.get_CH_pointID(q, opt_idxes);
}

void top_region(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes,
                unordered_map<int, vector<vector<double>>> &ret){
    int dim=square_vertexes[0].size();
//    for (int opt_idx:opt_idxes) {
//        for (int i = 0; i < dim; ++i) {
//            square_vertexes[i][i] = max(square_vertexes[i][i], PG[opt_idx][i]+SIDELEN);
//        }
//    }
//    while(square_vertexes.size()>dim+1){
//        square_vertexes.pop_back();
//    }
    square_vertexes.clear();
    square_vertexes.emplace_back(dim);
    for (int opt_idx:opt_idxes) {
        for (int i = 0; i < dim; ++i) {
            vector<double> tmp(PG[opt_idx], PG[opt_idx]+dim);
            for(auto &t:tmp){
                t+=SIDELEN;
                t=max(t, 0);
            }
            tmp[i]=0;
            square_vertexes.push_back(tmp);
        }
    }
    string s = to_string(dim) + " " + to_string(opt_idxes.size() + square_vertexes.size()) + " ";
    for(int opt_idx:opt_idxes){
        for (int i = 0; i <dim ; ++i) {
//            assert(opt_idx>=0 && opt_idx<=objCnt);
            if(PG[opt_idx][i]+SIDELEN < 1e-6){
                s += to_string(SIDELEN)+ " ";
            }else{
                s += to_string(PG[opt_idx][i]+SIDELEN) + " "; // in case of precision problem
            }
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
    // 性能瓶颈
    Qhull q(rbox, "QJ");
    // 性能瓶颈
    qhull_user qu;
    unordered_map<int, vector<vector<double>>> points;
    qu.get_neiFacetsNorm_of_point(q, opt_idxes, points);
    for (auto &pt:points) {
        ret[pt.first]=points_to_halfspace(pt.second);
    }
}



void build_qhull(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes, Qhull *q_ptr){
    int dim=square_vertexes[0].size();
    square_vertexes.clear();
    square_vertexes.emplace_back(dim);
    string s = to_string(dim) + " " + to_string(opt_idxes.size() + square_vertexes.size()) + " ";
    for(int opt_idx:opt_idxes){
        for (int i = 0; i <dim ; ++i) {
            if(PG[opt_idx][i]+SIDELEN < 1e-6){
                s += to_string(SIDELEN)+ " ";// in case of precision problem
            }else{
                s += to_string(PG[opt_idx][i]+SIDELEN) + " ";
            }
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
    q_ptr->runQhull(rbox, "QJ");
}

unordered_map<int, vector<vector<double>>> top_region(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes){
    //
    unordered_map<int, vector<vector<double>>> ret;
    top_region(opt_idxes, PG, square_vertexes, ret);
    return ret;
}

unordered_map<int, vector<vector<double>>> top_region2(const vector<int> &opt_idxes, float **PG, vector<vector<double>> &square_vertexes){
    //
    Qhull q;
    build_qhull(opt_idxes, PG, square_vertexes, &q);
    qhull_user qu;
    vector<int> ch=qu.get_CH_pointID(q, opt_idxes);
    unordered_map<int, vector<int>> reti;
    int dim=square_vertexes[0].size();
    qu.get_neiVT_of_VT(q, opt_idxes, reti);
    unordered_map<int, vector<vector<double>>> ret;
    for(int i:ch){
        auto iter=reti.find(i);
        if(iter!=reti.end()){
            ret[i]=vector<vector<double>>();
            for (int nei:iter->second) {
                vector<double> tmp(PG[nei], PG[nei]+dim);
                for (int j = 0; j <dim ; ++j) {
                    tmp[j]-=PG[i][j];
                }
                ret[i].push_back(tmp);
            }
        }
    }
    return ret;
}


inline void update_square_vertexes(vector<vector<double>> &square_vertexes, float *new_ele, int dim){
    for (int i=0;i<dim;++i){
        square_vertexes[i][i]=max(square_vertexes[i][i], new_ele[i]);
    }
}

class ch{
    unordered_set<int> rest;
    unordered_map<int, int> pdtid_layer;
    unordered_map<int, int> dominated_cnt;// only use for build k-convex-hull

    vector<vector<int>> chs;

    vector<int> EMPTY;
    float** pointSet;
    int d;
public:
    unordered_set<int> all;
    vector<int> rskyband;
    unordered_map<int, vector<int>> A_p;
    unordered_map<int, unordered_set<int>> do_map;
    unordered_map<int, unordered_set<int>> dominated_map;

    ch(vector<int> &idxes, float** &point_set, int dim){
        this->rskyband=idxes;
        this->rest=unordered_set<int>(idxes.begin(), idxes.end());
        this->all=unordered_set<int>(idxes.begin(), idxes.end());
        this->pointSet=point_set;
        this->d=dim;
        build_do_re(idxes, point_set, dim);
    }

    void fast_non_dominate_sort(
            const unordered_map<int, unordered_set<int>> &do_map,
            unordered_map<int, int>& dominated_cnt,
            const vector<int> &last_layer){
        for (int opt:last_layer) {
            auto iter=do_map.find(opt);
            if(iter!=do_map.end()){
                for(int dominated:iter->second){
                    auto cnt_iter=dominated_cnt.find(dominated);
                    if(cnt_iter!=dominated_cnt.end()){
                        cnt_iter->second-=1;
                    }
                }
            }
        }
    }

    void build_do_re(vector<int> &idxes, float** &point_set, int dim){
        for (int i:idxes){
            dominated_cnt[i]=0;
            do_map[i]=unordered_set<int>();
            dominated_map[i]=unordered_set<int>();
        }
        for (int ii = 0;ii<idxes.size();++ii) {
            int i=idxes[ii];
            for (int ji = ii+1; ji <idxes.size() ; ++ji) {
                int j=idxes[ji];
                if(v1_dominate_v2(point_set[i], point_set[j], dim)){
                    do_map[i].insert(j);
                    dominated_map[j].insert(i);
                    ++dominated_cnt[j];
                }else if(v1_dominate_v2(point_set[j], point_set[i], dim)){
                    do_map[j].insert(i);
                    dominated_map[i].insert(j);
                    ++dominated_cnt[i];
//                }else{// non-dominate
                }
            }
        }
    }

    const vector<int>& get_next_layer(){
        vector<vector<double>> square_vertexes(d+1, vector<double>(d));
        if(!chs.empty()){
            fast_non_dominate_sort(do_map, dominated_cnt, chs[chs.size()-1]);
        }
        vector<int> rest_v;
        for(int i:rest){
            auto iter=dominated_cnt.find(i);
            if(iter!=dominated_cnt.end() && iter->second<=0){
                rest_v.push_back(i);
            }
        }
        cout<<"no. of points to build convex hull: "<<rest_v.size()<<endl;
        vector<int> ch;
        if(rest_v.size()>=d+1){
            Qhull q;
            qhull_user qu;
            auto begin=chrono::steady_clock::now();
            build_qhull(rest_v, pointSet, square_vertexes, &q);
            auto end=chrono::steady_clock::now();
            chrono::duration<double> elapsed_seconds= end-begin;
            cout<<"finish build convex hull: "<<elapsed_seconds.count()<<endl;
            ch=qu.get_CH_pointID(q, rest_v);
            qu.get_neiVT_of_VT(q, rest_v, A_p);
        }else{
            for(int i:rest_v){
                vector<int> tmp;
                for(int j:rest_v){
                    if(i!=j){
                        tmp.push_back(j);
                    }
                }
                A_p[i]=tmp;
                ch.push_back(i);
            }
        }
        chs.push_back(ch);
        for (int idx:ch) {
            pdtid_layer[idx] =  chs.size();
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

    const vector<int>& get_neighbor_vertex(int option){
//        assert(option>=0 && option <=objCnt);
        auto lazy_get=A_p.find(option);
        if(lazy_get!=A_p.end()){
            return lazy_get->second;
        }else{
            return EMPTY;
        }
    }

    const vector<int>& get_layer(int which_layer){
        // layer is starting from 1
        while(chs.size()<which_layer && !rest.empty()){
            this->get_next_layer();
        }
        if(chs.size()<which_layer || which_layer<=0){
            return EMPTY;
        }
        return this->chs[which_layer-1];
    }

    ~ch(){
    }
};

bool region_overlap(vector<vector<double>> &r1, vector<vector<double>> &r2, vector<vector<double>> &r1_r2){
    return isFeasible(r1, r2, r1_r2);
}

float dist_region_w(vector<vector<double>> &region, vector<float> &w){
    float ret=qp_solver2(w, region);
    return ret;
}

double dist_region_w(vector<vector<double>> &r1,vector<vector<double>> &r2,  vector<float> &w){
    double ret=qp_solver2(w, r1, r2);
    return ret;
}

void topRegions(vector<vector<double>> &parent_region, const vector<int> &CH_upd, ch &ch_obj,
                multimap<double, region*> &id_radius, int deepest_layer, float **PG, int dim,
                const int k, vector<int> &topi, vector<float> &w, vector<int> &neighbors){
    //dfs
    if(topi.size()==k){
        double dist=dist_region_w(parent_region, w);
        if(dist!=INFINITY){
            region *r=new region(topi, parent_region);
            id_radius.emplace(dist, r);
        }
        if(id_radius.size()%1000==0){
            if(true){
                cout<< id_radius.size()<<"\n";
                cout<<"top:";
                for (int i:topi) {
                    cout<<i<<" ";
                }
                cout<<"\n";
                cout<< deepest_layer << " " << CH_upd.size()<< " " << neighbors.size() << "\n";
            }
        }
        return;
    }
    int d=dim;
    vector<vector<double>> square_vertexes(d+1, vector<double>(d));
    unordered_map<int, vector<vector<double>>> tops_region=top_region(CH_upd, PG, square_vertexes);
    for(auto &opt_r: tops_region){
        vector<vector<double>> new_region;
        if(region_overlap(parent_region, opt_r.second, new_region)){
            // update topi
            vector<int> topi1=topi;
//            assert(opt_r.first>=0 && opt_r.first<=objCnt);
            topi1.push_back(opt_r.first);
            // append neighbor, should be a set in nature
            vector<int> new_nei=neighbors;
            const vector<int> &aps=ch_obj.get_neighbor_vertex(opt_r.first); // TOREAD
            for(int ap:aps){
                new_nei.push_back(ap);
            }
            // judge whether need next CH layers
            int deeper=max(deepest_layer, ch_obj.get_option_layer(opt_r.first));
            unordered_set<int> update_s(new_nei.begin(), new_nei.end());
            const vector<int> &CH_mp1=ch_obj.get_layer(deeper+1);
            for(int next_layer_opt:CH_mp1){
                update_s.insert(next_layer_opt);
            }
            for (int i:topi1) {
                update_s.erase(i);
            }
            vector<int> update(update_s.begin(), update_s.end());
            topRegions(new_region, update, ch_obj,id_radius, deeper, PG, dim,k, topi1, w, new_nei);
        }
    }

}

vector<int> r_dominate_filter(const vector<int> &top_ip1_cdd,
                              const vector<vector<double>> &R, float **PG);

int topRegions_efficient(vector<vector<double>> &parent_region, ch &ch_obj,
                         multimap<double, region*> &id_radius, float **PG, int dim, int X,
                         const int k, vector<float> &w, unordered_set<int> &top1_calculated,
                         vector<pair<int, double>> &utk_option_ret,
                         vector<pair<double, region*>> &utk_cones_ret,
                         unordered_map<int, vector<vector<double>>> &top1_region){
    // init "top1_calculated" with top1 respecting to w
    auto begin=chrono::steady_clock::now();
    unordered_set<int> options;
    for (int i:top1_calculated) {
        options.insert(i);
        utk_option_ret.emplace_back(i, 0);
        cout << "radius: " << 0 << "\n";
    }
    int cnt=0;
    int rcnt=0;
    cout<<parent_region.size()<<endl;
    while(options.size() < X && !id_radius.empty()){
        ++rcnt;
        pair<double, region*> popped=*id_radius.begin(); // must deep copy because next line we use erase
        id_radius.erase(id_radius.begin());
        bool new_option=False;
        bool newOption=False;
        if(popped.second->topk.size()==1){
            // if it is top1, push its adjacent vertex
            const vector<int> &top1_adj=ch_obj.get_neighbor_vertex(popped.second->topk.front());
            for (int adj_opt:top1_adj) {
                if(top1_calculated.find(adj_opt)!=top1_calculated.end()){
                    continue;
                }
                top1_calculated.insert(adj_opt);
                auto iter = top1_region.find(adj_opt);
                if(iter==top1_region.end()){
                    continue; // precision problem occurs!
                }
                double dist=dist_region_w(parent_region,iter->second, w);
                if(dist!=INFINITY){
                    vector<int> tmp;
                    tmp.push_back(adj_opt);
                    vector<vector<double>> new_region(parent_region);
                    for (vector<double> &its_own_topr: iter->second) {
                        new_region.push_back(its_own_topr);
                    }
                    region *r=new region(tmp, new_region);
                    id_radius.emplace(dist, r);
                    cnt++;
                }
            }
        }
        if(popped.second->topk.size()==k){ // a region that don't need to be partitioned
            for (int opt:popped.second->topk) {
                auto iter=options.find(opt);
                if(iter==options.end()){// new option
                    options.insert(opt);
                    utk_option_ret.emplace_back(opt, popped.first);
                    cout << "radius: " << popped.first << "\n";
                    if(true){
                        vector<double> tmp(PG[opt], PG[opt]+dim);
                        cout<<options.size()<<": "<<popped.second->cone.size()<<"," << popped.first << ", "<<utk_cones_ret.size()<<
                            ", "<<rcnt <<"," << ch_obj.get_neighbor_vertex(popped.second->topk.front()).size()<<" # ";
                        cout<<popped.second->topk<<":"<<tmp<<"\n";
                    }
                    auto now = chrono::steady_clock::now();
                    chrono::duration<double> elapsed_seconds= now-begin;
                    cout << "time: " << options.size() << ", " << elapsed_seconds.count()<<"\n";
                    new_option=True;
                }
            }
            if(new_option){
                popped.second->setRadius(popped.first);
                utk_cones_ret.emplace_back(popped);
            }
        }
        else{
            unordered_set<int> ch_upd_s;
            int m=0;
            for(int top: popped.second->topk){
                for(int adj: ch_obj.get_neighbor_vertex(top)){
                    ch_upd_s.insert(adj);
                }
                m=max(ch_obj.get_option_layer(top), m);
            }
            for(int mp1_opt: ch_obj.get_layer(m+1)){
                ch_upd_s.insert(mp1_opt);
            }
            for (int top: popped.second->topk) {
                ch_upd_s.erase(top);
            }

            unordered_map<int, int> dominated_cnt;
            for(int i:ch_upd_s){
                dominated_cnt[i]=ch_obj.dominated_map[i].size();
            }
            for(int opt:ch_upd_s){
                auto iter=ch_obj.dominated_map.find(opt);
                if(iter!=ch_obj.dominated_map.end()){
                    for(int i:popped.second->topk){
                        if(iter->second.find(i)!=iter->second.end()){
                            --dominated_cnt[opt];
                        }
                    }
                }
            }
            vector<int> ch_upd;
            for (int i:ch_upd_s) {
                auto iter=dominated_cnt.find(i);
                if(iter!=dominated_cnt.end()&&iter->second<=0){
                    ch_upd.push_back(i);
                }
            }

//            vector<int> ch_upd(ch_upd_s.begin(), ch_upd_s.end());// from set to vector
            // TODO
            const int square_vertex_cnt=dim+1;
            vector<vector<double>> square_vertexes(1, vector<double>(dim));
            // 时间主要开销， 主要性能瓶颈, 求凸包以及每个凸包点的top region
//            cout<<"option num:"<<ch_upd.size()<<endl;
            unordered_map<int, vector<vector<double>>> tops_region=top_region2(ch_upd, PG, square_vertexes);

            // 时间主要开销， 主要性能瓶颈, 求凸包以及每个凸包点的top region
            for (int ch_upd_opt: ch_upd) {
                auto iter = tops_region.find(ch_upd_opt);
                if(iter!=tops_region.end()){
                    double dist=dist_region_w(popped.second->cone, iter->second, w);
                    if(dist!=INFINITY){
                        vector<int> tmp(popped.second->topk);
                        tmp.push_back(ch_upd_opt);
                        vector<vector<double>> new_region(popped.second->cone);
                        // TODO here is deep copy maybe faster if shadow copy
                        for (vector<double> &its_own_topr: iter->second) {
                            new_region.push_back(its_own_topr);
                        }
                        region *r=new region(tmp, new_region);
                        id_radius.emplace(dist, r);
                        cnt++;
                    }
                }
            }

        }

        if(popped.second->topk.size()!=k){
            delete(popped.second);
        }
        else if(!new_option && !newOption){
            delete(popped.second);   // to save mem
        }
    }
    cout<<"cnt: "<<cnt<<"\n";
    for (auto left:id_radius) {
        delete (left.second);
    }
    return cnt;
}


int topRegions_efficient2(vector<vector<double>> &parent_region, ch &ch_obj,
                          multimap<double, region*> &id_radius, float **PG, int dim, int X,
                          const int k, vector<float> &w, unordered_set<int> &top1_calculated,
                          vector<pair<int, double>> &utk_option_ret,
                          vector<pair<double, region*>> &utk_cones_ret,
                          unordered_map<int, vector<vector<double>>> &top1_region){
    // init "top1_calculated" with top1 respecting to w
    auto begin=chrono::steady_clock::now();
    unordered_set<int> options;
    for (int i:top1_calculated) {
        options.insert(i);
        utk_option_ret.emplace_back(i, 0);
        cout << "radius: " << 0 << "\n";
    }
    int cnt=0;
    int rcnt=0;
    cout<<parent_region.size()<<endl;
    while(options.size() < X && !id_radius.empty()){
        ++rcnt;
        pair<double, region*> popped=*id_radius.begin(); // must deep copy because next line we use erase
        id_radius.erase(id_radius.begin());
        bool new_option=False;
        bool newOption=False;
        if(popped.second->topk.size()==1){
            // if it is top1, push its adjacent vertex
            const vector<int> &top1_adj=ch_obj.get_neighbor_vertex(popped.second->topk.front());
            for (int adj_opt:top1_adj) {
                if(top1_calculated.find(adj_opt)!=top1_calculated.end()){
                    continue;
                }
                top1_calculated.insert(adj_opt);
                auto iter = top1_region.find(adj_opt);
                if(iter==top1_region.end()){
                    continue; // precision problem occurs!
                }
                double dist=dist_region_w(parent_region,iter->second, w);
                if(dist!=INFINITY){
                    vector<int> tmp;
                    tmp.push_back(adj_opt);
                    vector<vector<double>> new_region(parent_region);
                    for (vector<double> &its_own_topr: iter->second) {
                        new_region.push_back(its_own_topr);
                    }
                    region *r=new region(tmp, new_region);
                    id_radius.emplace(dist, r);
                    cnt++;
                }
            }
        }
        if(popped.second->topk.size()==k){ // a region that don't need to be partitioned
            for (int opt:popped.second->topk) {
                auto iter=options.find(opt);
                if(iter==options.end()){// new option
                    options.insert(opt);
                    utk_option_ret.emplace_back(opt, popped.first);
                    cout << "radius: " << popped.first << "\n";
                    if(true){
                        vector<double> tmp(PG[opt], PG[opt]+dim);
                        cout<<options.size()<<": "<<popped.second->cone.size()<<"," << popped.first << ", "<<utk_cones_ret.size()<<
                            ", "<<popped.second->topk.back() <<"," << ch_obj.get_neighbor_vertex(popped.second->topk.back()).size()<<" # ";
                        cout<<popped.second->topk<<":"<<tmp<<"\n";
                    }
                    auto now = chrono::steady_clock::now();
                    chrono::duration<double> elapsed_seconds= now-begin;
                    cout << "time: " << options.size() << ", " << elapsed_seconds.count()<<"\n";
                    new_option=True;
                }
            }
            if(new_option){
                popped.second->setRadius(popped.first);
                utk_cones_ret.emplace_back(popped);
            }
        }
        else{
            auto now1 = chrono::steady_clock::now();

            unordered_set<int> ch_upd_s;
            int m=0;
            for(int top: popped.second->topk){
                for(int adj: ch_obj.get_neighbor_vertex(top)){
                    ch_upd_s.insert(adj);
                }
                m=max(ch_obj.get_option_layer(top), m);
            }
            for(int mp1_opt: ch_obj.get_layer(m+1)){
                ch_upd_s.insert(mp1_opt);
            }
            for (int top: popped.second->topk) {
                ch_upd_s.erase(top);
            }

            unordered_map<int, int> dominated_cnt;
            for(int i:ch_upd_s){
                dominated_cnt[i]=ch_obj.dominated_map[i].size();
            }
            for(int opt:ch_upd_s){
                auto iter=ch_obj.dominated_map.find(opt);
                if(iter!=ch_obj.dominated_map.end()){
                    for(int i:popped.second->topk){
                        if(iter->second.find(i)!=iter->second.end()){
                            --dominated_cnt[opt];
                        }
                    }
                }
            }
            vector<int> ch_upd;
            for (int i:ch_upd_s) {
                auto iter=dominated_cnt.find(i);
                if(iter!=dominated_cnt.end()&&iter->second<=0){
                    ch_upd.push_back(i);
                }
            }

//            vector<int> EMPTY;
//            vector<double> wi;
//            qp_solver2(w, popped.second->cone, 0, EMPTY, PG, wi);
//            int top_opt=ch_upd[0];
//            double max_dot=0;
//            for(int ch_upd_opt: ch_upd){
//                double dot=0;
//                for (int i = 0; i <dim ; ++i) {
//                    dot+=(PG[ch_upd_opt][i]+SIDELEN)*wi[i];
//                }
//                if(dot>max_dot){
//                    max_dot=dot;
//                    top_opt=ch_upd_opt;
//                }
//            }
//            EMPTY.push_back(top_opt);
//            vector<double> UNUSED;
//            vector<int> new_ch_upd;
//            for(int ch_upd_opt: ch_upd){
//                double ldist=qp_solver2(w, popped.second->cone, ch_upd_opt, EMPTY, PG, UNUSED);
//                if(ldist!=INFINITY){
//                    new_ch_upd.emplace_back(ch_upd_opt);
//                }
//            }
            auto now2 = chrono::steady_clock::now();
//            vector<int> tc=vector<int>(ch_upd_s.begin(), ch_upd_s.end());
            vector<int> tc=vector<int>(ch_upd.begin(), ch_upd.end());
            vector<int> new_ch_upd=r_dominate_filter(tc, popped.second->cone, PG);
//            cout<<ch_upd.size()<<"-->"<<new_ch_upd.size()<<"\n";
            auto now3 = chrono::steady_clock::now();

            const int square_vertex_cnt=dim+1;
            vector<vector<double>> square_vertexes(1, vector<double>(dim));

            vector<pair<int, double>> topip1;
            for (int ch_upd_opt: new_ch_upd) {
                vector<int> cmp;
                for(int tmp_opt: new_ch_upd){
                    if(tmp_opt!=ch_upd_opt){
                        cmp.push_back(tmp_opt);
                    }
                }
                double dist=qp_solver2(w, popped.second->cone, ch_upd_opt, cmp, PG);
                if(dist!=INFINITY){
                    topip1.emplace_back(ch_upd_opt, dist);
                }
            }
            auto now4 = chrono::steady_clock::now();

            for(auto &opt_radius:topip1){
                vector<int> tmp(popped.second->topk);
                tmp.push_back(opt_radius.first);
                vector<vector<double>> new_region(popped.second->cone);
                for (auto& otherOpt_radius: topip1) {
                    if(otherOpt_radius.first==opt_radius.first){
                        continue;
                    }
                    vector<double> its_own_topr(dim);
                    for (int i = 0; i < dim; ++i) {
                        its_own_topr[i]=PG[otherOpt_radius.first][i]-PG[opt_radius.first][i];
                    }
                    new_region.push_back(its_own_topr);
                }
                region *r=new region(tmp, new_region);
                id_radius.emplace(opt_radius.second, r);
                cnt++;
            }
        }

        if(popped.second->topk.size()!=k){
            delete(popped.second);
        }
        else if(!new_option && !newOption){
            delete(popped.second);   // to save mem
        }
    }
    cout<<"cnt: "<<cnt<<"\n";
    for (auto left:id_radius) {
        delete (left.second);
    }
    return cnt;
}
#include <lp_lib.h>

class topi_region{
    double radius;
public:
    int opt_i; // for root, this value is -1
    int top_what;
    topi_region *parent;
    vector<int> child_constraint;
    unordered_set<int> candidate ;// TODO 2021/2/1
    unordered_set<int> fixed_candidate ;
    topi_region(topi_region *pr, int topi, int rank){
        if(pr){
            fixed_candidate=pr->fixed_candidate;
        }
        parent= pr;
        opt_i=topi;
        top_what=rank;
    }
    void set_child_constr(const vector<int> &child_constr){
        child_constraint=child_constr;
    }
    void set_radius(double r){
        radius=r;
    }
    vector<int> get_topi(){
        vector<int> ret;
        topi_region *pr=this;
        while(pr->parent){
            ret.push_back(pr->opt_i);
            pr=pr->parent;
        }
        return ret;
    }
    void lpModel(lprec* lp, int& dimen)
    {
        if (lp == nullptr)
        {
            fprintf(stderr, "Unable to create new LP model\n");
            exit(0);
        }
        // constraints in each dimension
        for (int di = 0; di < dimen; di++) {
            set_bounds(lp, di+1, 0.0, 1.0);
        }
    }
    void set_candidate(const unordered_set<int> &c, const vector<vector<double>>& H1, float **PG, int dim, int k){
        vector<double> this_w(dim);
        vector<float> UNUSED(dim);
        vector<int> EMPTY;
        feasible(UNUSED, H1, 0, EMPTY, PG, this_w);
        multimap<double, int, greater<double>> find_top_ip1_to_k;
        for(int opt:c){
            double dot=0;
            for (int i = 0; i <dim ; ++i) {
                dot+=(PG[opt][i]+SIDELEN)*this_w[i];
            }
            find_top_ip1_to_k.emplace(dot, opt);
        }
        int size=k-this->top_what;
        vector<int> topi=this->get_topi();
        unordered_set<int> top_1_to_k(topi.begin(), topi.end());
        while (top_1_to_k.size()<k) {
            top_1_to_k.insert(find_top_ip1_to_k.begin()->second);
            find_top_ip1_to_k.erase(find_top_ip1_to_k.begin());
        }
        for(int i:topi){
            top_1_to_k.erase(i);
        }
        assert(top_1_to_k.size()==size);
        vector<double> UNUSED2(dim);
        for(int opt:c){
            for(int cmp:top_1_to_k){
                vector<int> tmp;
                tmp.push_back(cmp);
                double ldist = feasible(UNUSED, H1, opt, tmp, PG, UNUSED2);
                if (ldist != INFINITY) {
                    candidate.insert(opt);
                    break;
                }
            }
        }
        for(int cmp:top_1_to_k){
            candidate.insert(cmp);
        }
    }
    double feasible(const vector<float>& w, const vector<vector<double>>& H1, int opt, vector<int>& cmp, float **PG,
                    vector<double> &retv){
//        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * m
//        ce0: m
//
//        CI: n * p
//        ci0: p
//
//        x: n
        quadprogpp::Matrix<double> G, CE, CI;
        quadprogpp::Vector<double> g0, ce0, ci0, x;
        int n=w.size(), m, p;
        G.resize(n, n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if(i==j){
                    G[i][j]=1.0;
                }else{
                    G[i][j]=0.0;
                }
            }
        }
        g0.resize(n);
        {
            for (int i = 0; i < n; i++){
                g0[i]=-w[i];
            }
        }


        m = 1;
        CE.resize(n, m);
        {
            for (int i = 0; i < n; i++)
                CE[i][0]=1.0;
        }
        ce0.resize(m);
        {
            for (int j = 0; j < m; j++)
                ce0[j]=-1.0;
        }


        p = H1.size()+cmp.size()+w.size();
        vector<int> opts;
        vector<int> cmps;
        int child_opt=opt_i;
        topi_region *pr=parent;
        while(pr!= nullptr){
            p+=pr->child_constraint.size()-1;
            for(int i:pr->child_constraint){
                if(i==child_opt){
                    continue;
                }
                opts.push_back(child_opt);
                cmps.push_back(i);
            }
            child_opt=pr->opt_i;
            pr=pr->parent;
        }
        CI.resize(n, p);
        {
            for (int i = 0; i < n; i++){
                for (int j = 0; j < H1.size(); j++){
                    CI[i][j]=-H1[j][i];
                }
                for (int j = H1.size(); j < H1.size()+cmp.size(); j++){
                    CI[i][j]=PG[opt][i]-PG[cmp[j-H1.size()]][i];
                }
                for (int j = H1.size()+cmp.size(); j < p-w.size(); j++){
                    CI[i][j]=PG[opts[j-H1.size()-cmp.size()]][i]-PG[cmps[j-H1.size()-cmp.size()]][i];
                }
                for (int j = p-w.size(); j < p; j++){
                    if(i==j+w.size()-p){
                        CI[i][j]=1.0;
                    }else{
                        CI[i][j]=0.0;
                    }
                }
            }
        }
        ci0.resize(p);
        {
            for (int j = 0; j < p; j++)
                ci0[j]=0;
        }
        x.resize(n);
        double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
        double ret=0;
        for (int k = 0; k < n; ++k) {
            ret+=(x[k]-w[k])*(x[k]-w[k]);
        }
        retv.resize(n);
        for (int k = 0; k < n; ++k) {
            retv[k]=x[k];
        }
        if(solver_ret == std::numeric_limits<double>::infinity()){
            return INFINITY;
        }
        return sqrt(ret);
    }
    bool linear_feasible(const vector<float>& w, int opt, vector<int>& cmp, float **PG){
        int dim=w.size();
        lprec *lp = make_lp(0, dim);
        set_verbose(lp, IMPORTANT);
        set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
        set_add_rowmode(lp, TRUE);
        lpModel(lp, dim);
        vector<double> obj(dim+1);
        set_obj_fn(lp, obj.data());
        vector<double *> store;
        int child_opt=opt_i;
        topi_region *pr=parent;
        while(pr!= nullptr){
            for(int other:pr->child_constraint){
                if(other==child_opt){
                    continue;
                }
                double *tmp=new double[dim+1];
                for (int i = 0; i <dim ; ++i) {
                    tmp[i+1]=PG[other][i]-PG[child_opt][i];
                }
                store.push_back(tmp);
                add_constraint(lp, tmp, LE, 0.0);
            }
            child_opt=pr->opt_i;
            pr=pr->parent;
        }
        for (int other:cmp) {
            double *tmp=new double[dim+1];
            for (int i = 0; i <dim ; ++i) {
                tmp[i+1]=PG[other][i]-PG[opt][i];
            }
            store.push_back(tmp);
            add_constraint(lp, tmp, LE, 0.0);
        }
        vector<double> ONE(dim+1, 1.0);
        add_constraint(lp, ONE.data(), EQ, 1.0);
        set_add_rowmode(lp, FALSE);
        set_timeout(lp, 1);
        int ret = solve(lp);
        delete_lp(lp);
        for(double *tmp: store){
            delete [] (tmp);
        }
        return ret<=1;
    }

};


vector<double> get_minmax_dimm1_3(const vector<vector<double>> &R) {
    //        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * m
//        ce0: m
//
//        CI: n * p
//        ci0: p
//
//        x: n
    int dim=R[0].size();
    double min_max[]={1, -1};
    vector<double> ret;
    for (int i1 = 0; i1 < dim-1; ++i1) {
        for (int j1 = 0; j1 <2 ; ++j1){
            vector<double> g0(dim);
            int n=dim;
            g0.resize(n);
            for (int k = 0; k < n; ++k){
                if(i1==k){
                    g0[k]=min_max[j1];
                }else{
                    g0[k]=0;
                }
            }
            double solver_ret;
            qp_solver qs(g0, R, solver_ret);
            ret.push_back(min_max[j1]*solver_ret);
        }
    }

    return ret;
}

vector<double> get_minmax_dimm1_2(const vector<vector<double>> &R) {
    //        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * m
//        ce0: m
//
//        CI: n * p
//        ci0: p
//
//        x: n
    int dim=R[0].size();
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=dim, m, p;
    double min_max[]={1, -1};
    vector<double> ret;

    for (int i1 = 0; i1 < dim-1; ++i1) {
        for (int j1 = 0; j1 <2 ; ++j1){
            g0.resize(n);

            for (int k = 0; k < n; ++k){
                if(i1==k){
                    g0[k]=min_max[j1];
                }else{
                    g0[k]=0;
                }
            }

            G.resize(n, n);
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    if(i==j){
                        G[i][j]=1e-8;
                    }else{
                        G[i][j]=0.0;
                    }
                }
            }

            m = 1;
            CE.resize(n, m);
            {
                for (int i = 0; i < n; i++)
                    CE[i][0]=1.0;
            }
            ce0.resize(m);
            {
                for (int j = 0; j < m; j++)
                    ce0[j]=-1.0;
            }


            p = R.size()+dim;
//            p = R.size();

            CI.resize(n, p);
            {
                for (int i = 0; i < n; i++){
                    for (int j = 0; j < R.size(); j++){
                        CI[i][j]=-R[j][i];
                    }
                    for (int j = R.size(); j < p; j++){
                        if(i==j-R.size()){
                            CI[i][j]=1.0;
                        }else{
                            CI[i][j]=0.0;
                        }
                    }
                }
            }
            ci0.resize(p);
            {
                for (int j = 0; j < p; j++)
                    ci0[j]=0;
            }
//            ci0.resize(R.size());
//            {
//                for (int j = 0; j < R.size(); j++)
//                    ci0[j]=0;
//            }
            x.resize(n);

            double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
            if(solver_ret==std::numeric_limits<double>::infinity()){
                ret.clear();
                return ret;
            }
            ret.push_back(min_max[j1]*solver_ret);
            vector<double> tmp(n);
            for (int k = 0; k < n; ++k) {
                tmp[k]=x[k];
            }
//            cout<<solver_ret<<":"<<tmp<<endl;
        }
    }
    vector<double> tmp(n);
    for (int k = 0; k < n; ++k) {
        tmp[k]=x[k];
    }
    return ret;
}

vector<double> get_minmax_dimm1(const vector<vector<double>> &R){
    vector<double> ret;
    if(!R.empty()){
        int dim=R[0].size();

        vector<double *> store;

        vector<double> ONE(dim+1, 1.0);

        double min_max[]={-1, 1};
        for (int i = 0; i < dim-1; ++i) {
            for (int j = 0; j <2 ; ++j) {
                vector<double> w(1+R.size()+1+dim);
                lprec *lp = make_lp(0, dim);
                lpModel(lp, dim);

                set_verbose(lp, IMPORTANT);
//                set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE );

                vector<double> obj(dim+1);
                obj[i+1]=min_max[j];
                set_obj_fn(lp, obj.data());
                set_add_rowmode(lp, TRUE);
                add_constraint(lp, ONE.data(), EQ, 1.0);
                for(const vector<double>& a_constrain: R){
                    double *tmp=new double[dim+1];
                    store.push_back(tmp);
                    for (int k = 0; k <dim ; ++k) {
                        tmp[k+1]=a_constrain[k];
                    }
                    add_constraint(lp, tmp, LE, 0.0);
                }
                set_add_rowmode(lp, FALSE);
                set_maxim(lp);
                set_timeout(lp, 1);
                solve(lp);
                get_primal_solution(lp, w.data());
                ret.push_back(w[1+R.size()+1+i]);
                delete_lp(lp);
            }
        }
        for(double *tmp: store){
            delete [] (tmp);
        }
    }
    return ret;
}

void rec_get_bound(const vector<double> &minmax,  int cur, vector<double> &v, vector<vector<double>> &ret){
    if(cur<minmax.size()){
        vector<double> v1=v;
        v1.push_back(minmax[cur]);
        vector<double> v2=v;
        v2.push_back(minmax[cur+1]);
        rec_get_bound(minmax, cur+2, v1, ret);
        rec_get_bound(minmax, cur+2, v2, ret);
    }else{
        double s=sum(v.begin(), v.end());
        v.push_back(1.0-s);
        ret.push_back(v);
    }
}

vector<vector<double>> get_region_bound(const vector<vector<double>> &R){
    auto now1 = chrono::steady_clock::now();
    vector<double> minmax = get_minmax_dimm1_2(R); // a 2*(dim-1) vector

    auto now2 = chrono::steady_clock::now();
    vector<double> EMPTY;
    vector<vector<double>> bounds;
    if(minmax.empty()){
        return bounds;
    }
    rec_get_bound(minmax, 0, EMPTY, bounds);
    auto now3 = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds1= now2-now1;
    chrono::duration<double> elapsed_seconds2= now3-now2;


//            cout << "t:"
//                    << elapsed_seconds1.count()<< ","
//                    << elapsed_seconds2.count()<< ","<<"\n";
    return bounds;
}


bool r_dominate(const vector<vector<double>>& ws, int opt1, int opt2, float **PG, int dim){
    double dot1, dot2;
    for(auto &w:ws){
        dot1=0;
        dot2=0;
        for (int i = 0; i < dim; ++i) {
            dot1+=PG[opt1][i]*w[i];
        }
        for (int i = 0; i < dim; ++i) {
            dot2+=PG[opt2][i]*w[i];
        }
        if(dot1<dot2){
            return false;
        }
    }
    return true;
}

vector<int> r_dominate_filter(const vector<int> &top_ip1_cdd,
                              const vector<vector<double>> &R, float **PG){
    // 1. use LP solver get min max of dim-1 dimension
    // 2. generate the 2^(d-1) points, marked as S
    // 3. for wi in S:
    // 4.     find top-k of wi, add them to top k set
    // 5. if top k set size is k:
    // 6.     prune this region, directly add top-k to result option set
    // 7. else:
    // 8.     r-dominate test, return the non-r-dominated options

    // r-dominated test
    // a. if it is r-dominated by over or equals to k-i options, prune it for this region
    // b. if it is r-domianted by over or equals to 1 option, prune it for current candidate of top-(i+1)

    // input
    // \para R: current region
    // \para top_ip1_cdd: top-(i+1) candidate
    vector<int> ret;
    if(R.empty()){
        return ret;
    }
    auto now1 = chrono::steady_clock::now();
    vector<vector<double>> ws=get_region_bound(R);
    if(ws.empty()){
        return ret;
    }
    auto now2 = chrono::steady_clock::now();

    int dim=R[0].size();
    for(int opt:top_ip1_cdd){
        bool flag=true;
        for(int cmp:top_ip1_cdd){
            if(opt!=cmp){
                if(r_dominate(ws, cmp, opt, PG, dim)){
                    flag=false;
                    break;
                }
            }
        }
        if(flag){
            ret.push_back(opt);
        }
    }
    return ret;
}

int topRegions_efficient3(vector<vector<double>> &parent_region, ch &ch_obj,
                          multimap<double, topi_region*> &id_radius, float **PG, int dim, int X,
                          const int k, vector<float> &w, unordered_set<int> &top1_calculated,
                          vector<pair<int, double>> &utk_option_ret,
                          vector<pair<double, topi_region*>> &utk_cones_ret,
                          unordered_map<int, vector<int>> &top1_region){
    // init "top1_calculated" with top1 respecting to w
    auto begin=chrono::steady_clock::now();
    unordered_set<int> options;
    for (int i:top1_calculated) {
        options.insert(i);
        utk_option_ret.emplace_back(i, 0);
        cout << "radius: " << 0 << "\n";
    }
    vector<topi_region *> all;
    int cnt=0;
    int rcnt=0;
    cout<<parent_region.size()<<endl;
    while(options.size() < X && !id_radius.empty()){
        ++rcnt;
        pair<double, topi_region*> popped=*id_radius.begin(); // must deep copy because next line we use erase
        id_radius.erase(id_radius.begin());
        bool new_option=False;
        bool newOption=False;
        if(popped.second->top_what==1){
            // if it is top1, push its adjacent vertex
            const vector<int> &top1_adj=ch_obj.get_neighbor_vertex(popped.second->opt_i);
            for (int adj_opt:top1_adj) {
                if(top1_calculated.find(adj_opt)!=top1_calculated.end()){
                    continue;
                }
                top1_calculated.insert(adj_opt);
                auto iter = top1_region.find(adj_opt);
                if(iter==top1_region.end()){
                    continue; // precision problem occurs!
                }
                double dist=qp_solver2(w, parent_region,adj_opt, iter->second, PG);
                if(dist!=INFINITY){
                    topi_region *r=new topi_region(popped.second->parent, adj_opt, 1);
                    all.push_back(r); // to delete
                    id_radius.emplace(dist, r);
                    cnt++;
                }
            }
        }
        vector<int> topk=popped.second->get_topi();
        for(int opt:topk){
            auto iter=options.find(opt);
            if(iter==options.end()) {// new option
                options.insert(opt);
                utk_option_ret.emplace_back(opt, popped.first);
                cout << "radius: " << popped.first << "\n";
                if(true){
                    vector<double> tmp(PG[opt], PG[opt]+dim);
                    cout<<options.size()<<": " << popped.first << ", "<<utk_cones_ret.size()<<
                        ", "<<rcnt <<"," << ch_obj.get_neighbor_vertex(opt).size()<<" # "<< popped.second->get_topi()<<"\n";
                }
                auto now = chrono::steady_clock::now();
                chrono::duration<double> elapsed_seconds= now-begin;
                cout << "time: " << options.size() << ", " << elapsed_seconds.count()<<"\n";
                new_option=True;
            }
        }
        if(new_option){
            popped.second->set_radius(popped.first);
        }
        if(popped.second->top_what==k){ // a region that don't need to be partitioned
            vector<int> topk=popped.second->get_topi();
            for(int opt:topk){
                auto iter=options.find(opt);
                if(iter==options.end()) {// new option
                    options.insert(opt);
                    utk_option_ret.emplace_back(opt, popped.first);
                    cout << "radius: " << popped.first << "\n";
                    if(true){
                        vector<double> tmp(PG[opt], PG[opt]+dim);
                        cout<<options.size()<<": " << popped.first << ", "<<utk_cones_ret.size()<<
                            ", "<<rcnt <<"," << ch_obj.get_neighbor_vertex(opt).size()<<" # "<< popped.second->get_topi()<<"\n";
                    }
                    auto now = chrono::steady_clock::now();
                    chrono::duration<double> elapsed_seconds= now-begin;
                    cout << "time: " << options.size() << ", " << elapsed_seconds.count()<<"\n";
                    new_option=True;
                }
            }
            if(new_option){
                popped.second->set_radius(popped.first);
            }
        }
        else{
            vector<int> topi=popped.second->get_topi();
            unordered_set<int> ch_upd_s;
            int m=0; // find maximal convex hull layer
            for(int top: topi){
                for(int adj: ch_obj.get_neighbor_vertex(top)){
                    ch_upd_s.insert(adj);
                }
                m=max(ch_obj.get_option_layer(top), m);
            }
            for(int mp1_opt: ch_obj.get_layer(m+1)){
                ch_upd_s.insert(mp1_opt);
            }
            for (int top: topi) {
                ch_upd_s.erase(top);
            }
            // TODO 这里其实可以直接测dominate
            unordered_map<int, int> dominated_cnt;// remove options dominated by other options than topi
            for(int i:ch_upd_s){
                dominated_cnt[i]=ch_obj.dominated_map[i].size();
            }
            for(int opt:ch_upd_s){
                auto iter=ch_obj.dominated_map.find(opt);
                if(iter!=ch_obj.dominated_map.end()){
                    for(int i:topi){
                        if(iter->second.find(i)!=iter->second.end()){
                            --dominated_cnt[opt];
                        }
                    }
                }
            }
            vector<int> ch_upd;
            for (int i:ch_upd_s) {
                auto iter=dominated_cnt.find(i);
                if(iter!=dominated_cnt.end()&&iter->second<=0){
                    ch_upd.push_back(i);
                }
            }

            popped.second->set_candidate(popped.second->parent->candidate, parent_region, PG, dim, k);
            vector<int> new_ch_upd;
            for(int ch_upd_opt: ch_upd){
                if(popped.second->candidate.find(ch_upd_opt)!=popped.second->candidate.end()){
                    new_ch_upd.push_back(ch_upd_opt);
                }
            }
//            cout<<popped.second->parent->candidate.size()<<"\n";
//            cout<<ch_upd.size()<<","<<popped.second->candidate.size()<<"\n";


            const int square_vertex_cnt=dim+1;
            vector<vector<double>> square_vertexes(1, vector<double>(dim));
            // 时间主要开销， 主要性能瓶颈, 求凸包以及每个凸包点的top region

            vector<pair<int, double>> topip1;
            vector<int> next_topi;
            vector<double> UNUSED(dim);

            for (int ch_upd_opt: new_ch_upd) {
                vector<int> cmp;
                for(int tmp_opt: new_ch_upd){
                    if(tmp_opt!=ch_upd_opt){
                        cmp.push_back(tmp_opt);
                    }
                }
                double dist=popped.second->feasible(w, parent_region, ch_upd_opt, cmp, PG, UNUSED);
                if(dist!=INFINITY){
                    topip1.emplace_back(ch_upd_opt, dist);
                    next_topi.push_back(ch_upd_opt);
                }
            }
//            cout<<"original size: "<<ch_upd.size()<<"\n";
//            cout<<"midterm size: "<<new_ch_upd.size()<<"\n";
//            cout<<"result size: "<<next_topi.size()<<"\n";

            for(auto &opt_radius:topip1){
                topi_region *r=new topi_region(popped.second, opt_radius.first, popped.second->top_what+1);
                all.push_back(r);
                id_radius.emplace(opt_radius.second, r);
                cnt++;
            }
            popped.second->child_constraint=next_topi;
        }
    }
    cout<<"cnt: "<<cnt<<"\n";
    for (auto left:all) {
        delete (left);
    }
    return cnt;
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
//    square_vertexes.emplace_back(dim, 1.0);

    string s = to_string(dim) + " " + to_string(points.size() + square_vertexes.size()) + " ";
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
    Qhull q(rbox, "QJ");
    qhull_user qu;
    return qu.get_cone_norms(q, points); // make sure the first of Qhull input is \vec{0}_{d}
}

void utk_basic(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
               vector<pair<int, double>> &utk_option_ret,
               vector<pair<vector<int>, vector<vector<double>>>> &utk_cones_ret){

    // 2 return value
    // 1. array of <option, topk_radius>
    // 2. array of <topk, region>
    // "apply k=1"
//    test_build_qhull();
//    return;
    auto begin = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, PointSet);
    pair<int, float> next={-1, INFINITY};
    //fetch top X options
    cout<< "begin fetch top X"<<endl;
    vector<int> CH_1_X_opt;
    while(CH_1_X_opt.size() < X){  // fetch top X
        next=get_next_obj.get_next();
        CH_1_X_opt.push_back(next.first);
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;

    cout<< elapsed_seconds.count() << " fetch top X finish\n";
    // qhull class in lib qhull
    const int square_vertex_cnt=dim+1;
    vector<vector<double>> square_vertexes(square_vertex_cnt, vector<double>(dim));

    // init qhull with top X options
    CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;

    cout<< elapsed_seconds.count() << " first time build qhull finish\n";

    // a 3-d example of square_vertex, square_vertex_cnt=4
    // point 0: (max(points[:, 0]), 0, 0)
    // point 1: (0, max(points[:, 1]), 0)
    // point 2: (0, 0, max(points[:, 1]))
    // point 3: \vec{0}

    // rho_star is computed as when CH_1.size()=X
    int cnt=0;
    while(CH_1_X_opt.size() < X){ // while(CH_1.size<X)
        while(CH_1_X_opt.size() < X){
            next=get_next_obj.get_next();
            if(next.second==INFINITY){
                break;
            }
            update_square_vertexes(square_vertexes, PointSet[next.first], dim);
            CH_1_X_opt.push_back(next.first);
            cnt++;
        }
        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
    }
    // for now, qhull_obj contains points 0~3 and convex hull vertexes
    float rho_star=next.second;

    cout << "init rho_star: "<<rho_star<<endl;
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish rho_star compute\n";
    // use known X version code to fetch rskyband options,
    // bear in mind such that we init \rho as \rho_star and X as INFINITY
    vector<pair<long int, float>> interval;
    computeRho(dim, k, INFINITY, w, *rtree, PointSet, interval, rho_star);
    vector<int> rskyband_CS;
    for (pair<long int, float> &p:interval) {
        rskyband_CS.push_back(p.first);
    }
    cout<< "rskyband size: "<<rskyband_CS.size()<< "\n";
    for (int i = 0; i < rskyband_CS.size(); ++i) {
        for (int j = 0; j < dim; ++j) {
            cout<< PointSet[rskyband_CS[i]][j] << ",";
        }
        cout << "\n";
    }
    ch ch_obj(rskyband_CS, PointSet, dim);
    const vector<int>& top1_idxes=ch_obj.get_layer(1);
    vector<vector<double>> tmp;

    double rho_star_d=rho_star; // change from float to double
    for (vector<c_float> &e:g_r_domain_vec) {
        for (int i = 0; i < dim; ++i) {
            cout << e[i] <<", ";
        }
        cout <<"\n";
    }
    cout <<"\n";

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    cout<< " begin generate domain" << endl;
//    vector<vector<double>> begin_region=points_to_halfspace(tmp);
    vector<vector<double>> begin_region;
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    for (int i = 0; i < begin_region.size(); ++i) {
        for (int j = 0; j < dim; ++j){
            cout << begin_region[i][j] << ", ";
        }
        cout<<"\n";
    }
    cout <<"\n";

    cout<< " end generate domain" << endl;
    multimap<double, region*> id_radius; // <radius, region>
    vector<int> init_topi;
    vector<int> init_neighbors;
    cout<< "starting recursively get top regions\n";
    topRegions(begin_region, CH_1_X_opt, ch_obj, id_radius, 0,  PointSet, dim,
               k, init_topi, w, init_neighbors);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    cout<< " finish recursively get top regions\n";

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

int utk_efficient3(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                   vector<pair<int, double>> &utk_option_ret,
                   vector<pair<double, region*>> &utk_cones_ret){

    // two return values
    // 1. array of <option, topk_radius>
    // 2. array of <topk, region>
    auto begin = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;

    // 1. fetch 1-skyband with m options
    // 2. continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimated \rho^* from step 2
    // 4. apply \rho^* to get r-k-skyband (this is for pruning irrelevent options)
    // 5. apply \rho^* to get initial domain (this is for pruning unnecessary regions)
    // 6. get the top regions of CH1's options
    // 7. init ORU heap with top1's top region
    // 8. apply ORU, that is recursively divide regions


    // 1. begin: fetch 1-skyband with m options
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, PointSet);
    pair<int, float> next={-1, INFINITY};
    cout<< "begin fetch CH1"<<endl;
    vector<int> CH_1_X_opt;
    int top1=0;
    bool top1f= false;
    while(CH_1_X_opt.size() < X){
        next=get_next_obj.get_next();
        cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
        CH_1_X_opt.push_back(next.first);
        if(!top1f){
            top1=CH_1_X_opt.back();
            top1f=true;
        }
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " fetch top-m finish\n";
    // 1. end: fetch 1-skyband with m options

    // a 3-d example of square_vertex, square_vertex_cnt=4
    // point 0: (max(points[:, 0]), 0, 0)
    // point 1: (0, max(points[:, 1]), 0)
    // point 2: (0, 0, max(points[:, 2]))
    // point 3: \vec{0}
    const int square_vertex_cnt=dim+1;
    vector<vector<double>> square_vertexes(1, vector<double>(dim));

    // qhull class in lib qhull
    // init qhull with top X options
    CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " first time build qhull finish\n";

    // rho_star is computed when CH_1.size()=X
    int cnt=0;
    // 2. begin: continue step 1 to fill CH1 until CH1 contains m options

    int minc=X, maxc=2*X, midc=X;
    while(CH_1_X_opt.size() < X){ // find the maxc
        midc=maxc; // make sure correct rho_star
        if(get_next_obj.interval.size()<maxc){
            while(get_next_obj.interval.size()<maxc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <maxc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();

        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()>=X){
            break;
        }else{
            minc=maxc;
            maxc*=2;
        }
    }
    while(CH_1_X_opt.size() != X){ // while(CH_1.size<X)
        midc=(maxc+minc)/2;
        if(get_next_obj.interval.size()<=midc){
            while(get_next_obj.interval.size()<midc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                update_square_vertexes(square_vertexes, PointSet[next.first], dim);
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <midc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();

        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()==X){
            break;
        }else if(CH_1_X_opt.size()>X){
            maxc=midc-1;
        }else{
            minc=midc+1;
        }
    }

    // 2. end: continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimate \rho^* from step 2
    float rho_star=get_next_obj.interval[min(midc-1, get_next_obj.interval.size()-1)].second;

    cout << "init rho_star: "<<rho_star<<endl;
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish rho_star compute\n";
    // use known X version code to fetch rskyband options,
    // bear in mind such that we init \rho as \rho_star and X as INFINITY
    // 4. begin: apply \rho^* to get r-k-skyband
    vector<pair<long int, float>> interval;
    computeRho(dim, k, INFINITY, w, *rtree, PointSet, interval, rho_star);
    vector<int> rskyband_CS;
    for (pair<long int, float> &p:interval) {
        rskyband_CS.push_back(p.first);
    }
    // 4. end: apply \rho^* to get r-k-skyband
    cout<< elapsed_seconds.count() << " rskyband size: "<<rskyband_CS.size()<< "\n";

    vector<vector<double>> tmp;
    // 5. begin: apply \rho^* to get initial domain

//    double rho_star_d=rho_star; // change from float to double
//
//    for (vector<c_float> &e:g_r_domain_vec) {
//        double atc_rho=rho_star_d;
////        double atc_rho=rho_star_d*sqrt(dim);
//
//        for (int i = 0; i < e.size(); ++i) {
//            if(e[i]<0){
//                atc_rho=min(atc_rho, -w[i]/e[i]); // in case of w[i] + \rho * e[i] <0 or >1
//            }
//        }
//        tmp.push_back(atc_rho*e+w);
//    }
//    now = chrono::steady_clock::now();
//    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " begin generate domain" << endl;
//    vector<vector<double>> begin_region=points_to_halfspace(tmp);
    vector<vector<double>> begin_region;

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() <<"finish generate domain \n";
    // 5. end: apply \rho^* to get initial domain

    // 6. begin: get the top regions of CH1's options
    cout << " begin find top1 region" << endl;
    vector<vector<double>> square_vertexes2(square_vertex_cnt, vector<double>(dim));
    // init convex hull object, this object can conveniently return CH_i given input i
    ch ch_obj(rskyband_CS, PointSet, dim);
    const vector<int>& top1_idxes=ch_obj.get_layer(1);// get CH1
    cout<<"finish finding CH1"<<endl;
    unordered_map<int, vector<int>> top1_region;
    for(int i:top1_idxes){
        auto iter=ch_obj.A_p.find(i);
        if(iter!=ch_obj.A_p.end()){
            top1_region[i]=iter->second;
        }
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish find top1 region" << endl;
    // 6. begin: get the top regions of CH1's options

    // 7. begin: init ORU heap with top1's top region
    auto iter=top1_region.find(top1);
    assert(iter != top1_region.end());
    unordered_set<int> top1_calculated;
    top1_calculated.insert(top1);
    topi_region *root=new topi_region(nullptr, -1, -1);
    root->set_child_constr(top1_idxes);
    root->candidate=ch_obj.all;

    topi_region *top1r=new topi_region(root, top1, 1);
    top1r->set_radius(0);

    multimap<double, topi_region*> id_radius; // <radius, region>
    id_radius.emplace(0, top1r);
    // 7. end: init ORU heap with top1's top region

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << "starting recursively get top regions\n";
    // 8. begin: apply ORU, that is recursively divide regions
    vector<pair<double, topi_region*>> utk_cones_ret_;
    int ret=topRegions_efficient3(begin_region, ch_obj, id_radius,  PointSet, dim, X,
                                  k,  w, top1_calculated, utk_option_ret, utk_cones_ret_, top1_region);
    delete(root);
    delete(top1r);
    // 8. end: apply ORU, that is recursively divide regions
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    cout<< " finish recursively get top regions\n";
    return ret;
}


int utk_efficient(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                  vector<pair<int, double>> &utk_option_ret,
                  vector<pair<double, region*>> &utk_cones_ret){

    // two return values
    // 1. array of <option, topk_radius>
    // 2. array of <topk, region>
    auto begin = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;

    // 1. fetch 1-skyband with m options
    // 2. continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimated \rho^* from step 2
    // 4. apply \rho^* to get r-k-skyband (this is for pruning irrelevent options)
    // 5. apply \rho^* to get initial domain (this is for pruning unnecessary regions)
    // 6. get the top regions of CH1's options
    // 7. init ORU heap with top1's top region
    // 8. apply ORU, that is recursively divide regions


    // 1. begin: fetch 1-skyband with m options
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, PointSet);
    pair<int, float> next={-1, INFINITY};
    cout<< "begin fetch CH1"<<endl;
    vector<int> CH_1_X_opt;
    int top1=0;
    bool top1f= false;
    while(CH_1_X_opt.size() < X){
        next=get_next_obj.get_next();
        cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
        CH_1_X_opt.push_back(next.first);
        if(!top1f){
            top1=CH_1_X_opt.back();
            top1f=true;
        }
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " fetch top-m finish\n";
    // 1. end: fetch 1-skyband with m options

    // a 3-d example of square_vertex, square_vertex_cnt=4
    // point 0: (max(points[:, 0]), 0, 0)
    // point 1: (0, max(points[:, 1]), 0)
    // point 2: (0, 0, max(points[:, 2]))
    // point 3: \vec{0}
    const int square_vertex_cnt=dim+1;
    vector<vector<double>> square_vertexes(1, vector<double>(dim));

    // qhull class in lib qhull
    // init qhull with top X options
    CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " first time build qhull finish\n";

    // rho_star is computed when CH_1.size()=X
    int cnt=0;
    // 2. begin: continue step 1 to fill CH1 until CH1 contains m options

    int minc=X, maxc=2*X, midc=X;
    while(CH_1_X_opt.size() < X){ // find the maxc
        midc=maxc; // make sure correct rho_star
        if(get_next_obj.interval.size()<maxc){
            while(get_next_obj.interval.size()<maxc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <maxc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();
        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()>=X){
            break;
        }else{
            minc=maxc;
            maxc*=2;
        }
    }
    while(CH_1_X_opt.size() != X){ // while(CH_1.size<X)
        midc=(maxc+minc)/2;
        if(get_next_obj.interval.size()<=midc){
            while(get_next_obj.interval.size()<midc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                update_square_vertexes(square_vertexes, PointSet[next.first], dim);
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <midc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();
        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()==X){
            break;
        }else if(CH_1_X_opt.size()>X){
            maxc=midc-1;
        }else{
            minc=midc+1;
        }
    }

    // 2. end: continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimate \rho^* from step 2
    float rho_star=get_next_obj.interval[min(midc-1, get_next_obj.interval.size()-1)].second;

    cout << "init rho_star: "<<rho_star<<endl;
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish rho_star compute\n";
    // use known X version code to fetch rskyband options,
    // bear in mind such that we init \rho as \rho_star and X as INFINITY
    // 4. begin: apply \rho^* to get r-k-skyband
    vector<pair<long int, float>> interval;
    computeRho(dim, k, INFINITY, w, *rtree, PointSet, interval, rho_star);
    vector<int> rskyband_CS;
    for (pair<long int, float> &p:interval) {
        rskyband_CS.push_back(p.first);
    }
    // 4. end: apply \rho^* to get r-k-skyband
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " rskyband size: "<<rskyband_CS.size()<< "\n";
    vector<vector<double>> tmp;
    // 5. begin: apply \rho^* to get initial domain
    cout<< elapsed_seconds.count() << " begin generate domain" << endl;
//    vector<vector<double>> begin_region=points_to_halfspace(tmp);
    vector<vector<double>> begin_region;

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() <<"finish generate domain \n";
    // 5. end: apply \rho^* to get initial domain

    // 6. begin: get the top regions of CH1's options
    cout << " begin find top1 region" << endl;
    multimap<double, region*> id_radius; // <radius, region>
    vector<vector<double>> square_vertexes2(square_vertex_cnt, vector<double>(dim));
    // init convex hull object, this object can conveniently return CH_i given input i
    ch ch_obj(rskyband_CS, PointSet, dim);
    const vector<int>& top1_idxes=ch_obj.get_layer(1);// get CH1
    cout<<"finish finding CH1"<<endl;
    unordered_map<int, vector<vector<double>>> top1_region;
    for(int i:top1_idxes){
        auto iter=ch_obj.A_p.find(i);
        if(iter!=ch_obj.A_p.end()){
            top1_region[i]=vector<vector<double>>();
            for (int nei:iter->second) {
                vector<double> tmp(PointSet[nei], PointSet[nei]+dim);
                for (int j = 0; j <dim ; ++j) {
                    tmp[j]-=PointSet[i][j];
                }
                top1_region[i].push_back(tmp);
            }
        }
    }
//    unordered_map<int, vector<vector<double>>> top1_region=top_region(top1_idxes, PointSet, square_vertexes2);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish find top1 region" << endl;
    // 6. begin: get the top regions of CH1's options

    // 7. begin: init ORU heap with top1's top region
    auto iter=top1_region.find(top1);
    assert(iter != top1_region.end());
    vector<int> topi;
    topi.push_back(top1);
    unordered_set<int> top1_calculated;
    top1_calculated.insert(top1);
    region *r=new region(topi, iter->second);
    id_radius.emplace(dist_region_w(iter->second, w), r);
    // 7. end: init ORU heap with top1's top region

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << "starting recursively get top regions\n";
    // 8. begin: apply ORU, that is recursively divide regions
    int ret=topRegions_efficient2(begin_region, ch_obj, id_radius,  PointSet, dim, X,
                                  k,  w, top1_calculated, utk_option_ret, utk_cones_ret, top1_region);
    // 8. end: apply ORU, that is recursively divide regions
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    cout<< " finish recursively get top regions\n";
    return ret;
}


int utk_efficient_cs3(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                      vector<pair<int, double>> &utk_option_ret,
                      vector<pair<double, region*>> &utk_cones_ret, double &rho_star){

    // two return values
    // 1. array of <option, topk_radius>
    // 2. array of <topk, region>
    auto begin = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;

    // 1. fetch 1-skyband with m options
    // 2. continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimated \rho^* from step 2
    // 4. apply \rho^* to get r-k-skyband (this is for pruning irrelevent options)
    // 5. apply \rho^* to get initial domain (this is for pruning unnecessary regions)
    // 6. get the top regions of CH1's options
    // 7. init ORU heap with top1's top region
    // 8. apply ORU, that is recursively divide regions


    // 1. begin: fetch 1-skyband with m options
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, PointSet);
    pair<int, float> next={-1, INFINITY};
    cout<< "begin fetch CH1"<<endl;
    vector<int> CH_1_X_opt;
    int top1=0;
    bool top1f= false;
    while(CH_1_X_opt.size() < X){
        next=get_next_obj.get_next();
        cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
        CH_1_X_opt.push_back(next.first);
        if(!top1f){
            top1=CH_1_X_opt.back();
            top1f=true;
        }
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " fetch top-m finish\n";
    // 1. end: fetch 1-skyband with m options

    // a 3-d example of square_vertex, square_vertex_cnt=4
    // point 0: (max(points[:, 0]), 0, 0)
    // point 1: (0, max(points[:, 1]), 0)
    // point 2: (0, 0, max(points[:, 2]))
    // point 3: \vec{0}
    const int square_vertex_cnt=dim+1;
    vector<vector<double>> square_vertexes(1, vector<double>(dim));

    // qhull class in lib qhull
    // init qhull with top X options
    CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " first time build qhull finish\n";

    // rho_star is computed when CH_1.size()=X
    int cnt=0;
    // 2. begin: continue step 1 to fill CH1 until CH1 contains m options

    int minc=X, maxc=2*X, midc=X;
    while(CH_1_X_opt.size() < X){ // find the maxc
        midc=maxc; // make sure correct rho_star
        if(get_next_obj.interval.size()<maxc){
            while(get_next_obj.interval.size()<maxc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <maxc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();

        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()>=X){
            break;
        }else{
            minc=maxc;
            maxc*=2;
        }
    }
    while(CH_1_X_opt.size() != X){ // while(CH_1.size<X)
        midc=(maxc+minc)/2;
        if(get_next_obj.interval.size()<=midc){
            while(get_next_obj.interval.size()<midc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                update_square_vertexes(square_vertexes, PointSet[next.first], dim);
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <midc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();

        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()==X){
            break;
        }else if(CH_1_X_opt.size()>X){
            maxc=midc-1;
        }else{
            minc=midc+1;
        }
    }

    // 2. end: continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimate \rho^* from step 2
    rho_star=get_next_obj.interval[min(midc-1, get_next_obj.interval.size()-1)].second;

    return 0;
}


vector<int> kskyband(float **PG, int cnt, int dim, int k){
    vector<long> idx(cnt);
    vector<int> ret;
    iota(idx.begin(), idx.end(), 1);
    for (int i = 1; i <=cnt ; ++i) {
//        bool dominatedByK(const int dimen, const float pt[], vector<long> &kskyband, float* PG[], int k);
        if(!dominatedByK_noSL(dim, PG[i], idx, PG, k+1)){
            ret.push_back(i);
        }
    }
    return ret;
}

int utk_efficient_anti(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                       vector<pair<int, double>> &utk_option_ret,
                       vector<pair<double, region*>> &utk_cones_ret){

    // A FASTER VERSION FOR anti data
    // 2 return value
    // 1. array of <option, topk_radius>
    // 2. array of <topk, region>
    // "apply k=1"
    auto begin = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    //fetch top X options
    cout<< "begin fetch top X"<<endl;
    assert(X>0);
    // qhull class in lib qhull
    // use known X version code to fetch rskyband options,
    // bear in mind such that we init \rho as \rho_star and X as INFINITY
    vector<int> rskyband_CS=kskyband(PointSet, objCnt, dim, k);
    cout<< "finish kskyband"<<endl;
    int top1=1;
    double dmax=0;
    for (int i = 1; i < objCnt; ++i) {
        double dot=0;
        for (int j = 0; j < dim; ++j) {
            dot+=(PointSet[i][j]+SIDELEN)*w[j];
        }
        if(dot>dmax){
            dmax=dot;
            top1=i;
        }
    }
    cout<< elapsed_seconds.count() << " rskyband size: "<<rskyband_CS.size()<< "\n";
    ch ch_obj(rskyband_CS, PointSet, dim);
    const vector<int>& top1_idxes=ch_obj.get_layer(1);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " begin generate domain" << endl;
    vector<vector<double>> begin_region;
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() <<"finish generate domain \n";
    cout << " begin find top1 region" << endl;
    multimap<double, region*> id_radius; // <radius, region>
    vector<vector<double>> square_vertexes2(dim+1, vector<double>(dim));
    cout<< "CH_1 size:" << top1_idxes.size()<<endl;

    unordered_map<int, vector<vector<double>>> top1_region=top_region(top1_idxes, PointSet, square_vertexes2);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish find top1 region" << endl;
    auto iter=top1_region.find(top1);
    assert(iter != top1_region.end());
    vector<int> topi;
    topi.push_back(top1);
    unordered_set<int> top1_calculated;
    top1_calculated.insert(top1);
    region *r=new region(topi, iter->second);
    id_radius.emplace(dist_region_w(iter->second, w), r);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << "starting recursively get top regions\n";

    int ret=topRegions_efficient(begin_region, ch_obj, id_radius,  PointSet, dim, X,
                                 k,  w, top1_calculated, utk_option_ret, utk_cones_ret, top1_region);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    cout<< " finish recursively get top regions\n";
    return ret;
}

double ord_m_t=0;
double rho_star_t=0;
double rskyband_t=0;
double ch_1_t=0;
double rec_oru_t=0;

double non_order_feasible(vector<int> &gt, vector<int> &le, float ** PG, vector<float> &w){
    //        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * m
//        ce0: m
//
//        CI: n * p
//        ci0: p
//
//        x: n
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=w.size(), m, p;
    G.resize(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if(i==j){
                G[i][j]=1.0;
            }else{
                G[i][j]=0.0;
            }
        }
    }
    g0.resize(n);
    {
        for (int i = 0; i < n; i++){
            g0[i]=-w[i];
        }
    }

    m = 1;
    CE.resize(n, m);
    {
        for (int i = 0; i < n; i++)
            CE[i][0]=1.0;
    }
    ce0.resize(m);
    {
        for (int j = 0; j < m; j++)
            ce0[j]=-1.0;
    }


    p = gt.size()*le.size()+w.size();
    CI.resize(n, p);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < p-w.size(); j++){
                int gtid=j/le.size();
                int leid=j%le.size();
                CI[i][j]=PG[gt[gtid]][i]-PG[le[leid]][i];
            }
            for (int j = p-w.size(); j < p; j++){
                if(i==j+w.size()-p){
                    CI[i][j]=1.0;
                }else{
                    CI[i][j]=0.0;
                }
            }
        }
    }
    ci0.resize(p);
    {
        for (int j = 0; j < p; j++)
            ci0[j]=0;
    }
    x.resize(n);
    double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

    if(solver_ret == std::numeric_limits<double>::infinity()){
        return INFINITY;
    }else{
        double ret=0;
        for (int k = 0; k < n; ++k) {
            ret+=(x[k]-w[k])*(x[k]-w[k]);
        }
        return sqrt(ret);
    }
}

vector<pair<int, double>> non_order_sensitive_ORU(int dim, float **PointSet, vector<float> &w, int X, int k, ch &ch_obj){
    // 1. get top-k
    // 2. cur_region={R made by top-k(non-order sensitive region)}
    // 3. while not heap.empty():
    // 4.     R=popped region
    // 5.     add all top-k in R to RET
    // 6.     cdd_replace=get all adj options of R's top-k
    // 7.     erase all i in cdd_replace such that already in RET
    // 8.     for i in cdd_replace:
    // 9.         for j in top-k:
    // 10.            replace j witj i in top-k, test feasible, if feasible push into heap

    // step 1 get top-k
    vector<int> direct_topk=computeTopK(dim, PointSet, ch_obj.rskyband, w, k);
    bool flag=False;
    while (!flag){
        flag=True;
        for(int i:direct_topk){
            if(ch_obj.get_option_layer(i)==-1){
                flag=False;
                continue;
            }
        }
        ch_obj.get_next_layer();
    }
    cout<<"finish build CH k"<<endl;

    // step 2
    set<int> ret;
    vector<pair<int, double>> retp;
    set<int> topks(direct_topk.begin(), direct_topk.end());
    multimap<double, set<int>> heap;// init heap
    heap.emplace(0, topks);
    set<set<int>> cb;
    cb.insert(topks);

    while(ret.size()<X && !heap.empty()){
        pair<double, set<int>> popped=*heap.begin(); // deap copy, later erase
        heap.erase(heap.begin());
        for(int id:popped.second){
            if(ret.find(id)==ret.end()){
                ret.insert(id);
                retp.emplace_back(id, popped.first);
                cout<<ret.size()<<":"<<id<<","<<popped.first<<endl;
            }
        }
        unordered_set<int> cdd_replace;
        int m=0;
        for(int top: popped.second){
            for(int adj: ch_obj.get_neighbor_vertex(top)){
                cdd_replace.insert(adj);
            }
            m=max(ch_obj.get_option_layer(top), m);
        }
        for(int mp1_opt: ch_obj.get_layer(m+1)){
            cdd_replace.insert(mp1_opt);
        }
        for(int f:ret){
            cdd_replace.erase(f);
        }

        unordered_map<int, int> dominated_cnt;
        for(int i:cdd_replace){
            dominated_cnt[i]=ch_obj.dominated_map[i].size();
        }
        for(int opt:cdd_replace){
            auto iter=ch_obj.dominated_map.find(opt);
            if(iter!=ch_obj.dominated_map.end()){
                for(int i:popped.second){
                    if(iter->second.find(i)!=iter->second.end()){
                        --dominated_cnt[opt];
                    }
                }
            }
        }
        vector<int> ch_upd;
        for (int i:cdd_replace) {
            auto iter=dominated_cnt.find(i);
            if(iter!=dominated_cnt.end()&&iter->second<=0){
                ch_upd.push_back(i);
            }
        }
        for (int i:ch_upd) {
            for(int j=0;j<popped.second.size();++j){
                vector<int> cdd(popped.second.begin(), popped.second.end());
                cdd[j]=i;
                set<int> cdds(cdd.begin(), cdd.end());
                if(cb.find(cdds)!=cb.end()){
                    continue;
                }
                cb.insert(cdds);
                // TODO test feasible, input is cdd, ch_upd-i
                vector<int> le;
                for(int id:ch_obj.rskyband){
                    if(cdds.find(id)==cdds.end()){
                        le.push_back(id);
                    }
                }
                double radius=non_order_feasible(cdd, le, PointSet, w);
                bool feasible=radius!=INFINITY;
                if(feasible){
                    heap.emplace(radius, cdds);
                }
            }
        }

    }
    return retp;
}

int non_order_utk_efficient(float **PointSet, int pt_cnt, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                            vector<pair<int, double>> &utk_option_ret,
                            vector<pair<double, region*>> &utk_cones_ret){

    // two return values
    // 1. array of <option, topk_radius>
    // 2. array of <topk, region>
    auto begin = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;

    // 1. fetch 1-skyband with m options
    // 2. continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimated \rho^* from step 2
    // 4. apply \rho^* to get r-k-skyband (this is for pruning irrelevent options)
    // 5. apply \rho^* to get initial domain (this is for pruning unnecessary regions)
    // 6. get the top regions of CH1's options
    // 7. init ORU heap with top1's top region
    // 8. apply ORU, that is recursively divide regions


    // 1. begin: fetch 1-skyband with m options
    if(X<=k){
        vector<int> input(pt_cnt);
        iota(input.begin(), input.end(), 1);
        vector<int> ret=computeTopK(dim, PointSet, input, w, k);
        for(int i:ret){
            utk_option_ret.emplace_back(i, 0.0);
        }
        return utk_option_ret.size();
    }
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, PointSet);
    pair<int, float> next={-1, INFINITY};
    cout<< "begin fetch CH1"<<endl;
    vector<int> CH_1_X_opt;
    int top1=0;
    bool top1f= false;
    while(CH_1_X_opt.size() < X){
        next=get_next_obj.get_next();
        cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
        CH_1_X_opt.push_back(next.first);
        if(!top1f){
            top1=CH_1_X_opt.back();
            top1f=true;
        }
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " fetch top-m finish\n";
    ord_m_t+=elapsed_seconds.count();// for stat time
    // 1. end: fetch 1-skyband with m options

    // a 3-d example of square_vertex, square_vertex_cnt=4
    // point 0: (max(points[:, 0]), 0, 0)
    // point 1: (0, max(points[:, 1]), 0)
    // point 2: (0, 0, max(points[:, 2]))
    // point 3: \vec{0}
    const int square_vertex_cnt=dim+1;
    vector<vector<double>> square_vertexes(1, vector<double>(dim));

    // qhull class in lib qhull
    // init qhull with top X options
    CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " first time build qhull finish\n";

    // rho_star is computed when CH_1.size()=X
    int cnt=0;
    // 2. begin: continue step 1 to fill CH1 until CH1 contains m options

    int minc=X, maxc=2*X, midc=X;
    while(CH_1_X_opt.size() < X){ // find the maxc
        midc=maxc; // make sure correct rho_star
        if(get_next_obj.interval.size()<maxc){
            while(get_next_obj.interval.size()<maxc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <maxc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();
        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()>=X){
            break;
        }else{
            minc=maxc;
            maxc*=2;
        }
    }
    while(CH_1_X_opt.size() != X && minc < maxc){ // while(CH_1.size<X)
        midc=(maxc+minc)/2;
        if(get_next_obj.interval.size()<=midc){
            while(get_next_obj.interval.size()<midc){
                next=get_next_obj.get_next();
                cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
                if(next.second==INFINITY){
                    break;
                }
                update_square_vertexes(square_vertexes, PointSet[next.first], dim);
                CH_1_X_opt.push_back(next.first);
                cnt++;
            }
        }else{
            CH_1_X_opt.clear();
            for (int j = 0; j <midc ; ++j) {
                CH_1_X_opt.push_back(get_next_obj.interval[j].first);
            }
        }

        auto lbegin = chrono::steady_clock::now();
        CH_1_X_opt=build_qhull(CH_1_X_opt, PointSet, square_vertexes);
        auto lnow = chrono::steady_clock::now();
        chrono::duration<double> lelapsed_seconds= lnow-lbegin;
        cout << lelapsed_seconds.count()<<endl;
        cout<< cnt<<" rebuild qhull finish "<< CH_1_X_opt.size() <<endl;
        if(next.second==INFINITY){
            break;
        }
        if(CH_1_X_opt.size()==X){
            break;
        }else if(CH_1_X_opt.size()>X){
            maxc=midc-1;
        }else{
            minc=midc+1;
        }
    }

    // 2. end: continue step 1 to fill CH1 until CH1 contains m options
    // 3. getting estimate \rho^* from step 2
    float rho_star=get_next_obj.interval[min(midc-1, get_next_obj.interval.size()-1)].second;

    cout << "init rho_star: "<<rho_star<<endl;
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish rho_star compute\n";
    rho_star_t+=elapsed_seconds.count(); // for stat time
    // use known X version code to fetch rskyband options,
    // bear in mind such that we init \rho as \rho_star and X as INFINITY
    // 4. begin: apply \rho^* to get r-k-skyband
    vector<pair<long int, float>> interval;
    computeRho(dim, k, INFINITY, w, *rtree, PointSet, interval, rho_star);
    vector<int> rskyband_CS;
    for (pair<long int, float> &p:interval) {
        rskyband_CS.push_back(p.first);
    }
    // 4. end: apply \rho^* to get r-k-skyband
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " rskyband size: "<<rskyband_CS.size()<< "\n";
    rskyband_t+=elapsed_seconds.count();
    vector<vector<double>> tmp;
    // 5. begin: apply \rho^* to get initial domain
    cout<< elapsed_seconds.count() << " begin generate domain" << endl;
//    vector<vector<double>> begin_region=points_to_halfspace(tmp);
    vector<vector<double>> begin_region;

    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() <<"finish generate domain \n";
    // 5. end: apply \rho^* to get initial domain

    // 6. begin: get the top regions of CH1's options
    cout << " begin find top1 region" << endl;
    multimap<double, region*> id_radius; // <radius, region>
    vector<vector<double>> square_vertexes2(square_vertex_cnt, vector<double>(dim));
    // init convex hull object, this object can conveniently return CH_i given input i
    ch ch_obj(rskyband_CS, PointSet, dim);
    const vector<int>& top1_idxes=ch_obj.get_layer(1);// get CH1
    cout<<"finish finding CH1"<<endl;

//    unordered_map<int, vector<vector<double>>> top1_region=top_region(top1_idxes, PointSet, square_vertexes2);
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count() << " finish find top1 region" << endl;
    ch_1_t+=elapsed_seconds.count();

    cout<< elapsed_seconds.count() << "starting recursively get top regions\n";
    // 8. begin: apply ORU, that is recursively divide regions
    utk_option_ret=non_order_sensitive_ORU(dim, PointSet, w, X, k, ch_obj);

    // 8. end: apply ORU, that is recursively divide regions
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    cout<< " finish recursively get top regions\n";
    rec_oru_t+=elapsed_seconds.count();
}
