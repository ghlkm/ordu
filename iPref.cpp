#include "iPref.h"
#include <lp_lib.h>
//#include "qp_solver.h"
#include "math_lib.h"
#include "qhull_user.h"
#include "lp_user.h"
#include <chrono>
#include <RboxPoints.h>
#include "qp_solver2.h"
#include "vector_operator.h"

extern int objCnt;

extern qp_solver qp;

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

float computeDis(const vector<float> &tmpHS, const vector<float> &userpref, double *&solution)
{
    float ret=qp.update_w_h_solve(userpref, tmpHS, solution);
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

float computeRho(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, unordered_map<long int, RtreeNode*> &a_ramTree, float* PG[],
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
            node = a_ramTree[pageID];
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
    return topK_dominate_radius.begin()->first;
}

vector<double> unknown_X_node::radius_LB_point(){
    return topK_dominate_radius.begin()->second;
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
        vector<double> EMPTY;
        topK_dominate_radius.emplace(INFINITY, EMPTY);
    }else{
        vector<float> tmpHS(data, data+dim);
        for (int i = 0; i < dim; ++i) {
            tmpHS[i]-=other_option[i];
        }
//        float tmpDis = computeDis(tmpHS, w);
        double*solution= nullptr;
        float tmpDis = computeDis(tmpHS, w, solution);
        vector<double> sv;
        if(solution!= nullptr){
            sv=vector<double>(solution, solution+dim);
        }
        topK_dominate_radius.emplace(tmpDis, sv);
    }
}

void unknown_X_node::update_radius_erase(const float* other_option, vector<float> &w){
    update_radius(other_option, w);
    topK_dominate_radius.erase(topK_dominate_radius.begin());
}


unknown_x_efficient::unknown_x_efficient(const int dim, int K, vector<float> &userPref,
        Rtree &aRtree, unordered_map<long int, RtreeNode*> &aRamTree, float **pg) : a_ramTree(aRamTree), a_rtree(aRtree) {
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
                if(C.begin()->first==INFINITY){
                    return {-1, INFINITY};
                }
                interval.emplace_back(C.begin()->second->page_id - MAXPAGEID, C.begin()->first);
                assert(!C.begin()->second->radius_LB_point().empty());
                drill_position.emplace_back(C.begin()->second->radius_LB_point());
//                cout<<"drill point:"<<C.begin()->second->radius_LB_point()<<endl;
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
            }
            else {
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
            RtreeNode *node = a_ramTree[pageID];
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
    int dim=square_vertexes[0].size();
    auto begin = chrono::steady_clock::now();
    square_vertexes.clear();
    square_vertexes.emplace_back(dim);
    vector<int> nd_opt_idxes=non_dominate_set(opt_idxes, PG, dim); // must non-dominated each other
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

extern ofstream myfile;

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
        cout << "time: 1," << i <<", 0"  << "\n";
        if(GEN_F1) myfile << "time: 1," << i <<", 0"  << "\n";
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
                    cout << "time: " << options.size() << ", " << opt <<"," << elapsed_seconds.count()<<"\n";
                    if(GEN_F1) myfile << "time: " << options.size() << ", " << opt <<"," << elapsed_seconds.count()<<"\n";
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

            vector<int> tc=vector<int>(ch_upd.begin(), ch_upd.end());
            vector<int> new_ch_upd=r_dominate_filter(tc, popped.second->cone, PG);
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
            delete(popped.second);
            popped.second=nullptr;
        }

        if(popped.second && popped.second->topk.size()!=k && !new_option && !newOption){
            delete(popped.second);   // to save mem
        }
    }
    cout<<"cnt: "<<cnt<<"\n";
    for (auto left:id_radius) {
        delete (left.second);
    }
    return cnt;
}

class topi_region{
    double radius;
public:
    int opt_i; // for root, this value is -1
    int top_what;
    topi_region *parent;
    vector<int> child_constraint;
    unordered_set<int> candidate ;//
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
    vector<int> get_rtopi(){
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
    void set_candidate(const unordered_set<int> &c, const vector<vector<double>>& UNSED3, float **PG, int dim, int k){
        vector<double> this_w(dim);
        vector<float> UNUSED(dim);
        vector<int> EMPTY;
//        feasible(UNUSED, H1, 0, EMPTY, PG, this_w);
        multimap<double, int, greater<double>> find_top_ip1_to_k;
        for(int opt:c){
            double dot=0;
            for (int i = 0; i <dim ; ++i) {
                dot+=(PG[opt][i]+SIDELEN)*this_w[i];
            }
            find_top_ip1_to_k.emplace(dot, opt);
        }
        int size=k-this->top_what;
        vector<int> topi= this->get_rtopi();
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
                double ldist = feasible(UNUSED, UNSED3, opt, tmp, PG, UNUSED2);
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


    double main_diagonal(float **PG, int dim, const vector<int> &constr, vector<double> &center_mbb_ret){
//        double main_diagonal(vector<int> &order_vals_greater, vector<int> &unordered_vals_less, float **PG, int dim);
        vector<int> rtopi= this->get_rtopi();
        vector<int> topi(rtopi.rbegin(), rtopi.rend());
//        remove print info below later
//        cout<<"***\n";
//        for(auto i: rtopi){
//            for (int j = 0; j < dim; ++j) {
//                cout<<PG[i][j]<<",";
//            }
//            cout<<"\n";
//        }
        return ::main_diagonal(topi, constr, PG, dim, center_mbb_ret);
    }

    double feasible(const vector<float>& w, const vector<vector<double>>& H1, int opt, vector<int>& cmp, float **PG,
                    vector<double> &retv){
          // TODO: there are huge space to optimize here
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

template<typename VVT>
vector<int> topk_single(int k,  VVT &P, vector<double>&w, double uHeat,
                        Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree){
    RtreeNode* node;
    priority_queue<pair<double, int>> heap;
    double bound=uHeat;
    heap.emplace(INFINITY, rtree_rt->m_memory.m_rootPageID);
    int dim=w.size();
    vector<double> tmp_v(dim, 0.0);
    vector<int> topk;
    double tmp_score;
    long pageID;

    while(!heap.empty() && topk.size()<k){
        tmp_score=heap.top().first;
        pageID=heap.top().second;
        heap.pop();
        if (pageID >= MAXPAGEID){ // an option
            topk.push_back(pageID-MAXPAGEID);
        }else{ // an inter node
            node = ramTree[pageID];
            if (node->isLeaf()){ // this inter node contains leaves
                for (int i = 0; i < node->m_usedspace; i++){
                    tmp_score=P[node->m_entry[i]->m_id] * w;
                    if(tmp_score>=bound){
                        // the usual topk bound is set to 0
                        // we here can also set it as MaxMin_k
                        heap.emplace(tmp_score, node->m_entry[i]->m_id + MAXPAGEID);
                    }
                }
            }
            else{
                for (int i = 0; i < node->m_usedspace; i++){
                    auto &tmp=node->m_entry[i]->m_hc.getUpper();
                    for (int j = 0; j < dim; j++){
                        tmp_v[j]=tmp[j];
                    }
                    tmp_score=tmp_v*w;
                    if (tmp_score>=bound){
                        heap.emplace(tmp_score, node->m_entry[i]->m_id);
                    }
                }
            }
        }
    }
    return topk;
}


inline vector<double> multiple(const vector<vector<double>> &vv, const vector<int>& v){
    vector<double> ret(vv[0].size());
    for (int i = 0; i < v.size(); ++i) {
        ret+=v[i]*vv[i];
    }
    return ret;
}

template<typename VVT>
void dfs_combination_fix_far_point(vector<int>& unit, const vector<vector<double>> &base_norm_delta, int idx, int left,
                                   unordered_map<int, double> &topk_maxm, const vector<float>&w, VVT& P, int dim, int k, int m,
                                   Rtree *rtree_rt, unordered_map<long int, RtreeNode *> &ramTree, const double length){
    // base, dim-1 elements, length=delta/(dim-1)
    // unit, dim-1 elements
    if(idx==base_norm_delta.size()-1){
        // last
        unit[idx]=left;
        vector<double> drill_position=multiple(base_norm_delta, unit)+w;
        // find topk
//        cout<<"begin to drill topk"<<endl;
        vector<int> topk=topk_single(k, P, drill_position, 0, rtree_rt, ramTree);
        for(int opt: topk){
            auto iter=topk_maxm.find(opt);
            if(iter==topk_maxm.end()){
                topk_maxm[opt]=length;
                cout<<topk_maxm.size()<<":"<<opt<<","<<length<<endl;
            }
        }
        unit[idx]=-left;
        vector<double> drill_position2=multiple(base_norm_delta, unit)+w;
        vector<int> topk2=topk_single(k, P, drill_position2, 0, rtree_rt, ramTree);
        for(int opt: topk2){
            auto iter=topk_maxm.find(opt);
            if(iter==topk_maxm.end()){
                topk_maxm[opt]=length;
                cout<<topk_maxm.size()<<":"<<opt<<","<<length<<endl;
            }
        }
    }else{
        unit[idx]=0;
        dfs_combination_fix_far_point(unit, base_norm_delta, idx+1, left, topk_maxm, w, P, dim, k, m, rtree_rt, ramTree, length);
        for (int i = 1; i <= left; ++i) {
            unit[idx]=i;
            dfs_combination_fix_far_point(unit, base_norm_delta, idx+1, left-i, topk_maxm, w, P, dim, k, m, rtree_rt, ramTree, length);
            unit[idx]=-i;
            dfs_combination_fix_far_point(unit, base_norm_delta, idx+1, left-i, topk_maxm, w, P, dim, k, m, rtree_rt, ramTree, length);
        }
    }
}

class box{
public:
//    vector<int> coor;
    vector<double> cntr;
//    explicit box(const vector<int> &co): coor(co){}
    explicit box(const vector<double> &co): cntr(co){}

    box(const box& other): cntr(other.cntr){}

    box& operator=(const box& other){
        this->cntr=other.cntr;
    }

//    inline int cal_L2_int_dist(){
////        double ret=0;
//        int j=0;
//        for(int i: coor){
//            j+=i*i;
//        }
//        return j;
//    }

    inline double cal_dist(vector<double> &w) const{
        vector<double> tmp=cntr-w;
        return sqrt(tmp*tmp);
    }

    inline double cal_dist(vector<float> &w) const{
        vector<double> tmp=cntr-w;
        return sqrt(tmp*tmp);
    }

//    vector<double> center(const vector<double> &base, const vector<vector<double>> &baseV){
//        vector<double> ret(base);
//        for (int i = 0; i < coor.size(); ++i) {
//            ret+=coor[i]*baseV[i];
//        }
//        return ret;
//    }
    const vector<double>& center() const{
        return cntr;
    }
};

bool add_box(multimap<double, box*> &heap, vector<int> &tmp, vector<float> &base,
             int dim_m1, vector<vector<double>> &e){

    vector<double> center(base.begin(), base.end());
    for (int i1 = 0; i1 < e.size(); ++i1) {
        center += tmp[i1] * e[i1];
    }
    // check whether outside the boundary
    // such as reduce out of border cell visited
    bool flag = true;
    for (double i1: center) {
        if (!(0.0 <= i1 && i1 <= 1.0)) {
            flag = false;
            break;
        }
    }
    // if not outside the boundary
    // TODO this is a complicated case need to discuss, border condition
    if (flag) {
        box *b=new box(center);
        heap.emplace(b->cal_dist(base), b);
        return true;
    }else{
        return false;
    }
}
void append_boxes(vector<float> &base, int n, multimap<double, box*> &heap,
                  vector<vector<double>> &e, int idx, int mid, int dim_m1, vector<int>& tmp){
    if(idx<mid){
        for (int i = -n; i <= n; ++i) {
            tmp[idx]=i;
            append_boxes(base, n, heap, e, idx+1, mid, dim_m1,tmp);
        }
    }else if(idx>=dim_m1){
        add_box(heap, tmp, base, dim_m1, e);
    }else if(idx!=mid){
        for (int i = -(n-1); i <= n-1; ++i) {
            tmp[idx]=i;
            append_boxes(base, n, heap, e, idx+1, mid, dim_m1,tmp);
        }
    }else{
        append_boxes(base, n, heap, e, idx+1, mid, dim_m1,tmp);
    }
}

uint32_t append_boxes(vector<float> &base, int n, multimap<double, box*> &heap,
                      vector<vector<double>> &e){
    uint32_t ret_cnt=0;
    int dim_m1=base.size()-1;
    vector<int> tmp(dim_m1);
    for (int i = 0; i <dim_m1 ; ++i) {
        size_t bsize=heap.size();
        tmp[i]=n;
        append_boxes(base, n, heap, e, 0, i, dim_m1,tmp);
        tmp[i]=-n;
        append_boxes(base, n, heap, e, 0, i, dim_m1,tmp);
        ret_cnt+=heap.size()-bsize;
    }
    return ret_cnt;
}

vector<pair<int, double>> drill_by_grid(vector<float>&w,float** P, int dim, int k, int m, double delta,
                                        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree, unordered_map<int, double> &p1_options){
    // for doru
    // progressively drilling box with side distance
    // e.g., drill all cells (i*e1, j*e2) that i,j <= l
    // util the top of minimal heap is greater equal to l
    // then drill all cell that i, j <= l+1
    vector<double> wd(w.begin(), w.end());
    vector<vector<c_float>> e = gen_r_domain_basevec(dim); // size with d-1
    for(auto &ee: e){
        for(c_float &num: ee){
            num*=delta/(dim-1);
        }
    }
    cout<<"debug: delta/2/(dim-1), "<<delta/2/(dim-1)<<endl;
    for(auto &&ee:e){
        cout<<ee<<"\n";
    }
    int n=1;
    // because the box is small, we use object in heap
    multimap<double, box*> heap;
    // drill the top-k options respecting to input preference w
//    heap.emplace(0.0, new box(wd));
    int ret=0;
    vector<int> gtopk =topk_single_extend(k, P, wd, 0, rtree, ram_Tree);
    for(int opt: gtopk){
        p1_options[opt]=0;
    }
    unordered_map<int, double> ret_options;
    multimap<double, uint, greater<>> distOptHeap;
    for(auto &id_dis: p1_options){
        ret_options.insert(id_dis);
        distOptHeap.emplace(id_dis.second, id_dis.first);
    }
    while (ret_options.size()>m){
        uint opt=distOptHeap.begin()->second;
        ret_options.erase(opt);
        distOptHeap.erase(distOptHeap.begin());
    }
//    assert(distOptHeap.size()==m);
//    double rho=min(p1_options.begin(), p1_options.end(), [](auto &a, auto &b){return a->second<b->second;})->second;
    while (distOptHeap.begin()->first > n * delta / (dim - 1)) {
        uint32_t new_boxes=append_boxes(w, n, heap, e); // adding (2*n)^(dim-1) - (2*(n-1))^(dim-1) boxes
        ret+=new_boxes;
        if(new_boxes==0){
            break;
        }
        while(!heap.empty() && heap.begin()->first<n*delta/(dim-1)){
            box *top=heap.begin()->second;
            const vector<double>& center=top->center();
            vector<int> topk=topk_single_extend(k, P, center, 0, rtree, ram_Tree);
            for(int opt:topk){
                auto iter=ret_options.find(opt);
                if(iter==ret_options.end()) {// new option
                    if(distOptHeap.begin()->first > heap.begin()->first) {
                        ret_options.emplace(opt, heap.begin()->first);
                        distOptHeap.emplace(heap.begin()->first, opt);
                        ret_options.erase(distOptHeap.begin()->second);
                        distOptHeap.erase(distOptHeap.begin());
//                        assert(distOptHeap.size()==m);
                    }
                }else{
                    if(iter->second > heap.begin()->first){
                        for(auto iter2=distOptHeap.begin(); iter2 != distOptHeap.end(); iter2++){
                            if(iter2->second==opt){
                                assert(iter2->first > heap.begin()->first);
                                distOptHeap.erase(iter2);
                                break;
                            }
                        }
                        distOptHeap.emplace(heap.begin()->first, opt);
                        iter->second=heap.begin()->first;
//                        assert(distOptHeap.size()==m);

                    }
                }
            }
            heap.erase(heap.begin());
            delete (top);
        }
        n++;
    }

    while(!heap.empty()){
        delete(heap.begin()->second);
        heap.erase(heap.begin());
    }

    cout<<"drilled boxes: "<<ret<<endl;
//    assert(distOptHeap.size()==m);
    vector<pair<int, double>> ret_idDis;
    for(auto iter=distOptHeap.rbegin(); iter != distOptHeap.rend(); ++iter){
        ret_idDis.emplace_back(iter->second, iter->first);
    }
//    assert(ret_idDis.size()==m);
    return ret_idDis;
}

vector<pair<int, double>> drill_by_grid(vector<float>&w,float** P, int dim, int k, int m, double delta,
                                        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree, unordered_map<int, double> &p1_options,
                                        doGraphPlus &dg){
    // for doru
    // progressively drilling box with side distance
    // e.g., drill all cells (i*e1, j*e2) that i,j <= l
    // util the top of minimal heap is greater equal to l
    // then drill all cell that i, j <= l+1
    vector<double> wd(w.begin(), w.end());
    vector<vector<c_float>> e = gen_r_domain_basevec(dim); // size with d-1
    for(auto &ee: e){
        for(c_float &num: ee){
            num*=delta/(dim-1);
        }
    }
    cout<<"debug: delta/2/(dim-1), "<<delta/2/(dim-1)<<endl;
    for(auto &&ee:e){
        cout<<ee<<"\n";
    }
    int n=1;
    // because the box is small, we use object in heap
    multimap<double, box*> heap;
    // drill the top-k options respecting to input preference w
//    heap.emplace(0.0, new box(wd));
    int ret=0;
    vector<int> gtopk =topk(P, objCnt+1, dim,k, w, dg);
    for(int opt: gtopk){
        p1_options[opt]=0;
    }
    unordered_map<int, double> ret_options;
    multimap<double, uint, greater<>> distOptHeap;
    for(auto &id_dis: p1_options){
        ret_options.insert(id_dis);
        distOptHeap.emplace(id_dis.second, id_dis.first);
    }
    while (ret_options.size()>m){
        uint opt=distOptHeap.begin()->second;
        ret_options.erase(opt);
        distOptHeap.erase(distOptHeap.begin());
    }
//    assert(distOptHeap.size()==m);
//    double rho=min(p1_options.begin(), p1_options.end(), [](auto &a, auto &b){return a->second<b->second;})->second;
    while (distOptHeap.begin()->first > n * delta / (dim - 1)) {
        uint32_t new_boxes=append_boxes(w, n, heap, e); // adding (2*n)^(dim-1) - (2*(n-1))^(dim-1) boxes
        ret+=new_boxes;
        if(new_boxes==0){
            break;
        }
        while(!heap.empty() && heap.begin()->first<n*delta/(dim-1)){
            box *top=heap.begin()->second;
            const vector<double>& center=top->center();
            vector<int> topk=topk_more(P,objCnt+1, dim, k, center, dg);
            for(int opt:topk){
                auto iter=ret_options.find(opt);
                if(iter==ret_options.end()) {// new option
                    if(distOptHeap.begin()->first > heap.begin()->first) {
                        ret_options.emplace(opt, heap.begin()->first);
                        distOptHeap.emplace(heap.begin()->first, opt);
                        ret_options.erase(distOptHeap.begin()->second);
                        distOptHeap.erase(distOptHeap.begin());
//                        assert(distOptHeap.size()==m);
                    }
                }else{
                    if(iter->second > heap.begin()->first){
                        for(auto iter2=distOptHeap.begin(); iter2 != distOptHeap.end(); iter2++){
                            if(iter2->second==opt){
                                assert(iter2->first > heap.begin()->first);
                                distOptHeap.erase(iter2);
                                break;
                            }
                        }
                        distOptHeap.emplace(heap.begin()->first, opt);
                        iter->second=heap.begin()->first;
//                        assert(distOptHeap.size()==m);

                    }
                }
            }
            heap.erase(heap.begin());
            delete (top);
        }
        n++;
    }

    while(!heap.empty()){
        delete(heap.begin()->second);
        heap.erase(heap.begin());
    }

    cout<<"drilled boxes: "<<ret<<endl;
//    assert(distOptHeap.size()==m);
    vector<pair<int, double>> ret_idDis;
    for(auto iter=distOptHeap.rbegin(); iter != distOptHeap.rend(); ++iter){
        ret_idDis.emplace_back(iter->second, iter->first);
    }
//    assert(ret_idDis.size()==m);
    return ret_idDis;
}


class fourQEle{
public:
    bool isVertex;
    vector<int> topk;
    vector<int> S_w;
    vector<int> S_NR;
    vector<double> closest_point;
    vector<vector<double>> vertexes;
    double bound;

    void setBound(float **P, const vector<int> &topkRec, int dim){
        bound=INFINITY;
        for(int r: topkRec){
            double score=0;
            for (int i = 0; i < dim; ++i) {
                score+=(P[r][i]+P[r][i+dim])/2*closest_point[i];
            }
            bound=min(bound, score);
        }
    }

    fourQEle(){
        isVertex=false;
    }

    explicit fourQEle(bool isV){
        isVertex=isV;
    }

    fourQEle(const vector<int> &r, vector<double>& feasiblePoint, bool isV){
        isVertex=true;
        closest_point=feasiblePoint;
        topk=r;
    }

//    fourQEle(const vector<int> &r, bool isV){
//        isVertex=isV;
//        topk=r;
//    }

    fourQEle(const vector<int> &r, const vector<int> &sw, const vector<int> &snr){
        topk=r;
        S_w=sw;
        S_NR=snr;
        isVertex=false;
    }

    fourQEle(const fourQEle& other)=default;

    inline bool outOfBoundary(const vector<double> &input)const{
        for (auto &val: input) {
            if(val<0 || val >1){
                return true;
            }
        }
        return false;
    }

    inline bool onBoundary(const vector<double> &input)const{
        for (auto &val: input) {
            if(abs(val)<=1e-5){
                return true;
            }
        }
        return false;
    }

    vector<vector<double>>& calGeometry(float **P, uint dim, vector<double> &w){
        auto begin = chrono::steady_clock::now();
        vector<double> innerPoint=(closest_point+0.001*(closest_point-w))*0.9;
//        vector<double> innerPoint=feasiblePoint(S_w, S_NR, P, dim);
//        cout << "debug: "<<S_w.size()<<"," <<S_NR.size()<<endl;
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
        act_ql_time+=elapsed_seconds.count();
        if(innerPoint.empty() || outOfBoundary(innerPoint)) return vertexes;
        // given halfspaces, return vertexes on \Sigma u_i = 1
        uint32_t halfspace_num=S_w.size()*S_NR.size()+1;
        realT *pointCoordinates=new realT[(dim+1)*halfspace_num];
        size_t counter=0;
        for(int l = 0; l < dim; ++l){
            pointCoordinates[counter]=1;
            ++counter;
        }
        pointCoordinates[counter]=-1;
        ++counter;
        for(int ri: S_w){
            for(int rj: S_NR){
                for (int l = 0; l < dim; ++l) {
                    assert(rj<objCnt && rj>0 && ri<objCnt && ri>0);
                    pointCoordinates[counter]=P[rj][l]-P[ri][l];
                    ++counter;
                }
                pointCoordinates[counter]=0;
                ++counter;
            }
        }
        begin = chrono::steady_clock::now();
        vertexes=halfspace2vertices(pointCoordinates, halfspace_num, innerPoint);
        now = chrono::steady_clock::now();
        elapsed_seconds= now-begin;
        act_qh_time+=elapsed_seconds.count();
        for(auto iter=vertexes.begin(); iter!=vertexes.end(); ){
            double s=sum(iter->begin(), iter->end());
            if(abs(s-1)>=1e-6 || outOfBoundary(*iter)){
                iter=vertexes.erase(iter);
            }else{
                ++iter;
            }
        }
        delete[](pointCoordinates);
        return vertexes;
    }

    vector<vector<double>>& calGeometry_open(float **P, uint dim, vector<double> &w){
        auto begin = chrono::steady_clock::now();
        vector<double> innerPoint=(closest_point+0.001*(closest_point-w))*0.9;
//        vector<double> innerPoint=feasiblePoint(S_w, S_NR, P, dim);
//        cout << "debug: "<<S_w.size()<<"," <<S_NR.size()<<endl;
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
        act_ql_time+=elapsed_seconds.count();
        if(innerPoint.empty() || outOfBoundary(innerPoint)) return vertexes;
        // given halfspaces, return vertexes on \Sigma u_i = 1
        uint32_t halfspace_num=S_w.size()*S_NR.size()+1;
        realT *pointCoordinates=new realT[(dim+1)*halfspace_num];
        size_t counter=0;
        for(int l = 0; l < dim; ++l){
            pointCoordinates[counter]=1;
            ++counter;
        }
        pointCoordinates[counter]=-1;
        ++counter;
        for(int ri: S_w){
            for(int rj: S_NR){
                for (int l = 0; l < dim; ++l) {
                    assert(rj<objCnt && rj>0 && ri<objCnt && ri>0);
                    pointCoordinates[counter]=P[rj][l]-P[ri][l];
                    ++counter;
                }
                pointCoordinates[counter]=0;
                ++counter;
            }
        }
        begin = chrono::steady_clock::now();
        vertexes=halfspace2vertices(pointCoordinates, halfspace_num, innerPoint);
        now = chrono::steady_clock::now();
        elapsed_seconds= now-begin;
        act_qh_time+=elapsed_seconds.count();
        for(auto iter=vertexes.begin(); iter!=vertexes.end(); ){
            double s=sum(iter->begin(), iter->end());
            if(abs(s-1)>=1e-6 || outOfBoundary(*iter) || onBoundary(*iter)){
                iter=vertexes.erase(iter);
            }else{
                ++iter;
            }
        }
        delete[](pointCoordinates);
        return vertexes;
    }


    vector<pair<uint, uint>> calGeometry2(float **P, uint dim, vector<double> &w){
        auto begin = chrono::steady_clock::now();
//        cout<<"debug 2: #"<<S_w <<" || "<<S_NR.size()<<endl;
        vector<pair<uint, uint>> ret;
//        vector<double> innerPoint;
//        qp_solver q=qp_solver(P, S_w, S_NR, dim, innerPoint);
//        if(innerPoint.empty()) return ret;
//        vector<double> innerPoint=feasiblePoint(S_w, S_NR, P, dim);
//        cout<<"debug "<<closest_point<<endl;
        vector<double> innerPoint=(closest_point+0.001*(closest_point-w))*0.9;
//        vector<double> innerPoint=feasiblePointReduce(S_w, S_NR, P, dim);
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
//        act_ql_time+=elapsed_seconds.count();
        if(innerPoint.empty()) return ret;
        // given halfspaces, return vertexes on \Sigma u_i = 1
        uint32_t halfspace_num=S_w.size()*S_NR.size()+1;
        realT *pointCoordinates=new realT[(dim+1)*halfspace_num];
        size_t counter=0;
        for(int ri: S_w){
            for(int rj: S_NR){
                for (int l = 0; l < dim; ++l) {
                    assert(rj<objCnt && rj>0 && ri<objCnt && ri>0);
                    pointCoordinates[counter]=P[rj][l]-P[ri][l];
                    ++counter;
                }
                pointCoordinates[counter]=0;
                ++counter;
            }
        }
        for(int l = 0; l < dim; ++l){
            pointCoordinates[counter]=1;
            ++counter;
        }
        pointCoordinates[counter]=-1;
        ++counter;
        begin = chrono::steady_clock::now();
        vector<uint> halfspace_idx=halfspace2nonRedundant(pointCoordinates, halfspace_num, innerPoint);
        delete[](pointCoordinates);
        now = chrono::steady_clock::now();
        elapsed_seconds= now-begin;
        act_qh_time+=elapsed_seconds.count();
        ret.reserve(halfspace_idx.size());
        for(uint idx: halfspace_idx){
            if(idx<halfspace_num-1) {
//                cout << "debug: fuck " << idx << ", " << halfspace_num << endl;
//                cout << "debug: fuck " << S_w[idx / S_NR.size()] << ", " << S_NR[idx % S_NR.size()] << endl;
                ret.emplace_back(S_w[idx / S_NR.size()], S_NR[idx % S_NR.size()]);
            }
        }
        return ret;
    }

};

bool same_kth(const vector<double> &u, vector<int> &records, float** P){
    if(records.size()<=1) return false;
    vector<double> scores;
    scores.reserve(records.size());
    for (int r: records) {
        scores.push_back(P[r]*u);
    }
    sort(scores.begin(),scores.end());
    return scores[1]-scores[0]<1e-5;
}

double isFeasibleFOUR(const vector<int> &topk, float** P,
fourQEle *R, int dim, fourQEle *Sv, vector<float> &w){
    R->S_w=no_dominate_set(topk, P, dim);
    double ret=non_order_feasible(R->S_w, Sv->S_NR, P, w, R->closest_point);
    if(ret!=INFINITY)  R->S_NR=Sv->S_NR;
//    auto now = chrono::steady_clock::now();
//    chrono::duration<double> elapsed_seconds= now-begin;
//    cout<<"debug, copy data and qp time: "<<elapsed_seconds.count()<<endl;
//    vector<vector<double>> A(R->S_w.size() * R->S_NR.size(), vector<double>(dim));
//    for (int m = 0; m < R->S_w.size(); ++m) {
//        for (int i = 0; i < R->S_NR.size(); ++i) {
//            for (int j = 0; j < dim; ++j) {
//                A[m*R->S_NR.size()+i][j]=P[R->S_NR[i]][j]-P[R->S_w[m]][j];
//            }
//        }
//    }
//    qp_solver qp_ser(w, A);
//    c_float *solution= nullptr;
//    double ret=qp_ser.qp_solve(solution);
//    if(solution!= nullptr) R->closest_point=vector<double>(solution, solution+dim);
    return ret;
}





double isFeasibleFOUR2(const vector<int> &topk, const vector<vector<int>>& doGraph, float** P,
                      fourQEle *R, int dim, fourQEle *Sv, vector<float> &w, pair<uint, uint>& kthRec){
    R->S_w=no_dominate_set(topk, P, dim);
    vector<int> R1mR2=doGraph[kthRec.first];
//    set_difference(doGraph[kthRec.first].begin(),  doGraph[kthRec.first].end(),
//                   doGraph[kthRec.second].begin(), doGraph[kthRec.second].end(),
//                   inserter(R1mR2, R1mR2.begin()));
    R1mR2.push_back(kthRec.second);
    sort(R1mR2.begin(), R1mR2.end());
    vector<int> R2mR1;
    set_difference(doGraph[kthRec.second].begin(),  doGraph[kthRec.second].end(),
                   doGraph[kthRec.first].begin(), doGraph[kthRec.first].end(),
                   back_inserter(R2mR1));
    R2mR1.push_back(kthRec.first);
    sort(R2mR1.begin(), R2mR1.end());
    vector<int> tmp;
    set_difference(Sv->S_NR.begin(), Sv->S_NR.end(),
                   R1mR2.begin(),  R1mR2.end(),
                   back_inserter( tmp));
    set_union(R2mR1.begin(), R2mR1.end(),
                   tmp.begin(), tmp.end(),
                   back_inserter(R->S_NR));
//    cout<<"debug: #"<<R->S_w <<" || "<<R->S_NR.size()<<endl;
    double ret=non_order_feasible(R->S_w, R->S_NR, P, w, R->closest_point);
    return ret;
}

vector<int> shrinkRecTopk(int k, float** P, vector<double> &w,
        Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree,
        vector<int> &retSw, vector<int> &retSnr){
    // topk return should be sorted
    // return topk record set, could be more than k records
    // return retSw, the shrink topk
    // return rtSnr, the shrink other records
    RtreeNode* node;
    priority_queue<pair<double, int>> heap;
    heap.emplace(INFINITY, lrtree->m_memory.m_rootPageID);
    int dim=w.size();
//    vector<double> tmp_v(dim, 0.0);
    vector<int> topk;
    double tmp_score;
    long pageID;
    double last_score=INFINITY;
//    bound=bound-1e-6; // for accuracy
    while(!heap.empty() && (topk.size()<k || heap.top().first>=last_score-1e-7)){ // 1e-6 for accuracy
        tmp_score=heap.top().first;
        if(topk.size()<=k) last_score=tmp_score;
        pageID=heap.top().second;
        heap.pop();
        if (pageID >= MAXPAGEID){ // an option
            topk.push_back(pageID-MAXPAGEID);
        }else{ // an inter node
            node = lramTree[pageID];
            if (node->isLeaf()){ // this inter node contains leaves
                for (int i = 0; i < node->m_usedspace; i++){
                    auto tmp=node->m_entry[i]->m_hc.getCenter();
                    tmp_score=0;
                    for (int j = 0; j < dim; ++j) {
                        tmp_score+=tmp[j]*w[j];
                    }
                    heap.emplace(tmp_score, node->m_entry[i]->m_id + MAXPAGEID);
                }
            }
            else{
                for (int i = 0; i < node->m_usedspace; i++){
                    auto &tmp=node->m_entry[i]->m_hc.getUpper();
                    tmp_score=0;
                    for (int j = 0; j < dim; j++){
                        tmp_score+=tmp[j]*w[j];
                    }
                    heap.emplace(tmp_score, node->m_entry[i]->m_id);
                }
            }
        }
    }
    retSw=no_dominate_set(topk, P, dim);
//    if(topk.size()>2*k){// abnormal
//        cout<<"debug: something might be incorrect here"<<endl;
//        cout<<w<<endl;
//        for(int r: topk){
//            cout<<dot(P[r], w, dim)<<", "<<r<<": ";
//            for (int i = 0; i < dim; ++i) {
//                cout<<P[r][i]<<", ";
//            }
//            cout<<endl;
//        }
//    }
//    while(!heap.empty()){
//        tmp_score=heap.top().first;
//        last_score=tmp_score;
//        pageID=heap.top().second;
//        heap.pop();
//        if (pageID >= MAXPAGEID){ // an option
//            bool flag=true;
//            int cur_r=pageID-MAXPAGEID;
//            for(int r: retSnr){
//                if(v1_dominate_v2(P[r], P[cur_r], dim)){
//                    flag=false;
//                    break;
//                }
//            }
//            if(flag) retSnr.push_back(cur_r);
//        }else{ // an inter node
//            node = lramTree[pageID];
//            if (node->isLeaf()){ // this inter node contains leaves
//                for (int i = 0; i < node->m_usedspace; i++){
//                    int cur_r=node->m_entry[i]->m_id;
//                    auto tmp=node->m_entry[i]->m_hc.getCenter();
//                    tmp_score=0;
//                    for (int j = 0; j < dim; ++j) {
//                        tmp_score+=tmp[j]*w[j];
//                    }
//                    bool flag=true;
//                    for(int r: retSnr){
//                        if(v1_dominate_v2(P[r], P[cur_r], dim)){
//                            flag=false;
//                            break;
//                        }
//                    }
//                    if(flag) heap.emplace(tmp_score, node->m_entry[i]->m_id + MAXPAGEID);
//                }
//            }
//            else{
//                for (int i = 0; i < node->m_usedspace; i++){
//                    auto &tmp=node->m_entry[i]->m_hc.getUpper();
//                    tmp_score=0;
//                    for (int j = 0; j < dim; j++){
//                        tmp_score+=tmp[j]*w[j];
//                    }
//                    bool flag=true;
//                    for(int r: retSnr){
//                        auto tmp2=vector<double>(dim);
//                        for (int j = 0; j < dim; ++j) {
//                            tmp2[j]=(P[r][j]+P[r][j+dim])/2.0;
//                        }
//                        if(v1_dominate_v2(tmp2, tmp, dim)){
//                            flag=false;
//                            break;
//                        }
//                    }
//                    if(flag) heap.emplace(tmp_score, node->m_entry[i]->m_id);
//                }
//            }
//        }
//    }
    return topk;
}

bool approximate_dominate(const float* r1, const float* r2, int dim){
    for (int i = 0; i < dim; ++i) {
        if(r2[i]-1e-7 >r1[i]){
            return false;
        }
    }
    return true;
}

vector<int> approximate_ksksyband(vector<int> &candidate, vector<int> &allR, float** P, int dim, int k){
    vector<int> ret;
    for(int r: candidate){
        int docnt=0;
        for (int rc: allR) {
            if(r==rc) continue;
            if(approximate_dominate(P[rc], P[r], dim)){
                docnt++;
                if(docnt>=k){
                    break;
                }
            }
        }
        if(docnt<k){
            ret.push_back(r);
        }
    }
    return ret;
}

vector<pair<int, double>> foru(vector<float>&w,float** P, int dim, int k, int m,
                               vector<vector<int>> &doGraph, vector<int> &kskybandR, vector<int> &OneSkybandR,
        Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree) {
    // A: k-skyband records and their dominated count
    act_qh_time=0;
    vector<double> wd(w.begin(), w.end());
    fourQEle Rw=fourQEle();
//    if(kskybandR.size()<200) {
//        Rw.topk = computeTopK_Extend(dim, P, kskybandR, w, k, Rw.S_w, Rw.S_NR);
//        sort(Rw.S_NR.begin(), Rw.S_NR.end());
//    }else {
        Rw.topk = shrinkRecTopk(k, P, wd, lrtree, lramTree, Rw.S_w, Rw.S_NR);
        vector<int> tmp;
        vector<int> tmp2;
        for(int r: Rw.S_w){
            set_union(doGraph[r].begin(), doGraph[r].end(),
                      tmp.begin(), tmp.end(), back_inserter(tmp2));
            tmp=tmp2;
            tmp2.clear();
        }
        set_union(OneSkybandR.begin(), OneSkybandR.end(),
                  tmp.begin(), tmp.end(), back_inserter(tmp2));
        tmp=Rw.topk;
        sort(tmp.begin(), tmp.end());
        set_difference(tmp2.begin(),  tmp2.end(),
                       tmp.begin(), tmp.end(),
                       back_inserter(Rw.S_NR));
//    }
    assert(Rw.topk.size()>=k);
    myHeap G=myHeap();
    for(int record: Rw.topk){
        G.insert(record, 0);
    }
    while(G.size()>m){
        G.pop();
    }
    while(G.size()<m){
        G.insert(-(int)G.size(), INFINITY);
    }
    assert(G.size()==m);
    sort(Rw.topk.begin(), Rw.topk.end()); // since Ts is a combination set, hence it needs a sort
    unordered_set<vector<int>, VectorHash<int>> Ts;
    Ts.insert(Rw.topk);

    unordered_set<myPoint<double>, myPointDoubleHash> Tv;
    unordered_set<myPoint<double>, myPointDoubleHash> To;

    // IMPORTANT: check every element that is popped from Q and "delete" them using their pointer
    multimap<double, fourQEle*> Q;
    Rw.closest_point=wd;
    auto begin = chrono::steady_clock::now();
    Rw.calGeometry(P, dim, wd);
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    cout<<"debug, Rw geo time: "<<elapsed_seconds.count()<<endl;
    for(auto &v: Rw.vertexes){
        myPoint<double> v_(v);
        if(Tv.find(v_)!=Tv.end()) continue;
        Tv.emplace(v_);
        fourQEle *Sv=new fourQEle(true);
        Sv->closest_point=v;
        double dist_v_w=dist(v, w);
        Q.insert({dist_v_w, Sv});
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<<"debug, Rw top-k time: "<<elapsed_seconds.count()<<endl;
    double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;
    int c1=0, c2=0, c3=0, c4=0;
    for (int j = 0; j < dim; ++j) {
        fourQEle *Sv=new fourQEle(true);
        Sv->closest_point=vector<double>(dim);
        Sv->closest_point[j]=1;
        double dist_v_w=dist(Sv->closest_point, w);
        Q.insert({dist_v_w, Sv});
    }
    while(!Q.empty() && G.topDist()>Q.begin()->first){
//        cout<<"top dist:"<<Q.begin()->first<<","<<Q.begin()->second->isVertex<<"#"<<Q.begin()->second->closest_point<<endl;
        if(Q.begin()->second->isVertex){
            begin = chrono::steady_clock::now();
            fourQEle* Sv=Q.begin()->second;
            double dist_v_w=Q.begin()->first;
            Q.erase(Q.begin());
            c1++;
            Sv->topk = shrinkRecTopk(k, P, Sv->closest_point, lrtree, lramTree, Sv->S_w, Sv->S_NR);
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t6+=elapsed_seconds.count();
//            int debug_s=0;
            for(int r: Sv->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_v_w<G.topDist()){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
                    G.pop();
//                    debug_s++;
                }else if(iter!=G.recordEnd() && iter->second > dist_v_w){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
//                    debug_s++;
                }
            }
//            if(debug_s>k){
//                cout<< "debug --------- ";
//                cout<< Sv->closest_point <<endl;
//                for(int r: Sv->topk){
//                    for (int i = 0; i < dim; ++i) {
//                        cout<< P[r][i]<<", ";
//                    }
//                    cout<<"\n";
//                }
//            }
            if(Sv->topk.size()<=k) {
                delete(Sv);
                continue;
            }
//            cout<<"debug: topk size: "<<Sv->topk.size()<<endl;
//            cout<<Sv->closest_point<<endl;
//            for (int j: Sv->topk) {
//                cout<<dot(P[j], Sv->closest_point, dim)<<"|";
//                for (int i = 0; i < dim; ++i) {
//                    cout<<P[j][i]<<", ";
//                }
//                cout<<endl;
//            }
            if(Sv->topk.size()>=2*k){// abnormal
                bool f=false;
                for (int i = 0; i < dim; ++i) {
                    if(Sv->closest_point[i]>1-1e-5){
                        f=true;
                    }
                }
                if(f){
                    delete(Sv);
                    continue;
                }else{
                    cout<<"too many top-k records, maybe bug in code"<<endl;
                    cout<<Sv->closest_point<<endl;
                    cout<<Sv->topk.size()<<endl;
                    for(int r: Sv->topk){
                        cout<<P[r]*Sv->closest_point<<": ";
                        for (int i = 0; i < dim; ++i) {
                            cout<<P[r][i]<<", ";
                        }
                        cout<<endl;
                    }
//                    cout<<"[";
//                    for(auto &t: Tv){
//                        cout<<t.data<<endl;
//                    }
//                    cout<<"]\n";
//                    exit(1);
                }
            }
            // -------------------------------------------
            double lastScore=P[Sv->topk.back()] * Sv->closest_point ;
            vector<int> fixedTopk;
            vector<int> unfixedTopk;
            for(auto iter=Sv->topk.rbegin(); iter!=Sv->topk.rend(); iter++){
                double tmpScore=P[*iter] * Sv->closest_point;
                if(abs(tmpScore-lastScore)<1e-6){
                    unfixedTopk.push_back(*iter);
                }else{
                    fixedTopk=vector<int>(iter, Sv->topk.rend());
                    break;
                }
            }
            unfixedTopk=approximate_ksksyband(unfixedTopk, Sv->topk, P, dim, k);
            // -------------------------------------------

            sort(Sv->topk.begin(), Sv->topk.end());
            if(Ts.find(Sv->topk)!=Ts.end() || fixedTopk.size()>=k) {
                delete(Sv);
                continue;
            }
            Ts.insert(Sv->topk);

            // select r from n do combination and check total dominate cnt (TDC)
            // return the the combination with maximal TDC
            // This will generate correct result for OSS-skyline, but a brute force solution
            // This practical fast since 1-skyline is small when d is small
//            begin=chrono::steady_clock::now();

            vector<int> tmp;
            vector<int> tmp2;
            for(int r: Sv->S_w){
                set_union(doGraph[r].begin(), doGraph[r].end(),
                          tmp.begin(), tmp.end(), back_inserter(tmp2));
                tmp=tmp2;
                tmp2.clear();
            }
            set_union(OneSkybandR.begin(), OneSkybandR.end(),
                      tmp.begin(), tmp.end(), back_inserter(tmp2));
            tmp=Sv->topk;
            sort(tmp.begin(), tmp.end());
            set_difference(tmp2.begin(),  tmp2.end(),
                           tmp.begin(), tmp.end(),
                           back_inserter(Sv->S_NR));

            int n=unfixedTopk.size();
            int r=k-fixedTopk.size();
            std::vector<bool> v(n);
//            cout<<"debug, (n, r): " <<n<<", "<<r<<endl;
            std::fill(v.begin(), v.begin() + r, true);
            do {
                vector<int> S;
//                S.reserve(k);
                S=fixedTopk;
                vector<int> S_com;
                S_com.reserve(n-r);
                for (int i = 0; i < n; ++i) {
                    if (v[i]) {
                        S.push_back(unfixedTopk[i]);
                    }else{
                        S_com.push_back(unfixedTopk[i]);
                    }
                }
                sort(S.begin(), S.end());
                if(Ts.find(S)!=Ts.end()) continue;
                c2++;
                Ts.insert(S);
                auto lbegin=chrono::steady_clock::now();
                auto *Rs=new fourQEle();
                for(int lessRec: S_com) Sv->S_NR.push_back(lessRec);
                double distRS=isFeasibleFOUR(S, P, Rs, dim, Sv, w);
                for(int lessRec: S_com) Sv->S_NR.pop_back();
                auto lnow = chrono::steady_clock::now();
                chrono::duration<double> lelapsed_seconds= lnow-lbegin;
                t7+=lelapsed_seconds.count();
                // Rs, we know Rs contains a point v
//                myPoint<double> tmp(Rs->closest_point);
//                if(distRS!=INFINITY && To.find(tmp)==To.end() && Tv.find(tmp)==Tv.end()){ // TODO 这个优化到底是不是正确的
                if(distRS!=INFINITY){
//                    To.insert(tmp);
                    Rs->topk=S;
                    // nearest point checking top-k
                    // delta point checking top-k
                    vector<double> innerPoint=(Rs->closest_point+0.001*(Rs->closest_point-w))*0.9;
                    if (!same_kth(Rs->closest_point, Rs->topk, P) && !same_kth(innerPoint, Rs->topk, P)){
                        Q.insert({distRS, Rs});
                    }else{
                        delete(Rs);
                    }
                }else{
                    delete(Rs);
                }
            } while (std::prev_permutation(v.begin(), v.end())); // forward next permutation
            delete(Sv);
        }
        else{
            begin=chrono::steady_clock::now();
            fourQEle* Rs=Q.begin()->second;
            double dist_Rs_w=Q.begin()->first;
            Q.erase(Q.begin());
//            int debug_s=0;
            for(int r: Rs->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_Rs_w<G.topDist()){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
                    G.pop();
//                    debug_s++;
                }else if(iter!=G.recordEnd() && iter->second > dist_Rs_w){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
//                    debug_s++;
                }
            }
//            if(debug_s>k){
//                cout<< "debug --------- ";
//                cout<< Rs->closest_point <<endl;
//                for(int r: Rs->topk){
//                    for (int i = 0; i < dim; ++i) {
//                        cout<< P[r][i]<<", ";
//                    }
//                    cout<<"\n";
//                }
//
//            }
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t3+=elapsed_seconds.count();
            begin = chrono::steady_clock::now();
            c3++;
            Rs->calGeometry(P, dim, wd);
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t4+=elapsed_seconds.count();
            begin = chrono::steady_clock::now();
            if(Rs->vertexes.size()<dim){
                c4++;
                delete(Rs);
                continue;
            }
            for(auto &v: Rs->vertexes){
                myPoint<double> v_(v);
                if(Tv.find(v_)!=Tv.end()) continue;
                Tv.emplace(v_);
                fourQEle *Sv=new fourQEle(true);
//                Sv->topk=shrinkRecTopk(k, P, v, lrtree, lramTree, A, Sv->S_w, Sv->S_NR);
//                sort(Sv->topk.begin(), Sv->topk.end());
//                if(Ts.find(Sv->topk)!=Ts.end()) {
//                    delete(Sv);
//                    continue;
//                }
//                Ts.insert(Sv->topk);
                Sv->closest_point=v;
                double dist_v_w=dist(v, w);
                Q.insert({dist_v_w, Sv});
            }
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t5+=elapsed_seconds.count();
            delete(Rs);
        }
    }
    cout<<"t1: "<<t1<<endl;
    cout<<"t2: "<<t2<<endl;
    cout<<"t3: "<<t3<<endl;
    cout<<"t4: "<<t4<<endl;
    cout<<"t5: "<<t5<<endl;
    cout<<"t6: "<<t6<<endl;
    cout<<"t7: "<<t7<<endl;
    cout<<"tq: "<<act_qp_time<<endl;
    cout<<"th: "<<act_qh_time<<endl;
    cout<<"tl: "<<act_ql_time<<endl;
    cout<<"c1: "<<c1<<endl;
    cout<<"c2: "<<c2<<endl;
    cout<<"c3: "<<c3<<endl;
    cout<<"c4: "<<c4<<endl;
    for(auto &e: Q){
        delete(e.second);
    }
    return G.report();
}

vector<pair<int, double>> foru_open(vector<float>&w,float** P, int dim, int k, int m,
                               vector<vector<int>> &doGraph, vector<int> &kskybandR, vector<int> &OneSkybandR,
                               Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree) {
    // A: k-skyband records and their dominated count
    act_qh_time=0;
    vector<double> wd(w.begin(), w.end());
    fourQEle Rw=fourQEle();
//    if(kskybandR.size()<200) {
//        Rw.topk = computeTopK_Extend(dim, P, kskybandR, w, k, Rw.S_w, Rw.S_NR);
//        sort(Rw.S_NR.begin(), Rw.S_NR.end());
//    }else {
    Rw.topk = shrinkRecTopk(k, P, wd, lrtree, lramTree, Rw.S_w, Rw.S_NR);
    vector<int> tmp;
    vector<int> tmp2;
    for(int r: Rw.S_w){
        set_union(doGraph[r].begin(), doGraph[r].end(),
                  tmp.begin(), tmp.end(), back_inserter(tmp2));
        tmp=tmp2;
        tmp2.clear();
    }
    set_union(OneSkybandR.begin(), OneSkybandR.end(),
              tmp.begin(), tmp.end(), back_inserter(tmp2));
    tmp=Rw.topk;
    sort(tmp.begin(), tmp.end());
    set_difference(tmp2.begin(),  tmp2.end(),
                   tmp.begin(), tmp.end(),
                   back_inserter(Rw.S_NR));
//    }
    assert(Rw.topk.size()>=k);
    myHeap G=myHeap();
    for(int record: Rw.topk){
        G.insert(record, 0);
    }
    while(G.size()>m){
        G.pop();
    }
    while(G.size()<m){
        G.insert(-(int)G.size(), INFINITY);
    }
    assert(G.size()==m);
    sort(Rw.topk.begin(), Rw.topk.end()); // since Ts is a combination set, hence it needs a sort
    unordered_set<vector<int>, VectorHash<int>> Ts;
    Ts.insert(Rw.topk);

    unordered_set<myPoint<double>, myPointDoubleHash> Tv;
    unordered_set<myPoint<double>, myPointDoubleHash> To;

    // IMPORTANT: check every element that is popped from Q and "delete" them using their pointer
    multimap<double, fourQEle*> Q;
    Rw.closest_point=wd;
    auto begin = chrono::steady_clock::now();
    Rw.calGeometry_open(P, dim, wd);
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    cout<<"debug, Rw geo time: "<<elapsed_seconds.count()<<endl;
    for(auto &v: Rw.vertexes){
        myPoint<double> v_(v);
        if(Tv.find(v_)!=Tv.end()) continue;
        Tv.emplace(v_);
        fourQEle *Sv=new fourQEle(true);
        Sv->closest_point=v;
        double dist_v_w=dist(v, w);
        Q.insert({dist_v_w, Sv});
    }
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<<"debug, Rw top-k time: "<<elapsed_seconds.count()<<endl;
    double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;
    int c1=0, c2=0, c3=0, c4=0;
    while(!Q.empty() && G.topDist()>Q.begin()->first){
//        cout<<"top dist:"<<Q.begin()->first<<","<<Q.begin()->second->isVertex<<"#"<<Q.begin()->second->closest_point<<endl;
        if(Q.begin()->second->isVertex){
            begin = chrono::steady_clock::now();
            fourQEle* Sv=Q.begin()->second;
            double dist_v_w=Q.begin()->first;
            Q.erase(Q.begin());
            c1++;
            Sv->topk = shrinkRecTopk(k, P, Sv->closest_point, lrtree, lramTree, Sv->S_w, Sv->S_NR);
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t6+=elapsed_seconds.count();
//            int debug_s=0;
            for(int r: Sv->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_v_w<G.topDist()){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
                    G.pop();
//                    debug_s++;
                }else if(iter!=G.recordEnd() && iter->second > dist_v_w){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
//                    debug_s++;
                }
            }
//            if(debug_s>k){
//                cout<< "debug --------- ";
//                cout<< Sv->closest_point <<endl;
//                for(int r: Sv->topk){
//                    for (int i = 0; i < dim; ++i) {
//                        cout<< P[r][i]<<", ";
//                    }
//                    cout<<"\n";
//                }
//            }
            if(Sv->topk.size()<=k) {
                delete(Sv);
                continue;
            }
//            cout<<"debug: topk size: "<<Sv->topk.size()<<endl;
//            cout<<Sv->closest_point<<endl;
//            for (int j: Sv->topk) {
//                cout<<dot(P[j], Sv->closest_point, dim)<<"|";
//                for (int i = 0; i < dim; ++i) {
//                    cout<<P[j][i]<<", ";
//                }
//                cout<<endl;
//            }
            if(Sv->topk.size()>=2*k){// abnormal
                bool f=false;
                for (int i = 0; i < dim; ++i) {
                    if(Sv->closest_point[i]>1-1e-5){
                        f=true;
                    }
                }
                if(f){
                    delete(Sv);
                    continue;
                }else{
                    cout<<"too many top-k records, maybe bug in code"<<endl;
                    cout<<Sv->closest_point<<endl;
                    cout<<Sv->topk.size()<<endl;
                    for(int r: Sv->topk){
                        cout<<P[r]*Sv->closest_point<<": ";
                        for (int i = 0; i < dim; ++i) {
                            cout<<P[r][i]<<", ";
                        }
                        cout<<endl;
                    }
//                    cout<<"[";
//                    for(auto &t: Tv){
//                        cout<<t.data<<endl;
//                    }
//                    cout<<"]\n";
//                    exit(1);
                }
            }
            // -------------------------------------------
            double lastScore=P[Sv->topk.back()] * Sv->closest_point ;
            vector<int> fixedTopk;
            vector<int> unfixedTopk;
            for(auto iter=Sv->topk.rbegin(); iter!=Sv->topk.rend(); iter++){
                double tmpScore=P[*iter] * Sv->closest_point;
                if(abs(tmpScore-lastScore)<1e-6){
                    unfixedTopk.push_back(*iter);
                }else{
                    fixedTopk=vector<int>(iter, Sv->topk.rend());
                    break;
                }
            }
            unfixedTopk=approximate_ksksyband(unfixedTopk, Sv->topk, P, dim, k);
            // -------------------------------------------

            sort(Sv->topk.begin(), Sv->topk.end());
            if(Ts.find(Sv->topk)!=Ts.end() || fixedTopk.size()>=k) {
                delete(Sv);
                continue;
            }
            Ts.insert(Sv->topk);

            // select r from n do combination and check total dominate cnt (TDC)
            // return the the combination with maximal TDC
            // This will generate correct result for OSS-skyline, but a brute force solution
            // This practical fast since 1-skyline is small when d is small
//            begin=chrono::steady_clock::now();

            vector<int> tmp;
            vector<int> tmp2;
            for(int r: Sv->S_w){
                set_union(doGraph[r].begin(), doGraph[r].end(),
                          tmp.begin(), tmp.end(), back_inserter(tmp2));
                tmp=tmp2;
                tmp2.clear();
            }
            set_union(OneSkybandR.begin(), OneSkybandR.end(),
                      tmp.begin(), tmp.end(), back_inserter(tmp2));
            tmp=Sv->topk;
            sort(tmp.begin(), tmp.end());
            set_difference(tmp2.begin(),  tmp2.end(),
                           tmp.begin(), tmp.end(),
                           back_inserter(Sv->S_NR));

            int n=unfixedTopk.size();
            int r=k-fixedTopk.size();
            std::vector<bool> v(n);
//            cout<<"debug, (n, r): " <<n<<", "<<r<<endl;
            std::fill(v.begin(), v.begin() + r, true);
            do {
                vector<int> S;
//                S.reserve(k);
                S=fixedTopk;
                vector<int> S_com;
                S_com.reserve(n-r);
                for (int i = 0; i < n; ++i) {
                    if (v[i]) {
                        S.push_back(unfixedTopk[i]);
                    }else{
                        S_com.push_back(unfixedTopk[i]);
                    }
                }
                sort(S.begin(), S.end());
                if(Ts.find(S)!=Ts.end()) continue;
                c2++;
                Ts.insert(S);
                auto lbegin=chrono::steady_clock::now();
                auto *Rs=new fourQEle();
                for(int lessRec: S_com) Sv->S_NR.push_back(lessRec);
                double distRS=isFeasibleFOUR(S, P, Rs, dim, Sv, w);
                for(int lessRec: S_com) Sv->S_NR.pop_back();
                auto lnow = chrono::steady_clock::now();
                chrono::duration<double> lelapsed_seconds= lnow-lbegin;
                t7+=lelapsed_seconds.count();
                // Rs, we know Rs contains a point v
//                myPoint<double> tmp(Rs->closest_point);
//                if(distRS!=INFINITY && To.find(tmp)==To.end() && Tv.find(tmp)==Tv.end()){ // TODO 这个优化到底是不是正确的
                if(distRS!=INFINITY){
//                    To.insert(tmp);
                    Rs->topk=S;
                    // nearest point checking top-k
                    // delta point checking top-k
                    vector<double> innerPoint=(Rs->closest_point+0.001*(Rs->closest_point-w))*0.9;
                    if (!same_kth(Rs->closest_point, Rs->topk, P) && !same_kth(innerPoint, Rs->topk, P)){
                        Q.insert({distRS, Rs});
                    }else{
                        delete(Rs);
                    }
                }else{
                    delete(Rs);
                }
            } while (std::prev_permutation(v.begin(), v.end())); // forward next permutation
            delete(Sv);
        }
        else{
            begin=chrono::steady_clock::now();
            fourQEle* Rs=Q.begin()->second;
            double dist_Rs_w=Q.begin()->first;
            Q.erase(Q.begin());
//            int debug_s=0;
            for(int r: Rs->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_Rs_w<G.topDist()){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
                    G.pop();
//                    debug_s++;
                }else if(iter!=G.recordEnd() && iter->second > dist_Rs_w){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
//                    debug_s++;
                }
            }
//            if(debug_s>k){
//                cout<< "debug --------- ";
//                cout<< Rs->closest_point <<endl;
//                for(int r: Rs->topk){
//                    for (int i = 0; i < dim; ++i) {
//                        cout<< P[r][i]<<", ";
//                    }
//                    cout<<"\n";
//                }
//
//            }
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t3+=elapsed_seconds.count();
            begin = chrono::steady_clock::now();
            c3++;
            Rs->calGeometry_open(P, dim, wd);
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t4+=elapsed_seconds.count();
            begin = chrono::steady_clock::now();
            if(Rs->vertexes.size()<dim){
                c4++;
                delete(Rs);
                continue;
            }
            for(auto &v: Rs->vertexes){
                myPoint<double> v_(v);
                if(Tv.find(v_)!=Tv.end()) continue;
                Tv.emplace(v_);
                fourQEle *Sv=new fourQEle(true);
//                Sv->topk=shrinkRecTopk(k, P, v, lrtree, lramTree, A, Sv->S_w, Sv->S_NR);
//                sort(Sv->topk.begin(), Sv->topk.end());
//                if(Ts.find(Sv->topk)!=Ts.end()) {
//                    delete(Sv);
//                    continue;
//                }
//                Ts.insert(Sv->topk);
                Sv->closest_point=v;
                double dist_v_w=dist(v, w);
                Q.insert({dist_v_w, Sv});
            }
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t5+=elapsed_seconds.count();
            delete(Rs);
        }
    }
    cout<<"t1: "<<t1<<endl;
    cout<<"t2: "<<t2<<endl;
    cout<<"t3: "<<t3<<endl;
    cout<<"t4: "<<t4<<endl;
    cout<<"t5: "<<t5<<endl;
    cout<<"t6: "<<t6<<endl;
    cout<<"t7: "<<t7<<endl;
    cout<<"tq: "<<act_qp_time<<endl;
    cout<<"th: "<<act_qh_time<<endl;
    cout<<"tl: "<<act_ql_time<<endl;
    cout<<"c1: "<<c1<<endl;
    cout<<"c2: "<<c2<<endl;
    cout<<"c3: "<<c3<<endl;
    cout<<"c4: "<<c4<<endl;
    for(auto &e: Q){
        delete(e.second);
    }
    return G.report();
}

vector<pair<int, double>> foru_dg(vector<float>&w,float** P, int dim, int k, int m, doGraphPlus &dg) {
    // A: k-skyband records and their dominated count
    vector<int> OneSkybandR=dg.oneSkyband;
    vector<vector<size_t>> doGraph=dg.outGraph;
    act_qh_time=0;
    auto begin = chrono::steady_clock::now();
    vector<double> wd(w.begin(), w.end());
    fourQEle Rw=fourQEle();
    Rw.topk=topk(P, objCnt+1, dim, k, w, dg);
    Rw.S_w=no_dominate_set(Rw.topk, P, dim);

    vector<int> tmp;
    vector<int> tmp2;
    for(int r: Rw.S_w){
        set_union(doGraph[r].begin(), doGraph[r].end(),
                  tmp.begin(), tmp.end(), back_inserter(tmp2));
        tmp=tmp2;
        tmp2.clear();
    }
    set_union(OneSkybandR.begin(), OneSkybandR.end(),
              tmp.begin(), tmp.end(), back_inserter(tmp2));
    tmp=Rw.topk;
    sort(tmp.begin(), tmp.end());
    set_difference(tmp2.begin(),  tmp2.end(),
                   tmp.begin(), tmp.end(),
                   back_inserter(Rw.S_NR));

    assert(Rw.topk.size()>=k);
    myHeap G=myHeap();
    for(int record: Rw.topk){
        G.insert(record, 0);
    }
    while(G.size()>m){
        G.pop();
    }
    while(G.size()<m){
        G.insert(-(int)G.size(), INFINITY);
    }
    assert(G.size()==m);
    sort(Rw.topk.begin(), Rw.topk.end()); // since Ts is a combination set, hence it needs a sort
    unordered_set<vector<int>, VectorHash<int>> Ts;
    Ts.insert(Rw.topk);

    unordered_set<myPoint<double>, myPointDoubleHash> Tv;
    unordered_set<myPoint<double>, myPointDoubleHash> To;

    // IMPORTANT: check every element that is popped from Q and "delete" them using their pointer
    multimap<double, fourQEle*> Q;
    Rw.closest_point=wd;
    Rw.calGeometry(P, dim, wd);
    for(auto &v: Rw.vertexes){
        myPoint<double> v_(v);
        if(Tv.find(v_)!=Tv.end()) continue;
        Tv.emplace(v_);
        fourQEle *Sv=new fourQEle(true);
        Sv->closest_point=v;
        double dist_v_w=dist(v, w);
        Q.insert({dist_v_w, Sv});
    }
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    cout<<"debug, Rw top-k time: "<<elapsed_seconds.count()<<endl;
    double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;
    int c1=0, c2=0, c3=0, c4=0;
    for (int j = 0; j < dim; ++j) {
        fourQEle *Sv=new fourQEle(true);
        Sv->closest_point=vector<double>(dim);
        Sv->closest_point[j]=1;
        double dist_v_w=dist(Sv->closest_point, w);
        Q.insert({dist_v_w, Sv});
    }
    while(!Q.empty() && G.topDist()>Q.begin()->first){
//        cout<<"top dist:"<<Q.begin()->first<<","<<Q.begin()->second->isVertex<<"#"<<Q.begin()->second->closest_point<<endl;
        if(Q.begin()->second->isVertex){
            begin = chrono::steady_clock::now();
            fourQEle* Sv=Q.begin()->second;
            double dist_v_w=Q.begin()->first;
            Q.erase(Q.begin());
            c1++;
            Sv->topk=topk_more(P, objCnt+1, dim, k, Sv->closest_point, dg);
            Sv->S_w=no_dominate_set(Sv->topk, P, dim);
            tmp.clear();
            tmp2.clear();
            for(int r: Sv->S_w){
                set_union(doGraph[r].begin(), doGraph[r].end(),
                          tmp.begin(), tmp.end(), back_inserter(tmp2));
                tmp=tmp2;
                tmp2.clear();
            }
            set_union(OneSkybandR.begin(), OneSkybandR.end(),
                      tmp.begin(), tmp.end(), back_inserter(tmp2));
            tmp=Sv->topk;
            sort(tmp.begin(), tmp.end());
            set_difference(tmp2.begin(),  tmp2.end(),
                           tmp.begin(), tmp.end(),
                           back_inserter(Sv->S_NR));
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t6+=elapsed_seconds.count();
            for(int r: Sv->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_v_w<G.topDist()){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
                    G.pop();
                }else if(iter!=G.recordEnd() && iter->second > dist_v_w){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
                }
            }
            if(Sv->topk.size()<=k) {
                delete(Sv);
                continue;
            }
            if(Sv->topk.size()>=2*k){// abnormal
                bool f=false;
                for (int i = 0; i < dim; ++i) {
                    if(Sv->closest_point[i]>1-1e-5){
                        f=true;
                    }
                }
                if(f){
                    delete(Sv);
                    continue;
                }else{
                    cerr<<"too many top-k records, maybe bug in code"<<endl;
                    cout<<Sv->closest_point<<endl;
                    cout<<Sv->topk.size()<<endl;
                    for(int r: Sv->topk){
                        cout<<P[r]*Sv->closest_point<<": ";
                        for (int i = 0; i < dim; ++i) {
                            cout<<P[r][i]<<", ";
                        }
                        cout<<endl;
                    }
//                    exit(1);
                }
            }
            // -------------------------------------------
            double lastScore=P[Sv->topk.back()] * Sv->closest_point ;
            vector<int> fixedTopk;
            vector<int> unfixedTopk;
            for(auto iter=Sv->topk.rbegin(); iter!=Sv->topk.rend(); iter++){
                double tmpScore=P[*iter] * Sv->closest_point;
                if(abs(tmpScore-lastScore)<1e-6){
                    unfixedTopk.push_back(*iter);
                }else{
                    fixedTopk=vector<int>(iter, Sv->topk.rend());
                    break;
                }
            }
            unfixedTopk=approximate_ksksyband(unfixedTopk, Sv->topk, P, dim, k);
            // -------------------------------------------

            sort(Sv->topk.begin(), Sv->topk.end());
            if(Ts.find(Sv->topk)!=Ts.end() || fixedTopk.size()>=k) {
                delete(Sv);
                continue;
            }
            Ts.insert(Sv->topk);

            // select r from n do combination and check total dominate cnt (TDC)
            // return the the combination with maximal TDC
            // This will generate correct result for OSS-skyline, but a brute force solution
            // This practical fast since 1-skyline is small when d is small
            int n=unfixedTopk.size();
            int r=k-fixedTopk.size();
            std::vector<bool> v(n);
            std::fill(v.begin(), v.begin() + r, true);
            do {
                vector<int> S;
//                S.reserve(k);
                S=fixedTopk;
                vector<int> S_com;
                S_com.reserve(n-r);
                for (int i = 0; i < n; ++i) {
                    if (v[i]) {
                        S.push_back(unfixedTopk[i]);
                    }else{
                        S_com.push_back(unfixedTopk[i]);
                    }
                }
                sort(S.begin(), S.end());
                if(Ts.find(S)!=Ts.end()) continue;
                c2++;
                Ts.insert(S);
                auto lbegin=chrono::steady_clock::now();
                auto *Rs=new fourQEle();
                for(int lessRec: S_com) Sv->S_NR.push_back(lessRec);
                double distRS=isFeasibleFOUR(S, P, Rs, dim, Sv, w);
                for(int lessRec: S_com) Sv->S_NR.pop_back();
                auto lnow = chrono::steady_clock::now();
                chrono::duration<double> lelapsed_seconds= lnow-lbegin;
                t7+=lelapsed_seconds.count();
                // Rs, we know Rs contains a point v
                myPoint<double> tmp(Rs->closest_point);
//                if(distRS!=INFINITY && To.find(tmp)==To.end() && Tv.find(tmp)==Tv.end()){ // TODO 这个优化到底是不是正确的
                if(distRS!=INFINITY){
//                    To.insert(tmp);
                    Rs->topk=S;
                    vector<double> innerPoint=(Rs->closest_point+0.001*(Rs->closest_point-w))*0.9;
                    if (!same_kth(Rs->closest_point, Rs->topk, P) && !same_kth(innerPoint, Rs->topk, P)){
                        Q.insert({distRS, Rs});
                    }else{
                        delete(Rs);
                    }
                }else{
                    delete(Rs);
                }
            } while (std::prev_permutation(v.begin(), v.end())); // forward next permutation
            delete(Sv);
        }
        else{
            begin=chrono::steady_clock::now();
            fourQEle* Rs=Q.begin()->second;
            double dist_Rs_w=Q.begin()->first;
            Q.erase(Q.begin());
            for(int r: Rs->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_Rs_w<G.topDist()){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
                    G.pop();
                }else if(iter!=G.recordEnd() && iter->second > dist_Rs_w){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
                }
            }
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t3+=elapsed_seconds.count();
            begin = chrono::steady_clock::now();
            c3++;
            Rs->calGeometry(P, dim, wd);
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t4+=elapsed_seconds.count();
            begin = chrono::steady_clock::now();
            if(Rs->vertexes.size()<dim){
                c4++;
                delete(Rs);
                continue;
            }
            for(auto &v: Rs->vertexes){
                myPoint<double> v_(v);
                if(Tv.find(v_)!=Tv.end()) continue;
                Tv.insert(v_);
                fourQEle *Sv=new fourQEle(true);
//                Sv->topk=shrinkRecTopk(k, P, v, lrtree, lramTree, A, Sv->S_w, Sv->S_NR);
//                sort(Sv->topk.begin(), Sv->topk.end());
//                if(Ts.find(Sv->topk)!=Ts.end()) {
//                    delete(Sv);
//                    continue;
//                }
//                Ts.insert(Sv->topk);
                Sv->closest_point=v;
                double dist_v_w=dist(v, w);
                Q.insert({dist_v_w, Sv});
            }
            now = chrono::steady_clock::now();
            elapsed_seconds= now-begin;
            t5+=elapsed_seconds.count();
            delete(Rs);
        }
    }
    cout<<"t1: "<<t1<<endl;
    cout<<"t2: "<<t2<<endl;
    cout<<"t3: "<<t3<<endl;
    cout<<"t4: "<<t4<<endl;
    cout<<"t5: "<<t5<<endl;
    cout<<"t6: "<<t6<<endl;
    cout<<"t7: "<<t7<<endl;
    cout<<"tq: "<<act_qp_time<<endl;
    cout<<"th: "<<act_qh_time<<endl;
    cout<<"tl: "<<act_ql_time<<endl;
    cout<<"c1: "<<c1<<endl;
    cout<<"c2: "<<c2<<endl;
    cout<<"c3: "<<c3<<endl;
    cout<<"c4: "<<c4<<endl;
    for(auto &e: Q){
        delete(e.second);
    }
    return G.report();
}


vector<vector<int>> buildDominateGraph(float** P, int dim, vector<int> &kskybandR, int objCnt){
    vector<vector<int>> ret(objCnt);
    for(int i: kskybandR){
        for(int j: kskybandR){
            if(i!=j) {
                if (v1_dominate_v2(P[i], P[j], dim)){
                    ret[i].emplace_back(j);
                }
            }
        }
        ret[i]=non_dominate_set(ret[i], P, dim);
        sort(ret[i].begin(), ret[i].end());
    }
    return ret;
}

vector<pair<int, double>> foru2(vector<float>&w,float** P, int dim, int k, int m,
                                vector<vector<int>> &doGraph, vector<int> &kskybandR, vector<int> &OneSkybandR,
                               Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree) {
    // A: k-skyband records and their dominated count
    double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;
    int c1=0, c2=0, c3=0, c4=0;
    act_qh_time=0;
    unordered_map<uint, uint> A;
    auto begin = chrono::steady_clock::now();
    vector<double> wd(w.begin(), w.end());
    fourQEle *Rw=new fourQEle();
    if(kskybandR.size()<200) {
        Rw->topk = computeTopK_Extend(dim, P, kskybandR, w, k, Rw->S_w, Rw->S_NR);
    }else {
        Rw->topk = shrinkRecTopk(k, P, wd, lrtree, lramTree, Rw->S_w, Rw->S_NR);
        vector<int> tmp=OneSkybandR;
        vector<int> tmp2;
        for(int r: Rw->S_w){
            set_union(doGraph[r].begin(), doGraph[r].end(),
                      tmp.begin(), tmp.end(), back_inserter(tmp2));
            tmp=tmp2;
            tmp2.clear();
        }
        tmp2=Rw->topk;
        sort(tmp2.begin(), tmp2.end());
        set_difference(tmp.begin(),  tmp.end(),
                       tmp2.begin(), tmp2.end(),
                       back_inserter(Rw->S_NR));
    }
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    t3+=elapsed_seconds.count();
    begin = chrono::steady_clock::now();
    sort(Rw->S_NR.begin(), Rw->S_NR.end());
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    t4+=elapsed_seconds.count();
    assert(Rw->topk.size()>=k);
    map<int, double> ret;
    for(int record: Rw->topk){
        ret.emplace(record, 0);
    }
    sort(Rw->topk.begin(), Rw->topk.end()); // since Ts is a combination set, hence it needs a sort
    unordered_set<vector<int>, VectorHash<int>> Ts;
    Ts.insert(Rw->topk);

    unordered_set<myPoint<double>, myPointDoubleHash> Tv;
    unordered_set<myPoint<double>, myPointDoubleHash> To;

    // IMPORTANT: check every element that is popped from Q and "delete" them using their pointer
    multimap<double, fourQEle*> Q;
    Rw->closest_point=wd;
    Q.emplace(0, Rw);

    while(!Q.empty() && ret.size()<m) {
//        cout<<"top dist:"<<Q.begin()->first<<","<<Q.begin()->second->isVertex<<"#"<<Q.begin()->second->closest_point<<endl;
        begin=chrono::steady_clock::now();
        fourQEle* Rt=Q.begin()->second;
        double dist_Rs_w=Q.begin()->first;
        Q.erase(Q.begin());
        for(int r: Rt->topk){
            auto iter=ret.find(r);
            if(iter==ret.end()){
                ret.emplace(r, dist_Rs_w);
            }else{
                iter->second=min(iter->second, dist_Rs_w);
            }
        }

        begin = chrono::steady_clock::now();
        vector<pair<uint, uint>> kthRecs = Rt->calGeometry2(P, dim, wd);
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
        t1+=elapsed_seconds.count();
        c1++;
        for (auto &v: kthRecs) {
            fourQEle *Rs = new fourQEle();
            Rs->topk = Rt->topk;
            for (auto &i: Rs->topk) {
                if (i == v.first) {
                    i = v.second;
                    break;
                }
            }
            sort(Rs->topk.begin(), Rs->topk.end());
            if (Ts.find(Rs->topk) != Ts.end()) {
                delete (Rs);
                continue;
            }
            Ts.insert(Rs->topk);
            begin = chrono::steady_clock::now();
            double distRS = isFeasibleFOUR2(Rs->topk, doGraph, P, Rs, dim, Rt, w, v);
            now = chrono::steady_clock::now();
            chrono::duration<double> elapsed_seconds= now-begin;
            t2+=elapsed_seconds.count();
            if(distRS==INFINITY){
                delete(Rs);
                continue;
            }
            c2++;
            Q.insert({distRS, Rs});
        }
        delete(Rt);
    }

    cout<<"t1: "<<t1<<endl;
    cout<<"t2: "<<t2<<endl;
    cout<<"t3: "<<t3<<endl;
    cout<<"t4: "<<t4<<endl;
    cout<<"t5: "<<t5<<endl;
    cout<<"t6: "<<t6<<endl;
    cout<<"t7: "<<t7<<endl;
    cout<<"tq: "<<act_qp_time<<endl;
    cout<<"th: "<<act_qh_time<<endl;
    cout<<"tl: "<<act_ql_time<<endl;
    cout<<"c1: "<<c1<<endl;
    cout<<"c2: "<<c2<<endl;
    cout<<"c3: "<<c3<<endl;
    cout<<"c4: "<<c4<<endl;
    for(auto &e: Q){
        delete(e.second);
    }
    vector<pair<int, double>> tmp(ret.begin(), ret.end());
    sort(tmp.begin(), tmp.end(), [](auto &a, auto &b){
       return a.second < b.second;
    });
    return tmp;
}

vector<pair<int, double>> forut(vector<float>&w,float** P, int dim, int k, int m,
                               unordered_map<uint, uint> &A) {
    // A: k-skyband records and their dominated count

//    template<typename INT, typename VV>
//    vector<INT> computeTopK(const int dim, VV &PG, vector<INT> &skyband, vector<float>& weight, int k);
//    template<typename VV>
//    vector<int> k_skyband(VV &P, const int &k);
    vector<int> skybandR=k_skyband(P, k, objCnt, dim);
    vector<double> wd(w.begin(), w.end());
    fourQEle Rw;
    Rw.topk =computeTopK_Extend(dim, P, skybandR, w, k, Rw.S_w, Rw.S_NR);
    sort(Rw.topk.begin(), Rw.topk.end()); // since Ts is a combination set, hence it needs a sort
    unordered_set<vector<int>, VectorHash<int>> Ts;
    Ts.insert(Rw.topk);
    myHeap G;
    for(int record: Rw.topk){
        G.insert(record, 0);
    }
    while(G.size()>m){
        G.pop();
    }
    while(G.size()<m){
        G.insert(-(int)G.size(), INFINITY);
    }

    unordered_set<myPoint<double>, myPointDoubleHash> Tv;
    unordered_set<myPoint<double>, myPointDoubleHash> To;

    // IMPORTANT: check every element that is popped from Q and "delete" them using their pointer
    multimap<double, fourQEle*> Q;
    Rw.closest_point=wd;
    Rw.calGeometry(P, dim, wd);
//    cout<<"debug, Rw.vertexes.size(): "<<Rw.vertexes.size()<<endl;
//    auto begin = chrono::steady_clock::now();
    for(auto &v: Rw.vertexes){
        Tv.emplace(v);
        fourQEle *Sv=new fourQEle(true);
//        Sv->topk=shrinkRecTopk(k, P, v, lrtree, lramTree, A, Sv->S_w, Sv->S_NR);
//        sort(Sv->topk.begin(), Sv->topk.end());
//        if(Ts.find(Sv->topk)!=Ts.end()) {
//            delete(Sv);
//            continue;
//        }
//        Ts.insert(Sv->topk);
        Sv->closest_point=v;
        double dist_v_w=dist(v, w);
        Q.insert({dist_v_w, Sv});
    }
//    auto now = chrono::steady_clock::now();
//    chrono::duration<double> elapsed_seconds= now-begin;
//    cout<<"debug, Rw top-k time: "<<elapsed_seconds.count()<<endl;
//    double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;
//    int c1=0, c2=0, c3=0, c4=0;
    while(!Q.empty() && G.topDist()>Q.begin()->first){
//        cout<<"top dist:"<<Q.begin()->first<<","<<Q.begin()->second->isVertex<<"#"<<Q.begin()->second->closest_point<<endl;
        if(Q.begin()->second->isVertex){
//            begin = chrono::steady_clock::now();
            fourQEle* Sv=Q.begin()->second;
            double dist_v_w=Q.begin()->first;
            Q.erase(Q.begin());
//            c1++;
            Sv->topk=computeTopK_Extend(dim, P, skybandR, Sv->closest_point, k, Sv->S_w, Sv->S_NR);
            for(int r: Sv->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_v_w<G.topDist()){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
                    G.pop();
                }else if(iter!=G.recordEnd() && iter->second > dist_v_w){
//                    cout<<Sv->isVertex<<","<<r<<", "<<dist_v_w<<endl;
                    G.insert(r, dist_v_w);
                }
            }
            if(Sv->topk.size()<=k) {
                delete(Sv);
                continue;
            }
            // -------------------------------------------
            double lastScore=P[Sv->topk.back()] * Sv->closest_point ;
            vector<int> fixedTopk;
            vector<int> unfixedTopk;
            for(auto iter=Sv->topk.rbegin(); iter!=Sv->topk.rend(); iter++){
                double tmpScore=P[*iter] * Sv->closest_point;
                if(abs(tmpScore-lastScore)<1e-6){
                    unfixedTopk.push_back(*iter);
                }else{
                    fixedTopk=vector<int>(iter, Sv->topk.rend());
                    break;
                }
            }
            // -------------------------------------------

            sort(Sv->topk.begin(), Sv->topk.end());
            if(Ts.find(Sv->topk)!=Ts.end() || fixedTopk.size()>=k) {
                delete(Sv);
                continue;
            }
            Ts.insert(Sv->topk);
//            now = chrono::steady_clock::now();
//            elapsed_seconds= now-begin;
//            t6+=elapsed_seconds.count();
            // select r from n do combination and check total dominate cnt (TDC)
            // return the the combination with maximal TDC
            // This will generate correct result for OSS-skyline, but a brute force solution
            // This practical fast since 1-skyline is small when d is small
//            begin=chrono::steady_clock::now();

            int n=unfixedTopk.size();
            int r=k-fixedTopk.size();
            assert(r>=1);
            std::vector<bool> v(n);
            std::fill(v.begin(), v.begin() + r, true);
            do {
                vector<int> S;
//                S.reserve(k);
                S=fixedTopk;
                vector<int> S_com;
                S_com.reserve(n-r);
                for (int i = 0; i < n; ++i) {
                    if (v[i]) {
                        S.push_back(unfixedTopk[i]);
                    }else{
                        S_com.push_back(unfixedTopk[i]);
                    }
                }
                sort(S.begin(), S.end());
                if(Ts.find(S)!=Ts.end()) continue;
//                c2++;
                Ts.insert(S);
//                auto begin = chrono::steady_clock::now();
//                auto lbegin=chrono::steady_clock::now();
                auto *Rs=new fourQEle();
                for(int lessRec: S_com) Sv->S_NR.push_back(lessRec);
                double distRS=isFeasibleFOUR(S, P, Rs, dim, Sv, w);
                for(int lessRec: S_com) Sv->S_NR.pop_back();
//                auto lnow = chrono::steady_clock::now();
//                chrono::duration<double> lelapsed_seconds= lnow-lbegin;
//                t7+=lelapsed_seconds.count();
//                auto now = chrono::steady_clock::now();
//                chrono::duration<double> elapsed_seconds= now-begin;
//                cout<<"debug, foru qp time: "<<Rs->S_w.size()<<", "<<Rs->S_NR.size()<<", "<<elapsed_seconds.count()<<endl;

                // Rs, we know Rs contains a point v
                myPoint<double> tmp(Rs->closest_point);
//                if(distRS!=INFINITY && To.find(tmp)==To.end() && Tv.find(tmp)==Tv.end()){ // TODO 这个优化到底是不是正确的
                if(distRS!=INFINITY){
//                    To.insert(tmp);
                    Rs->topk=S;
                    Q.insert({distRS, Rs});
                }else{
                    delete(Rs);
                }
            } while (std::prev_permutation(v.begin(), v.end())); // forward next permutation

//            int r=k;
//            int n=Sv->topk.size();
//            std::vector<bool> v(n);
//            std::fill(v.begin(), v.begin() + r, true);
//            do {
//                vector<int> S;
//                S.reserve(k);
//                S=fixedTopk;
//                vector<int> S_com;
//                S_com.reserve(n-r);
//                for (int i = 0; i < n; ++i) {
//                    if (v[i]) {
//                        S.push_back(Sv->topk[i]);
//                    }else{
//                        S_com.push_back(Sv->topk[i]);
//                    }
//                }
//                sort(S.begin(), S.end());
//                if(Ts.find(S)!=Ts.end()) continue;
//                c2++;
//                Ts.insert(S);
////                auto begin = chrono::steady_clock::now();
//                auto lbegin=chrono::steady_clock::now();
//                auto *Rs=new fourQEle();
//                for(int lessRec: S_com) Sv->S_NR.push_back(lessRec);
//                double distRS=isFeasibleFOUR(S, A, P, Rs, dim, lrtree, lramTree, Sv, w);
//                for(int lessRec: S_com) Sv->S_NR.pop_back();
//                auto lnow = chrono::steady_clock::now();
//                chrono::duration<double> lelapsed_seconds= lnow-lbegin;
//                t7+=lelapsed_seconds.count();
////                auto now = chrono::steady_clock::now();
////                chrono::duration<double> elapsed_seconds= now-begin;
////                cout<<"debug, foru qp time: "<<Rs->S_w.size()<<", "<<Rs->S_NR.size()<<", "<<elapsed_seconds.count()<<endl;
//
//                // Rs, we know Rs contains a point v
//                myPoint<double> tmp(Rs->closest_point);
////                if(distRS!=INFINITY && To.find(tmp)==To.end() && Tv.find(tmp)==Tv.end()){ // TODO 这个优化到底是不是正确的
//                if(distRS!=INFINITY){
//                    To.insert(tmp);
//                    Rs->topk=S;
//                    Q.insert({distRS, Rs});
//                }else{
//                    delete(Rs);
//                }
//            } while (std::prev_permutation(v.begin(), v.end())); // forward next permutation
//            now = chrono::steady_clock::now();
//            elapsed_seconds= now-begin;
//            t2+=elapsed_seconds.count();
//            cout<<"debug, qp count "<<debug_qp<<endl;
            delete(Sv);
        }
        else{
//            begin=chrono::steady_clock::now();
            fourQEle* Rs=Q.begin()->second;
            double dist_Rs_w=Q.begin()->first;
            Q.erase(Q.begin());
            for(int r: Rs->topk){
                auto iter=G.findRecord(r);
                if(iter==G.recordEnd() && dist_Rs_w<G.topDist()){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
                    G.pop();
                }else if(iter!=G.recordEnd() && iter->second > dist_Rs_w){
//                    cout<<Rs->isVertex<<","<<r<<", "<<dist_Rs_w<<endl;
                    G.insert(r, dist_Rs_w);
                }
            }
//            now = chrono::steady_clock::now();
//            elapsed_seconds= now-begin;
//            t3+=elapsed_seconds.count();
//            begin = chrono::steady_clock::now();
//            c3++;
            Rs->calGeometry(P, dim, wd);
//            now = chrono::steady_clock::now();
//            elapsed_seconds= now-begin;
//            t4+=elapsed_seconds.count();
//            begin = chrono::steady_clock::now();
            if(Rs->vertexes.size()<dim){
//                c4++;
                delete(Rs);
                continue;
            }
            for(auto &v: Rs->vertexes){
                myPoint<double> v_(v);
                if(Tv.find(v_)!=Tv.end()) continue;
                Tv.insert(v_);
                fourQEle *Sv=new fourQEle(true);
//                Sv->topk=shrinkRecTopk(k, P, v, lrtree, lramTree, A, Sv->S_w, Sv->S_NR);
//                sort(Sv->topk.begin(), Sv->topk.end());
//                if(Ts.find(Sv->topk)!=Ts.end()) {
//                    delete(Sv);
//                    continue;
//                }
//                Ts.insert(Sv->topk);
                Sv->closest_point=v;
                double dist_v_w=dist(v, w);
                Q.insert({dist_v_w, Sv});
            }
//            now = chrono::steady_clock::now();
//            elapsed_seconds= now-begin;
//            t5+=elapsed_seconds.count();
            delete(Rs);
        }
    }
//    cout<<"t1: "<<t1<<endl;
//    cout<<"t2: "<<t2<<endl;
//    cout<<"t3: "<<t3<<endl;
//    cout<<"t4: "<<t4<<endl;
//    cout<<"t5: "<<t5<<endl;
//    cout<<"t6: "<<t6<<endl;
//    cout<<"t7: "<<t7<<endl;
//    cout<<"tq: "<<act_qp_time<<endl;
//    cout<<"c1: "<<c1<<endl;
//    cout<<"c2: "<<c2<<endl;
//    cout<<"c3: "<<c3<<endl;
//    cout<<"c4: "<<c4<<endl;
    return G.report();
}

int topRegions_efficient3(vector<vector<double>> &parent_region, ch &ch_obj,
                          multimap<double, topi_region*> &id_radius, float **PG, int dim, int X,
                          const int k, vector<float> &w, unordered_set<int> &top1_calculated,
                          Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree,
                          vector<pair<int, double>> &utk_option_ret,
                          vector<pair<double, topi_region*>> &utk_cones_ret,
                          unordered_map<int, vector<int>> &top1_region, double delta){
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
        vector<int> topk= popped.second->get_rtopi();
        for(int opt:topk){
            auto iter=options.find(opt);
            if(iter==options.end()) {// new option
                options.insert(opt);
                utk_option_ret.emplace_back(opt, popped.first);
                cout << "radius: " << popped.first << "\n";
                if(true){
                    vector<double> tmp(PG[opt], PG[opt]+dim);
                    cout << options.size() << ": " << popped.first << ", " << utk_cones_ret.size() <<
                         ", " << rcnt << "," << ch_obj.get_neighbor_vertex(opt).size() << " # " << popped.second->get_rtopi() << "\n";
                }
                auto now = chrono::steady_clock::now();
                chrono::duration<double> elapsed_seconds= now-begin;
                cout << "time: " << options.size() << ", " << opt <<"," << elapsed_seconds.count()<<"\n";
                if(GEN_F1) myfile << "time: " << options.size() << ", " << opt <<"," << elapsed_seconds.count()<<"\n";
                new_option=True;
            }
        }
        if(new_option){
            popped.second->set_radius(popped.first);
        }
        if(popped.second->top_what==k){ // a region that don't need to be partitioned
            vector<int> topk= popped.second->get_rtopi();
            for(int opt:topk){
                auto iter=options.find(opt);
                if(iter==options.end()) {// new option
                    options.insert(opt);
                    utk_option_ret.emplace_back(opt, popped.first);
                    cout << "radius: " << popped.first << "\n";
                    if(true){
                        vector<double> tmp(PG[opt], PG[opt]+dim);
                        cout << options.size() << ": " << popped.first << ", " << utk_cones_ret.size() <<
                             ", " << rcnt << "," << ch_obj.get_neighbor_vertex(opt).size() << " # " << popped.second->get_rtopi() << "\n";
                    }
                    auto now = chrono::steady_clock::now();
                    chrono::duration<double> elapsed_seconds= now-begin;
                    cout << "time: " << options.size() << ", " << opt <<"," << elapsed_seconds.count()<<"\n";
                    if(GEN_F1) myfile << "time: " << options.size() << ", " << opt <<"," << elapsed_seconds.count()<<"\n";
                    new_option=True;
                }
            }
            if(new_option){
                popped.second->set_radius(popped.first);
            }
        }
        else{
            // if main diagonal less than delta
            if(delta>0) {
                vector<double> center_mbb(dim);
                double md;
                if(popped.second->top_what==1){
                    md=popped.second->main_diagonal(PG, dim, ch_obj.get_neighbor_vertex(popped.second->opt_i), center_mbb);
                }else{
                    md=popped.second->main_diagonal(PG, dim, popped.second->parent->child_constraint, center_mbb);
                }
                if ( md< delta) {
                    // center of mmb
//                    cout<<"cut\n";
                    vector<int> topk=topk_single(k, PG, center_mbb, 0, rtree, lramTree);
                    for(int opt:topk){
                        auto iter=options.find(opt);
                        if(iter==options.end()) {// new option
                            options.insert(opt);
                            utk_option_ret.emplace_back(opt, popped.first);
                            cout << "radius: " << popped.first << "\n";
                            if(true){
                                vector<double> tmp(PG[opt], PG[opt]+dim);
                                cout << options.size() << ": " << popped.first << ", " << utk_cones_ret.size() <<
                                     ", " << rcnt << "," << ch_obj.get_neighbor_vertex(opt).size() << " # " << topk << "\n";
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
                    continue;
                }
            }
            vector<int> topi= popped.second->get_rtopi();
            unordered_set<int> ch_upd_s;
            int m=0; // find ma[ximal convex hull layer
            // find all option adjacent to current top 1 to top i options
            for(int top: topi){
                for(int adj: ch_obj.get_neighbor_vertex(top)){
                    ch_upd_s.insert(adj);
                }
                m=max(ch_obj.get_option_layer(top), m);
            }
            // find the maximal layer to convex points
            for(int mp1_opt: ch_obj.get_layer(m+1)){
                ch_upd_s.insert(mp1_opt);
            }
            // erase options for top 1 to top i
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
            // fetch those options non-dominated by the top-1 to top-i
            vector<int> ch_upd;
            for (int i:ch_upd_s) {
                auto iter=dominated_cnt.find(i);
                if(iter!=dominated_cnt.end()&&iter->second<=0){
                    ch_upd.push_back(i);
                }
            }
            // parent_region is no longer used
            //
            vector<int> new_ch_upd(ch_upd);

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

int utk_efficient3(float **PointSet, int dim, vector<float> &w,
        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree,
        int X, int k, vector<pair<int, double>> &utk_option_ret,
                   vector<pair<double, region*>> &utk_cones_ret, double delta){

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
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, ram_Tree, PointSet);
    pair<int, float> next={-1, INFINITY};
    cout<< "begin fetch CH1"<<endl;
    vector<int> CH_1_X_opt;
    int top1=0;
    bool top1f= false;
    while(CH_1_X_opt.size() < X){
        next=get_next_obj.get_next();
        if(next.second==INFINITY){
            break;
        }
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
        if(minc>maxc){
            break;
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
    computeRho(dim, k, INFINITY, w, *rtree, ram_Tree, PointSet, interval, rho_star);
    vector<int> rskyband_CS;
    for (pair<long int, float> &p:interval) {
        rskyband_CS.push_back(p.first);
    }
    // 4. end: apply \rho^* to get r-k-skyband
    cout<< elapsed_seconds.count() << " rskyband size: "<<rskyband_CS.size()<< "\n";

    vector<vector<double>> tmp;
    // 5. begin: apply \rho^* to get initial domain

    cout<< elapsed_seconds.count() << " begin generate domain" << endl;
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
                                  k,  w, top1_calculated, rtree, ram_Tree, utk_option_ret, utk_cones_ret_, top1_region, delta);
    delete(root);
    delete(top1r);
    // 8. end: apply ORU, that is recursively divide regions
    now = chrono::steady_clock::now();
    elapsed_seconds= now-begin;
    cout<< elapsed_seconds.count();
    cout<< " finish recursively get top regions\n";
    return ret;
}


int utk_efficient(float **PointSet, int dim, vector<float> &w,
        Rtree* rtree, unordered_map<long int, RtreeNode*> &a_ramTree,
        int X, int k, vector<pair<int, double>> &utk_option_ret,
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
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, a_ramTree, PointSet);
    pair<int, float> next={-1, INFINITY};
    cout<< "begin fetch CH1"<<endl;
    vector<int> CH_1_X_opt;
    int top1=0;
    bool top1f= false;
    while(CH_1_X_opt.size() < X){
        next=get_next_obj.get_next();
        if(next.second==INFINITY){
            break;
        }
        cout<<get_next_obj.interval.size()<<" "<<next.second<<endl;
        if(next.first<0) break;
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
        if(minc>maxc){
            break;
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
    computeRho(dim, k, INFINITY, w, *rtree, a_ramTree, PointSet, interval, rho_star);
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


int utk_efficient_cs3(float **PointSet, int dim, vector<float> &w,
        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree,
        int X, int k, vector<pair<int, double>> &utk_option_ret,
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
    unknown_x_efficient get_next_obj(dim, 1, w, *rtree, ram_Tree, PointSet);
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
double act_qp_time=0;
double act_qh_time=0;
double act_ql_time=0;

double non_order_feasible(vector<int> &gt, vector<int> &le, float ** PG, vector<float> &w, vector<double> &retX){
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
    auto begin = chrono::steady_clock::now();
    double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    auto now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds= now-begin;
    act_qp_time+=elapsed_seconds.count();

    if(solver_ret == std::numeric_limits<double>::infinity()){
        return INFINITY;
    }else{
        retX.resize(x.size());
        for (int i = 0; i < n; ++i) {
            retX[i]=x[i];
        }
        double ret=0;
        for (int k = 0; k < n; ++k) {
            ret+=(x[k]-w[k])*(x[k]-w[k]);
        }
        return sqrt(ret);
    }
}



