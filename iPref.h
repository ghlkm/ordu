#ifndef _IPREF_H_
#define _IPREF_H_

#include "rtree_lib/header.h"
#include "rtree_lib/rtree.h"
#include "rtree_lib/rnode.h"
#include "rtree_lib/rentry.h"
#include "rtree_lib/virtualRNode.h"
#include "rtree_lib/filemem.h"
//#include "global.h"
#include "skyline.h"
#include "qhull_user.h"
#include "region.h"
#include "vector_operator.h"
#include "math_lib.h"
#define GEN_F1 false
//extern unordered_map<long int, RtreeNode*> ramTree;
extern Rtree* rtree;
extern double act_qp_time;
extern double act_qh_time;
extern double act_ql_time;

// compute the hyperplane of incomparable record pair: pi and pj
vector<float> computePairHP(const int dimen, float* PG[], long int pi, long int pj);

// compute the distance from w to hyperplane HS
float computeDis(const vector<float> &tmpHS, const vector<float> &userpref);

// optimized algorithm
float computeRho(const int dimen, const int k, const int X, vector<float>& userpref,
        Rtree& a_rtree, unordered_map<long int, RtreeNode*> &a_ramTree, float* PG[],
                 vector<pair<long int, float>> &interval, float radius = INFINITY);


// compute the radius in phase 2, optimized algorithm
float computeradius(const int k, const int dim, long int pi, vector<float>& userpref, vector<long int>& incompSet, float* PG[], float rho=INFINITY);


class unknown_X_node{
private:
    int dim;
public:
    multimap<float, vector<double>> topK_dominate_radius;// size with k
    int tau;// used in unknown efficient
    int needed_to_update; // need how many options to update this node's radius
    const float *data;  // the values
    int page_id;       // rtree pageid
    bool fetched;
    unknown_X_node(int d,  const float *dt, int pid);
    float radius_LB();
    vector<double> radius_LB_point();

    void init_radius(vector<long int> &fetched_options, float **PointSet, vector<float> &w, int k, float rho=INFINITY);

    inline bool update(int need_to_update) const;

    inline bool update() const;

    void update_radius(vector<long int>::iterator begin, vector<long int>::iterator end, float **PointSet, vector<float> &w, float rho=INFINITY);

    void update_radius(const float* other_option, vector<float> &w);

    void update_radius_erase(const float* other_option, vector<float> &w);
};

class unknown_x_efficient {
    typedef float RADIUS;  // is the inflection radius
    typedef int PAGE_ID;  // page id in rtree
    typedef int PG_IDX;  // idx in PG
    typedef long PG_IDX_L;  // idx in PG
    typedef float DOT;   // is the result of dot product with w
    typedef unknown_X_node *NODE_PTR;

    multimap<DOT, NODE_PTR, greater<float>> heap; //BBS max_heap, <w\cdot node, node>, if node is an option then node=id+MAXPAGEID
    unordered_set<NODE_PTR> S; // reflect which nodes and options in BBS
    multimap<RADIUS, NODE_PTR, less<float>> Q; // min_heap based on inflection radius, lazy update for S
    multimap<RADIUS, NODE_PTR, less<float>> C; // min_heap based on inflection radius, candidate list
    vector<PG_IDX_L> incompSet;
    float pt[MAXDIMEN];
    int dimen;
    int k;
    vector<float> userpref;
    Rtree &a_rtree;
    unordered_map<long int, RtreeNode*> &a_ramTree;
    float **PG;
public:
    vector<pair<PG_IDX, RADIUS>> interval; // return
    vector<vector<double>> drill_position;

//    unknown_x_efficient(const int dim, int K, vector<float> &userPref, Rtree &aRtree, float **pg);

    unknown_x_efficient(const int dim, int K, vector<float> &userPref,
            Rtree &aRtree, unordered_map<long int, RtreeNode*> &ram_Tree, float **pg);

    pair<int, float> get_next();
    ~unknown_x_efficient();
};


int utk_efficient(float **PointSet, int dim, vector<float> &w,
        Rtree* rtree, unordered_map<long int, RtreeNode*> &a_ramTree,
        int X, int k, vector<pair<int, double>> &utk_option_ret,
                  vector<pair<double, region*>> &utk_cones_ret);

vector<vector<double>> points_to_halfspace(vector<vector<double>> &points);

int utk_efficient_anti(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                       vector<pair<int, double>> &utk_option_ret,
                       vector<pair<double, region*>> &utk_cones_ret);

int utk_efficient_cs3(float **PointSet, int dim, vector<float> &w,
        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree,
        int X, int k, vector<pair<int, double>> &utk_option_ret,
                      vector<pair<double, region*>> &utk_cones_ret, double &rho_star);

int utk_efficient3(float **PointSet, int dim, vector<float> &w,
        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree,
        int X, int k, vector<pair<int, double>> &utk_option_ret,
                   vector<pair<double, region*>> &utk_cones_ret, double delta=-1);

int non_order_utk_efficient(float **PointSet, int pt_cnt, int dim, vector<float> &w,
                            Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree,
                            int X, int k, vector<pair<int, double>> &utk_option_ret,
                            vector<pair<double, region*>> &utk_cones_ret);




template<typename VVT>
vector<int> topk_single_extend(int k,  VVT &P, const vector<double>&w, double uHeat,
                               Rtree *lrtree_rt, unordered_map<long int, RtreeNode *> &lramTree){
    RtreeNode* node;
    priority_queue<pair<double, int>> heap;
    double bound=uHeat;
    heap.emplace(INFINITY, lrtree_rt->m_memory.m_rootPageID);
    int dim=w.size();
    vector<double> tmp_v(dim, 0.0);
    vector<int> topk;
    double tmp_score;
    long pageID;
    double last_score=INFINITY;
    while(!heap.empty() && (topk.size()<k || heap.top().first>=last_score-1e-6)){ // 1e-6 for accuracy
        tmp_score=heap.top().first;
        last_score=tmp_score;
        pageID=heap.top().second;
        heap.pop();
        if (pageID >= MAXPAGEID){ // an option
            topk.push_back(pageID-MAXPAGEID);
        }else{ // an inter node
            node = lramTree[pageID];
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

// I need a dominance graph top-k retrieval method
// before that, I need to build the dominance graph
// building a dominance graph is easy
// while building the dominance graph, I need to cut those indirect edges in dominance graph

// graph:
// each node should know its incoming edges

// store the dominance graph as an adjacency list
template<typename VVF>
vector<vector<size_t>> build_doInGraph(VVF& D, size_t Dsize, int dim, int k,
        Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree,
        vector<int> &oneSkyband){
    // build kskyband and dominance graph simultaneously
    RtreeNode* node;
    multimap<float, int, greater<>> heap;
    multimap<float, int, greater<>>::iterator heapIter;
    float pt[MAXDIMEN];
    float ORIGNIN[MAXDIMEN];
    float mindist;
    for (int i = 0; i < dim; i++)
        ORIGNIN[i] = 1;
    int pageID;
    float dist_tmp;
    heap.emplace(INFINITY, lrtree->m_memory.m_rootPageID);
    vector<vector<size_t>> inGraph(Dsize+1);
    vector<size_t> kskyband;
    while (!heap.empty()){
        heapIter = heap.begin();
        dist_tmp = heapIter->first;
        pageID = heapIter->second;
        heap.erase(heapIter);
        if (pageID > MAXPAGEID){
            for (int d = 0; d < dim; d++) {
                pt[d] = (D[pageID - MAXPAGEID][d] + D[pageID - MAXPAGEID][d + dim]) / 2;
            }
            assert(pageID - MAXPAGEID <inGraph.size());
            if (!dominatedByK(dim, pt, kskyband, D, k, inGraph[pageID-MAXPAGEID])){

                // prune indirect edges
                set<size_t> grandparent;
                for(size_t inVtx: inGraph[pageID-MAXPAGEID]){
//                    cout<<"debug: "<<inVtx<<endl;
                    assert(inVtx<inGraph.size());
                    grandparent.insert(inGraph[inVtx].begin(), inGraph[inVtx].end());
                }
                set<size_t> inVtx(inGraph[pageID-MAXPAGEID].begin(), inGraph[pageID-MAXPAGEID].end());
                inGraph[pageID-MAXPAGEID].clear();
                set_difference(inVtx.begin(),  inVtx.end(),
                               grandparent.begin(), grandparent.end(),
                               back_inserter(inGraph[pageID-MAXPAGEID]));

                if(inGraph[pageID-MAXPAGEID].empty()){
                    oneSkyband.push_back(pageID-MAXPAGEID);
                }
                kskyband.push_back(pageID - MAXPAGEID);
            }
        }else{
            node = lramTree[pageID];
            if (node->isLeaf()){
                for (int i = 0; i < node->m_usedspace; i++){
                    for (int d = 0; d < dim; d++){
                        pt[d] = (node->m_entry[i]->m_hc.getLower()[d] + node->m_entry[i]->m_hc.getUpper()[d])/2;
                    }
                    mindist = dot(pt, ORIGNIN, dim);
//                    mindist = minDist(pt, ORIGNIN, dim);
                    heap.emplace(mindist, node->m_entry[i]->m_id + MAXPAGEID);
                }
            }else{
                for (int i = 0; i < node->m_usedspace; i++){
                    for (int d = 0; d < dim; d++)
                        pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
                    if (!dominatedByK(dim, pt, kskyband, D, k)){
//                        mindist = minDist(pt, ORIGNIN, dim);
                        mindist = dot(pt, ORIGNIN, dim);
                        heap.emplace(mindist, node->m_entry[i]->m_id);
                    }
                }
            }
        }
    }
    return inGraph;
}

template<typename VVF>
vector<vector<size_t>> build_doInGraph2(VVF& D, size_t Dsize, int dim, int k,
                                       vector<int> &kskyband,
                                       vector<int> &oneSkyband){
    // build kskyband and dominance graph simultaneously
    RtreeNode* node;
    multimap<float, int, greater<>> heap;

    multimap<float, int, greater<>>::iterator heapIter;
    float pt[MAXDIMEN];
    float ORIGNIN[MAXDIMEN];
    float mindist;
    for (int i = 0; i < dim; i++)
        ORIGNIN[i] = 1;
    for(int r: kskyband){
        heap.insert({dot(ORIGNIN, D[r], dim), r});
    }
    kskyband.clear();
    int pageID;
    float dist_tmp;
    vector<vector<size_t>> inGraph(Dsize+1);
    while (!heap.empty()){
        heapIter = heap.begin();
        dist_tmp = heapIter->first;
        pageID = heapIter->second;
        heap.erase(heapIter);
        if (!dominatedByK(dim, D[pageID], kskyband, D, k, inGraph[pageID])){
            // prune indirect edges
            set<size_t> grandparent;
            for(size_t inVtx: inGraph[pageID]){
//                    cout<<"debug: "<<inVtx<<endl;
                assert(inVtx<inGraph.size());
                grandparent.insert(inGraph[inVtx].begin(), inGraph[inVtx].end());
            }
            set<size_t> inVtx(inGraph[pageID].begin(), inGraph[pageID].end());
            inGraph[pageID].clear();
            set_difference(inVtx.begin(),  inVtx.end(),
                           grandparent.begin(), grandparent.end(),
                           back_inserter(inGraph[pageID]));
//            inVtx.erase(grandparent.begin(), grandparent.end());
//            inGraph[pageID]=vector<size_t>(inVtx.begin(), inVtx.end());

            if(inGraph[pageID].empty()){
                oneSkyband.push_back(pageID);
            }
            kskyband.push_back(pageID);
        }
    }
    return inGraph;
}

class doGraphPlus{
public:
    vector<vector<size_t>> inGraph;
    vector<vector<size_t>> outGraph;
    vector<int> oneSkyband;
    doGraphPlus()=default;

    doGraphPlus(float** D, size_t Dsize, int dim, int k,
                Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree){
        inGraph=build_doInGraph(D, Dsize, dim, k, lrtree, lramTree, oneSkyband);
        inGraph2outGraph();
        sort(oneSkyband.begin(), oneSkyband.end());
        for(auto &i: inGraph){
            sort(i.begin(), i.end());
        }
        for(auto &i: outGraph){
            sort(i.begin(), i.end());
        }
    }

    void inGraph2outGraph(){
        outGraph.resize(inGraph.size());
        for (int outVtx = 0; outVtx < inGraph.size(); ++outVtx) {
            for (int inVtx: inGraph[outVtx]) {
                outGraph[inVtx].push_back(outVtx);
            }
        }

    }

};

class doGraphPlus2: public doGraphPlus{
public:
//    vector<vector<size_t>> inGraph;
//    vector<vector<size_t>> outGraph;
//    vector<int> oneSkyband;

    doGraphPlus2(float** D, size_t Dsize, int dim, int k,
                vector<int> &kskyband){
        inGraph=build_doInGraph2(D, Dsize, dim, k, kskyband, oneSkyband);
        inGraph2outGraph();
        sort(oneSkyband.begin(), oneSkyband.end());
        for(auto &i: inGraph){
            sort(i.begin(), i.end());
        }
        for(auto &i: outGraph){
            sort(i.begin(), i.end());
        }
    }

    void inGraph2outGraph(){
        outGraph.resize(inGraph.size());
        for (int outVtx = 0; outVtx < inGraph.size(); ++outVtx) {
            for (int inVtx: inGraph[outVtx]) {
                outGraph[inVtx].push_back(outVtx);
            }
        }

    }

    vector<vector<size_t>> multihop(size_t hop){
        vector<vector<size_t>> ret(hop);
        if(hop==0) return ret;
        for(size_t r: oneSkyband){
            ret[0].push_back(r);
        }
        for (int i = 0; i < hop-1; ++i) {
            for (size_t r: ret[i]) {
                ret[i+1].push_back(r);
            }
        }
        return ret;
    }

};

vector<pair<int, double>> drill_by_grid(vector<float>&w, float** P, int dim, int k, int m, double delta,
                                        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree, unordered_map<int, double> &p1_options);

vector<pair<int, double>> drill_by_grid(vector<float>&w, float** P, int dim, int k, int m, double delta,
                                        Rtree* rtree, unordered_map<long int, RtreeNode*> &ram_Tree, unordered_map<int, double> &p1_options,
                                        doGraphPlus &dg);
//template<typename VF1, typename VF2>
//inline float dot(VF1& v1, VF2& v2, size_t size){
//    float ret=0;
//    for (int i = 0; i < size; ++i) {
//        ret+=v1[i]*v2[i];
//    }
//    return ret;
//}

inline bool findValue(multimap<float, size_t> &heap, size_t v){
    for (auto iter=heap.begin(); iter!=heap.end(); ++iter) {
        if(iter->second==v) return true;
    }
    return false;
}


// a top-k retrieval algorithm using dominance graph
// fetch the non-dominance records' topk first
// then the records that only dominated by top-k non-dominated records
// the stop condition should be no new record update top-(k-1)
template<typename VVF, typename F>
vector<int> topk(VVF& D, size_t Dsize, int dim, int k, vector<F> &w, doGraphPlus &dg){
    vector<int>& oneSkyband=dg.oneSkyband;
    vector<vector<size_t>> &inGraph=dg.inGraph;
    vector<vector<size_t>> &outGragh=dg.outGraph;
    multimap<float, size_t> heap;
    if(oneSkyband.size()<k){
        for (int i = 0; i < oneSkyband.size(); ++i) {
            float score=dot(D[oneSkyband[i]], w, dim);
            heap.insert({score, oneSkyband[i]});
        }
        while(heap.size()<k){
            heap.insert({0, 1});
        }
    }else {
        for (int i = 0; i < k; ++i) {
            float score = dot(D[oneSkyband[i]], w, dim);
            heap.insert({score, oneSkyband[i]});
        }
    }
    for (int i = k; i < oneSkyband.size(); ++i) {
        float score=dot(D[oneSkyband[i]], w, dim);
        float tScore=heap.begin()->first;
        if(score>tScore){
            heap.erase(heap.begin());
            heap.insert({score, oneSkyband[i]});
        }
    }

    for (int i = 0; i < k-1; ++i) {
        auto iter=heap.rbegin();
        for (int j = 0; j < i; ++j) {
            iter++;
        }
        size_t inVtx=iter->second;
        float tScore=heap.begin()->first;

        for(size_t outVtx: outGragh[inVtx]){
            float score=dot(D[outVtx], w, dim);
            if(score>tScore){
                // see if this record is in heap
                if(!findValue(heap, outVtx)){
                    heap.erase(heap.begin());
                    heap.insert({score, outVtx});
                    tScore=heap.begin()->first;
                }
            }
        }
    }
    vector<int> ret;
    ret.reserve(k);
    for (auto iter=heap.rbegin(); iter!=heap.rend(); ++iter) {
        ret.push_back(iter->second);
    }
    return ret;
}

// For this version topk, I have to consider that retrievals more than k topk records.
template<typename VVF, typename F>
vector<int> topk_more(VVF& D, size_t Dsize, int dim, int k, const vector<F> &w, doGraphPlus &dg){
    vector<int>& oneSkyband=dg.oneSkyband;
    vector<vector<size_t>> &inGraph=dg.inGraph;
    vector<vector<size_t>> &outGragh=dg.outGraph;
    multimap<float, size_t> heap;
    if(oneSkyband.size()<k){
        for (int i = 0; i < oneSkyband.size(); ++i) {
            float score=dot(D[oneSkyband[i]], w, dim);
            heap.insert({score, oneSkyband[i]});
        }
        while(heap.size()<k){
            heap.insert({0, 1});
        }
    }else {
        for (int i = 0; i < k; ++i) {
            float score = dot(D[oneSkyband[i]], w, dim);
            heap.insert({score, oneSkyband[i]});
        }
    }
    for (int i = k; i < oneSkyband.size(); ++i) {
        float score=dot(D[oneSkyband[i]], w, dim);
        float tScore=heap.begin()->first;
        if(score>tScore){
            heap.erase(heap.begin());
            heap.insert({score, oneSkyband[i]});
        }else if(score>=tScore-1e-7){
            heap.insert({score, oneSkyband[i]});
        }
    }

    for (int i = 0; i < k-1; ++i) {
        auto iter=heap.rbegin();
        for (int j = 0; j < i; ++j) {
            iter++;
        }
        size_t inVtx=iter->second;
        float tScore=heap.begin()->first;

        for(size_t outVtx: outGragh[inVtx]){
            float score=dot(D[outVtx], w, dim);
            if(score>tScore){
                // see if this record is in heap
                if(!findValue(heap, outVtx)){
                    heap.erase(heap.begin());
                    heap.insert({score, outVtx});
                    tScore=heap.begin()->first;
                }
            }else if(score>=tScore-1e-7){
                if(!findValue(heap, outVtx)){
                    if(!findValue(heap, outVtx)) {
                        heap.insert({score, oneSkyband[i]});
                    }
                }
            }
        }
    }
    vector<int> ret;
    ret.reserve(k);
    float last_score;
    auto iter=heap.rbegin();
    for (int i=0; iter!=heap.rend() && i<k; ++iter, ++i) {
        ret.push_back(iter->second);
        last_score=iter->first;
    }
    while(iter!=heap.rend()){
        if(iter->first >= last_score-1e-7){
            ret.push_back(iter->second);
            ++iter;
        }else{
            break;
        }
    }
    return ret;
}

vector<pair<int, double>> foru(vector<float>&w,float** P, int dim, int k, int m,
                               vector<vector<int>> &doGraph, vector<int> &kskybandR, vector<int> &OneSkybandR,
                               Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree);

vector<pair<int, double>> foru_open(vector<float>&w,float** P, int dim, int k, int m,
                                    vector<vector<int>> &doGraph, vector<int> &kskybandR, vector<int> &OneSkybandR,
                                    Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree);

vector<pair<int, double>> foru2(vector<float>&w,float** P, int dim, int k, int m,
                               vector<vector<int>> &doGraph, vector<int> &kskybandR, vector<int> &OneSkybandR,
                               Rtree *lrtree, unordered_map<long int, RtreeNode *> &lramTree);

double non_order_feasible(vector<int> &gt, vector<int> &le, float ** PG, vector<float> &w, vector<double> &x);

vector<pair<int, double>> forut(vector<float>&w,float** P, int dim, int k, int m,
                             unordered_map<uint, uint> &A);

vector<vector<int>> buildDominateGraph(float** P, int dim, vector<int> &kskybandR, int objCnt);


class myHeap{
public:
    multimap<double, uint, greater<>> dist_id;
    unordered_map<uint, double> ids;
    void insert(int id, double dist){
        auto iter=ids.find(id);
        if(iter!=ids.end()){
            auto range=dist_id.equal_range(iter->second);
            for(auto i=range.first; i!=range.second;++i){
                if(i->second==id){
                    dist_id.erase(i);
                    break;
                }
            }
//            for(auto iter2=dist_id.begin();iter2!=dist_id.end(); ++iter2){
//                if(iter2->second==id){
//                    dist_id.erase(iter2);
//                    break;
//                }
//            }
            dist_id.emplace(dist, id);
            ids[id]=dist;
        }else{
            ids.emplace(id, dist);
            dist_id.emplace(dist, id);
        }
    }

    size_t size()const{
        return ids.size();
    }

    void pop(){
        ids.erase(dist_id.begin()->second);
        dist_id.erase(dist_id.begin());
    }

    _Rb_tree_iterator<pair<const double, unsigned int>> top(){
        return dist_id.begin();
    }

    unordered_map<uint, double>::iterator findRecord(int id){
        return ids.find(id);
    }

    unordered_map<uint, double>::iterator recordEnd(){
        return ids.end();
    }

    double topDist(){
        return dist_id.begin()->first;
    }

    vector<pair<int, double>> report()const{
        cout<<"$"<<dist_id.size()<<","<<ids.size()<<endl;
        vector<pair<int, double>> ret;
        ret.reserve(dist_id.size());
        for(auto iter=dist_id.rbegin(); iter!=dist_id.rend(); ++iter){
            ret.emplace_back(iter->second, iter->first);
        }
        return ret;
    }

};

vector<pair<int, double>> foru_dg(vector<float>&w,float** P, int dim, int k, int m, doGraphPlus &dg);

#endif