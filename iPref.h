#ifndef _IPREF_H_
#define _IPREF_H_

#include "header.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "virtualRNode.h"
#include "filemem.h"
#include "global.h"
#include "skyline.h"
#include "qhull_user.h"

extern unordered_map<long int, RtreeNode*> ramTree;


vector<int> computeTopK(const int dim, float* PG[], vector<int> &skyband, vector<float>& weight, int k);

//compute Top-K result set
vector<long int> computeTopK(const int dim, float* PG[], vector<long int> &skyband, vector<float>& weight, int k);

// Is pj dominate pi in traditional sense
bool IsPjdominatePi(const int dimen, float* PG[], long int pi, long int pj);

// test whether pj is incomparable with pi, and pj w > pi w
bool incomparableset(float* PG[], long int pi, long int pj, vector<float>& weight);

// compute the hyperplane of incomparable record pair: pi and pj
vector<float> computePairHP(const int dimen, float* PG[], long int pi, long int pj);

// compute the distance from w to hyperplane HS
float computeDis(const vector<float> &tmpHS, const vector<float> &userpref);
// sort pair
bool sortbysec(const pair<long int, float> &a, const pair<long int, float> &b);

// optimized algorithm
float computeRho(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[],
         vector<pair<long int, float>> &interval, float radius = INFINITY);

// unknown x baseline
float computeRho_unknownX_basic(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float** PG);

// unknown x efficient
float computeRho_unknownX_efficient(const int dimen, const int k, const int X, vector<float>& userpref, Rtree& a_rtree, float* PG[]);

// compute the radius in phase 2, optimized algorithm
float computeradius(const int k, const int dim, long int pi, vector<float>& userpref, vector<long int>& incompSet, float* PG[], float rho=INFINITY);




class unknown_X_node{
private:
    int dim;
    multiset<float> topK_dominate_radius;// size with k
public:
    int tau;// used in unknown efficient
    int needed_to_update;
    const float *data;
    int page_id;
    bool fetched;
    unknown_X_node(int d,  const float *dt, int pid);
    float radius_LB();

    void init_radius(vector<long int> &fetched_options, float **PointSet, vector<float> &w, int k, float rho=INFINITY);

    inline bool update(int need_to_update) const;

    inline bool update() const;

    void update_radius(vector<long int>::iterator begin, vector<long int>::iterator end, float **PointSet, vector<float> &w, float rho=INFINITY);

    void update_radius(const float* other_option, vector<float> &w);

    void update_radius_erase(const float* other_option, vector<float> &w);
};

class unknown_x_baseline{
    multimap<float, unknown_X_node*, greater<float>> heap; //BBS heap, <w\cdot node, node>, if node is an option then node=id+MAXPAGEID
    unordered_set<unknown_X_node*> S;
    vector<long int> incompSet;
    float pt[MAXDIMEN];
    vector<float> ones;
    vector<float> zeros;
    unknown_X_node *zeros_node;
    unknown_X_node *rt;
    vector<unknown_X_node*> to_dl; // todo change into intel pointer
    int k;
    vector<float> userpref;
    float** PG;
    int dimen;
    int next;
public:
    vector<pair<int, float>> interval; // return
    unknown_x_baseline(const int dim, const int K, vector<float>& userPref, Rtree& a_rtree, float** pg){
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
    pair<int, float> get_next(){
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
//                delete(LB.second);
//                    if (interval.size() == X) {
//                        break;
//                    }
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
    ~unknown_x_baseline(){
        for(unknown_X_node *node:to_dl){
            delete (node);
        }
    }


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
    float **PG;
public:
    vector<pair<PG_IDX, RADIUS>> interval; // return

    unknown_x_efficient(const int dim, int K, vector<float> &userPref, Rtree &aRtree, float **pg);

    pair<int, float> get_next();

    ~unknown_x_efficient();
};

void utk_basic(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
               vector<pair<int, double>> &utk_option_ret,
               vector<pair<vector<int>, vector<vector<double>>>> &utk_cones_ret);

vector<int> get_CH_pdtid(const vector<int> &pdt_ids, Qhull &q);


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


int utk_efficient(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                   vector<pair<int, double>> &utk_option_ret,
                   vector<pair<double, region*>> &utk_cones_ret);

vector<vector<double>> points_to_halfspace(vector<vector<double>> &points);

int utk_efficient_anti(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                       vector<pair<int, double>> &utk_option_ret,
                       vector<pair<double, region*>> &utk_cones_ret);
#endif