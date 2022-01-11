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

extern unordered_map<long int, RtreeNode*> ramTree;

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
    int needed_to_update; // need how many options to update this node's radius
    const float *data;  // the values
    int page_id;       // rtree pageid
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
    unordered_set<unknown_X_node*> S; // the nodes in rtree
    vector<long int> incompSet;      // a list store the options that poped from rtree
    float pt[MAXDIMEN]; //a tmp array
    vector<float> ones;
    vector<float> zeros;
    unknown_X_node *zeros_node;
    unknown_X_node *rt;
    vector<unknown_X_node*> to_dl; // todo change into intel pointer // store the heap val that to be deleted
    int k;
    vector<float> userpref;
    float** PG;
    int dimen;
    int next;
public:
    vector<pair<int, float>> interval; // return <option id, radius>
    unknown_x_baseline(const int dim, const int K, vector<float>& userPref, Rtree& a_rtree, float** pg);
    pair<int, float> get_next();
    ~unknown_x_baseline();
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



int utk_efficient(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                  vector<pair<int, double>> &utk_option_ret,
                  vector<pair<double, region*>> &utk_cones_ret);

vector<vector<double>> points_to_halfspace(vector<vector<double>> &points);

int utk_efficient_anti(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                       vector<pair<int, double>> &utk_option_ret,
                       vector<pair<double, region*>> &utk_cones_ret);

int utk_efficient_cs3(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                      vector<pair<int, double>> &utk_option_ret,
                      vector<pair<double, region*>> &utk_cones_ret, double &rho_star);

int utk_efficient3(float **PointSet, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                   vector<pair<int, double>> &utk_option_ret,
                   vector<pair<double, region*>> &utk_cones_ret);

int non_order_utk_efficient(float **PointSet, int pt_cnt, int dim, vector<float> &w, Rtree* rtree, int X, int k,
                            vector<pair<int, double>> &utk_option_ret,
                            vector<pair<double, region*>> &utk_cones_ret);
#endif