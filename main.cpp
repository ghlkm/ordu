
#include <chrono>
#include "rtree_lib/header.h"
#include "rtree_lib/hypercube.h"
#include "rtree_lib/rentry.h"
#include "rtree_lib/rnode.h"
#include "rtree_lib/rtree.h"
#include "rtree_lib/filemem.h"
#include "rtree_lib/tgs.h"
#include "rtree_lib/global.h"
#include "skyline.h"
#include "rtree_lib/param.h"
#include "iPref.h"
#include "qp_solver.h"
#include "math_lib.h"
#include "case_study.h"

// gobal variables
vector<vector<float>> HalfSpaces; // halfspace 
unordered_map<long int, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory
qp_solver qp;  //qp solver use for calculate dominate radius
vector<vector<c_float>> g_r_domain_vec; // use for r-dominate testing

int objCnt = 0; // # of data objects
double totalSpaceCost = 0.0; // space cost (MB)
clock_t at, ad;


int main(const int argc, const char** argv)
{
    cout.precision(6);
    cout << "ORD/ORU" << endl;
    clock_t at, ad;

    // parameter parser
    cout << "Parse Parameters" << endl;
    if (argc == 1)
    {
        helpmsg(argv[0]);
        return -1;
    }
//    product file format:
//      <id1  w_1,1  w_1,2  w_1,3  ...  w_1,d  w_1,d+1  w_1,d+2  w_1,d+3  ...  w_1,d+d>
//      <id2  w_2,1  w_2,2  w_2,3  ...  w_2,d  w_2,d+1  w_2,d+2  w_2,d+3  ...  w_2,d+d>
//      ...
//      <idn  w_n,1  w_n,2  w_n,3  ...  w_n,d  w_n,d+1  w_n,d+2  w_n,d+3  ...  w_n,d+d>

//    user file format:
//      <w_11 w_12 w_13 ... w_1d>
//      <w_21 w_22 w_23 ... w_2d>
//      ...
//      <w_n1 w_n2 w_n3 ... w_nd>

    const int k = atoi(Param::read(argc, argv, "-k", ""));
    int dim = atoi(Param::read(argc, argv, "-d", ""));
    const char* datafile = Param::read(argc, argv, "-f", "");  // option file name
    // index file for rtree
    const char* indexfile = Param::read(argc, argv, "-i", "");
    const int m = atoi(Param::read(argc, argv, "-m", ""));
    const char* methodName = Param::read(argc, argv, "-mt", "");
    int w_num = atoi(Param::read(argc, argv, "-w", "")); // the number of tested user weights
    int n=atoi(Param::read(argc, argv, "-n", "")); //  specific load how many options

    const char *w_file = Param::read(argc, argv, "-W", "");// user preference file

    vector<vector<float>> ws(w_num*3, vector<float>(dim));
    vector<int> ks(w_num*3);
    fstream fpdata;
    fpdata.open(w_file, ios::in);
    for (int i = 0; i < ws.size(); ++i) {
        ks[i]=k;
        for (int j = 0; j < ws[i].size(); ++j) {
            fpdata >> ws[i][j];
        }
    }
    fpdata.close(); // in case of memory leaking


    int resultSize = 0;
    // data loading
    cout << "Load data points from file" << endl;
    float** PointSet = new float*[MAXPTS + 1];
    RtreeNodeEntry** p = new RtreeNodeEntry*[MAXPTS];
    fpdata.open(datafile, ios::in);

    while (n--)
    {
        int id;
        float* cl = new float[dim];
        float* cu = new float[dim];
        fpdata >> id;
        if (fpdata.eof())
            break;

        PointSet[objCnt + 1] = new float[2 * dim];

        for (int d = 0; d < dim; d++)
        {
            fpdata >> cl[d];
            PointSet[objCnt + 1][d] = cl[d];
        }

        for (int d = 0; d < dim; d++)
        {
            fpdata >> cu[d];
            PointSet[objCnt + 1][d + dim] = cu[d];
        }

        Hypercube hc(dim, cl, cu);
        p[objCnt++] = new RtreeNodeEntry(id, hc);

        //log information
        if (objCnt % 1000 == 0)
            cout << ".";
        if (objCnt % 10000 == 0)
            cout << objCnt << " objects loaded" << endl;
    }

    double rawSize = dataSize(objCnt, dim);
    cout << "Total number of objects: " << objCnt << endl;
    totalSpaceCost += rawSize;

    // build rtree
    cout << "Bulkloading R-tree..." << endl;
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    FileMemory mem(PAGESIZE, indexfile, RtreeNodeEntry::fromMem, true);
    Rtree* rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild*0.3, (int)maxChild*0.3, p, objCnt, false);
    cout << "[Rtree build done]" << endl;

    // in-memory rtree
    cout << "cache R-tree into memory" << endl;
    rtreeRAM(*rtree, ramTree);
    totalSpaceCost += ramTree.size()*4096.00 / MB;
    cout << "Total R-tree size:" << totalSpaceCost << "MB" << endl;

    // aggregate rtree
    aggregateRecords(*rtree);
    cout << "[Aggregate Rtree done]" << endl;
    // init qp solver
    qp=qp_solver(dim);
    g_r_domain_vec=gen_r_domain_vec(dim);
    if (strcmp(methodName, "BB") == 0)
    {
        // ORD knownX baseline algorithm
        // (1) compute k-skyband set SK
        // (2) for each pi, suppose T_pi is the incomparable set (a) p_j is incomparable with p_i in SK, (2) p_j w > p_i w
        // (3) for each p in T_pi, compute the distance from point w to hyperplane H_{i,j} equation: d = \frac{a1 w1 +  a2 w2 + ... ad_1 wd-1+D}{\sqrt{a1^2+a2^2+ ... + ad-1^2}}
        // (4) compute the inflection point of pi: the k-th laregst value in step (3), denotes as inf_pi
        // (5) pi's rksykband interval: inf_pi to infinity
        // (6) the radius rho is the T-th minimum value in all pi in SK
        int k=ks[0];
        clock_t bat = clock();
        vector<long int> skyband;
        auto kbegin = chrono::steady_clock::now();
        kskyband(dim, *rtree, skyband, PointSet, k); // step (1)
        cout<< skyband.size() << endl;
        auto kend = chrono::steady_clock::now();
        at=clock();
        auto begin = chrono::steady_clock::now();

        clock_t bad = clock();
        cout << "Total time cost: " << fixed << (bad - bat) * 1.0 / (CLOCKS_PER_SEC) << " SEC " << endl;
        for (int wi = 0; wi < w_num; wi++)
        {
            // weight vector for testing, we should remove the redundant one
            auto wbegin = chrono::steady_clock::now();

            vector<float> w(ws[wi].begin(), ws[wi].end());

            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout<< w.back() <<endl;
            cout << "Testing k: " << k <<endl;
            //k-skyband
            vector<long int> topKRet = computeTopK(dim, PointSet, skyband, w, k);
            assert(topKRet.size() == k);
            vector<pair<long int, float>> interval;
            for (int i = 0; i < topKRet.size(); i++) {
                interval.emplace_back(topKRet[i], 0);
            }

            for (int ski = 0; ski < skyband.size(); ski++)
            {
                if (find(topKRet.begin(), topKRet.end(), skyband[ski]) != topKRet.end())
                {
                    continue;
                }
                multiset<float> radiusSKI;
                vector<long int> incompset;
                vector<long int> dominatorSet;
                int dominated_cnt=0;
                for (int pj = 0; pj < skyband.size(); pj++)
                {
                    if (ski == pj)
                        continue;
                    if(v1_dominate_v2(PointSet[skyband[pj]], PointSet[skyband[ski]], dim))
                    {
                        radiusSKI.insert(INFINITY);
                        dominated_cnt++;
                        if(dominated_cnt>=k)
                            break;
                    }
                    else if(!v1_dominate_v2(PointSet[skyband[ski]], PointSet[skyband[pj]], dim)
                    && dot(w, PointSet[skyband[ski]])< dot(w, PointSet[skyband[pj]]))// step (2)
                    {
                        incompset.push_back(skyband[pj]);
                    }
                }

                // here we need a function to compute the inflection radius of option pi
                if(dominated_cnt>=k){
                    interval.emplace_back(skyband[ski], INFINITY);
                    continue;
                }
                while(radiusSKI.size()>k){
                    radiusSKI.erase(radiusSKI.begin());
                }
                for (int inpi = 0; inpi < incompset.size(); inpi++)
                {
                    vector<float> tmpHS = computePairHP(dim, PointSet, skyband[ski], incompset[inpi]);
                    //compute the distance from w to hyperplane.
                    float tmpDis = computeDis(tmpHS, w);
                    radiusSKI.insert(tmpDis);
                    if(radiusSKI.size()>k){
                        radiusSKI.erase(radiusSKI.begin());
                        if(*radiusSKI.begin()==INFINITY){
                            break;
                        }
                    }
                }
                if(*radiusSKI.begin()==INFINITY){
                    interval.emplace_back(skyband[ski], INFINITY);
                }else{
                    assert(radiusSKI.size() >= k);
                    interval.emplace_back(skyband[ski], *radiusSKI.begin());
                }
            }
            sort(interval.begin(), interval.end(), sortbysec);
            for (auto i = 0; i < k; ++i) {
                interval[i].first = topKRet[i]; // although top-k is the same, we guarantee the correct top-k order
            }
            cout << "The inflection radius is: " << interval[m].second << endl;
            auto now = chrono::steady_clock::now();
            chrono::duration<double> elapsed_seconds= (now-wbegin)+(kend-kbegin);
            cout << "Total time cost: " << elapsed_seconds.count()<< " SEC " << endl;
        }
        ad = clock();
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= (now-begin)/w_num+(kend-kbegin);
        cout << "Total time cost: " << elapsed_seconds.count()<< " SEC " << endl;
    }
    if (strcmp(methodName, "OA") == 0)
    {
        // ORD knownX optimized algorithm
        // (1) revised BBS for w, compute top-k result, store in T directly
        // (2) for each pi (uncertain status), compute its inflection distance inf_pi, and push it into a max-heap Q (inf_pi as key)
        // (3) keep fetching until max-heap with size m-k+1, pop the top node out, and initalize p0 = inf_ptop
        // (4) test the rest candidates by r-dominanace, instead of BBS, thus fetching r-skyband options
        // (5) shrinks p0 quicly, simulating a scanning from infinity to 0
        // (6) fetching terminates none of options can be the r-skyband
        // (7) append Q to T, this is the final result.
        vector<pair<long int, float>> interval;
        at = clock();
        auto begin = chrono::steady_clock::now();
        for (int wi = 0; wi < w_num; wi++)
        {
            clock_t wat = clock();
            auto wbegin = chrono::steady_clock::now();

            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());

            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;
            float rho = computeRho(dim, k, m, w, *rtree, PointSet, interval);
            interval.clear();
            cout << "The inflection radius is: " << rho << endl;
            ad = clock();
            cout << "Total time cost: " << fixed << (ad - wat) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << "\n";
            auto now = chrono::steady_clock::now();
            chrono::duration<double> elapsed_seconds= now-wbegin;
            cout << elapsed_seconds.count() <<endl;
        }
        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << "\n";
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
        cout << elapsed_seconds.count()/w_num <<endl;
    }
    if (strcmp(methodName, "UB") == 0) // ORD unknown m basic
    {
        // (1) if get_next_time<=k, from BBS fetch topK
        // (2) for each element in BBS, calculate their inflection radius and store them into heap S (so we can know whose inflection radius is smallest)
        // (3) while the the top of S is not an option:
        //       (a) pop and push node A out and into BBS heap
        //       (b) if node A is an option:
        //             not update S
        //             marked A as fetched
        //           else
        //             update S to remove node
        //           update all nodes in S such that not fetched by BBS yet (use a flag FETCH to record)
        at = clock();
        for (int wi = 0; wi < w_num; wi++)
        {
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;

            float rho = computeRho_unknownX_basic(dim, k, m, w, *rtree, PointSet);
            cout << "The inflection radius is: " << rho << endl;
        }
        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
    }
    if (strcmp(methodName, "UB_GN") == 0) // ORD unknown m baseline get_next version
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(m);
        vector<float> avg_radius(m);
        for (int wi = 0; wi < w_num; wi++)
        {
            auto w_begin = chrono::steady_clock::now();
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;

            unknown_x_baseline obj(dim, k, w, *rtree, PointSet);
            for (int i = 0; i < m; ++i) {
                obj.get_next();
                auto now = chrono::steady_clock::now();
                chrono::duration<double> elapsed_seconds= now-w_begin;
                avg_time[i]+=elapsed_seconds.count();
            }
            for (int i = 0; i < m; ++i) {
                avg_radius[i]+=obj.interval[i].second;
            }
            float rho = obj.interval.back().second;
            cout << "The inflection radius is: " << rho << endl;
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        for (double time:avg_time) {
            cout<<time/w_num<<endl;
        }
        for (float radius:avg_radius) {
            cout<<radius/w_num<<endl;
        }
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
    }
    // incremental version, without known m
    // We do not have exact m, we need tell the user the radius rho and its corresponding T
    // It is similar to optimized algorithm, however, it computes incrementally, from rho = 0 to infinity, the size T is from k to k-skyband.
    if (strcmp(methodName, "UA") == 0) // unknown m efficient
    {
        at = clock();
        for (int wi = 0; wi < w_num; wi++)
        {
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;

            float rho = computeRho_unknownX_efficient(dim, k, m, w, *rtree, PointSet);
            cout << "The inflection radius is: " << rho << endl;
        }
        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
    }
    if (strcmp(methodName, "UA_GN") == 0) // unknown m efficient get_next version
    {
        vector<double> grank(m);

        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(m);
        vector<float> avg_radius(m);
        for (int wi = 0; wi < w_num; wi++)
        {
            auto w_begin = chrono::steady_clock::now();
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());

            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;

            unknown_x_efficient obj(dim, k, w, *rtree, PointSet);
            for (int i = 0; i < m; ++i) {
                obj.get_next();
                auto now = chrono::steady_clock::now();
                chrono::duration<double> elapsed_seconds= now-w_begin;
                avg_time[i]+=elapsed_seconds.count();
            }
            for (int i = 0; i < m; ++i) {
                avg_radius[i]+=obj.interval[i].second;
            }
            float rho = obj.interval.back().second;
            cout << "The inflection radius is: " << rho << endl;
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        for (float radius:avg_radius) {
            cout<<radius/w_num<<endl;
        }
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
    }
    if (strcmp(methodName, "UTK_BB") == 0) // ORU baseline
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(m);
        vector<float> avg_radius(m);
        for (int wi = 0; wi < w_num; wi++)
        {
            auto w_begin = chrono::steady_clock::now();
            int k=ks[wi];
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;
            vector<pair<int, double>> utk_option_ret;
            vector<pair<vector<int>, vector<vector<double>>>> utk_cones_ret;
            utk_basic(PointSet, dim, w, rtree, m, k, utk_option_ret, utk_cones_ret);

            double rho = utk_option_ret.back().second;
            cout << "The inflection radius is: " << rho << endl;
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
    }
    double avg_rt_cnt=0;
    double avg_rr_cnt=0;
    if (strcmp(methodName, "UTK_OA") == 0) // ORU efficient
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(m);
        vector<float> avg_radius(m, 0.0);
        for (int wi = 0; wi < w_num; wi++)
        {
            auto w_begin = chrono::steady_clock::now();
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;
            vector<pair<int, double>> utk_option_ret;
            vector<pair<double, region*>> utk_cones_ret;
            // the code commented below just for anti data
//            int generated_r_cnt= utk_efficient_anti(PointSet, dim, w, rtree, m, k, utk_option_ret,utk_cones_ret);
            int generated_r_cnt= utk_efficient(PointSet, dim, w, rtree, m, k, utk_option_ret, utk_cones_ret);
            avg_rt_cnt+=generated_r_cnt;
            avg_rr_cnt+=utk_cones_ret.size();
            cout<<"ret size: "<<utk_option_ret.size()<<"\n";
            for (int i = 0; i < avg_radius.size(); ++i) {
                if(i<utk_option_ret.size()){
                    avg_radius[i]+=utk_option_ret[i].second;
                }else{
                    avg_radius[i]+=1000000;
                }
            }
            for (float radius:avg_radius) {
                cout<<radius/(wi+1)<<endl;
            }
            double rho = utk_option_ret.back().second;
            for (pair<double, region*> &tmp: utk_cones_ret) {
                delete(tmp.second);
            }
            cout << "The inflection radius is: " << rho << endl;
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        cout << "Total generated regions: " << avg_rt_cnt/w_num<<endl;
        cout << "Total return regions: " << avg_rr_cnt/w_num<<endl;
        for (float radius:avg_radius) {
            cout<<radius/w_num<<endl;
        }
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
        cout<<elapsed_seconds.count();

        cout<<"output time stat:\n";
    }
    if (strcmp(methodName, "UTK_OA3") == 0) // ORU efficient
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(m);
        vector<float> avg_radius(m, 0.0);
        for (int wi = 0; wi < w_num; wi++)
        {
            auto w_begin = chrono::steady_clock::now();
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;
            vector<pair<int, double>> utk_option_ret;
            vector<pair<double, region*>> utk_cones_ret;
            // the code commented below just for anti data
//            int generated_r_cnt= utk_efficient_anti(PointSet, dim, w, rtree, m, k, utk_option_ret,utk_cones_ret);
            int generated_r_cnt= utk_efficient3(PointSet, dim, w, rtree, m, k, utk_option_ret, utk_cones_ret);
            avg_rt_cnt+=generated_r_cnt;
            avg_rr_cnt+=utk_cones_ret.size();
            cout<<"ret size: "<<utk_option_ret.size()<<"\n";
            for (int i = 0; i < avg_radius.size(); ++i) {
                if(i<utk_option_ret.size()){
                    avg_radius[i]+=utk_option_ret[i].second;
                }else{
                    avg_radius[i]+=1000000;
                }
            }
            for (float radius:avg_radius) {
                cout<<radius/(wi+1)<<endl;
            }
            double rho = utk_option_ret.back().second;
            for (pair<double, region*> &tmp: utk_cones_ret) {
                delete(tmp.second);
            }
            cout << "The inflection radius is: " << rho << endl;
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        cout << "Total generated regions: " << avg_rt_cnt/w_num<<endl;
        cout << "Total return regions: " << avg_rr_cnt/w_num<<endl;
        for (float radius:avg_radius) {
            cout<<radius/w_num<<endl;
        }
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
        cout<<elapsed_seconds.count();
    }
    if (strcmp(methodName, "CS") == 0) // case study
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        for (int wi = 0; wi < w_num; wi++)
        {
            auto w_begin = chrono::steady_clock::now();
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;
            vector<int> skyline_topm=OSS_skyline(objCnt, m, rtree, PointSet, dim);
            set<int> sky_topm_s(skyline_topm.begin(), skyline_topm.end());


            // fetch top-m directly
            vector<int> all_id(objCnt);
            iota(all_id.begin(), all_id.end(), 1);
//          computeTopK(const int dim, float* PG[], vector<int> &skyband, vector<float>& weight, int k)
            vector<int> direct_topm=computeTopK(dim, PointSet, all_id, w, m);
            set<int> d_topm_s(direct_topm.begin(), direct_topm.end());

            // fetch ORD top-m
            vector<pair<long int, float>> interval;
            float rho = computeRho(dim, k, m, w, *rtree, PointSet, interval);
            cout<<"rho:"<<rho<<endl;
            vector<int> ord_topm;
            for(pair<long int, float>&pi: interval){
                ord_topm.push_back(pi.first);
            }
            set<int> ord_topm_s(ord_topm.begin(), ord_topm.end());

            // fetch ORU top-m
            vector<pair<int, double>> utk_option_ret;
            vector<pair<double, region*>> utk_cones_ret;
//            int generated_r_cnt= utk_efficient_anti(PointSet, dim, w, rtree, m, k, utk_option_ret,utk_cones_ret);
            int generated_r_cnt= utk_efficient(PointSet, dim, w, rtree, m, k, utk_option_ret, utk_cones_ret);
            vector<int> oru_topm;
            for(pair<int, double>&pi: utk_option_ret){
                oru_topm.push_back(pi.first);
            }
            set<int> oru_topm_s(oru_topm.begin(), oru_topm.end());
            double rho2 = utk_option_ret.back().second;
//            cout<<"rho2:"<<rho2<<endl;

            cout << "1-skyband top m size: "<< sky_topm_s.size()<<endl;
            cout << "direct top m size: " << d_topm_s.size()<<endl;
            cout << "ORD top m size: " << ord_topm_s.size()<<endl;
            cout<< "ORU top m size: "<<oru_topm_s.size()<<endl;

//            assert(sky_topm_s.size()==m);
//            assert(d_topm_s.size()==m);
//            assert(ord_topm_s.size()==m);
//            assert(oru_topm_s.size()==m);

            // M1 (skyline, ord)
            // M2 (skyline, oru)
            // M3 (d-top-m, ord)
            // M4 (d-top-m, oru)
            int ma[]={10, 30, 50, 70, 90};
            vector<int> mv(ma, ma+5);
            // M1
            int mins=min(skyline_topm.size(), ord_topm.size());
            cout<<"M1\n";
            for (int i1:mv) {
                cout << jaccard(skyline_topm.begin(), skyline_topm.begin() + min(i1, mins), ord_topm.begin(),
                                ord_topm.begin() + min(i1, mins)) << "\n";
            }
            cout<<endl;
            // M2
            cout<<"M2\n";
            for (int i2:mv) {
                cout << jaccard(skyline_topm.begin(), skyline_topm.begin() + i2, oru_topm.begin(), oru_topm.begin() + i2) << "\n";
            }
            cout<<endl;
            // M3
            cout<<"M3\n";
            for (int i3:mv) {
                cout << jaccard(direct_topm.begin(), direct_topm.begin() + i3, ord_topm.begin(), ord_topm.begin() + i3) << "\n";
            }
            cout<<endl;
            // M4
            cout<<"M4\n";
            for (int i4:mv) {
                cout << jaccard(direct_topm.begin(), direct_topm.begin() + i4, oru_topm.begin(), oru_topm.begin() + i4) << "\n";
            }
            cout<<endl;
            cout<<"topm=[\n"<<direct_topm<<"\n]\n";
            cout<<"ord=[\n"<<ord_topm<<"\n]\n";
            cout<<"oru=[\n"<<oru_topm<<"\n]\n";
            cout<<"skyline=[\n"<<skyline_topm<<"\n]\n";
            for (pair<double, region*> &tmp: utk_cones_ret) {
                delete(tmp.second);
            }
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        cout << "Total generated regions: " << avg_rt_cnt/w_num<<endl;
        cout << "Total return regions: " << avg_rr_cnt/w_num<<endl;
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
    }

    if (strcmp(methodName, "CS2") == 0) // case study 2
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> jas; // contains jaccard results
        vector<double> pres; // contains precision results
        vector<double> recs; // contains recall results
        for (int wi = 0; wi < w_num; wi++)
        {
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;


            // fetch ORD top-m
            vector<pair<long int, float>> interval;
            float rho = computeRho(dim, k, m, w, *rtree, PointSet, interval);
            cout<<"rho:"<<rho<<endl;
            vector<int> ord_topm;
            for(pair<long int, float>&pi: interval){
                ord_topm.push_back(pi.first);
            }
            set<int> ord_topm_s(ord_topm.begin(), ord_topm.end());

            // fetch ORU top-m
            vector<pair<int, double>> utk_option_ret;
            vector<pair<double, region*>> utk_cones_ret;
//            int generated_r_cnt= utk_efficient_anti(PointSet, dim, w, rtree, m, k, utk_option_ret,utk_cones_ret);
            int generated_r_cnt= utk_efficient(PointSet, dim, w, rtree, m, k, utk_option_ret, utk_cones_ret);
            vector<int> oru_topm;
            for(pair<int, double>&pi: utk_option_ret){
                oru_topm.push_back(pi.first);
            }
            set<int> oru_topm_s(oru_topm.begin(), oru_topm.end());
            double rho2 = utk_option_ret.back().second;

            cout << "ORD top m size: " << ord_topm_s.size()<<endl;
            cout<< "ORU top m size: "<<oru_topm_s.size()<<endl;

            double ja=jaccard(oru_topm.begin(), oru_topm.end(), ord_topm.begin(),
                              ord_topm.end());
            jas.push_back(ja);
            cout<<"jaccard:";
            cout << ja << "\n";
            cout<<endl;

            double pre=precision(oru_topm.begin(), oru_topm.end(), ord_topm.begin(),
                                 ord_topm.end());
            pres.push_back(pre);
            cout<<"precision:";
            cout << pre << "\n";
            cout<<endl;

            double re=recall(oru_topm.begin(), oru_topm.end(), ord_topm.begin(),
                             ord_topm.end());
            recs.push_back(re);
            cout<<"recall:";
            cout << re << "\n";
            cout<<endl;

            cout<<endl;
            cout<<"ord=[\n"<<ord_topm<<"\n]\n";
            cout<<"oru=[\n"<<oru_topm<<"\n]\n";
            for(int i:oru_topm){
                for (int j = 0; j < dim; ++j) {
                    cout<<PointSet[i][j]+SIDELEN<<" ";
                }
                cout<<endl;
            }
            for (pair<double, region*> &tmp: utk_cones_ret) {
                delete(tmp.second);
            }

        }
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
        cout<<elapsed_seconds.count()<<endl;
//        ofstream myfile;
//        myfile.open ("cs2.txt", ios::out | ios::app | ios::binary);
//        for (int l = 0; l < argc; ++l) {
//            myfile<< argv[l]<<" ";
//        }
//        myfile<<endl;
//        myfile<<"avg jaccard: ";
//        myfile<< sum(jas.begin(), jas.end())/jas.size() <<"\n";
//        myfile<<"avg precision: ";
//        myfile<< sum(pres.begin(), pres.end())/pres.size() <<"\n";
//        myfile<<"avg recall: ";
//        myfile<< sum(recs.begin(), recs.end())/recs.size() <<"\n";
//
//        for(auto ja:jas){
//            myfile<<ja<<" ";
//        }
//        myfile<<"\n";
//
//        for(auto pre:pres){
//            myfile<<pre<<" ";
//        }
//        myfile<<"\n";
//
//        for(auto rec:recs){
//            myfile<<rec<<" ";
//        }
//        myfile<<"\n";
//
//        myfile<<endl;
//        myfile.close();
    }

    if (strcmp(methodName, "CS3") == 0) // case study3, see the change of rho_star with change of d
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> rss;
        for (int wi = 0; wi < w_num; wi++)
        {
            auto w_begin = chrono::steady_clock::now();
            int k=ks[wi];
            // weight vector for testing, we should remove the redundant one
            vector<float> w(ws[wi].begin(), ws[wi].end());
            cout << "Testing w: ";
            for (int di = 0; di < dim-1; di++)
            {
                cout << w[di] << ", ";
            }
            cout <<w.back()<< endl;
//            void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k)
            //get 1-skyband, fetch top m from them
            vector<pair<int, double>> utk_option_ret;
            vector<pair<double, region*>> utk_cones_ret;
            double rho_star;
            int generated_r_cnt= utk_efficient_cs3(PointSet, dim, w, rtree, m, k, utk_option_ret, utk_cones_ret, rho_star);
            rss.push_back(rho_star);
            for (pair<double, region*> &tmp: utk_cones_ret) {
                delete(tmp.second);
            }
        }
        cout<<"output rho_star begin"<<endl;
        for(double rho:rss){
            cout<<rho<<"\n";
        }
        cout<<"output rho_star end"<<endl;

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        cout << "Total generated regions: " << avg_rt_cnt/w_num<<endl;
        cout << "Total return regions: " << avg_rr_cnt/w_num<<endl;
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
    }
//    ofstream myfile;
//    myfile.open ("result_oru.txt", ios::out | ios::app | ios::binary);
//    myfile <<fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num)<<": ";
//    myfile << "(" << avg_rt_cnt/w_num << "," << avg_rr_cnt/w_num << ") ";
//    for (int l = 0; l < argc; ++l) {
//        myfile<< argv[l]<<" ";
//    }
//    myfile<<endl;
//    myfile.close();
    return 0;
}
