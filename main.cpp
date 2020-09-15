/*----------------------- iPref problem ---------------
----------------------DBGroup@SUSTech, 2020------------
-------------------- tangloner@gmail.com-------------*/

#include <chrono>
#include "header.h"
#include "hypercube.h"
#include "rentry.h"
#include "rnode.h"
#include "rtree.h"
#include "filemem.h"
#include "tgs.h"
#include "global.h"
#include "skyline.h"
#include "param.h"
#include "iPref.h"
#include "qp_solver.h"
#include "utk_math_lib.h"
/*the headers for osqp solver*/
//#include "qp_solver.h"
//#include "osqp.h"

// gobal variables
vector<vector<float>> HalfSpaces; // halfspace 
unordered_map<long int, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory
qp_solver qp;  //qp solver use for calculate dominate radius
vector<vector<c_float>> g_r_domain_vec; // use for r-dominate testing

int objCnt = 0; // # of data objects
double totalIO = 0; // # of IO access
double totalSpaceCost = 0.0; // space cost (MB)
double treeSpacecost = 0.0;
double Spacecost = 0.0;
clock_t at, ad;
//#include "lp_lib/lp_lib.h"
//int lp_solve_test(){
//    lprec *lp;
//
//    /* Create a new LP model */
//    lp = make_lp(0, 2);
//    if(lp == NULL) {
//        fprintf(stderr, "Unable to create new LP model\n");
//        return(1);
//    }
////    set_verbose(lp, IMPORTANT);
////    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
//    set_add_rowmode(lp, TRUE);
//    REAL r1[]={0, 1.0, 2.0};
//    add_constraint(lp, r1, GE, 0.0);
//    REAL r2[]={0, 1.0, 1.0};
//    add_constraint(lp, r2, GE, 1.0);
//    REAL r3[]={0, 0.0, -1.0};
//    add_constraint(lp, r3, GE, -1.0);
//    REAL r4[]={0, -1.0, 0.0};
//    add_constraint(lp, r4, GE, -1.0);
//    set_add_rowmode(lp, FALSE);
//    set_timeout(lp, 1);
//
//    cout << solve(lp) << endl;
//    int ccnt=get_Nrows(lp);
//    REAL *constr=new REAL[ccnt]();
//    get_constraints(lp, constr);
//    for (int i = 0; i <ccnt ; ++i) {
//        cout<<constr[i]<<endl;
//    }
//    delete_lp(lp);
//    delete [] (constr);
//    return(0);
//}
int main(const int argc, const char** argv)
{
	cout.precision(6);
	cout << "iPref Problem (Size-constrained R-kSkyband/UTK )" << endl;
	clock_t at, ad;

	// parameter parser
	cout << "Parse Parameters" << endl;
	if (argc == 1)
	{
		helpmsg(argv[0]);
		return -1;
	}
    const int k = atoi(Param::read(argc, argv, "-k", ""));
	int dim = atoi(Param::read(argc, argv, "-d", ""));
	const char* datafile = Param::read(argc, argv, "-f", "");
	const char* indexfile = Param::read(argc, argv, "-i", "");
	const int X = atoi(Param::read(argc, argv, "-X", ""));
	const char* methodName = Param::read(argc, argv, "-m", "");
    int w_num = atoi(Param::read(argc, argv, "-w", "")); // the number of tested user weights
    int n=atoi(Param::read(argc, argv, "-n", ""));
//    write file format:
//      <k1, w_11, w_12, w_13, ..., w_1d>
//      <k2, w_21, w_22, w_23, ..., w_2d>
//      ...
//      <kn, w_n1, w_n2, w_n3, ..., w_nd>
    const char *w_file = Param::read(argc, argv, "-W", "");

    vector<vector<float>> ws(w_num, vector<float>(dim));
    vector<int> ks(w_num);
	fstream fpdata;
    fpdata.open(w_file, ios::in);
    for (int i = 0; i < ws.size(); ++i) {
//        fpdata >> ks[i];
        ks[i]=k;
        for (int j = 0; j < ws[i].size(); ++j) {
            fpdata >> ws[i][j];
        }
    }
    fpdata.close();


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
    qp=qp_solver(dim);
    g_r_domain_vec=gen_r_domain_vec(dim);
    if (strcmp(methodName, "BB") == 0)
	{
		// baseline algorithm
		// (1) compute k-skyband set SK
		// (2) for each pi, suppose T_pi is the incomparable set (a) p_j is incomparable with p_i in SK, (2) p_j w > p_i w
		// (3) for each p in T_pi, compute the distance from point w to hyperplane H_{i,j} equation: d = \frac{a1 w1 +  a2 w2 + ... ad_1 wd-1+D}{\sqrt{a1^2+a2^2+ ... + ad-1^2}}
		// (4) compute the inflection point of pi: the k-th laregst value in step (3), denotes as inf_pi
		// (5) pi's rksykband interval: inf_pi to infinity
		// (6) the radius rho is the T-th minimum value in all pi in SK

		at = clock();
		for (int wi = 0; wi < w_num; wi++)
		{
			vector<long int> skyband;
			int k=ks[wi];
			kskyband(dim, *rtree, skyband, PointSet, k); // step (1)
            cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
			// weight vector for testing, we should remove the redundant one
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
					if (IsPjdominatePi(dim, PointSet, skyband[ski], skyband[pj]))
					{
						radiusSKI.insert(INFINITY);
                        dominated_cnt++;
						if(dominated_cnt>=k)
						    break;
					}
					else if (incomparableset(PointSet, skyband[ski], skyband[pj], w)) // step (2)
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
//			assert(interval.size() >= X );
			sort(interval.begin(), interval.end(), sortbysec);
			for (auto i = 0; i < k; ++i) {
				interval[i].first = topKRet[i];
			}

//			int cnt = 0;
//			for (auto i = interval.begin(); i != interval.end(); ++i)
//			{
//				cout << i->second << ", " << i->first << endl;
//				++cnt;
//				if (cnt > X)
//					break;
//			}
			cout << "The inflection radius is: " << interval[X].second << endl;
		}
		ad = clock();
		cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
	}
	if (strcmp(methodName, "OA") == 0)
	{
		// optimized algorithm
		// (1) revised BBS for w, compute top-k result, store in T directly
		// (2) for each pi (uncertain status), compute its inflection distance inf_pi, and push it into a max-heap Q (inf_pi as key)
		// (3) keep fetching until max-heap with size X-k+1, pop the top node out, and initalize p0 = inf_ptop
		// (4) test the rest candidates by r-dominanace, instead of BBS, thus fetching r-skyband options
		// (5) shrinks p0 quicly, simulating a scanning from infinity to 0
		// (6) fetching terminates none of options can be the r-skyband
		// (7) append Q to T, this is the final result.

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
            vector<pair<long int, float>> interval;
			float rho = computeRho(dim, k, X, w, *rtree, PointSet, interval);
			cout << "The inflection radius is: " << rho << endl;
		}
		ad = clock();
		cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
	}


	// inclremental version, without known X
	// We do not have exact X, we need tell the user the radius rho and its corresponding T
	// It is similar to optimized algorithm, however, it computes incrementally, from rho = 0 to infinity, the size T is from k to k-skyband.
    if (strcmp(methodName, "UB") == 0) // unknown X basic
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

            float rho = computeRho_unknownX_basic(dim, k, X, w, *rtree, PointSet);
            cout << "The inflection radius is: " << rho << endl;
        }
        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
    }

    // inclremental version, without known X
    // We do not have exact X, we need tell the user the radius rho and its corresponding T
    // It is similar to optimized algorithm, however, it computes incrementally, from rho = 0 to infinity, the size T is from k to k-skyband.
    if (strcmp(methodName, "UA") == 0) // unknown X efficient
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

            float rho = computeRho_unknownX_efficient(dim, k, X, w, *rtree, PointSet);
            cout << "The inflection radius is: " << rho << endl;
        }
        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
    }

    if (strcmp(methodName, "UA_GN") == 0) // unknown X efficient get_next version
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(X);
        vector<float> avg_radius(X);
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
            for (int i = 0; i < X; ++i) {
                obj.get_next();
                auto now = chrono::steady_clock::now();
                chrono::duration<double> elapsed_seconds= now-w_begin;
                avg_time[i]+=elapsed_seconds.count();
//                cout<<elapsed_seconds.count()<<"\n";
            }
            for (int i = 0; i < X; ++i) {
                avg_radius[i]+=obj.interval[i].second;
//                cout<<obj.interval[i].second<<"\n";
            }
//            pair<int, float> next=obj.get_next();
//            while(!(next.second==INFINITY)){
//                cout<<next.second<<"\n";
//                next=obj.get_next();
//            }
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
    if (strcmp(methodName, "UTK_BB") == 0) // utk baseline
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(X);
        vector<float> avg_radius(X);
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
            vector<pair<vector<int>, vector<vector<double>>>> utk_cones_ret;
            utk_basic(PointSet, dim, w, rtree, X, k, utk_option_ret, utk_cones_ret);

            double rho = utk_option_ret.back().second;
            cout << "The inflection radius is: " << rho << endl;
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
    }
    if (strcmp(methodName, "UTK_OA") == 0) // utk efficient
    {
        at = clock();
        auto begin = chrono::steady_clock::now();
        vector<double> avg_time(X);
        vector<float> avg_radius(X);
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
            utk_efficient(PointSet, dim, w, rtree, X, k, utk_option_ret,utk_cones_ret);
            double rho = utk_option_ret.back().second;
            cout << "The inflection radius is: " << rho << endl;
        }

        ad = clock();
        cout << "Total time cost: " << fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num) << " SEC " << endl;
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds= now-begin;
    }
    ofstream myfile;
    myfile.open ("result2.txt", ios::out | ios::app | ios::binary);
    myfile <<fixed << (ad - at) * 1.0 / (CLOCKS_PER_SEC*w_num)<<": ";
    for (int l = 0; l < argc; ++l) {
        myfile<< argv[l]<<" ";
    }
    myfile<<endl;
    myfile.close();
	return 0;
}
