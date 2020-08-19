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

// gobal variables
vector<vector<float>> HalfSpaces; // halfspace 
unordered_map<long int, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory

int objCnt = 0; // # of data objects
double totalIO = 0; // # of IO access
double totalSpaceCost = 0.0; // space cost (MB)
double treeSpacecost = 0.0;
double Spacecost = 0.0;

template<typename V>
inline V proj(const V& u, const V &v){
    /*
     * <u, v> / <u, u> * u
     */
    float dot=0;
    assert(u.size()==v.size());
    for (int i = 0; i < u.size(); ++i) {
        dot+=u[i]*v[i];
    }
    float u_mod2=0;
    for (int i = 0; i < u.size(); ++i) {
        u_mod2+=u[i]*u[i];
    }
    float length=dot/u_mod2;
    V ret(u);
    for(int i = 0; i < u.size(); ++i){
        ret[i]*=length;
    }
    return ret;
}

template<typename V>
float domin_r_ij2(const V &w, const V &h_ij) {
    //TODO this fun is in utk_math_lib, merge utk_math_lib later
    //TODO when vertical point is out of domain, use QP solver
    /*
     *     V ones(h_ij.size(), 1);
     *     V ones_n=ones-proj(h_ij, ones);
     *     V d=proj(h_ij, WA)+proj(ones_n, WA);
     *     return sqrt(d*d);
     */
    assert(h_ij.size()>=2); //dim should be no less than 2
    V ones(h_ij.size(), 1.0);
    V ones_proj_on_hij=proj(h_ij, ones); //proj(h_ij, ones)
    for (int i = 0; i < ones.size(); ++i) {
        ones[i]-=ones_proj_on_hij[i];   // V ones_n=ones-proj(h_ij, ones);
    }
    V point_of_hij(h_ij.size(), 0.0);// point A
    float sum=-h_ij[0]+h_ij[1];
    if(sum!=0){
        point_of_hij[0]= h_ij[1] / sum;
        point_of_hij[1]= -h_ij[0] / sum; // normalize A such that \Sigma A[i]=1
    }else{
        point_of_hij[0]=0.5;
        point_of_hij[1]=0.5;
    }
    // point A init finish
    for (int k = 0; k < w.size(); ++k) {
        point_of_hij[k]-=w[k]; //A-->WA
    }
    V A_proj_on_hij=proj(h_ij, point_of_hij); //proj(h_ij, WA)
    V A_proj_on_onesn=proj(ones, point_of_hij); //proj(ones_n, WA)
    for (int j = 0; j < A_proj_on_hij.size(); ++j) {
        A_proj_on_hij[j]+=A_proj_on_onesn[j]; //d=proj(h_ij, WA)+proj(ones_n, WA)
    }
    float WA_proj_on_hij_ones_mod2=0;
    for (int j = 0; j < A_proj_on_hij.size(); ++j) {
        WA_proj_on_hij_ones_mod2+= A_proj_on_hij[j] * A_proj_on_hij[j]; // d*d
    }
    return sqrt(WA_proj_on_hij_ones_mod2); //return sqrt(d*d)
}

int main(const int argc, const char** argv)
{
	cout.precision(4);
	cout << "iPref Problem (Size-constrained R-kSkyband/UTK )" << endl;
	clock_t at, ad;
	// parameter parser
	cout << "Parse Parameters" << endl;
	if (argc == 1)
	{
		helpmsg(argv[0]);
		return -1;
	}

	int dim = atoi(Param::read(argc, argv, "-d", ""));
	int k = atoi(Param::read(argc, argv, "-k", ""));
	const char* datafile = Param::read(argc, argv, "-f", "");
	const char* indexfile = Param::read(argc, argv, "-i", "");
	const int X = atoi(Param::read(argc, argv, "-X", ""));
	const char* methodName = Param::read(argc, argv, "-m", "");
    int w_num = atoi(Param::read(argc, argv, "-w", ""));
    int p_num = atoi(Param::read(argc, argv, "-p", ""));// how many option want to run
    const char *w_file = Param::read(argc, argv, "-W", "");
    std::vector<std::vector<float>> ws(w_num, std::vector<float>(dim));
    fstream fpdata;
    fpdata.open(w_file, ios::in);
    for (auto i = 0; i < ws.size(); ++i) {
        for (auto j = 0; j < ws[i].size(); ++j) {
            fpdata >> ws[i][j];
        }
    }
    fpdata.close();

    int resultSize = 0;
	
	
	// data loading
	cout << "Load data points from file" << endl;
	float** PointSet = new float*[MAXPTS + 1];
	RtreeNodeEntry** p = new RtreeNodeEntry*[MAXPTS];
//	fstream fpdata;
	fpdata.open(datafile, ios::in);
	while (p_num--)
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
	// aggregateRecords(*rtree);
	// cout << "[Aggregate Rtree done]" << endl;

	// at = clock();

	// baseline algorithm 
	// (1) compute k-skyband set SK
	// (2) for each pi, suppose T_pi is the incomparable set (a) p_j is incomparable with p_i in SK, (2) p_j w > p_i w
	// (3) for each p in T_pi, compute the distance from point w to hyperplane H_{i,j} equation: d = \frac{a1 w1 +  a2 w2 + ... ad_1 wd-1+D}{\sqrt{a1^2+a2^2+ ... + ad-1^2}}
	// (4) compute the inflection point of pi: the k-th laregst value in step (3), denotes as inf_pi
	// (5) pi's rksykband interval: inf_pi to infinity
	// (6) the radius rho is the T-th minimum value in all pi in SK
	
	if (strcmp(methodName, "BB") == 0)
	{
	    vector<long int> skyband;
	    kskyband(dim, *rtree, skyband, PointSet, k); // step (1)
	    cout << skyband.size() << endl;
        for(auto iteration=0;iteration<w_num;++iteration) {
            // TODO Actually, k-skyband don't needed for each iteration for diff weights
            // weight vector for testing
            vector<float> userpref(ws[iteration].begin(), ws[iteration].begin() + (ws[iteration].size() - 1));
            vector<float> w(ws[iteration].begin(), ws[iteration].end());
            cout<<"w: ";
            for (auto i:userpref) {
                cout << i << ", ";
            }
            cout <<w.back()<<endl;
            //k-skyband
            auto begin = chrono::steady_clock::now();
            vector<long int> topKRet = computeTopK(dim, PointSet, skyband, userpref, k);
            assert(topKRet.size() == k);
            vector<pair<long int, float>> interval;
            for (int i = 0; i < topKRet.size(); i++) {
                interval.emplace_back(topKRet[i], 0);
            }

            for (int ski = 0; ski < skyband.size(); ski++) {
                if (find(topKRet.begin(), topKRet.end(), skyband[ski]) != topKRet.end()) {
                    continue;
                }
                vector<float> radiusSKI;
                vector<long int> incompset;
                vector<long int> dominatorSet;
                for (int pj = 0; pj < skyband.size(); pj++) {
                    if (ski == pj) {
                        continue;
                    }
                    if(radiusSKI.size()>=k)// prune
                        continue;
                    if (IsPjdominatePi(dim, PointSet, skyband[ski], skyband[pj]))
                    {
                        radiusSKI.push_back(FLT_MAX);
                    }
                    else if (incomparableset(PointSet, skyband[ski], skyband[pj], userpref)) // step (2)
                    {
                        incompset.push_back(skyband[pj]);
                    }
                    else
                    {
//                        assert(find(topKRet.begin(), topKRet.end(), skyband[pj]) == topKRet.end());
                    }
                }
                // here we need a function to compute the inflection radius of option pi
                for (int inpi = 0; inpi < incompset.size(); inpi++) {
                    vector<float> HS(PointSet[skyband[ski]], PointSet[skyband[ski]] + dim);
                    for (auto i = 0; i < dim; ++i) {
                        HS[i] -= PointSet[incompset[inpi]][i];
                    }
                    float tmpdis = domin_r_ij2(w, HS);
                    radiusSKI.push_back(tmpdis);
                    //compute the hyperplane of pi and pj
//				vector<float> tmpHS = computePairHP(dim, PointSet, skyband[ski], incompset[inpi]);
//				//compute the distance from w to hyperplane.
//				float tmpDis = computeDis(tmpHS, userpref);
//				radiusSKI.push_back(tmpDis);
                }
                sort(radiusSKI.begin(), radiusSKI.end());
                assert(radiusSKI.size() >= k);
                if (radiusSKI.size() >= k) {
                    interval.emplace_back(skyband[ski], radiusSKI[radiusSKI.size() - k]);
                }
            }
            assert(interval.size() >= X);
            sort(interval.begin(), interval.end(), sortbysec);
            for (auto i=0;i<k;++i) {
                interval[i].first=topKRet[i];
            }
            auto end = chrono::steady_clock::now();
            chrono::duration<double> t = end - begin;
//            cout << "run time:" << t.count() << endl;
            int cnt = 0;
            for (auto i = interval.begin(); i != interval.end(); ++i) {
                cout << i->second << ", " << i->first << endl;
                ++cnt;
                if (cnt > X)
                    break;
            }
            cout<<"\n\n";
//            cout << "The inflection radius is: " << interval[X].second << "\n";
        }
	}
	
	


	

	// optimized algorithm
	// (1) revised BBS for w, compute top-k result, store in T directly
	// (2) for each pi (uncertain status), compute its inflection distance inf_pi, and push it into a max-heap Q (inf_pi as key)
	// (3) keep fetching until max-heap with size X-k+1, pop the top node out, and initalize p0 = inf_ptop
	// (4) test the rest candidates by r-dominanace, instead of BBS, thus fetching r-skyband options
	// (5) shrinks p0 quicly, simulating a scanning from infinity to 0
	// (6) fetching terminates none of options can be the r-skyband
	// (7) append Q to T, this is the final result.

	// inclremental version, without known X
	// We do not have exact X, we need tell the user the radius rho and its corresponding T
	// It is similar to optimized algorithm, however, it computes incrementally, from rho = 0 to infinity, the size T is from k to k-skyband.

	return 0;
}
