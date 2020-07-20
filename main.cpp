/*----------------------- iPref problem ---------------
----------------------DBGroup@SUSTech, 2020------------
-------------------- tangloner@gmail.com-------------*/

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

// gobal variables
vector<vector<float>> HalfSpaces; // halfspace 
unordered_map<long int, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory

int objCnt = 0; // # of data objects
double totalIO = 0; // # of IO access
double totalSpaceCost = 0.0; // space cost (MB)
double treeSpacecost = 0.0;
double Spacecost = 0.0;

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
	// removed in iPref problem
	//const float sigma = atof(Param::read(argc, argv, "-R", ""));
	//const char* methodName = Param::read(argc, argv, "-m", "");

	int resultSize = 0;
	
	
	// data loading
	cout << "Load data points from file" << endl;
	float** PointSet = new float*[MAXPTS + 1];
	RtreeNodeEntry** p = new RtreeNodeEntry*[MAXPTS];
	fstream fpdata;
	fpdata.open(datafile, ios::in);
	while (true)
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

	//k-skyband 
	vector<long int> skyband;
	kskyband(dim, *rtree, skyband, PointSet, k);
	cout << skyband.size() << endl;
	
	
	return 0;
}
