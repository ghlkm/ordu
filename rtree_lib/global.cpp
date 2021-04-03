#include "global.h"
using namespace std;
void helpmsg(const char* pgm)
{
	cout << "Suggested arguments:" << endl;
	cout << "> " << pgm << " ";
	cout << "-d 2 -f raw_data.txt -i index_file.idx -q query_data.txt -k 0.1 -m CTA" << endl;
	cout << "explanations:" << endl;
	cout << "-d: data dimensionality, e.g. 2" << endl;
	cout << "-f: a tab delimited file, format:" << endl;
	cout << "    id xmin ymin [zmin] xmax ymax [zmax]" << endl;
	cout << "-q: query file, format:" << endl;
	cout << "    id xmin ymin xmax ymax" << endl;
	cout << "-k: monochromatic reverse top-k" << endl;
	cout << "   >= 1: exact k" << endl;
	cout << "-i: index file" << endl;
	cout << "-m: method, e.g., CTA, CTA+, ACTA, ACTA*" << endl;
}

double dataSize(int objCnt, int dim)
{
	double tmpMemUsage = sizeof(float)* 2 * dim;
	return tmpMemUsage / MB;
}

void GenString(long int stringLen, long int HammingDist, long int curLen, long int start, vector<char>& hammstr, multimap<int, vector<char>>& binString)
{
	typedef multimap<int, vector<char>>::value_type VT;

	if (curLen < 0)
		return;

	for (long int i = start; i < stringLen; i++)
	{
		for (long int j = start; j < i; j++)
			hammstr.push_back('0');
		hammstr.push_back('1');
		GenString(stringLen, HammingDist, curLen - 1, i + 1, hammstr, binString);

		if (curLen == 0)
		{
			for (long int j = i + 1; j < stringLen; j++)
			{
				hammstr.push_back('0');
			}
			vector<char> tmphamstr = hammstr;
			binString.insert(VT(HammingDist, tmphamstr));

			for (long int j = i + 1; j < stringLen; j++)
			{
				hammstr.pop_back();
			}
		}

		hammstr.pop_back();
		for (long int j = start; j < i; j++)
		{
			hammstr.pop_back();
		}
	}
}


