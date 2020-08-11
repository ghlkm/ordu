#include "iPref.h"

vector<long int> computeTopK(const int dim, float* PG[], vector<long int> skyband, vector<float>& weight, int k)
{
	multimap<float, long int> heap;
	multimap<float, long int>::iterator heapIter;

	for (int i = 0; i < skyband.size(); i++)
	{
		float score = 0;
		float l2_sumWeight = 0;
		for (int d = 0; d < dim-1; d++)
		{
			score += (PG[skyband[i]][d] + PG[skyband[i]][d + dim]) / 2*weight[d];
            l2_sumWeight += weight[d]*weight[d];
		}
        score += (PG[skyband[i]][dim - 1] + PG[skyband[i]][dim - 1 + dim]) / 2 * sqrt(1 - l2_sumWeight);

		if (heap.size() < k)
		{
			heap.insert(multimap<float, long int>::value_type(score, skyband[i]));
			//could use emplace， e.g. “heap.emplace(score, skyband[i])”
		}
		else if (heap.size() == k && heap.begin()->first < score)
		{
			heapIter = heap.begin();
			heap.erase(heapIter);
			heap.insert(multimap<float, long int>::value_type(score, skyband[i]));
		}
	}

	vector<long int> topkRet;
	for (heapIter = heap.begin(); heapIter != heap.end(); heapIter++)
	{
		topkRet.push_back(heapIter->second);
	}
	return topkRet;
}

bool IsPjdominatePi(const int dimen, float* PG[], long int pi, long int pj)
{
	float pid;
	float pjd;
	int count = 0;
	for (int d = 0; d < dimen; d++)
	{
		pid = (PG[pi][d] + PG[pi][dimen + d]) / 2.0;
		pjd = (PG[pj][d] + PG[pj][dimen + d]) / 2.0;
		if (pjd >= pid)
			count++;
	}


	if (count == dimen)
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool incomparableset(float* PG[], long int pi, long int pj, vector<float>& weight)
{
	int dimen = weight.size() + 1;
	float l2_wd = 0;
	float spi = 0, spj = 0;
	int cpos = 0;
	int cneg = 0;

	float piv[MAXDIMEN];
	float pjv[MAXDIMEN];
	for (int d = 0; d < dimen; d++)
	{
		piv[d] = (PG[pi][d] + PG[pi][dimen + d]) / 2.0;
		pjv[d] = (PG[pj][d] + PG[pj][dimen + d]) / 2.0;
	}

	for (int i = 0; i < dimen - 1; i++)
	{
		spi += piv[i] * weight[i];
		spj += pjv[i] * weight[i];
        l2_wd += weight[i]*weight[i];

		if (piv[i] <= pjv[i])
		{
			cneg++;
		}
		else
		{
			cpos++;
		}
	}

	spi += sqrt(1 - l2_wd) * piv[dimen - 1];
	spj += sqrt(1 - l2_wd) * pjv[dimen - 1];

	if (piv[dimen - 1] <= pjv[dimen - 1])
	{
		cneg++;
	}
	else
	{
		cpos++;
	}

	if (spj >= spi && cneg != 0 && cpos != 0)
		return true;
	else
		return false;
}


vector<float> computePairHP(const int dimen, float* PG[], long int pi, long int pj)
{
	vector<float> retHS;
	float piv[MAXDIMEN];
	float pjv[MAXDIMEN];
	for (int d = 0; d < dimen; d++)
	{
		piv[d] = (PG[pi][d] + PG[pi][dimen + d]) / 2.0;
		pjv[d] = (PG[pj][d] + PG[pj][dimen + d]) / 2.0;
	}

	float pj_d = pjv[dimen - 1];
	float pi_d = piv[dimen - 1];
	for (int d = 0; d < dimen - 1; d++)
	{
		retHS.push_back((pjv[d] - pj_d) - (piv[d] - pi_d));
	}
	retHS.push_back(pi_d - pj_d);
	return retHS;
}


float computeDis(vector<float> tmpHS, vector<float> userpref)
{
	float normVec = 0;
	float dotRet = 0;
	for (int i = 0; i < tmpHS.size() - 1; i++)
	{
		dotRet += userpref[i] * tmpHS[i];
		normVec += tmpHS[i] * tmpHS[i];
	}
//  wrong
//	return abs(dotRet - tmpHS[tmpHS.size() - 1]) / normVec;

//  right for min distance of point to hyperplane
    return abs(dotRet - tmpHS[tmpHS.size() - 1]) / sqrt(normVec);
}

bool sortbysec(const pair<long int, float> &a, const pair<long int, float> &b)
{
	return (a.second < b.second);
}