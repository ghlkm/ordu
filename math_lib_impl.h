
#ifndef UTK_BV_UTK_MATH_LIB_IMPL_H
#define UTK_BV_UTK_MATH_LIB_IMPL_H

#include "vector_operator.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <algorithm>
#include "qp_solver.h"

template<typename V, typename VV>
double dist(const V &tmpT, const VV &PG){
    auto tmp=tmpT-PG;
    return sqrt(tmp*tmp);
}

template<typename V>
bool dominate(const V &v1, const V &v2) {
    /*
     * /tpara V vector, vector
     */
    assert(v1.size() == v2.size());
    return v1_dominate_v2(v1, v2, v1.size());
}

template<typename V, typename U>
bool v1_dominate_v2(const V &v1, const U &v2, size_t size) {
    /*
     * /tpara V array, pointer
     */
    for (auto i = 0; i < size; ++i) {
        if (v1[i] < v2[i]) {
            return false;
        }
    }
    return true;
}

template<typename V, typename VV>
double inflection_radius(const V &tmpT, const VV &PG, vector<c_float> &w,
                         int pdt_idx, int dimen, int k) {
    return inflection_radius(tmpT, PG, w, PG[pdt_idx], dimen, k);
}


template<typename V, typename VV, typename A>
double inflection_radius(const V &fetch, const VV &PG, vector<c_float> &w,
                         const A &pdt, int dimen, int k) {
    vector<c_float> radius;
    vector<c_float> h_ij(dimen, 0);
    for (auto &j:fetch) {
        if (v1_dominate_v2(PG[j], pdt, dimen)) {
            radius.push_back(INFINITY);
        } else {
            for (auto idx = 0; idx < dimen; ++idx) {
                h_ij[idx] = pdt[idx] - PG[j][idx];
            }
            radius.push_back(domin_r_ij3(w, h_ij));
        }
    }
    sort(radius.begin(), radius.end());
    return radius[radius.size() - k];
}

template<typename V, typename VV, typename A>
double inflection_radius(const V &fetch, const VV &PG, vector<c_float> &w,
                         const A &pdt, int dimen, int k, const double &rho,
                         multiset<double> &radius, unsigned int &cnt) {
    vector<c_float> h_ij(dimen, 0);
    auto iter = fetch.begin();
    cnt = 0;
    while (radius.size() < k && iter != fetch.end()) {
        ++cnt;
        if (v1_dominate_v2(PG[*iter], pdt, dimen)) {
            radius.insert(INFINITY);
        } else {
            for (auto idx = 0; idx < dimen; ++idx) {
                h_ij[idx] = pdt[idx] - PG[*iter][idx];
            }
            radius.insert(domin_r_ij3(w, h_ij));
        }
        ++iter;
    }
    for (; iter != fetch.end(); ++iter) {
        if (*radius.begin() >= rho) {
            break;
        } else {
            ++cnt;
            if (v1_dominate_v2(PG[*iter], pdt, dimen)) {
                radius.insert(INFINITY);
            } else {
                for (auto idx = 0; idx < dimen; ++idx) {
                    h_ij[idx] = pdt[idx] - PG[*iter][idx];
                }
                radius.insert(domin_r_ij3(w, h_ij));
            }
            radius.erase(radius.begin());
        }
    }
    return *radius.begin();
}

template<typename V, typename VV, typename A>
double inflection_radius(const V &fetch, const VV &PG, vector<c_float> &w,
                         const A &pdt, int dimen, int k, const double &rou) {
    /*
     * use rou to avoid unnecessary calculation
     */
    multiset<double> radius; // size with k
    unsigned int cnt;
    return inflection_radius(fetch, PG, w, pdt, dimen, k, rou, radius, cnt);
}

template<typename A, typename B>
double inflection_radius(vector<c_float> &w,
                         const A &the, const B &cmp, int dimen, int k) {
    if (v1_dominate_v2(cmp, the, dimen)) {
        return INFINITY;
    }
    vector<c_float> h_ij(dimen, 0);
    for (auto idx = 0; idx < dimen; ++idx) {
        h_ij[idx] = the[idx] - cmp[idx];
    }
    return domin_r_ij3(w, h_ij);
}


template<typename VV>
vector<int> k_skyband(VV &P, const int &k) {
    /*
     * /tpara vector<vector<double>> or vector<vector<double>>
     * the k-skyband contains thoes records that are dominated by fewer than k others
     */
    vector<int> do_cnt(P.size(), 0);
    vector<int> ret;
    for (auto i = 0; i < P.size(); ++i) {
        for (auto j = i + 1; j < P.size(); ++j) {
            if (do_cnt[i] >= k) {
                break;
            }
            if (dominate(P[i], P[j])) {
                ++do_cnt[j];
            } else if (dominate(P[j], P[i])) {
                ++do_cnt[i];
            }
        }
        if (do_cnt[i] < k) {
            ret.push_back(i);
        }
    }
    return ret;
}

template<typename VV>
vector<int> k_skyband(VV &P, const int &k, std::size_t p1, std::size_t p2) {
    /*
     * TODO: be careful to use this, begin from 1!!!!
     * /tpara vector<vector<double>> or vector<vector<double>>
     * the k-skyband contains thoes records that are dominated by fewer than k others
     */
    vector<int> do_cnt(p1, 0);
    vector<int> ret;
    for (auto i = 1; i <p1; ++i) {
        for (auto j = i + 1; j < p1; ++j) {
            if (do_cnt[i] >= k) {
                break;
            }
            if (v1_dominate_v2(P[i], P[j], p2)) {
                ++do_cnt[j];
            } else if (v1_dominate_v2(P[j], P[i], p2)) {
                ++do_cnt[i];
            }
        }
        if (do_cnt[i] < k) {
            ret.push_back(i);
        }
    }
    return ret;
}

template<typename ITER>
double sum(const ITER &begin, const ITER &end){
    double ret=0;
    for(auto i=begin;i!=end;++i){
        ret+=*i;
    }
    return ret;
}


extern qp_solver *qp_ptr;
template<typename V>
double domin_r_ij3(const V &w, const V &h_ij) {
    //TODO this can be lazy only update h_ij sometimes
    return qp_ptr->update_w_h_solve(w, h_ij);
}

template<typename V>
inline V proj(const V &u, const V &v) {
    return (u * v) / (u * u) * u;
}

template<typename V>
inline c_float vector_length(V &v) {
    c_float ret = c_float();
    for (auto &i:v) {
        ret += i * i;
    }
    return sqrt(ret);
}

template<typename VV>
VV gram_schmidt_process(const VV &input) {
    /*
     * be careful when:
     * e.g.
     *  input={a_0, a_1, a_2}
     *  if a_2 is a linear combination of {a_0, a_1},
     *   then ret[2] will be a vector of 0s
     */
    assert(!input.empty());
    VV ret(input.size());
    // begin generate orthogonal vectors
    for (int i = 0; i < input.size(); ++i) {
        ret[i] = input[i];
        for (int j = 0; j < i; ++j) {
            ret[i] -= proj(ret[j], input[i]);
        }
    }
    // begin normalize result
    for (int k = 0; k < ret.size(); ++k) {
        c_float norm_of_v = vector_length(ret[k]);
        for (auto &ele:ret[k]) {
            ele /= norm_of_v;
        }
    }
    return ret;
}

template<typename INTEGER>
vector<vector<c_float>> gen_r_domain_basevec(INTEGER dim) {
    /*
     * \tpara INTEGER bit, byte, char, short, int, long, size_t, unsigned
     */
    vector<vector<c_float>> u(dim);
    u[0] = vector<c_float>(dim, 1.0);
    for (int i = 1; i < dim; ++i) {
        u[i] = vector<c_float>(dim);
        u[i][i] = 1.0;
    }
    vector<vector<c_float>> e = gram_schmidt_process(u);
    e.erase(e.begin());
    return e;
}

template<typename VV, typename V>
void gen_r_domain_vec_dfs(V &&cur, int next,
                          VV &ret, VV &e) {
    if (next == e.size()) {
        ret.push_back(cur);
    } else {
        gen_r_domain_vec_dfs(cur + e[next], next + 1, ret, e);
        gen_r_domain_vec_dfs(cur - e[next], next + 1, ret, e);
    }
}

template<typename INTEGER>
vector<vector<c_float>> gen_r_domain_vec(INTEGER dim) {
    vector<vector<c_float>> e = gen_r_domain_basevec(dim); // size with d-1
    vector<vector<c_float>> ret; // size with 2^(d-1)
    vector<c_float> cur(dim);
    gen_r_domain_vec_dfs(cur, 0, ret, e);
    return ret;
}

extern vector<vector<c_float>> g_r_domain_vec;

template<typename VF, typename VI, typename VV, typename FLOAT>
bool isR_skyband(const VV &PG, const VI&vs, const VF &opt, const VF &w, const FLOAT &rho, int k) {
    int r_dominate_cnt=0;
    for (const int&v:vs) {
        if(v2_r_dominate_v1(opt, PG[v], w, g_r_domain_vec, rho)){
            ++r_dominate_cnt;
            if(r_dominate_cnt>=k){
                break;
            }
        }
    }
    return r_dominate_cnt<k;
}



template<typename V>
inline bool v1_v2_nonDominate(const V&v1, const V&v2, std::size_t size){
    return !v1_dominate_v2(v1, v2, size) && !v1_dominate_v2(v2, v1, size);
}

template<typename V>
inline bool v1_v2_nonDominate(const V&v1, const V&v2){
    return v1_v2_nonDominate(v1, v2, v1.size());
}



template<typename V, typename VV, typename FLOAT>
bool v2_r_dominate_v1(const V &v1, const V &v2, const V &w, const VV &r_domain_vec, const FLOAT &rho) {
    for (const V &v:r_domain_vec) {
        double atc_rho=rho;
        for (int i = 0; i < v.size(); ++i) {
            if(v[i]<0){
                atc_rho=min(atc_rho, -w[i]/v[i]); // in case of w[i] + \rho * v[i] <0 or >1
            }
        }
        V tmp_w = w + atc_rho * v;
        if (v1 * tmp_w < v2 * tmp_w) {
            return false;
        }
    }
    return true;
}

template<typename INT, typename VV>
vector<INT> computeTopK(const int dim, VV &PG, vector<INT> &skyband, vector<float>& weight, int k)
{
    multimap<float, INT> heap;
    for (int i = 0; i < skyband.size(); i++)
    {
        float score = 0;
        for (int d = 0; d < dim; d++)
        {
            score += PG[skyband[i]][d] *weight[d];
        }

        if (heap.size() < k)
        {
            heap.emplace(score, skyband[i]);
        }
        else if (heap.size() == k && heap.begin()->first < score)
        {
            heap.erase(heap.begin());
            heap.emplace(score, skyband[i]);
        }
    }

    vector<INT> topkRet;
    for (auto heapIter = heap.rbegin(); heapIter != heap.rend(); ++heapIter)
    {
        topkRet.push_back(heapIter->second);
    }
    return topkRet;
}

template<typename INT, typename VV, typename  FLOAT>
vector<INT> computeTopK_Extend(const int dim, VV &PG, vector<INT> &skyband, vector<FLOAT>& weight, int k,
        vector<INT> &retSw, vector<INT> &retSnr)
{
    vector<pair<INT, double>> rec;
    rec.reserve(skyband.size());
    for(INT r: skyband){
        double s=0;
        for (int d = 0; d < dim; d++){
            s += PG[r][d] *weight[d];
        }
        rec.emplace_back(r, s);
    }
    sort(rec.begin(), rec.end(), [](auto &a, auto& b){
        return a.second>b.second;
    });
    double kthScore=rec[k-1].second;
    int position =k;
    vector<INT> topkRet;
    for (int i = 0; i < k; ++i) {
        topkRet.push_back(rec[i].first);
    }
    while(rec.size()>position){
        if(abs(rec[position].second-kthScore)>1e-6){
            break;
        }
        topkRet.push_back(rec[position].first);
        position++;
    }

    retSw=no_dominate_set(topkRet, PG, dim);
    if(position<rec.size()){
        vector<int> tmp(rec.size()-position);
        for (int i = position; i <rec.size() ; ++i) {
            tmp[i-position]=rec[i].first;
        }
        retSnr=non_dominate_set(tmp, PG, dim);
    }
    return topkRet;
}

template<typename V1, typename V2>
inline double dot(V1 &v1, V2 &v2){
    return dot(v1, v2, v1.size());
}

template<typename V1, typename V2>
inline double dot(V1 &v1, V2 &v2, std::size_t size){
    double ret=0;
    for (int i = 0; i < size; ++i) {
        ret+=v1[i]*v2[i];
    }
    return ret;
}

template<typename FF>
vector<int> non_dominate_set(const vector<int> &candidate, const FF* P, std::size_t dim){
    vector<int> ret;
    for (int l: candidate) {
        bool flag=true;
        for (int i: candidate) {
            if(l!=i){
                if(v1_dominate_v2(P[i], P[l], dim)) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag){
            ret.push_back(l);
        }
    }
    return ret;
}

template<typename FF>
vector<int> no_dominate_set(const vector<int> &candidate, const FF* P, std::size_t dim){
    vector<int> ret;
    for (int l: candidate) {
        bool flag=true;
        for (int i: candidate) {
            if(l!=i){
                if(v1_dominate_v2(P[l], P[i], dim)) {
                    flag = false;
                    break;
                }
            }
        }
        if(flag){
            ret.push_back(l);
        }
    }
    return ret;
}

#endif //UTK_BV_UTK_MATH_LIB_IMPL_H

