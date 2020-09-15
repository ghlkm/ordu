//
// Created by 12859 on 2020/7/15.
//
#ifndef UTK_BV_UTK_MATH_LIB_IMPL_H
#define UTK_BV_UTK_MATH_LIB_IMPL_H

#include "utk_vector.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <algorithm>
#include "qp_solver.h"

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
                         const A &pdt, int dimen, int k, const double &rou,
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
        if (*radius.begin() >= rou) {
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

template<typename V>
double domin_r_ij(const V &w, const V &h_ij) {
    /*
     * \tpara V utk_vector
     *
     * return:
     *   r_ij
     */
    /*
     * n is one of the normal vector of hyperplane h_ij, actually the coefficients of h_ij is n
     *
     * cos(\theta) |w| = w \cdot n / |n|
     *
     * w_ij_ = w- n cos(\theta) |w|/|n| = w- w \cdot n \cdot n / |n|^2
     *
     * w_ij = \lambda w_ij_, \Sigma w_ij =1
     * --> w_ij = w_ij_ / \Sigma w_ij_
     *
     * inflection distance r_ij^2 = (w-w_ij) \cdot (w-w_ij)
     */

    // h_ij and w can be with dif direction because of
    // the resulting w_ij_=w-w*h_ij*h_ij/(h_ij*h_ij)
    V w_ij_ = w - w * h_ij * h_ij / (h_ij * h_ij); // dot product can not change operate order!!!
    double mod1 = 0;
    for (auto i = 0; i < w_ij_.size(); ++i) {
        mod1 += w_ij_[i];
    }
    for (auto i = 0; i < w_ij_.size(); ++i) {
        w_ij_[i] /= mod1;
    }
    // now w_ij=w_ij_
    V &&d = w - w_ij_;
    return d * d; // rou^2;
}


template<typename V>
double domin_r_ij2(const V &w, const V &h_ij) {
    // sample from h_ij
    assert(h_ij.size() >= 2);
    V A(h_ij.size(), 0.0);
    float sum = -h_ij[0] + h_ij[1];
    if (sum != 0) {
        A[0] = h_ij[1] / sum;
        A[1] = -h_ij[0] / sum;
    } else {
        A[0] = 0.5;
        A[1] = 0.5;
    }
    V WA = A - w;
    V ones(h_ij.size(), 1.0);
    V ones_n = ones - proj(h_ij, ones);
    V d = proj(h_ij, WA) + proj(ones_n, WA);
    return sqrt(d * d);
}

extern qp_solver *qp_ptr;
template<typename V>
double domin_r_ij3(const V &w, const V &h_ij) {
    return qp_ptr->solve_update_h(h_ij);
}
template<typename V>
inline V proj(const V &u, const V &v) {
    return (u * v) / (u * u) * u;
}

template<typename V>
inline c_float L2_norm(V &v) {
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
        c_float l2_norm = L2_norm(ret[k]);
        for (auto &ele:ret[k]) {
            ele /= l2_norm;
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
    return ret; // each element of it is with L2_norm=sqrt(d-1)
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

#endif //UTK_BV_UTK_MATH_LIB_IMPL_H

