//
// Created by likeming on 2020/8/24.
//

#ifndef TEST_QP_SOLVER_H
#define TEST_QP_SOLVER_H
#include  <numeric>
#include "osqp.h"
#include <vector>
#include <cassert>

using namespace std;


class qp_solver{
    c_int  dim;
    c_float w_L2;
    OSQPWorkspace *work;
    OSQPSettings  *settings;
    OSQPData      *data;
    vector<c_float> P_x;
    vector<c_int> P_i;
    vector<c_int> P_p;
    vector<c_float> A_x;
    vector<c_int> A_i;
    vector<c_int> A_p;
    vector<c_float> q;
    vector<c_float> l;
    vector<c_float> u;

    inline void destroy(){
        if (data) {
            if (data->A) c_free(data->A);
            if (data->P) c_free(data->P);
            c_free(data);
        }
        if (settings) c_free(settings);
    }

    template<typename FLOAT>
    inline void init_all(const vector<FLOAT>& w, const vector<FLOAT>& H);

    inline void init_all(const vector<float>& w, const vector<vector<double>>& H){
        /*
         * csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p)
         * in OSQP,
         * QP problem's "P" and "A" are represented as csc format sparse matrix
         * for more detail about transform dense_matrix into csc_matrix
         * see https://en.wikipedia.org/wiki/Sparse_matrix
         *
         * The OSQP example can be seen from here:
         * https://osqp.org/docs/examples/setup-and-solve.html
         */
        assert(!w.empty());
        dim=w.size();
        w_L2=0;
        for (float i : w) {
            w_L2+=i*i;
        }

        P_x=vector<c_float>(dim, 1.0);
        c_int P_nnz = dim;
        P_i=vector<c_int>(dim);
        iota(P_i.begin(), P_i.end(), 0); //{0, 1, 2, ..., dim-1}
        P_p=vector<c_int>(dim+1);
        iota(P_p.begin(), P_p.end(), 0); //{0, 1, 2, ..., dim-1}

        q=vector<c_float>(w.begin(), w.end());
        for (c_float &qi:q) {
            qi=-qi;
        }

        A_x=vector<c_float>((2+H.size())*dim);
        c_int A_x_idx=0, h_idx=0;
        while(A_x_idx<A_x.size()){
            A_x[A_x_idx++]=1.0;
            for (int i = 0; i < H.size(); ++i) {
                A_x[A_x_idx++]=H[i][h_idx];
            }
            h_idx++;
            A_x[A_x_idx++]=1.0;
        }
        c_int A_nnz = A_x.size();
        A_i=vector<c_int>(3*dim);
        c_int A_i_idx=0, h_iidx=0;
        while(A_i_idx<A_i.size()){
            A_i[A_i_idx++]=0;
            for (int i = 0; i < H.size(); ++i) {
                A_i[A_i_idx++]=i+1;
            }
            A_i[A_i_idx++]=h_iidx+1+H.size();
            ++h_iidx;
        }
        A_p=vector<c_int>(dim+1);
        c_int A_p_idx=0, h_pidx=0;
        while(A_p_idx+1<A_p.size()){
            A_p[A_p_idx+1]=A_p[A_p_idx]+2+H.size();
            ++A_p_idx;
        }

        l=vector<c_float>(dim+1+H.size(), 0);
        l[0]=1.0;
        for (int j = 0; j <H.size() ; ++j) {
            l[j+1]=-INFINITY;
        }
        u=vector<c_float>(dim+1+H.size(), 1.0);
        for (int j = 0; j <H.size() ; ++j) {
            u[j+1]=0;
        }

        settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        data     = (OSQPData *)c_malloc(sizeof(OSQPData));
        data->n = dim;
        data->m = dim+1+H.size(); //constrain number
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x.data(), P_i.data(), P_p.data());
        data->q = q.data();
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x.data(), A_i.data(), A_p.data());
        data->l = l.data();
        data->u = u.data();


        if (settings) {
            osqp_set_default_settings(settings);
            settings->verbose=false; // keep quiet
            settings->alpha = 1.0; // Change alpha parameter
        }
        // Setup workspace
        osqp_setup(&work, data, settings);

    }
public:
    qp_solver() = default;

    template<typename INT>
    qp_solver(INT dim);

    qp_solver(qp_solver &other){
        if(&other!=this){
            dim=other.dim;
            c_int P_nnz = dim;
            w_L2=other.w_L2;
            P_x=other.P_x;
            P_i=other.P_i;
            P_p=other.P_p;
            A_x=other.A_x;
            c_int A_nnz = A_x.size();
            A_i=other.A_i;
            A_p=other.A_p;
            q=other.q;
            l=other.l;
            u=other.u;
            settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
            data     = (OSQPData *)c_malloc(sizeof(OSQPData));
            data->n = dim;
            data->m = dim+2; //constrain number
            data->P = csc_matrix(data->n, data->n, P_nnz, P_x.data(), P_i.data(), P_p.data());
            data->q = q.data();
            data->A = csc_matrix(data->m, data->n, A_nnz, A_x.data(), A_i.data(), A_p.data());
            data->l = l.data();
            data->u = u.data();
            if (settings) {
                osqp_set_default_settings(settings);
                settings->verbose=false; // keep quiet
                settings->alpha = 1.0; // Change alpha parameter
            }
            // Setup workspace
            osqp_setup(&work, data, settings);
        }
    }

    template<typename FLOAT>
    explicit qp_solver(const vector<FLOAT>& w);

    template<typename FLOAT>
    qp_solver(const vector<FLOAT>& w, const vector<FLOAT>& h_ij);

    qp_solver(const vector<float>& w, const vector<vector<double>>& H){
        init_all(w, H);
    }

    template<typename FLOAT>
    inline void update_w(const vector<FLOAT> &w);

    template<typename FLOAT>
    inline void update_h(const vector<FLOAT> &h);

    template<typename FLOAT>
    inline FLOAT solve_update_h(const vector<FLOAT> &h);

    template<typename FLOAT>
    inline void update_w_h(const vector<FLOAT> &w, const vector<FLOAT>&h_ij);

    template<typename FLOAT>
    inline FLOAT update_w_h_solve(const vector<FLOAT> &w, const vector<FLOAT> &h);

    inline c_float solve(){
        osqp_solve(work);
        c_float ret=2*work->info->obj_val+w_L2;
        return ret>=0?sqrt(ret):0;
    }

    qp_solver& operator=(qp_solver &&other) noexcept {
        if(&other!=this){
            dim=other.dim;
            c_int P_nnz = dim;
            w_L2=other.w_L2;
            P_x=other.P_x;
            P_i=other.P_i;
            P_p=other.P_p;
            A_x=other.A_x;
            c_int A_nnz = A_x.size();
            A_i=other.A_i;
            A_p=other.A_p;
            q=other.q;
            l=other.l;
            u=other.u;
            settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
            data     = (OSQPData *)c_malloc(sizeof(OSQPData));
            data->n = dim;
            data->m = dim+2; //constrain number
            data->P = csc_matrix(data->n, data->n, P_nnz, P_x.data(), P_i.data(), P_p.data());
            data->q = q.data();
            data->A = csc_matrix(data->m, data->n, A_nnz, A_x.data(), A_i.data(), A_p.data());
            data->l = l.data();
            data->u = u.data();
            if (settings) {
                osqp_set_default_settings(settings);
                settings->verbose=false; // keep quiet
                settings->alpha = 1.0; // Change alpha parameter
            }
            // Setup workspace
            osqp_setup(&work, data, settings);
        }
        return *this;
    }

    ~qp_solver(){
        destroy();
    }
};

#include "qp_solver_impl.h"

#endif //TEST_QP_SOLVER_H
