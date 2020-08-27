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
    inline void init_all(const vector<FLOAT>& w, const vector<FLOAT>& h_ij);
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
