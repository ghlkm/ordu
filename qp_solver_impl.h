
#include  <numeric>
#include "osqp.h"


template<typename FLOAT>
inline void qp_solver::init_all(const vector<FLOAT>& w, const vector<FLOAT>& h_ij){
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
    for (FLOAT i : w) {
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

    A_x=vector<c_float>(3*dim);
    c_int A_x_idx=0, h_idx=0;
    while(A_x_idx<A_x.size()){
        A_x[A_x_idx++]=1.0;
        A_x[A_x_idx++]=h_ij[h_idx++];
        A_x[A_x_idx++]=1.0;
    }
    c_int A_nnz = A_x.size();
    A_i=vector<c_int>(3*dim);
    c_int A_i_idx=0, h_iidx=0;
    while(A_i_idx<A_i.size()){
        A_i[A_i_idx++]=0;
        A_i[A_i_idx++]=1;
        A_i[A_i_idx++]=h_iidx+2;
        ++h_iidx;
    }
    A_p=vector<c_int>(dim+1);
    c_int A_p_idx=0, h_pidx=0;
    while(A_p_idx+1<A_p.size()){
        A_p[A_p_idx+1]=A_p[A_p_idx]+3;
        ++A_p_idx;
    }

    l=vector<c_float>(dim+2, 0);
    l[0]=1.0;
    u=vector<c_float>(dim+2, INFINITY);
    u[0]=1.0;
    u[1]=0.0;

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

template<typename INT>
qp_solver::qp_solver(INT dim){
    vector<c_float> h_ij(dim);
    vector<c_float> w(dim);
    init_all(w, h_ij);
}

template<typename FLOAT>
qp_solver::qp_solver(const vector<FLOAT>& w){
    vector<FLOAT> h_ij(w.size());
    init_all(w, h_ij);
}

template<typename FLOAT>
qp_solver::qp_solver(const vector<FLOAT>& w, const vector<FLOAT>& h_ij){
    init_all(w, h_ij);
}



template<typename FLOAT>
inline void qp_solver::update_h(const vector<FLOAT>&h){
    c_int h_idx=1, h_ij_idx=0;
    while (h_idx<A_x.size()){
        A_x[h_idx]=h[h_ij_idx];
        h_idx+=3;
        h_ij_idx+=1;
    }
    osqp_update_A(work, A_x.data(), OSQP_NULL, A_x.size());
}


template<typename FLOAT>
inline FLOAT qp_solver::solve_update_h(const vector<FLOAT> &h){
    this->update_h(h);
    //the real dominate radius is
    // 2 * work->info->obj_val + w^T \cdot w
    FLOAT ret= this->qp_solve();
    return ret;
}

template<typename FLOAT>
inline void qp_solver::update_w(const vector<FLOAT>&new_w){
    if(new_w.size()==this->dim){
        q=vector<c_float>(new_w.begin(), new_w.end());
        for (c_float &qi:q) {
            qi=-qi;
        }
        w_L2=new_w*new_w;
        osqp_update_lin_cost(work, q.data());
    }else{
        destroy();
        vector<FLOAT> h_ij(new_w.size());
        init_all(new_w, h_ij);
    }
}

template<typename FLOAT>
inline void qp_solver::update_w_h(const vector<FLOAT> &w, const vector<FLOAT>&h_ij){
    this->update_w(w);
    this->update_h(h_ij);
}

template<typename FLOAT>
inline FLOAT qp_solver::update_w_h_solve(const vector<FLOAT> &w, const vector<FLOAT> &h){
//    assert(w.size()==dim);
//    assert(h.size()==dim);
    this->update_w_h(w, h);
    c_float ret= this->qp_solve();
    return (FLOAT)ret;
}

template<typename FLOAT>
inline FLOAT qp_solver::update_w_h_solve(const vector<FLOAT> &w, const vector<FLOAT> &h, c_float *&solution){
//    assert(w.size()==dim);
//    assert(h.size()==dim);
    this->update_w_h(w, h);
    c_float ret= this->qp_solve(solution);
    return (FLOAT)ret;
}









