

#include "qp_solver2.h"
#include "cassert"
#include <limits>
double qp_solver2(const vector<float>& w, const vector<vector<double>>& H){
//        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * p
//        ce0: p
//
//        CI: n * m
//        ci0: m
//
//        x: n
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=w.size(), m, p;
    G.resize(n, n);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if(i==j){
                    G[i][j]=1.0;
                }else{
                    G[i][j]=0.0;
                }
            }
        }
    }
    g0.resize(n);
    {
        for (int i = 0; i < n; i++){
            g0[i]=-w[i];
        }
    }
    m = 1;
    CE.resize(n, m);
    {
        for (int i = 0; i < n; i++)
            CE[i][0]=1.0;
    }
    ce0.resize(m);
    {
        for (int j = 0; j < m; j++)
            ce0[j]=-1.0;
    }
    p = H.size()+w.size();
    CI.resize(n, p);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < H.size(); j++){
                CI[i][j]=-H[j][i];
            }
            for (int j = H.size(); j < p; j++){
                if(i==j-H.size()){
                    CI[i][j]=1.0;
                }else{
                    CI[i][j]=0.0;
                }
            }
        }
    }
    ci0.resize(p);
    {
        for (int j = 0; j < p; j++)
            ci0[j]=0;
    }
    x.resize(n);
    double solverret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    double ret=0;
    for (int k = 0; k < n; ++k) {
        ret+=(x[k]-w[k])*(x[k]-w[k]);
    }
    if(solverret==std::numeric_limits<double>::infinity()){
        return INFINITY;
    }
    return sqrt(ret);
}


double qp_solver2(const vector<float>& w, const vector<vector<double>>& H1, const vector<vector<double>>& H2){
//        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * p
//        ce0: p
//
//        CI: n * m
//        ci0: m
//
//        x: n
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=w.size(), m, p;
    G.resize(n, n);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if(i==j){
                    G[i][j]=1.0;
                }else{
                    G[i][j]=0.0;
                }
            }
        }
    }
    g0.resize(n);
    {
        for (int i = 0; i < n; i++){
            g0[i]=-w[i];
        }
    }
    m = 1;
    CE.resize(n, m);
    {
        for (int i = 0; i < n; i++)
            CE[i][0]=1.0;
    }
    ce0.resize(m);
    {
        for (int j = 0; j < m; j++)
            ce0[j]=-1.0;
    }
    p = H1.size()+H2.size()+w.size();
    CI.resize(n, p);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < H1.size(); j++){
                CI[i][j]=-H1[j][i];
            }
            for (int j = H1.size(); j < H1.size()+H2.size(); j++){
                CI[i][j]=-H2[j-H1.size()][i];
            }
            for (int j = H1.size()+H2.size(); j < p; j++){
                if(i==j-H1.size()-H2.size()){
                    CI[i][j]=1.0;
                }else{
                    CI[i][j]=0.0;
                }
            }
        }
    }
    ci0.resize(p);
    {
        for (int j = 0; j < p; j++)
            ci0[j]=0;
    }
    x.resize(n);
    double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
//    assert(solver_ret!=std::numeric_limits<double>::infinity());
    double ret=0;
    for (int k = 0; k < n; ++k) {
        ret+=(x[k]-w[k])*(x[k]-w[k]);
    }
    if(solver_ret == std::numeric_limits<double>::infinity()){
        return INFINITY;
    }
    return sqrt(ret);
}

double qp_solver2(const vector<float>& w, const vector<vector<double>>& H1, int opt, const vector<int>& cmp, float **PG){
//        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * m
//        ce0: m
//
//        CI: n * p
//        ci0: p
//
//        x: n
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=w.size(), m, p;
    G.resize(n, n);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if(i==j){
                    G[i][j]=1.0;
                }else{
                    G[i][j]=0.0;
                }
            }
        }
    }
    g0.resize(n);
    {
        for (int i = 0; i < n; i++){
            g0[i]=-w[i];
        }
    }
    m = 1;
    CE.resize(n, m);
    {
        for (int i = 0; i < n; i++)
            CE[i][0]=1.0;
    }
    ce0.resize(m);
    {
        for (int j = 0; j < m; j++)
            ce0[j]=-1.0;
    }
    p = H1.size()+cmp.size()+w.size();
    CI.resize(n, p);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < H1.size(); j++){
                CI[i][j]=-H1[j][i];
            }
            for (int j = H1.size(); j < H1.size()+cmp.size(); j++){
                CI[i][j]=PG[opt][i]-PG[cmp[j-H1.size()]][i];
            }
            for (int j = H1.size()+cmp.size(); j < p; j++){
                if(i==j-H1.size()-cmp.size()){
                    CI[i][j]=1.0;
                }else{
                    CI[i][j]=0.0;
                }
            }
        }
    }
    ci0.resize(p);
    {
        for (int j = 0; j < p; j++)
            ci0[j]=0;
    }
    x.resize(n);
    double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    double ret=0;
    for (int k = 0; k < n; ++k) {
        ret+=(x[k]-w[k])*(x[k]-w[k]);
    }
    if(solver_ret == std::numeric_limits<double>::infinity()){
        return INFINITY;
    }
    return sqrt(ret);
}


double qp_solver2(const vector<float>& w, const vector<vector<double>>& H1, int opt, const vector<int>& cmp, float **PG, vector<double> &retv){
//        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * m
//        ce0: m
//
//        CI: n * p
//        ci0: p
//
//        x: n
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=w.size(), m, p;
    G.resize(n, n);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if(i==j){
                    G[i][j]=1.0;
                }else{
                    G[i][j]=0.0;
                }
            }
        }
    }
    g0.resize(n);
    {
        for (int i = 0; i < n; i++){
            g0[i]=-w[i];
        }
    }
    m = 1;
    CE.resize(n, m);
    {
        for (int i = 0; i < n; i++)
            CE[i][0]=1.0;
    }
    ce0.resize(m);
    {
        for (int j = 0; j < m; j++)
            ce0[j]=-1.0;
    }
    p = H1.size()+cmp.size()+w.size();
    CI.resize(n, p);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < H1.size(); j++){
                CI[i][j]=-H1[j][i];
            }
            for (int j = H1.size(); j < H1.size()+cmp.size(); j++){
                CI[i][j]=PG[opt][i]-PG[cmp[j-H1.size()]][i];
            }
            for (int j = H1.size()+cmp.size(); j < p; j++){
                if(i==j-H1.size()-cmp.size()){
                    CI[i][j]=1.0;
                }else{
                    CI[i][j]=0.0;
                }
            }
        }
    }
    ci0.resize(p);
    {
        for (int j = 0; j < p; j++)
            ci0[j]=0;
    }
    x.resize(n);
    double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    double ret=0;
    for (int k = 0; k < n; ++k) {
        ret+=(x[k]-w[k])*(x[k]-w[k]);
    }
    retv.resize(n);
    for (int k = 0; k < n; ++k) {
        retv[k]=x[k];
    }
    if(solver_ret == std::numeric_limits<double>::infinity()){
        return INFINITY;
    }
    return sqrt(ret);
}


vector<float> find_point_in_region(const vector<float>& w, const vector<vector<double>>& H){
//        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * p
//        ce0: p
//
//        CI: n * m
//        ci0: m
//
//        x: n
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=w.size(), m, p;
    G.resize(n, n);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if(i==j){
                    G[i][j]=1.0;
                }else{
                    G[i][j]=0.0;
                }
            }
        }
    }
    g0.resize(n);
    {
        for (int i = 0; i < n; i++){
            g0[i]=-w[i];
        }
    }
    m = 1;
    CE.resize(n, m);
    {
        for (int i = 0; i < n; i++)
            CE[i][0]=1.0;
    }
    ce0.resize(m);
    {
        for (int j = 0; j < m; j++)
            ce0[j]=-1.0;
    }
    p = H.size()+w.size();
    CI.resize(n, p);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < H.size(); j++){
                CI[i][j]=-H[j][i];
            }
            for (int j = H.size(); j < p; j++){
                if(i==j-H.size()){
                    CI[i][j]=1.0;
                }else{
                    CI[i][j]=0.0;
                }
            }
        }
    }
    ci0.resize(p);
    {
        for (int j = 0; j < p; j++)
            ci0[j]=0;
    }
    x.resize(n);
    double solverret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
//    assert(solverret!=std::numeric_limits<double>::infinity());
    vector<float> ret;
    if(solverret == std::numeric_limits<double>::infinity()){
        return ret;
    }
    for (int k = 0; k <n ; ++k) {
        ret.push_back(x[k]);
    }
    return ret;
}


vector<double> find_point_in_region(const vector<float>& w, const vector<vector<double>>& H1, const vector<vector<double>>& H2){
//        min 0.5 * x G x + g0 x
//        s.t.
//                CE^T x + ce0 = 0
//        CI^T x + ci0 >= 0
//        G: n * n
//        g0: n
//
//        CE: n * p
//        ce0: p
//
//        CI: n * m
//        ci0: m
//
//        x: n
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    int n=w.size(), m, p;
    G.resize(n, n);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if(i==j){
                    G[i][j]=1.0;
                }else{
                    G[i][j]=0.0;
                }
            }
        }
    }
    g0.resize(n);
    {
        for (int i = 0; i < n; i++){
            g0[i]=-w[i];
        }
    }
    m = 1;
    CE.resize(n, m);
    {
        for (int i = 0; i < n; i++)
            CE[i][0]=1.0;
    }
    ce0.resize(m);
    {
        for (int j = 0; j < m; j++)
            ce0[j]=-1.0;
    }
    p = H1.size()+H2.size()+w.size();
    CI.resize(n, p);
    {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < H1.size(); j++){
                CI[i][j]=-H1[j][i];
            }
            for (int j = H1.size(); j < H1.size()+H2.size(); j++){
                CI[i][j]=-H2[j-H1.size()][i];
            }
            for (int j = H1.size()+H2.size(); j < p; j++){
                if(i==j-H1.size()-H2.size()){
                    CI[i][j]=1.0;
                }else{
                    CI[i][j]=0.0;
                }
            }
        }
    }
    ci0.resize(p);
    {
        for (int j = 0; j < p; j++)
            ci0[j]=0;
    }
    x.resize(n);
    double solver_ret=solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    vector<double> ret;
    if(solver_ret == std::numeric_limits<double>::infinity()){
        return ret;
    }
    for (int k = 0; k <n ; ++k) {
        ret.push_back(x[k]);
    }
    return ret;
}