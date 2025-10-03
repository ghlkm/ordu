//
// Created by Keming Li on 2023/1/12.
//

#ifndef IPREF_TEST_H
#define IPREF_TEST_H
#include "qhull_user.h"
#include <vector>
double tri_area(double a, double b, double c){
    double p=(a+b+c)/2.0;
    return sqrt(p * (p-a) * (p-b) * (p-c));
}

void test_volume(){
    // in 3 dimension
//    double volumn_at_half_inter(const realT *pointCoordinates, int p_num, vector<float> &innerPoint);
    uint32_t dim=3;
    uint32_t halfspace_num=3+1;
    // three points that intersect with user preference:
    // A(0.3, 0.6, 0.1), B(0.2, 0.3, 0.5), C(0.4, 0.4, 0.2)
    // AB^2=0.26, AC^2=0.06, BC^2=0.14
    // area \appraximateEqual 0.0433
    // volume \apparaximateEqual 1/3 * \sqrt{3}
    // n1 * A =0, n1 * B =0
    // n1=\lambda_1 * (-27/13, 1, 3/13)
    // n2 * A =0, n2 * C =0
    // n2=\lambda_2 * (-4, 1, 6)
    // n3 * B =0, n3 * C =0
    // n3=\lambda_3 * (7/8, 1, -1/4)
    // an obvious interior point
    // D(0.3/2, 1.3/3/2, 0.8/3/2)
    // n1 * D < 0, n2 * D >0, n3*D>0
    // thereby, the halfspace expression of the three points and (0, 0, 0) will be:
    // {w| n1 * w <=0, -n2 * w<=0, -n3 * w <=0, \vec{1} * w<=-(-1)}
    realT *pointCoordinates=new realT[(dim+1)*halfspace_num];
    pointCoordinates[0]=-27.0/13;
    pointCoordinates[1]=1;
    pointCoordinates[2]=3.0/13;
    pointCoordinates[3]=0;
    pointCoordinates[4]=4;
    pointCoordinates[5]=-1;
    pointCoordinates[6]=-6;
    pointCoordinates[7]=0;
    pointCoordinates[8]=-7.0/8;
    pointCoordinates[9]=-1;
    pointCoordinates[10]=1.0/4;
    pointCoordinates[11]=0;
    pointCoordinates[12]=1;
    pointCoordinates[13]=1;
    pointCoordinates[14]=1;
    pointCoordinates[15]=-1;
    std::vector<float> innerPoint{0.3/2, 1.3/3/2, 0.8/3/2};
    qhull_user qu;

//    qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
//    qhT* qh = &qh_qh;
//    QHULL_LIB_CHECK
//    qh_zero(qh, NULL);
//    fflush(NULL);

    cout<<qu.volumn_at_half_inter(pointCoordinates, halfspace_num, innerPoint)<<endl;
    delete [] (pointCoordinates);
}


void test_halfspace(){
    // in 3 dimension
//    double volumn_at_half_inter(const realT *pointCoordinates, int p_num, vector<float> &innerPoint);
    uint32_t dim=3;
    uint32_t pnum_num=3+1;
    // three points that intersect with user preference:
    // A(0.3, 0.6, 0.1), B(0.2, 0.3, 0.5), C(0.4, 0.4, 0.2)
    // AB^2=0.26, AC^2=0.06, BC^2=0.14
    // area \appraximateEqual 0.0433
    // volume \apparaximateEqual 1/3 * \sqrt{3}
    // n1 * A =0, n1 * B =0
    // n1=\lambda_1 * (-27/13, 1, 3/13)
    // n2 * A =0, n2 * C =0
    // n2=\lambda_2 * (-4, 1, 6)
    // n3 * B =0, n3 * C =0
    // n3=\lambda_3 * (7/8, 1, -1/4)
    // an obvious interior point
    // D(0.3/2, 1.3/3/2, 0.8/3/2)
    // n1 * D > 0, n2 * D >0, n3*D>0
    // thereby, the halfspace expression of the three points and (0, 0, 0) will be:
    // {w| -n1 * w <=0, -n2 * w<=0, -n3 * w <=0, \vec{1} * w<=-(-1)}
    realT *pointCoordinates=new realT[(dim)*pnum_num];
    pointCoordinates[0]=0;
    pointCoordinates[1]=0;
    pointCoordinates[2]=0;
    pointCoordinates[3]=.3;
    pointCoordinates[4]=.6;
    pointCoordinates[5]=.1;
    pointCoordinates[6]=.2;
    pointCoordinates[7]=.3;
    pointCoordinates[8]=.5;
    pointCoordinates[9]=.4;
    pointCoordinates[10]=.4;
    pointCoordinates[11]=.2;
//    pointCoordinates[12]=1;
//    pointCoordinates[13]=1;
//    pointCoordinates[14]=1;
//    pointCoordinates[15]=-1;
    std::vector<float> innerPoint{0.3/2, 1.3/3/2, 0.8/3/2};
    qhull_user qu;
    Qhull q;
    std::stringstream output;
    q.runQhull("normals of square", dim, pnum_num, pointCoordinates, "QJ"); // halfspace intersect


    output.clear();
    q.outputQhull("n");
    cout<<"debug: after trying to get volume"<<endl;
    cout<<output.str()<<endl;
//    qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
//    qhT* qh = &qh_qh;
//    QHULL_LIB_CHECK
//    qh_zero(qh, NULL);
//    fflush(NULL);

    delete [] (pointCoordinates);
}


void test_nonRedundantHalfspace(){
    uint32_t p_num=6;
    std::vector<float> innerPoint{0.5, 0.5, 0.5};
    int dim=innerPoint.size();
    realT *pointCoordinates=new realT[(dim+1)*p_num];
    pointCoordinates[0]= -1;
    pointCoordinates[1]=0;
    pointCoordinates[2]=0;
    pointCoordinates[3]=0;

    pointCoordinates[4]=0;
    pointCoordinates[5]=-1;
    pointCoordinates[6]=0;
    pointCoordinates[7]=0;

    pointCoordinates[8]=0;
    pointCoordinates[9]=0;
    pointCoordinates[10]=-1;
    pointCoordinates[11]=0;

    pointCoordinates[12]=1;
    pointCoordinates[13]=0;
    pointCoordinates[14]=0;
    pointCoordinates[15]= -1;

    pointCoordinates[16]=0;
    pointCoordinates[17]=1;
    pointCoordinates[18]=0;
    pointCoordinates[19]=-1;

    pointCoordinates[20]=0;
    pointCoordinates[21]=0;
    pointCoordinates[22]=1;
    pointCoordinates[23]=-1;

    Qhull q;
    orgQhull::Coordinates feasible;
    for (float i:innerPoint) {
        feasible << i;
    }
    q.setFeasiblePoint(feasible);
    try {
        std::stringstream output;
        q.setOutputStream(&output);
        q.runQhull("normals of square", dim+1, p_num, pointCoordinates, "H Pp"); // halfspace intersect
        // see http://www.qhull.org/html/qh-optf.htm#FS

        output.clear();
        q.outputQhull("Fp");

        int out_dim, out_vtxCnt;
        output>>out_dim;
        output>>out_vtxCnt;
        vector<vector<double>> ret(out_vtxCnt, vector<double>(out_dim));
        for (int i = 0; i < out_vtxCnt; ++i) {
            for (int j = 0; j < out_dim; ++j) {
                output>>ret[i][j];
            }
        }
        output.clear();

        cout<<"----------------"<<endl;
        q.outputQhull("i");
        cout<<output.str()<<endl;

    } catch (std::exception &e) {// catch by ref
        cout<<e.what();
        vector<vector<double>> UNUSED_VVD;
    }
}

#include "osqp.h"

void test_osqp2(){
    int dim=2;
    float **P=new float*[3];
    P[0]=new float[2];
    P[0][0]=0.6; P[0][1]=0.6;
    P[1]=new float[2];
    P[1][0]=0.3; P[1][1]=0.7;
    P[2]=new float[2];
    P[2][0]=0.7; P[2][1]=0.3;
    vector<int> gt(1, 0);
    vector<int> le={1, 2};
    vector<double> ret;
    qp_solver q(P,gt, le, 2, ret);
    cout<<ret<<endl;
}

int test_osqp() {
    // Load problem data
    c_float P_x[3] = {};
    c_int P_nnz = 0;
    c_int P_i[3] = { };
    c_int P_p[3] = {};
    c_float q[2] = {1.0, 1.0, };
    c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };
    c_int A_nnz = 4;
    c_int A_i[4] = {0, 1, 0, 2, };
    c_int A_p[3] = {0, 2, 4, };
    c_float l[3] = {1.0, 0.0, 0.0, };
    c_float u[3] = {1.0, 0.7, 0.7, };
    c_int n = 2;
    c_int m = 3;

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define solver settings as default
    if (settings) {
        osqp_set_default_settings(settings);
        settings->alpha = 1.0; // Change alpha parameter
    }

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve Problem
    osqp_solve(work);

    // Cleanup
    osqp_cleanup(work);
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings) c_free(settings);

    return exitflag;
};


#endif //IPREF_TEST_H
