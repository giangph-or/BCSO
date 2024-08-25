#include <iostream>
#include <vector>
#include <gurobi_c++.h>
#include "Data.h"
#include "Param.h"
#include "McCormick.h"
#include <string>
#include <fstream>
#include <chrono>

using namespace std;

class CBMISOCP : public GRBCallback {
public:
    int cb_T;
    int cb_m;
    int cb_K;
    vector<vector<double>> cb_hl;
    vector<vector<vector<double>>> cb_gamma_h;
    vector<double> cb_b;
    vector<double> cb_L;
    vector<double> cb_U;
    GRBVar* w;
    GRBVar* y;
    GRBVar** z;
    CBMISOCP(GRBVar* w, GRBVar* y, GRBVar** z, int T, int m, int K, vector<double> b, vector<double> L, vector<double> U, vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h);
protected:
    void callback();
};

class BCMISOCPSolver {
public:
    Data data;
    Param param;
    McCormick mcCormick;
    double time_for_param = 0;
    double time_for_mc = 0;
    double time_for_solve = 0;
    double gap = 0;
    double obj_val_cplex = 0;
    double obj_val_true = 0;

    BCMISOCPSolver();

    BCMISOCPSolver(Data data, int K, double tol_lamda, int M);

    void solve(string output, int time_limit);
};