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

class ExpConeAssortment {
public:
    DataAssortment data;
    ParamAssortment param;
    double time_for_param = 0;
    double time_for_mc = 0;
    double time_for_solve = 0;
    double gap = 0;
    double obj_val_gurobi = 0;
    double obj_val_true = 0;


    ExpConeAssortment();

    ExpConeAssortment(DataAssortment data, int K, double tol_lamda, int M);

    void solve(string output, int time_limit);

    void solve_BC(string output, int time_limit);
};

class CBExpAssortment : public GRBCallback {
public:
    int T_;
    int m_;
    int K_;
    vector<double> b_;
    vector<double> U_;
    vector<double> L_;
    vector<vector<double>> hl_;
    vector<vector<vector<double>>> gamma_h_;
    GRBVar* y;
    GRBVar* n;
    GRBVar* d;
    GRBVar* f;
    GRBVar** z;
    CBExpAssortment(GRBVar* y, GRBVar* n, GRBVar* d, GRBVar* f, GRBVar** z,
        int T_, int m_, int K_, vector<double> b_, vector<double> U_, vector<double> L_,
        vector<vector<double>> hl_, vector<vector<vector<double>>> gamma_h_);
protected:
    void callback();
};

class ExpConeFacility {
public:
    DataFacility data;
    ParamFacility param;
    double time_for_param = 0;
    double time_for_mc = 0;
    double time_for_solve = 0;
    double gap = 0;
    double obj_val_gurobi = 0;
    double obj_val_true = 0;


    ExpConeFacility();

    ExpConeFacility(DataFacility data, int K, double tol_lamda, int M);

    void solve(string output, int time_limit);
};
