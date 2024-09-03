#include <vector>
#include <iostream>
#include <gurobi_c++.h>
#include "Data.h"
#include "Param.h"

using namespace std;

class McCormick {
public:
    vector<vector<double>> lb_w_y_1;
    vector<vector<double>> lb_w_y_0;
    vector<vector<vector<double>>> lb_w_z_1;
    vector<vector<vector<double>>> lb_w_z_0;

    vector<vector<double>> ub_w_y_1;
    vector<vector<double>> ub_w_y_0;
    vector<vector<vector<double>>> ub_w_z_1;
    vector<vector<vector<double>>> ub_w_z_0;

    McCormick();

    McCormick(DataAssortment data, ParamAssortment param, bool exist_y);

    McCormick(DataFacility data, ParamFacility param, bool exist_y);

    bool max_1w_y_1(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    bool max_1w_y_0(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    bool max_1w_z_1(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    bool max_1w_z_0(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    vector<vector<double>> cpt_lb_w_y_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
        vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h,
        vector<double> max_1w);

    vector<vector<double>> cpt_lb_w_y_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
        vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h,
        vector<double> max_1w);

    vector<vector<vector<double>>> cpt_lb_w_z_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
        vector<double> U, vector<vector<double>> hl,
        vector<vector<vector<double>>> gamma_h, vector<double> max_1w);

    vector<vector<vector<double>>> cpt_lb_w_z_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
        vector<double> U, vector<vector<double>> hl,
        vector<vector<vector<double>>> gamma_h, vector<double> max_1w);

    bool min_1w_y_1(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    bool min_1w_y_0(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    bool min_1w_z_1(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    bool min_1w_z_0(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
        vector<vector<double>> gamma_h_t, double& res);

    vector<vector<double>> cpt_ub_w_y_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
        vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h);

    vector<vector<double>> cpt_ub_w_y_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
        vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h);

    vector<vector<vector<double>>> cpt_ub_w_z_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
        vector<double> U, vector<vector<double>> hl,
        vector<vector<vector<double>>> gamma_h);

    vector<vector<vector<double>>> cpt_ub_w_z_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
        vector<double> U, vector<vector<double>> hl,
        vector<vector<vector<double>>> gamma_h);
};

#pragma once
