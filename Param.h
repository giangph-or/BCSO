#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <stdio.h>
#include "Data.h"


using namespace std;
class Param {
public:
    vector<vector<double>> hl;
    vector<vector<double>> gl;
    vector<vector<double>> gu;
    vector<vector<vector<double>>> gamma_h;
    vector<vector<vector<double>>> gamma_g;
    vector<double> lamda;
    vector<double> max_1w;
    int K;
    double tol_lamda;
    int M;
    vector<vector<double>> Lw;
    vector<vector<double>> Uw;
    vector<vector<double>> Lwlog;
    vector<vector<double>> Uwlog;

    Param();

    Param(Param* pParam);

    Param(Data data, int K, double tol_lamda, int M);

    double h_func(vector<vector<double>> eta, vector<vector<double>> kappa, int t, int i, double x);

    double g_func(vector<vector<double>> eta, vector<vector<double>> kappa, int t, int i, double x);

    vector<vector<double>> cpt_hl(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L);

    vector<vector<double>> cpt_gl(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L);

    vector<vector<double>> cpt_gu(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> U);

    vector<vector<vector<double>>> cpt_gamma_h(int T, int m, int K, vector<vector<double>> eta, vector<vector<double>> kappa,
        vector<double> L, vector<double> U);

    vector<vector<vector<double>>> cpt_gamma_g(int T, int m, int K, vector<vector<double>> eta, vector<vector<double>> kappa,
        vector<double> L, vector<double> U);

    vector<double> cpt_lamda(int T, int m, int K, double tol, vector<double> a, vector<double> b, vector<vector<double>> hl,
        vector<vector<double>> gl, vector<vector<double>> gu, vector<vector<vector<double>>> gamma_h,
        vector<vector<vector<double>>> gamma_g);

    vector<double> cpt_max_1w(int T, int m, int K, vector<double> b, vector<vector<double>> hl,
        vector<vector<vector<double>>> gamma_h);

    vector<vector<double>> cpt_Lw(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L);

    vector<vector<double>> cpt_Uw(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> U);

    vector<vector<double>> cpt_Lwlog(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L);

    vector<vector<double>> cpt_Uwlog(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> U);
};
#pragma once