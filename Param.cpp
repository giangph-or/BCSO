#include "Param.h"

using namespace std;

Param::Param() {

}

Param::Param(Param* pParam) {
    this->hl = pParam->hl;
    this->gl = pParam->gl;
    this->gu = pParam->gu;
    this->gamma_h = pParam->gamma_h;
    this->gamma_g = pParam->gamma_g;
    this->lamda = pParam->lamda;
    this->max_1w = pParam->max_1w;
    this->K = pParam->K;
    this->tol_lamda = pParam->tol_lamda;
    this->M = pParam->M;
}

Param::Param(Data data, int K, double tol_lamda, int M) {
    this->K = K;
    this->tol_lamda = tol_lamda;
    this->M = M;
    hl = cpt_hl(data.T, data.m, data.eta, data.kappa, data.L);
    gl = cpt_gl(data.T, data.m, data.eta, data.kappa, data.L);
    gu = cpt_gu(data.T, data.m, data.eta, data.kappa, data.U);
    gamma_h = cpt_gamma_h(data.T, data.m, K, data.eta, data.kappa, data.L, data.U);
    gamma_g = cpt_gamma_g(data.T, data.m, K, data.eta, data.kappa, data.L, data.U);
    lamda = cpt_lamda(data.T, data.m, K, tol_lamda, data.a, data.b, hl, gl, gu, gamma_h, gamma_g);
    max_1w = cpt_max_1w(data.T, data.m, K, data.b, hl, gamma_h);
    Lw = cpt_Lw(data.T, data.m, data.eta, data.kappa, data.U);
    Uw = cpt_Uw(data.T, data.m, data.eta, data.kappa, data.L);
    Lwlog = cpt_Lwlog(data.T, data.m, data.eta, data.kappa, data.U);
    Uwlog = cpt_Uwlog(data.T, data.m, data.eta, data.kappa, data.L);
}

double Param::h_func(vector<vector<double>> eta, vector<vector<double>> kappa, int t, int i, double x) {
    //negative eta
    double res = exp(eta[t][i] * x + kappa[t][i]);
    return res;
}

double Param::g_func(vector<vector<double>> eta, vector<vector<double>> kappa, int t, int i, double x) {
    double res = x * exp(eta[t][i] * x + kappa[t][i]);
    return res;
}

vector<vector<double>> Param::cpt_hl(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp;

        for (int iter_m = 0; iter_m < m; iter_m++) {
            tmp.push_back(h_func(eta, kappa, iter_t, iter_m, L[iter_m]));

        }
        res.push_back(tmp);
    }
    return res;
}

vector<vector<double>> Param::cpt_gl(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp;

        for (int iter_m = 0; iter_m < m; iter_m++) {
            tmp.push_back(g_func(eta, kappa, iter_t, iter_m, L[iter_m]));

        }
        res.push_back(tmp);
    }
    return res;
}

vector<vector<double>> Param::cpt_gu(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> U) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp;

        for (int iter_m = 0; iter_m < m; iter_m++) {
            tmp.push_back(g_func(eta, kappa, iter_t, iter_m, U[iter_m]));

        }
        res.push_back(tmp);
    }
    return res;
}

vector<vector<vector<double>>> Param::cpt_gamma_h(int T, int m, int K, vector<vector<double>> eta, vector<vector<double>> kappa,
    vector<double> L, vector<double> U) {
    vector<vector<vector<double>>> res;

    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<vector<double>> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            vector<double> res_t_m;
            res_t_m.push_back(0);
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                double x1 = L[iter_m] + (U[iter_m] - L[iter_m]) * iter_k / K;
                double x2 = L[iter_m] + (U[iter_m] - L[iter_m]) * (iter_k - 1) / K;
                double res = h_func(eta, kappa, iter_t, iter_m, x1) - h_func(eta, kappa, iter_t, iter_m, x2);
                res = res * K / (U[iter_m] - L[iter_m]);
                res_t_m.push_back(res);
            }
            res_t.push_back(res_t_m);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<vector<vector<double>>> Param::cpt_gamma_g(int T, int m, int K, vector<vector<double>> eta, vector<vector<double>> kappa,
    vector<double> L, vector<double> U) {
    vector<vector<vector<double>>> res;

    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<vector<double>> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            vector<double> res_t_m;
            res_t_m.push_back(0);
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                double x1 = L[iter_m] + (U[iter_m] - L[iter_m]) * iter_k / K;
                double x2 = L[iter_m] + (U[iter_m] - L[iter_m]) * (iter_k - 1) / K;
                double res = g_func(eta, kappa, iter_t, iter_m, x1) - g_func(eta, kappa, iter_t, iter_m, x2);
                res = res * K / (U[iter_m] - L[iter_m]);
                res_t_m.push_back(res);
            }
            res_t.push_back(res_t_m);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<double> Param::cpt_lamda(int T, int m, int K, double tol, vector<double> a, vector<double> b, vector<vector<double>> hl,
    vector<vector<double>> gl, vector<vector<double>> gu, vector<vector<vector<double>>> gamma_h,
    vector<vector<vector<double>>> gamma_g) {
    vector<double> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        double max = gl[iter_t][0] / hl[iter_t][0];

        for (int iter_m = 0; iter_m < m; iter_m++) {
            if (max < gl[iter_t][iter_m] / hl[iter_t][iter_m])
                max = gl[iter_t][iter_m] / hl[iter_t][iter_m];
        }

        double sum_gu = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sum_gu += gu[iter_t][iter_m];
        }
        if (max < (a[iter_t] + sum_gu) / b[iter_t])
            max = (a[iter_t] + sum_gu) / b[iter_t];

        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                if (max < gamma_g[iter_t][iter_m][iter_k] / gamma_h[iter_t][iter_m][iter_k])
                    max = gamma_g[iter_t][iter_m][iter_k] / gamma_h[iter_t][iter_m][iter_k];
            }
        }
        max += tol;
        res.push_back(max);
    }
    return res;
}

vector<double> Param::cpt_max_1w(int T, int m, int K, vector<double> b, vector<vector<double>> hl,
    vector<vector<vector<double>>> gamma_h) {
    vector<double> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        double tmp = b[iter_t];
        for (int iter_m = 0; iter_m < m; iter_m++) {
            tmp += hl[iter_t][iter_m];
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                tmp += gamma_h[iter_t][iter_m][iter_k];
            }
        }
        res.push_back(tmp);
    }
    return res;
}

vector<vector<double>> Param::cpt_Lw(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> U) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp;

        for (int iter_m = 0; iter_m < m; iter_m++) {
            //TODO: it must be U[iter_m]
            tmp.push_back(h_func(eta, kappa, iter_t, iter_m, U[iter_m]));

        }
        res.push_back(tmp);
    }
    return res;
}

vector<vector<double>> Param::cpt_Uw(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp;

        for (int iter_m = 0; iter_m < m; iter_m++) {
            //TODO: it must be L[iter_m]
            tmp.push_back(h_func(eta, kappa, iter_t, iter_m, L[iter_m]));

        }
        res.push_back(tmp);
    }
    return res;
}

vector<vector<double>> Param::cpt_Lwlog(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> L) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp;

        for (int iter_m = 0; iter_m < m; iter_m++) {
            tmp.push_back(eta[iter_t][iter_m] * L[iter_m] + kappa[iter_t][iter_m]);
        }
        res.push_back(tmp);
    }
    return res;
}

vector<vector<double>> Param::cpt_Uwlog(int T, int m, vector<vector<double>> eta, vector<vector<double>> kappa, vector<double> U) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp;

        for (int iter_m = 0; iter_m < m; iter_m++) {
            tmp.push_back(eta[iter_t][iter_m] * U[iter_m] + kappa[iter_t][iter_m]);
        }
        res.push_back(tmp);
    }
    return res;
}