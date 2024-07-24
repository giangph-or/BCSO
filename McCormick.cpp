#include "McCormick.h"

McCormick::McCormick() {

}

McCormick::McCormick(Data data, Param param, bool exist_y) {
    lb_w_z_1 = cpt_lb_w_z_1(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h, param.max_1w);
    lb_w_z_0 = cpt_lb_w_z_0(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h, param.max_1w);
    ub_w_z_1 = cpt_ub_w_z_1(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h);
    ub_w_z_0 = cpt_ub_w_z_0(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h);

    if (exist_y) {
        lb_w_y_1 = cpt_lb_w_y_1(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h, param.max_1w);
        lb_w_y_0 = cpt_lb_w_y_0(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h, param.max_1w);
        ub_w_y_1 = cpt_ub_w_y_1(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h);
        ub_w_y_0 = cpt_ub_w_y_0(data.T, data.m, data.W, param.K, param.M, data.b, data.L, data.U, param.hl, param.gamma_h);
    }
}

bool McCormick::max_1w_y_1(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint use z instead of s
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint y_i = 1
        model.addConstr(y[i] == 1);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MAXIMIZE);
        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);

        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

bool McCormick::max_1w_y_0(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint use z instead of s
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint y_i = 1
        model.addConstr(y[i] == 0);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MAXIMIZE);

        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);
        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

bool McCormick::max_1w_z_1(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint use z instead of s
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint z_ik = 1
        model.addConstr(z[i][k] == 1);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MAXIMIZE);
        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);


        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

bool McCormick::max_1w_z_0(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint use z instead of s
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint z_ik = 1
        model.addConstr(z[i][k] == 0);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MAXIMIZE);
        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);

        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

vector<vector<double>> McCormick::cpt_lb_w_y_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
    vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h,
    vector<double> max_1w) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            double tmp_value = 0;
            bool lb_w_t_y_1 = max_1w_y_1(m, K, W, M, iter_m, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
            if (lb_w_t_y_1)
                tmp_value = 1 / tmp_value;
            else
                tmp_value = 1 / max_1w[iter_t];

            res_t.push_back(tmp_value);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<vector<double>> McCormick::cpt_lb_w_y_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
    vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h,
    vector<double> max_1w) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            double tmp_value = 0;
            bool lb_w_t_y_0 = max_1w_y_0(m, K, W, M, iter_m, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
            if (lb_w_t_y_0)
                tmp_value = 1 / tmp_value;
            else
                tmp_value = 1 / max_1w[iter_t];

            res_t.push_back(tmp_value);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<vector<vector<double>>> McCormick::cpt_lb_w_z_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
    vector<double> U, vector<vector<double>> hl,
    vector<vector<vector<double>>> gamma_h, vector<double> max_1w) {
    vector<vector<vector<double>>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<vector<double>> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            vector<double> res_t_m;
            res_t_m.push_back(0);
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                double tmp_value = 0;
                bool lb_w_t_s_1 = max_1w_z_1(m, K, W, M, iter_m, iter_k, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
                if (lb_w_t_s_1)
                    tmp_value = 1 / tmp_value;
                else
                    tmp_value = 1 / max_1w[iter_t];
                res_t_m.push_back(tmp_value);
            }

            res_t.push_back(res_t_m);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<vector<vector<double>>> McCormick::cpt_lb_w_z_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
    vector<double> U, vector<vector<double>> hl,
    vector<vector<vector<double>>> gamma_h, vector<double> max_1w) {
    vector<vector<vector<double>>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<vector<double>> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            vector<double> res_t_m;
            res_t_m.push_back(0);
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                double tmp_value = 0;
                bool lb_w_t_s_0 = max_1w_z_0(m, K, W, M, iter_m, iter_k, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
                if (lb_w_t_s_0)
                    tmp_value = 1 / tmp_value;
                else
                    tmp_value = 1 / max_1w[iter_t];
                res_t_m.push_back(tmp_value);
            }

            res_t.push_back(res_t_m);
        }
        res.push_back(res_t);
    }
    return res;
}

bool McCormick::min_1w_y_1(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint s = zy
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint y_i = 1
        model.addConstr(y[i] == 1);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);

        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

bool McCormick::min_1w_y_0(int m, int K, double W, int M, int i, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint s = zy
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint y_i = 1
        model.addConstr(y[i] == 0);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);

        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

bool McCormick::min_1w_z_1(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint s = zy
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint z[i][k] = 1
        model.addConstr(z[i][k] == 1);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MINIMIZE);

        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);

        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

bool McCormick::min_1w_z_0(int m, int K, double W, int M, int i, int k, double b_t, vector<double> L, vector<double> U, vector<double> hl_t,
    vector<vector<double>> gamma_h_t, double& res) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);

        // z variables
        vector<vector<GRBVar>> z(m, vector<GRBVar>(K + 1));
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int k = 0; k <= K; k++) {
                z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // y variables
        vector<GRBVar> y(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        // x^ variables
        vector<GRBVar> x(m);
        for (int iter_m = 0; iter_m < m; iter_m++) {
            x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        // Constraint z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            for (int iter_k = 1; iter_k <= K - 1; iter_k++) {
                model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
            }
        }

        // Constraint s = zy
        for (int iter_m = 0; iter_m < m; iter_m++) {
            model.addConstr(y[iter_m] >= z[iter_m][1]);
        }

        // Constraint x^ and z
        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr sumZi = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                sumZi += z[iter_m][iter_k];
            }
            sumZi = sumZi * (U[iter_m] - L[iter_m]) / K;
            sumZi += L[iter_m];
            model.addConstr(x[iter_m] == sumZi);
        }

        // Constraint sum_Y <= M
        GRBLinExpr sumY = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumY += y[iter_m];
        }
        model.addConstr(sumY <= M);

        // Constraint sum_X <= W
        GRBLinExpr sumX = 0;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            sumX += x[iter_m];
        }
        model.addConstr(sumX <= W);

        // Constraint z[i][k] = 1
        model.addConstr(z[i][k] == 0);

        // Objective function
        GRBLinExpr obj = b_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            obj += y[iter_m] * hl_t[iter_m];
        }

        for (int iter_m = 0; iter_m < m; iter_m++) {
            GRBLinExpr obj_gamma = 0;
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                obj_gamma += gamma_h_t[iter_m][iter_k] * z[iter_m][iter_k];
            }
            obj += obj_gamma * (U[iter_m] - L[iter_m]) / K;
        }

        model.setObjective(obj, GRB_MINIMIZE);

        model.set(GRB_DoubleParam_TimeLimit, 60.0);
        model.set(GRB_IntParam_OutputFlag, 0);

        bool check = false;
        model.optimize();

        check = true;
        res = model.get(GRB_DoubleAttr_ObjVal);

        return check;
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}

vector<vector<double>> McCormick::cpt_ub_w_y_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
    vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            double tmp_value = 0;
            bool ub_w_t_y_1 = min_1w_y_1(m, K, W, M, iter_m, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
            if (ub_w_t_y_1)
                tmp_value = 1 / tmp_value;
            else
                tmp_value = 1 / b[iter_t];

            res_t.push_back(tmp_value);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<vector<double>> McCormick::cpt_ub_w_y_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L, vector<double> U,
    vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h) {
    vector<vector<double>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            double tmp_value = 0;
            bool ub_w_t_y_0 = min_1w_y_0(m, K, W, M, iter_m, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
            if (ub_w_t_y_0)
                tmp_value = 1 / tmp_value;
            else
                tmp_value = 1 / b[iter_t];

            res_t.push_back(tmp_value);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<vector<vector<double>>> McCormick::cpt_ub_w_z_1(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
    vector<double> U, vector<vector<double>> hl,
    vector<vector<vector<double>>> gamma_h) {
    vector<vector<vector<double>>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<vector<double>> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            vector<double> res_t_m;
            res_t_m.push_back(0);
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                double tmp_value = 0;
                bool ub_w_t_z_1 = min_1w_z_1(m, K, W, M, iter_m, iter_k, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
                if (ub_w_t_z_1)
                    tmp_value = 1 / tmp_value;
                else
                    tmp_value = 1 / b[iter_t];
                res_t_m.push_back(tmp_value);
            }

            res_t.push_back(res_t_m);
        }
        res.push_back(res_t);
    }
    return res;
}

vector<vector<vector<double>>> McCormick::cpt_ub_w_z_0(int T, int m, double W, int K, int M, vector<double> b, vector<double> L,
    vector<double> U, vector<vector<double>> hl,
    vector<vector<vector<double>>> gamma_h) {
    vector<vector<vector<double>>> res;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<vector<double>> res_t;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            vector<double> res_t_m;
            res_t_m.push_back(0);
            for (int iter_k = 1; iter_k <= K; iter_k++) {
                double tmp_value = 0;
                bool ub_w_t_z_0 = min_1w_z_0(m, K, W, M, iter_m, iter_k, b[iter_t], L, U, hl[iter_t], gamma_h[iter_t], tmp_value);
                if (ub_w_t_z_0)
                    tmp_value = 1 / tmp_value;
                else
                    tmp_value = 1 / b[iter_t];
                res_t_m.push_back(tmp_value);
            }

            res_t.push_back(res_t_m);
        }
        res.push_back(res_t);
    }
    return res;
}