#include "ApproxByGurobi.h"

ApproxGurobiSolver::ApproxGurobiSolver() {

}

ApproxGurobiSolver::ApproxGurobiSolver(Data data, int K, double tol_lamda, int M) {
    this->data = data;

    auto start = chrono::steady_clock::now(); //get start time
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;

    param = Param(data, K, tol_lamda, M);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    time_for_param = elapsed_seconds.count();

    //start = chrono::steady_clock::now(); //get start time
    //mcCormick = McCormick(data, param, true);
    //end = std::chrono::steady_clock::now();
    //elapsed_seconds = end - start;
    //time_for_mc = elapsed_seconds.count();
}

void ApproxGurobiSolver::solve(string output, int time_limit) {
    GRBEnv env = GRBEnv();

    // Create a model
    GRBModel model = GRBModel(env);

    // Create variables theta, v, u
    vector<GRBVar> theta(data.T);
    vector<GRBVar> v(data.T);
    vector<GRBVar> u(data.T);
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        theta[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "theta_" + to_string(iter_t));
        v[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "v_" + to_string(iter_t));
        u[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "u_" + to_string(iter_t));
    }

    // Create variables x,y
    vector<GRBVar> y(data.m);
    vector<GRBVar> x(data.m);
    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        x[iter_m] = model.addVar(data.L[iter_m], data.U[iter_m], 0.0, GRB_CONTINUOUS, "x_" + to_string(iter_m));
        y[iter_m] = model.addVar(0, 1, 0.0, GRB_BINARY, "y_" + to_string(iter_m));
    }

    // Variable s, w
    auto** w = new GRBVar * [data.T];
    auto** w_log = new GRBVar * [data.T];
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        w[iter_t] = reinterpret_cast<GRBVar*>(new GRBVar * [data.m]);
        w_log[iter_t] = reinterpret_cast<GRBVar*>(new GRBVar * [data.m]);
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            w[iter_t][iter_m] = model.addVar(param.Lw[iter_t][iter_m], param.Uw[iter_t][iter_m], 0.0, GRB_CONTINUOUS, "w_" + to_string(iter_t) + "_" + to_string(iter_m));
            w_log[iter_t][iter_m] = model.addVar(param.Lwlog[iter_t][iter_m], param.Uwlog[iter_t][iter_m], 0.0, GRB_CONTINUOUS, "w_log_" + to_string(iter_t) + "_" + to_string(iter_m));
        }
    }

    vector<vector<GRBVar>> s(data.T);
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            s[iter_t].push_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
    }

    // Constraint theta, v, and u
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        model.addQConstr(theta[iter_t] * v[iter_t] <= u[iter_t]);
    }

    // PWL for w
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            model.addConstr(w_log[iter_t][iter_m] == x[iter_m] * data.eta[iter_t][iter_m] + data.kappa[iter_t][iter_m]);
            model.addGenConstrExp(w_log[iter_t][iter_m], w[iter_t][iter_m]);
        }
    }

    // constr v and s
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        GRBLinExpr sum_s = 1;
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            sum_s += s[iter_t][iter_m];
        }
        model.addConstr(v[iter_t] == sum_s);
    }

    // constr u and x, s
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        GRBQuadExpr sum_xs;
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            sum_xs += x[iter_m] * s[iter_t][iter_m];
        }
        model.addQConstr(sum_xs == u[iter_t]);
    }

    // constr y and x
    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        model.addConstr(x[iter_m] <= data.U[iter_m] * y[iter_m] + data.L[iter_m]);
    }

    // constr y and s
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            model.addConstr(param.Lw[iter_t][iter_m] * y[iter_m] <= s[iter_t][iter_m]);
            model.addConstr(param.Uw[iter_t][iter_m] * y[iter_m] >= s[iter_t][iter_m]);
            model.addConstr(w[iter_t][iter_m] - param.Lw[iter_t][iter_m] * (1 - y[iter_m]) >= s[iter_t][iter_m]);
            model.addConstr(w[iter_t][iter_m] - param.Uw[iter_t][iter_m] * (1 - y[iter_m]) <= s[iter_t][iter_m]);
        }
    }

    // Constraint Ax + By <= D : sum_Y <= M
    GRBQuadExpr sum_y;
    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        sum_y += y[iter_m];
    }
    model.addConstr(sum_y <= param.M);

    //Constraint Ax + By <= D : sum_X <= W
    GRBQuadExpr sum_x = 0;
    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        sum_x += x[iter_m]* y[iter_m];
    }
    model.addQConstr(sum_x <= data.W);

    GRBLinExpr objective;
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        objective += theta[iter_t];
    }

    model.setObjective(objective, GRB_MAXIMIZE);
    //model.set(GRB_IntParam_Threads, 8);
    model.set(GRB_DoubleParam_TimeLimit, time_limit);
    model.set(GRB_DoubleParam_MIPGap, 1e-8);
    model.set(GRB_IntParam_MIQCPMethod, 1);
    model.set(GRB_IntParam_FuncPieces, 25);
    //model.set(GRB_IntParam_FuncPieces, -1);
    //model.set(GRB_DoubleParam_FuncPieceError, 1e-5);

    auto start = chrono::steady_clock::now(); //get start time
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;

    model.optimize();

    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    time_for_solve = elapsed_seconds.count();
    cout << "Status: " << model.get(GRB_IntAttr_Status) << endl;
    cout << "Objective value: " << setprecision(8) << model.get(GRB_DoubleAttr_ObjVal) << endl;

    double obj_val_gurobi = model.get(GRB_DoubleAttr_ObjVal);

    //get value x^
    vector<double> ansX(data.m, 0);
    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        ansX[iter_m] = x[iter_m].get(GRB_DoubleAttr_X);
    }

    //get value y
    vector<bool> ansY(data.m);
    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        ansY[iter_m] = y[iter_m].get(GRB_DoubleAttr_X) > 1 - model.get(GRB_DoubleParam_IntFeasTol);
    }

    gap = model.get(GRB_DoubleAttr_MIPGap);

    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        cout << "y[" << iter_m << "] = " << ansY[iter_m] << endl;
    }

    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        cout << "x[" << iter_m << "] = " << setprecision(5) << ansX[iter_m] << endl;
    }

    //Compute true obj val
    obj_val_true = 0;
    for (int iter_t = 0; iter_t < data.T; iter_t++) {
        double tmp1 = data.a[iter_t];
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            tmp1 += ansY[iter_m] * param.g_func(data.eta, data.kappa, iter_t, iter_m, ansX[iter_m]);
        }
        double tmp2 = data.b[iter_t];
        for (int iter_m = 0; iter_m < data.m; iter_m++) {
            tmp2 += ansY[iter_m] * param.h_func(data.eta, data.kappa, iter_t, iter_m, ansX[iter_m]);
        }
        obj_val_true += tmp1 / tmp2;
    }

    cout << "Solving, it took " << elapsed_seconds.count() << " seconds" << endl;
    ofstream report_results(output, ofstream::out);
    report_results << data.m << "," << data.W << "," << param.M
        << setprecision(8) << fixed << "," << obj_val_gurobi
        << "," << obj_val_true << "," << time_for_solve << "," << gap << endl;
    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        report_results << "y[" << iter_m << "] = " << ansY[iter_m] << endl;
    }

    for (int iter_m = 0; iter_m < data.m; iter_m++) {
        report_results << "x[" << iter_m << "] = " << setprecision(5) << ansX[iter_m] << endl;
    }
    report_results.close();
}


