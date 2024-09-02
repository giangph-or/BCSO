#include "BilinearSolver.h"

using namespace std;

BilinearSolverAssortment::BilinearSolverAssortment() {

}

BilinearSolverAssortment::BilinearSolverAssortment(DataAssortment data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamAssortment(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();

	start = chrono::steady_clock::now(); //get start time
	mcCormick = McCormickAssortment(data, param, true);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_mc = elapsed_seconds.count();
}

void BilinearSolverAssortment::solve(string output, int time_limit) {
	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		// w variables
		vector<GRBVar> w(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			w[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// z variables
		vector<vector<GRBVar>> z(data.m, vector<GRBVar>(param.K + 1));
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			for (int k = 0; k <= param.K; k++) {
				z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}

		// y variables
		vector<GRBVar> y(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		// x^ variables
		vector<GRBVar> x(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// Constraint w
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			GRBQuadExpr sum = data.b[iter_t] * w[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				sum += param.hl[iter_t][iter_m] * w[iter_t] * y[iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBQuadExpr sum_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					sum_gamma += param.gamma_h[iter_t][iter_m][iter_k] * w[iter_t] * z[iter_m][iter_k];
				}
				sum += sum_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
			model.addQConstr(sum <= 1);
		}

		// Constraint z
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			for (int iter_k = 1; iter_k <= param.K - 1; iter_k++) {
				model.addConstr(z[iter_m][iter_k] >= z[iter_m][iter_k + 1]);
			}
		}

		// Constraint y and z
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			model.addConstr(y[iter_m] >= z[iter_m][1]);
		}

		// Constraint x^ and z
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			GRBLinExpr sumZi = 0;
			for (int iter_k = 1; iter_k <= param.K; iter_k++) {
				sumZi += z[iter_m][iter_k];
			}
			sumZi = sumZi * (data.U[iter_m] - data.L[iter_m]) / param.K;
			sumZi += data.L[iter_m];
			model.addConstr(x[iter_m] == sumZi);
		}

		// Constraint Ax + By <= D : sum_Y <= M
		GRBLinExpr sumY = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sumY += y[iter_m];
		}
		model.addConstr(sumY <= param.M);

		//// Constraint Ax + By <= D : sum_X <= W
		//GRBLinExpr sumX = 0;
		//for (int iter_m = 0; iter_m < data.m; iter_m++) {
		//    sumX += x[iter_m];
		//}
		//model.addConstr(sumX <= data.W);

		// Constraint Ax + By <= D : sum_X <= W
		GRBQuadExpr sumX = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sumX += x[iter_m] * y[iter_m];
		}
		model.addQConstr(sumX <= data.W);

		//// constr y and x
		//for (int iter_m = 0; iter_m < data.m; iter_m++) {
		//    model.addConstr(x[iter_m] <= data.U[iter_m] * y[iter_m] + data.L[iter_m]);
		//}

		// Constraint McCormick for w*z
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					model.addQConstr(w[iter_t] * z[iter_m][iter_k] <= mcCormick.ub_w_z_1[iter_t][iter_m][iter_k] * z[iter_m][iter_k]);
					model.addQConstr(w[iter_t] * z[iter_m][iter_k] >= mcCormick.lb_w_z_1[iter_t][iter_m][iter_k] * z[iter_m][iter_k]);
					model.addQConstr(w[iter_t] * z[iter_m][iter_k] <= w[iter_t] - mcCormick.lb_w_z_0[iter_t][iter_m][iter_k] * (1 - z[iter_m][iter_k]));
					model.addQConstr(w[iter_t] * z[iter_m][iter_k] >= w[iter_t] - mcCormick.ub_w_z_0[iter_t][iter_m][iter_k] * (1 - z[iter_m][iter_k]));
				}
			}
		}

		// Constraint McCormick for w*y
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				model.addQConstr(w[iter_t] * y[iter_m] <= mcCormick.ub_w_y_1[iter_t][iter_m] * y[iter_m]);
				model.addQConstr(w[iter_t] * y[iter_m] >= mcCormick.lb_w_y_1[iter_t][iter_m] * y[iter_m]);
				model.addQConstr(w[iter_t] * y[iter_m] <= w[iter_t] - mcCormick.lb_w_y_0[iter_t][iter_m] * (1 - y[iter_m]));
				model.addQConstr(w[iter_t] * y[iter_m] >= w[iter_t] - mcCormick.ub_w_y_0[iter_t][iter_m] * (1 - y[iter_m]));
			}
		}

		// Objective function
		GRBQuadExpr obj = 0;
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			obj += data.a[iter_t] * w[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				obj += w[iter_t] * y[iter_m] * param.gl[iter_t][iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBQuadExpr obj_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					obj_gamma += param.gamma_g[iter_t][iter_m][iter_k] * w[iter_t] * z[iter_m][iter_k];
				}
				obj += obj_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
		}
		model.setObjective(obj, GRB_MAXIMIZE);

		// Parameters
		model.set(GRB_DoubleParam_TimeLimit, time_limit);
		//model.set(GRB_DoubleParam_MIPGap, 1e-8);
		model.set(GRB_IntParam_MIQCPMethod, 1);

		auto start = chrono::steady_clock::now(); //get start time
		auto end = chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;

		model.optimize();

		end = chrono::steady_clock::now();
		elapsed_seconds = end - start;
		time_for_solve = elapsed_seconds.count();
		cout << "Status: " << model.get(GRB_IntAttr_Status) << endl;
		cout << "Objective value: " << setprecision(8) << model.get(GRB_DoubleAttr_ObjVal) << endl;

		obj_val_gurobi = model.get(GRB_DoubleAttr_ObjVal);

		//get value x^
		vector<double> ansX(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			ansX[iter_m] = x[iter_m].get(GRB_DoubleAttr_X);
		}

		double sum_X = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sum_X += ansX[iter_m];
		}
		cout << "sum X = " << sum_X << endl;

		//get value y
		vector<bool> ansY(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			if (y[iter_m].get(GRB_DoubleAttr_X) > 0.5) ansY[iter_m] = 1;
		}

		//get value z
		vector<vector<bool>> ansZ;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			vector<bool> ansZi(param.K);
			for (int iter_k = 0; iter_k <= param.K; iter_k++) {
				ansZi[iter_k] = 0;
			}
			ansZ.push_back(ansZi);
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			for (int iter_k = 1; iter_k <= param.K; iter_k++) {
				if (z[iter_m][iter_k].get(GRB_DoubleAttr_X) > 0.5) ansZ[iter_m][iter_k] = 1;
				else ansZ[iter_m][iter_k] = 0;
			}
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

		gap = model.get(GRB_DoubleAttr_MIPGap);

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "y[" << iter_m << "] = " << ansY[iter_m] << endl;
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "x[" << iter_m << "] = " << setprecision(5) << ansX[iter_m] << endl;
		}

		cout << "Solving, it took " << elapsed_seconds.count() << " seconds" << endl;
		ofstream report_results(output, ofstream::out);
		report_results << data.m << "," << data.W << "," << param.M
			<< setprecision(8) << fixed << "," << -obj_val_gurobi
			<< "," << obj_val_true << "," << time_for_solve << "," << gap << endl;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			report_results << "y[" << iter_m << "] = " << ansY[iter_m] << endl;
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			report_results << "x[" << iter_m << "] = " << setprecision(5) << ansX[iter_m] << endl;
		}
		report_results.close();
	}
	catch (GRBException e) {
		cerr << "Error code = " << e.getErrorCode() << endl;
		cerr << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Exception during optimization" << endl;
	}
}

