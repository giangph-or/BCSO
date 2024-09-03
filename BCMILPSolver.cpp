#include "BCMILPSolver.h"

using namespace std;

//===========Assortment + Price===========//

BCMILPSolverAssortment::BCMILPSolverAssortment() {

}

BCMILPSolverAssortment::BCMILPSolverAssortment(DataAssortment data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamAssortment(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();

	start = chrono::steady_clock::now(); //get start time
	mcCormick = McCormick(data, param, true);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_mc = elapsed_seconds.count();
}

CBMILPAssortment::CBMILPAssortment(GRBVar* cb_w, GRBVar* cb_y, GRBVar** cb_z, int T, int m, int K, vector<double> b, vector<double> L, vector<double> U, vector<vector<double>> hl, vector<vector<vector<double>>> gamma_h) {
	w = cb_w;
	y = cb_y;
	z = cb_z;
	cb_T = T;
	cb_m = m;
	cb_K = K;
	cb_b = b;
	cb_L = L;
	cb_U = U;
	cb_hl = hl;
	cb_gamma_h = gamma_h;
}

void BCMILPSolverAssortment::solve(string output, int time_limit) {
	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		// w
		GRBVar* w;
		w = new GRBVar[data.T];
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			w[iter_t] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		}

		// z
		GRBVar** z;
		z = new GRBVar * [data.m];
		for (int iter_m = 0; iter_m < data.m; ++iter_m)
			z[iter_m] = new GRBVar[param.K + 1];
		for (int iter_m = 0; iter_m < data.m; ++iter_m)
			for (int k = 0; k < param.K + 1; ++k)
				z[iter_m][k] = model.addVar(0, 1, 0, GRB_BINARY);

		// y
		GRBVar* y;
		y = new GRBVar[data.m];
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			y[iter_m] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		// x^ variables
		vector<GRBVar> x(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// u variables
		vector<vector<vector<GRBVar>>> u(data.T, vector<vector<GRBVar>>(data.m, vector<GRBVar>(param.K + 1)));
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				for (int iter_k = 0; iter_k <= param.K; iter_k++) {
					u[iter_t][iter_m][iter_k] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
				}
			}
		}

		// v variables
		vector<vector<GRBVar>> v(data.T, vector<GRBVar>(data.m));
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				v[iter_t][iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
			}
		}

		// Constraint w
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			GRBLinExpr sum = data.b[iter_t] * w[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				sum += param.hl[iter_t][iter_m] * v[iter_t][iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBLinExpr sum_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					sum_gamma += param.gamma_h[iter_t][iter_m][iter_k] * u[iter_t][iter_m][iter_k];
				}
				sum += sum_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
			model.addConstr(sum <= 1);
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

		// Constraint McCormick for u
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					model.addConstr(u[iter_t][iter_m][iter_k] <= mcCormick.ub_w_z_1[iter_t][iter_m][iter_k] * z[iter_m][iter_k]);
					model.addConstr(u[iter_t][iter_m][iter_k] >= mcCormick.lb_w_z_1[iter_t][iter_m][iter_k] * z[iter_m][iter_k]);
					model.addConstr(u[iter_t][iter_m][iter_k] <= w[iter_t] - mcCormick.lb_w_z_0[iter_t][iter_m][iter_k] * (1 - z[iter_m][iter_k]));
					model.addConstr(u[iter_t][iter_m][iter_k] >= w[iter_t] - mcCormick.ub_w_z_0[iter_t][iter_m][iter_k] * (1 - z[iter_m][iter_k]));
				}
			}
		}

		// Constraint McCormick for v
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				model.addConstr(v[iter_t][iter_m] <= mcCormick.ub_w_y_1[iter_t][iter_m] * y[iter_m]);
				model.addConstr(v[iter_t][iter_m] >= mcCormick.lb_w_y_1[iter_t][iter_m] * y[iter_m]);
				model.addConstr(v[iter_t][iter_m] <= w[iter_t] - mcCormick.lb_w_y_0[iter_t][iter_m] * (1 - y[iter_m]));
				model.addConstr(v[iter_t][iter_m] >= w[iter_t] - mcCormick.ub_w_y_0[iter_t][iter_m] * (1 - y[iter_m]));
			}
		}

		// Objective function
		GRBLinExpr obj = 0;
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			obj += data.a[iter_t] * w[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				obj += v[iter_t][iter_m] * param.gl[iter_t][iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBLinExpr obj_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					obj_gamma += param.gamma_g[iter_t][iter_m][iter_k] * u[iter_t][iter_m][iter_k];
				}
				obj += obj_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
		}
		model.setObjective(obj, GRB_MAXIMIZE);

		// Parameters
		model.set(GRB_DoubleParam_TimeLimit, time_limit);
		//model.set(GRB_DoubleParam_MIPGap, 1e-8);
		model.set(GRB_IntParam_LazyConstraints, 1);
		model.set(GRB_IntParam_PreCrush, 1);
		model.set(GRB_IntParam_MIQCPMethod, 1);

		auto start = chrono::steady_clock::now(); //get start time
		auto end = chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;

		CBMILPAssortment cb = CBMILPAssortment(w, y, z, data.T, data.m, param.K, data.b, data.L, data.U, param.hl, param.gamma_h);
		model.setCallback(&cb);

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

void CBMILPAssortment::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double* initial_y = new double[cb_m];
			double** initial_z = new double* [cb_m];
			initial_y = getSolution(y, cb_m);
			for (int m = 0; m < cb_m; ++m) {
				initial_z[m] = getSolution(z[m], cb_K + 1);
			}

			vector<double> initial_denominator(cb_T, 0);
			for (int t = 0; t < cb_T; ++t) {
				initial_denominator[t] += cb_b[t];
				for (int m = 0; m < cb_m; ++m) {
					initial_denominator[t] += cb_hl[t][m] * initial_y[m];
					double gamma = 0;
					for (int k = 1; k <= cb_K; ++k)
						gamma += cb_gamma_h[t][m][k] * initial_z[m][k];
					initial_denominator[t] += gamma * (cb_U[m] - cb_L[m]) / cb_K;
				}
			}

			//Outer-cuts for w
			vector<double> partial_w(cb_T, 0);
			for (int t = 0; t < cb_T; ++t)
				partial_w[t] = 1 / initial_denominator[t];

			vector<vector<double>> subgradient_y(cb_T);
			for (int t = 0; t < cb_T; ++t)
				subgradient_y[t].resize(cb_m, 0);
			for (int t = 0; t < cb_T; ++t)
				for (int m = 0; m < cb_m; ++m)
					subgradient_y[t][m] -= cb_hl[t][m] / (initial_denominator[t] * initial_denominator[t]);

			vector<vector<vector<double>>> subgradient_z(cb_T);
			for (int t = 0; t < cb_T; ++t)
				subgradient_z[t].resize(cb_m);
			for (int t = 0; t < cb_T; ++t)
				for (int m = 0; m < cb_m; ++m)
					subgradient_z[t][m].resize(cb_K + 1, 0);
			for (int t = 0; t < cb_T; ++t)
				for (int m = 0; m < cb_m; ++m)
					for (int k = 1; k <= cb_K; ++k)
						subgradient_z[t][m][k] -= (cb_U[m] - cb_L[m]) / cb_K * cb_gamma_h[t][m][k] / (initial_denominator[t] * initial_denominator[t]);

			for (int t = 0; t < cb_T; ++t) {
				GRBLinExpr gradient;
				for (int m = 0; m < cb_m; ++m) {
					gradient += subgradient_y[t][m] * (y[m] - initial_y[m]);
					for (int k = 1; k <= cb_K; ++k)
						gradient += subgradient_z[t][m][k] * (z[m][k] - initial_z[m][k]);
				}
				addLazy(w[t] >= partial_w[t] + gradient);
				cout << "addLazy: " << t << endl;
			}

		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}

//===========Facility Location + Cost===========//

OASolverFacility::OASolverFacility() {

}

OASolverFacility::OASolverFacility(DataFacility data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamFacility(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();
}

void OASolverFacility::solve(string output, int time_limit) {
	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		// theta
		GRBVar* theta;
		theta = new GRBVar[data.T];
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			theta[iter_t] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		}

		// z
		GRBVar** z;
		z = new GRBVar * [data.m];
		for (int iter_m = 0; iter_m < data.m; ++iter_m)
			z[iter_m] = new GRBVar[param.K + 1];
		for (int iter_m = 0; iter_m < data.m; ++iter_m)
			for (int k = 0; k < param.K + 1; ++k)
				z[iter_m][k] = model.addVar(0, 1, 0, GRB_BINARY);

		// y
		GRBVar* y;
		y = new GRBVar[data.m];
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			y[iter_m] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		// x variables
		vector<GRBVar> x(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
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

		// Constraint x and z
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
		model.addQConstr(sumX <= data.C);

		//// constr y and x
		//for (int iter_m = 0; iter_m < data.m; iter_m++) {
		//    model.addConstr(x[iter_m] <= data.U[iter_m] * y[iter_m] + data.L[iter_m]);
		//}

		// Objective function
		GRBLinExpr obj = 0;
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			obj += data.a[iter_t] - data.a[iter_t] * theta[iter_t];
		}
		model.setObjective(obj, GRB_MAXIMIZE);

		auto start = chrono::steady_clock::now(); //get start time
		auto end = chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;

		double run_time = time_limit - elapsed_seconds.count();

		double stop_param = 1e-4;
		double best_obj = 0;
		vector<double> best_y(data.m, 0);
		vector<double> best_x(data.m, 0);
		vector<vector<double>> best_z(data.m);
		for (int iter_m = 0; iter_m < data.m; ++iter_m)
			best_z[iter_m].resize(param.K + 1, 0);

		vector<double> initial_y(data.m, 0);
		vector<double> initial_x(data.m, 0);
		vector<vector<double>> initial_z(data.m);
		for (int iter_m = 0; iter_m < data.m; ++iter_m)
			initial_z[iter_m].resize(param.K + 1, 0);

		obj_val_true = 0;
		obj_val_gurobi = 1;

		while (obj_val_true < obj_val_gurobi - stop_param) {
			vector<double> initial_denominator(data.T, 0);
			for (int t = 0; t < data.T; ++t) {
				initial_denominator[t] += data.b[t];
				for (int m = 0; m < data.m; ++m) {
					initial_denominator[t] += param.hl[t][m] * initial_y[m];
					double gamma = 0;
					for (int k = 1; k <= param.K; ++k) {
						gamma += param.gamma_h[t][m][k] * initial_z[m][k];
					}
					initial_denominator[t] += gamma * (data.U[m] - data.L[m]) / param.K;
				}
			}

			//Outer-cuts for theta
			vector<double> partial_theta(data.T, 0);
			for (int t = 0; t < data.T; ++t)
				partial_theta[t] = 1 / initial_denominator[t];

			vector<vector<double>> subgradient_y(data.T);
			for (int t = 0; t < data.T; ++t)
				subgradient_y[t].resize(data.m, 0);
			for (int t = 0; t < data.T; ++t)
				for (int m = 0; m < data.m; ++m)
					subgradient_y[t][m] -= param.hl[t][m] / (initial_denominator[t] * initial_denominator[t]);

			vector<vector<vector<double>>> subgradient_z(data.T);
			for (int t = 0; t < data.T; ++t)
				subgradient_z[t].resize(data.m);
			for (int t = 0; t < data.T; ++t)
				for (int m = 0; m < data.m; ++m)
					subgradient_z[t][m].resize(param.K + 1, 0);
			for (int t = 0; t < data.T; ++t)
				for (int m = 0; m < data.m; ++m)
					for (int k = 1; k <= param.K; ++k)
						subgradient_z[t][m][k] -= (data.U[m] - data.L[m]) / param.K * param.gamma_h[t][m][k] / (initial_denominator[t] * initial_denominator[t]);

			for (int t = 0; t < data.T; ++t) {
				GRBLinExpr gradient;
				for (int m = 0; m < data.m; ++m) {
					gradient += subgradient_y[t][m] * (y[m] - initial_y[m]);
					for (int k = 1; k <= param.K; ++k)
						gradient += subgradient_z[t][m][k] * (z[m][k] - initial_z[m][k]);
				}
				model.addConstr(theta[t] >= partial_theta[t] + gradient);
			}

			model.set(GRB_DoubleParam_TimeLimit, run_time);
			model.set(GRB_IntParam_MIQCPMethod, 1);

			model.optimize();

			if (model.get(GRB_IntAttr_SolCount) > 0) {
				obj_val_gurobi = model.get(GRB_DoubleAttr_ObjVal);

				for (int iter_m = 0; iter_m < data.m; iter_m++)
					initial_x[iter_m] = x[iter_m].get(GRB_DoubleAttr_X);

				for (int iter_m = 0; iter_m < data.m; iter_m++)
					if (y[iter_m].get(GRB_DoubleAttr_X) > 0.5)
						initial_y[iter_m] = 1;

				for (int iter_m = 0; iter_m < data.m; iter_m++)
					for (int iter_k = 1; iter_k <= param.K; iter_k++)
						if (z[iter_m][iter_k].get(GRB_DoubleAttr_X) > 0.5) 
							initial_z[iter_m][iter_k] = 1;
						else 
							initial_z[iter_m][iter_k] = 0;

				//Compute true obj val
				obj_val_true = 0;
				for (int iter_t = 0; iter_t < data.T; iter_t++) {
					double tmp1 = data.a[iter_t];
					double tmp2 = data.b[iter_t];
					for (int iter_m = 0; iter_m < data.m; iter_m++) {
						tmp2 += initial_y[iter_m] * param.h_func(data.eta, data.kappa, iter_t, iter_m, initial_x[iter_m]);
					}
					obj_val_true += tmp1 - tmp1 / tmp2;
				}
				if (obj_val_true >= obj_val_gurobi - stop_param) {
					best_obj = obj_val_gurobi;
					best_x = initial_x;
					best_y = initial_y;
					best_z = initial_z;
				}

				auto time_now = std::chrono::steady_clock::now(); //get now time
				std::chrono::duration<double> after_cut = time_now - start;

				if (after_cut.count() > time_limit) break;
				run_time = time_limit - after_cut.count();
			}
			else {
				auto end = chrono::steady_clock::now();
				chrono::duration<double> elapsed_seconds = end - start;
				time_for_solve = elapsed_seconds.count();
				break;
			}
		}

		end = chrono::steady_clock::now();
		elapsed_seconds = end - start;
		time_for_solve = elapsed_seconds.count();
		cout << "Objective value: " << setprecision(8) << best_obj << endl;

		//get value x^
		vector<double> ansX = best_x;
		double sum_X = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sum_X += ansX[iter_m];
		}
		cout << "sum X = " << sum_X << endl;

		//get value y
		vector<double> ansY = best_y;

		//get value z
		vector<vector<double>> ansZ = best_z;

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "y[" << iter_m << "] = " << ansY[iter_m] << endl;
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "x[" << iter_m << "] = " << setprecision(5) << ansX[iter_m] << endl;
		}

		cout << "Solving, it took " << elapsed_seconds.count() << " seconds" << endl;
		ofstream report_results(output, ofstream::out);
		report_results << data.m << "," << data.C << "," << param.M
			<< setprecision(8) << fixed << "," << best_obj
			<< "," << best_obj << "," << time_for_solve << endl;
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

