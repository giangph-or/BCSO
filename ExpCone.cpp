#include "ExpCone.h"

using namespace std;

//===========Assortment + Price===========//

ExpConeAssortment::ExpConeAssortment() {

}

ExpConeAssortment::ExpConeAssortment(DataAssortment data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamAssortment(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();
}

CBExpAssortment::CBExpAssortment(GRBVar* cb_y, GRBVar* cb_n, GRBVar* cb_d, GRBVar* cb_f, GRBVar** cb_z,
	int cb_T, int cb_m, int cb_K, vector<double> cb_b, vector<double> cb_U, vector<double> cb_L,
	vector<vector<double>> cb_hl, vector<vector<vector<double>>> cb_gamma_h) {
	y = cb_y;
	n = cb_n;
	d = cb_d;
	f = cb_f;
	z = cb_z;
	T_ = cb_T;
	m_ = cb_m;
	K_ = cb_K;
	b_ = cb_b;
	U_ = cb_U;
	L_ = cb_L;
	hl_ = cb_hl;
	gamma_h_ = cb_gamma_h;
}

void ExpConeAssortment::solve(string output, int time_limit) {
	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		// n variables: numerator
		vector<GRBVar> n(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			n[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// e^n variables: numerator
		vector<GRBVar> e_n(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			e_n[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// d variables: denominator
		vector<GRBVar> d(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			d[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// f variables: fraction
		vector<GRBVar> f(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			f[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
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

		// x variables
		vector<GRBVar> x(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			x[iter_m] = model.addVar(data.L[iter_m], data.U[iter_m], 0.0, GRB_CONTINUOUS);
		}

		//Constraint e^n and n
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			model.addGenConstrExp(n[iter_t], e_n[iter_t]);

		double max_U = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++)
			if (data.U[iter_m] > max_U)
				max_U = data.U[iter_m];

		// Constraint e^n >= ...
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			GRBLinExpr sum = max_U * data.b[iter_t] - data.a[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				sum += (max_U * param.hl[iter_t][iter_m] - param.gl[iter_t][iter_m]) * y[iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBLinExpr sum_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					sum_gamma += (max_U * param.gamma_h[iter_t][iter_m][iter_k] - param.gamma_g[iter_t][iter_m][iter_k]) * z[iter_m][iter_k];
				}
				sum += sum_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
			model.addConstr(sum == e_n[iter_t]);
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

		// Constraint Ax + By <= D : sum_X <= W
		GRBQuadExpr sumX = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sumX += x[iter_m] * y[iter_m];
		}
		model.addQConstr(sumX <= data.W);

		// Objective function
		GRBLinExpr obj = 0;
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			obj += f[iter_t];

		model.setObjective(obj, GRB_MINIMIZE);

		auto start = chrono::steady_clock::now(); //get start time
		auto end = chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;

		double run_time = time_limit - elapsed_seconds.count();

		double stop_param = 1e-5;
		double best_obj = 0;
		double best_obj_true = 0;
		vector<double> best_y(data.m, 0);
		vector<double> best_x(data.m, 0);

		vector<double> initial_y(data.m, 0);
		vector<double> initial_x(data.m, 0);

		vector<double> initial_d(data.T, 0);
		vector<double> initial_n(data.T, 0);
		//for (int iter_t = 0; iter_t < data.T; iter_t++) {
		//	initial_d[iter_t] = data.b[iter_t];
		//	initial_n[iter_t] = max_U * data.b[iter_t] - data.a[iter_t];
		//}

		obj_val_true = 1;
		obj_val_gurobi = 0;
		int iter = 1;
		while (obj_val_true > obj_val_gurobi + stop_param) {
			// Constraint d
			for (int iter_t = 0; iter_t < data.T; iter_t++) {
				GRBLinExpr sum = data.b[iter_t];
				for (int iter_m = 0; iter_m < data.m; iter_m++) {
					sum += param.hl[iter_t][iter_m] * y[iter_m];
				}
				for (int iter_m = 0; iter_m < data.m; iter_m++) {
					GRBLinExpr sum_gamma = 0;
					for (int iter_k = 1; iter_k <= param.K; iter_k++) {
						sum_gamma += param.gamma_h[iter_t][iter_m][iter_k] * z[iter_m][iter_k];
					}
					sum += sum_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
				}
				model.addConstr(exp(initial_d[iter_t]) + exp(initial_d[iter_t]) * (d[iter_t] - initial_d[iter_t]) <= sum);
			}

			// Constraint f
			for (int iter_t = 0; iter_t < data.T; iter_t++)
				model.addConstr(f[iter_t] >=
					exp(initial_n[iter_t] - initial_d[iter_t])
					+ exp(initial_n[iter_t] - initial_d[iter_t]) * (n[iter_t] - initial_n[iter_t])
					- exp(initial_n[iter_t] - initial_d[iter_t]) * (d[iter_t] - initial_d[iter_t])
				);

			model.set(GRB_IntParam_OutputFlag, 0);
			model.set(GRB_DoubleParam_TimeLimit, run_time);
			model.set(GRB_IntParam_MIQCPMethod, 1);
			model.set(GRB_IntParam_FuncPieces, 30);
			model.set(GRB_DoubleParam_MIPGap, 1e-4);

			model.optimize();

			if (model.get(GRB_IntAttr_SolCount) > 0) {
				obj_val_gurobi = model.get(GRB_DoubleAttr_ObjVal);

				for (int iter_m = 0; iter_m < data.m; iter_m++)
					initial_x[iter_m] = x[iter_m].get(GRB_DoubleAttr_X);

				for (int iter_m = 0; iter_m < data.m; iter_m++)
					if (y[iter_m].get(GRB_DoubleAttr_X) > 0.5)
						initial_y[iter_m] = 1;
					else
						initial_y[iter_m] = 0;

				for (int iter_t = 0; iter_t < data.T; iter_t++) {
					initial_n[iter_t] = n[iter_t].get(GRB_DoubleAttr_X);
					initial_d[iter_t] = d[iter_t].get(GRB_DoubleAttr_X);
				}

				//Compute true obj val
				obj_val_true = 0;
				for (int iter_t = 0; iter_t < data.T; iter_t++)
					obj_val_true += exp(initial_n[iter_t] - initial_d[iter_t]);
				

				if (obj_val_true <= obj_val_gurobi + stop_param) {
					best_obj = obj_val_gurobi;
					best_x = initial_x;
					best_y = initial_y;
					for (int iter_t = 0; iter_t < data.T; iter_t++) {
						double tmp1 = data.a[iter_t];
						for (int iter_m = 0; iter_m < data.m; iter_m++) {
							tmp1 += best_y[iter_m] * param.g_func(data.eta, data.kappa, iter_t, iter_m, best_x[iter_m]);
						}
						double tmp2 = data.b[iter_t];
						for (int iter_m = 0; iter_m < data.m; iter_m++) {
							tmp2 += best_y[iter_m] * param.h_func(data.eta, data.kappa, iter_t, iter_m, best_x[iter_m]);
						}
						best_obj_true += tmp1 / tmp2;
					}
				}

				auto time_now = std::chrono::steady_clock::now(); //get now time
				std::chrono::duration<double> after_cut = time_now - start;

				if (after_cut.count() > time_limit) break;
				run_time = time_limit - after_cut.count();

				cout << "Iteration " << iter << ":\n";
				cout << "Obj Gurobi: " << obj_val_gurobi << endl;
				cout << "Obj sub: " << obj_val_true << endl;
				cout << "Time: " << run_time << endl;
				cout << endl;
				iter++;
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
		cout << "Objective value: " << setprecision(8) << best_obj_true << endl;

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "y[" << iter_m << "] = " << best_y[iter_m] << endl;
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "x[" << iter_m << "] = " << setprecision(5) << best_x[iter_m] << endl;
		}

		cout << "Solving, it took " << elapsed_seconds.count() << " seconds" << endl;
		ofstream report_results(output, ofstream::out);
		report_results << data.m << "," << data.W << "," << param.M
			<< setprecision(8) << fixed << "," << -obj_val_gurobi
			<< "," << best_obj_true << "," << time_for_solve << endl;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			report_results << "y[" << iter_m << "] = " << best_y[iter_m] << endl;
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			report_results << "x[" << iter_m << "] = " << setprecision(5) << best_x[iter_m] << endl;
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

void ExpConeAssortment::solve_BC(string output, int time_limit) {
	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		// n variables: numerator
		GRBVar* n;
		n = new GRBVar[data.T];
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			n[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// e^n variables: numerator
		vector<GRBVar> e_n(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			e_n[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// d variables: denominator
		GRBVar* d;
		d = new GRBVar[data.T];
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			d[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// f variables: fraction
		GRBVar* f;
		f = new GRBVar[data.T];
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			f[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

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
			x[iter_m] = model.addVar(data.L[iter_m], data.U[iter_m], 0.0, GRB_CONTINUOUS);
		}

		//Constraint e^n and n
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			model.addGenConstrExp(n[iter_t], e_n[iter_t]);

		double max_U = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++)
			if (data.U[iter_m] > max_U)
				max_U = data.U[iter_m];

		// Constraint e^n >= ...
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			GRBLinExpr sum = max_U * data.b[iter_t] - data.a[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				sum += (max_U * param.hl[iter_t][iter_m] - param.gl[iter_t][iter_m]) * y[iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBLinExpr sum_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					sum_gamma += (max_U * param.gamma_h[iter_t][iter_m][iter_k] - param.gamma_g[iter_t][iter_m][iter_k]) * z[iter_m][iter_k];
				}
				sum += sum_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
			model.addConstr(sum == e_n[iter_t]);
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

		// Constraint Ax + By <= D : sum_X <= W
		GRBQuadExpr sumX = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sumX += x[iter_m] * y[iter_m];
		}
		model.addQConstr(sumX <= data.W);

		// Objective function
		GRBLinExpr obj = 0;
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			obj += f[iter_t];

		model.setObjective(obj, GRB_MINIMIZE);

		vector<double> initial_y(data.m, 0);
		vector<double> initial_x(data.m, 0);
		double obj_true = 0;

		//model.set(GRB_IntParam_OutputFlag, 0);
		model.set(GRB_DoubleParam_TimeLimit, 3600);
		model.set(GRB_IntParam_MIQCPMethod, 1);
		model.set(GRB_IntParam_FuncPieces, 100);
		//model.set(GRB_DoubleParam_MIPGap, 1e-4);
		model.set(GRB_IntParam_LazyConstraints, 1);

		auto start = chrono::steady_clock::now(); //get start time
		auto end = chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;

		CBExpAssortment cb = CBExpAssortment(y, n, d, f, z, data.T, data.m, param.K, data.b, data.U, data.L, param.hl, param.gamma_h);
		model.setCallback(&cb);

		model.optimize();

		if (model.get(GRB_IntAttr_SolCount) > 0) {
			obj_val_gurobi = model.get(GRB_DoubleAttr_ObjVal);

			for (int iter_m = 0; iter_m < data.m; iter_m++)
				initial_x[iter_m] = x[iter_m].get(GRB_DoubleAttr_X);

			for (int iter_m = 0; iter_m < data.m; iter_m++)
				if (y[iter_m].get(GRB_DoubleAttr_X) > 0.5)
					initial_y[iter_m] = 1;
				else
					initial_y[iter_m] = 0;

			for (int iter_t = 0; iter_t < data.T; iter_t++) {
				double tmp1 = data.a[iter_t];
				for (int iter_m = 0; iter_m < data.m; iter_m++) {
					tmp1 += initial_y[iter_m] * param.g_func(data.eta, data.kappa, iter_t, iter_m, initial_x[iter_m]);
				}
				double tmp2 = data.b[iter_t];
				for (int iter_m = 0; iter_m < data.m; iter_m++) {
					tmp2 += initial_y[iter_m] * param.h_func(data.eta, data.kappa, iter_t, iter_m, initial_x[iter_m]);
				}
				obj_true += tmp1 / tmp2;
			}
		}

		end = chrono::steady_clock::now();
		elapsed_seconds = end - start;
		time_for_solve = elapsed_seconds.count();
		cout << "Objective value: " << setprecision(8) << obj_true << endl;

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "y[" << iter_m << "] = " << initial_y[iter_m] << endl;
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			cout << "x[" << iter_m << "] = " << setprecision(5) << initial_x[iter_m] << endl;
		}

		cout << "Solving, it took " << elapsed_seconds.count() << " seconds" << endl;
		ofstream report_results(output, ofstream::out);
		report_results << data.m << "," << data.W << "," << param.M
			<< setprecision(8) << fixed << "," << -obj_val_gurobi
			<< "," << obj_true << "," << time_for_solve << "," << model.get(GRB_DoubleAttr_MIPGap) << endl;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			report_results << "y[" << iter_m << "] = " << initial_y[iter_m] << endl;
		}

		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			report_results << "x[" << iter_m << "] = " << setprecision(5) << initial_x[iter_m] << endl;
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

void CBExpAssortment::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double* initial_y = new double[m_];
			double** initial_z = new double* [m_];
			double* initial_n = new double[T_];
			double* initial_d = new double[T_];
			initial_y = getSolution(y, m_);
			initial_n = getSolution(n, T_);
			initial_d = getSolution(d, T_);
			for (int m = 0; m < m; ++m)
				initial_z[m] = getSolution(z[m], K_ + 1);

			// Constraint d
			for (int iter_t = 0; iter_t < T_; iter_t++) {
				GRBLinExpr sum = b_[iter_t];
				for (int iter_m = 0; iter_m < m_; iter_m++) {
					sum += hl_[iter_t][iter_m] * y[iter_m];
				}
				for (int iter_m = 0; iter_m < m_; iter_m++) {
					GRBLinExpr sum_gamma = 0;
					for (int iter_k = 1; iter_k <= K_; iter_k++) {
						sum_gamma += gamma_h_[iter_t][iter_m][iter_k] * z[iter_m][iter_k];
					}
					sum += sum_gamma * (U_[iter_m] - L_[iter_m]) / K_;
				}
				addLazy(exp(initial_d[iter_t]) + exp(initial_d[iter_t]) * (d[iter_t] - initial_d[iter_t]) <= sum);
			}

			// Constraint f
			for (int iter_t = 0; iter_t < T_; iter_t++)
				addLazy(f[iter_t] >=
					exp(initial_n[iter_t] - initial_d[iter_t])
					+ exp(initial_n[iter_t] - initial_d[iter_t]) * (n[iter_t] - initial_n[iter_t])
					- exp(initial_n[iter_t] - initial_d[iter_t]) * (d[iter_t] - initial_d[iter_t])
				);
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

ExpConeFacility::ExpConeFacility() {

}

ExpConeFacility::ExpConeFacility(DataFacility data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamFacility(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();
}

void ExpConeFacility::solve(string output, int time_limit) {
	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		// d variables: denominator
		vector<GRBVar> d(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			d[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// f variables: fraction
		vector<GRBVar> f(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			f[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// y variables
		vector<GRBVar> y(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			y[iter_m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}

		// x variables
		vector<GRBVar> x(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			x[iter_m] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		}

		// z variables
		vector<vector<GRBVar>> z(data.m, vector<GRBVar>(param.K + 1));
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			for (int k = 0; k <= param.K; k++) {
				z[iter_m][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}

		// Constraint d
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			GRBLinExpr sum = data.b[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				sum += param.hl[iter_t][iter_m] * y[iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBLinExpr sum_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					sum_gamma += param.gamma_h[iter_t][iter_m][iter_k] * z[iter_m][iter_k];
				}
				sum += sum_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
			model.addConstr(sum == d[iter_t]);
		}

		// f=a/d
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			model.addQConstr(f[iter_t] * d[iter_t] >= data.a[iter_t]);

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

		// Constraint Ax + By <= D : sum_X <= W
		GRBQuadExpr sumX = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sumX += x[iter_m] * y[iter_m];
		}
		model.addQConstr(sumX <= data.C);

		// Objective function
		GRBLinExpr obj = 0;
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			obj += data.a[iter_t] - f[iter_t];

		model.setObjective(obj, GRB_MAXIMIZE);

		// Parameters
		model.set(GRB_DoubleParam_TimeLimit, time_limit);
		//model.set(GRB_DoubleParam_MIPGap, 1e-8);
		model.set(GRB_IntParam_MIQCPMethod, 1);
		//model.set(GRB_IntParam_MIPFocus, 3);

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
			double tmp2 = data.b[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				tmp2 += ansY[iter_m] * param.h_func(data.eta, data.kappa, iter_t, iter_m, ansX[iter_m]);
			}
			obj_val_true += tmp1 - tmp1 / tmp2;
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
		report_results << data.m << "," << data.C << "," << param.M
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
	catch (GRBException e) {
		cerr << "Error code = " << e.getErrorCode() << endl;
		cerr << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Exception during optimization" << endl;
	}
}
