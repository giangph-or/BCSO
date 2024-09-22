#include "BilinearNew.h"

using namespace std;

//===========Assortment + Price===========//

BilinearNewAssortment::BilinearNewAssortment() {

}

BilinearNewAssortment::BilinearNewAssortment(DataAssortment data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamAssortment(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();
}

void BilinearNewAssortment::solve(string output, int time_limit) {
	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		// n variables: numerator
		vector<GRBVar> n(data.T);
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			n[iter_t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
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

		// x^ variables
		vector<GRBVar> x(data.m);
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			x[iter_m] = model.addVar(data.L[iter_m], data.U[iter_m], 0.0, GRB_CONTINUOUS);
		}

		// Constraint n
		for (int iter_t = 0; iter_t < data.T; iter_t++) {
			GRBLinExpr sum = data.a[iter_t];
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				sum += param.gl[iter_t][iter_m] * y[iter_m];
			}
			for (int iter_m = 0; iter_m < data.m; iter_m++) {
				GRBLinExpr sum_gamma = 0;
				for (int iter_k = 1; iter_k <= param.K; iter_k++) {
					sum_gamma += param.gamma_g[iter_t][iter_m][iter_k] * z[iter_m][iter_k];
				}
				sum += sum_gamma * (data.U[iter_m] - data.L[iter_m]) / param.K;
			}
			model.addConstr(sum == n[iter_t]);
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

		// f=n/d
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			model.addQConstr(f[iter_t] * d[iter_t] <= n[iter_t]);

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

		//// constr y and x
		//for (int iter_m = 0; iter_m < data.m; iter_m++) {
		//    model.addConstr(x[iter_m] <= data.U[iter_m] * y[iter_m] + data.L[iter_m]);
		//}

		//// Constraint Ax + By <= D : sum_Y <= M
		//GRBLinExpr sumY = 0;
		//for (int iter_m = 0; iter_m < data.m; iter_m++) {
		//	sumY += y[iter_m];
		//}
		//model.addConstr(sumY <= param.M);

		//// Constraint Ax + By <= D : sum_X <= W
		//GRBQuadExpr sumX = 0;
		//for (int iter_m = 0; iter_m < data.m; iter_m++) {
		//	sumX += x[iter_m] * y[iter_m];
		//}
		//model.addQConstr(sumX <= data.W);

		// Constraint Ax + By <= D : sum_X <= W
		GRBLinExpr sumX = 0;
		for (int iter_m = 0; iter_m < data.m; iter_m++) {
			sumX += x[iter_m];
		}
		model.addConstr(sumX <= data.W);

		// Objective function
		GRBLinExpr obj = 0;
		for (int iter_t = 0; iter_t < data.T; iter_t++)
			obj += f[iter_t];

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

//===========Facility Location + Cost===========//

BilinearNewFacility::BilinearNewFacility() {

}

BilinearNewFacility::BilinearNewFacility(DataFacility data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamFacility(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();
}

void BilinearNewFacility::solve(string output, int time_limit) {
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

		//// constr y and x
		//for (int iter_m = 0; iter_m < data.m; iter_m++) {
		//    model.addConstr(x[iter_m] <= data.U[iter_m] * y[iter_m] + data.L[iter_m]);
		//}

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
