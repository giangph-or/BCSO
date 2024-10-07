#include "ScipSolver.h"

ScipSolverAssortment::ScipSolverAssortment() {

}

ScipSolverAssortment::ScipSolverAssortment(DataAssortment data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamAssortment(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();
}

void ScipSolverAssortment::solve(string output, int time_limit) {
	SCIP* scip = NULL;

	SCIPcreate(&scip);
	SCIPincludeDefaultPlugins(scip);
	SCIPcreateProbBasic(scip, "BCSO");

	vector<SCIP_VAR*> x(data.m);
	vector<SCIP_VAR*> y(data.m);

	for (int i = 0; i < data.m; ++i) {
		string varName = "x" + std::to_string(i + 1);
		SCIPcreateVarBasic(scip, &x[i], varName.c_str(), data.L[i], data.U[i], 0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(scip, x[i]);
		varName = "y" + std::to_string(i + 1);
		SCIPcreateVarBasic(scip, &y[i], varName.c_str(), 0, 1, 0, SCIP_VARTYPE_BINARY);
		SCIPaddVar(scip, y[i]);
	}

	//sum of y_i <= M
	SCIP_CONS* ct_sum_y;
	SCIPcreateConsBasicLinear(scip, &ct_sum_y, "sum_y", 0, NULL, NULL, -SCIPinfinity(scip), param.M);
	for (int i = 0; i < data.m; ++i)
		SCIPaddCoefLinear(scip, ct_sum_y, y[i], 1.0);
	SCIPaddCons(scip, ct_sum_y);

	//sum of x_i*y_i <= W
	SCIP_EXPR** expr_x;
	expr_x = new SCIP_EXPR * [data.m];
	SCIP_EXPR** expr_y;
	expr_y = new SCIP_EXPR * [data.m];
	SCIP_EXPR** expr_xy;
	expr_xy = new SCIP_EXPR * [data.m];
	for (int i = 0; i < data.m; ++i) {
		SCIP_EXPR* xy[2];
		SCIPcreateExprVar(scip, &expr_x[i], x[i], NULL, NULL);
		SCIPcreateExprVar(scip, &expr_y[i], y[i], NULL, NULL);
		xy[0] = expr_x[i];
		xy[1] = expr_y[i];
		SCIPcreateExprProduct(scip, &expr_xy[i], 2, xy, 1.0, NULL, NULL);
	}
	SCIP_EXPR* sum_xy;
	SCIPcreateExprSum(scip, &sum_xy, data.m, expr_xy, NULL, 0.0, NULL, NULL);
	SCIP_CONS* ct_sum_xy;
	SCIPcreateConsBasicNonlinear(scip, &ct_sum_xy, "sum_xy", sum_xy, -SCIPinfinity(scip), data.W);
	SCIPaddCons(scip, ct_sum_xy);

	//create objective
	SCIP_EXPR* obj;
	SCIP_EXPR** frac;
	frac = new SCIP_EXPR * [data.T];
	for (int j = 0; j < data.T; ++j) {
		SCIP_EXPR* exprs[2];
		vector<SCIP_EXPR*> exp(data.m);
		vector<SCIP_EXPR*> power(data.m);
		SCIP_EXPR** product_xy;
		product_xy = new SCIP_EXPR * [data.m];
		SCIP_EXPR** product_y;
		product_y = new SCIP_EXPR * [data.m];
		for (int i = 0; i < data.m; ++i) {
			SCIPcreateExprSum(scip, &power[i], 1, &expr_x[i], &data.eta[j][i], data.kappa[j][i], NULL, NULL);
			SCIPcreateExprExp(scip, &exp[i], power[i], NULL, NULL);
			SCIP_EXPR* xye[3];
			xye[0] = expr_x[i];
			xye[1] = expr_y[i];
			xye[2] = exp[i];
			SCIPcreateExprProduct(scip, &product_xy[i], 3, xye, 1, NULL, NULL);
			SCIP_EXPR* ye[2];
			ye[0] = expr_y[i];
			ye[1] = exp[i];
			SCIPcreateExprProduct(scip, &product_y[i], 2, ye, 1, NULL, NULL);
		}
		SCIPcreateExprSum(scip, &exprs[0], data.m, product_xy, NULL, 0, NULL, NULL);
		SCIP_EXPR* denominator;
		SCIPcreateExprSum(scip, &denominator, data.m, product_y, NULL, 0, NULL, NULL);
		SCIP_EXPR* denominator1;
		SCIPcreateExprSum(scip, &denominator1, 1, &denominator, NULL, 1, NULL, NULL);
		SCIPcreateExprPow(scip, &exprs[1], denominator1, -1.0, NULL, NULL);
		SCIPcreateExprProduct(scip, &frac[j], 2, exprs, 1.0, NULL, NULL);
	}
	SCIPcreateExprSum(scip, &obj, data.T, frac, NULL, 0.0, NULL, NULL);

	SCIP_VAR* z;
	SCIPcreateVarBasic(scip, &z, "z", 0, SCIPinfinity(scip), 1, SCIP_VARTYPE_CONTINUOUS);
	SCIPaddVar(scip, z);
	SCIP_EXPR* expr_z;
	SCIPcreateExprVar(scip, &expr_z, z, NULL, NULL);
	SCIP_EXPR* objective[2];
	objective[0] = obj;
	objective[1] = expr_z;
	SCIP_Real coefs[2];
	coefs[0] = 1;
	coefs[1] = -1;
	SCIP_EXPR* sum_obj;
	SCIPcreateExprSum(scip, &sum_obj, 2, objective, coefs, 0.0, NULL, NULL);
	SCIP_CONS* ct_obj;
	SCIPcreateConsBasicNonlinear(scip, &ct_obj, "obj", sum_obj, 0, 0);
	SCIPaddCons(scip, ct_obj);

	SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE);
	SCIPsetRealParam(scip, "limits/time", 3600);
	SCIPsetRealParam(scip, "limits/gap", 1e-4);
	SCIPsetRealParam(scip, "parallel/maxnthreads", 64);
	SCIPsetRealParam(scip, "parallel/minnthreads", 8);
	//SCIPwriteOrigProblem(scip, "model.lp", nullptr, FALSE);

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	SCIPsolve(scip);

	end = chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_solve = elapsed_seconds.count();

	SCIP_SOL* bestsol = SCIPgetBestSol(scip);
	obj_val_scip = SCIPgetSolOrigObj(scip, bestsol);

	//get value x
	vector<double> ansX(data.m);
	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		ansX[iter_m] = SCIPgetSolVal(scip, bestsol, x[iter_m]);
	}

	//get value y
	vector<bool> ansY(data.m);
	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		if (SCIPgetSolVal(scip, bestsol, y[iter_m]) > 0.5) ansY[iter_m] = 1;
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

	gap = SCIPgetGap(scip);

	ofstream report_results(output, ofstream::out);
	report_results << data.m << "," << data.W << "," << param.M
		<< setprecision(8) << fixed << "," << obj_val_scip
		<< "," << obj_val_true << "," << time_for_solve << "," << gap << endl;
	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		report_results << "y[" << iter_m << "] = " << ansY[iter_m] << endl;
	}

	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		report_results << "x[" << iter_m << "] = " << setprecision(5) << ansX[iter_m] << endl;
	}
	report_results.close();
}

ScipSolverFacility::ScipSolverFacility() {

}

ScipSolverFacility::ScipSolverFacility(DataFacility data, int K, double tol_lamda, int M) {
	this->data = data;

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	param = ParamFacility(data, K, tol_lamda, M);
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_param = elapsed_seconds.count();
}

void ScipSolverFacility::solve(string output, int time_limit) {
	SCIP* scip = NULL;

	SCIPcreate(&scip);
	SCIPincludeDefaultPlugins(scip);
	SCIPcreateProbBasic(scip, "BCSO");

	vector<SCIP_VAR*> x(data.m);
	vector<SCIP_VAR*> y(data.m);

	for (int i = 0; i < data.m; ++i) {
		string varName = "x" + std::to_string(i + 1);
		SCIPcreateVarBasic(scip, &x[i], varName.c_str(), data.L[i], data.U[i], 0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(scip, x[i]);
		varName = "y" + std::to_string(i + 1);
		SCIPcreateVarBasic(scip, &y[i], varName.c_str(), 0, 1, 0, SCIP_VARTYPE_BINARY);
		SCIPaddVar(scip, y[i]);
	}

	//sum of y_i <= M
	SCIP_CONS* ct_sum_y;
	SCIPcreateConsBasicLinear(scip, &ct_sum_y, "sum_y", 0, NULL, NULL, -SCIPinfinity(scip), param.M);
	for (int i = 0; i < data.m; ++i)
		SCIPaddCoefLinear(scip, ct_sum_y, y[i], 1.0);
	SCIPaddCons(scip, ct_sum_y);

	//sum of x_i*y_i <= W
	SCIP_EXPR** expr_x;
	expr_x = new SCIP_EXPR * [data.m];
	SCIP_EXPR** expr_y;
	expr_y = new SCIP_EXPR * [data.m];
	SCIP_EXPR** expr_xy;
	expr_xy = new SCIP_EXPR * [data.m];
	for (int i = 0; i < data.m; ++i) {
		SCIP_EXPR* xy[2];
		SCIPcreateExprVar(scip, &expr_x[i], x[i], NULL, NULL);
		SCIPcreateExprVar(scip, &expr_y[i], y[i], NULL, NULL);
		xy[0] = expr_x[i];
		xy[1] = expr_y[i];
		SCIPcreateExprProduct(scip, &expr_xy[i], 2, xy, 1.0, NULL, NULL);
	}
	SCIP_EXPR* sum_xy;
	SCIPcreateExprSum(scip, &sum_xy, data.m, expr_xy, NULL, 0.0, NULL, NULL);
	SCIP_CONS* ct_sum_xy;
	SCIPcreateConsBasicNonlinear(scip, &ct_sum_xy, "sum_xy", sum_xy, -SCIPinfinity(scip), data.C);
	SCIPaddCons(scip, ct_sum_xy);

	//create objective
	SCIP_EXPR* obj;
	SCIP_EXPR** frac;
	frac = new SCIP_EXPR * [data.T];
	for (int j = 0; j < data.T; ++j) {
		SCIP_EXPR* exprs[2];
		vector<SCIP_EXPR*> exp(data.m);
		vector<SCIP_EXPR*> power(data.m);
		SCIP_EXPR** product_xy;
		product_xy = new SCIP_EXPR * [data.m];
		SCIP_EXPR** product_y;
		product_y = new SCIP_EXPR * [data.m];
		for (int i = 0; i < data.m; ++i) {
			SCIPcreateExprSum(scip, &power[i], 1, &expr_x[i], &data.eta[j][i], data.kappa[j][i], NULL, NULL);
			SCIPcreateExprExp(scip, &exp[i], power[i], NULL, NULL);
			SCIP_EXPR* xye[2];
			xye[0] = expr_y[i];
			xye[1] = exp[i];
			SCIPcreateExprProduct(scip, &product_xy[i], 2, xye, data.a[j], NULL, NULL);
			SCIP_EXPR* ye[2];
			ye[0] = expr_y[i];
			ye[1] = exp[i];
			SCIPcreateExprProduct(scip, &product_y[i], 2, ye, 1, NULL, NULL);
		}
		SCIPcreateExprSum(scip, &exprs[0], data.m, product_xy, NULL, 0, NULL, NULL);
		SCIP_EXPR* denominator;
		SCIPcreateExprSum(scip, &denominator, data.m, product_y, NULL, 0, NULL, NULL);
		SCIP_EXPR* denominator1;
		SCIPcreateExprSum(scip, &denominator1, 1, &denominator, NULL, data.b[j], NULL, NULL);
		SCIPcreateExprPow(scip, &exprs[1], denominator1, -1.0, NULL, NULL);
		SCIPcreateExprProduct(scip, &frac[j], 2, exprs, 1.0, NULL, NULL);
	}
	SCIPcreateExprSum(scip, &obj, data.T, frac, NULL, 0.0, NULL, NULL);

	SCIP_VAR* z;
	SCIPcreateVarBasic(scip, &z, "z", 0, SCIPinfinity(scip), 1, SCIP_VARTYPE_CONTINUOUS);
	SCIPaddVar(scip, z);
	SCIP_EXPR* expr_z;
	SCIPcreateExprVar(scip, &expr_z, z, NULL, NULL);
	SCIP_EXPR* objective[2];
	objective[0] = obj;
	objective[1] = expr_z;
	SCIP_Real coefs[2];
	coefs[0] = 1;
	coefs[1] = -1;
	SCIP_EXPR* sum_obj;
	SCIPcreateExprSum(scip, &sum_obj, 2, objective, coefs, 0.0, NULL, NULL);
	SCIP_CONS* ct_obj;
	SCIPcreateConsBasicNonlinear(scip, &ct_obj, "obj", sum_obj, 0, 0);
	SCIPaddCons(scip, ct_obj);

	SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE);
	SCIPsetRealParam(scip, "limits/time", 3600);
	SCIPsetRealParam(scip, "limits/gap", 1e-4);
	SCIPsetRealParam(scip, "parallel/maxnthreads", 64);
	SCIPsetRealParam(scip, "parallel/minnthreads", 8);
	//SCIPwriteOrigProblem(scip, "model.lp", nullptr, FALSE);

	auto start = chrono::steady_clock::now(); //get start time
	auto end = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

	SCIPsolve(scip);

	end = chrono::steady_clock::now();
	elapsed_seconds = end - start;
	time_for_solve = elapsed_seconds.count();

	SCIP_SOL* bestsol = SCIPgetBestSol(scip);
	obj_val_scip = SCIPgetSolOrigObj(scip, bestsol);

	//get value x
	vector<double> ansX(data.m);
	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		ansX[iter_m] = SCIPgetSolVal(scip, bestsol, x[iter_m]);
	}

	//get value y
	vector<bool> ansY(data.m);
	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		if (SCIPgetSolVal(scip, bestsol, y[iter_m]) > 0.5) ansY[iter_m] = 1;
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

	gap = SCIPgetGap(scip);

	ofstream report_results(output, ofstream::out);
	report_results << data.m << "," << data.C << "," << param.M
		<< setprecision(8) << fixed << "," << obj_val_scip
		<< "," << obj_val_true << "," << time_for_solve << "," << gap << endl;
	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		report_results << "y[" << iter_m << "] = " << ansY[iter_m] << endl;
	}

	for (int iter_m = 0; iter_m < data.m; iter_m++) {
		report_results << "x[" << iter_m << "] = " << setprecision(5) << ansX[iter_m] << endl;
	}
	report_results.close();
}