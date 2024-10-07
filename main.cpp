#include <string>
#include "Data.h"
#include "MISOCPSolver.h"
#include "MILPSolver.h"
#include "OAMILPSolver.h"
#include "Bilinear.h"
#include "ApproxByGurobi.h"
#include "ScipSolver.h"

int main(int argc, char* argv[]) {
	string problem = argv[1];
	//===========Assortment + Price===========//
	if (problem == "AP") {
		string instance_name = argv[2];
		int id = atoi(argv[3]);
		string input = "Data_AP//" + instance_name + "_" + to_string(id) + ".txt";
		int M = atoi(argv[4]);
		int K = 25;
		int time_limit = 3600;
		double tol_lamda = 0.01;
		string solver = argv[5];

		DataAssortment data;
		data.read_data(input);

		if (solver == "misocp") {
			string output = "out_AP//misocp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			MISOCPSolverAssortment misocpSolver(data, K, tol_lamda, M);
			misocpSolver.solve(output, time_limit);
		}

		if (solver == "milp") {
			string output = "out_AP//milp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			MILPSolverAssortment milpSolver(data, K, tol_lamda, M);
			milpSolver.solve(output, time_limit);
		}

		if (solver == "binew") {
			string output = "out_AP//binew//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			BilinearNewAssortment bilinearNew(data, K, tol_lamda, M);
			bilinearNew.solve(output, time_limit);
		}

		if (solver == "pwlGRB") {
			string output = "out_AP//pwlGurobi//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			ApproxGurobiSolverAssortment approxSolver(data, K, tol_lamda, M);
			approxSolver.solve(output, time_limit);
		}

		if (solver == "scip") {
			string output = "out_AP//scip//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			ScipSolverAssortment scipSolver(data, K, tol_lamda, M);
			scipSolver.solve(output, time_limit);
		}
	}

	//===========Facility Location + Cost===========//
	if (problem == "FLC") {
		string instance_name = argv[2];
		int id = atoi(argv[3]);
		string input = "Data_FLC//" + instance_name + "_" + to_string(id) + ".txt";
		int M = atoi(argv[4]);
		int K = 25;
		int time_limit = 3600;
		double tol_lamda = 0.01;
		string solver = argv[5];

		DataFacility data;
		data.read_data(input);

		if (solver == "binew") {
			string output = "out_FLC//binew//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			BilinearNewFacility bilinearNew(data, K, tol_lamda, M);
			bilinearNew.solve(output, time_limit);
		}

		if (solver == "pwlGRB") {
			string output = "out_FLC//pwlGurobi//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			ApproxGurobiSolverFacility approxSolver(data, K, tol_lamda, M);
			approxSolver.solve(output, time_limit);
		}

		if (solver == "milp") {
			string output = "out_FLC//milp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			MILPSolverFacility milpSolver(data, K, tol_lamda, M);
			milpSolver.solve(output, time_limit);
		}

		if (solver == "oa") {
			string output = "out_FLC//oa//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			OASolverFacility oaSolver(data, K, tol_lamda, M);
			oaSolver.solve(output, time_limit);
		}

		if (solver == "misocp") {
			string output = "out_FLC//misocp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			MISOCPSolverFacility misocpSolver(data, K, tol_lamda, M);
			misocpSolver.solve(output, time_limit);
		}

		if (solver == "scip") {
			string output = "out_FLC//scip//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
			ScipSolverFacility scipSolver(data, K, tol_lamda, M);
			scipSolver.solve(output, time_limit);
		}
	}
}

