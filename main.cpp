#include <string>
#include "Data.h"
#include "MISOCPSolver.h"
#include "MILPSolver.h"
#include "BCMISOCPSolver.h"
#include "BCMILPSolver.h"
#include "BilinearSolver.h"
#include "BilinearNew.h"
#include "BCBilinearSolver.h"
#include "ApproxByGurobi.h"

int main(int argc, char* argv[]) {
	string instance_name = argv[1];
	int id = atoi(argv[2]);
	string input = "Data_AP_2//" + instance_name + "_" + to_string(id) + ".txt";
	int M = atoi(argv[3]);
	int K = 25;
	int time_limit = 3600;
	double tol_lamda = 0.01;
	string solver = argv[4];

	DataAssortment data;
	data.read_data(input);

	if (solver == "misocp") {
		string output = "out//misocp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		MISOCPSolverAssortment misocpSolver(data, K, tol_lamda, M);
		misocpSolver.solve(output, time_limit);
	}

	if (solver == "bcmisocp") {
		string output = "out//bcmisocp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		BCMISOCPSolverAssortment bcmisocpSolver(data, K, tol_lamda, M);
		bcmisocpSolver.solve(output, time_limit);
	}

	if (solver == "milp") {
		string output = "out//milp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		MILPSolverAssortment milpSolver(data, K, tol_lamda, M);
		milpSolver.solve(output, time_limit);
	}

	if (solver == "bcmilp") {
		string output = "out//bcmilp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		BCMILPSolverAssortment bcmilpSolver(data, K, tol_lamda, M);
		bcmilpSolver.solve(output, time_limit);
	}

	if (solver == "bilinear") {
		string output = "out//bilinear//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		BilinearSolverAssortment bilinearSolver(data, K, tol_lamda, M);
		bilinearSolver.solve(output, time_limit);
	}

	if (solver == "bcbilinear") {
		string output = "out//bcbilinear//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		BCBilinearSolverAssortment bcbilinearSolver(data, K, tol_lamda, M);
		bcbilinearSolver.solve(output, time_limit);
	}

	if (solver == "binew") {
		string output = "out//binew//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		BilinearNewAssortment bilinearNew(data, K, tol_lamda, M);
		bilinearNew.solve(output, time_limit);
	}

	if (solver == "pwlGRB") {
		string output = "out//pwlGurobi//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		ApproxGurobiSolverAssortment approxSolver(data, K, tol_lamda, M);
		approxSolver.solve(output, time_limit);
	}
}

