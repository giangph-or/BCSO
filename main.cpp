#include <string>
#include "Data.h"
#include "MISOCPSolver.h"
#include "MILPSolver.h"
#include "BCMISOCPSolver.h"
#include "BCMILPSolver.h"
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

	Data data;
	data.read_data(input);

	if (solver == "misocp") {
		string output = "out//misocp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		MISOCPSolver misocpSolver(data, K, tol_lamda, M);
		misocpSolver.solve(output, time_limit);
	}

	if (solver == "bcmisocp") {
		string output = "out//bcmisocp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		BCMISOCPSolver bcmisocpSolver(data, K, tol_lamda, M);
		bcmisocpSolver.solve(output, time_limit);
	}

	if (solver == "milp") {
		string output = "out//milp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		MILPSolver milpSolver(data, K, tol_lamda, M);
		milpSolver.solve(output, time_limit);
	}

	if (solver == "bcmilp") {
		string output = "out//bcmilp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		BCMILPSolver bcmilpSolver(data, K, tol_lamda, M);
		bcmilpSolver.solve(output, time_limit);
	}

	if (solver == "pwlGRB") {
		string output = "out//pwlGurobi//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
		ApproxGurobiSolver approxSolver(data, K, tol_lamda, M);
		approxSolver.solve(output, time_limit);
	}
}

