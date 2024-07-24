#include <string>
#include "Data.h"
#include "MISOCPSolver.h"

int main(int argc, char* argv[]) {
	string instance_name = argv[1];
	int id = atoi(argv[2]);
	string input = "Data_AP_2//" + instance_name + "_" + to_string(id) + ".txt";
	int M = atoi(argv[3]);
	int K = 25;
	string output = "out//misocp//" + instance_name + "_" + to_string(id) + "_" + to_string(M) + "_" + to_string(K) + ".txt";
	int time_limit = 3600;
	double tol_lamda = 0.01;

	Data data;
	data.read_data(input);

	MISOCPSolver misocpSolver(data, K, tol_lamda, M);
	misocpSolver.solve(output, time_limit);
}

