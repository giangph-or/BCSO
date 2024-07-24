#include <iostream>
#include <vector>
#include <gurobi_c++.h>
#include "Data.h"
#include "Param.h"
#include "McCormick.h"
#include <string>
#include <fstream>
#include <chrono>

using namespace std;

class MISOCPSolver {
public:
    Data data;
    Param param;
    McCormick mcCormick;
    double time_for_param = 0;
    double time_for_mc = 0;
    double time_for_solve = 0;
    double gap = 0;
    double obj_val_cplex = 0;
    double obj_val_true = 0;

    MISOCPSolver();

    MISOCPSolver(Data data, int K, double tol_lamda, int M);

    void solve(string output, int time_limit);
};