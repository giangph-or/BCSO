#include <iostream>
#include <vector>
#include "Data.h"
#include "Param.h"
#include <string>
#include <fstream>
#include <chrono>
#include "scip/scip.h"
#include "scip/scip_expr.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_quadratic.h"

using namespace std;

class ScipSolverAssortment {
public:
    DataAssortment data;
    ParamAssortment param;
    double time_for_param = 0;
    double time_for_solve = 0;
    double gap = 0;
    double obj_val_scip = 0;
    double obj_val_true = 0;

    ScipSolverAssortment();

    ScipSolverAssortment(DataAssortment data, int K, double tol_lamda, int M);

    void solve(string output, int time_limit);
};

class ScipSolverFacility {
public:
    DataFacility data;
    ParamFacility param;
    double time_for_param = 0;
    double time_for_solve = 0;
    double gap = 0;
    double obj_val_cplex = 0;
    double obj_val_true = 0;

    ScipSolverFacility();

    ScipSolverFacility(DataFacility data, int K, double tol_lamda, int M);

    void solve(string output, int time_limit);
};
