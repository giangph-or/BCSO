#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <iomanip>

using namespace std;

class Data {
public:
    int T;
    int m;
    double W;
    vector<double> a;
    vector<double> b;
    vector<vector<double>> eta;
    vector<vector<double>> kappa;
    vector<double> L;
    vector<double> U;
    Data();
    void read_data(string data_path);
};
#pragma once