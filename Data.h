#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <iomanip>

using namespace std;

class DataAssortment {
public:
    //Assortment + Price
    int T;//#customers
    int m;//#products
    double W;
    vector<double> a;
    vector<double> b;
    vector<vector<double>> eta;
    vector<vector<double>> kappa;
    vector<double> L;
    vector<double> U;

    DataAssortment();
    void read_data(string data_path);
};

class DataFacility {
public:
    //Facility Location + Cost
    int T;//#zones
    int m;//#locations
    double C;//sum y*x <= C
    vector<double> a;//a=competitors' utilities in (0,1)
    vector<double> b;//b=1
    vector<vector<double>> eta;
    vector<vector<double>> kappa;
    vector<double> L;
    vector<double> U;

    DataFacility();
    void read_data(string data_path);
};
#pragma once