#include "Data.h"

DataAssortment::DataAssortment() {

}

void DataAssortment::read_data(string data_path) {
    fstream input;
    input.open(data_path);
    input >> T >> m >> W;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        double tmp;
        input >> tmp;
        a.push_back(tmp);
        input >> tmp;
        b.push_back(tmp);
    }

    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp_eta;
        vector<double> tmp_kappa;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            double tmp;
            input >> tmp;
            tmp_eta.push_back(tmp);
            input >> tmp;
            tmp_kappa.push_back(tmp);
        }
        eta.push_back(tmp_eta);
        kappa.push_back(tmp_kappa);
    }

    for (int iter_m = 0; iter_m < m; iter_m++) {
        double tmp;
        input >> tmp;
        L.push_back(tmp);
        input >> tmp;
        U.push_back(tmp);
    }
}

DataFacility::DataFacility() {

}

void DataFacility::read_data(string data_path) {
    fstream input;
    input.open(data_path);
    input >> T >> m >> C;
    for (int iter_t = 0; iter_t < T; iter_t++) {
        double tmp;
        input >> tmp;
        a.push_back(tmp);
        input >> tmp;
        b.push_back(tmp);
    }

    for (int iter_t = 0; iter_t < T; iter_t++) {
        vector<double> tmp_eta;
        vector<double> tmp_kappa;
        for (int iter_m = 0; iter_m < m; iter_m++) {
            double tmp;
            input >> tmp;
            tmp_eta.push_back(tmp);
            input >> tmp;
            tmp_kappa.push_back(tmp);
        }
        eta.push_back(tmp_eta);
        kappa.push_back(tmp_kappa);
    }

    for (int iter_m = 0; iter_m < m; iter_m++) {
        double tmp;
        input >> tmp;
        L.push_back(tmp);
        input >> tmp;
        U.push_back(tmp);
    }
}