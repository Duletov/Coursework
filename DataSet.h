//
// Created by user on 05.10.2018.
//

#ifndef LASSOREGRESSION_DATASET_H
#define LASSOREGRESSION_DATASET_H

#include <vector>

class DataSet {
public:
    std::vector<std::vector<double>> sample;
    std::vector<double> target;

    DataSet(double* vSignal, double* mDictionary, int nAtoms, int szSignal);
};


#endif //LASSOREGRESSION_DATASET_H
