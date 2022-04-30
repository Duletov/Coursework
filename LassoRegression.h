//
// Created by user on 08.10.2018.
//

#ifndef LASSOREGRESSION_LASSOREGRESSION_H
#define LASSOREGRESSION_LASSOREGRESSION_H

#include <vector>

class LassoRegression {
public:
    double **features;
    double *weights;
    double *target;

    unsigned long numberOfSamples;

    unsigned long numberOfFeatures;

    LassoRegression(std::vector<std::vector<double>> samples, std::vector<double> target);

    void predictions(double *result);

    void ro(double *results);

    double coordinateDescentStep(int weightIdx, double alpha);


    double *cyclicalCoordinateDescent(double tolerance, double alpha);


    void dumpWeightsToFile();

private:

    double **featuresMatrix(std::vector<std::vector<double>> samples);

    double **normalizeFeatures(double **matrix);

    double **emptyMatrix();

    double *initialWeights();

    double *targetAsArray(std::vector<double> target);

    void feature(int featureIdx, double *result);
};


#endif //LASSOREGRESSION_LASSOREGRESSION_H
