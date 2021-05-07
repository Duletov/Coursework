//
// Created by user on 05.10.2018.
//
#define _USE_MATH_DEFINES

#include <chrono>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

class DataSet{
	public:
		std::vector<std::vector<double>> sample;
		std::vector<double> target;
		
		DataSet(double* vSignal, double* mDictionary, int nAtoms, int szSignal) {
			std::vector<double> tmp;
			double t;
		    for (int i = 0; i < szSignal; i++) {
		        target.push_back(vSignal[i]);
		    }
			for(int i=0;i<nAtoms;i++){
				for(int j=0;j<szSignal;j++){
					tmp.push_back(mDictionary[szSignal*i+j]);
				}
				sample.push_back(tmp);
				tmp.clear();
			}
		}
};
