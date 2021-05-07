#include <iostream>
#include "DataSet.cpp"
#include "LassoRegression.h"
#include "matrix.h"
#include "Algorithm.cpp"
#include <iomanip>
#include <vector>

class LassoAlgorithm : public Algorithm {
	public:
		LassoAlgorithm(int nAtoms, int szSignal, int szTest) : Algorithm(nAtoms, szSignal, szTest) {}
		

		void RunAlgorithm(double* vSignal, double* rSignal, double* mDictionary, double* fullDictionary) override {
			DataSet dataSet(vSignal, mDictionary, nAtoms, szSignal);
    		LassoRegression *lasso = new LassoRegression(dataSet.sample, dataSet.target);
    		
    		double *vCoefficients = lasso->cyclicalCoordinateDescent(0.001, 0.01);
    		
		    printf("\nvCoefficients\n");    
		    for(int i=0;i<nAtoms;i++){
		    	printf("%f\n",vCoefficients[i]);
		    }
    		
    		double aSignal[szTest];
		    for (int i=0;i<szTest;i++){
				double temp = 0;
		    	for (int j=0;j<nAtoms;j++){
					temp += fullDictionary[j*szTest+i]*vCoefficients[j];
		    	}
		    	aSignal[i] = temp;
			}
			
			printf("\naSignal\n");    
		    for(int i=0;i<szTest;i++){
		    	printf("%f ",aSignal[i]);
		    }
		    
		    double max=0.0;
		    for(int i=2;i<szTest-2;i++){
		    	if(fabs(rSignal[i]-aSignal[i])>max){
		    		max=fabs(rSignal[i]-aSignal[i]);
				}
			}
			std::cout << std::endl << std::endl << "Diff " << std::setprecision(10) << max << std::endl;
			
			std::cout << std::endl << "rSignal" << std::endl;
			for(int i=0;i<szTest;i++){
				std::cout << rSignal[i] << ' ';
			}
			
			return;
		}
};
