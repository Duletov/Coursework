#include <random>
#include <set>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Algorithm.cpp"
# include "qr_solve.hpp"

class StOMPAlgorithm : public Algorithm {
	public:
		StOMPAlgorithm(int nAtoms, int szSignal, int szTest) : Algorithm(nAtoms, szSignal, szTest) {}
		
		void RunAlgorithm(double* vSignal, double* rSignal, double* mDictionary, double* fullDictionary) override {
			double tolerance = 0.0001, check;
			double residue[szSignal];
			for(int i=0;i<szSignal;i++)
				residue[i] = vSignal[i];
			double *M = new double[nAtoms];
			bool indices[nAtoms];
			for(int i=0;i<nAtoms;i++){
				indices[i] = false;
			}
			int StOmp_sz = 0, iteration = 0;
			double *s_StOMP;
			check = sqrt(real_inner_product(residue,residue,szSignal));
			//printf("%f ",check);
			trans_multyplication(residue, mDictionary, szSignal, nAtoms, M);
			/*for(int i=0;i<szSignal;i++)
				std::cout << residue[i] << ' ';
			std::cout << std::endl;
			
			for(int i=0;i<nAtoms;i++){
					for(int j=0;j<szSignal;j++){
					std::cout << mDictionary[i*szSignal+j] << ' ';
				}
				std::cout << std::endl;
			}
			
			std::cout << std::endl;
			for(int i=0;i<nAtoms;i++)
				std::cout << M[i] << ' ';
			std::cout << std::endl;*/
			double random = ((double) std::rand() / (RAND_MAX));
			double lambda = (1.0 + random) * check / sqrt(szSignal);
			//std::cout << lambda << std::endl;
			for(int i=0;i<nAtoms;i++){
					if(indices[i] == false)
						StOmp_sz++;
					indices[i] = true;
			}
			double *tmp_darr = new double[StOmp_sz*szSignal];
			int k=0;
			for(int i=0;i<nAtoms;i++){
				if(indices[i]){
					for(int j=0;j<szSignal;j++){
					tmp_darr[k*szSignal+j] = mDictionary[i*szSignal+j];
				}
				k++;
				}
			}
			s_StOMP = svd_solve(StOmp_sz,szSignal,tmp_darr,vSignal);
			for (int i=0;i<szSignal;i++){
				double temp = 0.0;
	    		for (int j=0;j<nAtoms;j++){
	    			temp += mDictionary[j*szSignal+i]*s_StOMP[j];
	    		}
	    		residue[i] = vSignal[i] - temp;
			}
			std::cout << "OK";
			delete [] tmp_darr;
			check = sqrt(real_inner_product(residue,residue,szSignal));
			printf("%f ",check);
			iteration++;
			delete [] M;
			
			printf("\nvCoefficients\n");    
		    for(int i=0;i<nAtoms;i++){
		    	//printf("%f\n",s_StOMP[i]);
		    }
		    
		    double aSignal[szTest];
		    for (int i=0;i<szTest;i++){
				double temp = 0;
		    	for (int j=0;j<nAtoms;j++){
					temp += fullDictionary[j*szTest+i]*s_StOMP[j];
		    	}
		    	aSignal[i] = temp;
			}
		  	
			printf("\naSignal\n");    
		    for(int i=0;i<szTest;i++){
		    	//printf("%f ",aSignal[i]);
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
				//std::cout << rSignal[i] << ' ';
			}
		}
	
	
	
	private:
		double real_inner_product(double *v1, double *v2, int szVector)
		{
		  double sum = 0;
		  for ( int i = 0; i < szVector; i++)
		  {
		    sum += v1[i]*v2[i];    
		  }
		  
		return sum;
		}
		
		void trans_multyplication(double *vector, double *matrix, int m, int n, double *returnVector)
		{
		    double t;
		    for (int i=0;i<n;i++)
			{
				double temp = 0.0;
		    	for (int j=0;j<m;j++)
		    	{
		      		temp += fabs(matrix[i*m+j]*vector[j]);
		    	}
		    	returnVector[i] = temp;
		  	}
		}
};
