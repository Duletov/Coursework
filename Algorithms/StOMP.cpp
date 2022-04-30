#include <random>
#include <set>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "Algorithm.cpp"
# include "qr_solve.hpp"

class StOMPAlgorithm : public Algorithm {
	public:
		StOMPAlgorithm(int nAtoms, int szSignal, int szTest) : Algorithm(nAtoms, szSignal, szTest) {}
		
		void RunAlgorithm(double* vSignal, double* rSignal, double* mDictionary, double* fullDictionary) override {
			double tolerance = 0.1, check;
			double residue[szSignal];
			for(int i=0;i<szSignal;i++)
				residue[i] = vSignal[i];
			double *M = new double[nAtoms];
			bool indices[nAtoms];
			for(int i=0;i<nAtoms;i++){
				indices[i] = false;
			}
			int StOmp_sz = 0, iteration = 0, k = 0;
			double *s_StOMP;
			check = sqrt(real_inner_product(residue,residue,szSignal));
			while(iteration<2 && check>tolerance){
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
				*/
				std::cout << std::endl;
				for(int i=0;i<nAtoms;i++)
					std::cout << M[i] << ' ';
				std::cout << std::endl;
				double random = rand() / RAND_MAX;
				double lambda = (2.0 + random) * check / sqrt(szSignal);
				std::cout << lambda << ' ';
				for(int i=0;i<nAtoms;i++){
						if(indices[i] == false && M[i]>=lambda){
							StOmp_sz++;
							indices[i] = true;
							//std::cout << i << ' ';
						}
				}
				//std::cout << std::endl;
				double *tmp_darr = new double[StOmp_sz*szSignal];
				k=0;
				for(int i=0;i<nAtoms;i++){
					if(indices[i]){
						for(int j=0;j<szSignal;j++){
							tmp_darr[k+j*StOmp_sz] = mDictionary[i*szSignal+j];
						}
					k++;
					}
				}
				s_StOMP = svd_solve(szSignal, StOmp_sz, tmp_darr, vSignal);
				
				std::cout << StOmp_sz << ' ';
				for (int i=0;i<szSignal;i++){
					double temp = 0.0;
					k=0;
		    		for (int j=0;j<nAtoms;j++){
		    			if(indices[j]){
		    				temp += mDictionary[j*szSignal+i]*s_StOMP[k];
		    				k++;
						}
		    		}
		    		residue[i] = vSignal[i] - temp;
				}
				delete [] tmp_darr;
				check = sqrt(real_inner_product(residue,residue,szSignal));
				std::cout << check  << std::endl;
				iteration++;
			}
			delete [] M;
			printf("\nvCoefficients\n"); 
			k = 0;
			
		    for(int i=0;i<nAtoms;i++){
		    	if(indices[i]){
		    		std::cout << i << ' ' << s_StOMP[k] << std::endl;
		    		k++;
				}
		    }
		    
		    double aSignal[szTest];
		    for (int i=0;i<szTest;i++){
				double temp = 0;
		    	k=0;
	    		for (int j=0;j<nAtoms;j++){
	    			if(indices[j]){
	    				temp += fullDictionary[j*szTest+i]*s_StOMP[k];
	    				k++;
					}
				}
		    	aSignal[i] = temp;
	    	}
			/*
		  	std::ofstream fout;
			fout.open("imggg.txt");
			*/
			printf("\naSignal\n");    
		    for(int i=0;i<szTest;i++){
		    	//std::cout << aSignal[i] << ' ';
		    }
		    //fout.close();
		    
		    double max=0.0;
		    for(int i=2;i<szTest-2;i++){
		    	//if(fabs(rSignal[i]-aSignal[i])>max){
		    	//	max=fabs(rSignal[i]-aSignal[i]);
		    	//	gde = i;
				//}
				max += (rSignal[i]-aSignal[i]) * (rSignal[i]-aSignal[i]);
			}
			std::cout << std::endl << std::endl << "Diff " << std::setprecision(10) << sqrt(max/szTest) << std::endl;
			
			
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
