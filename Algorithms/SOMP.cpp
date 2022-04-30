#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "Algorithm.cpp"

class SplineOMPAlgorithm : public Algorithm {
	public:
		SplineOMPAlgorithm(int nAtoms, int szSignal, int szTest) : Algorithm(nAtoms, szSignal, szTest) {}
		

		void RunAlgorithm(double* vSignal, double* rSignal, double* mDictionary, double* fullDictionary) override {
			double *mNewDictionary = new double[szSignal*nAtoms];
			double *mOrthogonalDictionary = new double[szSignal*szSignal];
			double *mBiorthogonal = new double[szSignal*szSignal];
			double *vCoefficients = new double[nAtoms];
		    int *iOldDictionary = new int[nAtoms];
		    int chosen;
		    double tolerance1 = 0.0001;
		    
			omp_set_num_threads(omp_get_num_procs());
		    
		    /* Perform deep copy of the original dictionary matrix */
		    for ( int i = 0; i < nAtoms; i++ )
		    {
		    	#pragma omp for
		        for ( int j = 0; j < szSignal; j++ )
		        {
		            mNewDictionary[(i*szSignal) + j] = mDictionary[(i*szSignal) + j];
		        }
		    }
		
		    chosen = (omp(vSignal,mNewDictionary,tolerance1,szSignal, nAtoms,
		        iOldDictionary,mOrthogonalDictionary,mBiorthogonal,vCoefficients));
		  	
			printf("%d\n", chosen);
		    printf("\nvCoefficients\n");    
		    for(int i=0;i<chosen;i++){
		    	printf("%f %d\n",vCoefficients[i], iOldDictionary[i]);
		    }
		    double aSignal[szTest];
		    for (int i=0;i<szTest;i++){
				double temp = 0;
				#pragma omp for
		    	for (int j=0;j<chosen;j++){
					temp += fullDictionary[iOldDictionary[j]*szTest+i]*vCoefficients[j];
		    	}
		    	aSignal[i] = temp;
//		    	std::cout << temp << " ";
			}
			
			printf("\naSignal\n");    
		    /*for(int i=0;i<szTest;i++){
		    	fout << aSignal[i] << ' ';
		    }
		    fout.close();*/
		    
		    double max=0.0;
		    #pragma omp for
		    for(int i=2;i<szTest-2;i++){
		    	//if(fabs(rSignal[i]-aSignal[i])>max){
		    	//	max=fabs(rSignal[i]-aSignal[i]);
		    	//	gde = i;
				//}
				max += (rSignal[i]-aSignal[i]) * (rSignal[i]-aSignal[i]);
			}
			std::cout << std::endl << std::endl << "Diff " << std::setprecision(10) << sqrt(max/szTest) << std::endl;
			
		
			
//			std::cout << std::endl << "rSignal" << std::endl;
//			for(int i=0;i<szTest;i++){
//				std::cout << rSignal[i] << ' ';
//			}
		  	
		  	
		  	//delete [] mNewDictionary;
			//delete [] mOrthogonalDictionary;
			//delete [] mBiorthogonal;
			//delete [] vCoefficients;
		    //delete [] iOldDictionary;
		    
			return;
		}
	
	
	private:
		double real_inner_product(double *v1, double *v2, int szVector)
		{
		  double sum = 0;
		  #pragma omp for
		  for ( int i = 0; i < szVector; i++)
		  {
		    sum += v1[i]*v2[i];    
		  }
		  
		return sum;
		}
		
		
		double vectorNorm(double *first, int N) {
			return sqrt(real_inner_product(first, first, N));
		}
		
		
		double* trans_multyplication(double *vector, double *matrix, double *returnVector, int m, int n)
		{
		    double t;
		    
		    for (int i=0;i<n;i++)
			{
				double temp = 0;
				#pragma omp for
		    	for (int j=0;j<m;j++)
		    	{
		      		temp += matrix[i*m+j]*vector[j];
		    	}
		    	returnVector[i] = temp;
		  	}
		    return returnVector;
		}
		
		void choose_atom(double *residue, double *mUnselectedAtoms, int m, int n, std::vector<int> &chosenAtoms)
		{
		    int leftIndex = 0, rightIndex = 0, maxIndex = 0;
		    double localMaxValue = 0, maxValue = 0, tolerance2 = 1e-4; // Should be global constant
		    double cc[n];
			double *cc1 = trans_multyplication(residue, mUnselectedAtoms, cc, m, n);
		    /* MATLAB CODE
		     * [max_c,q]=max(cc);
		     */
		    #pragma omp for
		    
		    
		    for ( int rightIndex = 0; rightIndex < n; rightIndex++)
		    {
		        if(fabs(cc1[rightIndex]) > localMaxValue){
	                maxIndex = rightIndex;
	                localMaxValue = fabs(cc1[rightIndex]);
	                if(rightIndex-leftIndex>2 && fabs(cc1[leftIndex]) > tolerance2){
	                	chosenAtoms.push_back(leftIndex);
	                	leftIndex = rightIndex;
	                	if(cc1[leftIndex]>maxValue)
							maxValue = cc1[leftIndex];
					}
	            }
		        else{
		        	if(rightIndex-leftIndex>2){
	                	chosenAtoms.push_back(maxIndex);
	                	leftIndex = rightIndex;
	                	maxIndex = rightIndex;
	                	localMaxValue = cc1[rightIndex];
	                	if(cc1[maxIndex]>maxValue)
							maxValue = cc1[maxIndex];
					}
				}
		    }
		    /* MATLAB CODE
		     * if max_c<tol2 
		     *   fprintf('%s stopped, max(|<f,q>|/||q||)<= tol2=%g.\n',name,tol2);
		     *   break;
		     * end
		     */
			
		    if (maxValue < tolerance2)
		    {
		    	chosenAtoms.clear();
		    }
		}
		
		
		void reorthognalize(double *mOrthogonalDictionary,int szSignal,int k, int nRepetitions)
		{
		    double alpha;
		
		     /* MATLAB CODE
		      * for l=1:nRepetitions
		      *     for p=1:k-1
		      *              Q(:,k)=Q(:,k)-(Q(:,p)'*Q(:,k))*Q(:,p);
		      *      end
		      * end
		     */
		    for ( int l = 0; l < nRepetitions; l++ )
		    {
		        if (k > 0)
		        {            
		            for ( int i = 0; i < k ; i++)
		            {
		                alpha = real_inner_product(&mOrthogonalDictionary[(i*szSignal)],&mOrthogonalDictionary[(k*szSignal)],szSignal);
		                #pragma omp for
		                for ( int j = 0; j < szSignal; j++)
		                {
		                    mOrthogonalDictionary[(k*szSignal) + j] = mOrthogonalDictionary[(k*szSignal)+j] - alpha*mOrthogonalDictionary[(i*szSignal) + j];
		                }
		            }
		        }
		    }
		}
		
		void calc_biorthogonal(double *mBiorthogonal, double *vNewAtom, double *vOrthogonalAtom, double normAtom, int szSignal, int k)
		{
		    /* Compute biorthogonal functions beta from 1 to k-1
		     * MATLAB CODE
		     * beta=beta-Q(:,k)*(D(:,k)'*beta)/nork;
		     */
		
		    double alpha;
		    
		    if ( k > 0 )
		    {
		        for ( int j = 0; j < k ; j++ )
		        {
		            alpha = real_inner_product(vNewAtom,&mBiorthogonal[(j*szSignal)],szSignal)/normAtom;
					#pragma omp for
		            for ( int i = 0; i < szSignal; i++ )
		            {		 	  
		                mBiorthogonal[(j*szSignal)+i] = mBiorthogonal[(j*szSignal)+i] - alpha*vOrthogonalAtom[i];        
		            }
		        }
		    }
		    
		    /* Calculate the kth biorthogonal function
		     * MATLAB CODE
		     * beta(:,k)=Q(:,k)/nork; % kth biorthogonal function
		     */        
		    #pragma omp for
		    for ( int i = 0; i< szSignal; i++ )
		    {
		        mBiorthogonal[(k*szSignal)+i] = vOrthogonalAtom[i]/normAtom;
		    } 
		            
		    
		}
		
		void calc_residue(double *residue, double *vSignal, double *mOrthogonalDictionary, int szSignal)
		{
		    /* Calculate the residue
		     * MATLAB CODE
		     * Re=Re-f*Q(:,k)*Q(:,k)';
		     */
		
		    double alpha;
		    
		    alpha = real_inner_product(mOrthogonalDictionary,vSignal,szSignal); 
			#pragma omp for   
		    for ( int i = 0; i < szSignal; i++)
		    {
		       residue[i] = residue[i] - alpha*mOrthogonalDictionary[i];
		    }
		    
		}
		
		void swap_elements(double *swapA, double *swapB, int nRows)
		{
		    double swappedElement;
		    #pragma omp for
		    for ( int i = 0; i < nRows; i++ )
		    {
		        swappedElement = swapB[i];
		        swapB[i] = swapA[i];
		        swapA[i] = swappedElement;
		    }
		}
		
		void copy_elements(double *array1, double *array2, int nRows)
		{
		    for ( int i = 0; i < nRows; i++ )
		    {
		        array1[i] = array2[i];
		    }
		}
		        
		    
		double normalize(double *atom, int szSignal)
		{   
		    double normAtom = sqrt(real_inner_product(atom,atom,szSignal));
		    #pragma omp for
		    for ( int i = 0; i < szSignal; i++)
		    {
		        atom[i] = atom[i]/normAtom;
		    }
		    return normAtom;
		}
		
		int omp(   
		            double *vSignal,
		            double *mNewDictionary,
		            double tolerance1,
		            int szSignal,
		            int nAtoms,
		            int *iOldDictionary,
		            double *mOrthogonalDictionary,
		            double *mBiorthogonal,
		            double *vCoefficients           
		        )
		{
		    /* Declare constants */
		    const double tolerance2 = 1e-10;
		    
		    
		    /* Declare variables */
		    int k = -1, i, nIterations, iNewAtom, old, cAtom;
		    std::vector<int> iChosenAtom;
		    double normKthAtom, normresidue, swap;
		    
		    
		    /* Declare dynamic arrays remember to delete after use */
		    double residue[szSignal];
		    
		    /* Perform deep copy of the original signal */
		    copy_elements(residue,vSignal,szSignal);
		
		   
			/* Populate the index array - remember MATLAB index's start at 1 not 0*/
			#pragma omp for
		    for ( int i = 0; i < nAtoms; i++ )
		    {
		        iOldDictionary[i] = i;
		    }    
		
		    
		    /* The number of iterations should be less the min of the nAtoms
		     * and the szSignal. */
		    
		
		    nIterations = std::min(nAtoms, szSignal);
			nIterations = 2;
		    
		
		    
		    /*************************************************************/
		    /********************* Start of main routine *****************/
		    /*************************************************************/
		    
			
		    for (int y = 0; y < nIterations; y++)
		    {
		    	k++;
		        iNewAtom = k*szSignal;
		        
		        /* Choose the next atom from the unselected atoms. */ 
		        choose_atom(residue,&mNewDictionary[iNewAtom],szSignal,(nAtoms - k), iChosenAtom);
		        
		        /* Stopping criterion (coefficient)
		         * MATLAB CODE
		         * if max_c<tol2
		         *      fprintf('%s stopped, max(|<f,q>|/||q||)<= tol2=%g.\n',name,tol2);
		         *      break;
		         * end
		         */
		        if ( iChosenAtom.size() == 0 )
		        {
		        	printf("No Good Atoms");
		            break;            
		        }
		        old = k;
		        k--;
		        for(size_t p=0;p<iChosenAtom.size();p++){
		        	k++;
			        iNewAtom = k*szSignal;
			        cAtom = iChosenAtom[p];
		        	
			        /* Add k
			         * because we want the index in the mNewDictionary array,
			         * not the index in the smaller array we pass to choose_atom.
			         */
			        cAtom += old;
			
			        
			        /* Populate the matrices with the new atoms
			         * MATLAB CODE
			         * if q~=k
			         *      Q(:,[k q])=Q(:,[q k]); % swap k-th and q-th columns
			         *      D(:,[k q])=D(:,[q k]);
			         *      Di([k q])=Di([q k]);
			         * end
			         */        
			        int temp = iOldDictionary[k];
					iOldDictionary[k] = iOldDictionary[cAtom];
					iOldDictionary[cAtom] = temp;
			        swap_elements(&mNewDictionary[iNewAtom],&mNewDictionary[(cAtom*szSignal)],szSignal);
			        auto result1 = std::find(begin(iChosenAtom), end(iChosenAtom), k);
			        if(result1 != end(iChosenAtom))
			        	*result1 = cAtom;
			        /* Pick the kth atom as we have swapped the chosen one */
			        copy_elements(&mOrthogonalDictionary[iNewAtom],&mNewDictionary[iNewAtom],szSignal);
			    }
			    
		        /* Re-orthogonalization of Q(:,k)  w.r.t Q(:,1),..., Q(:,k-1) */
		        reorthognalize(mOrthogonalDictionary,szSignal,k,2);
		                
		        /* Normalize atom Q(:,k)
		         * MATLAB CODE
		         * nork=norm(Q(:,k));
		         * Q(:,k)=Q(:,k)/nork; %normalization
		         */
				for(int j=old;j<=k;j++){
					normKthAtom = normalize(&mOrthogonalDictionary[j*szSignal],szSignal);
					
					/* Compute biorthogonal functions beta*/
		        	calc_biorthogonal(mBiorthogonal, &mNewDictionary[iNewAtom], &mOrthogonalDictionary[iNewAtom], normKthAtom, szSignal, j);
					
					/* Calculate the residue */
		        	calc_residue(residue, vSignal, &mOrthogonalDictionary[iNewAtom], szSignal);
				}   
		        /*for(int i=0;i<k;i++){
			    	for(int j=0;j<szSignal;j++){
			    		std::cout << mBiorthogonal[i*szSignal+j] << " ";
					}
					std::cout << std::endl;
				}*/
		        /* Calculate the norm of the residue */
		        normresidue = sqrt(real_inner_product(residue,residue,szSignal));
		        
		        /*for(int i=0;i<szSignal;i++){
		        	printf("%f ", residue[i]);
				}
				printf(" %f\n", normresidue); */
		        
		        /* Stopping criteria (distance)
		         * MATLAB CODE - We use the residue check with omp if this is the same
		         * if (norm(f'-D(:,1:k)*(f*beta)')*sqrt(delta) < tol) && (tol~= 0)break;end;
		         */
		         
		        if ( normresidue < tolerance1 && tolerance1 != 0 )
		        {
		            /* Break so will not increment k before exiting the loop */
		            k += 1;
		            printf("%d %f\n", k, real_inner_product(residue,residue,szSignal));
		            break;                          
		        }
		    	printf("%d %f\n", k, real_inner_product(residue,residue,szSignal));
		    	iChosenAtom.clear();
		    }
		    /*Calculate the coefficients
		    * MATLAB CODE
		    * c=f*beta;
		    */
		    #pragma omp for
		    for ( i = 0; i <k; i++ )
		    {
		    	/*printf("%f ",vSignal);*/
		    	/*printf("%f \n",&mBiorthogonal[(i*szSignal)]);*/
		        vCoefficients[i] = real_inner_product(vSignal,&mBiorthogonal[(i*szSignal)],szSignal);
		    }
		    
		    return k;
		}
};
