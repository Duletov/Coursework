/*
 * Implements omp in c++, only accepts real arguments,
 * stops for a fixed tol, must have 3 inputs and five 
 * outputs.
 * MATLAB CALL
 *  [mNewDictionary,iOldDictionary,mOrthogonalDictionary,...
 *      mBiorthogonal,vCoefficients]=omp(vSignal,...
 *      mDictionary,tolerance1);
 */
 
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double real_inner_product(double *v1, double *v2, int szVector)
{
  double sum = 0;
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
    	for (int j=0;j<m;j++)
    	{
      		temp += matrix[i*m+j]*vector[j];
    	}
    	returnVector[i] = temp;
  	}
    return returnVector;
}

int choose_atom(double *residue, double *mUnselectedAtoms, int m, int n)
{
    int maxIndex = 0;
    double maxValue = 0, tolerance2 = 1e-10; // Should be global constant
    double cc[n];
	double *cc1 = trans_multyplication(residue, mUnselectedAtoms, cc, m, n);
    /* MATLAB CODE
     * [max_c,q]=max(cc);
     */
    
    for ( int i = 0; i < n; i++)
    {
         if ( abs(cc1[i]) > maxValue )
            {
                maxValue = abs(cc1[i]);
                maxIndex = i;
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
        maxIndex = -1;
    }

    /* Free the memory */

    return maxIndex;
  
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
    /* If picking more than 150 atoms may be quicker to use BLAS and the matrix vector multiplication */
    if ( k > 0 )
    {
        for ( int j = 0; j < k ; j++ )
        {
            alpha = real_inner_product(vNewAtom,&mBiorthogonal[(j*szSignal)],szSignal)/normAtom;

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
    for ( int i = 0; i < szSignal; i++)
    {
       residue[i] = residue[i] - alpha*mOrthogonalDictionary[i];
    }
    
}

void swap_elements(double *swapA, double *swapB, int nRows)
{
    double swappedElement;
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
    int k, i, iChosenAtom, nIterations, iNewAtom;
    double normKthAtom, normresidue, swap;
    
    
    /* Declare dynamic arrays remember to delete after use */
    double residue[szSignal];
    
    /* Perform deep copy of the original signal */
    copy_elements(residue,vSignal,szSignal);

   
	/* Populate the index array - remember MATLAB index's start at 1 not 0*/         
    for ( int i = 0; i < nAtoms; i++ )
    {
        iOldDictionary[i] = i;
    }    

    
    /* The number of iterations should be less the min of the nAtoms
     * and the szSignal. */
    

        nIterations = nAtoms;

    

    
    /*************************************************************/
    /********************* Start of main routine *****************/
    /*************************************************************/
    

    for ( k = 0; k < nIterations; k++)
    {    
        iNewAtom = k*szSignal;
        
        /* Choose the next atom from the unselected atoms. */ 
        iChosenAtom = choose_atom(residue,&mNewDictionary[iNewAtom],szSignal,(nAtoms - k));
        
        /* Stopping criterion (coefficient)
         * MATLAB CODE
         * if max_c<tol2
         *      fprintf('%s stopped, max(|<f,q>|/||q||)<= tol2=%g.\n',name,tol2);
         *      break;
         * end
         */
        if ( iChosenAtom == -1 )
        {
        	printf("No Good Atoms");
            break;            
        }
         
        /* Add k
         * because we want the index in the mNewDictionary array,
         * not the index in the smaller array we pass to choose_atom.
         */  
         iChosenAtom += k;

        
        /* Populate the matrices with the new atoms
         * MATLAB CODE
         * if q~=k
         *      Q(:,[k q])=Q(:,[q k]); % swap k-th and q-th columns
         *      D(:,[k q])=D(:,[q k]);
         *      Di([k q])=Di([q k]);
         * end
         */        
        int temp = iOldDictionary[k];
		iOldDictionary[k] = iOldDictionary[iChosenAtom];
		iOldDictionary[iChosenAtom] = temp;
        swap_elements(&mNewDictionary[iNewAtom],&mNewDictionary[(iChosenAtom*szSignal)],szSignal);
        /* Pick the kth atom as we have swapped the chosen one */
        copy_elements(&mOrthogonalDictionary[iNewAtom],&mNewDictionary[iNewAtom],szSignal);
        
        
        /* Re-orthogonalization of Q(:,k)  w.r.t Q(:,1),..., Q(:,k-1) */
        reorthognalize(mOrthogonalDictionary,szSignal,k,2);
                
        
        /* Normalize atom Q(:,k)
         * MATLAB CODE
         * nork=norm(Q(:,k));
         * Q(:,k)=Q(:,k)/nork; %normalization
         */      
        normKthAtom = normalize(&mOrthogonalDictionary[iNewAtom],szSignal);
         
        
        /* Compute biorthogonal functions beta*/
        calc_biorthogonal(mBiorthogonal, &mNewDictionary[iNewAtom], &mOrthogonalDictionary[iNewAtom], normKthAtom, szSignal, k);
         
        
        /* Calculate the residue */
        calc_residue(residue, vSignal, &mOrthogonalDictionary[iNewAtom], szSignal);
   
        
        /* Calculate the norm of the residue */
        normresidue = sqrt(real_inner_product(residue,residue,szSignal));
        
        for(int i=0;i<szSignal;i++){
        	printf("%f ", residue[i]);
		}
		printf(" %f\n", normresidue); 
        
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
    }
    
    /*Calculate the coefficients
     * MATLAB CODE
     * c=f*beta;
     */
    for ( i = 0; i <k; i++ )
    {
    	/*printf("%f ",vSignal);*/
    	/*printf("%f \n",&mBiorthogonal[(i*szSignal)]);*/
        vCoefficients[i] = real_inner_product(vSignal,&mBiorthogonal[(i*szSignal)],szSignal);
    }
    
    return k;
}



void initialization(double *vSignal, double *mDictionary, double tolerance1, int szSignal, int nAtoms, int szTest, double *rSignal)
{
    double mNewDictionary[szSignal*nAtoms], mOrthogonalDictionary[szSignal*szSignal], mBiorthogonal[szSignal*szSignal], vCoefficients[nAtoms];
    int iOldDictionary[nAtoms];
    int chosen;
    
    /* Perform deep copy of the original dictionary matrix */
    for ( int i = 0; i < nAtoms; i++ )
    {
        for ( int j = 0; j < szSignal; j++ )
        {
            mNewDictionary[(i*szSignal) + j] = mDictionary[(i*szSignal) + j];
        }
    }

    chosen = (omp(vSignal,mNewDictionary,tolerance1,szSignal, nAtoms,
        iOldDictionary,mOrthogonalDictionary,mBiorthogonal,vCoefficients));
        
    double fullDict[nAtoms*szTest];
	string fname = "d" + to_string(nAtoms) + to_string(szTest) + ".txt";
	ifstream in(fname);
	if (in.is_open()){
		for(int i=0;i<200;i++){
			in >> fullDict[i];
		}
	}
	in.close();	
    
    printf("\nDictFull\n");    
    for(int i=0;i<nAtoms;i++){
    	for(int j=0;j<szTest;j++){
    		printf("%0.5f ", fullDict[i*szTest+j]);
		}
		printf("\n");
    }
    printf("end\n");
    
	printf("%d\n", chosen);
    printf("\nvCoefficients\n");    
    for(int i=0;i<chosen;i++){
    	printf("%f %d\n",vCoefficients[i], iOldDictionary[i]);
    }
    
    double aSignal[szTest];
    for (int i=0;i<szTest;i++)
	{
		double temp = 0;
    	for (int j=0;j<chosen;j++)
    	{
      		temp += fullDict[iOldDictionary[j]*szTest+i]*vCoefficients[j];
    	}
    	aSignal[i] = temp;
  	}
  	
  	printf("\naSignal\n");    
    for(int i=0;i<szTest;i++){
    	printf("%f ",aSignal[i]);
    }
    
    double max=0.0;
    for(int i=2;i<szTest-2;i++){
    	if(abs(rSignal[i]-aSignal[i])>max){
    		max=abs(rSignal[i]-aSignal[i]);
		}
	}
	printf("\n\nDiff - %.15f\n", max); 
  	
	return;
}

int main(int argc, char *argv[]){
	int nAtoms, szSignal, szTest, i,j;
	nAtoms = stoi(argv[1]);
	szSignal = stoi(argv[2]);
	szTest = stoi(argv[3]);
	double vSignal[szSignal];
	double rSignal[szTest];
	double mCosines[szSignal*nAtoms];
	double tol = 0.0000001;
	
	for(int i=0;i<szSignal;i++){
		vSignal[i] = (i/double(szSignal))*(i/double(szSignal));
	}
	
	for(int i=0;i<szTest;i++){
		rSignal[i] = (i/double(szTest))*(i/double(szTest));
	}
	
	string fname = "d" + to_string(nAtoms) + to_string(szSignal) + ".txt";
	ifstream in(fname);
	if (in.is_open()){
		for(int i=0;i<szSignal*nAtoms;i++){
			in >> mCosines[i];
		}
	}
	in.close();
	
	printf("\nDict\n");    
    for(int i=0;i<nAtoms;i++){
    	for(int j=0;j<szSignal;j++){
    		printf("%0.5f ", mCosines[i*szSignal+j]);
		}
		printf("\n");
	}
    printf("end\n");
    
    
	initialization(vSignal, mCosines, tol, szSignal, nAtoms, szTest, rSignal);
	
	
	printf("\nrSignal\n");    
    for(int i=0;i<szTest;i++){
    	printf("%f ",rSignal[i]);
    }
}
