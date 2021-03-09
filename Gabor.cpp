#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>

double real_inner_product(double *v1, double *v2, int szVector)
{
  double sum = 0;
  for ( int i = 0; i < szVector; i++)
  {
    sum += v1[i]*v2[i];    
  }
  
  return sum;
}

int main(){
	int i, j, N, M, scale, timeShift;
	double frequency, phase, normAtom, tmp, n, m;
	N = 10;
	M = 100;
	n = double(N);
	m = double(M);
	scale = 2;
	timeShift = 0;
	frequency = 0.01;
	phase = 0;
	
	double atom[N][M];
	for(i = 0; i<N; i++){
		for(j = 0; j<M; j++){
			atom[i][j] = (1/sqrt(scale)) * exp(-M_PI*((j/m-i/n)*(j/m-i/n))/(scale*scale)) * cos(2*M_PI*frequency*(j/M-i/N)+phase);
			printf("%f\n", atom[i][j]);
		}
	}
	for(i = 0; i<N; i++){
		for(j = 0; j<M; j++){
    		printf("%f\n", atom[i][j]);
		}
	}
	for(j = 0; j<N; j++){
		normAtom = sqrt(real_inner_product(atom[j],atom[j],M));
		for ( i = 0; i < M; i++)
    	{
        	atom[j][i] = atom[j][i]/normAtom;
    	}
	}
    
}
