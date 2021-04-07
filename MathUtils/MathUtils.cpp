#include <cmath>
#include "MathUtils.h"

double get_inner_product(double *vector1, double *vector2, int vectorSize) {
  
  double inner_product = 0;
  for (int i = 0; i < vectorSize; i++)
  {
    inner_product += vector1[i] * vector2[i];    
  }
  
  return inner_product;
}

double get_vector_norm(double *vector, int vector_size) {

	return sqrt(get_inner_product(vector, vector, vector_size));
}

void normalize_vector(double *vector, int vector_size) {   
    
	double vectorNorm = get_vector_norm(vector, vector_size);
    
    for (int i = 0; i < vector_size; i++) {
        vector[i] = 1.0 * vector[i] / vectorNorm;
    }
}

double* multiply(double *vector, double *matrix, int matrix_size1, int matrix_size2, double *result_vector)
{
    for (int i = 0; i < matrix_size2; i++) {
		
		double element = 0;
    	for (int j = 0; j < matrix_size1; j++) {
      		element += matrix[i * matrix_size1 + j] * vector[j];
    	}
    	result_vector[i] = element;
  	}
  	
    return result_vector;
}
