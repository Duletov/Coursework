class MathUtils
{
	public:
		
		void normalize_vector(double *vector, int vector_size);
		double get_vector_norm(double *vector, int vector_size);
		
		double* multiply(double *vector, double *matrix, int matrix_size1, int matrix_size2, double *result_vector);

		double get_inner_product(double *vector1, double *vector2, int vectorSize);
};
