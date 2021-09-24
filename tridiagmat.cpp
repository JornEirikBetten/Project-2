#include "tridiagmat.hpp"



mat TriDiagMat::initialize_from_scalars(int N, double a, double d, double e) {
	mat A = mat(N, N, fill::eye)*d; 
	// fill the first row (row index 0)
	A(0,1) = e; 
	
	// loop that fills rows 2 to n-1 (row indices 1 to n-2)
	for (int i=1; i<N-1; i++) {
		A(i, i+1) = e; 
		A(i, i-1) = a; 
	}
	// fill last row (row index n-1) 
	A(N-1, N-2) = a;  

	return A; 
}


mat TriDiagMat::initialize_symmetric_from_scalars(int N, double a, double d) {
	return initialize_from_scalars(N, a, d, a); 
}

mat TriDiagMat::initialize_from_vectors(int N, const vec& a, const vec& d, const vec& e) {
	// start from identity matrix
	mat A = mat(N,N, fill::eye)*d; 
	//fill first row
	A(0,1) = a(0); 
	//loop that fills rows 2 to n-1 (row indices 1 to n-2)
	for (int i=1; i<N-1; i++) {
		A(i, i-1) = a(i); 
		A(i, i+1) = e(i); 
	}
	// fill last row
	A(N-1, N-2) = e(N-2); 
	// returning matrix
	return A; 

}

mat TriDiagMat::initialize_dense_from_scalars(int N, double a, double d) {
	mat A = ones(N, N)*a; 
	for (int i=0; i<N; i++) {
		A(i,i) = d; 
	}
	return A; 
}

