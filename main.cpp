#include <armadillo>
#include <iostream>
#include <cmath>
#include "tridiagmat.hpp"
#include "jacobi.hpp"

using namespace arma; 
using namespace std; 


int main(int argc, char const *argv[]) 
{
	int N = atoi(argv[1]); 
	double h = (double) 1./((double) N); 
	double h2 = h*h; 
	double diagonal = 2.*h2; 
	double subsupdiagonal = -1.*h2; 
	TriDiagMat myInit; 
	mat A = myInit.initialize_symmetric_from_scalars(N, subsupdiagonal, diagonal);  
	vec eigval = vec(N); 
	mat eigvecs = mat(N,N, fill::eye); 
	vec eigval1; 
	mat eigvecs1;   
	eig_sym(eigval1, eigvecs1, A); 
	cout << eigval1 << endl; 
	Jacobi mySolver; 
	double epsilon = pow(10.0, -8.0); 
	double max_offdiag; 
	int maximum_iterations = pow(10.0, 8.0); 
	int k, l; 
	int iterations = 0; 
	bool converged; 
	double mod = mySolver.max_offdiag_symmetric(A, k, l); 
	mySolver.jacobi_eigensolver(A, epsilon, eigval, eigvecs, maximum_iterations, iterations, converged); 
	cout << A << endl; 
	cout << eigval << endl; 
	cout << iterations << endl; 
	
	
	return 0; 
}
