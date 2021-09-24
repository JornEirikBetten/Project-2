#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include "tridiagmat.hpp"
#include "jacobi.hpp"

using namespace arma; 
using namespace std; 


int main(int argc, char const *argv[]) 
{
	int N = atoi(argv[1]);
	string filename = "eigensolve_N_" + to_string(N) + ".csv";  
	double h = (double) 1./((double) N); 
	double h2 = h*h; 
	double diagonal = 2.*h2; 
	double subsupdiagonal = -1.*h2; 
	TriDiagMat myInit; 
	//mat A = myInit.initialize_symmetric_from_scalars(N, subsupdiagonal, diagonal);  
	mat A = mat(N,N).randn(); 
	A = symmatu(A); 
	vec eigval = vec(N); 
	mat eigvecs = mat(N,N, fill::eye); 
	vec eigval1; 
	mat eigvecs1;   
	eig_sym(eigval1, eigvecs1, A); 
	//cout << eigval1 << endl; 
	//cout << eigvecs1 << endl; 
	Jacobi mySolver; 
	double epsilon = pow(10.0, -8.0); 
	double max_offdiag; 
	int maximum_iterations = pow(10.0, 8.0); 
	int k, l; 
	int iterations = 0; 
	bool converged; 
	double mod = mySolver.max_offdiag_symmetric(A, k, l); 
	mySolver.jacobi_eigensolver(A, epsilon, eigval, eigvecs, maximum_iterations, iterations, converged);
	mat normalized_eigvecs = eigvecs/norm(eigvecs); 
	cout << norm(eigvecs) << endl; 
	mySolver.write_to_csv(normalized_eigvecs, eigval, filename); 
	//cout << eigval << endl; 
	cout << "The number of iterations (similarity transformations) it took to diagonalize the matrix with dimensionality " << to_string(N) << "x" << to_string(N) << " is " << to_string(iterations) << endl; 
	//sort(eigval.begin(), eigval.end()); 
	//cout << eigval << endl; 
	
	
	return 0; 
}
