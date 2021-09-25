#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include "tridiagmat.hpp"

using namespace arma; 
using namespace std; 

int main(int argc, char const* argv[]) {
	int N = 6; 
	double h; 
	h = 1.0/((double) N); 
	double h2 = h*h; 
	double diagonal = 2.0*h2; 
	double subsupdiagonal = -1.*h2; 
	TriDiagMat myInit; 
	mat A = myInit.initialize_symmetric_from_scalars(N, subsupdiagonal, diagonal);
	vec eigvals; 
	mat eigvecs;
	eig_sym(eigvals, eigvecs, A); 
	cout << "Started out with 6x6-matrix A: " << endl; 
	cout << A << endl; 
	cout << "Used arma::eig_sym to calculate eigenvalues and eigenvectors." << endl; 
	cout << "Eigenvalues of the matrix A:\n"; 
	for (int i=0; i<N; i++) {
		cout << eigvals(i) << endl; 
	} 
	cout << "Matrix containing eigenvectors of A:" << endl; 
	cout << eigvecs << endl; 
	
	return 0; 


}
