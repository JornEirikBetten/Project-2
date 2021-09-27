#include <iostream>
#include <armadillo>
#include "jacobi.hpp"
#include "tridiagmat.hpp"


using namespace std; 
using namespace arma; 


int main(int argc, char const* argv[]) {

	int N = 4; 
	mat A = mat(N,N, fill::eye); 
	A(0, 3) = 0.5; 
	A(1, 2) = -0.7; 
	A(2,1) = -0.7; 
	A(3,0) = 0.5; 
	cout << "Matrix A looks like: " << endl; 
	cout << A << endl; 
	Jacobi solver; 
	double max_value; 
	int k; 
	int l; 
	max_value = solver.max_offdiag_symmetric(A, k, l); 
	cout << "The largest non-diagonal element of A is element (1, 2)," << endl; 
	cout << "which we want to be equal to (k,l), and has value 0.7. " << endl; 
	cout << "Max-value given from max_offdiag_symmetric: " << to_string(max_value) << endl; 
	cout << "Element k updated from max_offdiag_symmetric: " << to_string(k) << endl; 
	cout << "Element l updated from max_offdiag_symmetric: " << to_string(l) << endl; 
	cout << "If max-value = 0.7, k=1 and l=2, we are good, and function works in this case." << endl; 
}
