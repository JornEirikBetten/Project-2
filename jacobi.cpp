#include "jacobi.hpp"
#include <iostream>

using namespace arma; 
using namespace std;
 
 
 // Finding maximum offdiagonal element of matrix A
double Jacobi::max_offdiag_symmetric(const mat& A, int& k, int& l) {
	// finds size
	int N = A.n_rows; 
	// checks if square
	if (A.is_square()==true && N>1) {
		k = 0; 
		l = 1; 
		// First off diagonal element is set to maxval
		double maxval = abs(A(k,l)); 
		// looping through all elements in the upper triangle of A
		// to see if they are bigger than maxval. If bigger, change
		// k and l their new appropriate numbers. 
		for (int i=0; i<N; i++) {
			for (int j=i+1; j<N; j++) {
				if (abs(A(i,j)) > maxval) {
					maxval = A(i,j);
					k = i; 
					l = j; 
				}
			}
		}
	return maxval; 
	
	}
	else {
	 cout << "Matrix is either not square or N isn't larger than 1"; 
	}

}
// Rotation algorithm
void Jacobi::jacobi_rotate(mat& A, mat& R, int k, int l) {
	int N = A.n_rows;
	double tau, t, c, s; 
	double a_kk, a_kl, a_ll; 
	a_kk = A(k,k); 
	a_kl = A(k,l); 
	a_ll = A(l,l); 
	tau = (a_ll-a_kk)/2*a_kl; 
	if (tau>=0) {
		t = 1./(tau+sqrt(1+tau*tau)); 
	}
	else {
		t = -1./(-tau+sqrt(1+tau*tau)); 
	}
	c = 1./sqrt(1+t*t); 	// calculating cos(theta)
	s = c*t; 		// calculating sin(theta)
	
	// Updating A-matrix at diagonal and maximapoints
	A(k,k) = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s; 
	A(l,l) = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s; 
	A(k,l) = 0; 
	A(l,k) = 0; 
	for (int i = 0; i<N; i++) {
		if (i==k) {
			continue; 
		}
		if (i==l) {
			continue; 
		} 
		double a_ik,a_il; 
		a_ik = A(i,k); 
		a_il = A(i,l); 
		A(i,k) = a_ik - a_il;
		A(i,l) = a_il + a_ik; 
		A(k,i) = A(i,k); 
		A(l,i) = A(i,l);  
	}
	// Updating rotation matrix R
	for (int i=0; i<N; i++) {
		double r_ik, r_il; 
		r_ik = R(i,k); 
		r_il = R(i,l); 
		R(i,k) = r_ik*c - r_il*s; 
		R(i,l) = r_il*c + r_ik*s; 
	} 


}

void Jacobi::jacobi_eigensolver(mat& A, double eps, vec& eigenvalues, mat& eigenvectors, const int maxiter, int& iterations, bool& converged) {
	int k, l; 
	int N = A.n_rows; 
	double max_offdiag = max_offdiag_symmetric(A,k,l); 
	while (eps<max_offdiag) {
		if (iterations > maxiter) {
			cout << "Not converging before " << to_string(maxiter) << " iterations has passed"; 
			converged = false; 
			break; 
		}
		jacobi_rotate(A, eigenvectors, k, l);  
		max_offdiag = max_offdiag_symmetric(A, k, l);
		iterations += 1; 
	}
	converged = true; 
	for (int j=0; j<N; j++) {
		eigenvalues(j) = A(j,j); 
	}

}
