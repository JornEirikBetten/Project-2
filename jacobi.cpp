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
				if (abs(A(i,j)) >= maxval) {
					maxval = abs(A(i,j));
					k = i; 
					l = j; 
				}
			}
		}
	return maxval; 
	
	}
	else {
	 cout << "Matrix is either not square or N isn't larger than 1";
	 return 0;  
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
	if (A(k,l) != 0) {
		tau = (a_ll-a_kk)/(2.0*a_kl); 
		if (tau>0) {
			t = 1./(tau+sqrt(1+tau*tau)); 
		}
		else {
			t = -1./(-tau+sqrt(1+tau*tau)); 
			
		}
		c = 1./sqrt(1+t*t); 	// calculating cos(theta)
		s = c*t; 
	}
	else {
		c = 1.0; 	// calculating cos(theta)
		s = 0.0; 	
	}
	
	
	// Updating A-matrix at diagonal and maximapoints
	A(k,k) = a_kk*c*c - 2.0*a_kl*c*s + a_ll*s*s; 
	A(l,l) = a_ll*c*c + 2.0*a_kl*c*s + a_kk*s*s; 
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
		A(i,k) = a_ik*c - a_il*s;
		A(i,l) = a_il*c + a_ik*s; 
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

void Jacobi::swap_columns(mat& A, int col1, int col2) {
	int N = A.n_rows; 
	double temp; 
	if (col1<N && col2<N && col1 != col2) {
		for (int i=0; i<N; i++) {
			temp = A(i, col1); 
			A(i,col1) = A(i,col2);
			A(i,col2) = temp;  
		}
	
	}

}

void Jacobi::sort_cols(mat& A, vec& eigenvalues) {
	int N = A.n_rows; 
	double val; 
	int count;  
	for (int i=0; i<N; i++) {
		count = 0; 
		for (int j=0; j<N; j++) {
			if (j==i) {
				continue; 
			}
			if (eigenvalues(i) > eigenvalues(j)) {
				count += 1; 
			}
		}
		swap_columns(A, i, count); 
		
	
	}
	sort(eigenvalues.begin(), eigenvalues.end()); 


}

void Jacobi::write_to_csv(mat& eigvec, vec& eigval, string filename) {
	int N = eigvec.n_rows; 
	ofstream file; 
	file.open(filename); 
	file << "x,"; 
	file.precision(4); 
	file.scientific; 
	for (int i=0; i<N; i++) {
		file << eigval(i) << ",";
	}
	file << "\n"; 
	for (int i=0; i<N; i++) {
		file << "0.0000e+0,"; 
	}
	file << "0.0000e+0\n"; 
	for (int i=0; i<N; i++) {
		file << float(i+1)/(float(N)+1.0) << ","; 
		for (int j=0; j<N; j++) {
			file << eigvec(i,j) << ","; 
		}
		file << "\n"; 
	
	}
	file << "1.0000e+0, "; 
	for (int i=0; i<N-1; i++) {
		file << "0.0000e+0,"; 
	}
	file << "0.0000e+0\n"; 

}

