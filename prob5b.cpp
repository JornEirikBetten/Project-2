#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include "tridiagmat.hpp"
#include "jacobi.hpp"

using namespace arma; 
using namespace std; 

vec analytical_eigenvalues(int N) {
	vec eigval = vec(N); 
	double h = 1.0/((double) N); 
	double h2 = h*h; 
	for (int j=1; j<(N+1); j++) {
		eigval(j-1) = 2.0*h2*(1-cos((double) j*M_PI/((double) N +1.0))); 
	}
	return eigval; 

}

mat analytical_eigenvectors(int N) {
	mat A = mat(N,N); 
	for (int j=1; j<(N+1); j++) {
		for (int i=1; i<(N+1); i++) {
			double val; 
			val = sin((double) j*(double) i*M_PI/((double) N + 1.0)); 
			A(j-1, i-1) = val; 
		}
	}
	return A; 
	

}

int main(int argc, char const* argv[]) {
	int N = 6; //which is exactly how much the size of the matrix increases b
	double h; 
	h = 1.0/((double) N); 
	double h2 = h*h; 
	double diagonal = 2.0*h2; 
	double subsupdiagonal = -1.*h2; 
	TriDiagMat myInit; 
	// intitializing tridiagonal matrix. Also making matrices and vectors
	// that will contain analytical and calculated solutions. 
	mat A = myInit.initialize_symmetric_from_scalars(N, subsupdiagonal, diagonal); 
	cout << "Started out with 6x6 matrix A: " << endl; 
	cout << A << endl; 
	vec calculated_eigvals = vec(N); 
	mat calculated_eigvecs = mat(N,N, fill::eye);
	vec analytical_eigvals; 
	mat analytical_eigvecs;
	analytical_eigvals = analytical_eigenvalues(N); // Analytical eigenvalues
	analytical_eigvecs = analytical_eigenvectors(N); // Analytical eigenvectors
	// Using the eigensolver using similarity transformations
	Jacobi mySolver; // Initiating Jacobi class
	double epsilon = pow(10.0, -8.0); // defining tolerance
	int maximum_iterations = pow(10.0, 8.0); // maximum iterations, in case it won't converge
	int iterations = 0; // defining and setting the number of iterations to 0
	bool converged;  // bool set to true when it has converged
	mySolver.jacobi_eigensolver(A, epsilon, calculated_eigvals, calculated_eigvecs, maximum_iterations, iterations, converged); // solving eigenvalue-system, using jacobi_eigensolver (Jacobi.cpp)
	//Ordering calculated eigenvalues and eigenvectors so that the eigenvalues and eigenvectors are ordered in an ascending order (eigenvaluewise) 
	double templ1 = calculated_eigvals(0); 
	double templ4 = calculated_eigvals(1); 
	double templ0 = calculated_eigvals(2); 
	double templ5 = calculated_eigvals(3); 
	double templ2 = calculated_eigvals(4); 
	double templ3 = calculated_eigvals(5); 
	calculated_eigvals(0) = templ0; 
	calculated_eigvals(1) = templ1; 
	calculated_eigvals(2) = templ2; 
	calculated_eigvals(3) = templ3; 
	calculated_eigvals(4) = templ4; 
	calculated_eigvals(5) = templ5; 
	vec tempev0 = vec(N); 
	vec tempev1 = vec(N); 
	vec tempev2 = vec(N); 
	vec tempev3 = vec(N); 
	vec tempev4 = vec(N); 
	vec tempev5 = vec(N); 

	for (int i=0; i<N; i++) {
		tempev1(i) = calculated_eigvecs(i,0); 
		tempev4(i) = calculated_eigvecs(i,1); 
		tempev0(i) = calculated_eigvecs(i,2);
		tempev5(i) = calculated_eigvecs(i,3); 
		tempev2(i) = calculated_eigvecs(i,4); 
		tempev3(i) = calculated_eigvecs(i,5);  
	}
	for (int i=0; i<N;i++) {
		calculated_eigvecs(i,0) = tempev0(i); 
		calculated_eigvecs(i,1) = tempev1(i); 
		calculated_eigvecs(i,2) = tempev2(i); 
		calculated_eigvecs(i,3) = tempev3(i); 
		calculated_eigvecs(i,4) = tempev4(i); 
		calculated_eigvecs(i,5) = tempev5(i); 
	
	}
	// Printing output
	cout << "After all rotations, 6x6-matrix SAS^T: " << endl; 
	cout << A << endl; 
	cout << "Used jacobi_eigensolver to calculate eigenvalues and eigenvectors." << endl; 
	cout << "Calculated eigenvalues of the matrix A:\n"; 
	for (int i=0; i<N; i++) {
		cout << calculated_eigvals(i) << endl; 
	} 
	cout << "Matrix containing calculated eigenvectors of A:" << endl; 
	// outputting scaled eigenvectors, so that they are of the same scale as the
	// analytical eigenvectors. 
	cout << calculated_eigvecs*(analytical_eigvecs(1,1)/calculated_eigvecs(1,1)) << endl; 
	cout << "Analytical eigenvalues of A: " << endl; 
	for (int j=0; j<N; j++) {
		cout << analytical_eigvals(j) << endl; 
		}
	cout << "Matrix containing analytical eigenvectors of A: " << endl; 
	cout << analytical_eigvecs << endl; 
	
	vec difference_eigenvalues = abs(analytical_eigvals-calculated_eigvals); 
	mat difference_eigenvectors = abs(abs(analytical_eigvecs)-abs(calculated_eigvecs*(analytical_eigvecs(1,1)/calculated_eigvecs(1,1)))); 
	double tolerance = pow(10.0,-8.0); 
	cout << "The difference between the eigenvectors are shown in the matrix below: " << endl;
	cout << difference_eigenvectors << endl; 
	// Will print out to terminal if differences between the calculated
	// and analytical eigenvalues or eigenvectors are greater than the tolerance
	for (int i=0; i<N; i++) {
		if (difference_eigenvalues(i)>tolerance) {
			cout << "The difference between the calculated and analytical eigenvalues " << endl; 
			cout << "is larger than tolerance: " << to_string(tolerance) << endl; 
		}
		for (int j=0; j<N; j++) {
			if (difference_eigenvectors(i,j) > tolerance) {
				cout << "The difference in between the calculated and analytical " << endl; 
				cout << "eigenvectors is larger than the tolerance: " << to_string(tolerance) << endl; 
			}
		}
	
	}
	
	return 0; 


}
