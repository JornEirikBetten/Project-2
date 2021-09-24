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
	int N = atoi(argv[1]); // getting the number of datapoints from terminal
	string filename = "eigensolve_N_" + to_string(N) + ".csv";  // initializing name of output file
	double h = (double) 1./((double) N); // finding stepsize
	double h2 = h*h; // stepsize*stepsize
	double diagonal = 2.*h2; //calculating diagonal
	double subsupdiagonal = -1.*h2; //calculating super and subdiagonal
	TriDiagMat myInit; // using TriDiagMat as matrix-initiator
	mat A = myInit.initialize_symmetric_from_scalars(N, subsupdiagonal, diagonal);  // initiating tridiag matrix
	//mat A = mat(N,N).randn(); //used to check dense matrix
	//A = symmatu(A);  //used to check dense matrix
	vec eigval = vec(N); // vector to store eigenvalues
	mat eigvecs = mat(N,N, fill::eye); // vector to store eigenvectors
	//vec eigval1; //checking eigsym for eigenvalues and eigenvectors
	//mat eigvecs1;   
	//eig_sym(eigval1, eigvecs1, A); 
	//cout << eigval1 << endl; 
	//cout << eigvecs1 << endl; 
	
	Jacobi mySolver;  // Using the Jacobi-class as solver
	double epsilon = pow(10.0, -8.0); // defining tolerance
	int maximum_iterations = pow(10.0, 8.0); // maximum iterations, in case it won't converge
	int iterations = 0; // defining and setting the number of iterations to 0
	bool converged;  // bool set to true when it has converged
	mySolver.jacobi_eigensolver(A, epsilon, eigval, eigvecs, maximum_iterations, iterations, converged); // solving eigenvalue-system, using jacobi_eigensolver (Jacobi.cpp)
	eigvecs = normalise(eigvecs); // Normalizing eigenvectors
	mySolver.write_to_csv(eigvecs, eigval, filename); //writing eigenvectors and eigenvalues to file
	cout << "The number of iterations (similarity transformations) it took to diagonalize the matrix with dimensionality " << to_string(N) << "x" << to_string(N) << " is " << to_string(iterations) << endl; // printing out message with specifics about how code ran. 
	//sort(eigval.begin(), eigval.end()); // to see which eigenvectors were the first three
	//cout << eigval << endl; 
	
	
	return 0; 
}
