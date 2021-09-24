#ifndef TRIDIAGMAT_HPP
#define TRIDIAGMAT_HPP

#include <armadillo> 
#include <string>
#include <fstream>
#include <cmath>

using namespace arma; 
using namespace std; 


class TriDiagMat {
protected: 
	
public: 
	mat initialize_from_scalars(int N, double a, double d, double e); 
	mat initialize_symmetric_from_scalars(int N, double a, double d); 
	mat initialize_from_vectors(int N, const vec& a, const vec& d, const vec& e); 
	mat initialize_dense_from_scalars(int N, double a, double d); 
}; 







#endif
