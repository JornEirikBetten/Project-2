#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <armadillo>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>

using namespace arma; 
using namespace std; 

class Jacobi {
protected: 

public: 
	double max_offdiag_symmetric(const mat& A, int& k, int& l); 
	void jacobi_rotate(mat& A, mat& R, int k, int l); 
	void jacobi_eigensolver(mat&, double, vec&, mat&, const int, int&, bool&); 
	void sort_cols(mat&, vec&); 
	void swap_columns(mat&, int, int); 
	void write_to_csv(mat&, vec&, string filename); 
};

#endif
