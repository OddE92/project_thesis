#include "calculate_eigenvalues.h"

#include <vector>

#include <armadillo>

using std::vector;


int calculate_eigenvalues_3x3_sym(vector<double> &inVector, int startMatrix, vector<double> &outVector){

    //Assumes a 3x3 symmetric matrix, taking the top half of the matrix, including the diagonal
    //It also assumes it gets the full vector with all values, not only the 3x3 matrix you want
    //Thus startMatrix is also passed to find the start of the matrix.
    //This function is as dumb as they come, and assumes the theory is correct.

    arma::mat A(3, 3, arma::fill::zeros);

    A(0,0) = inVector[startMatrix + 0];   A(1,0) = inVector[startMatrix + 1];   A(2,0) = inVector[startMatrix + 2];
    A(0,1) = inVector[startMatrix + 1];   A(1,1) = inVector[startMatrix + 3];   A(2,1) = inVector[startMatrix + 4];
    A(0,2) = inVector[startMatrix + 2];   A(1,2) = inVector[startMatrix + 4];   A(2,2) = inVector[startMatrix + 5];

    arma::vec b(3);

    arma::eig_sym(b, A);

    for(int i = 0; i < 3; i++){
        outVector[i] = b(i);
    }
    
    return 0;
}

