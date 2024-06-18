
#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions


// using namespace std;

#include "src/matrices.cpp"
#include "src/eigen.cpp"

int main(){

    matrix M{.x=4, .y=4};
    M = identitymatrix(4,4);
    M.elements[0][2] = 4.8;
    M.randfill(100);
    M.forcesymmetry();
    std::cout << "Original matrix:" << std::endl;
    M.print();
    matrix original = M;


   

    std::cout << "Diagonalizing" << std::endl;
    matrix eigenvectors = cycle(M);
    std::cout << "Diagonalized matrix:" << std::endl;
    M.print();
    std::cout << "Eigenvectors:" << std::endl;
    eigenvectors.print();
    std::cout << "VDV^T:" << std::endl;
    matrix VDV = matmult(eigenvectors, matmult(M,transpose(eigenvectors)));
    VDV.print();
    std::cout << "VDV is = Matrix: " << samematrix(VDV,original) << std::endl;
    std::cout << "AVA is = D: " << samematrix(matmult(transpose(eigenvectors), matmult(original,eigenvectors)),M) << std::endl;
    std::cout << "V^TV is = I: " << samematrix(matmult(transpose(eigenvectors),eigenvectors),identitymatrix(M.x, M.x)) << std::endl;
    std::cout << "VV^T is = I: " << samematrix(matmult(eigenvectors,transpose(eigenvectors)),identitymatrix(M.x, M.x)) << std::endl;





    return 0;
}


