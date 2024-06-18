
#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions


// using namespace std;

#include "src/matrices.cpp"
#include "src/eigen.cpp"

int main(int argc, char *argv[]){
    // Making a file to output to
    std::ofstream outputting;
    outputting.open("output/eigen.txt");
    

    //////////////////////////////
    // This part is for diagonalization of a random square matrix

    matrix M{.x=8, .y=8};
    //M = identitymatrix(4,4);
    //M.elements[0][2] = 4.8;
    M.randfill(100);
    M.forcesymmetry();
    outputting << "Original matrix:" << std::endl;
    M.print(outputting);
    matrix original = M;


   

    outputting << "Diagonalizing" << std::endl;
    matrix eigenvectors = cycle(M);
    outputting << "Diagonalized matrix:" << std::endl;
    M.print(outputting);
    outputting << "Eigenvectors:" << std::endl;
    eigenvectors.print(outputting);
    outputting << "VDV^T:" << std::endl;
    matrix VDV = matmult(eigenvectors, matmult(M,transpose(eigenvectors)));
    VDV.print(outputting);
    outputting << "VDV is = Matrix: " << samematrix(VDV,original) << std::endl;
    outputting << "AVA is = D: " << samematrix(matmult(transpose(eigenvectors), matmult(original,eigenvectors)),M) << std::endl;
    outputting << "V^TV is = I: " << samematrix(matmult(transpose(eigenvectors),eigenvectors),identitymatrix(M.x, M.x)) << std::endl;
    outputting << "VV^T is = I: " << samematrix(matmult(eigenvectors,transpose(eigenvectors)),identitymatrix(M.x, M.x)) << std::endl;




    outputting.close(); // Closing the output file

    return 0;
}


