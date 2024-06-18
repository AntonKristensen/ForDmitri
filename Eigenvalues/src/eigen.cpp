#include <iostream> // For printing
#include <ctime> // For generating a starting seed for RNG
#include <vector> // For basic vectors
#include <cmath>


//#include "matrices.cpp" // My linear algebra classes etc

// Right multiplication by J
void timesJ(matrix &M, int p, int q, double theta){
    double cosine = cos(theta);
    double sine = sin(theta);
    for (int i=0; i<M.x; i++){
        double aip = M.elements[i][p];
        double aiq = M.elements[i][q];

        M.elements[i][p] = cosine * aip - sine * aiq;
        M.elements[i][q] = sine * aip + cosine * aiq;
    }
}

// Left multiplication by J
void Jtimes(matrix &M, int p, int q, double theta){
    double cosine = cos(theta);
    double sine = sin(theta);
    for (int j=0; j<M.x; j++){
        double apj = M.elements[p][j];
        double aqj = M.elements[q][j];

        M.elements[p][j] = cosine * apj - sine * aqj;
        M.elements[q][j] = sine * apj + cosine * aqj;
    }
}

// Function that cycles through a matrix until it has converged below some threshold
matrix cycle(matrix &M){
    matrix V = identitymatrix(M.x, M.y);
    if (M.x != M.y){
        std::cout << "Mega Oof, the dimensions don't match! >:[" << std::endl;
    }
    bool changed = false; // Boolean that controls the loop
    int loopcounter = 0;

    do{ // Loop to do the EVD
        changed = false;

        //std::cout << "Loop number " << loopcounter << std::endl; 
        loopcounter++;
        for (int p=0; p<M.x-1; p++){
            for (int q=p+1; q<M.x; q++){
                // Making variables
                double apq = M.elements[p][q];
                double app = M.elements[p][p];
                double aqq = M.elements[q][q];
                double theta = 0.5 * atan2(2 * apq, aqq-app);
                double cosine = cos(theta);
                double sine = sin(theta);

                // New values of the diagonal, to be able to see how much the diagonal changes, so the program can check convergence
                double napp = cosine * cosine * app - 2 * sine * cosine * apq + sine * sine * aqq;
                double naqq = sine * sine * app + 2 * sine * cosine * apq + cosine * cosine * aqq;

                //std::cout << p << "," << q << " with " << abs(napp - app) << " and " << abs(naqq - aqq) << "; theta= " << theta << std::endl;
                if (napp =! app || naqq != aqq){
                    changed = true;
                    timesJ(M, p, q, theta);
                    Jtimes(M, p, q, theta);
                    timesJ(V,p,q, theta);
                    //std::cout << "The Matrix is now:" << std::endl; M.print();
                }
            }
        }
    }while(changed && loopcounter < 100);

    return V;
}