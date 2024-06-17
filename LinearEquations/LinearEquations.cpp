
#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions


// using namespace std;

#include "src/matrices.cpp"



int main() {
    srand(time(0)); // Seeding the RNG
    std::cout << "Hello Dmitri!" << std::endl;


    //////////////////////////////
    // I make a tall matrix and put random stuff in it
    matrix tall(7, 5);
    tall.randfill(100);
    //std::cout << "Matrix is:" << std::endl;
    //tall.print();

    // I decompose
    std::cout << "Decomposing tall matrix:" << std::endl;
    QR talldecomp = decomp(tall);
    std::cout << "Decomposed!" << std::endl;

    // Matrices to check
    matrix QQ = matmult(transpose(talldecomp.Q), talldecomp.Q); // Q^T Q
    matrix QnR = matmult(talldecomp.Q, talldecomp.R); // Q R

    // Printing matrices out explicitly if I wanna look at them
    /* do
    std::cout << "R is:" << std::endl;
    talldecomp.R.print();
    std::cout << "Q is:" << std::endl;
    talldecomp.Q.print();

    std::cout << "Q^T is:" << std::endl;
    matrix qposed = transpose(talldecomp.Q);
    qposed.print();
    std::cout << "Q^T Q is:" << std::endl;
    
    QQ.print();
    std::cout << "Q R is:" << std::endl;
    matrix QnR = matmult(talldecomp.Q, talldecomp.R);
    QnR.print();
    std::cout << "Matrix is:" << std::endl;
    tall.print();
    */



    // Here I check that the matrices are as they should be
    std::cout << "R is upper triangular: " << istriangular(talldecomp.R) << std::endl;
    std::cout << "Q^T Q is I: " << samematrix(QQ,identitymatrix(QQ.x, QQ.y)) << std::endl;
    std::cout << "QR is = Matrix: " << samematrix(tall,QnR) << std::endl;

    


    /////////////////////////////
    // Now for making a square matrix and solving it

    // Making a square matrix and putting random stuff in it
    int squaredimension = 5;
    matrix square(squaredimension,squaredimension);
    square.randfill(100);
    QR squaredecomp = decomp(square); // QR decomposition of the matrix



    // Making a vector to set up equation system
    std::vector<double> b;
    b.resize(square.x);
    randvec(b,100);
    std::cout << "Solving square matrix equation:" << std::endl;
    std::vector<double> solution = squaredecomp.solve(b); // Solving the equation
    std::cout << "Solved!" << std::endl;
    
    /*
    std::cout << "b is: " << std::endl;
    vectorprint(b);
    std::cout << "testing: \nx=" << std::endl;
    vectorprint(solution);

    std::cout << "b " << std::endl;
    vectorprint(b);

    std::cout << "Q^Tb = c " << std::endl;
    vectorprint(matvecmult(transpose(squaredecomp.Q), b));

    std::cout << "Rx " << std::endl;
    vectorprint(matvecmult(squaredecomp.R, solution));

    std::cout << "Ax " << std::endl;
    vectorprint(matvecmult(m, solution));
    
    std::cout << "QRx " << std::endl;
    vectorprint(matvecmult(squaredecomp.Q, matvecmult(squaredecomp.R,solution)));
    */
 
    std::cout << "Rx is Q^Tb: " << samevector(matvecmult(squaredecomp.R, solution) ,matvecmult(transpose(squaredecomp.Q),b)) << std::endl;
    std::cout << "Ax is b: " << samevector(matvecmult(square, solution), b) << std::endl;






    ///////////////////////////////////
    // Matrix inversion

    // Making a matrix that is the inverse of the square matrix from above
    std::cout << "Inverting:" << std::endl;
    matrix inverted = squaredecomp.inverse();
    
    // Testing if A^-1 A = I
    std::cout << "A^-1 A = I: " << samematrix(matmult(inverted, square),identitymatrix(square.x, square.y)) << std::endl;



    //////////////////////////////////
    // Timing tests
    clock_t start, stop;
    int dimensionality;
    int end = 100;

    // Vector to hold data, both axes for plotting
    std::vector<double> dimensions;
    dimensions.resize(end);

    std::vector<double> timeresults;
    timeresults.resize(end);
    



    for (dimensionality=1; dimensionality<=end; dimensionality++){
        // Make matrix
        matrix timetester(dimensionality,dimensionality);
        timetester.randfill(100);

        // Then time the decomposition
        start = clock(); 
        decomp(timetester); // QR decomposition of the matrix
        stop = clock();

        // Write it into a vector that holds the results
        dimensions[dimensionality] = dimensionality;
        timeresults[dimensionality] = (double) (stop-start) / double(CLOCKS_PER_SEC);
    }
    

    std::ofstream myfile;
    myfile.open ("output/timingoutput.txt");
    for (int i=0; i<timeresults.size(); i++){
        myfile << dimensions[i] << ", " <<timeresults[i] << std::endl;
    }
    myfile.close();

    return 0;
}
