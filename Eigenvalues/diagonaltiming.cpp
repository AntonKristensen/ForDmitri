

#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions

#include "src/matrices.cpp"
#include "src/eigen.cpp"
  

int main(int argc, char *argv[]){  // Pass the size of the matrix as a parameter

    int n = std::stoi(argv[1]);

    matrix M{.x=n, .y=n};
    M.randfill(100);
    M.forcesymmetry();
    clock_t start, stop;
    start = clock();
    cycle(M);
    stop = clock();

    double time = (double) (stop-start) / double(CLOCKS_PER_SEC);
    
    std::cout << n << ", " << time << std::endl;




    return 0;
}