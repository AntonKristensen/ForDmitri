#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values

#include "src/matrices.cpp"
#include "src/diffeq.cpp"
#include "src/root.cpp"





// I use bubblesort to sort my d_i in ascending order. I already had it implemented in my "interpolations.cpp".
void bubblesort(std::vector<double> &x){ // Yeahyeah not the fastest but it's not like my vectors are big
    bool changed = false;
    double oldx;
    do{
        changed = false;
        for (int i=0; i<x.size()-1; i++){
            if (x[i]>=x[i+1]){
                oldx = x[i];
                x[i] = x[i+1];
                x[i+1] = oldx;
                changed = true;
            }
        }
    }while(changed);
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){  
    srand(time(0)); // Seeding the RNG

    // The dimension of the matrix problem
    int N = std::stoi(argv[1]); // The dimension is given as argument in the makefile

    
    // The diagonal matrix doesn't need elements on the off-diagonal, so I can just represent it by a vector
    std::vector<double> d;
    d.resize(N);
    randvec(d,1); // Filling the vector with random values between -1000 and 1000
    bubblesort(d);  // Sorting in ascending order
    // If you really need to start from some non-ascending order diagonal matrix, then you can do QR decomposition and work with the matrix A' = Q^T D Q + Q^T u u^T Q, since QR decomposition will put the diagonals in ascending order. This also works for any random (symmetric) matrix.

    // The vector which updates the diagonal matrix as D -> D + σ u u^T
    std::vector<double> u;
    u.resize(N);
    randvec(u,1); // Filling the vector with random values between -1000 and 1000

    // The value of σ. It is randomly generated, uniformly distributed from -1 to 1
    double sigma = -1;//((double)rand()/ (double)RAND_MAX - 0.5) * 2;

    // Making the secular equation as a lambda expression. It captures the d and u vectors and σ from above.
    auto secular = [d, u, sigma](std::vector<double> params){ // For beauty's sake it should take a double as input, but my root finding algorithm needs it to take a vector, so it's just a 1 element list. I could make an overloaded version of the root finding algorithm if I had infinite free time, but I'm going to CERN tomorrow so I don't feel like it right now.
        std::vector<double> result = {1};
        for (int i=0; i<d.size(); i++){
            result[0] += sigma * u[i] * u[i] / (d[i] - params[0]); // Fill in the equation
        }
        return result;
    };

    std::vector<double> roots;
    roots.resize(N);
    double guess;

    // Timing tests
    clock_t start, stop;
    start = clock();
    // Running the actual root finding
    for (int i=0; i<N; i++){
        // These are the two cases, that split up the ranges that the roots can be in, depending on the sign of sigma.
        // I'm guessing for the root to be in the middle of the range
        if (sigma >= 0){
            if (i == N-1){
                //guess = (2 * d[i] + sigma * dot(u, u))/2;
                guess = d[i] + sigma * dot(u, u); // Might as well guess at the top
                roots[i] = newton(secular, {guess})[0];
            }else {
                guess = (d[i] + d[i+1])/2;
                roots[i] = newton(secular, {guess})[0];
                if (roots[i] < d[i] || roots[i] > d[i+1]){ 
                    // It turns out that if the root is not found by guessing in the middle, then the root is very close to d[i+1], so I guess near that

                    //std::cout << "Fake root at i= " << i  <<  ". root - guess = " << roots[i] - guess << ". Guessed: " << guess << std::endl;
                    //std::cout << secular({d[i] + (d[i+1] - d[i])/100000})[0] << ", " << secular({d[i] + 999*(d[i+1] - d[i])/1000})[0] << std::endl;
                    roots[i] = newton(secular, {d[i] + 999*(d[i+1] - d[i])/1000}, 0.001)[0];
                    //std::cout << "Found actual root! at: " << roots[i] << ", with value: " << secular({roots[i]})[0] << std::endl;
                }
            }
        }else if (sigma <= 0){
            if (i == 0){
                guess = (2 * d[i] + sigma * dot(u, u))/2;
            }else {
                guess = (d[i] + d[i-1])/2;
            }
            roots[i] = newton(secular, {guess})[0];
            if (roots[i] < d[i-1] || roots[i] > d[i]){ // Same here. If the root is hard to find then it's because it's close to d[i]
                roots[i] = newton(secular, {d[i-1] + (d[i] - d[i-1])/1000}, 0.001)[0];
            }
        }
        
    }
    stop = clock();
    double time = (double) (stop-start) / double(CLOCKS_PER_SEC);

    // Writes out the time it took
    std::ofstream out;
    out.open("output/timing.txt", std::ofstream::app);
    out << N << ", " << time << std::endl;
    out.close();

    return 0;
}