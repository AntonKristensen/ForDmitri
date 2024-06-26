#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values

#include "src/matrices.cpp"
#include "src/interpolations.cpp"
#include "src/integration.cpp"


// Plain Monte Carlo returns a 2-element vector, where the first is the estimate of the integral, and the second is the estimate of the error
// The input is a function which returns a double. Its input is a single vector that holds all the double-valued parameters. This is so that I can make N-dimensional integrals.
std::vector<double> plainmc(std::function<double(std::vector<double>)> f, std::vector<double> start, std::vector<double> stop, int N){
    std::vector<double> result;
    result.resize(2);

    int dim = start.size();
    double V = 1;
    for (int i=0; i<dim; i++){
        V *= stop[i] - start[i]; // This makes a rectangular N-dimensional hypervolume
    }

    double sum = 0;
    double sum2 = 0;

    std::vector<double> x;
    x.resize(dim);

    // This little algorithm does not save the points, which is more memory efficient, but makes it harder to plot the points' distribution :(
    for (int i=0; i<N; i++){
        for (int j=0; j<dim; j++){
            x[j] = start[j] + (stop[j] - start[j]) * (double) rand() / (double)RAND_MAX; // Make the point in N-space, uniformly distributed in the hypervolume 
        }
        double fx = f(x);
        sum += fx; // Just the sum of the evaluated functions
        sum2 += fx*fx; // Sum of squares for error estimation
    }

    result[0] = V * sum / N;
    result[1] = V * sqrt(sum2 / N - sum / N * sum / N) / sqrt(N);

    return result;
}


// Functions for making Coput and Halton distributed numbers
double corput(int n, int b){
    
    double q = 0;
    double bk = (double) 1 / b;
    while(n>0){
        q += (double) (n % b) * bk;
        n /= b;
        bk /= b;
    }
    return q;
}


std::vector<double> halton(int n, int d){
    std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    std::vector<double> result;
    result.resize(d);
    for (int i=0; i<d; i++){
        result[i] = corput(n, primes[i]);
    }
    return result;
}


std::vector<double> lattice(int n, int d, std::vector<double> alphas){
    std::vector<double> result;
    result.resize(d);

    double irrelevantpointer[1];
    for (int i=0; i<d; i++){
        result[i] = modf((double) n * alphas[i], irrelevantpointer);
    }
    
    return result;
}


// Quasi-random Monte Carlo returns a 2-element vector, where the first is the estimate of the integral, and the second is the estimate of the error
// The input is a function which returns a double. Its input is a single vector that holds all the double-valued parameters. This is so that I can make N-dimensional integrals.
std::vector<double> quasimc(std::function<double(std::vector<double>)> f, std::vector<double> start, std::vector<double> stop, int N){
    std::vector<double> result;
    result.resize(2);

    int dim = start.size();
    double V = 1;
    for (int i=0; i<dim; i++){
        V *= stop[i] - start[i]; // This makes a rectangular N-dimensional hypervolume
    }

    double sum1 = 0;
    double sum2 = 0;

    std::vector<double> x1;
    x1.resize(dim);
    std::vector<double> x2;
    x2.resize(dim);

    std::vector<double> alphas; // Vector with fractional numbers for the lattice quasi-random generator
    alphas.resize(dim);
    for (int i=0; i<dim; i++){
        alphas[i] = (double) rand() / (double)RAND_MAX;
    }

    // This little algorithm does not save the points, which is more memory efficient, but makes it harder to plot the points' distribution :(
    for (int i=0; i<N; i++){
        x1 = halton(i, dim); // Value from the Halton
        double fx1 = f(x1);

        x2 = lattice(i, dim, alphas);
        double fx2 = f(x2);

        sum1 += fx1;
        sum2 += fx2; 
    }

    result[0] = V * (sum1 + sum2)/2 / N; // Returns the average of the two methods
    result[1] = V * fabs(sum2-sum1)/2 /N; // Estimates the error to be the difference between the two methods

    return result;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions

// Unit N-Sphere
double S(std::vector<double> point){
    if (dot(point,point) <= 1){ // Checks if the point is within the N-sphere
        return 1;
    }
    else {
        return 0;
    }
}



// Quick uniform gaussian function, works for any dimension
double gauss(std::vector<double> point){
    double value = 1;
    for (int i=0; i<point.size(); i++){
        value *= exp(- point[i] * point[i]); // It's just a sum of exponentials
    }
    return value;
}


// That fancy very singular function
double singular(std::vector<double> point){
    return 0.032251534433 / (1 - cos(point[0]) * cos(point[1]) * cos(point[2])); // 0.032251534433 is 1/pi^3
}











//////////////////////////////////////////////////////////


int main(int argc, char *argv[]){  
    srand(time(0)); // Seeding the RNG


    // Printing the results of MC integration on a circle as a function of number of points
    std::ofstream circleout;
    circleout.open("output/circle.txt");
    std::vector<double> circle;
    std::vector<double> quasi;
    for (int i=100; i<=10000; i += 1 + (int) sqrt(i)){ // This iteration is just so the points are spaced farther and farther apart the bigger the N is
        circle = plainmc(S, {-1, -1}, {1, 1}, i);
        circleout << i << ", " << circle[0] << ", " << circle[1];
        // Also doing it with the quasi random method for comparison
        quasi = quasimc(S, {-1, -1}, {1, 1}, i);
        circleout << ", " << quasi[0] << ", " << quasi[1]<< std::endl;
    }
    circleout.close();


    std::vector<double> gamma = plainmc(singular, {0, 0, 0}, {3.1415965, 3.1415965, 3.1415965}, 100000);
    std::cout << "Should be 1.3932039: " << gamma[0] << ", error: " << gamma[1] << std::endl;

    // Trying the 
    std::vector<double> quasigamma = quasimc(S, {-1, -1, -1}, {1, 1, 1}, 100000);
    std::cout << "Should be 4/3 pi (4.18879): " << quasigamma[0] << ", error: " << quasigamma[1] << std::endl;

    
}