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


double erf(double z){
    if(z<0){
        return - erf(-z);
    }
    else if(z >= 0 && z <= 1){
        return 2/sqrt(3.1415926535897932384626433832795) * integrate([](double x){return exp(-x * x);}, 0, z);
    }
    else if(z>1){
        return 1 - 2/sqrt(3.1415926535897932384626433832795) * integrate([z](double t){return exp(-(z+(1-t)/t) * (z+(1-t)/t)) /t /t;}, 0, 1);
    }
    else {
        return NAN;
    }
}


int main(int argc, char *argv[]){  

    
    // The small functions that need to be integrated
    std::cout << "Testing my integration on the first four functions:" << std::endl;
    int called = calls;
    double number1 = integrate([](double x){return sqrt(x);}, 0, 1);
    std::cout << "Should be 2/3: " << number1 << " #calls: " << calls - called << std::endl;
    
    called = calls;
    double number2 = integrate([](double x){return 1/sqrt(x);}, 0, 1);
    std::cout << "Should be 2: " << number2 << " #calls: " << calls - called << std::endl;
    int invsqrtcalls = calls - called;
    
    called = calls;
    double number3 = integrate([](double x){return 4*sqrt(1-x*x);}, 0, 1);
    std::cout << "Should be pi: " << number3 << " #calls: " << calls - called << std::endl;
    
    called = calls;
    double number4 = integrate([](double x){return log(x)/sqrt(x);}, 0, 1);
    std::cout << "Should be -4: " << number4 << " #calls: " << calls - called << std::endl;
    int loginvsqrtcalls = calls - called;

    // Making a linspace with points where the errorfunction is calculated on
    std::vector<double> xs = linspace(-5, 5, 0.1);
    std::vector<double> errs;
    errs.resize(xs.size());
    for (int i=0; i< xs.size(); i++){
        errs[i] = erf(xs[i]);
    }
    twovectorwrite(xs, errs, "output/erf.txt");

    // Try the Clenshaw-Curtis variable transformation:
    std::cout << "Testing the Clenshaw-Curtis transformation:" << std::endl;
    called = calls;
    double number5 = cctransform([](double x){return 1/sqrt(x);}, 0, 1);
    std::cout << "Should be 2: " << number5 << " #calls: " << calls - called << " compared to the non transformed's: " << invsqrtcalls << std::endl;

    called = calls;
    double number6 = cctransform([](double x){return log(x)/sqrt(x);}, 0, 1);
    std::cout << "Should be -4: " << number6 << " #calls: " << calls - called << " compared to the non transformed's: " << loginvsqrtcalls << std::endl;

    // Now with error estimate
    std::cout << "Now with error estimate!" << std::endl;
    called = calls;
    std::vector<double> errortest = integratewitherror([](double x){return 1/sqrt(x);}, 0, 1);
    std::cout << "Should be 2: " << errortest[0] << ", error: " << errortest[1] << ", #calls: " << calls - called << std::endl;

    // Indefinite integrals!
    std::cout << "Testing infinite limits:" << std::endl;
    called = calls;
    std::vector<double> topinf = infiniteintegral([](double x){return 1/(x*x);}, 1, std::numeric_limits<double>::infinity()); // f(x) = 1/x^2
    std::cout << "Should be 1: " << topinf[0] << ", error: " << topinf[1] << ", #calls: " << calls - called << std::endl;

    called = calls;
    std::vector<double> botinf = infiniteintegral([](double x){return exp(x);}, -std::numeric_limits<double>::infinity(), 1); // f(x) = exp(x)
    std::cout << "Should be e: " << botinf[0] << ", error: " << botinf[1] << ", #calls: " << calls - called << std::endl;

    called = calls;
    std::vector<double> bothinf = infiniteintegral([](double x){return exp(-x*x);}, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()); // f(x) = exp(-x^2)
    std::cout << "Should be sqrt(pi) ≈ 1.77245: " << bothinf[0] << ", error: " << bothinf[1] << ", #calls: " << calls - called << std::endl;

    return 0;
}





