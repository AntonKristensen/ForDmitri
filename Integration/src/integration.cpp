#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values

int calls = 0; // Super hack way to keep track of how many times integrate() has been called.
double integrate(std::function<double(double)> f, double start, double stop, double delta=0.001, double epsilon=0.001, double f2=NAN, double f3 = NAN){
    calls++;

    double h = stop-start;

    if (std::isnan(f2)){
        f2 = f(start + 2*h/6);
        f3 = f(start + 4*h/6);
    }

    double f1 = f(start + h/6);
    double f4 = f(start + 5*h/6);

    double Q = (2*f1 + f2 + f3 + 2*f4)/6 * (stop - start);
    double q = (f1 + f2 + f3 + f4)/4 * (stop - start);
    double error = fabs(Q - q);

    if (error <= delta + epsilon *fabs(Q)){
        return Q;
    }
    else {
        return integrate(f, start, (start + stop)/2, delta / sqrt(2), epsilon, f1, f2) + integrate(f, (start + stop)/2, stop, delta / sqrt(2), epsilon, f3, f4);
    }
}



double cctransform(std::function<double(double)> f, double start, double stop, double delta=0.001, double epsilon=0.001){
    // Lambda expressions are really convenient for stuff like this!
    return integrate([start, stop, f](double t){return f((start+stop)/2 + (stop-start) / 2*cos(t)) * sin(t) * (stop-start)/2;}, 0, 3.1415926535897932384626433832795, delta, epsilon);
}


std::vector<double> integratewitherror(std::function<double(double)> f, double start, double stop, double delta=0.001, double epsilon=0.001, double f2=NAN, double f3 = NAN){
    calls++;

    std::vector<double> result;
    result.resize(2);

    double h = stop-start;

    if (std::isnan(f2)){
        f2 = f(start + 2*h/6);
        f3 = f(start + 4*h/6);
    }

    double f1 = f(start + h/6);
    double f4 = f(start + 5*h/6);

    double Q = (2*f1 + f2 + f3 + 2*f4)/6 * (stop - start);
    double q = (f1 + f2 + f3 + f4)/4 * (stop - start);
    double error = fabs(Q - q);

    if (error <= delta + epsilon *fabs(Q)){
        result[0] = Q;
        result[1] = error;
        return result;
    }
    else {
        return add(integratewitherror(f, start, (start + stop)/2, delta / sqrt(2), epsilon, f1, f2), integratewitherror(f, (start + stop)/2, stop, delta / sqrt(2), epsilon, f3, f4));
    }
}

std::vector<double> infiniteintegral(std::function<double(double)> f, double start, double stop, double delta=0.001, double epsilon=0.001, double f2=NAN, double f3 = NAN){

    // The three cases for integrals with infinite limits, which need to be transformed away
    if (start == -std::numeric_limits<double>::infinity() && stop == std::numeric_limits<double>::infinity()){
        return integratewitherror([f](double t){return f(t / (1-t*t)) * (1 + t*t) / ((1-t*t)*(1-t*t));}, -1, 1, delta, epsilon);
    }
    else if(start == -std::numeric_limits<double>::infinity()){
        return integratewitherror([f, stop](double t){return f(stop + t/(1+t)) * 1/((1+t)*(1+t));}, -1, 0, delta, epsilon);
    }else if(stop == std::numeric_limits<double>::infinity()){
        return integratewitherror([f, start](double t){return f(start + t/(1-t)) * 1/((1-t)*(1-t));}, 0, 1, delta, epsilon);
    }else {
        return integratewitherror(f, start, stop, delta, epsilon);
    }


}

