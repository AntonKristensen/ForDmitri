
///////


#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers

#include "src/matrices.cpp"
#include "src/interpolations.cpp"
#include "src/diffeq.cpp"


// Convenient function for writing two vectors into a file so Gnuplot can plot them
void twovectorwrite(std::vector<double> xs, std::vector<double> ys, std::string filename){
    for (int i=0; i<xs.size(); i++){
        std::ofstream out;
        out.open(filename);
        for (int i=0; i<xs.size(); i++){
            out << xs[i] << ", " << ys[i] << std::endl;
        }
        out.close();
    }
}



//////////////////////////////////////
// Various differential equations:



std::vector<double> harmonic(double x, std::vector<double> y){
    // Simple harmonic motion
    // d^2/dt^2(r) = -r

    std::vector<double> result;
    result.resize(y.size());

    result[0] = y[1]; // d/dt(y1) = y2
    result[1] = -y[0]; // d/dt(y2) = -y1

    return result;
}



std::vector<double> damped(double x, std::vector<double> y){
    // Damped harmonic motion
    // d^2/dt^2(r) = -r - 0.2 * d/dt(r)

    std::vector<double> result;
    result.resize(y.size());

    result[0] = y[1]; // d/dt(y1) = y2
    result[1] = -y[0] - 0.2 * y[1]; // d/dt(y2) = -y1 - y2

    return result;
}


std::vector<double> driven(double x, std::vector<double> y){
    // Damped harmonic motion
    // d^2/dt^2(r) = -r - 0.2 * d/dt(r) + sin(t)

    std::vector<double> result;
    result.resize(y.size());

    result[0] = y[1]; // d/dt(y1) = y2
    result[1] = -y[0] - 0.2 * y[1] + sin(x); // d/dt(y2) = -y1 - y2

    return result;
}



std::vector<double> newton(double x, std::vector<double> y){
    // Newtonian orbital motion
    // d^2/dt^2(u) + u = 1

    std::vector<double> result;
    result.resize(y.size());

    result[0] = y[1]; // d/dt(y1) = y2
    result[1] = 1 - y[0]; // d/dt(y2) = 1 - y1

    return result;
}

std::vector<double> albert(double x, std::vector<double> y){
    // Relativistc orbital motion
    // d^2/dt^2(u) + u = 1 + e * u^2

    std::vector<double> result;
    result.resize(y.size());

    result[0] = y[1]; // d/dt(y1) = y2
    result[1] = 1 - y[0] + 0.01 * y[0] * y[0]; // d/dt(y2) = 1 - y1

    return result;
}


///////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){  

    // Showing that my ODE solver works by solving some harmonic differential equations
    oderesult simple = driver(harmonic, 2*3.14159 * 5, {1, 0});
    simple.write("output/harmonic.txt");

    oderesult damp = driver(damped, 2*3.14159 * 5, {1, 0});
    damp.write("output/damped.txt");

    oderesult driv = driver(driven, 2*3.14159 * 5, {1, 0});
    driv.write("output/driven.txt");

    double range = 2*3.14159 * 5;
    linterp intie = interpode(harmonic, range, {1, 1}, 1, 1, 1); // Making h, acc and eps big to maybe show the interpolation a bit
    std::vector<double> interpolationpoints = linspace(0,range,0.1);
    twovectorwrite(interpolationpoints, evalinterpode(intie,interpolationpoints), "output/interpolated.txt");

    // Now for orbital motion!
    double circlerange = 2*3.14159 * 0.9;
    linterp circle = interpode(newton, circlerange, {1, 0}); // This should give a circle
    std::vector<double> circleinterpolationpoints = linspace(0,circlerange,0.1);
    twovectorwrite(circleinterpolationpoints, evalinterpode(circle,circleinterpolationpoints), "output/circular.txt");
    
    double ellipserange = 2*3.14159 * 0.9;
    linterp ellipse = interpode(newton, ellipserange, {1, -0.5}); // This should give an ellipse
    std::vector<double> ellipseinterpolationpoints = linspace(0,ellipserange,0.1);
    twovectorwrite(ellipseinterpolationpoints, evalinterpode(ellipse,ellipseinterpolationpoints), "output/ellipse.txt");

    double relrange = 2*3.14159 * 10;
    linterp rel = interpode(albert, relrange, {1, -0.5}); // This should give an ellipse
    std::vector<double> relinterpolationpoints = linspace(0,relrange,0.1);
    twovectorwrite(relinterpolationpoints, evalinterpode(rel,relinterpolationpoints), "output/relativistic.txt");

    return 0;
}
