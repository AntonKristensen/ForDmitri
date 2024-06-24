
///////


#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers

#include "src/matrices.cpp"


std::vector<std::vector<double>> stepper(std::function<std::vector<double>(double, std::vector<double>)> f, double x, std::vector<double> y, double h){
    std::vector<double> k0 = f(x,y);
    std::vector<double> k1 = f(x + 0.5*h, add(y, scale(k0, 0.5 * h))); // 0.5h is midpoint method
    std::vector<double> yh = add(y, scale(k1, h));
    std::vector<double> ey = scale(add(k1, scale(k0,-1)), h); // Error estimate

    std::vector<std::vector<double>> result;
    result.resize(2);
    result[0] = yh;
    result[1] = ey;
    return result;
}


// Struct to hold the results of an ODE solution
struct oderesult{
    std::vector<double> x = {};
    std::vector<std::vector<double>> y = {};

    // Convenient function for writing the results into a file
    void write(std::string file){
        std::ofstream out;
        out.open(file);
        for (int i=0; i<x.size(); i++){
            out << x[i] << ", ";
            for (int j=0; j<y[i].size(); j++){
                out << y[i][j] << ", ";
            }
            out << std::endl;
        }
        out.close();
    }
};


// The driver function
oderesult driver(std::function<std::vector<double>(double, std::vector<double>)> f, double b, std::vector<double> y0, double h = 0.01, double acc = 0.001, double eps = 0.001){
    double x = 0; // I have defined the interval b - a = b, so a = 0, so that the driver always starts at t=0 because I really like that
    std::vector<double> y = y0;

    //extra = 2;
    //lengthguess = (int) extra * b / h; // To avoid constantly reallocating new memory each step, I allocate too much memory (guessing a factor above), and then at the end I cut down to size.

    // Making output
    std::vector<double> xs;
    //xs.resize(lengthguess);
    //xs[0] = x;
    xs.push_back(x);

    std::vector<std::vector<double>> ys;
    //ys.resize(lengthguess);
    //ys[0] = y;
    ys.push_back(y);

    do{
        // To make sure the last step ends exactly at b
        if (x + h > b){
            h = b - x;
        }

        std::vector<std::vector<double>> step = stepper(f, x, y, h);

        double tolerance = (acc + eps * norm(step[0])) * sqrt(h / b);
        double error = norm(step[1]);

        if (error <= tolerance){
            x += h;
            y = step[0];
            xs.push_back(x);
            ys.push_back(y);
        }

        h *= std::min(pow(tolerance/error, 0.25) * 0.95, 2.0); // Stepsize can max get 2 times bigger than previous step. Hmm I could consider making the power and safety and max multiplication into variables

    }while(x < b);

    oderesult result{.x = xs, .y = ys}; 
    return result;

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




///////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){  

    // Showing that my ODE solver works by solving some harmonic differential equations
    oderesult simple = driver(harmonic, 2*3.14159 * 5, {1, 0});
    simple.write("output/harmonic.txt");

    oderesult damp = driver(damped, 2*3.14159 * 5, {1, 0});
    damp.write("output/damped.txt");

    oderesult driv = driver(driven, 2*3.14159 * 5, {1, 0});
    driv.write("output/driven.txt");


    return 0;
}
