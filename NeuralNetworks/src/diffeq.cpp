///////


#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers




// This only takes a function which takes a double as the first input and a vector (which has y[0] and y[1]) as the second input. This can only solve very few diffeqs, so if I wanted it to do more,
// I would rewrite it so that any function just has a list or a struct as input, which then has variable size depending on the specific function, which the stepper and driver function would take into account
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
oderesult driver(std::function<std::vector<double>(double, std::vector<double>)> f, double a, double b, std::vector<double> y0, double h = 0.01, double acc = 0.001, double eps = 0.001){
    double x = a; //
    std::vector<double> y = y0;


    // Making output
    std::vector<double> xs;
    xs.push_back(x);

    std::vector<std::vector<double>> ys;
    ys.push_back(y);

    do{
        // To make sure the last step ends exactly at b
        if (x + h > b){
            h = b - x;
        }

        std::vector<std::vector<double>> step = stepper(f, x, y, h);

        double tolerance = (acc + eps * norm(step[0])) * sqrt(h / (b-a));
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

/*
// A function the returns an interpolated list of the ode solution result
// I have made it such that you pass the interpode() function a list of x points where you want the solution evaluated in. It then solves the ODE once, and then returns a list of interpolated y_i's.
linterp interpode(std::function<std::vector<double>(double, std::vector<double>)> f, double b, std::vector<double> y0, double h = 0.1, double acc = 0.01, double eps = 0.01){
    oderesult ode = driver(f, b, y0, h, acc, eps);

    // Filling a vector with the y values
    std::vector<double> ys;
    ys.resize(ode.x.size());
    for (int i=0; i<ode.x.size(); i++){
        ys[i] = ode.y[i][0];
    }

    linterp interpolation{ode.x, ys};
    return interpolation;
}

// Convenient function for evaluating an interpolation from the function above at some points xs
std::vector<double> evalinterpode(linterp interpolation, std::vector<double> xs){
    std::vector<double> ys;
    ys.resize(xs.size());

    for (int i=0; i<xs.size(); i++){
        ys[i] = interpolation.eval(xs[i]);
    }
    
    return ys;
}
*/