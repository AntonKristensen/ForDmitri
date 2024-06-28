#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values

#include "src/matrices.cpp"
#include "src/diffeq.cpp"





// Function for finding the jacobian
matrix jacobian(std::function<std::vector<double>(std::vector<double>)> f, std::vector<double> x, std::vector<double> fx = {}, std::vector<double> dx = {}){
    if (dx.size() == 0){
        dx.resize(x.size());
        for (int i=0; i<x.size(); i++){
            dx[i] = std::max(fabs(x[i]), 1.0) * pow(2, -26);
        }
    }

    if (fx.size() == 0){
        std::vector<double> temp = f(x); // A little dimensional hacking :p
        fx.resize(temp.size());
        fx = temp;
    }

    matrix J{fx.size(), x.size()};

    for (int j=0; j<x.size(); j++){
        x[j] += dx[j];
        std::vector<double> df = add(f(x), scale(fx, -1)); // Finite difference, "f(x) - fx"
        for (int i=0; i<fx.size(); i++){
            //std::cout << fx.size() << std::endl;
            J.elements[i][j] = df[i] / dx[j];
        }
        x[j] -= dx[j];
    }
    return J;
}


// Newton method for finding roots
std::vector<double> newton(std::function<std::vector<double>(std::vector<double>)> f, std::vector<double> start, double acc = 0.01, std::vector<double> dx = {}){
    std::vector<double> x = start;
    std::vector<double> fx = f(x);
    std::vector<double> z;
    std::vector<double> fz;

    do{
        if(norm(fx) < acc){
            break;
        }
        matrix J = jacobian(f, x, fx, dx);
        QR decomped = decomp(J);

        std::vector<double> Dx = decomped.solve(scale(fx, -1));
        double lambda = 1;
        do{
            z = add(x, scale(Dx, lambda));
            fz = f(z);


            if (norm(fz) < (1 - lambda/2) * norm(fx)){
                break;
            }
            //if (lambda < pow(2,-10)){
            //    break; // Stop if lambda gets below double precision
            //}
            lambda /= 2;       
        }while(true);
        x = z;
        fx = fz;

    }while(true);
    return x;
}






/////////////////////////////////////////////////////////////////
// Functions

// Making a test function whose i'th output is x_i^i
std::vector<double> testfunc(std::vector<double> params){
    std::vector<double> result;
    result.resize(params.size());
    
    for (int i=0; i<params.size(); i++){
        result[i] = pow(params[i],(i+1));
    }
    return result;
}

// f(x,y,z) = {x, y, z}
std::vector<double> trivialfunc(std::vector<double> params){
    return params;
}


// f_i = x_i^2 - 1
std::vector<double> sombrerofunc(std::vector<double> params){
    std::vector<double> result;
    result.resize(params.size());
    for (int i=0; i<params.size(); i++){
        result[i] = pow(params[i],2) - 1;
    }
    return result;
}

// The gradient of the Rosenbrock valley function
std::vector<double> rosenbrockfunc(std::vector<double> params){
    std::vector<double> result;
    result.resize(params.size());
    result[0] = 2 * (1 - params[0]) * (-1) + 100 * 2 * (params[1] - params[0] * params[0]) * (- 2 * params[0]);
    result[1] = 0 + 100 * 2 * (params[1] - params[0] * params[0]) * 1;
    return result;
}

// The gradient of the Himmelblau function
std::vector<double> himmelblaufunc(std::vector<double> params){
    std::vector<double> result;
    result.resize(params.size());
    result[0] = 2 * (params[0] * params[0] + params[1] - 11) * 2 * params[0] + 2 * (params[0] + params[1] * params[1] -7) * 1;
    result[1] = 2 * (params[0] * params[0] + params[1] - 11) * 1 + 2 * (params[0] + params[1] * params[1] -7) * 2 * params[1];
    return result;
}




////////////////////
// For Schrödinger equation

// Global variables for convenience
double rmin = 0.01;
double rmax = 8;
double odeacc = 0.001;
double odeeps = 0.001;

// Returns the result of an ODE run on the schrodinger equation
oderesult schrodinger(double E = -1){
    auto schrodingerequation = [E](double x, std::vector<double> y){
        std::vector<double> result;
        result.resize(y.size());

        result[0] = y[1]; // 
        result[1] = -2 * (E + 1/x) * y[0]; // The schrodinger equation
        return result;
    };


    double starty = rmin - rmin * rmin;
    double startv = 1 - 2 * rmin;

    oderesult schrode = driver(schrodingerequation, rmin, rmax, {starty, startv}, odeacc, odeeps);
    

    return schrode;
}


    // Returns f(rmax) of an ODE run on the schrodinger equation
    std::vector<double> frmax(std::vector<double> params){
        std::vector<double> result;
        result.resize(1);
        double E = params[0];

        oderesult schrode = schrodinger(E);

        result[0] = schrode.y[schrode.x.size()-1][0];
        return result;
    };






/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){  
   
    /*
    std::cout << "Let's find a root for fi = xi^2 -1:" << std::endl;
    std::vector<double> sombreroroots = newton(sombrerofunc, {0,0});
    std::cout << "Root found at:" << std::endl;
    vectorprint(sombreroroots);
    std::cout << "With fi's being:" << std::endl;
    vectorprint(sombrerofunc(sombreroroots));
    */


    
    std::cout << "Let's find a root for that Rosenbrock function:" << std::endl;
    std::vector<double> rosenroots = newton(rosenbrockfunc, {0,0});
    std::cout << "Root found at:" << std::endl;
    vectorprint(rosenroots);
    std::cout << "With fi's being:" << std::endl;
    vectorprint(rosenbrockfunc(rosenroots));
    
    std::cout << "Let's find a root for that Himmelblau function:" << std::endl;
    std::vector<double> himmelroots = newton(himmelblaufunc, {0,0});
    std::cout << "Root found at:" << std::endl;
    vectorprint(himmelroots);
    std::cout << "With fi's being:" << std::endl;
    vectorprint(himmelblaufunc(himmelroots));


    std::cout << "Now for the Schrödinger equation:" << std::endl;
    std::vector<double> schrodingerroots = newton(frmax, {-0.8});
    std::cout << "Energy of the ground state is: " << schrodingerroots[0] << std::endl;

    oderesult groundstate = schrodinger(schrodingerroots[0]);
    groundstate.write("output/schrodinger.txt");

    // Convergence tests
    rmin = 1;
    std::vector<double> rmins = {};
    std::vector<double> ermins = {};
    do{
        rmins.push_back(rmin); // This appends, so it has to reallocate every iteration, but whatever, I'm pretty sure it's the ODE and root finder that takes the performance
        ermins.push_back(newton(frmax, {-0.8})[0]);
        rmin -= 0.01;
    }while(rmin > 0.01);
    rmin = 0.01;
    twovectorwrite(rmins, ermins, "output/rmin.txt");

    rmax = 0.1;
    std::vector<double> rmaxs = {};
    std::vector<double> ermaxs = {};
    do{
        rmaxs.push_back(rmax); // This appends, so it has to reallocate every iteration, but whatever, I'm pretty sure it's the ODE and root finder that takes the performance
        ermaxs.push_back(newton(frmax, {-0.8})[0]);
        rmax += 0.1;
    }while(rmax < 8);
    rmax = 8;
    twovectorwrite(rmaxs, ermaxs, "output/rmax.txt");

    odeacc = 0.1;
    std::vector<double> accs = {};
    std::vector<double> eaccs = {};
    do{
        accs.push_back(odeacc); // This appends, so it has to reallocate every iteration, but whatever, I'm pretty sure it's the ODE and root finder that takes the performance
        eaccs.push_back(newton(frmax, {-0.8})[0]);
        odeacc -= 0.001;
    }while(odeacc > 0.001);
    odeacc = 0.001;
    twovectorwrite(accs, eaccs, "output/acc.txt");

    odeeps = 0.1;
    std::vector<double> epss = {};
    std::vector<double> eepss = {};
    do{
        epss.push_back(odeeps); // This appends, so it has to reallocate every iteration, but whatever, I'm pretty sure it's the ODE and root finder that takes the performance
        eepss.push_back(newton(frmax, {-0.8})[0]);
        odeeps -= 0.001;
    }while(odeeps > 0.001);
    odeeps = 0.001;
    twovectorwrite(epss, eepss, "output/eps.txt");


    return 0;
}