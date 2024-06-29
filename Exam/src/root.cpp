#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values


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
    int count = 0;
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
            if (lambda < pow(2,-10)){; //std::cout << "small lambda!" << std::endl;
                break; // Stop if lambda gets below certain precision
            }
            lambda /= 2;       
        }while(true);
        x = z;
        fx = fz;
    if (count > 50){std::cout << "Oof the rootfinder algorithm didn't find a root :(" << std::endl; break;} // 50 is a decent ish value for this specific problem here. Normally it takes less than 10 to find a root
    count++;
    }while(true);
    return x;
}


