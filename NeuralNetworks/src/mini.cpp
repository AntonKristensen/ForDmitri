#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values


/*
// Super special function for finding the gradient, only works for cost function of neural networks where the activation functions are gaussian wavelets
std::vector<double> gradient(std::function<double(std::vector<double>)> f, std::vector<double> x){
    std::vector<double> result;
    result.resize(x.size());

    // Notice that x (which is the vector with all the parameters of the neural network) comes sizes of multiples of 3
    // The first third is the a_i's, the middle third is the b_i's and the last third is the w_i's
    for (int i=0; i<x.size(); i++){
        // I do the differentiation myself, but the derivative depends on whether the parameter is an a, b or w:
        // The derivative of t*exp(-t²) is = exp(-t²) - 2t² exp(-t²)
        // The variable t is of course t = (x-a)/b
        // The value of x is captured from when newton is called inside the nn class where the cost function captures stuff
        if (i < x.size()/3){
            result[i] = 3.14;
        }else if(i < 2 * x.size() / 3){
            result[i] = 3.14;
        }else {
            result[i] = 3.14;
        }
    }

    return result;
}
*/


// Function for finding the gradient
std::vector<double> gradient(std::function<double(std::vector<double>)> f, std::vector<double> x){
    std::vector<double> result;
    result.resize(x.size());

    double fx = f(x);

    for (int i=0; i<x.size(); i++){
        double dx = std::max(fabs(x[i]), 1.0) * pow(2, -26);

        x[i] += dx;
        result[i] = (f(x) - fx) / dx;
        x[i] -= dx;

    }

    return result;
}





// Function for finding the Hessian
matrix hessian(std::function<double(std::vector<double>)> f, std::vector<double> x){
    matrix result{x.size(), x.size()};

    std::vector<double> gradx = gradient(f, x); 

    for (int j=0; j<x.size(); j++){

        double dx = std::max(fabs(x[j]), 1.0) * pow(2, -13);

        x[j] += dx;
        std::vector<double> dgradx = add(gradient(f, x), scale(gradx, -1));
        for (int i=0; i< x.size(); i++){
            result.elements[i][j] = dgradx[i] / dx;
        }
        x[j] -= dx;
    }

    return result;
}



// Newton's method for finding minimum of a function
std::vector<double> newton(std::function<double(std::vector<double>)> f, std::vector<double> x, double acc = 0.001){ // A bit bad naming it the same as the method from the "roots" homework, so need to beware if I want to use them both in the same program some time!

    int n=0; // Keeps track of how many iterations has been done
    do{
        std::vector<double> gradx = gradient(f, x);
        if (norm(gradx) < acc){
            break;
        }

        matrix H = hessian(f, x);

        QR decomped = decomp(H);
        std::vector<double> dx = decomped.solve(scale(gradx, -1));

        double lambda = 1;
        double fx = f(x);

        do{
            if (f(add(x, scale(dx, lambda))) < fx){
                break;
            }
            if (lambda < pow(2, -26)){
                break; // If lambda gets mega small then do the step anyways
            }
            lambda /= 2;

        }while(true);
        x =add(x, scale(dx, lambda));
        if (n > 50){
            std::cout << "Oh no it didn't converge! :|" << std::endl;
            break; // Breaks if it has been going on for too long
        }
        n++;
    }while(true);
    std::cout << "Steps taken: " << n << std::endl;
    return x;
}


