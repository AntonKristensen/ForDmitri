#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values

#include "src/matrices.cpp"
#include "src/diffeq.cpp"


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
        if (n > 1000){
            std::cout << "Oh no it didn't converge! :|" << std::endl;
            break; // Breaks if it has been going on for too long
        }
        n++;
    }while(true);
    std::cout << "Steps taken: " << n << std::endl;
    return x;
}











/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Various functions


double rosenbrock(std::vector<double> params){
    return pow(1 - params[0], 2) + 100 * pow(params[1] - pow(params[0], 2), 2);
}


double himmelblau(std::vector<double> params){
    return pow(params[0] * params[0] + params[1] - 11, 2) + pow(params[0] + params[1] * params[1] - 7, 2);
}













/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){  

    std::cout << "Testing:" << std::endl;
    std::cout << "Value of Rosenbrock(0,0) = " << std::endl;
    double rtest = rosenbrock({0,0});
    std::cout << rtest << std::endl;

    std::cout << "Gradient of Rosenbrock(0,0) = " << std::endl;
    std::vector<double> rgrad = gradient(rosenbrock, {0,0});
    vectorprint(rgrad);

    std::cout << "Hessian of Rosenbrock(0,0) = " << std::endl;
    matrix rhess = hessian(rosenbrock, {0,0});
    rhess.print();

    std::cout << "Minimum of Rosenbrock(0,0): " << std::endl;
    std::vector<double> rosenmin = newton(rosenbrock, {0,0});
    vectorprint(rosenmin);
    std::cout << "Rosenbrock(minimum) = " << rosenbrock(rosenmin) << std::endl;

    std::cout << "Value of Himmelblau(0,0) = " << std::endl;
    double htest = himmelblau({0,0});
    std::cout << htest << std::endl;

    std::cout << "Gradient of Himmelblau(0,0) = " << std::endl;
    std::vector<double> hgrad = gradient(himmelblau, {0,0});
    vectorprint(hgrad);

    std::cout << "Hessian of Himmelblau(0,0) = " << std::endl;
    matrix hhess = hessian(himmelblau, {0,0});
    hhess.print();

    std::cout << "Minimum of Himmelblau(0,0): " << std::endl;
    std::vector<double> himmelmin = newton(himmelblau, {0,0});
    vectorprint(himmelmin);
    std::cout << "Himmelblau(minimum) = " << himmelblau(himmelmin) << std::endl;

    std::cout << "This isn't a correct minimum of the Himmelblau function." << std::endl;
    std::cout << "My algorithm does, however, find a correct minimum if the starting point is better:" << std::endl;

    std::cout << "Minimum of Himmelblau(2,3): " << std::endl;
    std::vector<double> chimmelmin = newton(himmelblau, {2,3});
    vectorprint(chimmelmin);
    std::cout << "Himmelblau(minimum) = " << himmelblau(chimmelmin) << std::endl;





    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Higgs!

    //Loading data
    std::vector<double> es;
    std::vector<double> sigs;
    std::vector<double> esigs;

    std::ifstream data("src/data.csv");
    std::string line;
    size_t start = 0;
    while (std::getline(data, line)) {
        size_t end;
        es.push_back(std::stod(line.substr(start), &end));
        sigs.push_back(std::stod(line.substr(start + end), &end));
        esigs.push_back(std::stod(line.substr(start + end), &end));
    }
    data.close();


    auto breitwignerfit = [es, sigs, esigs](std::vector<double> params){
        double result = 0;
        for (int i=0; i<es.size(); i++){
            result += pow((params[2] / (pow(es[i] - params[0],2) + params[1] * params[1] / 4) - sigs[i]) / esigs[i], 2);
        }
        //std::cout << result << std::endl;
        return result;
    };

    std::vector<double> higgsmin = newton(breitwignerfit,{125, 5, 10});
    std::cout << "D = " << breitwignerfit(higgsmin) << std::endl;
    std::cout << "With parameters m,Γ,A = ";
    vectorprint(higgsmin);
    std::ofstream fit;
    fit.open("output/fitresult.txt");
    fit << "m = " << higgsmin[0] << "\ngamma = " << higgsmin[1] << "\na = " << higgsmin[2] << std::endl;
    fit.close();


    std::cout << "Guess by inspection: D = " << breitwignerfit({126, 3, 14}) << std::endl;
    std::cout << "With parameters m,Γ,A = ";
    vectorprint({126, 3, 14});
    std::cout << "So maybe my minimization algorithm isn't the best at finding such difficult minima :(" << std::endl;

    return 0;
}