#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions
#include <functional>   // convenient way to use std::function as datatype so I don't have to actually think about pointers
#include <limits>       // To allow for infinite values

#include "src/matrices.cpp"
#include "src/diffeq.cpp"
#include "src/mini.cpp"



class nn{
    public:
        int n; // Number of neutrons. There's just the single layer
        std::function<double(double)> f = [](double x){return x * exp(-x*x);}; // Using Gaussian wavelets as activation function
        std::vector<double> p; 

        // constructor
        nn(int neurons){
            n = neurons;
            p.resize(n*3); // just a list of (a_i's, b_i's, w_i's)
            randvec(p,1); // Initializes all parameters to be random between -1 and 1
        }

        double response(double point){
            double result = 0;

            for (int i=0; i<n; i++){
                result += f((point - p[i])/p[n+i]) * p[2*n + i]; // The one-dimensional output is just the sum of the responses*weights
            }
            return result;
        }

        
     

        // A function that trains the network to some data x,y
        void train(std::vector<double> x, std::vector<double> y){

            // Define a cost function that can then be minimized by my minimization routine. If I had lots of time I might have made it so that the cost function is either a parameter of the nn class, or an input to the train() function    
            std::function<double(std::vector<double>)> cost = [xs=x, ys=y, N=n, func=f](std::vector<double> params){
                double result = 0;

                for (int j=0; j<xs.size(); j++){
                    double respons = 0;
                    for (int i=0; i<N; i++){
                        respons += func((xs[j] - params[i])/params[N+i]) * params[2*N + i]; // The one-dimensional output is just the sum of the responses*weights
                        }
                    result += pow(respons - ys[j] ,2); 
                }

                return result;
            };
                
                std::cout << cost(p) << std::endl;
            p = newton(cost, p);            
        }

};






int main(){

    std::vector<double> testx = linspace(-1, 1, 0.01);
    std::vector<double> testy;
    testy.resize(testx.size());

    std::function<double(double)> testfunc = [](double x){return  sin(x);}; // Testing on a gaussian

    for (int i=0; i<testx.size(); i++){
        testy[i] = testfunc(testx[i]); // Filling in the tabulated values of the testfunction
    }


    nn network{10};

    network.train(testx, testy);

    double point = 0.221;
    vectorprint(network.p);
    std::cout << "The functions value is: " << testfunc(point) << ", and the NN interpolation is: " << network.response(point) << std::endl;
    


    return 0;
}


