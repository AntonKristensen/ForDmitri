#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions

#include "src/matrices.cpp"


int main(int argc, char *argv[]){  


    // Loading data from experiment
    std::vector<double> times = {1,  2,  3, 4, 6, 9,   10,  13,  15};
    std::vector<double> activity = {117,100,88,72,53,29.5,25.2,15.2,11.1};
    std::vector<double> sigma = {6,5,4,4,4,3,3,2,2};

    std::vector<double> lnactivity; // Making a vector that has the logarithm of the activities, because that's how we get the fit model to be a linear one
    lnactivity.resize(times.size());
    for (int i=0; i<times.size(); i++){
        lnactivity[i] = log(activity[i]); 
    }


    // Setting up the least squares matrix
    matrix M{.x=times.size(), .y=2}; // As many rows as measurements, as many columns as parameters in the fitting model.
    std::vector<double> ones;
    ones.resize(times.size());
    vectorfill(ones, 1.0); 
    M.setcol(0,ones); // The first parameter is constant in x, so just the coefficient times 1.0 
    M.setcol(1,times); // The second parameter is linear in x
    
    //M.print();
    QR decomped = decomp(M);
    std::vector<double> solution = decomped.solve(lnactivity); // Using the logarithm of the activities

    std::ofstream datas;
    datas.open("output/datas.txt");
    for (int i=0; i<times.size(); i++){
        datas << times[i] << ", " << activity[i] << ", " << sigma[i] << std::endl;
    }
    datas.close();

    std::ofstream fitresult;
    fitresult.open("output/fitresult.txt");
    // This is a super hack way of storing the result from the fit in a way that Gnuplot can use easily
    fitresult << "a = " << exp(solution[0]) <<  std::endl; // a is the exponential of ln(a), and ln(a) is the linear parameter that we fitted
    fitresult << "lambda = " << - solution[1] <<  std::endl; // lambda is = -c2
    fitresult.close();



    return 0;
}