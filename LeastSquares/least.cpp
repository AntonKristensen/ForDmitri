#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions

#include "src/matrices.cpp"
#include "src/eigen.cpp"

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
    vectorprint(lnactivity);

    std::vector<double> lnsigma; // Similarly making a vector that has the uncertainties of the logarithm
    lnsigma.resize(times.size());
    for (int i=0; i<times.size(); i++){
        lnsigma[i] = sigma[i] / activity[i]; 
    }

    // Least squares that takes uncertainty in the data into account has Aik = fk(xi)/sigmayi and bi = yi / sigmayi

    // Setting up the b vector

    std::vector<double> b = lnactivity;
    for (int i=0; i<times.size(); i++){
        b[i] /= lnsigma[i];
    }


    // Setting up the least squares matrix
    matrix M{.x=times.size(), .y=2}; // As many rows as measurements, as many columns as parameters in the fitting model.
    std::vector<double> ones;
    ones.resize(times.size());
    vectorfill(ones, 1.0); 
    M.setcol(0,ones); // The first parameter is constant in x, so just the coefficient times 1.0 
    M.setcol(1,times); // The second parameter is linear in x
    // To take the uncertainties into account every i'th row in A must be scaled with sigmai
    for (int i=0; i<times.size(); i++){
        M.elements[i][0] /= lnsigma[i];
        M.elements[i][1] /= lnsigma[i];
    }
    
    
    //M.print();
    QR decomped = decomp(M);
    std::vector<double> solution = decomped.solve(b); // Using the logarithm of the activities

    // Outputting the data into a file
    std::ofstream datas;
    datas.open("output/datas.txt");
    for (int i=0; i<times.size(); i++){
        datas << times[i] << ", " << activity[i] << ", " << sigma[i] << ", " << lnactivity[i] << ", " << lnsigma[i] << std::endl;
    }
    datas.close();

    // Outputting the result of the fit so it can be plotted (is done in the file "gnuplotter" when called by "make plot")
    std::ofstream fitresult;
    fitresult.open("output/fitresult.txt");
    // This is a super hack way of storing the result from the fit in a way that Gnuplot can use easily
    fitresult << "a = " << exp(solution[0]) <<  std::endl; // a is the exponential of ln(a), and ln(a) is the linear parameter that we fitted
    fitresult << "lambda = " << - solution[1] <<  std::endl; // lambda is = -c2
    fitresult.close();


    // The covariance matrix of A is (A^T A)^-1, so I use the QR decomposition on A^t A to find the inverse of my matrix
    matrix ATA = matmult(transpose(M), M);
    QR ATAdecomped = decomp(ATA);
    matrix covariance = ATAdecomped.inverse();
    //covariance.print();
    fitresult.open("output/fitresult.txt", std::ofstream::app);
    fitresult << "sigma_ln(a) = " << sqrt(covariance.elements[0][0]) <<  std::endl; // This is the uncertainty in the fitting coefficient
    fitresult << "sigma_a = " << sqrt(covariance.elements[0][0]) * exp(solution[0]) <<  std::endl; // By using error propagation, f(a) = e^a: sigma_a = sqrt([e^a * sigma_f(a)]^2)
    fitresult << "sigma_(minuslambda) = " << sqrt(covariance.elements[1][1]) <<  std::endl; // 
    fitresult << "sigma_lambda = " << sqrt(covariance.elements[1][1]) <<  std::endl; // This one is just the f(lambda) = -lambda, so the minus sign goes away when squaring in the error propagation formula
    fitresult.close();
    
    // Calculating the half-life if it is of interest
    double halflife = -log(2) / solution[1];
    double sigmahalflife = log(2) / (solution[1] * solution[1]) * sqrt(covariance.elements[1][1]); // From propagation
    std::cout << "The half-life is " << halflife << " ± " << sigmahalflife << " days" << std::endl;
    std::cout << "Compared to 3.631(2) days from NIST this is " << fabs(halflife - 3.631) / sigmahalflife << "σ away, so not really indicative of the same value, but honestly not too bad for the time" << std::endl;
    return 0;
}