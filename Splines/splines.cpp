#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions

#include "src/matrices.cpp"
#include "src/eigen.cpp"
#include "src/interpolations.cpp"



int main(int argc, char *argv[]){  
    srand(time(0)); // Seeding the RNG

    // Making a vector of random length
    double xrange = 10;
    int minsize = 10;
    int maxsize = 200;
    int len;
    do{
        len = ((double)rand()/ (double)RAND_MAX) * maxsize ;
    }while(len <= minsize); // Makes sure the length isn't TOO small
    
    std::vector<double> x, y; x.resize(len); y.resize(len);
    randvec(x,xrange);
    randvec(y,0.1); // The 0.1 adds a bit of uniformly distributed 10% noise
    for (int i=0; i<len; i++){
        y[i] = sin(x[i]) + y[i];
    }

    std::ofstream data;
    data.open("output/data.txt");
    for (int i=0; i<len; i++){
        data << x[i] << ", " << y[i] << std::endl;
    }
    data.close();
    

    //////////////////////////////////////////////////////////////////////////////
    // Let's do some linear interpolation!


    // Creating an instance of my linear interpolation class on the given data vectors x and y
    linterp linnie{x,y};
    

    // Making a linspace vector of xpoints that I want to evaluate the interpolation in
    std::vector<double> inx;
    inx.resize(1000-2);
    vectorfill(inx, linnie.x[0]); // Starts out at the minimum value, but will get added to
    double range = linnie.x[linnie.x.size()-1] - linnie.x[0];
    for (int i=0; i<inx.size(); i++){
        inx[i] += ((double)i+1)/(inx.size()+1) * range; // Makes sure I dont try to interpolate farther than the data allows
    }
    
    // Putting evaluations into file for plotting
    std::ofstream lin;
    lin.open("output/lin.txt");
    for (int i=0; i<inx.size(); i++){
        lin << inx[i] << ", " << linnie.eval(inx[i]) << std::endl;
    }
    lin.close();

    // Putting integrals into file for plotting
    std::ofstream linint;
    linint.open("output/linint.txt");
    for (int i=0; i<inx.size(); i++){
        linint << inx[i] << ", " << linnie.inteval(inx[i]) << std::endl;
    }
    linint.close();


    ////////////////////////////////
    // Quadratic spline now

    qinterp qinnie{x,y};

    // Putting evaluations into file for plotting
    std::ofstream qin;
    qin.open("output/qin.txt");
    for (int i=0; i<inx.size(); i++){
        qin << inx[i] << ", " << qinnie.eval(inx[i]) << std::endl;
    }
    qin.close();

    // Putting integrals into file for plotting
    std::ofstream qinint;
    qinint.open("output/qinint.txt");
    for (int i=0; i<inx.size(); i++){
        qinint << inx[i] << ", " << qinnie.inteval(inx[i]) << std::endl;
    }
    qinint.close();


    // Putting differentials into file for plotting
    std::ofstream dinint;
    dinint.open("output/dinint.txt");
    for (int i=0; i<inx.size(); i++){
        dinint << inx[i] << ", " << qinnie.difeval(inx[i]) << std::endl;
    }
    dinint.close();



    return 0;
}