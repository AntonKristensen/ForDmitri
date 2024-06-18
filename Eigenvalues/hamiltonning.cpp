#include <iostream>     // For printing to terminal for debugging etc
#include <fstream>      // For writing to files
#include <ctime>        // For timing
#include <vector>       // I want the std::vector<double> since it knows its own length
#include <cmath>        // fancy math functions

#include "src/matrices.cpp"
#include "src/eigen.cpp"
  

int main(int argc, char *argv[]){  
    // Making a file to output to
    std::ofstream outputting;
    outputting.open("output/hamil.txt");

    // Now for the Schrödinger equation

    // The variables come from calling the program with them, so like ./program.exe 10 0.1
    double rmax = std::stod(argv[1]); // The 0th index is always just the program's own name
    double dr = std::stod(argv[2]); // From the command line you input strings, so the stod() makes them into doubles

    int len = (int) (rmax/dr); // Making an integer that is the size of the physical space in r
    std::vector<double> space; // A vector that has the space coordinates, for convenient plotting later on
    space.resize(len);
    for (int i=0; i<len; i++){
        space[i] = dr * (i+1);
    }

    // Making the Hamiltonian. It might be more performant to make a seperate submatrices by some fancy routine and then add them.
    matrix H{.x=len, .y=len};
    H.fill(0);
    for (int i=0; i<len; i++){
        for (int j=0; j<len; j++){
            if (i==j){ // The diagonal elements
                H.elements[i][j] = +1 / dr  / dr  -1/space[i] ; // The i+1 is to make sure there is no 1/0, so the space physically starts at r=dr, not r=0.
            }
            if (j==i+1 || i==j+1){ // The off-diagonal elements
                H.elements[i][j] = -0.5 / dr  / dr ; 
            }
        }
    }

    // Making a copy of the matrix for later use
    matrix original = H;

    // Diagonalizing the Hamiltonian
    matrix wavefunctions = cycle(H);
    std::vector<double> energies; // A vector that will contain the energies. I can use it to find the lowest eigenvalues and later on plot the corresponding wavefunctions
    energies.resize(len);
    for (int i=0; i<len; i++){
        energies[i] = H.elements[i][i]; // Filling in the eigenvalues
    }
    //std::vector<double> groundstate = wavefunctions.col(0);




    // Making an output file stream where I can put values etc and stuff
    std::ofstream energyoutput;
    energyoutput.open("output/energies.txt", std::ofstream::app); // Appending mode so it doesn't overwrite the values from last run
    energyoutput << energies[0] << std::endl;
    std::cout << "Energies: " << energies[0] << " and " << energies[1] << " and " << energies[2] <<   std::endl;
    //std::cout << "Groundstate: " << wavefunctions.col(0)[0] << " and " << wavefunctions.col(0)[1] << " and " << wavefunctions.col(0)[2] << " and " << wavefunctions.col(0)[3] << " and " << wavefunctions.col(0)[4] <<   std::endl;
    //std::cout << "1st excited: " << wavefunctions.col(1)[0] << " and " << wavefunctions.col(1)[1] << " and " << wavefunctions.col(1)[2] << " and " << wavefunctions.col(1)[3] << " and " << wavefunctions.col(1)[4] <<   std::endl;
    energyoutput.close();

    std::ofstream eigenfunctions;
    eigenfunctions.open("output/eigenfunctions.txt");
    for (int i=0; i<len; i++){
        eigenfunctions << space[i] << ", " << wavefunctions.col(0)[i] << ", " << wavefunctions.col(1)[i]  << ", " << wavefunctions.col(2)[i]  << ", " << wavefunctions.col(3)[i] << std::endl;
    }
    eigenfunctions.close();

    std::ofstream states;
    states.open("output/states.txt");
    for (int i=0; i<len; i++){
        states << space[i] << ", " << wavefunctions.col(0)[i] *space[i] << ", " << wavefunctions.col(1)[i] *space[i] << ", " << wavefunctions.col(2)[i] *space[i] << ", " << wavefunctions.col(3)[i] *space[i] << std::endl;
    }
    states.close();

    outputting.close(); // Closing the file

    return 0;
}