#include <iostream> // For printing
#include <fstream>      // For writing to files
#include <ctime> // For generating a starting seed for RNG
#include <vector> // For basic vectors
#include <cmath>



class matrix {
    public:
        int x;
        int y;
        std::vector<std::vector<double>> elements;

        matrix(int n, int m) : x(n), y(m) {
            elements.resize(x);
            for (int i=0; i<x;i++){
                elements[i].resize(y);
            }
        }

        // Fill with some value
        void fill(double value){
            for (int i=0; i<x;i++){
                for (int j=0; j<y;j++){
                    elements[i][j] = value;
                }
            }
        }
        // Fill with random value between -value and value (default = 1)
        void randfill(double value=1){
            for (int i=0; i<x;i++){
                for (int j=0; j<y;j++){
                    elements[i][j] = ((double)rand()/ (double)RAND_MAX - 0.5) * value * 2;
                }
            }
        }
        // Print out the matrix
        void print(std::ostream& output = std::cout){
            for (int i=0; i<x;i++){
                output << "|";
                for (int j=0; j<y;j++){
                    output << elements[i][j] << ", ";
                }
                output<<"\b \b"; // This removes the extra space
                output<<"\b \b"; // This removes the extra comma
                output << "|" << std::endl;
            }
        }

        // C is row major order, but I want to also be able to easily work on the columns
        std::vector<double> col(int j){
            std::vector<double> jcol;
            jcol.resize(x);
            for (int i=0; i<x; i++){
                jcol[i] = elements[i][j];
            }
            return jcol;
        }
        void setcol(int j, std::vector<double> colvector){
            for (int i=0; i<x; i++){
                elements[i][j] = colvector[i];
            }
        }

        // Function that overrides the bottom left triangle with the values from upper right, so to force the matrix to be symmetric
        void forcesymmetry(){
            for (int i=1; i<x;i++){
                for (int j=0; j<i;j++){
                    elements[i][j] = elements[j][i];
                }
            }
        }
};






////////////////////////////////////////////////////////////////////
// Vector functions

void vectorfill(std::vector<double> &vec, double value=1){
    for (int i=0; i<vec.size();i++){
        vec[i] = value;
    }
}


void randvec(std::vector<double> &vec, double value=1){
    for (int i=0; i<vec.size();i++){
        vec[i] = ((double)rand()/ (double)RAND_MAX - 0.5) * value * 2;
    }
}

std::vector<double> scale(std::vector<double> vec, double scalar){
    for (int i=0; i<vec.size(); i++){
       vec[i] = vec[i] * scalar;
    }
    return vec;
}

std::vector<double> add(std::vector<double> vec1, std::vector<double> vec2, std::ostream& output = std::cout){
    std::vector<double> result;
    if (vec1.size() == vec2.size()) {
        result.resize(vec1.size()); 
        for (int i=0; i<vec1.size(); i++){
            result[i] = vec1[i] + vec2[i];
        }
        return result;
    }
    else{
        output << "Mega Oof, the dimensions don't match! >:[" << std::endl;
        return result;
    }
}

// Norm of vector
double norm(std::vector<double> vec){
    double sum = 0;
    for (int i=0; i<vec.size(); i++){
        sum += pow(vec[i],2);
    }
    return sqrt(sum);
}

// Dot product between two vectors
double dot(std::vector<double> vec1, std::vector<double> vec2, std::ostream& output = std::cout){
    double sum = 0;
    if (vec1.size() == vec2.size()) {
        for (int i=0; i<vec1.size(); i++){
            sum += vec1[i] * vec2[i];
        }
        return sum;
    }
    else{
        output << "Mega Oof, the dimensions don't match! >:[" << std::endl;
        return -1;
    }
}

// Convenient function for printing vectors
void vectorprint(std::vector<double> vec, std::ostream& output = std::cout){
    output << "|";
    for (int i=0; i<vec.size(); i++){
        output << vec[i] << ", ";
    }
    output<<"\b \b"; // This removes the extra space
    output<<"\b \b"; // This removes the extra comma
    output << "|" << std::endl;
}

// Convenient function for writing two vectors into a file so Gnuplot can plot them
void twovectorwrite(std::vector<double> xs, std::vector<double> ys, std::string filename){
    for (int i=0; i<xs.size(); i++){
        std::ofstream out;
        out.open(filename);
        for (int i=0; i<xs.size(); i++){
            out << xs[i] << ", " << ys[i] << std::endl;
        }
        out.close();
    }
}


// Check if two vectors are roughly equal
bool samevector(std::vector<double> v1, std::vector<double> v2, double threshold = 1/10000000000){
    bool sameness = true;
    if (v1.size() == v2.size()){
        for (int i=0; i<v1.size(); i++){
            if (abs(v1[i] - v2[i]) > threshold){
                sameness = false;
            }
        }
    }
    else{
        sameness = false;
    }
    return sameness;
}


// Neat little function for making a linspace
std::vector<double> linspace(double start, double stop, double interval){
    int len = (int) ((stop - start) / interval);
    std::vector<double> xs;
    xs.resize(len+1);
    for (int i=0; i<len+1; i++){
        xs[i] = start + (i) * interval;
    }
    return xs;
}



// A function for transposing matrices
matrix transpose(matrix M){
    matrix transposed(M.y, M.x);
    for (int i=0; i<M.y; i++){
        for (int j=0; j<M.x; j++){
            transposed.elements[i][j] = M.elements[j][i];
        }
    }
    return transposed;
}

// Multiplying a matrix with a vector
std::vector<double> matvecmult(matrix M, std::vector<double> vec, std::ostream& output = std::cout){
    std::vector<double> result;
    result.resize(M.x);
    if (M.y == vec.size()){
        for (int i=0; i<M.x; i++){
            for (int j=0; j<M.y; j++){
                result[i] += M.elements[i][j] * vec[j];
            }
        }
    }
    else{
        output << "Mega Oof, the dimensions for matrix-vector don't match! >:[" << std::endl;
    }
    return result; // Will be empty if dimensions are bad
}

// Multiplying matrices
matrix matmult(matrix M1, matrix M2, std::ostream& output = std::cout){
    matrix result(M1.x, M2.y);
    if (M1.y == M2.x){ // This checks dimensions real quick
        result.fill(0.0);

        for (int i=0; i<M1.x; i++){
            for (int j=0; j<M2.y; j++){
                for (int k=0; k<M1.y; k++){
                    result.elements[i][j] = result.elements[i][j] + M1.elements[i][k] * M2.elements[k][j];
                }
            }
        }
    }
    else {
        output << "Mega Oof, the dimensions for the matrix multiplication don't match! >:[" << std::endl;
    }
    
    return result;  
}


// This one checks if matrices are equal
bool samematrix(matrix M1, matrix M2, double threshold=1 / 10000000000 /* The default threshold is 10^-10 */){
    bool sameness = false;
    if (M1.x == M2.x && M1.y == M2.y){
        sameness = true;
        for (int i=0; i<M1.x; i++){
            for (int j=0; j<M1.y; j++){
                if (abs(M1.elements[i][j] - M2.elements[i][j]) > threshold ){
                    sameness = false;
                }
            }
        }
    }
    return sameness;
}


// This one checks if a matrix is upper triangular
bool istriangular(matrix M, double threshold=1 / 10000000000 /* The default threshold is 10^-10 */){
    bool triangular = true;
    for (int i=1; i<M.x; i++){
        for (int j=0; j<i; j++){
            if (abs(M.elements[i][j]) > threshold ){
                triangular = false;
            }
        }
    }
    return triangular;
} // The method by only checking within a threshold is to handle floating point errors etc. It's a dumb method, so if the matrices only contain VERY small numbers then it's kind of not the best way to go about it

// This function makes an identity matrix of size n,m
matrix identitymatrix(int n,int m){
    matrix ident{.x = n, .y=m};
    ident.fill(0.0);
    for (int i=0; i<n; i++){
        ident.elements[i][i] = 1;
    }
    return ident;
}


// Make a class that holds both Q and R matrices
class QR{
    public:
        matrix Q;
        matrix R;


        // Making a function to solve linear systems of equations
        std::vector<double> solve(std::vector<double> vec, std::ostream& output = std::cout){
            if (Q.x != Q.y){
            //    output << "Careful, the matrix is not square!" << std::endl;
            }
            std::vector<double> c = matvecmult(transpose(Q), vec);
            std::vector<double> y;
            y.resize(R.y);
            y[R.y-1] = vec[R.y-1] / R.elements[R.y-1][R.y-1];
            for (int i=c.size()-1; i>=0; i--){
                double sum=0;
                for (int k=i+1; k<R.x; k++){
                    sum += R.elements[i][k] * y[k];
                }
                y[i] = (c[i] - sum) / R.elements[i][i];
            }
            return y;
        }

        // Function for finding the determinant of the original matrix, which is just = det(R)
        double det(){
            double result = 0;
            for (int i=0; i<R.x; i++){
                result += R.elements[i][i];
            }
            return result;
        }

        // Function for finding the inverse of the original matrix
        matrix inverse(std::ostream& output = std::cout){
            matrix result{.x=Q.y, .y=Q.x};

            if (Q.x != Q.y){
                //output << "Careful, the matrix is not square!" << std::endl;
            }
            for (int i=0; i<Q.x; i++){
                // Making the ei vectors which are zero except at the i'th index
                std::vector<double> ei;
                ei.resize(Q.y);
                vectorfill(ei,0.0);
                ei[i] = 1;

                // Solving the equation to get a column vector
                std::vector<double> xi = solve(ei);
                
                // Setting that column vector into the result matrix
                result.setcol(i,xi);
            }
            return result;
        }


};

// This will do the decompositions
QR decomp(matrix M){
    QR decomped {
        .Q = M,
        .R = matrix(M.y, M.y)
    }; //Initialize an instance of the QR class

    decomped.R.fill(0.0); //Make the R matrix have values of 0 everywhere

    for (int i=0; i<M.y; i++){
        decomped.R.elements[i][i] = norm(decomped.Q.col(i));
        decomped.Q.setcol(i, scale(decomped.Q.col(i) , 1/decomped.R.elements[i][i]));

        for (int j=i+1; j<M.y; j++){
            decomped.R.elements[i][j] = dot(decomped.Q.col(i), decomped.Q.col(j));

            decomped.Q.setcol(j, add(decomped.Q.col(j) , scale(decomped.Q.col(i), - decomped.R.elements[i][j])) );
        }
    }
    return decomped;


}

