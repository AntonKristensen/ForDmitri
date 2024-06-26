// Making a class for linear interpolation
class linterp{
    public:
        std::vector<double> x;
        std::vector<double> y;

        // Constructor
        linterp(std::vector<double> a, std::vector<double>  b) : x(a), y(b) {
            bubblesort(); // Sorts the vectors
        }



    // Binary search only works on sorted lists, so I make a function that sorts the x in ascending order, and also sorts y in that same order
    void bubblesort(){ // Yeahyeah not the fastest but it's not like my vectors are big
        bool changed = false;
        double oldx;
        double oldy;
        do{
            changed = false;
            for (int i=0; i<x.size()-1; i++){
                if (x[i]>x[i+1]){
                    oldx = x[i]; oldy = y[i];
                    x[i] = x[i+1]; y[i] = y[i+1];
                    x[i+1] = oldx; y[i+1] = oldy;
                    changed = true;
                }
            }
        }while(changed);
    }
 
    // Method that finds the two indices of x that the interpolation point is between
    int binarysearch(double point){
        int index = 0;
        if (point < x[0] || point > x[x.size()-1]){
            std::cout << "Oopsie, the point is outside the binsearch data range! :[" << std::endl;
        }
        else{ 
            int j=x.size()-1;
            while (j-index > 1){
                int mid=(index+j)/2;
		        if(point>x[mid]){
                    index=mid;
                }
                else{ 
                    j=mid;
                }
            }
        }
        return index;
    }


    // Method for actually evaluating a point
    double eval(double point){
        double value;
        int index = binarysearch(point);
        double dx=x[index+1]-x[index];
        if (dx == 0){ // If it's actually just evaluated at a known point then just return the known point
            value = y[index];
        }
        else{ // If evaluated at other point, do the algorithm etc
            double dy=y[index+1]-y[index];
            value = y[index]+dy/dx*(point-x[index]); // Yeah it's one extra allocation, but if the compiler is good the eval() will get optimized into one calculation with no extraflorous allocations
        }
        return value; 
    }
    double inteval(double point){
        double value = 0;
        int index = binarysearch(point);
        for (int i=0; i<=index-1; i++){ // This loop runs all index up to the last interval. It calculates piecewise integrals from x[i] to x[i+1]
            double dx=x[i+1]-x[i];
            double dy=y[i+1]-y[i];
            value += y[i]*(x[i+1] - x[i])+dy/dx * (x[i+1]-x[i]) * (x[i+1]-x[i]) / 2; // The integral is the sum of all piecewise integrals, up to "point"
        }
        // This runs the last interval, which is from x[index] to the point where the integral needs to be evaluated in
        double dx=x[index+1]-x[index];
        if (dx == 0){ // If it's actually to a point then there is no integration distance
            // Do nothing
        }
        else{
            double dy=y[index+1]-y[index];
        value += y[index]*(point - x[index])+dy/dx * (point-x[index]) * (point-x[index]) / 2;
        }
        return value; 
    }
};


// Making a class for quadratic interpolation
class qinterp{
    public:
        std::vector<double> x;
        std::vector<double> y;

        // Variables for holding coefficients
        std::vector<double> p;
        std::vector<double> b;
        std::vector<double> dx;
        std::vector<double> cup;
        std::vector<double> cdown;
        std::vector<double> c;
        // Yeah it holds a decent amount of memory, but I don't really think it is that worrysome (scales O(N)). In worst case I can always just make a deallocator method to free up the memory stored

        // Constructor
        qinterp(std::vector<double> a, std::vector<double>  b) : x(a), y(b) {
            bubblesort(); // Sorts the vectors so x is in ascending order, necessary for binary search


            // Variables for holding coefficients 
            p.resize(x.size()-1);
            this->b.resize(x.size()-1); // Need "this->" since I want to store the values. These variables are all just copies otherwise, so "this->" makes it point to the actual
            dx.resize(x.size()-1);
            cup.resize(x.size()-1);
            cdown.resize(x.size()-1);
            this->c.resize(x.size()-1);


            //std::cout << "length" << c.size() <<  std::endl;

            // Finding dx's and p's
            for (int i=0; i<dx.size(); i++){
                dx[i] = x[i+1] - x[i];
                p[i] = dx[i] / (y[i+1] - y[i]);
            }
            // Calculating c upwards
            cup[0] = 0; cup[cup.size()-1] = 0; // Setting the first and last elements to zero () effectively making the interpolation there linear.
            for (int i=0; i<cup.size()-2; i++){ //std::cout << "calculating cup value at index " << i+1 <<  std::endl;
                cup[i+1] = (p[i+1] - p[i] - cup[i] * dx[i]) / dx[i+1];
            }
            // And downwards
            cdown[0] = 0; cdown[cdown.size()-1] = 0; // Setting the first and last elements to zero () effectively making the interpolation there linear.
            for (int i=cdown.size()-2; i>0; i--){ //std::cout << "calculating cdown value at index " << i <<  std::endl;
                cdown[i] = (p[i+1] - p[i] - cup[i+1] * dx[i+1]) / dx[i];
            }
            // Averaging the two c's
            for (int i=0; i<cup.size(); i++){
                this->c[i] = (cup[i] + cdown[i])/2; //std::cout << c[i] << std::endl;
                // Also calculating the b's
                this->b[i] = p[i] - c[i] * dx[i];
            }
        }



        // Binary search only works on sorted lists, so I make a function that sorts the x in ascending order, and also sorts y in that same order
        void bubblesort(){ // Yeahyeah not the fastest but it's not like my vectors are big
            bool changed = false;
            double oldx;
            double oldy;
            do{
                changed = false;
                for (int i=0; i<x.size()-1; i++){
                    if (x[i]>x[i+1]){
                        oldx = x[i]; oldy = y[i];
                        x[i] = x[i+1]; y[i] = y[i+1];
                        x[i+1] = oldx; y[i+1] = oldy;
                        changed = true;
                    }
                }
            }while(changed);
        }
    
        // Method that finds the two indices of x that the interpolation point is between
        int binarysearch(double point){
            int index = 0;
            if (point < x[0] || point > x[x.size()-1]){
                std::cout << "Oopsie, the point is outside the binsearch data range! :[" << std::endl;
            }
            else{ 
                int j=x.size()-1;
                while (j-index > 1){
                    int mid=(index+j)/2;
                    if(point>x[mid]){
                        index=mid;
                    }
                    else{ 
                        j=mid;
                    }
                }
            }
            return index;
        }


        // Method for actually evaluating a point
    double eval(double point){
            double value;
            int index = binarysearch(point);
            if (dx[index] == 0){ // If it's actually just evaluated at a known point then just return the known point
                value = y[index];
            }
            else{ // If evaluated at other point, do the algorithm etc
                value = y[index] + b[index] * (point - x[index]) + c[index] * (point - x[index]) * (point - x[index]); // Yeah it's one extra allocation, but if the compiler is good the eval() will get optimized into one calculation with no extraflorous allocations
            }
            return value; 
    }

        // The integration
        double inteval(double point){ // Not written yet
            double value = 0;
            int index = binarysearch(point);
            for (int i=0; i<=index-1; i++){ // This loop runs all index up to the last interval. It calculates piecewise integrals from x[i] to x[i+1]
                value += y[i] * (x[i+1] - x[i]) + 0.5 * b[i] * (x[i+1] - x[i]) * (x[i+1] - x[i]) + c[i] * (x[i+1] - x[i]) * (x[i+1] - x[i]) * (x[i+1] - x[i]) / 3;
                //value += y[i]*(x[i+1] - x[i])+dy/dx * (x[i+1]-x[i]) * (x[i+1]-x[i]) / 2; // The integral is the sum of all piecewise integrals, up to "point"
            }
            // This runs the last interval, which is from x[index] to the point where the integral needs to be evaluated in
            double dx=x[index+1]-x[index];
            if (dx == 0){ // If it's actually to a point then there is no integration distance
                // Do nothing
            }
            else{
                value += y[index] * (point - x[index]) + 0.5 * b[index] * (point - x[index]) * (point - x[index]) + c[index] * (point - x[index]) * (point - x[index]) * (point - x[index]) / 3;
                //value += y[index]*(point - x[index])+dy/dx * (point-x[index]) * (point-x[index]) / 2;
            }
            return value; 
        }

        // The differentiation
        double difeval(double point){
            double value;
            int index = binarysearch(point);
            /*
            if (dx[index] == 0){ // If it's actually just evaluated at a known point then just return the known point
                value = y[index];
            }
            else{ // If evaluated at other point, do the algorithm etc
                value =  b[index]  + c[index] * (point - x[index]) * 2; // 
            }
            */
            value =  b[index]  + c[index] * (point - x[index]) * 2;
            return value; 
        }
};