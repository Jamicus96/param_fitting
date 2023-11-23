// header guard:
#ifndef fitVar
#define fitVar

// include
#include <iostream>


class FitVar {
    private:
        char* parname;
        double value = 0.0, vprior = 0.0, verr = 0.0, verr_copy = 0.0, vlow = 0.0, vhigh = 0.0;

        bool holdConstant = false;
        unsigned int parIdx = 0;

    public:
        // Constructors
        FitVar() {};
        FitVar(char* Parname, double Value, double Verr, double Vlow, double Vhigh, bool HoldConstant = false);

        // Member function
        char* name() {return parname};
        double& val() {return value};
        double& prior() {return vprior};
        double& err() {return verr};
        double& min() {return vlow};
        double& max() {return vhigh};
        void HoldConstant(bool isTrue);
        bool isConstant() {return holdConstant;};
        unsigned int& ParIdx() {return parIdx;};

        // Destructor
        ~FitVar() {delete parname;};
};

//end header guard
#endif