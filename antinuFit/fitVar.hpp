// header guard:
#ifndef fitVar
#define fitVar

// include
#include <iostream>
#include <iostream>


class FitVar {
    private:
        std::string parname;
        double value = 0.0, vprior = 0.0, verr = 0.0, verr_copy = 0.0, vlow = 0.0, vhigh = 0.0;

        bool holdConstant = false;
        unsigned int parIdx = 0;

    public:
        // Constructors
        FitVar() {};
        FitVar(const FitVar& fit) {
            value = fit.value; vprior = fit.vprior; verr = fit.verr; verr_copy = fit.verr_copy;
            vlow = fit.vlow; vhigh = fit.vhigh; holdConstant = fit.holdConstant; parIdx = fit.parIdx;
        };
        FitVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool HoldConstant = false);

        // Member function
        std::string& name() {return parname;};
        double& val() {return value;};
        double& prior() {return vprior;};
        double& err() {return verr;};
        double& min() {return vlow;};
        double& max() {return vhigh;};
        void HoldConstant(bool isTrue);
        bool isConstant() {return holdConstant;};
        unsigned int& ParIdx() {return parIdx;};

        // Destructor
        ~FitVar() {};
};


FitVar::FitVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool HoldConstant) {
    parname = Parname;
    value = Value;
    verr_copy = verr;  // Stored error info, in case is verr is set to 0 temporarily, to hold it constant
    vlow = Vlow;
    vhigh = Vhigh;

    holdConstant = HoldConstant;
    if (!holdConstant) verr = Verr;
}

void FitVar::HoldConstant(bool isTrue) {
    if (isTrue) {
        holdConstant = true;  // This will tell us that it's held constant
        verr = 0.0;  // This will tell Minuit to hold it constant
    } else {
        holdConstant = false;
        verr = verr_copy;
    }
}


//end header guard
#endif