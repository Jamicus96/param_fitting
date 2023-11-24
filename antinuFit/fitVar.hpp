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

//end header guard
#endif