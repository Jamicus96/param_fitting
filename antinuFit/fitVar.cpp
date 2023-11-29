#include "fitVar.hpp"



FitVar::FitVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool HoldConstant) {
    parname = Parname;
    value = Value;
    vprior = value;
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