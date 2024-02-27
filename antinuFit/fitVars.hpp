// header guard:
#ifndef fitVars
#define fitVars

// include
#include <iostream>
#include <TVectorD.h>


class FitVars : Fitter {
    private:
        static std::vector<std::string> parnames;
        static std::vector<double> values;
        static std::vector<double> vpriors;
        static std::vector<double> verrs;
        static std::vector<double> verr_copies;
        static std::vector<double> vlows;
        static std::vector<double> vhighs;
        static std::vector<bool> holdConstants;

    public:
        // Constructors
        FitVars() {numVars = 0;};

        // Member function
        static void AddVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool HoldConstant = false) {
            ++numVars;
            parnames.push_back(Parname);
            values.push_back(Value);
            vpriors.push_back(Value);
            verrs.push_back(Verr);
            verr_copies.push_back(Verr);
            vlows.push_back(Vlow);
            vhighs.push_back(Vhigh);
            holdConstants.push_back(HoldConstant);

            if (HoldConstant) verrs.at(numVars-1) = 0.0;
        };

        static unsigned int findIdx(const std::string Parname) {
            for (unsigned int varIdx = 0; varIdx < numVars; ++varIdx) {
                if (Parname == parnames.at(varIdx)) return varIdx;
            }
            std::cout << "[FitVars::findIdx]: Parameter '" << Parname << "' not found!" << std::endl;
            exit(1);
            return 0;
        };

        static void HoldConstant(const unsigned int idx, const bool isTrue) {
            if (isTrue) {
                holdConstant.at(idx) = true;  // This will tell us that it's held constant
                verrs.at(idx) = 0.0;  // This will tell Minuit to hold it constant
            } else {
                holdConstant.at(idx) = false;
                verrs.at(idx) = verr_copy.at(idx);
            }
        };
        static unsigned int GetNumVars() {return numVars;};
        static std::string& name(unsigned int idx) {return parnames.at(idx);};
        static double& val(unsigned int idx) {return values.at(idx);};
        static double& prior(unsigned int idx) {return vpriors.at(idx);};
        static double& err(unsigned int idx) {return verrs.at(idx);};
        static double& min(unsigned int idx) {return vlows.at(idx);};
        static double& max(unsigned int idx) {return vhighs.at(idx);};
        static bool isConstant(unsigned int idx) {return holdConstants.at(idx);};

        static void HoldConstant(const std::string Parname, const bool isTrue) {HoldConstant(findIdx(Parname), isTrue);};
        static double& val(std::string Parname) {return val(findIdx(Parname));};
        static double& prior(std::string Parname) {return prior(findIdx(Parname));};
        static double& err(std::string Parname) {return err(findIdx(Parname));};
        static double& min(std::string Parname) {return min(findIdx(Parname));};
        static double& max(std::string Parname) {return max(findIdx(Parname));};
        static bool isConstant(std::string Parname) {return isConstant(findIdx(Parname));};

        static void GetVarValues(Double_t* p) {
            #ifdef SUPER_DEBUG
                std::cout << "[FitVars::GetVarValues]: numVars = " << numVars << std::endl;
            #endif
            for (unsigned int i = 0; i < numVars; ++i) {
                if (holdConstants.at(i)) continue;  // provided by user
                values.at(i) = p[i];  // provided by Minuit

                #ifdef SUPER_DEBUG
                    std::cout << "[FitVars::GetVarValues]: Vars.at(" << i << "): name = " << Vars.at(i)->name()
                              << ", val = " << Vars.at(i)->val() << ", prior = " << Vars.at(i)->prior() << ", err = " << Vars.at(i)->err()
                              << ", min = " << Vars.at(i)->min() << ", max = " << Vars.at(i)->max() << ", parIdx = " << Vars.at(i)->ParIdx() << std::endl;
                    std::cout << "[FitVars::GetVarValues]: vars.at(" << i << ") = " << vars.at(i) << std::endl;
                #endif
            }
        };

        // Destructor
        ~FitVar() {};
};

//end header guard
#endif