// header guard:
#ifndef fitVars
#define fitVars

// include
#include <iostream>
#include <TVectorD.h>


class FitVars {
    private:
        static FitVars *FitVarInstance_;

        unsigned int numVars;
        std::vector<std::string> parnames;
        std::vector<double> values;
        std::vector<double> vpriors;
        std::vector<double> verrs;
        std::vector<double> verr_copies;
        std::vector<double> vlows;
        std::vector<double> vhighs;
        std::vector<bool> holdConstants;
        std::vector<bool> constraineds;
    
    protected:
        // Constructors/desctructors
        FitVars() : numVars(0) {}
        ~FitVars() {}

    public:
        // Should not be cloneable.
        FitVars(FitVars &other) = delete;

        // Should not be assignable.
        void operator=(const FitVars &) = delete;
        /**
         * This is the static method that controls the access to the FitVars
         * instance. On the first run, it creates a FitVars object and places it
         * into the static field. On subsequent runs, it returns the client existing
         * object stored in the static field.
         */
        static FitVars *GetInstance();

        // Member function
        void AddVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool holdConstant = false, bool isConstrained = true) {
            ++numVars;
            parnames.push_back(Parname);
            values.push_back(Value);
            vpriors.push_back(Value);
            verrs.push_back(Verr);
            verr_copies.push_back(Verr);
            vlows.push_back(Vlow);
            vhighs.push_back(Vhigh);
            holdConstants.push_back(holdConstant);
            constraineds.push_back(isConstrained);
        }

        void AddVar_unconstrained(std::string Parname, double Value, double Vlow, double Vhigh) {
            ++numVars;
            parnames.push_back(Parname);
            values.push_back(Value);
            vpriors.push_back(Value);
            verrs.push_back(0.25 * (Vhigh - Vlow));
            verr_copies.push_back(0.25 * (Vhigh - Vlow));
            vlows.push_back(Vlow);
            vhighs.push_back(Vhigh);
            holdConstants.push_back(false);
            constraineds.push_back(false);
        }

        void resetVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool holdConstant = false, bool isConstrained = true) {
            unsigned int idx = findIdx(Parname);

            values.at(idx) = Value;
            vpriors.at(idx) = Value;
            verrs.at(idx) = Verr;
            verr_copies.at(idx) = Verr;
            vlows.at(idx) = Vlow;
            vhighs.at(idx) = Vhigh;
            holdConstants.at(idx) = holdConstant;
            constraineds.at(idx) = isConstrained;
        }

        unsigned int findIdx(const std::string Parname) {
            for (unsigned int varIdx = 0; varIdx < numVars; ++varIdx) {
                if (Parname == parnames.at(varIdx)) return varIdx;
            }
            std::cout << "[FitVars::findIdx]: Parameter '" << Parname << "' not found!" << std::endl;
            exit(1);
            return 0;
        }

        unsigned int GetNumVars() {return numVars;}
        std::string& name(unsigned int idx) {return parnames.at(idx);}
        double& val(unsigned int idx) {return values.at(idx);}
        double& prior(unsigned int idx) {return vpriors.at(idx);}
        double& err(unsigned int idx) {return verrs.at(idx);}
        double& min(unsigned int idx) {return vlows.at(idx);}
        double& max(unsigned int idx) {return vhighs.at(idx);}
        bool IsConstant(unsigned int idx) {return holdConstants.at(idx);}
        bool IsConstrained(unsigned int idx) {return constraineds.at(idx);}
        void HoldConstant(unsigned int idx, bool flag) {holdConstants.at(idx) = flag;}
        void Constrain(unsigned int idx, bool flag) {constraineds.at(idx) = flag;}

        double& val(std::string Parname) {return val(findIdx(Parname));}
        double& prior(std::string Parname) {return prior(findIdx(Parname));}
        double& err(std::string Parname) {return err(findIdx(Parname));}
        double& min(std::string Parname) {return min(findIdx(Parname));}
        double& max(std::string Parname) {return max(findIdx(Parname));}
        bool IsConstant(std::string Parname) {return IsConstant(findIdx(Parname));}
        bool IsConstrained(std::string Parname) {return IsConstrained(findIdx(Parname));}
        void HoldConstant(std::string Parname, bool flag) {HoldConstant(findIdx(Parname), flag);}
        void Constrain(std::string Parname, bool flag) {Constrain(findIdx(Parname), flag);}

        void GetVarValues(Double_t* p) {
            #ifdef SUPER_DEBUG
                std::cout << "[FitVars::GetVarValues]: numVars = " << numVars << std::endl;
            #endif
            for (unsigned int i = 0; i < numVars; ++i) {
                if (holdConstants.at(i)) continue;  // provided by user
                values.at(i) = p[i];  // provided by Minuit

                #ifdef SUPER_DEBUG
                    std::cout << "[FitVars::GetVarValues]: Vars.at(" << i << "): name = " << parnames.at(i)
                              << ", val = " << values.at(i) << ", prior = " << vpriors.at(i) << ", err = " << verrs.at(i)
                              << ", min = " << vlows.at(i) << ", max = " << vhighs.at(i) << std::endl;
                #endif
            }
        }
};

/**
 * Static methods should be defined outside the class.
 */
FitVars* FitVars::FitVarInstance_{nullptr};

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
FitVars *FitVars::GetInstance() {
    if (FitVarInstance_ == nullptr) {
        FitVarInstance_ = new FitVars();
    }
    return FitVarInstance_;
}

//end header guard
#endif