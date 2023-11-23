// header guard:
#ifndef geoNu_model
#define geoNu_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "model.hpp"
#include "fitVar.hpp"


class geoNu: public Model {
    private:
        TH1D* hist;
        double hist_integral;
        double survival_prob;
        bool computed_survival_prob = false

    public:
        // Constructors
        geoNu(FitVar* vNorm, FitVar* vS_12_2, FitVar* vS_13_2, TH1D* Hist) {
            Vars.push_back(vS_12_2); Vars.push_back(vS_13_2);
            Vars.push_back(vNorm);
            vars.resize(3);
            hist = Hist;
            hist_integral = hist->Integral();
        };

        // Member function
        // Geo-nu survival probability: averaged over baseline -> No E-depence, only depends on mixing angles.
        void geoNu_survival_prob() {
            const double fSSqrTheta12 = vars.at(0);
            const double fSSqrTheta13 = vars.at(1);
            survival_prob = fSSqrTheta13*fSSqrTheta13 + (1. - fSSqrTheta13)*(1. - fSSqrTheta13) * (1. - 2. * fSSqrTheta12 * (1. - fSSqrTheta12));
        }
        void compute_spec() {
            // If the oscillation parameters are constant and the survival prob was already computed, can skip this step!
            if (!(Vars.at(vS_12_2)->isConstant() && Vars.at(iS_13_2)->isConstant() && computed_survival_prob)) {
                this->geoNu_survival_prob();
                computed_survival_prob = true;
            }

            /* INSERT ENERGY SCALING AND SMEARING HERE */

            model_spec->Add(hist, survival_prob * vars.at(3) / hist_integral);
        };
    
        // Destructor
        ~geoNu() {delete hist;};
};


//end header guard
#endif