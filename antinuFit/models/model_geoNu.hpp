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
        bool computed_survival_prob;

    public:
        // Constructors
        geoNu(const geoNu& mod) : Model(mod) {
            hist = mod.hist; hist_integral = mod.hist_integral; survival_prob = mod.survival_prob;
            computed_survival_prob = mod.computed_survival_prob;
        };
        geoNu(FitVar* vNorm, FitVar* vS_12_2, FitVar* vS_13_2, TH1D* Hist) {
            Vars.push_back(vS_12_2); Vars.push_back(vS_13_2);
            Vars.push_back(vNorm);
            vars.resize(3);
            hist = Hist;
            hist_integral = hist->Integral();

            for (unsigned int i = 0; i < Vars.size(); ++i) {
                vars.at(i) = Vars.at(i)->val();
            }
            computed_survival_prob = false;
        };

        // Member functions
        // Geo-nu survival probability: averaged over baseline -> No E-depence, only depends on mixing angles.
        void geoNu_survival_prob() {
            const double fSSqrTheta12 = vars.at(0);
            const double fSSqrTheta13 = vars.at(1);
            survival_prob = fSSqrTheta13*fSSqrTheta13 + (1. - fSSqrTheta13)*(1. - fSSqrTheta13) * (1. - 2. * fSSqrTheta12 * (1. - fSSqrTheta12));
        };
        void hold_osc_params_const(bool isTrue) {
            if (!isTrue) {
                computed_survival_prob = false;
            } else if (isTrue && Vars.at(0)->isConstant() && Vars.at(1)->isConstant()) {
                this->geoNu_survival_prob();
                computed_survival_prob = true;
            } else {
                std::cout << "[geoNu] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
            }
        };
        void compute_spec() {
            // If the oscillation parameters are constant and the survival prob was already computed, can skip this step!
            if (!(Vars.at(0)->isConstant() && Vars.at(1)->isConstant() && computed_survival_prob)) {
                this->geoNu_survival_prob();
            }

            /* INSERT ENERGY SCALING AND SMEARING HERE */

            model_spec->Add(hist, survival_prob * vars.at(3) / hist_integral);
        };
        TH1D* GetOscHist() {
            TH1D* rescaled_osc_hist = (TH1D*)(hist->Clone());
            rescaled_osc_hist->Add(hist, survival_prob * vars.at(3) / hist_integral);
            return rescaled_osc_hist;
        };

        // Destructor
        ~geoNu() {delete hist; for (auto p : Vars) {delete p;} Vars.clear();};
};


// friend void hold_osc_params_const(geoNu& Mod, bool isTrue) {
//     if (!isTrue) {
//         Mod.GetCompSurvProbBool() = false;
//     } else if (isTrue && Mod.GetVars().at(0)->isConstant() && Mod.GetVars().at(1)->isConstant()) {
//         Mod.geoNu_survival_prob();
//         Mod.GetCompSurvProbBool() = true;
//     } else {
//         std::cout << "[geoNu] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
//     }
// };


//end header guard
#endif