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
        geoNu(const geoNu& mod);
        void operator = (const geoNu& mod);
        geoNu(FitVar* vNorm, FitVar* vS_12_2, FitVar* vS_13_2, TH1D* Hist);

        // Member functions
        // Geo-nu survival probability: averaged over baseline -> No E-depence, only depends on mixing angles.
        void geoNu_survival_prob();
        void hold_osc_params_const(bool isTrue);
        void compute_spec();
        TH1D* GetOscHist();

        // Destructor
        ~geoNu();
};

geoNu::geoNu(const geoNu& mod) {
    Vars = mod.Vars; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec;
    hist = mod.hist; hist_integral = mod.hist_integral; survival_prob = mod.survival_prob;
    computed_survival_prob = mod.computed_survival_prob;
}

void geoNu::operator = (const geoNu& mod) {
    Vars = mod.Vars; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec;
    hist = mod.hist; hist_integral = mod.hist_integral; survival_prob = mod.survival_prob;
    computed_survival_prob = mod.computed_survival_prob;
}

geoNu::geoNu(FitVar* vNorm, FitVar* vS_12_2, FitVar* vS_13_2, TH1D* Hist) {
    Vars.push_back(vS_12_2); Vars.push_back(vS_13_2);
    Vars.push_back(vNorm);
    vars.resize(3);
    hist = Hist;
    hist_integral = hist->Integral();

    for (unsigned int i = 0; i < Vars.size(); ++i) {
        vars.at(i) = Vars.at(i)->val();
    }
    computed_survival_prob = false;
}

void geoNu::geoNu_survival_prob() {
    const double fSSqrTheta12 = vars.at(0);
    const double fSSqrTheta13 = vars.at(1);
    survival_prob = fSSqrTheta13*fSSqrTheta13 + (1. - fSSqrTheta13)*(1. - fSSqrTheta13) * (1. - 2. * fSSqrTheta12 * (1. - fSSqrTheta12));
}

void geoNu::hold_osc_params_const(bool isTrue) {
    if (!isTrue) {
        computed_survival_prob = false;
    } else if (isTrue && Vars.at(0)->isConstant() && Vars.at(1)->isConstant()) {
        this->geoNu_survival_prob();
        computed_survival_prob = true;
    } else {
        std::cout << "[geoNu] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
    }
}

void geoNu::compute_spec() {
    // If the oscillation parameters are constant and the survival prob was already computed, can skip this step!
    if (!(Vars.at(0)->isConstant() && Vars.at(1)->isConstant() && computed_survival_prob)) {
        this->geoNu_survival_prob();
    }

    /* INSERT ENERGY SCALING AND SMEARING HERE */

    model_spec->Add(hist, survival_prob * vars.at(3) / hist_integral);
}

TH1D* geoNu::GetOscHist() {
    TH1D* rescaled_osc_hist = (TH1D*)(hist->Clone());
    rescaled_osc_hist->Add(hist, survival_prob * vars.at(3) / hist_integral);
    return rescaled_osc_hist;
}

// Destructor
geoNu::~geoNu() {
    delete hist;
    for (auto p : Vars) {delete p;}
    Vars.clear();
}


//end header guard
#endif