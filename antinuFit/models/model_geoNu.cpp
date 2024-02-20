#include "model_geoNu.hpp"

geoNu::geoNu(const geoNu& mod) {
    Vars = mod.Vars; E_systs = mod.E_systs; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec; model_spec_sys = mod.model_spec_sys;
    histTh = mod.histTh; histU = mod.histU;
    Th_integral = mod.Th_integral; U_integral = mod.U_integral;
    survival_prob = mod.survival_prob;
    computed_survival_prob = mod.computed_survival_prob; E_systs = mod.E_systs;
}

void geoNu::operator = (const geoNu& mod) {
    Vars = mod.Vars; E_systs = mod.E_systs; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec; model_spec_sys = mod.model_spec_sys;
    histTh = mod.histTh; histU = mod.histU;
    Th_integral = mod.Th_integral; U_integral = mod.U_integral;
    survival_prob = mod.survival_prob;
    computed_survival_prob = mod.computed_survival_prob; E_systs = mod.E_systs;
}

geoNu::geoNu(FitVar* NormTh, FitVar* NormU, FitVar* vS_12_2, FitVar* vS_13_2, Esys* E_syst, TH1D* Hist_Th, TH1D* Hist_U) {
    Vars.push_back(vS_12_2); Vars.push_back(vS_13_2);
    Vars.push_back(NormTh); Vars.push_back(NormU);
    E_systs.push_back(E_syst);
    histTh = Hist_Th;
    histU = Hist_U;
    Th_integral = histTh->Integral();
    U_integral = histU->Integral();

    numVars = Vars.size();
    vars.resize(numVars);
    for (unsigned int i = 0; i < Vars.size(); ++i) {
        vars.at(i) = Vars.at(i)->val();
    }
    computed_survival_prob = false;

    model_spec = (TH1D*)(Hist->Clone());
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
        vars.at(0) = Vars.at(0)->val();
        vars.at(1) = Vars.at(1)->val();
        this->geoNu_survival_prob();
        computed_survival_prob = true;
    } else {
        std::cout << "[geoNu] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
    }
}

void geoNu::compute_spec(Double_t* p) {
    this->GetVarValues(p);
    model_spec->Reset("ICES");  // empty it before re-computing it
    model_spec_sys->Reset("ICES");  // empty it before re-computing it

    // If the oscillation parameters are constant and the survival prob was already computed, can skip this step!
    if (!(Vars.at(0)->isConstant() && Vars.at(1)->isConstant() && computed_survival_prob)) {
        this->geoNu_survival_prob();
    }

    model_spec->Add(histTh, survival_prob * vars.at(3) / Th_integral);
    model_spec->Add(histU, survival_prob * vars.at(4) / U_integral);

    // Apply energy systematics
    E_systs.at(0)->apply_systematics(model_spec, model_spec_sys);
}

void geoNu::GetOscHists(TH1D* rescaled_osc_hist_Th, TH1D* rescaled_osc_hist_U) {
    rescaled_osc_hist_Th = (TH1D*)(histTh->Clone());
    rescaled_osc_hist_Th->Reset("ICES");
    rescaled_osc_hist_Th->Add(histTh, survival_prob * vars.at(3) / Th_integral);

    rescaled_osc_hist_U = (TH1D*)(histU->Clone());
    rescaled_osc_hist_U->Reset("ICES");
    rescaled_osc_hist_U->Add(histU, survival_prob * vars.at(4) / U_integral);
}

// Destructor
geoNu::~geoNu() {
    // delete hist;
    // for (auto p : Vars) {delete p;}
    Vars.clear();
}
