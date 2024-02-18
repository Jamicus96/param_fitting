#include "model_alphaN.hpp"


alphaN::alphaN(const alphaN& mod) {
    Vars = mod.Vars; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec;
    hist_ProtontR = mod.hist_ProtontR; hist_C12Scatter = mod.hist_C12Scatter;
    hist_O16Deex = mod.hist_O16Deex; Integral_hist_ProtontR = mod.Integral_hist_ProtontR;
    Integral_hist_C12Scatter = mod.Integral_hist_C12Scatter; Integral_hist_O16Deex = mod.Integral_hist_O16Deex;
    E_systs = mod.E_systs;
}

void alphaN::operator = (const alphaN& mod) {
    Vars = mod.Vars; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec;
    hist_ProtontR = mod.hist_ProtontR; hist_C12Scatter = mod.hist_C12Scatter;
    hist_O16Deex = mod.hist_O16Deex; Integral_hist_ProtontR = mod.Integral_hist_ProtontR;
    Integral_hist_C12Scatter = mod.Integral_hist_C12Scatter; Integral_hist_O16Deex = mod.Integral_hist_O16Deex;
    E_systs = mod.E_systs;
}

alphaN::alphaN(FitVar* NormProtonR, FitVar* NormC12Scatter, FitVar* NormO16Deex, Esys* E_syst, Esys* E_syst_proton, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {
    Vars.push_back(NormProtonR);
    Vars.push_back(NormC12Scatter);
    Vars.push_back(NormO16Deex);
    E_systs.push_back(E_syst); iEsys = 0;
    E_systs.push_back(E_syst_proton); iEsysP = 1;

    hist_ProtontR = Hist_ProtontR;
    Integral_hist_ProtontR = hist_ProtontR->Integral();
    hist_C12Scatter = Hist_C12Scatter;
    Integral_hist_C12Scatter = hist_C12Scatter->Integral();
    hist_O16Deex = Hist_O16Deex;
    Integral_hist_O16Deex = hist_O16Deex->Integral();

    numVars = Vars.size();
    vars.resize(numVars);
    for (unsigned int i = 0; i < Vars.size(); ++i) {
        vars.at(i) = Vars.at(i)->val();
    }

    model_spec = (TH1D*)(Hist_ProtontR->Clone());
}

// Member function
void alphaN::compute_spec(Double_t* p) {
    this->GetVarValues(p);
    model_spec->Reset("ICES");  // empty it before re-computing it
    model_spec_sys->Reset("ICES");  // empty it before re-computing it

    // Add proton recoild to spectrum
    model_spec->Add(hist_ProtontR, vars.at(0) / Integral_hist_ProtontR);

    // Apply Proton energy systematics
    E_systs.at(iEsysP)->apply_systematics(p, model_spec, model_spec_sys);

    // Add the rest to spectrum
    model_spec->Add(hist_C12Scatter, vars.at(1) / Integral_hist_C12Scatter);
    model_spec->Add(hist_O16Deex, vars.at(2) / Integral_hist_O16Deex);

    // Apply Normal energy systematics (adds stuff to model_spec_sys, doesn't reset it)
    E_systs.at(iEsys)->apply_systematics(p, model_spec, model_spec_sys);
}

// Destructor
alphaN::~alphaN() {
    // delete hist_ProtontR;
    // delete hist_C12Scatter;
    // delete hist_O16Deex;
    // for (auto p : Vars) {delete p;}
    Vars.clear();
}