#include "model_alphaN.hpp"


alphaN::alphaN(const alphaN& mod) {
    Vars = mod.Vars; E_systs = mod.E_systs; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec; model_spec_sys = mod.model_spec_sys;
    hist_ProtontR = mod.hist_ProtontR; hist_C12Scatter = mod.hist_C12Scatter;
    hist_O16Deex = mod.hist_O16Deex; Integral_hist_GS = mod.Integral_hist_GS;
    Integral_hist_ES = mod.Integral_hist_ES;
    E_systs = mod.E_systs; model_Proton = mod.model_Proton;
    iEsys = mod.iEsys; iEsysP = mod.iEsysP;
}

void alphaN::operator = (const alphaN& mod) {
    Vars = mod.Vars; E_systs = mod.E_systs; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec; model_spec_sys = mod.model_spec_sys;
    hist_ProtontR = mod.hist_ProtontR; hist_C12Scatter = mod.hist_C12Scatter;
    hist_O16Deex = mod.hist_O16Deex; Integral_hist_GS = mod.Integral_hist_GS;
    Integral_hist_ES = mod.Integral_hist_ES;
    E_systs = mod.E_systs; model_Proton = mod.model_Proton;
    iEsys = mod.iEsys; iEsysP = mod.iEsysP;
}

/**
 * @brief Construct a new alpha N::alpha N object
 * 
 * @param NormGS  Normalisation for events from ground state neutrons (proton recoil + C12 scatter)
 * @param NormES  Normalisation for events from excited state neutrons (O16 de-excitation)
 * @param E_syst  Energy systematics applied to all events
 * @param E_syst_proton  Addition energy systematics applied to proton recoil events
 * @param Hist_ProtontR 
 * @param Hist_C12Scatter 
 * @param Hist_O16Deex 
 */
alphaN::alphaN(FitVar* NormGS, FitVar* NormES, Esys* E_syst, Esys* E_syst_proton, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {
    Vars.push_back(NormGS);
    Vars.push_back(NormES);
    E_systs.push_back(E_syst); iEsys = 0;
    E_systs.push_back(E_syst_proton); iEsysP = 1;

    hist_ProtontR = Hist_ProtontR;
    Integral_hist_GS = hist_ProtontR->Integral();
    hist_C12Scatter = Hist_C12Scatter;
    Integral_hist_GS += hist_C12Scatter->Integral();
    hist_O16Deex = Hist_O16Deex;
    Integral_hist_ES = hist_O16Deex->Integral();

    numVars = Vars.size();
    vars.resize(numVars);
    for (unsigned int i = 0; i < Vars.size(); ++i) {
        vars.at(i) = Vars.at(i)->val();
    }

    model_spec = (TH1D*)(Hist_ProtontR->Clone());
    model_Proton = (TH1D*)(Hist_ProtontR->Clone());
    model_spec_sys = (TH1D*)(Hist_ProtontR->Clone());
}

// Member function
void alphaN::compute_spec(Double_t* p) {
    this->GetVarValues(p);
    model_spec->Reset("ICES");  // empty it before re-computing it
    model_Proton->Reset("ICES");  // empty it before re-computing it
    model_spec_sys->Reset("ICES");  // empty it before re-computing it

    // Add the non proton-recoil spectra to spectrum
    model_spec->Add(hist_C12Scatter, vars.at(0) / Integral_hist_GS);
    model_spec->Add(hist_O16Deex, vars.at(1) / Integral_hist_ES);

    // Apply Normal (beta) energy systematics (adds stuff to model_spec_sys, doesn't reset it)
    E_systs.at(iEsys)->apply_systematics(model_spec, model_spec_sys);

    // Add proton recoil to spectrum, and its own spectrum, to apply proton systematics separately
    model_Proton->Add(hist_ProtontR, vars.at(0) / Integral_hist_GS);
    // model_spec->Add(model_Proton);

    // Apply Proton energy systematics
    E_systs.at(iEsysP)->apply_systematics(model_Proton, model_spec_sys);
}

void alphaN::Spectra(std::vector<TH1D*>& hists) {
    TH1D* temp_hist = (TH1D*)(hist_ProtontR->Clone("temp_hist"));
    
    temp_hist->Reset("ICES");
    temp_hist->Add(hist_ProtontR, vars.at(0) / Integral_hist_GS);
    hists.push_back((TH1D*)(hist_ProtontR->Clone("model_alphaN_PR")));
    hists.at(hists.size()-1)->Reset("ICES");
    E_systs.at(iEsysP)->apply_systematics(temp_hist, hists.at(hists.size()-1));

    temp_hist->Reset("ICES");
    temp_hist->Add(hist_C12Scatter, vars.at(0) / Integral_hist_GS);
    hists.push_back((TH1D*)(hist_C12Scatter->Clone("model_alphaN_C12")));
    hists.at(hists.size()-1)->Reset("ICES");
    E_systs.at(iEsys)->apply_systematics(temp_hist, hists.at(hists.size()-1));

    temp_hist->Reset("ICES");
    temp_hist->Add(hist_O16Deex, vars.at(1) / Integral_hist_ES);
    hists.push_back((TH1D*)(hist_O16Deex->Clone("model_alphaN_O16")));
    hists.at(hists.size()-1)->Reset("ICES");
    E_systs.at(iEsys)->apply_systematics(temp_hist, hists.at(hists.size()-1));
}

// Destructor
alphaN::~alphaN() {
    // delete hist_ProtontR;
    // delete hist_C12Scatter;
    // delete hist_O16Deex;
    // for (auto p : Vars) {delete p;}
    Vars.clear();
}