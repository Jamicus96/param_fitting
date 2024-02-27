#include "model_alphaN.hpp"


FitVars Fitter::Vars;
Esys Fitter::Esysts;


alphaN::alphaN(const unsigned int NormGS_idx, const unsigned int NormES_idx, const unsigned int E_syst_idx, const unsigned int E_syst_proton_idx, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {

    iNormGS = NormGS_idx; iNormES = NormES_idx; iEsys = E_syst_idx; iEsysP = E_syst_proton_idx;

    hist_ProtontR = Hist_ProtontR;
    Integral_hist_GS = hist_ProtontR->Integral();
    hist_C12Scatter = Hist_C12Scatter;
    Integral_hist_GS += hist_C12Scatter->Integral();
    hist_O16Deex = Hist_O16Deex;
    Integral_hist_ES = hist_O16Deex->Integral();

    model_Proton = (TH1D*)(Hist_ProtontR->Clone("alphaN::temp_model_PR"));

    this->AddModel("alphaN", Hist_ProtontR);
}

// Member function
void alphaN::compute_spec() {
    models_noEsys.at(modIdx)->Reset("ICES");  // empty it before re-computing it
    model_Proton->Reset("ICES");  // empty it before re-computing it
    models_Esys.at(modIdx)->Reset("ICES");  // empty it before re-computing it

    // Add the non proton-recoil spectra to spectrum
    models_noEsys.at(modIdx)->Add(hist_C12Scatter, Vars.val(iNormGS) / Integral_hist_GS);
    models_noEsys.at(modIdx)->Add(hist_O16Deex, Vars.val(iNormES) / Integral_hist_ES);

    // Apply Normal (beta) energy systematics (adds stuff to model_spec_sys, doesn't reset it)
    Esysts.apply_systematics(iEsys, models_noEsys.at(modIdx), models_Esys.at(modIdx));

    // Add proton recoil to spectrum, and its own spectrum, to apply proton systematics separately
    model_Proton->Add(hist_ProtontR, Vars.val(iNormGS) / Integral_hist_GS);
    // models_noEsys.at(modIdx)->Add(model_Proton);

    // Apply Proton energy systematics
    Esysts.apply_systematics(iEsysP, model_Proton, models_Esys.at(modIdx));
}

void alphaN::Spectra(std::vector<TH1D*>& hists) {
    TH1D* temp_hist = (TH1D*)(hist_ProtontR->Clone("temp_hist"));
    
    temp_hist->Reset("ICES");
    temp_hist->Add(hist_ProtontR, Vars.val(iNormGS) / Integral_hist_GS);
    hists.push_back((TH1D*)(hist_ProtontR->Clone("model_alphaN_PR")));
    hists.at(hists.size()-1)->Reset("ICES");
    Esysts.apply_systematics(iEsysP, temp_hist, hists.at(hists.size()-1));

    temp_hist->Reset("ICES");
    temp_hist->Add(hist_C12Scatter, Vars.val(iNormGS) / Integral_hist_GS);
    hists.push_back((TH1D*)(hist_C12Scatter->Clone("model_alphaN_C12")));
    hists.at(hists.size()-1)->Reset("ICES");
    Esysts.apply_systematics(iEsys, temp_hist, hists.at(hists.size()-1));

    temp_hist->Reset("ICES");
    temp_hist->Add(hist_O16Deex, Vars.val(iNormES) / Integral_hist_ES);
    hists.push_back((TH1D*)(hist_O16Deex->Clone("model_alphaN_O16")));
    hists.at(hists.size()-1)->Reset("ICES");
    Esysts.apply_systematics(iEsys, temp_hist, hists.at(hists.size()-1));
}
