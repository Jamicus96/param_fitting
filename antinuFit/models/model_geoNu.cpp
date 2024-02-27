#include "model_geoNu.hpp"



geoNu::geoNu(const unsigned int NormTh_idx, const unsigned int NormU_idx, const unsigned int vS_12_2_idx,
      const unsigned int vS_13_2_idx, const unsigned int E_syst_idx, TH1D* Hist_Th, TH1D* Hist_U) {

    ModName = "geoNu";
    iNormTh = NormTh_idx; iNormU = NormU_idx; iS_12_2 = vS_12_2_idx; iS_13_2 = vS_13_2_idx; iE_syst = E_syst_idx;
    histTh = Hist_Th; histTh->SetName("geoNu::histTh");
    histU = Hist_U; histU->SetName("geoNu::histU");
    Th_integral = histTh->Integral();
    U_integral = histU->Integral();

    computed_survival_prob = false;

    this->AddModel(ModName, histTh);
}

void geoNu::geoNu_survival_prob() {
    survival_prob = Vars.val(iS_13_2)*Vars.val(iS_13_2) + (1. - Vars.val(iS_13_2))*(1. - Vars.val(iS_13_2)) * (1. - 2. * Vars.val(iS_12_2) * (1. - Vars.val(iS_12_2)));
}

void geoNu::hold_osc_params_const(const bool isTrue) {
    if (!isTrue) {
        computed_survival_prob = false;
    } else if (isTrue && Vars.isConstant(iS_12_2) && Vars.isConstant(iS_13_2)) {
        this->geoNu_survival_prob();
        computed_survival_prob = true;
    } else {
        std::cout << "[geoNu] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
    }
}

void geoNu::compute_spec() {
    model_noEsys->Reset("ICES");  // empty it before re-computing it
    model_Esys->Reset("ICES");  // empty it before re-computing it

    // If the oscillation parameters are constant and the survival prob was already computed, can skip this step!
    if (!(Vars.isConstant(iS_12_2) && Vars.isConstant(iS_13_2) && computed_survival_prob)) {
        this->geoNu_survival_prob();
    }

    model_noEsys->Add(histTh, survival_prob * Vars.val(iNormTh) / Th_integral);
    model_noEsys->Add(histU, survival_prob * Vars.val(iNormU) / U_integral);

    // Apply energy systematics
    Esysts.apply_systematics(iE_syst, model_noEsys, model_Esys);
}

void geoNu::Spectra(std::vector<TH1D*>& hists) {
    TH1D* temp_hist = (TH1D*)(histTh->Clone("temp_hist"));
    
    temp_hist->Reset("ICES");
    temp_hist->Add(histTh, survival_prob * Vars.val(iNormTh) / Th_integral);
    hists.push_back((TH1D*)(histTh->Clone("model_geoNu_Th")));
    hists.at(hists.size()-1)->Reset("ICES");
    Esysts.apply_systematics(iE_syst, temp_hist, hists.at(hists.size()-1));

    temp_hist->Reset("ICES");
    temp_hist->Add(histU, survival_prob * Vars.val(iNormU) / U_integral);
    hists.push_back((TH1D*)(histU->Clone("model_geoNu_U")));
    hists.at(hists.size()-1)->Reset("ICES");
    Esysts.apply_systematics(iE_syst, temp_hist, hists.at(hists.size()-1));
}

