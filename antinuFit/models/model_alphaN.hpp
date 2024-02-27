// header guard:
#ifndef alphaN_model
#define alphaN_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "model.hpp"
#include "fitVars.hpp"


class alphaN: public Model {
    private:
        unsigned int iNormGS, iNormES, iEsys, iEsysP;

        TH1D* hist_ProtontR;
        TH1D* hist_C12Scatter;
        TH1D* hist_O16Deex;
        double Integral_hist_GS, Integral_hist_ES;

        TH1D* model_Proton;

    public:
        // Constructors
        alphaN(const unsigned int NormGS_idx, const unsigned int NormES_idx, const unsigned int E_syst_idx,
               const unsigned int E_syst_proton_idx, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex);
        
        alphaN(const std::string NormGS_name, const std::string NormES_name, const std::string E_syst_name,
               const std::string E_syst_proton_name, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {
            alphaN(Vars.findIdx(NormGS_name), Vars.findIdx(NormES_name), Esysts.findIdx(E_syst_name),
                   Esysts.findIdx(E_syst_proton_name), Hist_ProtontR, Hist_C12Scatter, Hist_O16Deex);
        };

        // Member function
        void compute_spec();
        void Spectra(std::vector<TH1D*>& hists);
    
        // Destructor
        ~alphaN() {delete hist_ProtontR; delete hist_C12Scatter; delete hist_O16Deex; delete model_Proton;};
};

//end header guard
#endif