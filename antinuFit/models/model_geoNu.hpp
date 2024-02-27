// header guard:
#ifndef geoNu_model
#define geoNu_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "model.hpp"
#include "fitVars.hpp"
#include "E_systematics.hpp"


class geoNu: public Model {
    private:
        unsigned int iNormTh, iNormU, iS_12_2, iS_13_2, iE_syst;
        TH1D *histTh, *histU;
        double Th_integral, U_integral;
        double survival_prob;
        bool computed_survival_prob;

    public:
        // Constructors
        geoNu(const unsigned int NormTh_idx, const unsigned int NormU_idx, const unsigned int vS_12_2_idx, const unsigned int vS_13_2_idx, const unsigned int E_syst_idx, TH1D* Hist_Th, TH1D* Hist_U);
        geoNu(const std::string NormTh_name, const std::string NormU_name, const std::string vS_12_2_name, const std::string vS_13_2_name, const std::string E_syst_name, TH1D* Hist_Th, TH1D* Hist_U) {
            geoNu(Vars.findIdx(NormTh_name), Vars.findIdx(NormU_name), Vars.findIdx(vS_12_2_name), Vars.findIdx(vS_13_2_name), Esysts.findIdx(E_syst_name), Hist_Th, Hist_U);
        }

        // Member functions
        // Geo-nu survival probability: averaged over baseline -> No E-depence, only depends on mixing angles.
        void geoNu_survival_prob();
        TH1D* GetHit_Th() {return histTh;};
        TH1D* GetHit_U() {return histU;};
        void hold_osc_params_const(const bool isTrue);
        void compute_spec();
        void Spectra(std::vector<TH1D*>& hists);

        // Destructor
        ~geoNu() {delete histTh; delete histU;};
};

//end header guard
#endif