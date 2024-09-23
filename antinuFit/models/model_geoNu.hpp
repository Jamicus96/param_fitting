// header guard:
#ifndef geoNu_model
#define geoNu_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "fitVars.hpp"
#include "E_systematics.hpp"


class geoNu {
    private:
        static geoNu *geoNuInstance_;

        TH1D* model_noEsys;
        TH1D* model_Esys;
        bool isInit;

        unsigned int iNorm, iUThRatio, iS_12_2, iS_13_2, iE_syst;
        TH1D *histTh, *histU;
        double Th_integral, U_integral;
        double survival_prob;
        bool computed_survival_prob;
        unsigned int iMinBin, iMaxBin;
    
    protected:
        // Constructors/desctructors
        geoNu() : isInit(false) {}
        ~geoNu() {}

    public:
        // Should not be cloneable.
        geoNu(geoNu &other) = delete;

        // Should not be assignable.
        void operator=(const geoNu&) = delete;
        /**
         * This is the method that controls the access to the geoNu
         * instance. On the first run, it creates a geoNu object and places it
         * into the field. On subsequent runs, it returns the client existing
         * object stored in the field.
         */
        static geoNu *GetInstance();

        // Initialisers
        void InitGeoNu(const unsigned int Norm_idx, const unsigned int UThRatio_idx, const unsigned int vS_12_2_idx, const unsigned int vS_13_2_idx, const unsigned int E_syst_idx, TH1D* Hist_Th, TH1D* Hist_U) {
            iNorm = Norm_idx; iUThRatio = UThRatio_idx; iS_12_2 = vS_12_2_idx; iS_13_2 = vS_13_2_idx; iE_syst = E_syst_idx;
            iMinBin = 1; iMaxBin = Hist_Th->GetXaxis()->GetNbins();
            histTh = Hist_Th; histTh->SetName("geoNu::histTh");
            histU = Hist_U; histU->SetName("geoNu::histU");
            Th_integral = histTh->Integral(iMinBin, iMaxBin);
            U_integral = histU->Integral(iMinBin, iMaxBin);

            computed_survival_prob = false;

            model_noEsys = (TH1D*)(histTh->Clone("geoNu::model_noEsys"));
            model_Esys = (TH1D*)(histTh->Clone("geoNu::model_Esys"));

            isInit = true;
        }
        
        void InitGeoNu(const std::string Norm_name, const std::string UThRatio_name, const std::string vS_12_2_name, const std::string vS_13_2_name, const std::string E_syst_name, TH1D* Hist_Th, TH1D* Hist_U) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            InitGeoNu(Vars->findIdx(Norm_name), Vars->findIdx(UThRatio_name), Vars->findIdx(vS_12_2_name), Vars->findIdx(vS_13_2_name), Esysts->findIdx(E_syst_name), Hist_Th, Hist_U);
        }

        // Member functions
        // Geo-nu survival probability: averaged over baseline -> No E-depence, only depends on mixing angles.
        void geoNu_survival_prob() {
            FitVars* Vars = FitVars::GetInstance();
            survival_prob = Vars->val(iS_13_2)*Vars->val(iS_13_2) + (1. - Vars->val(iS_13_2))*(1. - Vars->val(iS_13_2)) * (1. - 2. * Vars->val(iS_12_2) * (1. - Vars->val(iS_12_2)));
        }

        TH1D* GetHist_Th() {return histTh;}
        TH1D* GetHist_U() {return histU;}

        void hold_osc_params_const(const bool isTrue) {
            FitVars* Vars = FitVars::GetInstance();
            if (!isTrue) {
                computed_survival_prob = false;
            } else if (isTrue && Vars->IsConstant(iS_12_2) && Vars->IsConstant(iS_13_2)) {
                geoNu_survival_prob();
                computed_survival_prob = true;
            } else {
                std::cout << "[geoNu] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
            }
        }

        void compute_spec() {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            model_noEsys->Reset("ICES");  // empty it before re-computing it
            model_Esys->Reset("ICES");  // empty it before re-computing it

            // If the oscillation parameters are constant and the survival prob was already computed, can skip this step!
            if (!(Vars->IsConstant(iS_12_2) && Vars->IsConstant(iS_13_2) && computed_survival_prob)) {
                geoNu_survival_prob();
            }
            double N_Th = Vars->val(iNorm) / (1.0 + Vars->val(iUThRatio));
            double N_U = Vars->val(iNorm) * Vars->val(iUThRatio) / (1.0 + Vars->val(iUThRatio));

            model_noEsys->Add(histTh, survival_prob * N_Th / Th_integral);
            model_noEsys->Add(histU, survival_prob * N_U / U_integral);

            // Apply energy systematics
            Esysts->apply_systematics(iE_syst, model_noEsys, model_Esys);
        }
        
        void Spectra(std::vector<TH1D*>& hists) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            TH1D* temp_hist = (TH1D*)(histTh->Clone("temp_hist"));

            double N_Th = Vars->val(iNorm) / (1.0 + Vars->val(iUThRatio));
            double N_U = Vars->val(iNorm) * Vars->val(iUThRatio) / (1.0 + Vars->val(iUThRatio));
            
            temp_hist->Reset("ICES");
            temp_hist->Add(histTh, survival_prob * N_Th / Th_integral);
            hists.push_back((TH1D*)(histTh->Clone("model_geoNu_Th")));
            hists.at(hists.size()-1)->Reset("ICES");
            Esysts->apply_systematics(iE_syst, temp_hist, hists.at(hists.size()-1));

            temp_hist->Reset("ICES");
            temp_hist->Add(histU, survival_prob * N_U / U_integral);
            hists.push_back((TH1D*)(histU->Clone("model_geoNu_U")));
            hists.at(hists.size()-1)->Reset("ICES");
            Esysts->apply_systematics(iE_syst, temp_hist, hists.at(hists.size()-1));
        }

        TH1D* GetModelNoEsys() {return model_noEsys;}
        TH1D* GetModelEsys() {return model_Esys;}

        bool IsInit() {return isInit;}
        void SetBinLims(const unsigned int MinBin, const unsigned int MaxBin) {
            iMinBin = MinBin;
            iMaxBin = MaxBin;

            Th_integral = histTh->Integral(iMinBin, iMaxBin);
            U_integral = histU->Integral(iMinBin, iMaxBin);
        }
};

// Static methods should be defined outside the class.
geoNu* geoNu::geoNuInstance_{nullptr};

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
geoNu* geoNu::GetInstance() {
    if (geoNuInstance_ == nullptr) {
        geoNuInstance_ = new geoNu();
    }
    return geoNuInstance_;
}

//end header guard
#endif