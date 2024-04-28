// header guard:
#ifndef alphaN_model
#define alphaN_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "fitVars.hpp"


class alphaN {
    private:
        static alphaN *alphaNInstance_;

        TH1D* model_noEsys;
        TH1D* model_Esys;
        bool isInit;

        unsigned int iNormGS, iNormES, iEsys, iEsysP;
        TH1D* hist_ProtontR;
        TH1D* hist_C12Scatter;
        TH1D* hist_O16Deex;
        double Integral_hist_GS, Integral_hist_ES;
        TH1D* model_Proton;
        unsigned int iMinBin, iMaxBin;
    
    protected:
        // Constructors/desctructors
        alphaN() : isInit(false) {}
        ~alphaN() {}

    public:
        // Should not be cloneable.
        alphaN(alphaN &other) = delete;

        // Should not be assignable.
        void operator=(const alphaN&) = delete;
        /**
         * This is the method that controls the access to the alphaN
         * instance. On the first run, it creates a alphaN object and places it
         * into the field. On subsequent runs, it returns the client existing
         * object stored in the field.
         */
        static alphaN *GetInstance();

        // Initialisers
        void InitAlphaN(const unsigned int NormGS_idx, const unsigned int NormES_idx, const unsigned int E_syst_idx,
               const unsigned int E_syst_proton_idx, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {

            iNormGS = NormGS_idx; iNormES = NormES_idx; iEsys = E_syst_idx; iEsysP = E_syst_proton_idx;
            iMinBin = 1; iMaxBin = Hist_ProtontR->GetXaxis()->GetNbins();

            hist_ProtontR = Hist_ProtontR;
            Integral_hist_GS = hist_ProtontR->Integral(iMinBin, iMaxBin);
            hist_C12Scatter = Hist_C12Scatter;
            Integral_hist_GS += hist_C12Scatter->Integral(iMinBin, iMaxBin);
            hist_O16Deex = Hist_O16Deex;
            Integral_hist_ES = hist_O16Deex->Integral(iMinBin, iMaxBin);

            model_Proton = (TH1D*)(hist_ProtontR->Clone("alphaN::temp_model_PR"));
            model_noEsys = (TH1D*)(hist_ProtontR->Clone("alphaN::model_noEsys"));
            model_Esys = (TH1D*)(hist_ProtontR->Clone("alphaN::model_Esys"));

            isInit = true;
        }
        
        void InitAlphaN(const std::string NormGS_name, const std::string NormES_name, const std::string E_syst_name,
               const std::string E_syst_proton_name, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            InitAlphaN(Vars->findIdx(NormGS_name), Vars->findIdx(NormES_name), Esysts->findIdx(E_syst_name),
                   Esysts->findIdx(E_syst_proton_name), Hist_ProtontR, Hist_C12Scatter, Hist_O16Deex);
        };

        // Member function
        void compute_spec() {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            model_noEsys->Reset("ICES");  // empty it before re-computing it
            model_Proton->Reset("ICES");  // empty it before re-computing it
            model_Esys->Reset("ICES");  // empty it before re-computing it

            #ifdef SUPER_DEBUG
                std::cout << "[alphaN::compute_spec]: hist_C12Scatter->Integral(iMinBin, iMaxBin) = " << hist_C12Scatter->Integral(iMinBin, iMaxBin) << std::endl;
                std::cout << "[alphaN::compute_spec]: hist_O16Deex->Integral(iMinBin, iMaxBin) = " << hist_O16Deex->Integral(iMinBin, iMaxBin) << std::endl;
                std::cout << "[alphaN::compute_spec]: hist_ProtontR->Integral(iMinBin, iMaxBin) = " << hist_ProtontR->Integral(iMinBin, iMaxBin) << std::endl;

                std::cout << "[alphaN::compute_spec]: Vars->val(iNormGS) = " << Vars->val(iNormGS) << std::endl;
                std::cout << "[alphaN::compute_spec]: Vars->val(iNormES) = " << Vars->val(iNormES) << std::endl;

                std::cout << "[alphaN::compute_spec]: Integral_hist_GS = " << Integral_hist_GS << std::endl;
                std::cout << "[alphaN::compute_spec]: Integral_hist_ES = " << Integral_hist_ES << std::endl;
            #endif

            // Add proton recoil to spectrum, and apply extra proton systematics separately
            model_Proton->Add(hist_ProtontR, Vars->val(iNormGS) / Integral_hist_GS);
            Esysts->apply_systematics(iEsysP, model_Proton, model_noEsys);

            // Add the non proton-recoil spectra to spectrum
            model_noEsys->Add(hist_C12Scatter, Vars->val(iNormGS) / Integral_hist_GS);
            model_noEsys->Add(hist_O16Deex, Vars->val(iNormES) / Integral_hist_ES);

            // Apply Normal (beta) energy systematics to all (adds stuff to model_Esys, doesn't reset it)
            Esysts->apply_systematics(iEsys, model_noEsys, model_Esys);  
            #ifdef SUPER_DEBUG
                std::cout << "[alphaN::compute_spec]: model_Esys->Integral(iMinBin, iMaxBin) = " << model_Esys->Integral(iMinBin, iMaxBin) << std::endl;
            #endif
        }

        void Spectra(std::vector<TH1D*>& hists) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            TH1D* temp_hist = (TH1D*)(hist_ProtontR->Clone("temp_hist"));
            
            temp_hist->Reset("ICES");
            temp_hist->Add(hist_ProtontR, Vars->val(iNormGS) / Integral_hist_GS);
            hists.push_back((TH1D*)(hist_ProtontR->Clone("model_alphaN_PR")));
            hists.at(hists.size()-1)->Reset("ICES");
            Esysts->apply_systematics(iEsysP, temp_hist, hists.at(hists.size()-1));

            temp_hist->Reset("ICES");
            temp_hist->Add(hist_C12Scatter, Vars->val(iNormGS) / Integral_hist_GS);
            hists.push_back((TH1D*)(hist_C12Scatter->Clone("model_alphaN_C12")));
            hists.at(hists.size()-1)->Reset("ICES");
            Esysts->apply_systematics(iEsys, temp_hist, hists.at(hists.size()-1));

            temp_hist->Reset("ICES");
            temp_hist->Add(hist_O16Deex, Vars->val(iNormES) / Integral_hist_ES);
            hists.push_back((TH1D*)(hist_O16Deex->Clone("model_alphaN_O16")));
            hists.at(hists.size()-1)->Reset("ICES");
            Esysts->apply_systematics(iEsys, temp_hist, hists.at(hists.size()-1));
        }
    
        TH1D* GetModelNoEsys() {return model_noEsys;}
        TH1D* GetModelEsys() {return model_Esys;}

        bool IsInit() {return isInit;}
        void SetBinLims(const unsigned int MinBin, const unsigned int MaxBin) {
            iMinBin = MinBin;
            iMaxBin = MaxBin;

            Integral_hist_GS = hist_ProtontR->Integral(iMinBin, iMaxBin);
            Integral_hist_GS += hist_C12Scatter->Integral(iMinBin, iMaxBin);
            Integral_hist_ES = hist_O16Deex->Integral(iMinBin, iMaxBin);
        }
};

// Static methods should be defined outside the class.
alphaN* alphaN::alphaNInstance_{nullptr};

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
alphaN* alphaN::GetInstance() {
    if (alphaNInstance_ == nullptr) {
        alphaNInstance_ = new alphaN();
    }
    return alphaNInstance_;
}

//end header guard
#endif