// header guard:
#ifndef model
#define model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "fitVars.hpp"
#include "E_systematics.hpp"


class Model {
    private:
        static Model *ModelInstance_;

        TH1D* model_noEsys;
        TH1D* model_Esys;
        bool isInit;
        std::string ModName;

        unsigned int iNorm, iE_syst;
        TH1D *hist;
        double integral;
        unsigned int iMinBin, iMaxBin;
    
    protected:
        // Constructors/desctructors
        Model() : isInit(false) {}
        ~Model() {}

    public:
        // Should not be cloneable.
        Model(Model &other) = delete;

        // Should not be assignable.
        void operator=(const Model&) = delete;
        /**
         * This is the method that controls the access to the geoNu
         * instance. On the first run, it creates a geoNu object and places it
         * into the field. On subsequent runs, it returns the client existing
         * object stored in the field.
         */
        static Model *GetInstance();

        // Initialisers
        void InitModel(const unsigned int Norm_idx, const unsigned int E_syst_idx, TH1D* Hist, const std::string name) {
            iNorm = Norm_idx; iE_syst = E_syst_idx;
            ModName = name;
            iMinBin = 1; iMaxBin = Hist->GetXaxis()->GetNbins();
            hist = Hist; hist->SetName((ModName + "::hist").c_str());
            integral = hist->Integral(iMinBin, iMaxBin);

            model_noEsys = (TH1D*)(hist->Clone("Model::model_noEsys"));
            model_Esys = (TH1D*)(hist->Clone("Model::model_Esys"));

            isInit = true;
        }
        
        void InitModel(const std::string Norm_name, const std::string E_syst_name, TH1D* Hist, const std::string name) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            InitModel(Vars->findIdx(Norm_name), Esysts->findIdx(E_syst_name), Hist, name);
        }

        TH1D* GetHit() {return hist;}

        void compute_spec() {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            model_noEsys->Reset("ICES");  // empty it before re-computing it
            model_Esys->Reset("ICES");  // empty it before re-computing it

            model_noEsys->Add(hist, Vars->val(iNorm) / integral);

            // Apply energy systematics
            Esysts->apply_systematics(iE_syst, model_noEsys, model_Esys);
        }
        
        void Spectra(std::vector<TH1D*>& hists) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            TH1D* temp_hist = (TH1D*)(hist->Clone("temp_hist"));
            
            temp_hist->Reset("ICES");
            temp_hist->Add(hist, Vars->val(iNorm) / integral);
            hists.push_back((TH1D*)(hist->Clone(ModName.c_str())));
            hists.at(hists.size()-1)->Reset("ICES");
            Esysts->apply_systematics(iE_syst, temp_hist, hists.at(hists.size()-1));
        }

        TH1D* GetModelNoEsys() {return model_noEsys;}
        TH1D* GetModelEsys() {return model_Esys;}

        bool IsInit() {return isInit;}
        void SetBinLims(const unsigned int MinBin, const unsigned int MaxBin) {
            iMinBin = MinBin;
            iMaxBin = MaxBin;

            integral = hist->Integral(iMinBin, iMaxBin);
        }
};

// Static methods should be defined outside the class.
Model* Model::ModelInstance_{nullptr};

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
Model* Model::GetInstance() {
    if (ModelInstance_ == nullptr) {
        ModelInstance_ = new Model();
    }
    return ModelInstance_;
}

//end header guard
#endif