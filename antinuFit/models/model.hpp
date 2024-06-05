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

        unsigned int numMods;
        std::vector<TH1D*> model_noEsys;
        std::vector<TH1D*> model_Esys;
        std::vector<std::string> names;

        std::vector<unsigned int> iNorm, iE_syst;
        std::vector<TH1D*> hist;
        std::vector<double> integral;
        std::vector<unsigned int> iMinBin, iMaxBin;
    
    protected:
        // Constructors/desctructors
        Model() : numMods(0) {}
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
        void AddModel(const unsigned int Norm_idx, const unsigned int E_syst_idx, TH1D* Hist, const std::string name) {
            ++numMods;

            iNorm.push_back(Norm_idx); iE_syst.push_back(E_syst_idx);
            names.push_back(name);
            iMinBin.push_back(1); iMaxBin.push_back(Hist->GetXaxis()->GetNbins());
            hist.push_back(Hist); hist.at(numMods-1)->SetName((name + "::hist").c_str());
            integral.push_back(hist.at(numMods-1)->Integral(iMinBin.at(numMods-1), iMaxBin.at(numMods-1)));

            model_noEsys.push_back((TH1D*)(hist.at(numMods-1)->Clone("Model::model_noEsys")));
            model_Esys.push_back((TH1D*)(hist.at(numMods-1)->Clone("Model::model_Esys")));
        }
        
        void AddModel(const std::string Norm_name, const std::string E_syst_name, TH1D* Hist, const std::string name) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            AddModel(Vars->findIdx(Norm_name), Esysts->findIdx(E_syst_name), Hist, name);
        }

        unsigned int findIdx(const std::string ModName) {
            for (unsigned int ModIdx = 0; ModIdx < numMods; ++ModIdx) {
                if (ModName == names.at(ModIdx)) return ModIdx;
            }
            std::cout << "[Model::findIdx]: Model '" << ModName << "' not found!" << std::endl;
            exit(1);
            return 0;
        }

        unsigned int GetNumMods() {return numMods;}

        TH1D* GetHist(const unsigned int idx) {return hist.at(idx);}

        void compute_spec(const unsigned int idx) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            model_noEsys.at(idx)->Reset("ICES");  // empty it before re-computing it
            model_Esys.at(idx)->Reset("ICES");  // empty it before re-computing it

            model_noEsys.at(idx)->Add(hist.at(idx), Vars->val(iNorm.at(idx)) / integral.at(idx));

            // Apply energy systematics
            Esysts->apply_systematics(iE_syst.at(idx), model_noEsys.at(idx), model_Esys.at(idx));
        }
        
        void Spectra(const unsigned int idx, std::vector<TH1D*>& hists) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            TH1D* temp_hist = (TH1D*)(hist.at(idx)->Clone("temp_hist"));
            
            temp_hist->Reset("ICES");
            temp_hist->Add(hist.at(idx), Vars->val(iNorm.at(idx)) / integral.at(idx));
            hists.push_back((TH1D*)(hist.at(idx)->Clone(names.at(idx).c_str())));
            hists.at(hists.size()-1)->Reset("ICES");
            Esysts->apply_systematics(iE_syst.at(idx), temp_hist, hists.at(hists.size()-1));
        }

        TH1D* GetModelNoEsys(const unsigned int idx) {return model_noEsys.at(idx);}
        TH1D* GetModelEsys(const unsigned int idx) {return model_Esys.at(idx);}

        void SetBinLims(const unsigned int idx, const unsigned int MinBin, const unsigned int MaxBin) {
            iMinBin.at(idx) = MinBin;
            iMaxBin.at(idx) = MaxBin;

            integral.at(idx) = hist.at(idx)->Integral(iMinBin.at(idx), iMaxBin.at(idx));
        }

        TH1D* GetHist(const std::string ModName) {return GetHist(findIdx(ModName));}
        void compute_spec(const std::string ModName) {compute_spec(findIdx(ModName));}
        void Spectra(const std::string ModName, std::vector<TH1D*>& hists) {Spectra(findIdx(ModName), hists);}
        TH1D* GetModelNoEsys(const std::string ModName) {return GetModelNoEsys(findIdx(ModName));}
        TH1D* GetModelEsys(const std::string ModName) {return GetModelEsys(findIdx(ModName));}
        void SetBinLims(const std::string ModName, const unsigned int MinBin, const unsigned int MaxBin) {SetBinLims(findIdx(ModName), MinBin, MaxBin);}
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