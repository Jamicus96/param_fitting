// header guard:
#ifndef model
#define model

// include
#include <TVectorD.h>
#include <iostream>
#include <TH1.h>
#include "fitVars.hpp"
#include "E_systematics.hpp"

// #define SUPER_DEBUG

class Model {
    private:

    protected:
        std::string ModName;
        TH1D* model_noEsys;
        TH1D* model_Esys;
        FitVars Vars;
        Esys Esysts;

    public:
        // Constructors
        Model() {};

        // Member function
        void AddModel(std::string modName, TH1D* templateHist) {
            model_noEsys = (TH1D*)(templateHist->Clone((modName + "::model_noEsys").c_str()));
            model_Esys = (TH1D*)(templateHist->Clone((modName + "::model_Esys").c_str()));
        };
        virtual void compute_spec() {};
        virtual void hold_osc_params_const(bool isTrue) {};
        virtual void Spectra(std::vector<TH1D*>& hists) {};
        TH1D* GetModelNoEsys() {return model_noEsys;};
        TH1D* GetModelEsys() {return model_noEsys;};
        std::string& GetName() {return ModName;};

        // Destructor
       virtual ~Model() {delete model_noEsys; delete model_Esys;};
};

//end header guard
#endif