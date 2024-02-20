// header guard:
#ifndef model
#define model

// include
#include <TVectorD.h>
#include <iostream>
#include <TH1.h>
#include "fitVar.hpp"
#include "E_systematics.hpp"

// #define SUPER_DEBUG

class Model {
    private:

    protected:
        std::vector<FitVar*> Vars;  // sub-classes provide these
        std::vector<Esys*> E_systs;  // sub-classes provide these
        std::vector<double> vars;  // but sub-classes use values from these
        unsigned int numVars;
        TH1D* model_spec;
        TH1D* model_spec_sys;

        void GetVarValues(Double_t* p) {
            #ifdef SUPER_DEBUG
                std::cout << "[Model::GetVarValues]: numVars = " << numVars << std::endl;
            #endif
            for (unsigned int i = 0; i < numVars; ++i) {
                if (Vars.at(i)->isConstant()) vars.at(i) = Vars.at(i)->val();  // provided by user
                else vars.at(i) = p[Vars.at(i)->ParIdx()];  // provided by Minuit

                #ifdef SUPER_DEBUG
                    std::cout << "[Model::GetVarValues]: Vars.at(" << i << "): name = " << Vars.at(i)->name()
                              << ", val = " << Vars.at(i)->val() << ", prior = " << Vars.at(i)->prior() << ", err = " << Vars.at(i)->err()
                              << ", min = " << Vars.at(i)->min() << ", max = " << Vars.at(i)->max() << ", parIdx = " << Vars.at(i)->ParIdx() << std::endl;
                    std::cout << "[Model::GetVarValues]: vars.at(" << i << ") = " << vars.at(i) << std::endl;
                #endif
            }

            for (unsigned int i = 0; i < E_systs.size(); ++i) {
                #ifdef SUPER_DEBUG
                    std::cout << "[Model::GetVarValues]: E_systs.at(" << i << "):" std::endl;
                #endif
                E_systs.at(i)->GetVarValues(p);
            }
        };

    public:
        // Constructors
        Model(const Model& mod);
        Model();

        // Member function
        virtual void compute_spec(Double_t* p);
        virtual void hold_osc_params_const(bool isTrue);
        TH1D* Spectrum();

        virtual void operator = (const Model& mod);

        // Destructor
        virtual ~Model();
};

//end header guard
#endif