// header guard:
#ifndef model
#define model

// include
#include <TVectorD.h>
#include <iostream>
#include <TH1.h>
#include "fitVar.hpp"


class Model {
    private:
        void GetVarValues(Double_t* p) {
            for (unsigned int i = 0; i < numVars; ++i) {
                if (Vars.at(i)->isConstant()) vars.at(i) = Vars.at(i)->val();  // provided by user
                else vars.at(i) = p[Vars.at(i)->ParIdx()];  // provided by Minuit
            }
        };

    protected:
        std::vector<FitVar*> Vars;  // sub-classes provide these
        std::vector<double> vars;  // but sub-classes use values from these
        unsigned int numVars;
        TH1D* model_spec;

    public:
        // Constructors
        Model(const Model& mod);
        Model();

        // Member function
        void compute_spec(Double_t* p);

        virtual void compute_spec();
        virtual void hold_osc_params_const(bool isTrue);
        TH1D* Spectrum();

        virtual void operator = (const Model& mod);

        // Destructor
        virtual ~Model();
};

//end header guard
#endif