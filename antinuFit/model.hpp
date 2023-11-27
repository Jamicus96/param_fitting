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


Model::Model(const Model& mod) {
    Vars = mod.Vars;
    vars = mod.vars;
    numVars = mod.numVars;
    model_spec = mod.model_spec;
}

Model::Model() {}

void Model::compute_spec(Double_t* p) {
    GetVarValues(p);
    model_spec->Reset("ICES");  // empty it before re-computing it
    compute_spec();
}

void Model::compute_spec() {}
void Model::hold_osc_params_const(bool isTrue) {}

TH1D* Model::Spectrum() {
    return model_spec;
}

void Model::operator = (const Model& mod) {
    Vars = mod.Vars;
    vars = mod.vars;
    numVars = mod.numVars;
    model_spec = mod.model_spec;
}

// Destructor
Model::~Model() {
    for (auto p : Vars) {delete p;} Vars.clear();
}

//end header guard
#endif