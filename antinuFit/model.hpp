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
        Model(const Model& mod) {Vars = mod.Vars; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec;};
        Model(/* FitVar* Norm, FitVar* Var1, FitVar* Var2, const std::vector<double>& X_vals */) {/* numX = x_vals.size(); model_spec.resize(numX) */};

        // Member function
        void compute_spec(Double_t* p) {
            GetVarValues(p);
            model_spec->Reset("ICES");  // empty it before re-computing it
            compute_spec();
        };

        virtual void compute_spec() {/* model_spec.at(0) = vars.at(1); */}
        virtual void hold_osc_params_const(bool isTrue) {};
        TH1D* Spectrum() {return model_spec;};
        // Destructor
        ~Model() {for (auto p : Vars) {delete p;} Vars.clear();};
};

//end header guard
#endif