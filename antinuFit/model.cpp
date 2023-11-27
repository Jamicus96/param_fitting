#include "model.hpp"

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