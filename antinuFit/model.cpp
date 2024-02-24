#include "model.hpp"

Model::Model(const Model& mod) {
    Vars = mod.Vars;
    E_systs = mod.E_systs;
    vars = mod.vars;
    numVars = mod.numVars;
    model_spec = mod.model_spec;
    model_spec_sys = mod.model_spec_sys;
}

Model::Model() {}

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