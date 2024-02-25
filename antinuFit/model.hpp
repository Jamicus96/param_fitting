// header guard:
#ifndef model
#define model

// include
#include <TVectorD.h>
#include <iostream>
#include <TH1.h>
#include "fitter.hpp"

// #define SUPER_DEBUG

class Model : Fitter {
    private:

    protected:
        unsigned int model_idx;

    public:
        // Constructors
        Model(const Model& mod) {model_idx = mod.model_idx;};
        Model() {};
        Model(unsigned int ModIdx) {model_idx = ModIdx;};

        // Member function
        void SetModIdx(unsigned int ModIdx) {model_idx = ModIdx;};
        virtual void compute_spec(Double_t* p) {};
        virtual void hold_osc_params_const(bool isTrue) {};
        void Spectra(std::vector<TH1D*>& hists) {};

        virtual void operator = (const Model& mod) {model_idx = mod.model_idx;};

        // Destructor
        virtual ~Model() {};
};

//end header guard
#endif