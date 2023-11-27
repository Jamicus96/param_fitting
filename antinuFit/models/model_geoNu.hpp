// header guard:
#ifndef geoNu_model
#define geoNu_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "model.hpp"
#include "fitVar.hpp"


class geoNu: public Model {
    private:
        TH1D* hist;
        double hist_integral;
        double survival_prob;
        bool computed_survival_prob;

    public:
        // Constructors
        geoNu(const geoNu& mod);
        void operator = (const geoNu& mod);
        geoNu(FitVar* vNorm, FitVar* vS_12_2, FitVar* vS_13_2, TH1D* Hist);

        // Member functions
        // Geo-nu survival probability: averaged over baseline -> No E-depence, only depends on mixing angles.
        void geoNu_survival_prob();
        void hold_osc_params_const(bool isTrue);
        void compute_spec();
        TH1D* GetOscHist();

        // Destructor
        ~geoNu();
};

//end header guard
#endif