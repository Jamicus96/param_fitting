// header guard:
#ifndef alphaN_model
#define alphaN_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "model.hpp"
#include "fitVar.hpp"


class alphaN: public Model {
    private:
        TH1D* hist_ProtontR;
        TH1D* hist_C12Scatter;
        TH1D* hist_O16Deex;
        double Integral_hist_ProtontR, Integral_hist_C12Scatter, Integral_hist_O16Deex;

    public:
        // Constructors
        alphaN(const alphaN& mod) : Model(mod) {
            hist_ProtontR = mod.hist_ProtontR; hist_C12Scatter = mod.hist_C12Scatter;
            hist_O16Deex = mod.hist_O16Deex; Integral_hist_ProtontR = mod.Integral_hist_ProtontR;
            Integral_hist_C12Scatter = mod.Integral_hist_C12Scatter; Integral_hist_O16Deex = mod.Integral_hist_O16Deex;
        };
        alphaN(FitVar* NormProtonR, FitVar* NormC12Scatter, FitVar* NormO16Deex, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {
            Vars.push_back(NormProtonR);
            Vars.push_back(NormC12Scatter);
            Vars.push_back(NormO16Deex);
            vars.resize(3);

            hist_ProtontR = Hist_ProtontR;
            Integral_hist_ProtontR = hist_ProtontR->Integral();
            hist_C12Scatter = Hist_C12Scatter;
            Integral_hist_C12Scatter = hist_C12Scatter->Integral();
            hist_O16Deex = Hist_O16Deex;
            Integral_hist_O16Deex = hist_O16Deex->Integral();

            for (unsigned int i = 0; i < Vars.size(); ++i) {
                vars.at(i) = Vars.at(i)->val();
            }
        };

        // Member function
        void compute_spec() {
            /* INSERT ENERGY SCALING AND SMEARING HERE */

            model_spec->Add(hist_ProtontR, vars.at(0) / Integral_hist_ProtontR);
            model_spec->Add(hist_C12Scatter, vars.at(1) / Integral_hist_C12Scatter);
            model_spec->Add(hist_O16Deex, vars.at(2) / Integral_hist_O16Deex);
        };
    
        // Destructor
        ~alphaN(){delete hist_ProtontR; delete hist_C12Scatter; delete hist_O16Deex; for (auto p : Vars) {delete p;} Vars.clear();};
};


//end header guard
#endif