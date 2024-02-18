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
        unsigned int iEsys, iEsysP;

    public:
        // Constructors
        alphaN(const alphaN& mod);
        void operator = (const alphaN& mod);
        alphaN(FitVar* NormProtonR, FitVar* NormC12Scatter, FitVar* NormO16Deex, Esys* E_syst, Esys* E_syst_proton, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex);

        // Member function
        void compute_spec(Double_t* p);
    
        // Destructor
        ~alphaN();
};

//end header guard
#endif