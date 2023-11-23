// header guard:
#ifndef fitter
#define fitter

// include
#include <TVirtualFitter.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TH1.h>
#include <iostream>
#include <algorithm>
#include "fitVar.hpp"
#include "model.hpp"


class Fitter {
    private:
        TVirtualFitter* minuit;
        unsigned int maxNparams = 15;
        double arglist[2];

        std::vector<FitVar*> variables;
        Double_t* var_bestFits;
        std::vector<double> var_bestFit_errs;
        unsigned int numVars = 0;

        std::vector<Model*> models;
        unsigned int numModels = 0;

        TH1D* data, tot_fitModel;
        unsigned int numBins = 0;

    public:
        // Constructors
        Fitter(const unsigned int MaxNparams = 15);
        Fitter(TH1D* Data, std::vector<FitVar*>& Variables, std::vector<Model*>& Models, const unsigned int MaxNparams = 15);
        Fitter(TH1D* Data, const unsigned int MaxNparams = 15);

        // Member function
        double fit_models();
        void fit_models();
        Double_t fitFunc(Double_t* p);
        void fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t) {fval = this->fitFunc(p);};

        double GetFitParameter(const unsigned int parIdx) {return minuit->GetParameter(parIdx)};
        double GetFitParError(const unsigned int parIdx) {return minuit->GetParError(parIdx)};
        double GetFitCovarianceMatrixElement(const unsigned int i, const unsigned int j) {return minuit->GetCovarianceMatrixElement(i, j)};

        double ExtendedConstrainedLogLikelihood(Double_t* p);

        // Destructor
        ~Fitter() {delete arglist; delete var_bestFits;}
};

//end header guard
#endif


/* ~~~~~ Saving other potentially useful commands here ~~~~~ */

// unsigned int NtotParams = minuit->GetNumberTotalParameters();
// unsigned int NfreeParams = minuit->GetNumberFreeParameters();

// double chi2, edm, errdef;
// int nvpar, nparx;
// minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

// void TFitter::GetConfidenceIntervals(Int_t n, Int_t ndim, const Double_t *x, Double_t *ci, Double_t cl)