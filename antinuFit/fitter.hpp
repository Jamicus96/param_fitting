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
        static TVirtualFitter* minuit;
        static unsigned int maxNparams;
        static double arglist[2];

        static std::vector<FitVar*> variables;
        static Double_t* var_bestFits;
        static std::vector<double> var_bestFit_errs;
        static unsigned int numVars;

        static std::vector<Model*> models;
        static unsigned int numModels;

        static TH1D* data, *tot_fitModel;
        static unsigned int numBins;

    public:
        // Constructors
        Fitter(const Fitter& fit) {
            minuit = fit.minuit; maxNparams = fit.maxNparams; arglist[0] = fit.arglist[0]; arglist[1] = fit.arglist[1];
            variables = fit.variables; var_bestFits = fit.var_bestFits; var_bestFit_errs = fit.var_bestFit_errs;
            numVars = fit.numVars; models = fit.models; numModels = fit.numModels; data = fit.data;
            tot_fitModel = fit.tot_fitModel; numBins = fit.numBins;
        };
        Fitter(const unsigned int MaxNparams = 15);
        Fitter(TH1D* Data, std::vector<FitVar*>& Variables, std::vector<Model*>& Models, const unsigned int MaxNparams = 15);
        Fitter(TH1D* Data, const unsigned int MaxNparams = 15);

        // Member function
        static double fit_models();
        static Double_t fitFunc(Double_t* p);
        // (Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag*/)
        static void fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t) {fval = fitFunc(p);};
        static void SetFunc() {minuit->SetFCN(fitFunc_minuit_format);};

        static double GetFitParameter(const unsigned int parIdx) {return minuit->GetParameter(parIdx);};
        static double GetFitParError(const unsigned int parIdx) {return minuit->GetParError(parIdx);};
        static double GetFitCovarianceMatrixElement(const unsigned int i, const unsigned int j) {return minuit->GetCovarianceMatrixElement(i, j);};

        static double ExtendedConstrainedLogLikelihood(Double_t* p);

        // Destructor
        ~Fitter() {delete var_bestFits;}
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