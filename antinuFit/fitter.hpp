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
        inline static TVirtualFitter* minuit;
        inline static double arglist[2];

        inline static std::vector<FitVar*> variables;
        inline static Double_t* var_bestFits;
        inline static double* var_bestFit_errs;
        inline static unsigned int numVars;

        inline static std::vector<Model*> models;
        inline static unsigned int numModels;

        inline static TH1D* data, *tot_fitModel;
        inline static unsigned int numBins;

    public:
        // Constructors
        // Fitter(const Fitter& fit) {
        //     minuit = fit.minuit; maxNparams = fit.maxNparams; arglist[0] = fit.arglist[0]; arglist[1] = fit.arglist[1];
        //     variables = fit.variables; var_bestFits = fit.var_bestFits; var_bestFit_errs = fit.var_bestFit_errs;
        //     numVars = fit.numVars; models = fit.models; numModels = fit.numModels; data = fit.data;
        //     tot_fitModel = fit.tot_fitModel; numBins = fit.numBins;
        // };
        Fitter();
        Fitter(TH1D* Data, std::vector<FitVar*>& Variables, std::vector<Model*>& Models);
        Fitter(TH1D* Data);

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
        ~Fitter() {};
};


Fitter::Fitter(TH1D* Data, std::vector<FitVar*>& Variables, std::vector<Model*>& Models) : Fitter::Fitter(Data) {

    numVars = Variables.size();
    variables.resize(numVars);
    var_bestFit_errs = new double[numVars];

    for (unsigned int i = 0; i < numVars; ++i) {
        variables.at(i) = Variables.at(i);
        // Sets the parameter
        minuit->SetParameter(i, variables.at(i)->name().c_str(), variables.at(i)->prior(), variables.at(i)->err(), variables.at(i)->min(), variables.at(i)->max());
        // Binds the parameter's detail's addresses to fitVar's, so that changing the fitVar object changes this automatically
        minuit->GetParameter(i, const_cast<char*>(variables.at(i)->name().c_str()), variables.at(i)->prior(), variables.at(i)->err(), variables.at(i)->min(), variables.at(i)->max());
        variables.at(i)->ParIdx() = i;
    }

    numModels = Models.size();
    models.resize(numModels);
    for (unsigned int i = 0; i < numModels; ++i) {
        models.at(i) = Models.at(i);
    }
}

Fitter::Fitter(TH1D* Data) : Fitter::Fitter() {

    data = Data;
    numBins = data->GetXaxis()->GetNbins();
    tot_fitModel = (TH1D*)(data->Clone());
    tot_fitModel->SetName("tot_fit_model");
    tot_fitModel->SetTitle("Total Fit Spectrum");
}

Fitter::Fitter() {

    var_bestFits = new Double_t[numVars];

    //The default minimizer is Minuit, you can also try Minuit2
    TVirtualFitter::SetDefaultFitter("Minuit2");
    // TVirtualFitter::SetDefaultFitter("Minuit");
    minuit = TVirtualFitter::Fitter(0, numVars);
    
    arglist[0] = 0;
    // set print level
    minuit->ExecuteCommand("SET PRINT", arglist, 2);

    // minimiser settings
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance

    // (void (*fcn)(Int_t&, Double_t*, Double_t& f, Double_t*, Int_t))
    SetFunc();
}

Double_t Fitter::fitFunc(Double_t* p) {
    // Reset total model spectrum to zero
    tot_fitModel->Reset("ICES");

    // Compute the total spectrum from all the models, using parameters p
    for (unsigned int iModel = 0; iModel < numModels; ++iModel) {
        models.at(iModel)->compute_spec(p);
        tot_fitModel->Add(models.at(iModel)->Spectrum());
    }

    // Compute test statistic
    return ExtendedConstrainedLogLikelihood(p);
}

double Fitter::ExtendedConstrainedLogLikelihood(Double_t* p) {
    // Model PDFs have already been scaled by their respective norms, and added together
    double logL = - tot_fitModel->Integral();
    for (unsigned int ibin = 1; ibin < numBins+1; ++ibin) {
        logL += data->GetBinContent(ibin) * log(tot_fitModel->GetBinContent(ibin));
    }
    // Add in constraints
    for (unsigned int iVar = 0; iVar < numVars; ++ iVar) {
        // If variable is being held constant, assume it's equal to its prior most likely value
        if (!(variables.at(iVar)->isConstant())) {
            // Add gaussian constraint
            logL += 0.5 * (p[iVar] - variables.at(iVar)->prior())*(p[iVar] - variables.at(iVar)->prior()) / (variables.at(iVar)->err()*variables.at(iVar)->err());
        }
    }
    return logL;
}

double Fitter::fit_models() {
    // minimize
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    //get result
    return fitFunc(var_bestFits);
}


//end header guard
#endif


/* ~~~~~ Saving other potentially useful commands here ~~~~~ */

// unsigned int NtotParams = minuit->GetNumberTotalParameters();
// unsigned int NfreeParams = minuit->GetNumberFreeParameters();

// double chi2, edm, errdef;
// int nvpar, nparx;
// minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

// void TFitter::GetConfidenceIntervals(Int_t n, Int_t ndim, const Double_t *x, Double_t *ci, Double_t cl)