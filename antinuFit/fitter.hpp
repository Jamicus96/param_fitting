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
#include "fitVars.hpp"
#include "model_Reactor.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"


class Fitter {
    private:
        static Fitter *FitterInstance_;

        static TVirtualFitter* minuit;
        static double arglist[2];

        static TH1D *data, *tot_fitModel;
        static unsigned int numBins;
        static bool dataIsSet;

        // For minuit
        static double minfuncOut, edm, errdef;
        static int nvpar, nparx;

    
    protected:
        // Constructors/desctructors
        Fitter() {
            #ifdef antinuDEBUG
                std::cout << "[Fitter::Fitter]: Creating singleton." << std::endl;
            #endif

            FitVars* Vars = FitVars::GetInstance();
            if (Vars->GetNumVars() == 0) {
                std::cout << "Initialising miniut with zero variables! Must set up all variables first." << std::endl;
                exit(1);
            }

            //The default minimizer is Minuit, you can also try Minuit2
            // TVirtualFitter::SetDefaultFitter("Minuit2");
            TVirtualFitter::SetDefaultFitter("Minuit");
            minuit = TVirtualFitter::Fitter(0, Vars->GetNumVars());
            
            arglist[0] = 0;
            // set print level
            minuit->ExecuteCommand("SET PRINT", arglist, 2);

            // minimiser settings
            arglist[0] = 5000; // number of function calls
            arglist[1] = 0.01; // tolerance

            // (void (*fcn)(Int_t&, Double_t*, Double_t& f, Double_t*, Int_t))
            SetFunc();

            for (unsigned int iVar = 0; iVar < Vars->GetNumVars(); ++iVar) {
                // Sets the parameter
                minuit->SetParameter(iVar, Vars->name(iVar).c_str(), Vars->prior(iVar), Vars->err(iVar), Vars->min(iVar), Vars->max(iVar));
                // Binds the parameter's detail's addresses to fitVar's, so that changing the fitVar object changes this automatically
                minuit->GetParameter(iVar, const_cast<char*>(Vars->name(iVar).c_str()), Vars->prior(iVar), Vars->err(iVar), Vars->min(iVar), Vars->max(iVar));
            }
        }
        ~Fitter() {}

    public:
        // Should not be cloneable.
        Fitter(Fitter &other) = delete;

        // Should not be assignable.
        void operator=(const Fitter&) = delete;
        /**
         * This is the method that controls the access to the Fitter
         * instance. On the first run, it creates a Fitter object and places it
         * into the field. On subsequent runs, it returns the client existing
         * object stored in the field.
         */
        static Fitter *GetInstance();

        // Fitting functions
        static double fit_models();
        static Double_t fitFunc(Double_t* p);

        // (Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag*/)
        static void fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t);
        static void SetFunc();

        double GetFitParameter(const unsigned int parIdx) {return minuit->GetParameter(parIdx);}
        double GetFitParError(const unsigned int parIdx) {return minuit->GetParError(parIdx);}
        double GetFitCovarianceMatrixElement(const unsigned int i, const unsigned int j) {return minuit->GetCovarianceMatrixElement(i, j);}
        static double ExtendedConstrainedLogLikelihood();

        // Get info functions
        static void GetAllSpectra(std::vector<TH1D*>& hists);
        static TH1D* DataHist();
        static void SetData(TH1D* Data);
};

// Static methods should be defined outside the class.
Fitter* Fitter::FitterInstance_{nullptr};

TVirtualFitter* Fitter::minuit;
double Fitter::arglist[2];
TH1D *Fitter::data, *Fitter::tot_fitModel;
unsigned int Fitter::numBins;
double Fitter::minfuncOut, Fitter::edm, Fitter::errdef;
int Fitter::nvpar, Fitter::nparx;
bool Fitter::dataIsSet = false;

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
Fitter* Fitter::GetInstance() {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::GetInstance]: Creating Instance." << std::endl;
    #endif
    if (FitterInstance_ == nullptr) {
        FitterInstance_ = new Fitter();
    }
    return FitterInstance_;
}

void Fitter::fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t) {
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::fitFunc_minuit_format]: Running fitFunc_minuit_format()." << std::endl;
    #endif
    fval = fitFunc(p);
}

void Fitter::SetFunc() {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::SetFunc]: Bound minuit to fitFunc_minuit_format()." << std::endl;
    #endif
    minuit->SetFCN(fitFunc_minuit_format);
}

double Fitter::fit_models() {
    if (!dataIsSet) {
        std::cout << "Attempting to perform fit with no data!" << std::endl;
        exit(1);
    }

    // minimize
    #ifdef antinuDEBUG
        FitVars* Vars = FitVars::GetInstance();
        std::cout << "[Fitter::fit_models]: Running MIGRAD!" << std::endl;
        std::cout << "[Fitter::fit_models]: s_12^2 = " << Vars->val("sinsqrtheta12") << std::endl;
        std::cout << "[Fitter::fit_models]: Dm21^2 = " << Vars->val("deltamsqr21") << std::endl;
    #endif

    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    //get result
    minuit->GetStats(minfuncOut, edm, errdef, nvpar, nparx);
    std::cout << "[Fitter::fit_models]: minfuncOut = " << minfuncOut << std::endl;
    return minfuncOut;
}

Double_t Fitter::fitFunc(Double_t* p) {
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::fitFunc]: Running fitFunc()." << std::endl;
    #endif
    // Reset total model spectrum to zero
    tot_fitModel->Reset("ICES");

    // Compute the total spectrum from all the models, using parameters p
    FitVars* Vars = FitVars::GetInstance();
    Reactor* ReactorMod = Reactor::GetInstance();
    alphaN* alphaNMod = alphaN::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();

    if (!(ReactorMod->IsInit() || alphaNMod->IsInit() || geoNuMod->IsInit())) {
        std::cout << "Attempting to perform fit with no models!" << std::endl;
        exit(1);
    }

    Vars->GetVarValues(p);

    if (ReactorMod->IsInit()) {
        ReactorMod->compute_spec();
        tot_fitModel->Add(ReactorMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added ReactorMod, tot_fitModel integral = " << tot_fitModel->Integral() << std::endl;
        #endif
    }
    if (alphaNMod->IsInit()) {
        alphaNMod->compute_spec();
        tot_fitModel->Add(alphaNMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added alphaNMod, tot_fitModel integral = " << tot_fitModel->Integral() << std::endl;
        #endif
    }
    if (geoNuMod->IsInit()) {
        geoNuMod->compute_spec();
        tot_fitModel->Add(geoNuMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added geoNuMod, tot_fitModel integral = " << tot_fitModel->Integral() << std::endl;
        #endif
    }

    // Compute test statistic (minus, since minuit minizes, but we want max)
    return - ExtendedConstrainedLogLikelihood();
}

double Fitter::ExtendedConstrainedLogLikelihood() {
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: Running ExtendedConstrainedLogLikelihood()." << std::endl;
    #endif
    // Model PDFs have already been scaled by their respective norms, and added together
    // Assume data and models have the same binning (could generalise at some point)
    double logL = - tot_fitModel->Integral();
    for (unsigned int ibin = 1; ibin < numBins+1; ++ibin) {
        if (tot_fitModel->GetBinContent(ibin) > 0) logL += data->GetBinContent(ibin) * log(tot_fitModel->GetBinContent(ibin));
    }
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: +PDFs logL = " << logL << std::endl;
    #endif

    // Add in constraints
    FitVars* Vars = FitVars::GetInstance();
    for (unsigned int iVar = 0; iVar < Vars->GetNumVars(); ++ iVar) {
        // If variable is being held constant, assume it's equal to its prior most likely value
        if (!Vars->isConstant(iVar)) {
            // Add gaussian constraint
            logL -= 0.5 * (Vars->val(iVar) - Vars->prior(iVar))*(Vars->val(iVar) - Vars->prior(iVar)) / (Vars->err(iVar)*Vars->err(iVar));
        }
    }
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: +constraints logL = " << logL << std::endl;
    #endif

    return logL;
}

void Fitter::GetAllSpectra(std::vector<TH1D*>& hists) {

    TH1D* total_hist;

    Reactor* ReactorMod = Reactor::GetInstance();
    alphaN* alphaNMod = alphaN::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();

    if (!(ReactorMod->IsInit() || alphaNMod->IsInit() || geoNuMod->IsInit())) {
        std::cout << "Attempting to perform fit with no models!" << std::endl;
        exit(1);
    } else if (ReactorMod->IsInit()) {
        total_hist = (TH1D*)(ReactorMod->GetModelEsys()->Clone("Total_Model_Spectrum"));
    } else if (alphaNMod->IsInit()) {
        total_hist = (TH1D*)(alphaNMod->GetModelEsys()->Clone("Total_Model_Spectrum"));
    } else {
        total_hist = (TH1D*)(geoNuMod->GetModelEsys()->Clone("Total_Model_Spectrum"));
    }
    total_hist->SetTitle("Total Model Spectrum");

    if (ReactorMod->IsInit()) {
        ReactorMod->compute_spec();

        // Add total spectra
        hists.push_back((TH1D*)(ReactorMod->GetModelEsys()->Clone()));
        hists.at(hists.size()-1)->Add(ReactorMod->GetModelEsys());

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra
        ReactorMod->Spectra(hists);
    }

    if (alphaNMod->IsInit()) {
        alphaNMod->compute_spec();

        // Add total spectra
        hists.push_back((TH1D*)(alphaNMod->GetModelEsys()->Clone()));
        hists.at(hists.size()-1)->Add(alphaNMod->GetModelEsys());

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra
        alphaNMod->Spectra(hists);
    }

    if (geoNuMod->IsInit()) {
        geoNuMod->compute_spec();

        // Add total spectra
        hists.push_back((TH1D*)(geoNuMod->GetModelEsys()->Clone()));
        hists.at(hists.size()-1)->Add(geoNuMod->GetModelEsys());

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra
        geoNuMod->Spectra(hists);
    }

    hists.push_back(total_hist);
}

TH1D* Fitter::DataHist() {return data;}

void Fitter::SetData(TH1D* Data) {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::SetData]: Setting Data." << std::endl;
    #endif
    data = Data;
    numBins = data->GetXaxis()->GetNbins();
    tot_fitModel = (TH1D*)(data->Clone());
    tot_fitModel->SetName("tot_fit_model");
    tot_fitModel->SetTitle("Total Fit Spectrum");
    dataIsSet = true;
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