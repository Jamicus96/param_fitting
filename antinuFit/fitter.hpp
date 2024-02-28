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
#include "E_systematics.hpp"
#include "model_Reactor.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"


class Fitter {
    private:
        static Fitter *FitterInstance_;

        static TVirtualFitter* minuit;
        double arglist[2];

        static TH1D *data, *tot_fitModel;
        static unsigned int numBins;

        // For minuit
        double minfuncOut, edm, errdef;
        int nvpar, nparx;

    
    protected:
        // Constructors/desctructors
        Fitter() {}
        ~Fitter() {}

        bool dataIsSet;
        bool isInit;
        // FitVars* Vars = FitVars::GetInstance();
        Esys* Esysts = Esys::GetInstance();
        static unsigned int numMods;

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
        double fit_models() {
            if (!dataIsSet || numMods == 0 || isInit == 0) {
                std::cout << "Attempting to perform fit with no models or no data, or no initialisation!" << std::endl;
                exit(1);
            }

            // minimize
            #ifdef DEBUG
                std::cout << "[Fitter::fit_models]: Running MIGRAD!" << std::endl;
                std::cout << "[Fitter::fit_models]: s_12^2 = " << variables.at(variables.size() - 2)->val() << std::endl;
                std::cout << "[Fitter::fit_models]: Dm21^2 = " << variables.at(variables.size() - 4)->val() << std::endl;
            #endif

            minuit->ExecuteCommand("MIGRAD", arglist, 2);

            //get result
            minuit->GetStats(minfuncOut, edm, errdef, nvpar, nparx);
            std::cout << "[Fitter::fit_models]: minfuncOut = " << minfuncOut << std::endl;
            return minfuncOut;
        }

        static Double_t fitFunc(Double_t* p);

        // (Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag*/)
        static void fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t);
        void SetFunc() {minuit->SetFCN(fitFunc_minuit_format);}

        double GetFitParameter(const unsigned int parIdx) {return minuit->GetParameter(parIdx);}
        double GetFitParError(const unsigned int parIdx) {return minuit->GetParError(parIdx);}
        double GetFitCovarianceMatrixElement(const unsigned int i, const unsigned int j) {return minuit->GetCovarianceMatrixElement(i, j);}

        static double ExtendedConstrainedLogLikelihood();

        // Get info functions
        void GetAllSpectra(std::vector<TH1D*>& hists) {

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


        TH1D* DataHist() {return data;}

        void SetData(TH1D* Data) {
            data = Data;
            numBins = data->GetXaxis()->GetNbins();
            tot_fitModel = (TH1D*)(data->Clone());
            tot_fitModel->SetName("tot_fit_model");
            tot_fitModel->SetTitle("Total Fit Spectrum");
            dataIsSet = true;
        }

        void Initialise() {
            //The default minimizer is Minuit, you can also try Minuit2
            // TVirtualFitter::SetDefaultFitter("Minuit2");
            TVirtualFitter::SetDefaultFitter("Minuit");
            FitVars* Vars = FitVars::GetInstance();
            minuit = TVirtualFitter::Fitter(0, Vars->GetNumVars());
            
            arglist[0] = 0;
            // set print level
            minuit->ExecuteCommand("SET PRINT", arglist, 2);

            // minimiser settings
            arglist[0] = 5000; // number of function calls
            arglist[1] = 0.01; // tolerance

            // (void (*fcn)(Int_t&, Double_t*, Double_t& f, Double_t*, Int_t))
            SetFunc();

            isInit = true;
        }
};

// Static methods should be defined outside the class.
Fitter* Fitter::FitterInstance_{nullptr};

TVirtualFitter* Fitter::minuit;
TH1D *Fitter::data, *Fitter::tot_fitModel;
unsigned int Fitter::numBins;
unsigned int Fitter::numMods;

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
Fitter* Fitter::GetInstance() {
    if (FitterInstance_ == nullptr) {
        FitterInstance_ = new Fitter();
    }
    return FitterInstance_;
}

void Fitter::fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t) {fval = fitFunc(p);}

Double_t Fitter::fitFunc(Double_t* p) {
    // Reset total model spectrum to zero
    tot_fitModel->Reset("ICES");

    // Compute the total spectrum from all the models, using parameters p
    FitVars* Vars = FitVars::GetInstance();
    Vars->GetVarValues(p);

    Reactor* ReactorMod = Reactor::GetInstance();
    if (ReactorMod->IsInit()) {
        ReactorMod->compute_spec();
        tot_fitModel->Add(ReactorMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added ReactorMod, tot_fitModel integral = " << tot_fitModel->Integral() << std::endl;
        #endif
    }

    alphaN* alphaNMod = alphaN::GetInstance();
    if (alphaNMod->IsInit()) {
        alphaNMod->compute_spec();
        tot_fitModel->Add(alphaNMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added alphaNMod, tot_fitModel integral = " << tot_fitModel->Integral() << std::endl;
        #endif
    }

    geoNu* geoNuMod = geoNu::GetInstance();
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



//end header guard
#endif


/* ~~~~~ Saving other potentially useful commands here ~~~~~ */

// unsigned int NtotParams = minuit->GetNumberTotalParameters();
// unsigned int NfreeParams = minuit->GetNumberFreeParameters();

// double chi2, edm, errdef;
// int nvpar, nparx;
// minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

// void TFitter::GetConfidenceIntervals(Int_t n, Int_t ndim, const Double_t *x, Double_t *ci, Double_t cl)