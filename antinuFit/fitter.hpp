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
        static double arglist[2];

        static Double_t* var_bestFits;
        static double* var_bestFit_errs;
        static unsigned int numVars;

        static std::vector<Esys*> Esysts;
        static unsigned int numEsysts;

        static std::vector<Model*> models;
        static unsigned int numModels;

        static TH1D data, tot_fitModel;
        static unsigned int numBins;

        // For minuit
        static double minfuncOut, edm, errdef;
        static int nvpar, nparx;

    protected:
        static FitVars Vars;
        static std::vector<TH1D*> hists;

    public:
        // Constructors
        Fitter();
        Fitter(TH1D& Data) {data = Data;};

        // Constructing functions function
        void AddVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool HoldConstant = false) {
            Vars.AddVar(Parname, Value, Verr, Vlow, Vhigh, HoldConstant);
        }
        void AddEsys(const double kB, const unsigned int linScale_idx, const unsigned int kBp_idx, const unsigned int sigPerRootE_idx) {
            Esysts.push_back(new Esys(numEsysts, kB, linScale_idx, kBp_idx, sigPerRootE_idx));
            ++numEsysts;
        };
        void AddEsys(const double kB, const std::string linScale_name, const std::string kBp_name, const std::string sigPerRootE_name) {
            Esysts.push_back(new Esys(numEsysts, kB, Vars.findIdx(linScale_name), Vars.findIdx(kBp_name), Vars.findIdx(sigPerRootE_name)));
            ++numEsysts;
        };

        void AddModel() {models.push_back(new Model(numModels)); ++numModels;};
        void AddReactorModel(const unsigned int Dm21_2_idx, const unsigned int Dm32_2_idx, const unsigned int s12_2_idx, const unsigned int s13_2_idx,
                             const std::vector<unsigned int>& norms_idx, const unsigned int totNorm_idx, const unsigned int Esys_idx,
                             const std::vector<TH1D*>& Reactor_hists, TH2D* E_conv_hist, const std::vector<std::string>& Reactor_names, RAT::DB* DB) {
            models.push_back(new Reactor(numModels, Dm21_2_idx, Dm32_2_idx, s12_2_idx, s13_2_idx, norms_idx, totNorm_idx, Esys_idx, Reactor_hists, E_conv_hist, Reactor_names, DB));
            ++numModels;
        };

        // Fitting functions
        static double fit_models();
        static Double_t fitFunc(Double_t* p);
        // (Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag*/)
        static void fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t) {fval = fitFunc(p);};
        static void SetFunc() {minuit->SetFCN(fitFunc_minuit_format);};

        static double GetFitParameter(const unsigned int parIdx) {return minuit->GetParameter(parIdx);};
        static double GetFitParError(const unsigned int parIdx) {return minuit->GetParError(parIdx);};
        static double GetFitCovarianceMatrixElement(const unsigned int i, const unsigned int j) {return minuit->GetCovarianceMatrixElement(i, j);};

        static double ExtendedConstrainedLogLikelihood(Double_t* p);

        // Get info functions
        void GetAllSpectra(std::vector<TH1D*>& hists);
        static FitVars& GetVars() {return Vars;};
        static std::vector<Model*>& GetModels() {return models;};
        static TH1D& GetData() {return data;};
        static std::vector<TH1D*>& GetHists() {return hists;};

        // Destructor
        ~Fitter() {};
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