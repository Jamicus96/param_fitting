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
#include "model.hpp"
#include "model_Reactor.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"

class Fitter {
    private:
        static TVirtualFitter* minuit;
        static double arglist[2];

        static TH1D* data, tot_fitModel;
        static unsigned int numBins;

        // For minuit
        static double minfuncOut, edm, errdef;
        static int nvpar, nparx;

    protected:
        static FitVars Vars;
        static Esys Esysts;
        static std::vector<Model*> Mods;
        static unsigned int numMods;

    public:
        // Constructors
        Fitter();
        Fitter(TH1D* Data) {data = Data;};

        // Adding models
        void AddReactorMod(const std::string Dm21_2_name, const std::string Dm32_2_name, const std::string s12_2_name, const std::string s13_2_name,
                           const std::vector<std::string>& norms_names, const std::string totNorm_name, const std::string Esys_name,
                           const std::vector<TH1D*>& Reactor_hists, TH2D* E_conv_hist, const std::vector<std::string>& Reactor_names, RAT::DB* DB) {
            Mods.push_back(new Reactor(Dm21_2_name, Dm32_2_name, s12_2_name, s13_2_name, norms_names, totNorm_name, Esys_name, Reactor_hists, E_conv_hist, Reactor_names, DB));
            ++numMods;
        };
        void AddAlphaNMod(const std::string NormGS_name, const std::string NormES_name, const std::string E_syst_name,
                          const std::string E_syst_proton_name, TH1D* Hist_ProtontR, TH1D* Hist_C12Scatter, TH1D* Hist_O16Deex) {
            Mods.push_back(new alphaN(NormGS_name, NormES_name, E_syst_name, E_syst_proton_name, Hist_ProtontR, Hist_C12Scatter, Hist_O16Deex));
            ++numMods;
        };
        void AddGeoNuMod(const std::string NormTh_name, const std::string NormU_name, const std::string vS_12_2_name, const std::string vS_13_2_name, const std::string E_syst_name, TH1D* Hist_Th, TH1D* Hist_U) {
            Mods.push_back(new geoNu(NormTh_name, NormU_name, vS_12_2_name, vS_13_2_name, E_syst_name, Hist_Th, Hist_U));
            ++numMods;
        }

        // Fitting functions
        static double fit_models();
        static Double_t fitFunc(Double_t* p);
        // (Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag*/)
        static void fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t) {fval = fitFunc(p);};
        static void SetFunc() {minuit->SetFCN(fitFunc_minuit_format);};

        static double GetFitParameter(const unsigned int parIdx) {return minuit->GetParameter(parIdx);};
        static double GetFitParError(const unsigned int parIdx) {return minuit->GetParError(parIdx);};
        static double GetFitCovarianceMatrixElement(const unsigned int i, const unsigned int j) {return minuit->GetCovarianceMatrixElement(i, j);};

        static double ExtendedConstrainedLogLikelihood();

        // Get info functions
        void GetAllSpectra(std::vector<TH1D*>& hists);
        static FitVars& GetVars() {return Vars;};
        static Esys& GetEsysts() {return Esysts;};
        static std::vector<Model*>& GetModels() {return Mods;};
        static TH1D* DataHist() {return data;};
        static void SetData(TH1D* Data) {data = Data;};

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