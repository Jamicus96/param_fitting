// header guard:
#ifndef fitting_utils
#define fitting_utils

// include
#include <iostream>
#include <TH1.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
// #include "RooDataSet.h"
#include <RooClassFactory.h>
#include <RooTFnBinding.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooConstVar.h>

using namespace RooFit;


// double fDmSqr21;  // \Delta m_{21}^2 [MeV^2]
// double fDmSqr32;  // \Delta m_{32}^2 [MeV^2]
// double fSSqrTheta12;  // s_{12}^2
// double fSSqrTheta13;  // s_{13}^2


class PDFspec {
    private:
        // Initial oscillation parameters that only depend on
        // \Delta m_{21}^2, \Delta m_{32}^2, s_{12}^2, s_{13}^2 and electron density
        double H_ee_vac;
        double a0_vac;
        double a1_vac;
        double Y_ee_vac;
        double alpha;

        // Computed oscillation paramters depending only on E [MeV], but not on L [km]
        std::vector<double> eigen;  // Antinu propagation Hamiltonian eigenvalues
        std::vector<double> X_mat;  // Values from matrix governing oscillation

        // Raw reactor information
        std::vector<double>* baselines;  // Baselines [km]
        std::vector<TH1D*>* reactor_hists;  // Un-oscillated reactor IBD spectra
        unsigned int num_reactors;
        unsigned int hists_Nbins;

        // PDFs, after oscillation (RooFit)
        TH1D* tot_reactor_hist;
        RooDataHist* tempData_react;
        RooHistPdf* reactor_PDF;

        RooDataHist* tempData_alpha;
        RooHistPdf* alphaN_PDF;

        // PDF model and fit result
        RooAddPdf* model;
        RooFitResult *result;

        // Other RooFit paramters
        RooRealVar e;  // Energy (MeV)
        RooRealVar reactor_frac;
        RooRealVar alphaN_frac;

    public:
        // Constructors
        PDFspec(TH1D* AlphaN_hist, std::vector<TH1D*>* Reactor_hists, std::vector<double>* Baselines);

        // Member function
        void compute_oscillation_constants(const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13);
        void re_compute_consts(const double E);
        double survival_prob(const double E, const double L);

        void compute_tot_reactor_spec(const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13);
        double ML_fit(RooDataHist& dataHist);

};

//end header guard
#endif