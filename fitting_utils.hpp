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


class reactorINFO {
    private:
        // Initial oscillation parameters that only depend on
        // Dm_21^2, Dm_32^2, s_12^2, s_13^2 and electron density
        double fDmSqr21;
        double fDmSqr32;
        double fSSqrTheta12;
        double fSSqrTheta13;

        double H_ee_vac;
        double a0_vac;
        double a1_vac;
        double Y_ee_vac;
        double alpha;

        // Computed oscillation paramters depending only on E [MeV], but not on L [km]
        std::vector<double> eigen;  // Antinu propagation Hamiltonian eigenvalues
        std::vector<double> X_mat;  // Values from matrix governing oscillation

        // Raw reactor information
        std::vector<double> baselines;  // Baselines [km]
        std::vector<TH1D*> reactor_hists;  // Un-oscillated reactor IBD spectra
        unsigned int num_reactors;
        unsigned int hists_Nbins;

        // Oscillated reactor hist
        std::vector<TH1D*> osc_reactor_hist;

    public:
        // Constructors
        reactorINFO(std::vector<TH1D*>& Reactor_hists, std::vector<double>& Baselines);
        reactorINFO(std::vector<TH1D*>& Reactor_hists, std::vector<double>& Baselines, const double DmSqr21, const double DmSqr32, const double SSqrTheta12, const double SSqrTheta13);
        reactorINFO(std::vector<TH1D*>& Reactor_hists, std::vector<double>& Baselines, const double DmSqr32, const double SSqrTheta13);

        // Destructor
        ~reactorINFO() {};

        // Member function
        double& Dm21_2() {return fDmSqr21;};
        double& Dm32_2() {return fDmSqr32;};
        double& s12_2() {return fSSqrTheta12;};
        double& s13_2() {return fSSqrTheta13;};
        const std::vector<TH1D*>& Get_osc_reactor_specs() {return osc_reactor_hist;};

        void compute_oscillation_constants();
        void re_compute_consts(const double E);
        double survival_prob(const double E, const double L);

        void compute_osc_reactor_spec();
};

//end header guard
#endif