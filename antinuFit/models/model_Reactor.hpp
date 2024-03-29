// header guard:
#ifndef reactor_model
#define reactor_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <RAT/DB.hh>
#include "model.hpp"
#include "fitVar.hpp"
#include <TMath.h>


class Reactor: public Model {
    private:
        // Indices pointing to variables
        unsigned int iDm_21_2, iDm_32_2, iS_12_2, iS_13_2;
        std::vector<unsigned int> iNorms;

        // Initial oscillation parameters that only depend on
        // Dm_21^2, Dm_32^2, s_12^2, s_13^2 and electron density
        double H_ee_vac;
        double a0_vac;
        double a1_vac;
        double Y_ee_vac;

        // Define electron density Ne of the crust, based on 2.7g/cm3 mass density, and <N/A> = 0.5
        static constexpr double alpha = - 2.535e-31 * 8.13e23;  // conversion factor in eV2/MeV * Ne = 8.13e23

        // Computed oscillation paramters depending only on E [MeV], but not on L [km]
        std::vector<double> eigen;  // Antinu propagation Hamiltonian eigenvalues
        std::vector<double> X_mat;  // Values from matrix governing oscillation

        std::vector<double> baselines;  // Baselines [km]
        std::vector<TH1D*> reactor_hists;  // Un-oscillated reactor IBD spectra
        unsigned int num_reactors;
        unsigned int hists_Nbins;
        TH2D* E_conv;  // Conversion from E_e to E_nu (normalised for each E_e bin)

        // Oscillated reactor hists
        // Names of reactors whose norms can float independently, last one is always 'WORLD':
        std::vector<std::string> reactor_names; // example: {"BRUCE", "DARLINGTON", "PICKERING", "WORLD"}  
        std::vector<unsigned int> reactor_idx;  // Maps each reactor_hists to the idx it should be assigned to in reactor_names
        std::vector<TH1D*> osc_hists;
        std::vector<double> unosc_hist_ints;  // Integral un un-oscillated histograms (computed once)
        bool computed_osc_specs = false;

        RAT::DB* db;

    public:
        // Constructors
        Reactor(const Reactor& mod);
        void operator = (const Reactor& mod);
        Reactor(FitVar* vDm_21_2, FitVar* vDm_32_2, FitVar* vS_12_2, FitVar* vS_13_2, const std::vector<TH1D*>& Reactor_hists,
                TH2D* E_conv_hist, std::vector<std::string>& Reactor_names, const std::vector<FitVar*>& Norms, RAT::DB* DB);

        // Member function
        void compute_unosc_integrals();
        void compute_osc_specs();
        void compute_oscillation_constants();
        void re_compute_consts(const double& E);
        double survival_prob(const double& E, const double& L);
        void compute_baselines();
        void hold_osc_params_const(bool isTrue);

        void compute_spec(Double_t* p);

        std::vector<std::string>& GetReactorNames();
        void GetOscReactorHists(std::vector<TH1D*>& rescaled_osc_hists);
    
        // Destructor
        ~Reactor();
};

std::vector<std::string> SplitString(std::string str);
TVector3 LLAtoECEF(const double& longitude, const double& latitude, const double& altitude);
double GetReactorDistanceLLA(const double& longitude, const double& latitude, const double &altitude);

//end header guard
#endif