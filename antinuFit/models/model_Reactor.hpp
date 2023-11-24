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
        unsigned int numVars, iDm_21_2, iDm_32_2, iS_12_2, iS_13_2;
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
        Reactor(const Reactor& mod) : Model(mod) {
            numVars = mod.numVars; iDm_21_2 = mod.iDm_21_2; iDm_32_2 = mod.iDm_32_2; iS_12_2 = mod.iS_12_2; iS_13_2 = mod.iS_13_2;
            iNorms = mod.iNorms; H_ee_vac = mod.H_ee_vac; a0_vac = mod.a0_vac; a1_vac = mod.a1_vac; Y_ee_vac = mod.Y_ee_vac;
            eigen = mod.eigen; X_mat = mod.X_mat; baselines = mod.baselines; reactor_hists = mod.reactor_hists; num_reactors = mod.num_reactors;
            hists_Nbins = mod.hists_Nbins; E_conv = mod.E_conv; reactor_names = mod.reactor_names; reactor_idx = mod.reactor_idx;
            osc_hists = mod.osc_hists; unosc_hist_ints = mod.unosc_hist_ints; computed_osc_specs = mod.computed_osc_specs; db = mod.db;
        };
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

        void compute_spec();

        std::vector<std::string> GetReactorNames() {return reactor_names;};
        std::vector<TH1D*>& GetOscReactorHists();
    
        // Destructor
        ~Reactor() {
            for (auto p : reactor_hists) {delete p;}
            reactor_hists.clear();
            for (auto p : osc_hists) {delete p;}
            osc_hists.clear();
            delete db;
            for (auto p : Vars) {delete p;} Vars.clear();
        };
};

/* ~~~~~~~~~~~~~~~~ OTHER (NON-MEMBER) FUNCTIONS ~~~~~~~~~~~~~~~~ */


std::vector<std::string> SplitString(std::string str);

TVector3 LLAtoECEF(const double& longitude, const double& latitude, const double& altitude);

double GetReactorDistanceLLA(const double& longitude, const double& latitude, const double &altitude);



//end header guard
#endif