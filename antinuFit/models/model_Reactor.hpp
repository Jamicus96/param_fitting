// header guard:
#ifndef reactor_model
#define reactor_model

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <RAT/DB.hh>
#include "fitter.hpp"
#include "model.hpp"
#include "fitVars.hpp"
#include "E_systematics.hpp"

class Reactor : public Model {
    private:
        std::string ModName = "Reactor";
        // Indices pointing to variables
        unsigned int iDm_21_2, iDm_32_2, iS_12_2, iS_13_2, iTotNorm;
        std::vector<unsigned int> iNorms;
        unsigned int iEsys;

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
        Reactor(const unsigned int Dm21_2_idx, const unsigned int Dm32_2_idx, const unsigned int s12_2_idx, const unsigned int s13_2_idx,
                const std::vector<unsigned int>& norms_idx, const unsigned int totNorm_idx, const unsigned int Esys_idx,
                const std::vector<TH1D*>& Reactor_hists, TH2D* E_conv_hist, const std::vector<std::string>& Reactor_names, RAT::DB* DB);
        
        Reactor(const std::string Dm21_2_name, const std::string Dm32_2_name, const std::string s12_2_name, const std::string s13_2_name,
                const std::vector<std::string>& norms_names, const std::string totNorm_name, const std::string Esys_name,
                const std::vector<TH1D*>& Reactor_hists, TH2D* E_conv_hist, const std::vector<std::string>& Reactor_names, RAT::DB* DB) {
            std::vector<unsigned int> norms_idx;
            for (unsigned int i = 0; i < norms_names.size(); ++i) norms_idx.push_back(Vars.findIdx(norms_names.at(i)));
            Reactor(Vars.findIdx(Dm21_2_name), Vars.findIdx(Dm32_2_name), Vars.findIdx(s12_2_name), Vars.findIdx(s13_2_name), norms_idx,
                    Vars.findIdx(totNorm_name), Esysts.findIdx(Esys_name), Reactor_hists, E_conv_hist, Reactor_names, DB);
        };

        // Member function
        void compute_unosc_integrals();
        void compute_osc_specs();
        void compute_oscillation_constants();
        void re_compute_consts(const double& E);
        double survival_prob(const double& E, const double& L);
        void compute_baselines();
        void hold_osc_params_const(const bool isTrue);

        void compute_spec();

        std::vector<std::string>& GetReactorNames() {return reactor_names;};
        void Spectra(std::vector<TH1D*>& hists);
    
        // Destructor
        ~Reactor();
};

std::vector<std::string> SplitString(std::string str);
TVector3 LLAtoECEF(const double& longitude, const double& latitude, const double& altitude);
double GetReactorDistanceLLA(const double& longitude, const double& latitude, const double &altitude);

//end header guard
#endif