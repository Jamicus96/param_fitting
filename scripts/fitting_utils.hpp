// header guard:
#ifndef fitting_utils
#define fitting_utils

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <RAT/DB.hh>

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
        double N_IBD;  // number of unoscilalted reactor events
        double IBD_err;  // fractional error in N_IBD
        std::vector<double> baselines;  // Baselines [km]
        std::vector<TH1D*> reactor_hists;  // Un-oscillated reactor IBD spectra
        unsigned int num_reactors;
        unsigned int hists_Nbins;
        TH2D E_conv;  // Conversion from E_e to E_nu (normalised for each E_e bin)

        // Oscillated reactor hists
        // Names of reactors whose norms can float independently, last one is always 'WORLD':
        std::vector<std::string> reactor_names = {"BRUCE", "DARLINGTON", "PICKERING", "WORLD"};  
        std::vector<unsigned int> reactor_idx;  // Maps each reactor_hists to the idx it should be assigned to in reactor_names
        std::vector<TH1D*> osc_hists;
        double tot_hist_int;
        // std::map<std::string, unsigned int> osc_hists_map = {{"BRUCE", 0}, {"DARLINGTON", 1}, {"PICKERING", 2}, {"WORLD", 3}};
        std::vector<double> norms;

    public:
        // Constructors
        reactorINFO(std::vector<TH1D*>& Reactor_hists, const double N_IBDs, const double IBD_errs, TH2D& E_conv_hist);
        reactorINFO(std::vector<TH1D*>& Reactor_hists, const double N_IBDs, const double IBD_errs, TH2D& E_conv_hist, const double DmSqr21, const double DmSqr32, const double SSqrTheta12, const double SSqrTheta13);
        reactorINFO(std::vector<TH1D*>& Reactor_hists, const double N_IBDs, const double IBD_errs, TH2D& E_conv_hist, const double DmSqr32, const double SSqrTheta13);

        // Member function
        double& Dm21_2() {return fDmSqr21;};
        const double& Dm21_2() const {return fDmSqr21;};
        double& Dm32_2() {return fDmSqr32;};
        const double& Dm32_2() const {return fDmSqr32;};
        double& s12_2() {return fSSqrTheta12;};
        const double& s12_2() const {return fSSqrTheta12;};
        double& s13_2() {return fSSqrTheta13;};
        const double& s13_2() const {return fSSqrTheta13;};
        double& N_IBDs() {return N_IBD;};
        const double& N_IBDs() const {return N_IBD;};
        double& IBDs_err() {return IBD_err;};
        const double& IBDs_err() const {return IBD_err;};

        void compute_baselines(RAT::DB* db);

        std::vector<TH1D*>& Get_osc_reactor_specs() {return osc_hists;};
        const std::vector<TH1D*>& Get_osc_reactor_specs() const {return osc_hists;};
        std::vector<double>& Get_osc_reactor_norms() {return norms;};
        const std::vector<double>& Get_osc_reactor_norms() const {return norms;};

        void compute_oscillation_constants();
        void re_compute_consts(const double E);
        double survival_prob(const double E, const double L);

        void compute_osc_reactor_spec();

        // Geo-nu survival probability: averaged over baseline -> No E-depence, only depends on mixing angles.
        double geoNu_survival_prob() {return fSSqrTheta13*fSSqrTheta13 + (1. - fSSqrTheta13)*(1. - fSSqrTheta13) * (1. - 2. * fSSqrTheta12 * (1. - fSSqrTheta12));};

        // Destructor
        ~reactorINFO() {
            for (auto p : reactor_hists) {delete p;}
            reactor_hists.clear();
            for (auto p : osc_hists) {delete p;}
            osc_hists.clear();
        };
};

std::vector<std::string> SplitString(std::string str);
TVector3 LLAtoECEF(double longitude, double latitude, double altitude);
double GetReactorDistanceLLA(const double &longitude, const double&latitude, const double &altitude);

//end header guard
#endif