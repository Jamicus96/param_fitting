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
#include "fitVars.hpp"
#include "E_systematics.hpp"


std::vector<std::string> SplitString(std::string str);
TVector3 LLAtoECEF(const double& longitude, const double& latitude, const double& altitude);
double GetReactorDistanceLLA(const double& longitude, const double& latitude, const double &altitude);

class Reactor {
    private:
        static Reactor *ReactorInstance_;

        TH1D *model_noEsys, *model_Esys;
        double model_noEsys_integral;
        bool isInit;

        // Indices pointing to variables
        unsigned int iDm_21_2, iDm_32_2, iS_12_2, iS_13_2, iNorm;
        unsigned int iEsys;
        bool computed_osc_specs;

        // Initial oscillation parameters that only depend on
        // Dm_21^2, Dm_32^2, s_12^2, s_13^2 and electron density
        double H_ee_vac;
        double a0_vac;
        double a1_vac;
        double Y_ee_vac;

        // Define electron density Ne of the crust, based on 2.7g/cm3 mass density, and <N/A> = 0.5
        static constexpr double alpha = - 2.535e-31 * 8.13e23;  // -2*sqrt(2)*G_F*Ne [eV2/MeV], for Ne = 8.13e23 [cm^-3]

        // Computed oscillation paramters depending only on E [MeV], but not on L [km]
        std::vector<double> eigen;  // Antinu propagation Hamiltonian eigenvalues
        std::vector<double> X_mat;  // Values from matrix governing oscillation

        unsigned int hists_Nbins;
        TH2D* E_conv;  // Conversion from E_e to E_nu (normalised for each E_e bin)
        TH1D *PWR_promptE_hist, *PWR_Enu_hist, *PHWR_Enu_hist;

        double fPWR_promptE_frac;
        std::vector<double> fPWR_Enu_fracs, fPHWR_Enu_fracs, fPWR_Enu_baselines, fPHWR_Enu_baselines;
        unsigned int num_PWR, num_PHWR;

        double av_survival_prob;
        unsigned int iMinBin, iMaxBin;

        RAT::DB* db;
    
    protected:
        // Constructors/desctructors
        Reactor() : isInit(false), computed_osc_specs(false) {}
        ~Reactor() {}

    public:
        // Should not be cloneable.
        Reactor(Reactor &other) = delete;

        // Should not be assignable.
        void operator=(const Reactor&) = delete;
        /**
         * This is the method that controls the access to the Reactor
         * instance. On the first run, it creates a Reactor object and places it
         * into the field. On subsequent runs, it returns the client existing
         * object stored in the field.
         */
        static Reactor *GetInstance();

        // Initialisers
        void InitReactor(const unsigned int Dm21_2_idx, const unsigned int Dm32_2_idx, const unsigned int s12_2_idx, const unsigned int s13_2_idx,
                const unsigned int Norm_idx, const unsigned int Esys_idx,
                const std::vector<TH1D*> Reactor_hists, TH2D* E_conv_hist, const Double_t PWR_promptE_frac, const std::vector<Double_t>& PWR_Enu_fracs,
                const std::vector<Double_t>& PHWR_Enu_fracs, const std::vector<Double_t>& PWR_Enu_baselines, const std::vector<Double_t>& PHWR_Enu_baselines, RAT::DB* DB) {
                        
            if (PWR_Enu_fracs.size() != PWR_Enu_baselines.size() || PHWR_Enu_fracs.size() != PHWR_Enu_baselines.size()) {
                std::cout << "[Reactor::InitReactor] ERROR: flux fractions and baselines should have the same size!" << std::endl;
                exit(1);
            }

            #ifdef antinuDEBUG
                FitVars* Vars = FitVars::GetInstance();
                std::cout << "[Reactor::InitReactor]: Vars->GetNumVars() = " << Vars->GetNumVars() << std::endl;
            #endif

            // Save variables
            iDm_21_2 = Dm21_2_idx;
            iDm_32_2 = Dm32_2_idx;
            iS_12_2 = s12_2_idx;
            iS_13_2 = s13_2_idx;
            iEsys = Esys_idx;
            iNorm = Norm_idx;
            db = DB;

            #ifdef antinuDEBUG
                std::cout << "[Reactor::InitReactor]: iDm_21_2 = " << iDm_21_2 << ", iDm_32_2 = " << iDm_32_2 << ", iS_12_2 = " << iS_12_2 << ", iS_13_2 = " << iS_13_2 << ", iEsys = " << iEsys << ", iNorm = " << iNorm << std::endl;
            #endif

            // Save hist variables
            E_conv = (TH2D*)(E_conv_hist->Clone("Reactor::E_conv"));
            E_conv->Reset("ICES"); E_conv->Add(E_conv_hist);

            bool bPWR_promptE = false, bPWR_Enu = false, bPHWR_promptE = false, bPHWR_Enu = false;
            for (unsigned int i = 0; i < Reactor_hists.size(); ++i) {
                if (Reactor_hists.at(i)->GetName() == "PWR_promptE") {
                    PWR_promptE_hist = (TH1D*)(Reactor_hists.at(i)->Clone("Reactor::PWR_promptE"));
                    PWR_promptE_hist->Reset("ICES"); PWR_promptE_hist->Add(Reactor_hists.at(i));
                    bPWR_promptE = true;
                } else if (Reactor_hists.at(i)->GetName() == "PWR_Enu") {
                    PWR_Enu_hist = (TH1D*)(Reactor_hists.at(i)->Clone("Reactor::PWR_Enu"));
                    PWR_Enu_hist->Reset("ICES"); PWR_Enu_hist->Add(Reactor_hists.at(i));
                    bPWR_Enu = true;
                } else if (Reactor_hists.at(i)->GetName() == "PHWR_Enu") {
                    PHWR_Enu_hist = (TH1D*)(Reactor_hists.at(i)->Clone("Reactor::PHWR_Enu"));
                    PHWR_Enu_hist->Reset("ICES"); PHWR_Enu_hist->Add(Reactor_hists.at(i));
                    bPHWR_Enu = true;
                }
            }
            if (!(bPWR_promptE && bPWR_Enu && bPHWR_promptE && bPHWR_Enu)) {
                std::cout << "[Reactor::InitReactor]: ERROR: PDF missing!" << std::endl;
                exit(1);
            }

            hists_Nbins = PWR_promptE_hist->GetXaxis()->GetNbins();

            SetBinLims(1, PWR_promptE_hist->GetXaxis()->GetNbins());

            fPWR_promptE_frac = PWR_promptE_frac;
            num_PWR = PWR_Enu_fracs.size();
            num_PHWR = PHWR_Enu_fracs.size();
            for (unsigned int i = 0; i < num_PWR; ++i) {
                fPWR_Enu_fracs.push_back(PWR_Enu_fracs.at(i));
                fPWR_Enu_baselines.push_back(PWR_Enu_baselines.at(i));
            }
            for (unsigned int i = 0; i < num_PHWR; ++i) {
                fPHWR_Enu_fracs.push_back(PHWR_Enu_fracs.at(i));
                fPHWR_Enu_baselines.push_back(PHWR_Enu_baselines.at(i));
            }

            // Set up empty vectors
            eigen.resize(3);
            X_mat.resize(3);

            model_noEsys = (TH1D*)(PWR_promptE_hist->Clone("Reactor::model_noEsys"));
            model_Esys = (TH1D*)(PWR_promptE_hist->Clone("Reactor::model_Esys"));

            #ifdef antinuDEBUG
                std::cout << "[Reactor::InitReactor]: iDm_21_2 = " << iDm_21_2 << ", iDm_32_2 = " << iDm_32_2 << ", iS_12_2 = " << iS_12_2 << ", iS_13_2 = " << iS_13_2 << ", iEsys = " << iEsys << ", iNorm = " << iNorm << std::endl;
            #endif

            isInit = true;
        }
        
        void InitReactor(const std::string Dm21_2_name, const std::string Dm32_2_name, const std::string s12_2_name, const std::string s13_2_name,
                const std::string Norm_name, const std::string Esys_name, const std::vector<TH1D*>& Reactor_hists, TH2D* E_conv_hist,
                const Double_t PWR_promptE_frac, const std::vector<Double_t>& PWR_Enu_fracs, const std::vector<Double_t>& PHWR_Enu_fracs,
                const std::vector<Double_t>& PWR_Enu_baselines, const std::vector<Double_t>& PHWR_Enu_baselines, RAT::DB* DB) {

            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();

            #ifdef antinuDEBUG
                std::cout << "[Reactor::InitReactor]: Vars->GetNumVars() = " << Vars->GetNumVars() << std::endl;
            #endif

            InitReactor(Vars->findIdx(Dm21_2_name), Vars->findIdx(Dm32_2_name), Vars->findIdx(s12_2_name), Vars->findIdx(s13_2_name),
                    Vars->findIdx(Norm_name), Esysts->findIdx(Esys_name), Reactor_hists, E_conv_hist, PWR_promptE_frac, PWR_Enu_fracs,
                    PHWR_Enu_fracs, PWR_Enu_baselines, PHWR_Enu_baselines, DB);
        };

        void compute_spec() {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            model_Esys->Reset("ICES");  // empty it before re-computing it

            // If the oscillation constants are being held constant, and the oscillated spectra have already been computed, can skip this expensive step!
            if (!(Vars->IsConstant(iDm_21_2) && Vars->IsConstant(iDm_32_2) && Vars->IsConstant(iS_12_2) && Vars->IsConstant(iS_13_2) && computed_osc_specs)) {
                #ifdef SUPER_DEBUG
                    std::cout << "[Reactor::compute_spec]: computing oscillation specs" << std::endl;
                #endif
                compute_osc_specs();
            }

            model_noEsys->Scale(Vars->val(iNorm) / model_noEsys_integral);
            
            // Apply energy systematics
            Esysts->apply_systematics(iEsys, model_noEsys, model_Esys);

            #ifdef antinuDEBUG
                std::cout << "[Reactor::compute_spec]: Vars->val(iNorm) = " << Vars->val(iNorm) << std::endl;
                std::cout << "[Reactor::compute_spec]: model_noEsys->Integral(iMinBin, iMaxBin) = " << model_noEsys->Integral(iMinBin, iMaxBin) << std::endl;
                std::cout << "[Reactor::compute_spec]: model_Esys->Integral(iMinBin, iMaxBin) = " << model_Esys->Integral(iMinBin, iMaxBin) << std::endl;
            #endif
        }

        void compute_osc_specs() {
            // Compute oscillation constants
            compute_oscillation_constants();

            model_noEsys->Reset("ICES");  // empty it before re-computing it

            // Assume all the histograms have the same E binning
            double Enu, tot_weight, weight;
            std::string origin_reactor;
            #ifdef antinuDEBUG
                std::cout << "[Reactor::compute_osc_specs]: hists_Nbins = " << hists_Nbins << std::endl;
                std::cout << "[Reactor::compute_osc_specs]: E_conv->GetXaxis()->GetNbins() = " << E_conv->GetXaxis()->GetNbins() << std::endl;
                std::cout << "[Reactor::compute_osc_specs]: E_conv->GetYaxis()->GetNbins() = " << E_conv->GetYaxis()->GetNbins() << std::endl;
            #endif

            // Add constribution from distant reactors (average scaling from oscillation)
            model_noEsys->Add(PWR_promptE_hist, av_survival_prob * fPWR_promptE_frac);

            // Loop over event energy bins (prompt energy)
            for (unsigned int iBin = 1; iBin <= hists_Nbins; ++iBin) {
                // Loop over corresponding E_nu bins to apply oscillation and energy conversion
                tot_weight = 0;
                for (unsigned int iEnu = 1; iEnu <= E_conv->GetXaxis()->GetNbins(); ++iEnu) {
                    weight = E_conv->GetBinContent(iEnu, iBin);

                    // Most of the histogram is empty, so only bother if > 0
                    if (weight > 0.0) {
                        // Compute associated antinu energy, and re-compute oscillation parameters
                        Enu = E_conv->GetXaxis()->GetBinCenter(iEnu);
                        re_compute_consts(Enu);

                        // Loop over all PWRs
                        double PWR_weight = 0;
                        for (unsigned int iPWR = 0; iPWR < num_PWR; ++iPWR) {
                            // Add value of bin from reactor core to appropriate hist, scaled by survival probability
                            PWR_weight += survival_prob(Enu, fPWR_Enu_baselines.at(iPWR)) * fPWR_Enu_fracs.at(iPWR);
                        }
                        PWR_weight *= PWR_Enu_hist->GetBinContent(iEnu);

                        // Loop over all PHWRs
                        double PHWR_weight = 0;
                        for (unsigned int iPHWR = 0; iPHWR < num_PHWR; ++iPHWR) {
                            // Add value of bin from reactor core to appropriate hist, scaled by survival probability
                            PHWR_weight += survival_prob(Enu, fPHWR_Enu_baselines.at(iPHWR)) * fPHWR_Enu_fracs.at(iPHWR);
                        }
                        PHWR_weight *= PHWR_Enu_hist->GetBinContent(iEnu);

                        tot_weight += weight * (PWR_weight + PHWR_weight);
                    }
                }

                // Add total oscillated weight to the bin, multiplied by the total normalisation
                model_noEsys->AddBinContent(iBin, tot_weight);
                
                #ifdef SUPER_DEBUG
                    double prompt_E = model_noEsys->GetXaxis()->GetBinCenter(iBin);
                    std::cout << "[Reactor::compute_osc_specs]: iBin = " << iBin << ", prompt_E = " << prompt_E << ", tot_weight = " << tot_weight << std::endl;
                #endif
            }

            model_noEsys_integral = model_noEsys->Integral(iMinBin, iMaxBin);

            #ifdef antinuDEBUG
                std::cout << "[Reactor::compute_osc_specs]: model_noEsys->Integral(iMinBin, iMaxBin) = " << model_noEsys->Integral(iMinBin, iMaxBin) << std::endl;
            #endif
        }

        void compute_oscillation_constants() {
            FitVars* Vars = FitVars::GetInstance();
            #ifdef antinuDEBUG
                std::cout << "[Reactor::compute_oscillation_constants]: Computing oscillation constants" << std::endl;
                std::cout << "[Reactor::compute_oscillation_constants]: iDm_21_2 = " << iDm_21_2 << ", iDm_32_2 = " << iDm_32_2 << ", iS_12_2 = " << iS_12_2 << ", iS_13_2 = " << iS_13_2 << std::endl;
            #endif

            // Upacks variables
            double fDmSqr21 = Vars->val(iDm_21_2);
            double fDmSqr32 = Vars->val(iDm_32_2);
            double fSSqrTheta12 = Vars->val(iS_12_2);
            double fSSqrTheta13 = Vars->val(iS_13_2);
            #ifdef antinuDEBUG
                std::cout << "[Reactor::compute_oscillation_constants]: fDmSqr21 = " << fDmSqr21 << ", fDmSqr32 = " << fDmSqr32
                        << ", fSSqrTheta12 = " << fSSqrTheta12 << ", fSSqrTheta13 = " << fSSqrTheta13 << std::endl;
            #endif

            // Do calculation
            double fDmSqr31 = fDmSqr32 + fDmSqr21;

            H_ee_vac = fDmSqr21 * (fSSqrTheta12 * (1 - fSSqrTheta13) - (1.0/3.0)) + fDmSqr31 * (fSSqrTheta13 - (1.0/3.0));
            a0_vac = (fDmSqr21*fDmSqr21*fDmSqr21 + fDmSqr31*fDmSqr31*fDmSqr31) / 27.0
                     - (fDmSqr21*fDmSqr21 * fDmSqr31 + fDmSqr21 * fDmSqr31*fDmSqr31) / 18.0;
            a1_vac = (fDmSqr21*fDmSqr21 + fDmSqr31*fDmSqr31 - fDmSqr21 * fDmSqr31) / 9.0;
            Y_ee_vac = (fDmSqr21*fDmSqr21 * (fSSqrTheta12 * (1 - fSSqrTheta13) - (1.0/3.0)) + fDmSqr31*fDmSqr31 * (fSSqrTheta13 - (1.0/3.0))
                        + 2.0 * fDmSqr21 * fDmSqr31 * ((1 - fSSqrTheta12) * (1 - fSSqrTheta13) - (1.0/3.0))) / 3.0;
            
            #ifdef antinuDEBUG
                std::cout << "[Reactor::compute_oscillation_constants]: H_ee_vac = " << H_ee_vac << ", a0_vac = " << a0_vac << ", a1_vac = " << a1_vac << ", Y_ee_vac = " << Y_ee_vac << std::endl;
            #endif
        }

        void re_compute_consts(const double& E) {
            double A_CC = alpha * E; // for A_CC in [eV^2] and nuE in [MeV]
            
            // Compute new values for H_ee, Y, a0 and a1
            double alpha_1 = H_ee_vac * A_CC / 3.0 + A_CC*A_CC / 9.0;

            double a0 = a0_vac + 0.5 * Y_ee_vac * A_CC + H_ee_vac * A_CC*A_CC / 6.0 + A_CC*A_CC*A_CC / 27.0;
            double a1 = a1_vac + alpha_1;
            double Y_ee = Y_ee_vac + 2.0 * alpha_1;
            double H_ee = H_ee_vac + (2.0/3.0) * A_CC;

            // Get eigenvalues of H, and constants X and theta
            double sqrt_a1 = sqrt(a1);
            double arcCos = (1.0/3.0) * acos(a0 / (a1 * sqrt_a1));

            for (int i = 0; i < 3; ++i) {
                eigen.at(i) = 2.0 * sqrt_a1 * cos(arcCos - (2.0/3.0) * M_PI * i);
                X_mat.at(i) = (1.0 + (eigen.at(i) * H_ee + Y_ee) / (eigen.at(i)*eigen.at(i) - a1)) / 3.0;
            }

            #ifdef SUPER_DEBUG
                std::cout << "[Reactor::re_compute_consts]: A_CC = " << A_CC << ", H_ee = " << H_ee << ", a0 = " << a0 << ", a1 = " << a1 << ", Y_ee = " << Y_ee << std::endl;
                std::cout << "[Reactor::re_compute_consts]: eigen.at(0) = " << eigen.at(0) << ", eigen.at(1) = " << eigen.at(1) << ", eigen.at(2) = " << eigen.at(2) << ", X_mat.at(0) = "
                          << X_mat.at(0) << ", X_mat.at(1) = " << X_mat.at(1) << ", X_mat.at(2) = " << X_mat.at(2) << std::endl;
            #endif
        }

        double survival_prob(const double& E, const double& L) {

            const double scale = 1.267E3 * L / E; // for E in [MeV] and L in [km]

            const double s_10 = sin(scale * (eigen.at(1) - eigen.at(0)));
            const double s_20 = sin(scale * (eigen.at(2) - eigen.at(0)));
            const double s_21 = sin(scale * (eigen.at(2) - eigen.at(1)));

            // Compute probability
            double prob = 1.0 - 4.0 * (X_mat.at(1)*X_mat.at(0)*s_10*s_10 + X_mat.at(2)*X_mat.at(0)*s_20*s_20 + X_mat.at(2)*X_mat.at(1)*s_21*s_21);

            #ifdef SUPER_DEBUG
                std::cout << "[Reactor::survival_prob]: E = " << E << ", L = " << L << ", s_10 = " << s_10 << ", s_20 = " << s_20 << ", s_21 = " << s_21 << ", prob = " << prob << std::endl;
            #endif

            return prob;
        }

        // Geo-nu survival probability, used for most distant reactors: averaged over baseline -> No E-depence, only depends on mixing angles.
        void geoNu_survival_prob() {
            FitVars* Vars = FitVars::GetInstance();
            av_survival_prob = Vars->val(iS_13_2)*Vars->val(iS_13_2) + (1. - Vars->val(iS_13_2))*(1. - Vars->val(iS_13_2)) * (1. - 2. * Vars->val(iS_12_2) * (1. - Vars->val(iS_12_2)));
        }

        void hold_osc_params_const(const bool isTrue) {
            FitVars* Vars = FitVars::GetInstance();
            #ifdef antinuDEBUG
                std::cout << "[Reactor::hold_osc_params_const]: Vars->GetNumVars() = " << Vars->GetNumVars() << std::endl;
                std::cout << "[Reactor::hold_osc_params_const]: iDm_21_2 = " << iDm_21_2 << ", iDm_32_2 = " << iDm_32_2 << ", iS_12_2 = " << iS_12_2 << ", iS_13_2 = " << iS_13_2 << std::endl;
            #endif
            if (!isTrue) {
                #ifdef antinuDEBUG
                    std::cout << "[Reactor::hold_osc_params_const]: NOT holding oscillation parameters constant" << std::endl;
                #endif
                computed_osc_specs = false;
            } else if (Vars->IsConstant(iDm_21_2) && Vars->IsConstant(iDm_32_2) && Vars->IsConstant(iS_12_2) && Vars->IsConstant(iS_13_2)) {
                #ifdef antinuDEBUG
                    std::cout << "[Reactor::hold_osc_params_const]: Holding oscillation parameters constant" << std::endl;
                #endif
                geoNu_survival_prob();
                compute_osc_specs();
                computed_osc_specs = true;
            } else {
                std::cout << "[Reactor] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
            }
        }
    
        TH1D* GetModelNoEsys() {return model_noEsys;}
        TH1D* GetModelEsys() {return model_Esys;}

        bool IsInit() {return isInit;}
        void SetBinLims(const unsigned int MinBin, const unsigned int MaxBin) {
            iMinBin = MinBin;
            iMaxBin = MaxBin;

            // re-normalise 1D PDFs
            PWR_promptE_hist->Scale(1. / PWR_promptE_hist->Integral(iMinBin, iMaxBin));
            PWR_Enu_hist->Scale(1. / PWR_Enu_hist->Integral(iMinBin, iMaxBin));
            PHWR_Enu_hist->Scale(1. / PHWR_Enu_hist->Integral(iMinBin, iMaxBin));

            // re-normalise 2D PDF
            double integ;
            for (unsigned int iEnu = 1; iEnu <= E_conv->GetXaxis()->GetNbins(); ++iEnu) {
                integ = E_conv->Integral(iEnu, iEnu, iMinBin, iMaxBin);
                for (unsigned int iBin = 0; iBin <= E_conv->GetYaxis()->GetNbins(); ++iBin) {
                    E_conv->SetBinContent(iEnu, iBin, E_conv->GetBinContent(iEnu, iBin) / integ);
                }
            }
        }
};

// Static methods should be defined outside the class.
Reactor* Reactor::ReactorInstance_{nullptr};

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
Reactor* Reactor::GetInstance() {
    if (ReactorInstance_ == nullptr) {
        ReactorInstance_ = new Reactor();
    }
    return ReactorInstance_;
}


//end header guard
#endif