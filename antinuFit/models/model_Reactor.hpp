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

        TH1D* model_noEsys;
        TH1D* model_Esys;
        bool isInit;

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
    
    protected:
        // Constructors/desctructors
        Reactor() : isInit(false) {}
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
                const std::vector<unsigned int>& norms_idx, const unsigned int totNorm_idx, const unsigned int Esys_idx,
                const std::vector<TH1D*>& Reactor_hists, TH2D* E_conv_hist, const std::vector<std::string>& Reactor_names, RAT::DB* DB) {
                        
            if (Reactor_names.size() != norms_idx.size()) {
                std::cout << "[Reactor::InitReactor] ERROR: Reactor_names and Norms should have the same size!" << std::endl;
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
            reactor_names = Reactor_names;
            iNorms = norms_idx;
            iTotNorm = totNorm_idx;
            db = DB;

            #ifdef antinuDEBUG
                std::cout << "[Reactor::InitReactor]: iDm_21_2 = " << iDm_21_2 << ", iDm_32_2 = " << iDm_32_2 << ", iS_12_2 = " << iS_12_2 << ", iS_13_2 = " << iS_13_2 << ", iEsys = " << iEsys << ", iTotNorm = " << iTotNorm << std::endl;
            #endif

            // Save hist variables
            reactor_hists = Reactor_hists;
            num_reactors = reactor_hists.size();
            hists_Nbins = reactor_hists.at(0)->GetXaxis()->GetNbins();
            E_conv = (TH2D*)(E_conv_hist->Clone("Reactor::E_conv"));
            E_conv->Reset("ICES"); E_conv->Add(E_conv_hist);
            baselines.resize(num_reactors);

            // Set up empty vectors
            eigen.resize(3);
            X_mat.resize(3);

            // Create histogram to sum oscillated reactor events, based off reactor_names
            for (unsigned int i = 0; i < reactor_names.size(); ++i) {
                osc_hists.push_back((TH1D*)(reactor_hists.at(0)->Clone(("Reactor::osc_hists::" + reactor_names.at(i)).c_str())));
                osc_hists.at(i)->SetTitle((reactor_names.at(i) + " spectrum").c_str());
            }
            unosc_hist_ints.resize(reactor_names.size());
            
            // Assign each element in reactor_hists to the appropriate element in osc_hists/reactor_names
            std::string origin_reactor;
            bool reactor_not_found;
            for (unsigned int i = 0; i < num_reactors; ++i) {
                origin_reactor = SplitString(reactor_hists.at(i)->GetName())[0];

                // See if reactor is in reactor_names (except last element 'WORLD')
                reactor_not_found = true;
                for (unsigned int j = 0; j < (reactor_names.size()-1); ++j) {
                    if (origin_reactor == reactor_names.at(j)) {
                        reactor_idx.push_back(j);
                        reactor_not_found = false;
                        break;
                    }
                }
                // Reactor was not in the list, so assign it to 'WORLD' (last element in reactor_names)
                if (reactor_not_found) reactor_idx.push_back(reactor_names.size() - 1);
            }

            // Compute baselines
            compute_baselines();

            // Compute unoscillated histogram integrals for the different sources in Reactor_names
            compute_unosc_integrals();

            model_noEsys = (TH1D*)(reactor_hists.at(0)->Clone("Reactor::model_noEsys"));
            model_Esys = (TH1D*)(reactor_hists.at(0)->Clone("Reactor::model_Esys"));

            #ifdef antinuDEBUG
                std::cout << "[Reactor::InitReactor]: iDm_21_2 = " << iDm_21_2 << ", iDm_32_2 = " << iDm_32_2 << ", iS_12_2 = " << iS_12_2 << ", iS_13_2 = " << iS_13_2 << ", iEsys = " << iEsys << ", iTotNorm = " << iTotNorm << std::endl;
            #endif

            isInit = true;
        }
        
        void InitReactor(const std::string Dm21_2_name, const std::string Dm32_2_name, const std::string s12_2_name, const std::string s13_2_name,
                const std::vector<std::string>& norms_names, const std::string totNorm_name, const std::string Esys_name,
                const std::vector<TH1D*>& Reactor_hists, TH2D* E_conv_hist, const std::vector<std::string>& Reactor_names, RAT::DB* DB) {

            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();

            std::vector<unsigned int> norms_idx;
            for (unsigned int i = 0; i < norms_names.size(); ++i) norms_idx.push_back(Vars->findIdx(norms_names.at(i)));
            #ifdef antinuDEBUG
                std::cout << "[Reactor::InitReactor]: Vars->GetNumVars() = " << Vars->GetNumVars() << std::endl;
            #endif

            InitReactor(Vars->findIdx(Dm21_2_name), Vars->findIdx(Dm32_2_name), Vars->findIdx(s12_2_name), Vars->findIdx(s13_2_name), norms_idx,
                    Vars->findIdx(totNorm_name), Esysts->findIdx(Esys_name), Reactor_hists, E_conv_hist, Reactor_names, DB);
        };

        // Member function
        void compute_unosc_integrals() {
            for (unsigned int i = 0; i < unosc_hist_ints.size(); ++i) {
                unosc_hist_ints.at(i) = 0;
            }
            for (unsigned int i = 0; i < num_reactors; ++i) {
                unosc_hist_ints.at(reactor_idx.at(i)) += reactor_hists.at(i)->Integral();
            }
        }

        void compute_osc_specs() {
            // Compute oscillation constants
            compute_oscillation_constants();

            // Reset total reactor histograms (i.e. empty them)
            for (unsigned int i = 0; i < osc_hists.size(); ++i) {
                #ifdef antinuDEBUG
                    std::cout << "[Reactor::compute_osc_specs]: resetting osc_hists.at(" << i << ")" << std::endl;
                #endif
                osc_hists.at(i)->Reset("ICES");
            }

            // Assume all the histograms have the same E binning
            double E;
            double weighted_av_P;
            double weight;
            std::string origin_reactor;
            #ifdef antinuDEBUG
                std::cout << "[Reactor::compute_osc_specs]: hists_Nbins = " << hists_Nbins << std::endl;
                std::cout << "[Reactor::compute_osc_specs]: E_conv->GetXaxis()->GetNbins() = " << E_conv->GetXaxis()->GetNbins() << std::endl;
                std::cout << "[Reactor::compute_osc_specs]: E_conv->GetYaxis()->GetNbins() = " << E_conv->GetYaxis()->GetNbins() << std::endl;
            #endif

            for (unsigned int i = 1; i <= hists_Nbins; ++i) {
                E = reactor_hists.at(0)->GetXaxis()->GetBinCenter(i);
                re_compute_consts(E);
                for (unsigned int j = 0; j < num_reactors; ++j) {
                    // Integrate survival prob over E_nu, weigther by E_nu(E_e) distribution -> weigthed average survival prob
                    weighted_av_P = 0.0;
                    for (unsigned int k = 1; k <= E_conv->GetYaxis()->GetNbins(); ++k) {
                        weight = E_conv->GetBinContent(i, k);
                        if (weight > 0.0) weighted_av_P += weight * survival_prob(E, baselines.at(j));
                    }

                    // Add value of bin from reactor core to appropriate hist, scaled by survival probability
                    osc_hists.at(reactor_idx.at(j))->AddBinContent(i, weighted_av_P * reactor_hists.at(j)->GetBinContent(i));
                }
            }
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
            a0_vac = - (2.0/27.0) * (fDmSqr21*fDmSqr21*fDmSqr21 + fDmSqr31*fDmSqr31*fDmSqr31)
                                + (1.0/9.0) * (fDmSqr21*fDmSqr21 * fDmSqr31 + fDmSqr21 * fDmSqr31*fDmSqr31);
            a1_vac = (1.0/3.0) * (fDmSqr21 * fDmSqr31 - fDmSqr21*fDmSqr21 - fDmSqr31*fDmSqr31);
            Y_ee_vac = (2.0/3.0) * a1_vac + H_ee_vac*H_ee_vac + (1 - fSSqrTheta13) * (fDmSqr21*fDmSqr21 * fSSqrTheta12 * (1 + fSSqrTheta12 * (fSSqrTheta13 - 1))
                                + fDmSqr31*fDmSqr31 * fSSqrTheta13 - 2.0 * fDmSqr21 * fDmSqr31 * fSSqrTheta12 * fSSqrTheta13);;
        }

        void re_compute_consts(const double& E) {
            const double A_CC = alpha * E; // for A_CC in [eV^2] and nuE in [MeV]
            
            // Compute new values for H_ee, Y, a0 and a1
            const double alpha_1 = H_ee_vac * A_CC + (1.0/3.0) * A_CC*A_CC;

            const double a0 = a0_vac - Y_ee_vac * A_CC - (1.0/3.0) * H_ee_vac * A_CC*A_CC - (2.0/27.0) * A_CC*A_CC*A_CC;
            const double a1 = a1_vac - alpha_1;
            const double Y_ee = Y_ee_vac + (2.0/3.0) * alpha_1;
            const double H_ee = H_ee_vac + (2.0/3.0) * A_CC;

            // Get eigenvalues of H, and constants X and theta
            const double arcCos = (1.0/3.0) * acos(1.5 * (a0/a1) * sqrt(- 3.0 / a1));
            const double preFact = 2.0 * sqrt(- a1 / 3.0);

            for (int i = 0; i < 3; ++i) {
                eigen.at(i) = preFact * cos(arcCos - (2.0/3.0) * M_PI * i);
                X_mat.at(i) = (1.0/3.0) + (eigen.at(i) * H_ee + Y_ee) / (3.0 * eigen.at(i)*eigen.at(i) + a1);
            }
        }

        double survival_prob(const double& E, const double& L) {

            const double scale = 1.267e3 * L / E; // for E in [MeV] and L in [km]

            const double s_10 = sin(scale * (eigen.at(1) - eigen.at(0)));
            const double s_20 = sin(scale * (eigen.at(2) - eigen.at(0)));
            const double s_21 = sin(scale * (eigen.at(2) - eigen.at(1)));

            // Compute probability
            return 4.0 * (X_mat.at(1)*X_mat.at(0)*s_10*s_10 + X_mat.at(2)*X_mat.at(0)*s_20*s_20 + X_mat.at(2)*X_mat.at(1)*s_21*s_21);
        }

        void compute_baselines() {

            RAT::DBLinkPtr linkdb;
            std::vector<std::string> originReactorVect;
            std::vector<Double_t> fLatitude;
            std::vector<Double_t> fLongitute;
            std::vector<Double_t> fAltitude;
            for (unsigned int n = 0; n < num_reactors; ++n) {
                originReactorVect = SplitString(reactor_hists.at(n)->GetName());
                linkdb = db->GetLink("REACTOR",originReactorVect[0]);
                fLatitude  = linkdb->GetDArray("latitude");
                fLongitute = linkdb->GetDArray("longitude");
                fAltitude = linkdb->GetDArray("altitude");

                baselines.at(n) = GetReactorDistanceLLA(fLongitute[std::stoi(originReactorVect[1])], fLatitude[std::stoi(originReactorVect[1])], fAltitude[std::stoi(originReactorVect[1])]);
            }
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
            } else if (Vars->isConstant(iDm_21_2) && Vars->isConstant(iDm_32_2) && Vars->isConstant(iS_12_2) && Vars->isConstant(iS_13_2)) {
                #ifdef antinuDEBUG
                    std::cout << "[Reactor::hold_osc_params_const]: Holding oscillation parameters constant" << std::endl;
                #endif
                compute_osc_specs();
                computed_osc_specs = true;
            } else {
                std::cout << "[Reactor] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
            }
        }

        void compute_spec() {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            model_noEsys->Reset("ICES");  // empty it before re-computing it
            model_Esys->Reset("ICES");  // empty it before re-computing it

            // If the oscillation constants are being held constant, and the oscillated spectra have already been computed, can skip this expensive step!
            if (!(Vars->isConstant(iDm_21_2) && Vars->isConstant(iDm_32_2) && Vars->isConstant(iS_12_2) && Vars->isConstant(iS_13_2) && computed_osc_specs)) {
                #ifdef SUPER_DEBUG
                    std::cout << "[Reactor::compute_spec]: computing oscillation specs" << std::endl;
                #endif
                compute_osc_specs();
            }
            
            for (unsigned int i = 0; i < osc_hists.size(); ++i) {
                #ifdef SUPER_DEBUG
                    std::cout << "[Reactor::compute_spec]: Adding oscillation spectrum from " << osc_hists.at(i)->GetName() << std::endl;
                #endif
                model_noEsys->Add(osc_hists.at(i), Vars->val(iTotNorm) * Vars->val(iNorms.at(i)) / unosc_hist_ints.at(i));
            }

            // Apply energy systematics
            Esysts->apply_systematics(iEsys, model_noEsys, model_Esys);
        }

        std::vector<std::string>& GetReactorNames() {return reactor_names;}

        void Spectra(std::vector<TH1D*>& hists) {
            FitVars* Vars = FitVars::GetInstance();
            Esys* Esysts = Esys::GetInstance();
            TH1D* temp_hist = (TH1D*)(osc_hists.at(0)->Clone("temp_hist"));

            for (unsigned int i = 0; i < osc_hists.size(); ++i) {
                #ifdef SUPER_DEBUG
                    std::cout << "[Reactor::GetOscReactorHists]: filling hists.at(" << i << ")" << std::endl;
                #endif
                temp_hist->Reset("ICES");
                temp_hist->Add(osc_hists.at(i), Vars->val(iNorms.at(i)) / unosc_hist_ints.at(i));
                hists.push_back((TH1D*)(osc_hists.at(i)->Clone(("model_" + reactor_names.at(i)).c_str())));
                hists.at(hists.size()-1)->Reset("ICES");
                Esysts->apply_systematics(iEsys, temp_hist, hists.at(hists.size()-1));
            }
        }

    
        TH1D* GetModelNoEsys() {return model_noEsys;}
        TH1D* GetModelEsys() {return model_Esys;}

        bool IsInit() {return isInit;}
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


/**
 * @brief Split reactor name strings for getting baselines
 * 
 * @param str 
 * @return std::vector<std::string> 
 */
std::vector<std::string> SplitString(std::string str) {
    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end); //each word of string now in vector
    std::vector<std::string> info;
    std::string dummy = "";

    for(int i=0;i<tokens.size();i++){ //combine back to reactor name, core number
        if(i==0){
            dummy = tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else if(i!=tokens.size()-1){
            dummy += tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else{
            info.push_back(dummy);
            info.push_back(tokens.at(i));
        }
    }

    return info;
}

/**
 * @brief Coordinate conversion time, to get baselines
 * 
 * @param longitude 
 * @param latitude 
 * @param altitude 
 * @return TVector3 
 */
TVector3 LLAtoECEF(const double& longitude, const double& latitude, const double& altitude) {
  // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
  double toRad = TMath::Pi()/180.;
  double Earthradius = 6378137.0; //Radius of the Earth (in meters)
  double f = 1./298.257223563; //Flattening factor WGS84 Model
  double L, rs, x, y, z;
  L = atan( pow((1. - f),2)*tan(latitude*toRad))*180./TMath::Pi();
  rs = sqrt( pow(Earthradius,2)/(1. + (1./pow((1. - f),2) - 1.)*pow(sin(L*toRad),2)));
  x = (rs*cos(L*toRad)*cos(longitude*toRad) + altitude*cos(latitude*toRad)*cos(longitude*toRad))/1000; // in km
  y = (rs*cos(L*toRad)*sin(longitude*toRad) + altitude*cos(latitude*toRad)*sin(longitude*toRad))/1000; // in km
  z = (rs*sin(L*toRad) + altitude*sin(latitude*toRad))/1000; // in km

  TVector3 ECEF = TVector3(x,y,z);

  return ECEF;
}

/**
 * @brief Compute reactor baseline [km]
 * 
 * @param longitude 
 * @param latitude 
 * @param altitude 
 * @return double 
 */
double GetReactorDistanceLLA(const double& longitude, const double& latitude, const double &altitude) {
    const TVector3 SNO_ECEF_coord_ = TVector3(672.87,-4347.18,4600.51);
    double dist = (LLAtoECEF(longitude, latitude, altitude) - SNO_ECEF_coord_).Mag();
  return dist;
}



//end header guard
#endif