#include "model_Reactor.hpp"


Reactor::Reactor(const Reactor& mod) {
    Vars = mod.Vars; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec;
    numVars = mod.numVars; iDm_21_2 = mod.iDm_21_2; iDm_32_2 = mod.iDm_32_2; iS_12_2 = mod.iS_12_2; iS_13_2 = mod.iS_13_2;
    iNorms = mod.iNorms; H_ee_vac = mod.H_ee_vac; a0_vac = mod.a0_vac; a1_vac = mod.a1_vac; Y_ee_vac = mod.Y_ee_vac;
    eigen = mod.eigen; X_mat = mod.X_mat; baselines = mod.baselines; reactor_hists = mod.reactor_hists; num_reactors = mod.num_reactors;
    hists_Nbins = mod.hists_Nbins; E_conv = mod.E_conv; reactor_names = mod.reactor_names; reactor_idx = mod.reactor_idx;
    osc_hists = mod.osc_hists; unosc_hist_ints = mod.unosc_hist_ints; computed_osc_specs = mod.computed_osc_specs; db = mod.db;
}
void Reactor::operator = (const Reactor& mod) {
    Vars = mod.Vars; vars = mod.vars; numVars = mod.numVars; model_spec = mod.model_spec;
    numVars = mod.numVars; iDm_21_2 = mod.iDm_21_2; iDm_32_2 = mod.iDm_32_2; iS_12_2 = mod.iS_12_2; iS_13_2 = mod.iS_13_2;
    iNorms = mod.iNorms; H_ee_vac = mod.H_ee_vac; a0_vac = mod.a0_vac; a1_vac = mod.a1_vac; Y_ee_vac = mod.Y_ee_vac;
    eigen = mod.eigen; X_mat = mod.X_mat; baselines = mod.baselines; reactor_hists = mod.reactor_hists; num_reactors = mod.num_reactors;
    hists_Nbins = mod.hists_Nbins; E_conv = mod.E_conv; reactor_names = mod.reactor_names; reactor_idx = mod.reactor_idx;
    osc_hists = mod.osc_hists; unosc_hist_ints = mod.unosc_hist_ints; computed_osc_specs = mod.computed_osc_specs; db = mod.db;
}

Reactor::Reactor(FitVar* vDm_21_2, FitVar* vDm_32_2, FitVar* vS_12_2, FitVar* vS_13_2, const std::vector<TH1D*>& Reactor_hists,
                TH2D* E_conv_hist, std::vector<std::string>& Reactor_names, const std::vector<FitVar*>& Norms, Esys* E_syst, RAT::DB* DB) {

    // Save variables
    Vars.push_back(vDm_21_2); iDm_21_2 = 0;
    Vars.push_back(vDm_32_2); iDm_32_2 = 1;
    Vars.push_back(vS_12_2); iS_12_2 = 2;
    Vars.push_back(vS_13_2); iS_13_2 = 3;
    E_systs.push_back(E_syst);

    if (Reactor_names.size() != Norms.size()) {
        std::cout << "[Reactor] ERROR: Reactor_names and Norms should have the same size!" << std::endl;
        exit(1);
    }
    reactor_names = Reactor_names;

    for (unsigned int i = 0; i < Norms.size(); ++i) {
        Vars.push_back(Norms.at(i));
        iNorms.push_back(4 + i);
    }
    numVars = Vars.size();

    vars.resize(numVars);
    for (unsigned int i = 0; i < numVars; ++i) {
        vars.at(i) = Vars.at(i)->val();
    }

    db = DB;

    // Save hist variables
    reactor_hists = Reactor_hists;
    num_reactors = reactor_hists.size();
    hists_Nbins = reactor_hists.at(0)->GetXaxis()->GetNbins();
    E_conv = (TH2D*)(E_conv_hist->Clone());

    baselines.resize(num_reactors);

    // Set up empty vectors
    eigen.resize(3);
    X_mat.resize(3);

    // Create histogram to sum oscillated reactor events, based off reactor_names
    model_spec = (TH1D*)(reactor_hists.at(0)->Clone());
    for (unsigned int i = 0; i < reactor_names.size(); ++i) {
        osc_hists.push_back((TH1D*)(reactor_hists.at(0)->Clone()));

        osc_hists.at(i)->SetName(reactor_names.at(i).c_str());
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
    Reactor::compute_baselines();

    // Compute unoscillated histogram integrals for the different sources in Reactor_names
    this->compute_unosc_integrals();
}

void Reactor::compute_unosc_integrals() {
    for (unsigned int i = 0; i < unosc_hist_ints.size(); ++i) {
        unosc_hist_ints.at(i) = 0;
    }
    for (unsigned int i = 0; i < num_reactors; ++i) {
        unosc_hist_ints.at(reactor_idx.at(i)) += reactor_hists.at(i)->Integral();
    }
}

/**
 * @brief Set oscillation parameters that don't depend on energy or baseline
 * 
 */
void Reactor::compute_oscillation_constants() {
    #ifdef DEBUG
        std::cout << "[Reactor::compute_oscillation_constants]: Computing oscillation constants" << std::endl;
    #endif

    // Upacks variables
    const double fDmSqr21 = vars.at(iDm_21_2);
    const double fDmSqr32 = vars.at(iDm_32_2);
    const double fSSqrTheta12 = vars.at(iS_12_2);
    const double fSSqrTheta13 = vars.at(iS_13_2);
    #ifdef DEBUG
        std::cout << "[Reactor::compute_oscillation_constants]: fDmSqr21 = " << fDmSqr21 << ", fDmSqr32 = " << fDmSqr32
                  << ", fSSqrTheta12 = " << fSSqrTheta12 << ", fSSqrTheta13 = " << fSSqrTheta13 << std::endl;
    #endif

    // Do calculation
    const double fDmSqr31 = fDmSqr32 + fDmSqr21;

    H_ee_vac = fDmSqr21 * (fSSqrTheta12 * (1 - fSSqrTheta13) - (1.0/3.0)) + fDmSqr31 * (fSSqrTheta13 - (1.0/3.0));
    a0_vac = - (2.0/27.0) * (fDmSqr21*fDmSqr21*fDmSqr21 + fDmSqr31*fDmSqr31*fDmSqr31)
                          + (1.0/9.0) * (fDmSqr21*fDmSqr21 * fDmSqr31 + fDmSqr21 * fDmSqr31*fDmSqr31);
    a1_vac = (1.0/3.0) * (fDmSqr21 * fDmSqr31 - fDmSqr21*fDmSqr21 - fDmSqr31*fDmSqr31);
    Y_ee_vac = (2.0/3.0) * a1_vac + H_ee_vac*H_ee_vac + (1 - fSSqrTheta13) * (fDmSqr21*fDmSqr21 * fSSqrTheta12 * (1 + fSSqrTheta12 * (fSSqrTheta13 - 1))
                        + fDmSqr31*fDmSqr31 * fSSqrTheta13 - 2.0 * fDmSqr21 * fDmSqr31 * fSSqrTheta12 * fSSqrTheta13);;
}


/**
 * @brief Compute oscillation parameters depending on Energy [MeV] but not baseline yet
 * 
 * @param E Energy [MeV]
 */
void Reactor::re_compute_consts(const double& E) {
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

/**
 * @brief Survival probability for energy (MeV) and baseline (km), with oscillation parameters pre-computed.
 * 
 * @param E Energy (MeV)
 * @param L Baseline (km)
 * @return double 
 */
double Reactor::survival_prob(const double& E, const double& L) {

    const double scale = 1.267e3 * L / E; // for E in [MeV] and L in [km]

    const double s_10 = sin(scale * (eigen.at(1) - eigen.at(0)));
    const double s_20 = sin(scale * (eigen.at(2) - eigen.at(0)));
    const double s_21 = sin(scale * (eigen.at(2) - eigen.at(1)));

    // Compute probability
    return 4.0 * (X_mat.at(1)*X_mat.at(0)*s_10*s_10 + X_mat.at(2)*X_mat.at(0)*s_20*s_20 + X_mat.at(2)*X_mat.at(1)*s_21*s_21);
}

/**
 * @brief Get vector of reactor baslines [km] from vector of reactor histograms
 * 
 * @param db 
 * @param reactor_hists 
 * @param baselines
 */
void Reactor::compute_baselines() {

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


void Reactor::compute_osc_specs() {
    // Compute oscillation constants
    this->compute_oscillation_constants();

    // Reset total reactor histograms (i.e. empty them)
    for (unsigned int i = 0; i < osc_hists.size(); ++i) {
        #ifdef DEBUG
            std::cout << "[Reactor::compute_osc_specs]: resetting osc_hists.at(" << i << ")" << std::endl;
        #endif
        osc_hists.at(i)->Reset("ICES");
    }

    // Assume all the histograms have the same E binning
    double E;
    double weighted_av_P;
    double weight;
    std::string origin_reactor;
    for (unsigned int i = 1; i <= hists_Nbins; ++i) {
        E = reactor_hists.at(0)->GetXaxis()->GetBinCenter(i);
        this->re_compute_consts(E);
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


void Reactor::hold_osc_params_const(bool isTrue) {
    if (!isTrue) {
        #ifdef DEBUG
            std::cout << "[Reactor::hold_osc_params_const]: NOT holding oscillation parameters constant" << std::endl;
        #endif
        computed_osc_specs = false;
    } else if (Vars.at(iDm_21_2)->isConstant() && Vars.at(iDm_32_2)->isConstant() && Vars.at(iS_12_2)->isConstant() && Vars.at(iS_13_2)->isConstant()) {
        #ifdef DEBUG
            std::cout << "[Reactor::hold_osc_params_const]: Holding oscillation parameters constant" << std::endl;
        #endif
        vars.at(iDm_21_2) = Vars.at(iDm_21_2)->val();
        vars.at(iDm_32_2) = Vars.at(iDm_32_2)->val();
        vars.at(iS_12_2) = Vars.at(iS_12_2)->val();
        vars.at(iS_13_2) = Vars.at(iS_13_2)->val();
        this->compute_osc_specs();
        computed_osc_specs = true;
    } else {
        std::cout << "[Reactor] WARNING: Oscillation constants not held constant for fitting, since they were not set to constants before hand." << std::endl;
    }
}

void Reactor::compute_spec(Double_t* p) {
    this->GetVarValues(p);
    model_spec->Reset("ICES");  // empty it before re-computing it
    model_spec_sys->Reset("ICES");  // empty it before re-computing it

    // If the oscillation constants are being held constant, and the oscillated spectra have already been computed, can skip this expensive step!
    if (!(Vars.at(iDm_21_2)->isConstant() && Vars.at(iDm_32_2)->isConstant() && Vars.at(iS_12_2)->isConstant() && Vars.at(iS_13_2)->isConstant() && computed_osc_specs)) {
        #ifdef SUPER_DEBUG
            std::cout << "[Reactor::compute_spec]: computing oscillation specs" << std::endl;
        #endif
        this->compute_osc_specs();
    }
    
    for (unsigned int i = 0; i < osc_hists.size(); ++i) {
        #ifdef SUPER_DEBUG
            std::cout << "[Reactor::compute_spec]: Adding oscillation spectrum from " << osc_hists.at(i)->GetName() << std::endl;
        #endif
        model_spec->Add(osc_hists.at(i), vars.at(iNorms.at(i)) / unosc_hist_ints.at(i));
    }

    // Apply energy systematics
    E_systs.at(0)->apply_systematics(p, model_spec, model_spec_sys);
}


void Reactor::GetOscReactorHists(std::vector<TH1D*>& rescaled_osc_hists) {
    for (unsigned int i = 0; i < osc_hists.size(); ++i) {
        #ifdef SUPER_DEBUG
            std::cout << "[Reactor::GetOscReactorHists]: filling rescaled_osc_hists.at(" << i << ")" << std::endl;
        #endif
        rescaled_osc_hists.push_back((TH1D*)(osc_hists.at(i)->Clone()));
        rescaled_osc_hists.at(rescaled_osc_hists.size() - 1)->Reset("ICES");
        rescaled_osc_hists.at(rescaled_osc_hists.size() - 1)->Add(osc_hists.at(i), Vars.at(iNorms.at(i))->val() / unosc_hist_ints.at(i));
    }
}

std::vector<std::string>& Reactor::GetReactorNames() {return reactor_names;}

Reactor::~Reactor() {
    // for (auto p : reactor_hists) {delete p;}
    reactor_hists.clear();
    for (auto p : osc_hists) {delete p;}
    osc_hists.clear();
    // delete db;
    // for (auto p : Vars) {delete p;}
    Vars.clear();
}


/* ~~~~~~~~~~~~~~~~ OTHER (NON-MEMBER) FUNCTIONS ~~~~~~~~~~~~~~~~ */


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
  static double toRad = TMath::Pi()/180.;
  static double Earthradius = 6378137.0; //Radius of the Earth (in meters)
  static double f = 1./298.257223563; //Flattening factor WGS84 Model
  static double L, rs, x, y, z;
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

