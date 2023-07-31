#include "fitting_utils.hpp"



/**
 * @brief Construct a new reactorINFO::reactorINFO object. Only builds basic data structure.
 * Needs oscillation parameters to produce useful things.
 * 
 * @param Reactor_frac 
 * @param AlphaN_hist 
 * @param reactor_hists 
 * @param L 
 */
reactorINFO::reactorINFO(std::vector<TH1D*>& Reactor_hists, const double N_IBDs, const double IBD_errs, TH2D& E_conv_hist) {

    // Define electron density Ne of the crust, based on 2.7g/cm3 mass density, and <N/A> = 0.5
    alpha = - 2.535e-31 * 8.13e23;  // conversion factor in eV2/MeV * Ne = 8.13e23

    // Set up empty vectors
    for (unsigned int i = 0; i < 3; ++i) {
        eigen.push_back(0.0);
        X_mat.push_back(0.0);
    }

    // Assign values
    reactor_hists = Reactor_hists;
    num_reactors = reactor_hists.size();
    hists_Nbins = reactor_hists.at(0)->GetXaxis()->GetNbins();
    N_IBD = N_IBDs;
    IBD_err = IBD_errs;
    E_conv = *(TH2D*)(E_conv_hist.Clone());
    std::cout << "Set up E conversion histogram:" << std::endl;
    std::cout << "NbinsX = " << E_conv.GetXaxis()->GetNbins() << ", NbinsY = " << E_conv.GetYaxis()->GetNbins() << std::endl;

    tot_hist_int = 0.0;
    for (unsigned int i = 0; i < num_reactors; ++i) {
        tot_hist_int += reactor_hists.at(i)->Integral();
        baselines.push_back(0.);
    }

    // Create histogram to sum oscillated reactor events
    for (unsigned int i = 0; i < 4; ++i) {
        osc_hists.push_back((TH1D*)(reactor_hists.at(0)->Clone()));
        norms.push_back(0.);
    }
    osc_hists.at(0)->SetName("BRUCE");      osc_hists.at(0)->SetTitle("BRUCE spectrum");
    osc_hists.at(1)->SetName("DARLINGTON"); osc_hists.at(1)->SetTitle("DARLINGTON spectrum");
    osc_hists.at(2)->SetName("PICKERING");  osc_hists.at(2)->SetTitle("PICKERING spectrum");
    osc_hists.at(3)->SetName("WORLD");      osc_hists.at(3)->SetTitle("WORLD spectrum");
}

/**
 * @brief Construct a new reactorINFO::reactorINFO object. Only builds basic data structure.
 * Needs oscillation parameters to produce useful things.
 * 
 * @param Reactor_hists 
 * @param Baselines 
 * @param DmSqr21 
 * @param DmSqr32 
 * @param SSqrTheta12 
 * @param SSqrTheta13 
 */
reactorINFO::reactorINFO(std::vector<TH1D*>& Reactor_hists, const double N_IBDs, const double IBD_errs, TH2D& E_conv_hist, const double DmSqr21, const double DmSqr32,
                        const double SSqrTheta12, const double SSqrTheta13) : reactorINFO::reactorINFO(Reactor_hists, N_IBDs, IBD_errs, E_conv_hist) {

    this->Dm21_2() = DmSqr21;
    this->Dm32_2() = DmSqr32;
    this->s12_2() = SSqrTheta12;
    this->s13_2() = SSqrTheta13;
}

/**
 * @brief Construct a new reactorINFO::reactorINFO object. Only builds basic data structure.
 * Needs oscillation parameters to produce useful things.
 * 
 * @param Reactor_hists 
 * @param Baselines 
 * @param DmSqr32 
 * @param SSqrTheta13 
 */
reactorINFO::reactorINFO(std::vector<TH1D*>& Reactor_hists, const double N_IBDs, const double IBD_errs, TH2D& E_conv_hist, const double DmSqr32, const double SSqrTheta13)
                        : reactorINFO::reactorINFO(Reactor_hists, N_IBDs, IBD_errs, E_conv_hist) {

    this->Dm32_2() = DmSqr32;
    this->s13_2() = SSqrTheta13;
}


/**
 * @brief Set oscillation parameters that don't depend on energy or baseline
 * 
 */
void reactorINFO::compute_oscillation_constants() {

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
void reactorINFO::re_compute_consts(const double E) {
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
double reactorINFO::survival_prob(const double E, const double L) {

    const double scale = 1.267e3 * L / E; // for E in [MeV] and L in [km]

    const double s_10 = sin(scale * (eigen.at(1) - eigen.at(0)));
    const double s_20 = sin(scale * (eigen.at(2) - eigen.at(0)));
    const double s_21 = sin(scale * (eigen.at(2) - eigen.at(1)));

    // Compute probability
    return 4.0 * (X_mat.at(1)*X_mat.at(0)*s_10*s_10 + X_mat.at(2)*X_mat.at(0)*s_20*s_20 + X_mat.at(2)*X_mat.at(1)*s_21*s_21);
}


/**
 * @brief Use oscillation parameters and reactor information to compute oscillated histograms
 * 
 */
void reactorINFO::compute_osc_reactor_spec() {
    
    // Compute oscillation constants
    this->compute_oscillation_constants();

    // Reset total reactor histograms (i.e. empty them)
    for (unsigned int i = 0; i < osc_hists.size(); ++i) {
        osc_hists.at(i)->Reset("ICES");
        norms.at(i) = 0;
    }

    // Assume all the histograms have the same E binning
    double E;
    unsigned int hist_idx;
    double weighted_av_P;
    double weight;
    std::string origin_reactor;
    for (unsigned int i = 1; i <= hists_Nbins; ++i) {
        E = reactor_hists.at(0)->GetXaxis()->GetBinCenter(i);
        this->re_compute_consts(E);
        for (unsigned int j = 0; j < num_reactors; ++j) {
            origin_reactor = SplitString(reactor_hists.at(j)->GetName())[0];
            // Find which oscillated histogram to add to
            if (origin_reactor == "BRUCE")           hist_idx = 0;
            else if (origin_reactor == "DARLINGTON") hist_idx = 1;
            else if (origin_reactor == "PICKERING")  hist_idx = 2;
            else                                     hist_idx = 3;

            // Integrate survival prob over E_nu, weigther by E_nu(E_e) distribution -> weigthed average survival prob
            weighted_av_P = 0.0;
            for (unsigned int k = 1; k <= E_conv.GetYaxis()->GetNbins(); ++k) {
                weight = E_conv.GetBinContent(i, k);
                if (weight > 0.0) weighted_av_P += weight * survival_prob(E, baselines.at(j));
            }

            // Add value of bin from reactor core to appropriate hist, scaled by survival probability
            osc_hists.at(hist_idx)->AddBinContent(i, weighted_av_P * reactor_hists.at(j)->GetBinContent(i));
        }
    }

    // Compute hist norms, rescaling using N_IBD
    for (unsigned int i = 0; i < osc_hists.size(); ++i) {
        norms.at(i) = osc_hists.at(i)->Integral() * N_IBD / tot_hist_int;
    }
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
TVector3 LLAtoECEF(double longitude, double latitude, double altitude) {
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
double GetReactorDistanceLLA(const double &longitude, const double&latitude, const double &altitude) {
    const TVector3 SNO_ECEF_coord_ = TVector3(672.87,-4347.18,4600.51);
    double dist = (LLAtoECEF(longitude, latitude,altitude) - SNO_ECEF_coord_).Mag();
  return dist;
}

/**
 * @brief Get vector of reactor baslines [km] from vector of reactor histograms
 * 
 * @param db 
 * @param reactor_hists 
 * @param baselines
 */
void reactorINFO::compute_baselines(RAT::DB* db) {

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
