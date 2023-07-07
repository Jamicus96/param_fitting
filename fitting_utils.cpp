#include "fitting_utils.hpp"



/**
 * @brief Construct a new PDFspec::PDFspec object. Only builds basic data structure.
 * Needs oscillation parameters to produce useful things.
 * 
 * @param Reactor_frac 
 * @param AlphaN_hist 
 * @param reactor_hists 
 * @param L 
 */
PDFspec::PDFspec(TH1D* AlphaN_hist, std::vector<TH1D*>* Reactor_hists, std::vector<double>* Baselines) {

    // Define electron density Ne of the crust, based on 2.7g/cm3 mass density, and <N/A> = 0.5
    alpha = - 2.535e-31 * 8.13e23;  // conversion factor in eV2/MeV * Ne = 8.13e23

    // Set up empty vectors
    for (unsigned int i = 0; i < 3; ++i) {
        eigen.push_back(0.0);
        X_mat.push_back(0.0);
    }

    // Assign vectors
    reactor_hists = Reactor_hists;
    baselines = Baselines;
    num_reactors = reactor_hists->size();
    if (baselines->size() != num_reactors) {
        std::cout << "ERROR: Number of baselines (" << baselines->size() << ") not equal to number of reactors (" << num_reactors << ")!" << std::endl;
        exit(1);
    }
    hists_Nbins = reactor_hists->at(0)->GetXaxis()->GetNbins();

    // Create histogram to sum oscillated reactor events
    tot_reactor_hist = (TH1D*)(reactor_hists->at(0)->Clone());

    //declare observables
    e = RooRealVar("E", "energy", 0.9, 8);
    reactor_frac = RooRealVar("reactor_frac", "reactorIBD fraction", 0.4, 0.8);
    alphaN_frac = RooRealVar("alphaN_frac", "alpha-n fraction", 0.2, 0.6);

    // Create alphaN PDF
    tempData_alpha = new RooDataHist("tempData", "temporary data", e, AlphaN_hist);
    alphaN_PDF = new RooHistPdf("alphaN_PDF", "alphaN PDF", e, *tempData_alpha);
}


/**
 * @brief Set oscillation parameters that don't depend on energy or baseline
 * 
 * @param fDmSqr21 Dm_21^2
 * @param fDmSqr32 Dm_32^2
 * @param fSSqrTheta12 s_12^2
 * @param fSSqrTheta13 s_13^2
 */
void PDFspec::compute_oscillation_constants(const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13) {

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
void PDFspec::re_compute_consts(const double E) {
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
double PDFspec::survival_prob(const double E, const double L) {

    const double scale = 1.267e3 * L / E; // for E in [MeV] and L in [km]

    const double s_10 = sin(scale * (eigen.at(1) - eigen.at(0)));
    const double s_20 = sin(scale * (eigen.at(2) - eigen.at(0)));
    const double s_21 = sin(scale * (eigen.at(2) - eigen.at(1)));

    // Compute probability
    return 4.0 * (X_mat.at(1)*X_mat.at(0)*s_10*s_10 + X_mat.at(2)*X_mat.at(0)*s_20*s_20 + X_mat.at(2)*X_mat.at(1)*s_21*s_21);
}


/**
 * @brief Give PDFspec oscillation parameters, and it will compute the total oscillated reactor spectrum
 * 
 * @param fDmSqr21 
 * @param fDmSqr32 
 * @param fSSqrTheta12 
 * @param fSSqrTheta13 
 */
void PDFspec::compute_tot_reactor_spec(const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13) {
    
    // Compute oscillation constants
    this->compute_oscillation_constants(fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);

    // Reset total reactor histogram (i.e. empty it)
    tot_reactor_hist->Reset("ICES");

    // Assume all the histograms have the same E binning
    double E;
    for (unsigned int i = 1; i < hists_Nbins; ++i) {
        E = reactor_hists->at(0)->GetXaxis()->GetBinCenter(i);
        this->re_compute_consts(E);
        for (unsigned int j = 0; j < num_reactors; ++j) {
            // Add value of bin from reactor core to total, scaled by survival probability
            tot_reactor_hist->AddBinContent(i, survival_prob(E, baselines->at(j)) * reactor_hists->at(j)->GetBinContent(i));
        }
    }

    // Create reactor PDF
    tempData_react = new RooDataHist("tempData", "temporary data", e, tot_reactor_hist);
    reactor_PDF = new RooHistPdf("reactor_PDF", "reactor PDF", e, *tempData_react);
}

/**
 * @brief Fit PDFs to data, output minimum log likelihood
 * 
 * @param dataHist 
 * @return double 
 */
double PDFspec::ML_fit(RooDataHist& dataHist) {

    // Add PDFs together
    *model = RooAddPdf("model", "r+a", RooArgList(*reactor_PDF, *alphaN_PDF), RooArgList(reactor_frac, alphaN_frac));

    // Fit to data
    result = model->fitTo(dataHist, Extended(true), PrintLevel(-1), SumW2Error(kFALSE), Save());

    // Get results
    return result->minNll();
}


/**
 * @brief Destroy the PDFspec::PDFspec object
 * 
 */
PDFspec::~PDFspec() {
    delete(baselines);
    delete(reactor_hists);
    delete(tot_reactor_hist);
    delete(tempData_react);
    delete(reactor_PDF);
    delete(tempData_alpha);
    delete(alphaN_PDF);
    delete(model);
    delete(result);
}
