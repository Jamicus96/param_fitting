#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooConstVar.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TVector3.h>
#include <TTree.h>
#include <sstream>
#include <RAT/DB.hh>
#include <TRandom3.h>

using namespace RooFit;


/**
 * @brief Computes list of constants needed for neutrino oscillation computations.
 * 
 * @param fDmSqr21 
 * @param fDmSqr32 
 * @param fSSqrTheta12 
 * @param fSSqrTheta13 
 * @return std::vector<double> const 
 */
std::vector<double> compute_oscillation_constants(const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13) {

    const double fDmSqr31 = fDmSqr32 + fDmSqr21;

    const double H_ee_vac = fDmSqr21 * (fSSqrTheta12 * (1 - fSSqrTheta13) - (1.0/3.0)) + fDmSqr31 * (fSSqrTheta13 - (1.0/3.0));
    const double H_neq2 = (1 - fSSqrTheta13) * (fDmSqr21*fDmSqr21 * fSSqrTheta12 * (1 + fSSqrTheta12 * (fSSqrTheta13 - 1))
                          + fDmSqr31*fDmSqr31 * fSSqrTheta13 - 2.0 * fDmSqr21 * fDmSqr31 * fSSqrTheta12 * fSSqrTheta13);

    const double a0_vac = - (2.0/27.0) * (fDmSqr21*fDmSqr21*fDmSqr21 + fDmSqr31*fDmSqr31*fDmSqr31)
                          + (1.0/9.0) * (fDmSqr21*fDmSqr21 * fDmSqr31 + fDmSqr21 * fDmSqr31*fDmSqr31);
    const double a1_vac = (1.0/3.0) * (fDmSqr21 * fDmSqr31 - fDmSqr21*fDmSqr21 - fDmSqr31*fDmSqr31);
    const double Y_ee_vac = (2.0/3.0) * a1_vac + H_ee_vac*H_ee_vac + H_neq2;

    // Define electron density of the crust, based on 2.7g/cm3 mass density, and <N/A> = 0.5
    const double Ne = 8.13e23;
    const double alpha = - 2.535e-31 * Ne;  // conversion factor in eV2/MeV

    return {H_ee_vac, a0_vac, a1_vac, Y_ee_vac, alpha};
}

/**
 * @brief Computes survival probability for an electron antineutrino, in curst with constant matter density.
 * 
 * @param E  antineutrino recon energy (MeV)
 * @param L  baseline (km)
 * @param params = {H_ee_vac, a0_vac, a1_vac, Y_ee_vac, alpha}
 * @return double 
 */
double survival_prob(const double E, const double L, const std::vector<double>& params) {

    double scale = 1.267e3 * L / E; // for E in [MeV] and L in [km]
    double A_CC = params[4] * L; // for A_CC in [eV^2] and nuE in [MeV]
    
    // Compute new values for H_ee, Y, a0 and a1 (make sure and Y are updated after their use by others)
    double alpha_1 = params[0] * A_CC + (1.0/3.0) * A_CC*A_CC;

    double a0 = params[1] - params[3] * A_CC - (1.0/3.0) * params[0] * A_CC*A_CC - (2.0/27.0) * A_CC*A_CC*A_CC;
    double a1 = params[2] - alpha_1;
    double Y_ee = params[3] + (2.0/3.0) * alpha_1;
    double H_ee = params[0] + (2.0/3.0) * A_CC;

    // Get eigenvalues of H, and constants X and theta
    double arcCos = (1.0/3.0) * acos(1.5 * (a0/a1) * sqrt(- 3.0 / a1));
    double preFact = 2.0 * sqrt(- a1 / 3.0);

    double eigen[3];
    double X[3];
    for(int i=0; i<3; ++i){
        eigen[i] = preFact * cos(arcCos - (2.0/3.0) * M_PI * i);
        X[i] = (1.0/3.0) + (eigen[i] * H_ee + Y_ee) / (3.0 * eigen[i]*eigen[i] + a1);
    }

    double s_10 = sin(scale * (eigen[1] - eigen[0]));
    double s_20 = sin(scale * (eigen[2] - eigen[0]));
    double s_21 = sin(scale * (eigen[2] - eigen[1]));

    // Compute probability
    return 4.0 * (X[1]*X[0]*s_10*s_10 + X[2]*X[0]*s_20*s_20 + X[2]*X[1]*s_21*s_21);
}


int main(int argv, char** argc) {

    // Get oscillation paramters from ratdb
    // RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    // Double_t fDmSqr21 = linkdb->GetD("deltamsqr21");
    // Double_t fDmSqr32 = linkdb->GetD("deltamsqr32");
    // Double_t fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
    // Double_t fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    return 0;
}