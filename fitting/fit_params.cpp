#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <algorithm>
// #include <TVector3.h>
#include <TTree.h>
#include <sstream>
#include <RAT/DB.hh>
#include <TRandom3.h>
#include <TKey.h>
#include <TObject.h>
#include <TList.h>
#include <TVectorD.h>
#include "Math/DistFunc.h"

#include "fitter.hpp"
#include "fitVars.hpp"
#include "E_systematics.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"
#include "model_Reactor.hpp"
#include "fitting_utils.hpp"


void Fit_spectra(const std::vector<std::vector<double>>& var_params, TH2D* minllHist, std::vector<TH2D*>& par_hists, const std::vector<unsigned int>& start_idx, const bool verbose);
std::vector<std::vector<double>> make_var_param_vals(const double Dm21_min, const double Dm21_max, const unsigned int Dm21_nSteps, const double Theta12_min,
                                                     const double Theta12_max, const unsigned int Theta12_nSteps);



int main(int argv, char** argc) {

    /* ~~~~~~~~ READING IN ~~~~~~~~ */

    // file args
    std::string PDFs_address = argc[1];
    std::string data_ntuple_address = argc[2];
    std::string out_address = argc[3];

    bool use_Azimov = std::stoi(argc[4]);

    // 2d hist limit args
    double Dm21_lower = std::atof(argc[5]);
    double Dm21_upper = std::atof(argc[6]);
    double Theta12_lower = std::atof(argc[7]);  // degrees
    double Theta12_upper = std::atof(argc[8]);  // degrees
    unsigned int N_bins = std::atoi(argc[9]);

    // variable paramters limit args
    double Dm21_min = std::atof(argc[10]);
    double Dm21_max = std::atof(argc[11]);
    double Theta12_min = std::atof(argc[12]);  // degrees
    double Theta12_max = std::atof(argc[13]);  // degrees

    unsigned int start_idx_Dm21 = std::atoi(argc[14]);
    unsigned int start_idx_theta = std::atoi(argc[15]);

    unsigned int Dm21_nSteps = std::atoi(argc[16]);
    unsigned int Theta12_nSteps = std::atoi(argc[17]);

    bool verbose = std::stoi(argc[18]);

    if (verbose) {
        std::cout << "N_IBD = " << N_IBD << " ± (" << IBD_err_indiv*100. << "% (individual), " << IBD_err_tot*100. << "% (total))." << std::endl;
        std::cout << "N_alphaN = " << N_alphaN << " ± (" << alphaN_err_GS*100. << "% (GS), " << alphaN_err_ES*100. << "% (ES))." << std::endl;
        std::cout << "N_geoNu = " << N_geoNu << " ± " << geoNu_err*100. << "%." << std::endl;
        std::cout << "linScale = 1 ± " << linScale_err << "." << std::endl;
        std::cout << "Beta: kB = " << kB << " ± " << kB_err << "." << std::endl;
        std::cout << "Proton: kB = " << kB_P << " ± " << kB_err_P << "." << std::endl;
        std::cout << "sigPerSqrtE = 0 + " << sigPerSqrtE << "." << std::endl;
        std::cout << "fit params: Dm_21^2 € [" << Dm21_min << ", " << Dm21_max << "] (#" << Dm21_nSteps << "), " << "theta_12^2 € [" << Theta12_min << ", " << Theta12_max << "] (#" << Theta12_nSteps << ")." << std::endl;
        std::cout << "hist params: Dm_21^2 lims € {" << Dm21_lower << ", " << Dm21_upper << "}, " << "theta_12^2 lims € {" << Theta12_lower << ", " << Theta12_upper << "}, Nbins = " << N_bins << "." << std::endl;
    }

    // Get DB
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DB* db = RAT::DB::Get();
    db->LoadDefaults();

    // Get oscillation constants
    std::cout << "Getting oscillation parameters..." << std::endl;
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    const double fDmSqr32 = linkdb->GetD("deltamsqr32");
    const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
    const double fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    /* ~~~~~~~~ FITTING ~~~~~~~~ */

    // Create fitter object
    if (use_Azimov) {
        create_fitter(PDFs_address, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13, db);
    } else {
        TFile* DataFile = TFile::Open(data_ntuple_address.c_str());
        TTree *DataInfo = (TTree *) DataFile->Get("output");
        create_fitter(PDFs_address, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13, db, DataInfo);
    }

    // Make list of Dm_21^2 and s_12^2 values to iterate over
    std::vector<std::vector<double>> var_params = make_var_param_vals(Dm21_min, Dm21_max, Dm21_nSteps, Theta12_min, Theta12_max, Theta12_nSteps);
    std::vector<unsigned int> start_idx = {start_idx_theta, start_idx_Dm21};

    // Set up 2d log-likelihood histogram
    TH2D* minllHist = new TH2D("minllHist", "minimised likelihood values", N_bins, Theta12_lower, Theta12_upper, N_bins, Dm21_lower, Dm21_upper);
    minllHist->GetXaxis()->SetTitle("Theta_12");
    minllHist->GetYaxis()->SetTitle("Delta m_21^2");
    minllHist->GetYaxis()->SetTitleOffset(1);  // Move x-axis label further away from the axis (0 is default)
    minllHist->SetStats(0);  // Remove stats box

    // Set up the same 2d hist for each parameter in the fit
    std::vector<TH2D*> par_hists;
    setup_param_hists(minllHist, par_hists);

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;
    Fit_spectra(var_params, minllHist, par_hists, start_idx, verbose);


    /* ~~~~~~~~ OUTPUT ~~~~~~~~ */

    // Write hist to file and close
    TFile *outroot = new TFile(out_address.c_str(), "RECREATE");
    minllHist->Write();
    for (unsigned int i = 0; i < par_hists.size(); ++i) par_hists.at(i)->Write();
    outroot->Write();
    outroot->Close();
    delete(outroot);

    return 0;
}


/**
 * @brief Takes reactor and alpha-n data, and fits data to it for a range of Dm_21^2 and theta_12 values
 * 
 * @param spectrum  reactorINFO object that contains all the relevent reactor info
 * @param alphaN_hists  alphaN PDFs (hists)
 * @param N_alphaN  number of alphaN expected events
 * @param alphaN_err  fractional error of N_alphaN
 * @param geoNu_hists  geoNu PDFs (hists)
 * @param N_geoNu  number of geoNu expected events
 * @param geoNu_err  fractional error of N_geoNu
 * @param data  data to fit model to 
 * @param var_params  Dm21^2 and s12^2 values to iterate over
 * @param minllHist  2D histogram that min log-likelihood values are dumped in
 */
void Fit_spectra(const std::vector<std::vector<double>>& var_params, TH2D* minllHist, std::vector<TH2D*>& par_hists, const std::vector<unsigned int>& start_idx, const bool verbose) {

    Fitter* antinuFitter = Fitter::GetInstance();
    FitVars* Vars = FitVars::GetInstance();
    Reactor* ReactorMod = Reactor::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();

    // Unpack
    std::vector<double> sinTheta12 = var_params.at(0);
    std::vector<double> Dm21 = var_params.at(1);
    std::cout << "sinTheta12.size() = " << sinTheta12.size() << std::endl;
    std::cout << "Dm21.size() = " << Dm21.size() << std::endl;

    /* ~~~~~  Compute best fit Log likelihood for each set of parameters (fraction of alpha-n vs reactor IBD events is fit in each loop)  ~~~~~ */
    double value;
    std::cout << "Looping over oscillation parameters..." << std::endl;
    for (unsigned int i = 0; i < sinTheta12.size(); ++i) {
        // Set sinTheta12 value
        Vars->val("sinsqrtheta12") = sinTheta12.at(i);
        if (verbose) std::cout << "s_12^2 = " << Vars->val("sinsqrtheta12") << std::endl;
        geoNuMod->hold_osc_params_const(true);  // This will also pre-compute the geo-nu survival prob ahead of time
        for (unsigned int j = 0; j < Dm21.size(); ++j) {
            // Set Dm21 value
            Vars->val("deltamsqr21") = Dm21.at(j);
            if (verbose) std::cout << "Dm_21^2 = " << Vars->val("deltamsqr21") << std::endl;
            ReactorMod->hold_osc_params_const(true); // This will also compute oscillated reactor specs

            // Perform fit and record log-likelihood
            value = antinuFitter->fit_models();
            minllHist->SetBinContent(start_idx.at(0) + i + 1, start_idx.at(1) + j + 1, value);
            std::cout << "minLL = " << value << std::endl;

            // Also record best fit values of all parameters
            for (unsigned int iVar = 0; iVar < Vars->GetNumVars(); ++iVar) {
                value = antinuFitter->GetFitParameter(iVar);
                par_hists.at(iVar)->SetBinContent(start_idx.at(0) + i + 1, start_idx.at(1) + j + 1, value);
                std::cout << par_hists.at(iVar)->GetName() << " = " << value << std::endl;
            }
        }
    }
}


/**
 * @brief Produces list of Dm_21^2 and s_12^2 values (packaged together) based on Dm_21^2 and theta_12 limits.
 * 
 * @param Dm21_min [MeV^2]
 * @param Dm21_max [MeV^2]
 * @param Theta12_min [degrees]
 * @param Theta12_max [degrees]
 * @param N_steps same number of steps in both directions
 * @return const std::vector<std::vector<double>>& = {s_12^2, Dm_21^2}
 */
std::vector<std::vector<double>> make_var_param_vals(const double Dm21_min, const double Dm21_max, const unsigned int Dm21_nSteps, const double Theta12_min, const double Theta12_max, const unsigned int Theta12_nSteps) {

    // Dm21^2
    std::vector<double> Dm21;
    if (Dm21_nSteps == 0) {
        std::cout << "ERROR: Cannot have zero steps! (Dm21_nSteps = 0)" <<std::endl;
        exit(1);
    } else if (Dm21_nSteps == 1) {
        Dm21.push_back(0.5 * (Dm21_max + Dm21_min));
        std::cout << "Dm21.at(0) = " << Dm21.at(0) << std::endl;
    } else {
        double Dm21_step = (Dm21_max - Dm21_min) / (double)(Dm21_nSteps - 1);
        for (unsigned int n = 0; n < Dm21_nSteps; ++n) {
            Dm21.push_back(Dm21_min + (double)n * Dm21_step);
            std::cout << "Dm21.at(" << n << ") = " << Dm21.at(n) << std::endl;
        }
    }
    std::cout << "Dm21.size() = " << Dm21.size() << std::endl;

    // s_12^2
    std::vector<double> sinTheta12;
    if (Theta12_nSteps == 0) {
        std::cout << "ERROR: Cannot have zero steps! (Theta12_nSteps = 0)" <<std::endl;
        exit(1);
    } else if (Theta12_nSteps == 1) {
        sinTheta12.push_back(0.5 * (Theta12_max + Theta12_min));
        std::cout << "sinTheta12.at(0) = " << sinTheta12.at(0) << std::endl;
    } else {
        double Theta12_step = (Theta12_max - Theta12_min) / (double)(Theta12_nSteps - 1);
        for (unsigned int n = 0; n < Theta12_nSteps; ++n) {
            sinTheta12.push_back(pow(sin((Theta12_min + (double)n * Theta12_step)  * TMath::Pi() / 180.), 2));
            std::cout << "sinTheta12.at(" << n << ") = " << sinTheta12.at(n) << std::endl;
        }
    }
    std::cout << "sinTheta12.size() = " << sinTheta12.size() << std::endl;

    // Package them together
    return {sinTheta12, Dm21};
}
