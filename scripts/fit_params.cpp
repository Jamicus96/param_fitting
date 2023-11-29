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
#include "fitVar.hpp"
#include "model.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"
#include "model_Reactor.hpp"


void Fit_spectra(Fitter& antinuFitter, FitVar& vDm21_2, FitVar& vS_12_2, Model& geoNuMod, Model& ReactorMod, const std::vector<std::vector<double>>& var_params, TH2D* minllHist, const std::vector<unsigned int>& start_idx, const bool verbose);
std::vector<double>& compute_unosc_fracs(const std::vector<TH1D*>& Reactor_hists, const std::vector<std::string>& Reactor_names, std::vector<double>& hist_fracs);
std::vector<std::vector<double>> make_var_param_vals(const double Dm21_min, const double Dm21_max, const unsigned int Dm21_nSteps, const double Theta12_min, const double Theta12_max, const unsigned int Theta12_nSteps);
void read_hists_from_file(std::string file_address, std::vector<TH1D*>& reactor_hists, std::vector<TH1D*>& alphaN_hists, std::vector<TH1D*>& geoNu_hists, TH2D& E_conv);


int main(int argv, char** argc) {
    // file args
    std::string PDFs_address = argc[1];
    std::string out_address = argc[2];

    // Number of IBD (un-oscillated) and alpha-n events estimate (maily ratio matters)
    double N_IBD = std::atof(argc[3]);
    double IBD_err = std::atof(argc[4]);  // fractional error in N_IBD
    double N_alphaN = std::atof(argc[5]);
    double alphaN_err = std::atof(argc[6]);  // fractional error in N_alphaN
    double N_geoNu = std::atof(argc[7]);
    double geoNu_err = std::atof(argc[8]);  // fractional error in N_alphaN

    // 2d hist limit args
    double Dm21_lower = std::atof(argc[9]);
    double Dm21_upper = std::atof(argc[10]);
    double Theta12_lower = std::atof(argc[11]);  // degrees
    double Theta12_upper = std::atof(argc[12]);  // degrees
    unsigned int N_bins = std::atoi(argc[13]);

    // variable paramters limit args
    double Dm21_min = std::atof(argc[14]);
    double Dm21_max = std::atof(argc[15]);
    double Theta12_min = std::atof(argc[16]);  // degrees
    double Theta12_max = std::atof(argc[17]);  // degrees

    unsigned int start_idx_Dm21 = std::atoi(argc[18]);
    unsigned int start_idx_theta = std::atoi(argc[19]);

    unsigned int Dm21_nSteps = std::atoi(argc[20]);
    unsigned int Theta12_nSteps = std::atoi(argc[21]);

    bool verbose = std::stoi(argc[22]);

    if (verbose) {
        // Print input variables
        std::cout << "N_IBD = " << N_IBD << " ± " << IBD_err*100. << "%, " << "N_alphaN = " << N_alphaN << " ± " << alphaN_err*100. << "%, " << "N_geoNu = " << N_geoNu << " ± " << geoNu_err*100. << "%." << std::endl;
        std::cout << "fit params: Dm_21^2 € [" << Dm21_min << ", " << Dm21_max << "] (#" << Dm21_nSteps << "), " << "theta_12^2 € [" << Theta12_min << ", " << Theta12_max << "] (#" << Theta12_nSteps << ")." << std::endl;
        std::cout << "hist params: Dm_21^2 lims € {" << Dm21_lower << ", " << Dm21_upper << "}, " << "theta_12^2 lims € {" << Theta12_lower << ", " << Theta12_upper << "}, Nbins = " << N_bins << "." << std::endl;
    }

    // Read in file
    std::cout << "Reading in hists from file..." << std::endl;
    std::vector<TH1D*> reactor_hists;
    std::vector<TH1D*> alphaN_hists;
    std::vector<TH1D*> geoNu_hists;
    TH2D E_conv;
    read_hists_from_file(PDFs_address, reactor_hists, alphaN_hists, geoNu_hists, E_conv);

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

    // Create oscillation variables
    FitVar vDm21_2("deltamsqr21", fDmSqr21, 0, fDmSqr21, fDmSqr21);
    FitVar vDm32_2("deltamsqr32", fDmSqr32, 0, fDmSqr32, fDmSqr32);
    FitVar vS_12_2("sinsqrtheta12", fSSqrTheta12, 0, fSSqrTheta12, fSSqrTheta12);
    FitVar vS_13_2("sinsqrtheta13", fSSqrTheta13, 0, fSSqrTheta13, fSSqrTheta13);

    vDm21_2.HoldConstant(true);
    vDm32_2.HoldConstant(true);
    vS_12_2.HoldConstant(true);
    vS_13_2.HoldConstant(true);

    // Create geo-nu norm variable (allow to vary by ±3 sigma), and model
    double geoNuNorm_min = (1.0 - 3.0 * geoNu_err) * N_geoNu;
    if (geoNuNorm_min < 0.0) geoNuNorm_min = 0.0;
    FitVar geoNuNorm("geoNuNorm", N_geoNu, geoNu_err * N_geoNu, geoNuNorm_min, (1.0 + 3.0 * geoNu_err) * N_geoNu);

    // geoNu geoNuMod(&geoNuNorm, &vS_12_2, &vS_13_2, geoNu_hists.at(0));
    geoNu geoNuMod(&geoNuNorm, &vS_12_2, &vS_13_2, alphaN_hists.at(0));
    geoNuMod.hold_osc_params_const(true);  // This will also pre-compute the survival prob ahead of time

    // Create alpha-n norm variables (allow to vary by ±3 sigma), and model
    // Set the three norms to a third of total each for now (testing)
    double N_alphaN_third = N_alphaN/3.0;
    double alphaNNorm_min = (1.0 - 3.0 * alphaN_err) * N_alphaN_third;
    if (alphaNNorm_min < 0.0) alphaNNorm_min = 0.0;
    FitVar alphaNNorm_1("alphaNNorm_1", N_alphaN_third, alphaN_err * N_alphaN_third, alphaNNorm_min, (1.0 + 3.0 * alphaN_err) * N_alphaN_third);
    FitVar alphaNNorm_2("alphaNNorm_2", N_alphaN_third, alphaN_err * N_alphaN_third, alphaNNorm_min, (1.0 + 3.0 * alphaN_err) * N_alphaN_third);
    FitVar alphaNNorm_3("alphaNNorm_3", N_alphaN_third, alphaN_err * N_alphaN_third, alphaNNorm_min, (1.0 + 3.0 * alphaN_err) * N_alphaN_third);

    // alphaN alphaNMod(&alphaNNorm_1, &alphaNNorm_2, &alphaNNorm_3, alphaN_hists.at(0), alphaN_hists.at(1), alphaN_hists.at(2));
    alphaN alphaNMod(&alphaNNorm_1, &alphaNNorm_2, &alphaNNorm_3, alphaN_hists.at(0), alphaN_hists.at(0), alphaN_hists.at(0));

    // Create reactor norm variables (allow to vary by ±3 sigma), and model
    std::vector<std::string> reactor_names = {"BRUCE", "DARLINGTON", "PICKERING", "WORLD"};

    std::vector<double> hist_fracs;
    compute_unosc_fracs(reactor_hists, reactor_names, hist_fracs);
    std::vector<FitVar*> ReactorNorms;
    double reacNorm;
    double reacLowLim;
    for (unsigned int i = 0; i < hist_fracs.size(); ++i) {
        reacNorm = hist_fracs.at(i) * N_IBD;
        reacLowLim = (1.0 - 3.0 * IBD_err) * reacNorm;
        if (reacLowLim < 0.0) reacLowLim = 0.0;
        ReactorNorms.push_back(new FitVar(reactor_names.at(i).c_str(), reacNorm, IBD_err * reacNorm, reacLowLim, (1.0 + 3.0 * IBD_err) * reacNorm));
    }

    Reactor ReactorMod(&vDm21_2, &vDm32_2, &vS_12_2, &vS_13_2, reactor_hists, &E_conv, reactor_names, ReactorNorms, db);
    ReactorMod.hold_osc_params_const(true); // This will also compute oscillated reactor specs

    // Make fake dataset out of PDF hists (same function called to make PDFs)
    std::cout << "Creating fake dataset..." << std::endl;
    std::vector<TH1D*> osc_hists;
    ReactorMod.GetOscReactorHists(osc_hists);
    TH1D* data = (TH1D*)osc_hists.at(0)->Clone();
    data->Reset("ICES");
    std::cout << "data integral = " << data->Integral() << std::endl;

    // Add reactor spectra (already re-sclaed)
    for (unsigned int i = 0; i < osc_hists.size(); ++i) {
        std::cout << osc_hists.at(i)->GetName() << " integral = " << osc_hists.at(i)->Integral() << std::endl;
        data->Add(osc_hists.at(i));  // Add reactor events
    }
    std::cout << "data integral = " << data->Integral() << std::endl;

    // Add alpha-n spectra
    data->Add(alphaN_hists.at(0), alphaNNorm_1.val() / alphaN_hists.at(0)->Integral());
    // data->Add(alphaN_hists.at(1), alphaNNorm_2.val() / alphaN_hists.at(1)->Integral());
    // data->Add(alphaN_hists.at(2), alphaNNorm_3.val() / alphaN_hists.at(2)->Integral());
    data->Add(alphaN_hists.at(0), alphaNNorm_2.val() / alphaN_hists.at(0)->Integral());
    data->Add(alphaN_hists.at(0), alphaNNorm_3.val() / alphaN_hists.at(0)->Integral());

    // Add geo-nu spectrum
    // data->Add(geoNuMod.GetOscHist());

    std::cout << "data integral = " << data->Integral() << std::endl;

    // Package variables and models together, and pass to fitter (just use ReactorNorms vector)
    ReactorNorms.push_back(&alphaNNorm_1); ReactorNorms.push_back(&alphaNNorm_2); ReactorNorms.push_back(&alphaNNorm_3);
    // ReactorNorms.push_back(&geoNuNorm);
    ReactorNorms.push_back(&vDm21_2); ReactorNorms.push_back(&vDm32_2);
    ReactorNorms.push_back(&vS_12_2); ReactorNorms.push_back(&vS_13_2);

    std::vector<Model*> models = {&alphaNMod, &ReactorMod}; //, &geoNuMod};

    Fitter antinuFitter(data, ReactorNorms, models);

    // Make list of Dm_21^2 and s_12^2 values to iterate over
    std::vector<std::vector<double>> var_params = make_var_param_vals(Dm21_min, Dm21_max, Dm21_nSteps, Theta12_min, Theta12_max, Theta12_nSteps);
    std::vector<unsigned int> start_idx = {start_idx_theta, start_idx_Dm21};

    // Set up 2d log-likelihood histogram
    TH2D* minllHist = new TH2D("minllHist", "minimised likelihood values", N_bins, Theta12_lower, Theta12_upper, N_bins, Dm21_lower, Dm21_upper);
    minllHist->GetXaxis()->SetTitle("Theta_12");
    minllHist->GetYaxis()->SetTitle("Delta m_21^2");
    minllHist->GetYaxis()->SetTitleOffset(1);  // Move x-axis label further away from the axis (0 is default)
    minllHist->SetStats(0);  // Remove stats box

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;
    Fit_spectra(antinuFitter, vDm21_2, vS_12_2, geoNuMod, ReactorMod, var_params, minllHist, start_idx, verbose);

    // Write hist to file and close
    TFile *outroot = new TFile(out_address.c_str(), "RECREATE");
    minllHist->Write();
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
void Fit_spectra(Fitter& antinuFitter, FitVar& vDm21_2, FitVar& vS_12_2, Model& geoNuMod, Model& ReactorMod, const std::vector<std::vector<double>>& var_params, TH2D* minllHist, const std::vector<unsigned int>& start_idx, const bool verbose) {

    // Unpack
    std::vector<double> sinTheta12 = var_params.at(0);
    std::vector<double> Dm21 = var_params.at(1);

    /* ~~~~~  Compute best fit Log likelihood for each set of parameters (fraction of alpha-n vs reactor IBD events is fit in each loop)  ~~~~~ */
    std::cout << "Looping over oscillation parameters..." << std::endl;
    for (unsigned int i = 0; i < sinTheta12.size(); ++i) {
        // std::cout << "i = " << i << std::endl;
        vS_12_2.val() = sinTheta12.at(i);
        if (verbose) std::cout << "s_12^2 = " << vS_12_2.val() << std::endl;
        // geoNuMod.hold_osc_params_const(true);  // This will also pre-compute the survival prob ahead of time
        for (unsigned int j = 0; j < Dm21.size(); ++j) {
            // std::cout << "j = " << j << std::endl;
            vDm21_2.val() = Dm21.at(j);
            ReactorMod.hold_osc_params_const(true); // This will also compute oscillated reactor specs
            if (verbose) std::cout << "Dm_21^2 = " << vDm21_2.val() << std::endl;
            minllHist->SetBinContent(start_idx.at(0) + i + 1, start_idx.at(1) + j + 1, antinuFitter.fit_models());
        }
    }
}

std::vector<double>& compute_unosc_fracs(const std::vector<TH1D*>& Reactor_hists, const std::vector<std::string>& Reactor_names, std::vector<double>& hist_fracs) {

    double tot_int;
    for (unsigned int i = 0; i < Reactor_names.size(); ++i) {
        hist_fracs.push_back(0.0);
    }

    std::string origin_reactor;
    bool reactor_not_found;
    double integral;
    for (unsigned int i = 0; i < Reactor_hists.size(); ++i) {
        origin_reactor = SplitString(Reactor_hists.at(i)->GetName())[0];
        integral = Reactor_hists.at(i)->Integral();
        // See if reactor is in Reactor_names (except last element 'WORLD')
        reactor_not_found = true;
        for (unsigned int j = 0; j < (Reactor_names.size()-1); ++j) {
            if (origin_reactor == Reactor_names.at(j)) {
                hist_fracs.at(j) += integral;
                reactor_not_found = false;
                break;
            }
        }
        // Reactor was not in the list, so assign it to 'WORLD' (last element in Reactor_names)
        if (reactor_not_found) hist_fracs.at(hist_fracs.size()-1) += integral;

        tot_int += integral;
    }

    // Convert integrals to fractions of total integral
    for (unsigned int i = 0; i < Reactor_names.size(); ++i) {
        hist_fracs.at(i) /= tot_int;
    }

    return hist_fracs;
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
    double Dm21_step = (Dm21_max - Dm21_min) / (double)Dm21_nSteps;
    for (unsigned int n = 0; n < Dm21_nSteps; ++n) {
        Dm21.push_back(Dm21_min + (double)n * Dm21_step);
    }

    // s_12^2
    std::vector<double> sinTheta12;
    double Theta12_step = (Theta12_max - Theta12_min) / (double)(Theta12_nSteps - 1);
    for (unsigned int n = 0; n < Theta12_nSteps; ++n) {
        sinTheta12.push_back(pow(sin((Theta12_min + (double)n * Theta12_step)  * TMath::Pi() / 180.), 2));
    }

    // Package them together
    return {sinTheta12, Dm21};
}

/**
 * @brief Lists all histograms from root file into vectors of reactor, alphaN, geoNu and E_conv hists
 * 
 * @param file_address 
 * @param reactor_hists 
 * @param alphaN_hists
 * @param geoNu_hists
 * @param E_conv
 */
void read_hists_from_file(std::string file_address, std::vector<TH1D*>& reactor_hists, std::vector<TH1D*>& alphaN_hists, std::vector<TH1D*>& geoNu_hists, TH2D& E_conv) {

    TFile *fin = TFile::Open(file_address.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file." << std::endl;
        exit(1);
    }
    
    // Iterate through list of objects in root file
    TList* list = fin->GetListOfKeys() ;
    if (!list) {std::cout << "No keys found in file\n" << std::endl; exit(1);}
    TIter next(list);
    TObject* obj;
    TKey* key;
    std::string name;

    // Go through list of histograms and add them to temp list if they are included in name list
    while((key = (TKey*)next())){
        obj = key->ReadObj() ;
        if (obj->InheritsFrom(TH1::Class())) {
            // Check which histogram in file matches name
            name = obj->GetName();
            if (name == "geoNu") {
                geoNu_hists.push_back((TH1D*)obj);
            // } else if (name == "alphaN_1" || name == "alphaN_2" || name == "alphaN_3") {
            } else if (name == "alphaN") {
                alphaN_hists.push_back((TH1D*)obj);
            } else if (name == "E_conversion") {
                continue;
            } else {
                reactor_hists.push_back((TH1D*)obj);
            }
        }
    }

    // Get 2D hist
    E_conv = *(TH2D*)fin->Get("E_conversion");

    // Error handling
    if (reactor_hists.size() == 0) {
        std::cout << "ERROR: No reactor IBD histograms!" << std::endl;
        exit(1);
    }
    if (alphaN_hists.size() == 0) {
        std::cout << "ERROR: No alpha-n histograms!" << std::endl;
        exit(1);
    }
    // if (geoNu_hists.size() == 0) {
    //     std::cout << "ERROR: No geo-nu histograms!" << std::endl;
    //     exit(1);
    // }
}
