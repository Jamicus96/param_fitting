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
#include "fitting_utils.hpp"

#include <RooRealVar.h>
#include <RooDataHist.h>
// #include "RooDataSet.h"
#include <RooClassFactory.h>
#include <RooTFnBinding.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooConstVar.h>

using namespace RooFit;

void Fit_spectra(reactorINFO& spectrum, TH1D& alphaN_hist, const double N_alphaN, const double alphaN_err, TH1D* data, const std::vector<std::vector<double>>& var_params, TH2D* minllHist, const bool verbose);
std::vector<std::vector<double>> make_var_param_vals(const double Dm21_min, const double Dm21_max, const unsigned int Dm21_nSteps, const double Theta12_min, const double Theta12_max, const unsigned int Theta12_nSteps);
double ML_fit(reactorINFO& spectrum, RooHistPdf& alphaN_PDF, RooRealVar& norm_alphaN, RooGaussian& constraint_alphaN, RooDataHist& dataHist, const RooRealVar& E, const RooRealVar& reactor_frac, const RooRealVar& alphaN_frac);
void read_hists_from_file(std::string file_address, std::vector<TH1D*>& reactor_hists, TH1D& alphaN_hist);


int main(int argv, char** argc) {
    // file args
    std::string PDFs_address = argc[1];
    std::string out_address = argc[2];

    // Number of IBD (un-oscillated) and alpha-n events estimate (maily ratio matters)
    double N_IBD = atof(argc[3]);
    double IBD_err = atof(argc[4]);  // fractional error in N_IBD
    double N_alphaN = atof(argc[5]);
    double alphaN_err = atof(argc[6]);  // fractional error in N_alphaN

    // 2d hist limit args
    double Dm21_lower = atof(argc[7]);
    double Dm21_upper = atof(argc[8]);
    double Theta12_lower = atof(argc[9]);  // degrees
    double Theta12_upper = atof(argc[10]);  // degrees
    unsigned int N_bins = atoi(argc[11]);

    // variable paramters limit args
    double Dm21_min = atof(argc[12]);
    double Dm21_max = atof(argc[13]);
    double Theta12_min = atof(argc[14]);  // degrees
    double Theta12_max = atof(argc[15]);  // degrees
    unsigned int Dm21_nSteps = atoi(argc[16]);
    unsigned int Theta12_nSteps = atoi(argc[17]);

    bool verbose = std::stoi(argc[18]);

    if (verbose) {
        // Print input variables
        std::cout << "N_IBD = " << N_IBD << " ± " << IBD_err*100. << "%, " << "N_alphaN = " << N_alphaN << " ± " << alphaN_err*100. << "%." << std::endl;
        std::cout << "fit params: Dm_21^2 € [" << Dm21_min << ", " << Dm21_max << "] (#" << Dm21_nSteps << "), " << "theta_12^2 € [" << Theta12_min << ", " << Theta12_max << "] (#" << Theta12_nSteps << ")." << std::endl;
        std::cout << "hist params: Dm_21^2 lims € {" << Dm21_lower << ", " << Dm21_upper << "}, " << "theta_12^2 lims € {" << Theta12_lower << ", " << Theta12_upper << "}, Nbins = " << N_bins << "." << std::endl;
    }

    // Read in file
    std::cout << "Reading in hists from file..." << std::endl;
    std::vector<TH1D*> reactor_hists;
    TH1D alphaN_hist;
    read_hists_from_file(PDFs_address, reactor_hists, alphaN_hist);

    // Get baselines
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DB* db = RAT::DB::Get();
    db->LoadDefaults();
    std::cout << "Getting baselines..." << std::endl;

    // Get oscillation constants
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    const double fDmSqr32 = linkdb->GetD("deltamsqr32");
    const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
    const double fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    // Create PDF fitting object
    reactorINFO spectrum(reactor_hists, N_IBD, IBD_err, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);
    spectrum.compute_baselines(db);

    // Make fake dataset out of PDF hists (same function called to make PDFs)
    std::cout << "Creating fake dataset..." << std::endl;
    spectrum.compute_osc_reactor_spec();
    std::vector<TH1D*> osc_hists = spectrum.Get_osc_reactor_specs();
    std::vector<double> norms = spectrum.Get_osc_reactor_norms();
    TH1D* data = (TH1D*)osc_hists.at(0)->Clone();

    // Rescale each hist to have them integrate to the correct total number of events given by N_IBD (+oscillation) and N_alphaN
    for (unsigned int i; i < osc_hists.size(); ++i) {
        data->Add(osc_hists.at(i), norms.at(i) / osc_hists.at(i)->Integral());  // Add reactor events
    }
    data->Add(&alphaN_hist, N_alphaN / alphaN_hist.Integral());  // Add alpha-n events

    // Make list of Dm_21^2 and s_12^2 values to iterate over
    std::vector<std::vector<double>> var_params = make_var_param_vals(Dm21_min, Dm21_max, Dm21_nSteps, Theta12_min, Theta12_max, Theta12_nSteps);

    // Set up 2d log-likelihood histogram
    TH2D* minllHist = new TH2D("minllHist", "minimised likelihood values", N_bins, Theta12_lower, Theta12_upper, N_bins, Dm21_lower, Dm21_upper);
    minllHist->GetXaxis()->SetTitle("Theta_12");
    minllHist->GetYaxis()->SetTitle("Delta m_21^2");
    minllHist->GetYaxis()->SetTitleOffset(1);  // Move x-axis label further away from the axis (0 is default)
    minllHist->SetStats(0);  // Remove stats box

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;
    Fit_spectra(spectrum, alphaN_hist, N_alphaN, alphaN_err, data, var_params, minllHist, verbose);

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
 * @param spectrum  reactorINFO object that contains all the relevent reactor
 * @param alphaN_hist  alphaN pdf
 * @param N_alphaN  number of alphaN events
 * @param alphaN_err  fractional error of N_alphaN
 * @param data  data to fit model to 
 * @param var_params  Dm21^2 and s12^2 values to iterate over
 * @param minllHist  2D histogram that min log-likelihood values are dumped in
 */
void Fit_spectra(reactorINFO& spectrum, TH1D& alphaN_hist, const double N_alphaN, const double alphaN_err, TH1D* data, const std::vector<std::vector<double>>& var_params, TH2D* minllHist, const bool verbose) {

    // Unpack
    std::vector<double> sinTheta12 = var_params.at(0);
    std::vector<double> Dm21 = var_params.at(1);

    // Declare observables
    RooRealVar E("E", "energy", 0.9, 8);
    RooRealVar reactor_frac("reactor_frac", "reactorIBD fraction", 0.4, 0.8);
    RooRealVar alphaN_frac("alphaN_frac", "alpha-n fraction", 0.2, 0.6);

    // Make data hist (data that must be fit)
    RooDataHist dataHist("dataHist", "data hist", E, data);

    // Create alphaN PDF and norm + norm constraint
    RooDataHist* tempData_alphaN = new RooDataHist("tempData", "temporary data", E, &alphaN_hist);
    RooHistPdf* alphaN_PDF = new RooHistPdf("alphaN_PDF", "alphaN PDF", E, *tempData_alphaN);

    double lower_bound = (1.0 - 3.0 * alphaN_err) * N_alphaN;
    if (lower_bound < 0.) lower_bound = 0.;
    RooRealVar norm_alphaN("alpa-n Norm", "temporary norm", N_alphaN, lower_bound, (1.0 + 3.0 * alphaN_err) * N_alphaN); // Allow variable to vary between ±3 sigmas
    RooGaussian constraint_alphaN("alpha-n Norm Const", "temp constraint", norm_alphaN, RooConst(N_alphaN), RooConst(alphaN_err * N_alphaN));  // var, mean, sigma

    // Compute best fit Log likelihood for each set of paramters (fraction of alpha-n vs reactor IBD events is fit in each loop)
    double MLL;
    std::cout << "Looping over oscillation parameters..." << std::endl;
    for (unsigned int i = 0; i < sinTheta12.size(); ++i) {
        // std::cout << "i = " << i << std::endl;
        spectrum.s12_2() = sinTheta12.at(i);
        if (verbose) std::cout << "s_12^2 = " << spectrum.s12_2() << std::endl;
        for (unsigned int j = 0; j < Dm21.size(); ++j) {
            // std::cout << "j = " << j << std::endl;
            spectrum.Dm21_2() = Dm21.at(j);
            if (verbose) std::cout << "Dm_21^2 = " << spectrum.Dm21_2() << std::endl;
            minllHist->SetBinContent(i+1, j+1, ML_fit(spectrum, *alphaN_PDF, norm_alphaN, constraint_alphaN, dataHist, E, reactor_frac, alphaN_frac));
        }
    }

    // clean up
    delete(tempData_alphaN);
    delete(alphaN_PDF);
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

double ML_fit(reactorINFO& spectrum, RooHistPdf& alphaN_PDF, RooRealVar& norm_alphaN, RooGaussian& constraint_alphaN, RooDataHist& dataHist, const RooRealVar& E, const RooRealVar& reactor_frac, const RooRealVar& alphaN_frac) {

    // Compute oscillated reactor spectra
    spectrum.compute_osc_reactor_spec();

    // Declare argument lists for reactor PDFs
    RooArgList components = RooArgList();
    RooArgList coeffs = RooArgList();
    RooArgList constraints = RooArgList();

    RooDataHist* tempData;
    RooHistPdf* tempPDF;
    std::string name;
    std::string norm_name;
    std::string const_name;
    double lower_bound;
    for (unsigned int i = 0; i < spectrum.Get_osc_reactor_specs().size(); ++i) {
        // Make argument names
        name = spectrum.Get_osc_reactor_specs().at(i)->GetName();
        norm_name = name + "_norm";
        const_name = norm_name + "_const";

        // Add PDF with associated norm and norm constraint
        tempData = new RooDataHist("tempData", "temporary data", E, spectrum.Get_osc_reactor_specs().at(i));
        tempPDF = new RooHistPdf(name.c_str(), "temporary PDF", E, *tempData);

        lower_bound = (1.0 - 3.0 * spectrum.IBDs_err()) * spectrum.Get_osc_reactor_norms().at(i);
        if (lower_bound < 0.) lower_bound = 0.;
        RooRealVar normTemp(norm_name.c_str(), "temporary norm", spectrum.Get_osc_reactor_norms().at(i), lower_bound, (1.0 + 3.0 * spectrum.IBDs_err()) * spectrum.Get_osc_reactor_norms().at(i)); // Allow variable to vary between ±3 sigmas
        RooGaussian tempConstraint(const_name.c_str(), "temp constraint", normTemp, RooConst(spectrum.Get_osc_reactor_norms().at(i)), RooConst(spectrum.IBDs_err() * spectrum.Get_osc_reactor_norms().at(i)));  // var, mean, sigma

        components.addClone(*tempPDF);
        coeffs.addClone(normTemp);
        constraints.addClone(tempConstraint);
    }

    // Do the same for alpha-n
    components.addClone(alphaN_PDF);
    coeffs.addClone(norm_alphaN);
    constraints.addClone(constraint_alphaN);

    // Add PDFs together
    RooAddPdf model("model", "r+a", components,  coeffs);

    // Fit to data
    double result = model.fitTo(dataHist, Extended(true), PrintLevel(-1), SumW2Error(kFALSE), Save())->minNll();

    // clean up
    delete(tempData);
    delete(tempPDF);

    //return result
    return result;
}


/**
 * @brief Lists all TH1D histograms from root file into a vector of reactor hists and an alphaN hist
 * 
 * @param file_address 
 * @param reactor_hists 
 * @param alphaN_hist 
 */
void read_hists_from_file(std::string file_address, std::vector<TH1D*>& reactor_hists, TH1D& alphaN_hist) {

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
    bool alphaN_hist_found, reactor_hist_found = false;
    while((key = (TKey*)next())){
        obj = key->ReadObj() ;
        if(obj->InheritsFrom(TH1::Class())){
            // Check which histogram in file matches name
            name = obj->GetName();
            if (name == "alphaN") {
                alphaN_hist = *(TH1D*)obj;
                alphaN_hist_found = true;
            } else {
                reactor_hists.push_back((TH1D*)obj);
                reactor_hist_found = true;
            }
        }
    }

    // Error handling
    if (!reactor_hist_found) {
        std::cout << "ERROR: No reactor IBD histograms!" << std::endl;
        exit(1);
    }
    if (!alphaN_hist_found) {
        std::cout << "ERROR: No alpha-n histograms!" << std::endl;
        exit(1);
    }
}
