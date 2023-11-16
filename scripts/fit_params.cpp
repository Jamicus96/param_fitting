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

void Fit_spectra(reactorINFO& spectrum, std::vector<TH1D*>& alphaN_hists, const double N_alphaN, const double alphaN_err, std::vector<TH1D*>& geoNu_hists, const double  N_geoNu, const double  geoNu_err, TH1D* data, const std::vector<std::vector<double>>& var_params, TH2D* minllHist, const std::vector<unsigned int>& start_idx, const bool verbose);
std::vector<std::vector<double>> make_var_param_vals(const double Dm21_min, const double Dm21_max, const unsigned int Dm21_nSteps, const double Theta12_min, const double Theta12_max, const unsigned int Theta12_nSteps);
RooDataHist* Make_PDFs_norms_constraints(std::vector<TH1D*>& hists, const RooRealVar& E, std::vector<double>& norms, double frac_err, std::vector<RooHistPdf*>& RooPDFs, std::vector<RooRealVar>& RooNorms, std::vector<RooGaussian>& RooConstraints);
double ML_fit(reactorINFO& spectrum, std::vector<RooHistPdf*>& bckgnd_PDFs, std::vector<RooRealVar>& bckgnd_norms, std::vector<RooGaussian>& bckgnd_constraints, RooDataHist& dataHist, const RooRealVar& E);
void read_hists_from_file(std::string file_address, std::vector<TH1D*>& reactor_hists, std::vector<TH1D*>& alphaN_hists, std::vector<TH1D*>& geoNu_hists, TH2D& E_conv);


int main(int argv, char** argc) {
    // file args
    std::string PDFs_address = argc[1];
    std::string out_address = argc[2];

    // Number of IBD (un-oscillated) and alpha-n events estimate (maily ratio matters)
    double N_IBD = atof(argc[3]);
    double IBD_err = atof(argc[4]);  // fractional error in N_IBD
    double N_alphaN = atof(argc[5]);
    double alphaN_err = atof(argc[6]);  // fractional error in N_alphaN
    double N_geoNu = atof(argc[7]);
    double geoNu_err = atof(argc[8]);  // fractional error in N_alphaN

    // 2d hist limit args
    double Dm21_lower = atof(argc[9]);
    double Dm21_upper = atof(argc[10]);
    double Theta12_lower = atof(argc[11]);  // degrees
    double Theta12_upper = atof(argc[12]);  // degrees
    unsigned int N_bins = atoi(argc[13]);

    // variable paramters limit args
    double Dm21_min = atof(argc[14]);
    double Dm21_max = atof(argc[15]);
    double Theta12_min = atof(argc[16]);  // degrees
    double Theta12_max = atof(argc[17]);  // degrees

    unsigned int start_idx_Dm21 = atoi(argc[18]);
    unsigned int start_idx_theta = atoi(argc[19]);

    unsigned int Dm21_nSteps = atoi(argc[20]);
    unsigned int Theta12_nSteps = atoi(argc[21]);

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
    reactorINFO spectrum(reactor_hists, N_IBD, IBD_err, E_conv, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);
    spectrum.compute_baselines(db);

    // Make fake dataset out of PDF hists (same function called to make PDFs)
    std::cout << "Creating fake dataset..." << std::endl;
    spectrum.compute_osc_reactor_spec();
    std::vector<TH1D*> osc_hists = spectrum.Get_osc_reactor_specs();
    std::vector<double> norms = spectrum.Get_osc_reactor_norms();
    TH1D* data = (TH1D*)osc_hists.at(0)->Clone();
    data->Reset("ICES");

    // Rescale each hist to have them integrate to the correct total number of events given by N_IBD (+oscillation) and N_alphaN
    for (unsigned int i; i < osc_hists.size(); ++i) {
        data->Add(osc_hists.at(i), norms.at(i) / osc_hists.at(i)->Integral());  // Add reactor events
    }
    data->Add(alphaN_hists.at(0), N_alphaN / alphaN_hists.at(0)->Integral());  // Add alpha-n events

    std::cout << "data integral = " << data->Integral() << std::endl;

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
    Fit_spectra(spectrum, alphaN_hists, N_alphaN, alphaN_err, geoNu_hists, N_geoNu, geoNu_err, data, var_params, minllHist, start_idx, verbose);

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
void Fit_spectra(reactorINFO& spectrum, std::vector<TH1D*>& alphaN_hists, const double N_alphaN, const double alphaN_err, std::vector<TH1D*>& geoNu_hists, const double  N_geoNu, const double  geoNu_err, TH1D* data, const std::vector<std::vector<double>>& var_params, TH2D* minllHist, const std::vector<unsigned int>& start_idx, const bool verbose) {

    // Unpack
    std::vector<double> sinTheta12 = var_params.at(0);
    std::vector<double> Dm21 = var_params.at(1);

    // Declare observable
    RooRealVar E("E", "energy", 0.9, 8);

    // Make data hist (data that must be fit)
    RooDataHist dataHist("dataHist", "data hist", E, data);


    /* ~~~~~~~~~~  Create list of background PDFs, norms and constraints (alpha-n and goe-nu)  ~~~~~~~~~~ */
    std::vector<RooHistPdf*> bckgnd_PDFs;
    std::vector<RooRealVar> bckgnd_norms;
    std::vector<RooGaussian> bckgnd_constraints;

    // alpha-n
    double tot_alphaN_int = 0.0;
    for (unsigned int i = 0; i < alphaN_hists.size(); ++i) {
        tot_alphaN_int += alphaN_hists.at(i)->Integral();
    }
    std::vector<double> alphaN_norms;
    for (unsigned int i = 0; i < alphaN_hists.size(); ++i) {
        alphaN_norms.push_back(N_alphaN * alphaN_hists.at(i)->Integral() / tot_alphaN_int);  // get norm of each PDF from total number of alphaNs
    }
    RooDataHist* tempData = Make_PDFs_norms_constraints(alphaN_hists, E, alphaN_norms, alphaN_err, bckgnd_PDFs, bckgnd_norms, bckgnd_constraints);
    
    // Geo-nu
    tempData = new RooDataHist("tempData", "temporary data", E, geoNu_hists.at(0));
    bckgnd_PDFs.push_back(new RooHistPdf("geoNu", "temporary PDF", E, *tempData));
    // Add placeholder norms and constraints (will update them inside loop, below)
    unsigned int geoNu_idx = bckgnd_PDFs.size() - 1;
    bckgnd_norms.push_back(RooRealVar("geoNu", "temporary norm", N_geoNu, 0.99 * N_geoNu, 1.01 * N_geoNu));
    bckgnd_constraints.push_back(RooGaussian("geoNu", "temp constraint", bckgnd_norms.at(geoNu_idx), RooConst(N_geoNu), RooConst(0.01 * N_geoNu)));


    /* ~~~~~  Compute best fit Log likelihood for each set of parameters (fraction of alpha-n vs reactor IBD events is fit in each loop)  ~~~~~ */
    double MLL;
    double new_N_geoNu;
    double lower_bound;
    std::cout << "Looping over oscillation parameters..." << std::endl;
    for (unsigned int i = 0; i < sinTheta12.size(); ++i) {
        // std::cout << "i = " << i << std::endl;
        spectrum.s12_2() = sinTheta12.at(i);
        if (verbose) std::cout << "s_12^2 = " << spectrum.s12_2() << std::endl;

        // Compute geo-nu survival prob, and update re-scaled norm and constraint in background lists
        new_N_geoNu = N_geoNu * spectrum.geoNu_survival_prob();
        lower_bound = (1.0 - 3.0 * geoNu_err) * new_N_geoNu;
        if (lower_bound < 0.) lower_bound = 0.;
        bckgnd_norms.at(geoNu_idx) = RooRealVar("geoNu", "temporary norm", new_N_geoNu, lower_bound, (1.0 + 3.0 * geoNu_err) * new_N_geoNu); // Allow variable to vary between ±3 sigmas
        bckgnd_constraints.at(geoNu_idx) = RooGaussian("geoNu", "temp constraint", bckgnd_norms.at(geoNu_idx), RooConst(new_N_geoNu), RooConst(geoNu_err * new_N_geoNu));  // var, mean, sigma

        for (unsigned int j = 0; j < Dm21.size(); ++j) {
            // std::cout << "j = " << j << std::endl;
            spectrum.Dm21_2() = Dm21.at(j);
            if (verbose) std::cout << "Dm_21^2 = " << spectrum.Dm21_2() << std::endl;
            minllHist->SetBinContent(start_idx.at(0) + i + 1, start_idx.at(1) + j + 1, ML_fit(spectrum, bckgnd_PDFs, bckgnd_norms, bckgnd_constraints, dataHist, E));
        }
    }

    // clean up
    delete(tempData);
    for (auto p : bckgnd_PDFs) {delete p;} 
    bckgnd_PDFs.clear();
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
 * @brief Package histogram PDFs and total norm predictions (+errs) into the appropriate RooFit objects.
 * Outputs RooDataHist* object, to be deleted after fitting to avoid memory leaks.
 * 
 * @param hists 
 * @param norms 
 * @param frac_err 
 * @param bckgnd_PDFs 
 * @param bckgnd_norms 
 * @param bckgnd_constraints 
 * @return RooDataHist* 
 */
RooDataHist* Make_PDFs_norms_constraints(std::vector<TH1D*>& hists, const RooRealVar& E, std::vector<double>& norms, double frac_err, std::vector<RooHistPdf*>& RooPDFs, std::vector<RooRealVar>& RooNorms, std::vector<RooGaussian>& RooConstraints) {
    RooDataHist* tempData;
    double lower_bound;
    std::string name;
    std::string norm_name;
    std::string const_name;
    for (unsigned int i = 0; i < hists.size(); ++i) {
        // Make argument names
        name = hists.at(i)->GetName();
        norm_name = name + "_norm";
        const_name = norm_name + "_const";

        // Add PDF with associated norm and norm constraint
        tempData = new RooDataHist("tempData", "temporary data", E, hists.at(i));
        RooPDFs.push_back(new RooHistPdf(name.c_str(), "temporary PDF", E, *tempData));

        lower_bound = (1.0 - 3.0 * frac_err) * norms.at(i);
        if (lower_bound < 0.) lower_bound = 0.;
        RooNorms.push_back(RooRealVar(norm_name.c_str(), "temporary norm", norms.at(i), lower_bound, (1.0 + 3.0 * frac_err) * norms.at(i))); // Allow variable to vary between ±3 sigmas
        RooConstraints.push_back(RooGaussian(const_name.c_str(), "temp constraint", RooNorms.at(RooNorms.size()-1), RooConst(norms.at(i)), RooConst(frac_err * norms.at(i))));  // var, mean, sigma
    }
    return tempData;
}

double ML_fit(reactorINFO& spectrum, std::vector<RooHistPdf*>& bckgnd_PDFs, std::vector<RooRealVar>& bckgnd_norms, std::vector<RooGaussian>& bckgnd_constraints, RooDataHist& dataHist, const RooRealVar& E) {

    // Compute oscillated reactor spectra
    spectrum.compute_osc_reactor_spec();

    // Declare argument lists for reactor and background PDFs norms and constraints
    RooArgList components = RooArgList();
    RooArgList coeffs = RooArgList();
    RooArgList constraints = RooArgList();

    // Add reactor PDFs
    std::vector<RooHistPdf*> reactor_PDFs;
    std::vector<RooRealVar> reactor_norms;
    std::vector<RooGaussian> reactor_constraints;
    RooDataHist* tempData = Make_PDFs_norms_constraints(spectrum.Get_osc_reactor_specs(), E, spectrum.Get_osc_reactor_norms(), spectrum.IBDs_err(), reactor_PDFs, reactor_norms, reactor_constraints);
    for (unsigned int i = 0; i < reactor_PDFs.size(); ++i) {
        components.addClone(*(reactor_PDFs.at(i)));
        coeffs.addClone(reactor_norms.at(i));
        constraints.addClone(reactor_constraints.at(i));
    }

    // Do the same for backgrounds
    for (unsigned int i = 0; i < bckgnd_PDFs.size(); ++i) {
        components.addClone(*(bckgnd_PDFs.at(i)));
        coeffs.addClone(bckgnd_norms.at(i));
        constraints.addClone(bckgnd_constraints.at(i));
    }

    // Add PDFs together
    RooAddPdf model("model", "r+a", components,  coeffs);

    // Fit to data
    double result = model.fitTo(dataHist, Extended(true), PrintLevel(-1), SumW2Error(kFALSE), Save())->minNll();

    // clean up
    delete(tempData);
    for (auto p : reactor_PDFs) {delete p;} 
    reactor_PDFs.clear();

    //return result
    return result;
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
            } else if (name == "alphaN_1" || name == "alphaN_2" || name == "alphaN_3") {
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
    if (geoNu_hists.size() == 0) {
        std::cout << "ERROR: No geo-nu histograms!" << std::endl;
        exit(1);
    }
}
