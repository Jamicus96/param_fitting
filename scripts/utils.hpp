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
#include <map>

#include "fitter.hpp"
#include "fitVars.hpp"
#include "E_systematics.hpp"
#include "model.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"
#include "model_Reactor.hpp"


void create_fitter(std::string PDFs_address, double Dm21_2, double Dm32_2, double s_12_2, double s_13_2, RAT::DB* db);
void compute_hist_fracs(const std::vector<TH1D*>& hists, const std::vector<std::string>& hist_names, std::vector<double>& hist_fracs, std::vector<unsigned int>& hist_idx);
void compute_reac_unosc_fracs(const std::vector<TH1D*>& Reactor_hists, const std::vector<std::string>& Reactor_names, std::vector<double>& hist_fracs);
void read_hists_from_file(std::string file_address, std::vector<TH1D*>& reactor_hists, std::vector<TH1D*>& alphaN_hists, std::vector<TH1D*>& geoNu_hists, TH2D& E_conv);


/* ~~~~~~~~ CONSTRAINED PARAMETERS ~~~~~~~~ */

double N_IBD = 122.0;           // Total number of expected reactor IBDs (un-oscillated)
double IBD_err_indiv = 0.032;   // fractional error in N_IBD for each individual reactor PDF
double IBD_err_tot = 0.03;      // fractional error in N_IBD for total reactor IBDs

double N_alphaN = 50.0;         // Total number of expected alpha-n
double alphaN_err_GS = 0.3;     // fractional error in N_alphaN for ground state neutrons (PR + C12)
double alphaN_err_ES = 1.0;     // fractional error in N_alphaN for excited state neutrons (O16)

double N_geoNu = 5.0;           // Total number of expected geo-nu IBDs (un-oscillated)
double geoNu_err = 0.2;         // fractional error in N_geoNu for individual Th and U spectra

double linScale_err = 0.011;    // Error in linear scaling (scaling = 1) (not fractional)
double kB = 0.074;              // Birk's constant for betas
double kB_err = 0.004;          // Error in kB (not fractional)
double sigPerSqrtE = 0.042;     // smearing sigma = sigPerSqrtE * sqrt(E)

double linScale_err_P = 0.011;  // Error linear scaling for proton recoils (scaling = 1) (not fractional)
double kB_P = 0.078;            // Birk's constant for protons
double kB_err_P = 0.004;        // Error in kB_P for proton recoils (not fractional)
// Proton recoil uses the same smearing as the rest

/* ~~~~~~~~ SET UP FITTER ~~~~~~~~ */

/**
 * @brief Create a fitter object.
 * models = {alphaNMod, ReactorMod, geoNuMod}
 * Vars = {BruceMorm, DarlingtonNorm, PickeringNorm, WorldNorm, totNorm,
 *         alphaNNorm_GS, ReactorNorms, geoNuNormTh, geoNuNormU, Dm21_2
 *         Dm32_2, S_12_2, S_13_2, linScale, kBp, linScale_P, kBp_P, smearSigmaSqrtE}
 * 
 * @param PDFs_address 
 * @param Dm21_2 
 * @param Dm32_2 
 * @param s_12_2 
 * @param s_13_2 
 * @param db 
 */
void create_fitter(std::string PDFs_address, double Dm21_2, double Dm32_2, double s_12_2, double s_13_2, RAT::DB* db) {
    // Read in file
    std::cout << "Reading in hists from file..." << std::endl;
    std::vector<TH1D*> reactor_hists;
    std::vector<TH1D*> alphaN_hists;
    std::vector<TH1D*> geoNu_hists;
    TH2D E_conv;
    read_hists_from_file(PDFs_address, reactor_hists, alphaN_hists, geoNu_hists, E_conv);

    /* ~~~~~~~~ INITIALISE ~~~~~~~~ */

    Fitter antinuFitter = Fitter();

    /* ~~~~~~~~ OSCILLATION CONSTANTS ~~~~~~~~ */

    // Add oscillation variables
    antinuFitter.GetVars().AddVar("deltamsqr21", Dm21_2, 0, Dm21_2, Dm21_2, true);
    antinuFitter.GetVars().AddVar("deltamsqr32", Dm32_2, 0, Dm32_2, Dm32_2, true);
    antinuFitter.GetVars().AddVar("sinsqrtheta12", s_12_2, 0, s_12_2, s_12_2, true);
    antinuFitter.GetVars().AddVar("sinsqrtheta13", s_13_2, 0, s_13_2, s_13_2, true);


    /* ~~~~~~~~ ENERGY SYSTEMATICS ~~~~~~~~ */

    // Add beta variables
    double linScale_min = 1.0 - 3.0 * linScale_err;
    double kBp_min = kB - 3.0 * kB_err;
    if (linScale_min < 0.0) linScale_min = 0.0;
    if (kBp_min < 0.0) kBp_min = 0.0;
    antinuFitter.GetVars().AddVar("linScale", 1.0, linScale_err, linScale_min, 1.0 + 3.0 * linScale_err);
    antinuFitter.GetVars().AddVar("kBp", kB, kB_err, kBp_min, kB + 3.0 * kB_err);

    // Add proton variables (linear scaling is the same, but independent)
    double kBp_P_min = kB_P - 3.0 * kB_err_P;  // same error, but different Birk's constant
    if (kBp_P_min < 0.0) kBp_P_min = 0.0;
    antinuFitter.GetVars().AddVar("linScale_P", 1.0, linScale_err, linScale_min, 1.0 + 3.0 * linScale_err);
    antinuFitter.GetVars().AddVar("kBp_P", kB_P, kB_err_P, kBp_P_min, kB_P + 3.0 * kB_err_P);

    // Add smearing variable (same for all, allow to be negative, absolute value taken in calculations)
    antinuFitter.GetVars().AddVar("sigPerSqrtE", 0.0, sigPerSqrtE, -3.0 * sigPerSqrtE, 3.0 * sigPerSqrtE);

    // Add energy systamtics objects (beta and proton), and link them to variables above
    antinuFitter.GetEsysts().AddEsys("EsysBeta", kB, "linScale", "kBp", "sigPerSqrtE");
    antinuFitter.GetEsysts().AddEsys("EsysProton", kB, "linScale_P", "kBp_P", "sigPerSqrtE");


    /* ~~~~~~~~ GEO-NU ~~~~~~~~ */

    // Get fraction of (unoscillated) geo-nu flux coming from each histogram (Th and U, in order)
    std::vector<std::string> geoNu_names = {"geoNu_Th", "geoNu_U"};
    std::vector<double> geoNu_hist_fracs;
    std::vector<unsigned int> geoNu_hist_idx;
    compute_hist_fracs(geoNu_hists, geoNu_names, geoNu_hist_fracs, geoNu_hist_idx);

    // Create geo-nu norm variables (allow to vary by ±3 sigma), and model
    double geoNuNorm_frac_min = 1.0 - 3.0 * geoNu_err;
    if (geoNuNorm_frac_min < 0.0) geoNuNorm_frac_min = 0.0;
    
    antinuFitter.GetVars().AddVar("geoNuNorm_Th", N_geoNu * geoNu_hist_fracs.at(0), geoNu_err * N_geoNu * geoNu_hist_fracs.at(0),
                                  geoNuNorm_frac_min * N_geoNu * geoNu_hist_fracs.at(0), (1.0 + 3.0 * geoNu_err) * N_geoNu * geoNu_hist_fracs.at(0));
    antinuFitter.GetVars().AddVar("geoNuNorm_U", N_geoNu * geoNu_hist_fracs.at(1), geoNu_err * N_geoNu * geoNu_hist_fracs.at(1),
                                  geoNuNorm_frac_min * N_geoNu * geoNu_hist_fracs.at(1), (1.0 + 3.0 * geoNu_err) * N_geoNu * geoNu_hist_fracs.at(1));

    // Add geo-nu model, linking it to approproate variables and E-systematics defined above
    antinuFitter.AddGeoNuMod("geoNuNorm_Th", "geoNuNorm_U", "sinsqrtheta12", "sinsqrtheta13", "EsysBeta", geoNu_hists.at(geoNu_hist_idx.at(0)), geoNu_hists.at(geoNu_hist_idx.at(1)));
    antinuFitter.GetModels().at(0)->hold_osc_params_const(true);  // This will also pre-compute the survival prob ahead of time


    /* ~~~~~~~~ ALPHA-N ~~~~~~~~ */

    // Get fraction of alpha-n flux coming from each histogram (in order: proton recoil, C12 scatter, O16 de-excitation)
    std::vector<std::string> alphaN_names = {"alphaN_PR", "alphaN_C12", "alphaN_O16"};
    std::vector<double> alphaN_hist_fracs;
    std::vector<unsigned int> alphaN_hist_idx;
    compute_hist_fracs(alphaN_hists, alphaN_names, alphaN_hist_fracs, alphaN_hist_idx);
    double GS_frac = alphaN_hist_fracs.at(0) + alphaN_hist_fracs.at(1);
    double ES_frac = alphaN_hist_fracs.at(2);

    // Create alpha-n norm variables (allow to vary by ±3 sigma), and model
    double GS_Norm_min = (1.0 - 3.0 * alphaN_err_GS) * N_alphaN * GS_frac;
    if (GS_Norm_min < 0.0) GS_Norm_min = 0.0;
    double ES_Norm_min = (1.0 - 3.0 * alphaN_err_ES) * N_alphaN * ES_frac;
    if (ES_Norm_min < 0.0) ES_Norm_min = 0.0;

    antinuFitter.GetVars().AddVar("alphaNNorm_GS", N_alphaN * GS_frac, alphaN_err_GS * N_alphaN * GS_frac, GS_Norm_min, (1.0 + 3.0 * alphaN_err_GS) * N_alphaN * GS_frac);
    antinuFitter.GetVars().AddVar("alphaNNorm_ES", N_alphaN * ES_frac, alphaN_err_ES * N_alphaN * ES_frac, ES_Norm_min, (1.0 + 3.0 * alphaN_err_ES) * N_alphaN * ES_frac);

    for (unsigned int i = 0; i < 3; ++i) {
        std::cout << "alphaN_hists.at(" << i << ")->GetName() = " << alphaN_hists.at(i)->GetName() << std::endl;
        std::cout << "alphaN_hists.at(" << i << ")->Integral()  = " << alphaN_hists.at(i)->Integral() << std::endl;
    }

    // Add alpha-n model, linking it to approproate variables and E-systematics defined above
    antinuFitter.AddAlphaNMod("alphaNNorm_GS", "alphaNNorm_ES", "EsysBeta", "EsysProton", alphaN_hists.at(alphaN_hist_idx.at(0)), alphaN_hists.at(alphaN_hist_idx.at(1)), alphaN_hists.at(alphaN_hist_idx.at(2)));

    /* ~~~~~~~~ REACTOR-NU ~~~~~~~~ */

    // Create reactor norm variables (allow to vary by ±3 sigma), and model
    std::vector<std::string> reactor_names = {"BRUCE", "DARLINGTON", "PICKERING", "WORLD"};

    std::vector<double> reac_hist_fracs;
    compute_reac_unosc_fracs(reactor_hists, reactor_names, reac_hist_fracs);
    std::vector<std::string> ReactorNorms_VarNames;
    double reacNorm;
    double reacLowLim;
    for (unsigned int i = 0; i < reac_hist_fracs.size(); ++i) {
        // Independent scaling factor of reactor PDFs (with norms = number of expected events)
        reacNorm = reac_hist_fracs.at(i) * N_IBD;
        reacLowLim = (1.0 - 3.0 * IBD_err_indiv) * reacNorm;
        if (reacLowLim < 0.0) reacLowLim = 0.0;
        antinuFitter.GetVars().AddVar("reactorNorm_" + reactor_names.at(i), reacNorm, IBD_err_indiv * reacNorm, reacLowLim, (1.0 + 3.0 * IBD_err_indiv) * reacNorm);
        ReactorNorms_VarNames.push_back("reactorNorm_" + reactor_names.at(i));
    }
    // Extra overall normalisation (= 1), to add shared unceetainties 
    reacLowLim = 1.0 - 3.0 * IBD_err_tot;
    if (reacLowLim < 0.0) reacLowLim = 0.0;
    antinuFitter.GetVars().AddVar("reactorNorm_tot", 1.0, IBD_err_tot, reacLowLim, 1.0 + 3.0 * IBD_err_tot);

    // Add reactor model, linking it to approproate variables and E-systematics defined above
    antinuFitter.AddReactorMod("deltamsqr21", "deltamsqr32", "sinsqrtheta12", "sinsqrtheta13", ReactorNorms_VarNames,
                               "reactorNorm_tot", "EsysBeta", reactor_hists, &E_conv, reactor_names, db);
    antinuFitter.GetModels().at(2)->hold_osc_params_const(true); // This will also compute oscillated reactor specs


    /* ~~~~~~~~ AZIMOV DATASET ~~~~~~~~ */

    // Make fake dataset out of PDF hists (same function called to make PDFs)
    std::cout << "Creating fake dataset..." << std::endl;
    std::vector<TH1D*> hists;
    antinuFitter.GetAllSpectra(hists);
    
    TH1D* data = (TH1D*)hists.at(0)->Clone("data");
    data->SetTitle("Azimov Dataset");
    data->Reset("ICES");
    std::cout << "data integral = " << data->Integral() << std::endl;

    // Add spectra (already re-sclaed)
    for (unsigned int i = 0; i < hists.size(); ++i) {
        std::cout << hists.at(i)->GetName() << " integral = " << hists.at(i)->Integral() << std::endl;
        data->Add(hists.at(i));  // Add reactor events
    }
    std::cout << "data integral = " << data->Integral() << std::endl;

    // Set bins outside "real data" cuts to zero
    double dE = data->GetBinCenter(2) - data->GetBinCenter(1);
    double bin_centre, bin_top, bin_bottom;
    for (unsigned int i = 1; i < data->GetXaxis()->GetNbins() + 1; ++i) {
        bin_centre = data->GetBinCenter(i);
        bin_top = bin_centre + 0.5 * dE;
        bin_bottom = bin_centre - 0.5 * dE;

        if (bin_bottom >= 0.7 && bin_top <= 8.0) {
            continue;
        } else if (bin_top <= 0.7 || bin_bottom >= 0.8) {
            data->SetBinContent(i, 0.0);
        } else if (bin_bottom < 0.7) {
            data->SetBinContent(i, data->GetBinContent(i) * (0.7 - bin_bottom) / dE);
        } else {
            data->SetBinContent(i, data->GetBinContent(i) * (bin_top - 8.0) / dE);
        }
    }

    std::cout << "data integral = " << data->Integral() << std::endl;

    // Add Azimov dataset as data
    antinuFitter.SetData(data);
}


/* ~~~~~~~~ OTHER SHARED UTILS ~~~~~~~~ */

void compute_hist_fracs(const std::vector<TH1D*>& hists, const std::vector<std::string>& hist_names, std::vector<double>& hist_fracs, std::vector<unsigned int>& hist_idx) {

    if (hists.size() != hist_names.size()) {
        std::cout << "Error! hists and hist_names are not the same size: " << hists.size() << ", and " << hist_names.size() << std::endl;
        exit(1);
    }

    hist_fracs.resize(hist_names.size());
    hist_idx.resize(hist_names.size());

    double tot_int;
    double integral;
    bool hist_not_found;
    std::string name;
    for (unsigned int i = 0; i < hists.size(); ++i) {
        name = hists.at(i)->GetName();
        integral = hists.at(i)->Integral();
        // See if hist is in hist_names
        hist_not_found = true;
        for (unsigned int j = 0; j < hist_names.size(); ++j) {
            if (name == hist_names.at(j)) {
                hist_fracs.at(j) = integral;
                hist_idx.at(j) = i;
                hist_not_found = false;
                break;
            }
        }
        if (hist_not_found) {
            std::cout << "Error! Histogram " << name << " not found" << std::endl;
            exit(1);
        }

        tot_int += integral;
    }

    // Convert integrals to fractions of total integral
    for (unsigned int i = 0; i < hist_fracs.size(); ++i) {
        hist_fracs.at(i) /= tot_int;
    }
}

void compute_reac_unosc_fracs(const std::vector<TH1D*>& Reactor_hists, const std::vector<std::string>& Reactor_names, std::vector<double>& hist_fracs) {

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
            if (name == "geoNu_Th" || name == "geoNu_U") {
                geoNu_hists.push_back((TH1D*)obj);
            } else if (name == "alphaN_PR" || name == "alphaN_C12" || name == "alphaN_O16") {
                alphaN_hists.push_back((TH1D*)obj);
            } else if (name == "E_conversion") {
                continue;
            } else {
                reactor_hists.push_back((TH1D*)obj);
            }
        }
    }

    // Get 2D hist
    E_conv = *(TH2D*)(fin->Get("E_conversion"));

    // Error handling
    if (reactor_hists.size() == 0) {
        std::cout << "ERROR: No reactor IBD histograms!" << std::endl;
        exit(1);
    }
    if (alphaN_hists.size() < 3) {
        std::cout << "ERROR: Not enough alpha-n histograms! Only " << alphaN_hists.size() << std::endl;
        exit(1);
    }
    if (geoNu_hists.size() < 2) {
        std::cout << "ERROR: No geo-nu histograms! Only " << geoNu_hists.size() << std::endl;
        exit(1);
    }
}

