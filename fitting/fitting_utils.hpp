// header guard:
#ifndef fitting_utils
#define fitting_utils

// include
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
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"
#include "model_Reactor.hpp"
#include "model.hpp"
#include "cutting_utils.hpp"


void create_fitter(std::string PDFs_address, double Dm21_2, double Dm32_2, double s_12_2, double s_13_2, RAT::DB* db, const bool useAzimovData = true);
void compute_hist_fracs(const std::vector<TH1D*>& hists, const std::vector<std::string>& hist_names, std::vector<double>& hist_fracs, std::vector<unsigned int>& hist_idx, const unsigned int min_bin, const unsigned int max_bin);
void read_hists_from_file(TFile* fin, std::vector<TH1D*>& reactor_hists, std::vector<TH1D*>& alphaN_hists, std::vector<TH1D*>& geoNu_hists, std::vector<TH1D*>& Accidental_hists, TH2D& E_conv);

/* ~~~~~~~~ CONSTRAINED PARAMETERS ~~~~~~~~ */

double N_IBD = 52.2;            // Total number of expected reactor IBDs (at 30000 times rate, Raw entries: 2391909, scaled entries: 1566636, ratio: 0.654973077989171)
// double N_IBD = 52.2 * 0.96;     // Classifier cut
double IBD_err_indiv = 0;
double IBD_err = 0.03;      // fractional error in N_IBD for total reactor IBDs

double N_alphaN = 18.2;         // Total number of expected alpha-n
// double N_alphaN = 18.2 * 0.22;  // Classifier cut
double alphaN_err_PR = 0;     // fractional error in N_alphaN for PR, on top of GS uncertainty
double alphaN_err_GS = 0.3;     // fractional error in N_alphaN for ground state neutrons (PR & 12C)
double alphaN_err_ES = 1.0;     // fractional error in N_alphaN for excited state neutrons (O16)

double N_geoNu = 12.5;          // Total number of expected geo-nu IBDs (un-oscillated, 72% cut efficiency)
// double N_geoNu = 12.5 * 0.89;   // Classifier cut
// double geoNu_err = 1.0;         // fractional error in N_geoNu for individual Th and U spectra
double geoNuUThRatio = 3.7;
double geoNuRatio_err = 0.35;  // fractional error

double N_acc = 0.3;               // Total number of expected accidentals coincidence events (no err)
double N_side = 1.1;               // Total number of expected sideband events
double sideband_err = 1.0;      // fractional error

double linScale_err = 0.011;    // Error in linear scaling (scaling = 1) (not fractional)
double kB = 0.074;              // Birk's constant for betas
double kB_err = 0.004;          // Error in kB (not fractional)
double sigPerSqrtE = 0.042;     // smearing sigma = sigPerSqrtE * sqrt(E)

double linScale_err_P = 0.011;  // Error linear scaling for proton recoils (scaling = 1) (not fractional)
double kB_P = 0.078;            // Birk's constant for protons
double kB_err_P = 0;            // Error in kB_P for proton recoils (not fractional)
// Proton recoil uses the same smearing as the rest

/* ~~~~~~~~ SET UP FITTER ~~~~~~~~ */

void create_fitter(std::string PDFs_address, double Dm21_2, double Dm32_2, double s_12_2, double s_13_2, RAT::DB* db, TTree* Data) {
    create_fitter(PDFs_address, Dm21_2, Dm32_2, s_12_2, s_13_2, db, false);
    Fitter* antinuFitter = Fitter::GetInstance();
    antinuFitter->SetData(Data);
}

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
void create_fitter(std::string PDFs_address, const double Dm21_2, const double Dm32_2, const double s_12_2, const double s_13_2, RAT::DB* db, const bool useAzimovData) {
    // Read in file
    std::cout << "Reading in hists from file..." << std::endl;

    TFile *fin = TFile::Open(PDFs_address.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file." << std::endl;
        exit(1);
    }
    // Check for errors
    TList* list = fin->GetListOfKeys() ;
    if (!list) {std::cout << "No keys found in file\n" << std::endl; exit(1);}

    // Get hists
    std::vector<TH1D*> reactor_hists;
    std::vector<TH1D*> alphaN_hists;
    std::vector<TH1D*> geoNu_hists;
    std::vector<TH1D*> Accidental_hists;
    TH2D E_conv;

    read_hists_from_file(fin, reactor_hists, alphaN_hists, geoNu_hists, Accidental_hists, E_conv);

    // Get other reactor info
    TTree *inTree = (TTree*) fin->Get("reactor_vals");

    Double_t PWR_promptE_frac;
    std::vector<Double_t> PWR_Enu_fracs, PWR_Enu_baselines;
    std::vector<Double_t> PHWR_Enu_fracs, PHWR_Enu_baselines;

    inTree->SetBranchAddress("PWR_promptE_frac", &PWR_promptE_frac);
    inTree->SetBranchAddress("PWR_Enu_fracs", &PWR_Enu_fracs);
    inTree->SetBranchAddress("PWR_Enu_baselines", &PWR_Enu_baselines);
    inTree->SetBranchAddress("PHWR_Enu_fracs", &PHWR_Enu_fracs);
    inTree->SetBranchAddress("PHWR_Enu_baselines", &PHWR_Enu_baselines);

    // Find data bin limits
    double bin_centre;
    unsigned int min_bin = reactor_hists.at(0)->GetXaxis()->GetNbins();
    unsigned int max_bin = 0;
    for (unsigned int i = 1; i < reactor_hists.at(0)->GetXaxis()->GetNbins() + 1; ++i) {
        bin_centre = reactor_hists.at(0)->GetBinCenter(i);
        if (bin_centre > IBD_MIN_PROMPT_E && bin_centre < IBD_MAX_PROMPT_E) {
            if (i < min_bin) min_bin = i;
            if (i > max_bin) max_bin = i;
        }
    }

    /* ~~~~~~~~ INITIALISE ~~~~~~~~ */

    FitVars* Vars = FitVars::GetInstance();
    Esys* Esysts = Esys::GetInstance();
    Reactor* ReactorMod = Reactor::GetInstance();
    alphaN* alphaNMod = alphaN::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();
    Model* BasicMods = Model::GetInstance();

    /* ~~~~~~~~ OSCILLATION CONSTANTS ~~~~~~~~ */

    // Add oscillation variables
    Vars->AddVar("deltamsqr21", Dm21_2, 0, Dm21_2, Dm21_2, true, false);
    Vars->AddVar("deltamsqr32", Dm32_2, 0, Dm32_2, Dm32_2, true, false);
    Vars->AddVar("sinsqrtheta12", s_12_2, 0, s_12_2, s_12_2, true, false);
    Vars->AddVar("sinsqrtheta13", s_13_2, 0, s_13_2, s_13_2, true, false);


    /* ~~~~~~~~ ENERGY SYSTEMATICS ~~~~~~~~ */

    // Add beta variables
    double linScale_min = 1.0 - 3.0 * linScale_err;
    double kBp_min = kB - 3.0 * kB_err;
    if (linScale_min < 0.0) linScale_min = 0.0;
    if (kBp_min < 0.0) kBp_min = 0.0;
    Vars->AddVar("linScale", 1.0, linScale_err, linScale_min, 1.0 + 3.0 * linScale_err);
    Vars->AddVar("kBp", kB, kB_err, kBp_min, kB + 3.0 * kB_err);

    // Add extra proton variables (linear scaling is the same, but independent)
    Vars->AddVar("linScale_P", 1.0, linScale_err, linScale_min, 1.0 + 3.0 * linScale_err);
    if (kB_err_P == 0) {
        // No extra proton non-linear scaling
        Vars->AddVar("kBp_P", kB_P, 0.0, kB_P, kB_P, true, false);
    } else {
        double kBp_P_min = kB_P - 3.0 * kB_err_P;  // same error, but different Birk's constant
        if (kBp_P_min < 0.0) kBp_P_min = 0.0;
        Vars->AddVar("kBp_P", kB_P, kB_err_P, kBp_P_min, kB_P + 3.0 * kB_err_P);
    }

    // Add smearing variable (same for all, allow to be negative, absolute value taken in calculations)
    Vars->AddVar("sigPerSqrtE", 0.0, sigPerSqrtE, -3.0 * sigPerSqrtE, 3.0 * sigPerSqrtE);

    // Add energy systamtics objects (beta and proton), and link them to variables above
    Esysts->AddEsys("EsysBeta", kB, "linScale", "kBp", "sigPerSqrtE");
    Esysts->AddEsys("EsysProton", kB_P, "linScale_P", "kBp_P", "sigPerSqrtE");


    /* ~~~~~~~~ GEO-NU ~~~~~~~~ */

    // Get fraction of (unoscillated) geo-nu flux coming from each histogram (Th and U, in order)
    std::vector<std::string> geoNu_names = {"geoNu_Th", "geoNu_U"};
    std::vector<double> geoNu_hist_fracs;
    std::vector<unsigned int> geoNu_hist_idx;
    compute_hist_fracs(geoNu_hists, geoNu_names, geoNu_hist_fracs, geoNu_hist_idx, min_bin, max_bin);

    // Create geo-nu norm variables (allow to vary by ±3 sigma), and model
    double geoNuUThRatio_min = (1.0 - 3.0 * geoNuRatio_err) * geoNuUThRatio;
    if (geoNuUThRatio_min < 0) geoNuUThRatio_min = 0;
    Vars->AddVar("geoNuUThRatio", geoNuUThRatio, geoNuRatio_err * geoNuUThRatio, geoNuUThRatio_min, (1.0 + 3.0 * geoNuRatio_err) * geoNuUThRatio);
    Vars->AddVar_unconstrained("geoNuNorm", N_geoNu, 0.0, 3.0 * N_geoNu);

    // Add geo-nu model, linking it to approproate variables and E-systematics defined above
    geoNuMod->InitGeoNu("geoNuNorm", "geoNuUThRatio", "sinsqrtheta12", "sinsqrtheta13", "EsysBeta", geoNu_hists.at(geoNu_hist_idx.at(0)), geoNu_hists.at(geoNu_hist_idx.at(1)));
    geoNuMod->hold_osc_params_const(true);  // This will also pre-compute the survival prob ahead of time


    /* ~~~~~~~~ ALPHA-N ~~~~~~~~ */

    // Get fraction of alpha-n flux coming from each histogram (in order: proton recoil, C12 scatter, O16 de-excitation)
    std::vector<std::string> alphaN_names = {"alphaN_PR", "alphaN_C12", "alphaN_O16"};
    std::vector<double> alphaN_hist_fracs;
    std::vector<unsigned int> alphaN_hist_idx;
    compute_hist_fracs(alphaN_hists, alphaN_names, alphaN_hist_fracs, alphaN_hist_idx, min_bin, max_bin);
    double GS_frac = alphaN_hist_fracs.at(0) + alphaN_hist_fracs.at(1);
    double ES_frac = alphaN_hist_fracs.at(2);

    // Create alpha-n norm variables (allow to vary by ±3 sigma), and model
    double PR_frac_min = (1.0 - 3.0 * alphaN_err_PR);
    if (PR_frac_min < 0.0) PR_frac_min = 0.0;
    double GS_Norm_min = (1.0 - 3.0 * alphaN_err_GS) * N_alphaN * GS_frac;
    if (GS_Norm_min < 0.0) GS_Norm_min = 0.0;
    double ES_Norm_min = (1.0 - 3.0 * alphaN_err_ES) * N_alphaN * ES_frac;
    if (ES_Norm_min < 0.0) ES_Norm_min = 0.0;

    Vars->AddVar("alphaNNorm_PR", 1.0, alphaN_err_PR, PR_frac_min, 1.0 + 3.0 * alphaN_err_PR);
    Vars->AddVar("alphaNNorm_GS", N_alphaN * GS_frac, alphaN_err_GS * N_alphaN * GS_frac, GS_Norm_min, (1.0 + 3.0 * alphaN_err_GS) * N_alphaN * GS_frac);
    Vars->AddVar("alphaNNorm_ES", N_alphaN * ES_frac, alphaN_err_ES * N_alphaN * ES_frac, ES_Norm_min, (1.0 + 3.0 * alphaN_err_ES) * N_alphaN * ES_frac);

    for (unsigned int i = 0; i < 3; ++i) {
        std::cout << "alphaN_hists.at(" << i << ")->GetName() = " << alphaN_hists.at(i)->GetName() << std::endl;
        std::cout << "alphaN_hists.at(" << i << ")->Integral(min_bin, max_bin)  = " << alphaN_hists.at(i)->Integral(min_bin, max_bin) << std::endl;
    }

    // Add alpha-n model, linking it to approproate variables and E-systematics defined above
    alphaNMod->InitAlphaN("alphaNNorm_PR", "alphaNNorm_GS", "alphaNNorm_ES", "EsysBeta", "EsysProton", alphaN_hists.at(alphaN_hist_idx.at(0)), alphaN_hists.at(alphaN_hist_idx.at(1)), alphaN_hists.at(alphaN_hist_idx.at(2)));

    /* ~~~~~~~~ REACTOR-NU ~~~~~~~~ */

    // Create reactor norm variables (allow to vary by ±3 sigma), and model
    std::vector<std::string> reactor_names = {"BRUCE", "DARLINGTON", "PICKERING", "WORLD"};

    // Extra overall normalisation (= 1), to add shared unceetainties 
    double reacLowLim = (1.0 - 3.0 * IBD_err) * N_IBD;
    if (reacLowLim < 0.0) reacLowLim = 0.0;
    Vars->AddVar("reactorNorm", N_IBD, IBD_err * N_IBD, reacLowLim, (1.0 + 3.0 * IBD_err) * N_IBD);

    // Add reactor model, linking it to approproate variables and E-systematics defined above
    ReactorMod->InitReactor("deltamsqr21", "deltamsqr32", "sinsqrtheta12", "sinsqrtheta13",
                            "reactorNorm", "EsysBeta", reactor_hists, &E_conv, PWR_promptE_frac,
                            PWR_Enu_fracs, PHWR_Enu_fracs, PWR_Enu_baselines, PHWR_Enu_baselines, db);
    ReactorMod->hold_osc_params_const(true); // This will also compute oscillated reactor specs


    /* ~~~~~~~~ ACCIDENTALS ~~~~~~~~ */
    
    Vars->AddVar("AccidentalsNorm", N_acc, 0, N_acc, N_acc, true, false); // no normalisation error, since it is data driven
    Esysts->AddEsys_trivial("trivial");  // no energy systematics either, for the same reason
    // Add accidentals model, linking it to approproate variables and E-systematics defined above
    BasicMods->AddModel("AccidentalsNorm", "trivial", Accidental_hists.at(0), "Accidentals");

    /* ~~~~~~~~ INITIALISE FITTER ~~~~~~~~ */

    Fitter* antinuFitter = Fitter::GetInstance();
    antinuFitter->SetBinLims(min_bin, max_bin);  // should do this before making Azimov data


    /* ~~~~~~~~ AZIMOV DATASET ~~~~~~~~ */

    if (useAzimovData) {
        // Make fake dataset out of PDF hists (same function called to make PDFs)
        std::cout << "Creating fake dataset..." << std::endl;

        ReactorMod->compute_spec();
        alphaNMod->compute_spec();
        geoNuMod->compute_spec();

        TH1D* data = (TH1D*)(ReactorMod->GetModelEsys()->Clone("data"));
        data->SetTitle("Azimov Dataset");
        data->Reset("ICES");
        
        std::cout << "data integral = " << data->Integral(min_bin, max_bin) << std::endl;
        data->Add(ReactorMod->GetModelEsys());
        std::cout << "data integral = " << data->Integral(min_bin, max_bin) << std::endl;
        data->Add(alphaNMod->GetModelEsys());
        std::cout << "data integral = " << data->Integral(min_bin, max_bin) << std::endl;
        data->Add(geoNuMod->GetModelEsys());
        std::cout << "data integral = " << data->Integral(min_bin, max_bin) << std::endl;

        for (unsigned int iMod = 0; iMod < BasicMods->GetNumMods(); ++iMod) {
            BasicMods->compute_spec(iMod);
            data->Add(BasicMods->GetModelEsys(iMod));
            std::cout << "data integral = " << data->Integral(min_bin, max_bin) << std::endl;
        }

        // Set bins outside "real data" cuts to zero
        for (unsigned int i = 1; i < min_bin; ++i) data->SetBinContent(i, 0.0);
        for (unsigned int i = max_bin + 1; i < data->GetXaxis()->GetNbins() + 1; ++i) data->SetBinContent(i, 0.0);

        std::cout << "data integral = " << data->Integral(min_bin, max_bin) << std::endl;

        // Add Azimov dataset as data
        antinuFitter->SetData(data);
    }
}


/* ~~~~~~~~ OTHER SHARED UTILS ~~~~~~~~ */

void compute_hist_fracs(const std::vector<TH1D*>& hists, const std::vector<std::string>& hist_names, std::vector<double>& hist_fracs, std::vector<unsigned int>& hist_idx, const unsigned int min_bin, const unsigned int max_bin) {

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
        integral = hists.at(i)->Integral(min_bin, max_bin);
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

/**
 * @brief Lists all histograms from root file into vectors of reactor, alphaN, geoNu and E_conv hists
 * 
 * @param file_address 
 * @param reactor_hists 
 * @param alphaN_hists
 * @param geoNu_hists
 * @param E_conv
 */
void read_hists_from_file(TFile* fin, std::vector<TH1D*>& reactor_hists, std::vector<TH1D*>& alphaN_hists, std::vector<TH1D*>& geoNu_hists, std::vector<TH1D*>& Accidental_hists, TH2D& E_conv) {

    // Get 1-D hists
    geoNu_hists.push_back((TH1D*)fin->Get("geoNu_Th"));
    geoNu_hists.push_back((TH1D*)fin->Get("geoNu_U"));

    alphaN_hists.push_back((TH1D*)fin->Get("alphaN_PR"));
    alphaN_hists.push_back((TH1D*)fin->Get("alphaN_C12"));
    alphaN_hists.push_back((TH1D*)fin->Get("alphaN_O16"));

    Accidental_hists.push_back((TH1D*)fin->Get("Accidental"));

    reactor_hists.push_back((TH1D*)fin->Get("PWR_promptE"));
    reactor_hists.push_back((TH1D*)fin->Get("PWR_Enu"));
    reactor_hists.push_back((TH1D*)fin->Get("PHWR_Enu"));

    // Get 2D hist
    E_conv = *(TH2D*)(fin->Get("E_conversion"));
}

void setup_param_hists(TH2D* minllHist, std::vector<TH2D*>& par_hists) {
    FitVars* Vars = FitVars::GetInstance();

    for (unsigned int iVar = 0; iVar < Vars->GetNumVars(); ++iVar) {
        par_hists.push_back((TH2D*)minllHist->Clone((Vars->name(iVar)).c_str()));
        par_hists.at(iVar)->SetTitle((Vars->name(iVar)).c_str());
        par_hists.at(iVar)->Reset("ICES");
    }
}

//end header guard
#endif