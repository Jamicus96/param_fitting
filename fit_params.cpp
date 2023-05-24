#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TVector3.h>
#include <TTree.h>
#include <sstream>
#include <RAT/DB.hh>
#include <TRandom3.h>
#include <TKey.h>
#include <TObject.h>
#include <TList.h>
#include <TVectorD.h>
#include "Math/DistFunc.h"

using namespace RooFit;


std::vector<std::string> SplitString(std::string str){
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

double GetReactorDistanceLLA(const double &longitude, const double&latitude, const double &altitude) {
    const TVector3 SNO_ECEF_coord_ = TVector3(672.87,-4347.18,4600.51);
    double dist = (LLAtoECEF(longitude, latitude,altitude) - SNO_ECEF_coord_).Mag();
  return dist;
}

std::vector<double> get_baselines(RAT::DB *db, std::vector<TH1D*> hists) {

    std::vector<double> baselines;
    RAT::DBLinkPtr linkdb;
    std::vector<std::string> originReactorVect;
    std::vector<Double_t> fLatitude;
    std::vector<Double_t> fLongitute;
    std::vector<Double_t> fAltitude;
    for (unsigned int n = 0; n < hists.size(); ++n) {
        originReactorVect = SplitString(hists.at(n)->GetName());
        linkdb = db->GetLink("REACTOR",originReactorVect[0]);
        fLatitude  = linkdb->GetDArray("latitude");
        fLongitute = linkdb->GetDArray("longitude");
        fAltitude = linkdb->GetDArray("altitude");

        baselines.push_back(GetReactorDistanceLLA(fLongitute[std::stoi(originReactorVect[1])], fLatitude[std::stoi(originReactorVect[1])], fAltitude[std::stoi(originReactorVect[1])]));
    }

    return baselines;
}

/**
 * @brief Lists all TH1D histograms from root file into a vector.
 * 
 * @param file_address 
 * @return std::vector<TH1D*> 
 */
std::vector<std::vector<TH1D*>> read_hists_from_file(std::string file_address) {

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
    std::vector<TH1D*> reactor_hists;
    std::vector<TH1D*> alphaN_hists;

    // Go through list of histograms and add them to temp list if they are included in name list
    while((key = (TKey*)next())){
        obj = key->ReadObj() ;
        if(obj->InheritsFrom(TH1::Class())){
            // Check which histogram in file matches name
            name = obj->GetName();
            if (name == "alphaN") {
                alphaN_hists.push_back((TH1D*)obj);
            } else {
                reactor_hists.push_back((TH1D*)obj);
            }
        }
    }

    // Error handling
    if (reactor_hists.size() == 0) {
        std::cout << "ERROR: No reactor IBD histograms!" << std::endl;
        exit(1);
    }
    if (alphaN_hists.size() == 0) {
        std::cout << "ERROR: No alpha-n histograms!" << std::endl;
        exit(1);
    }

    return {reactor_hists, alphaN_hists};
}

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

    // Define electron density Ne of the crust, based on 2.7g/cm3 mass density, and <N/A> = 0.5
    const double alpha = - 2.535e-31 * 8.13e23;  // conversion factor in eV2/MeV * Ne = 8.13e23

    return {H_ee_vac, a0_vac, a1_vac, Y_ee_vac, alpha};
}

/**
 * @brief Re-compute oscillation constants to take into account matter effects.
 * Essencially computes everything possible in oscillations that don't requite the baseline.
 * 
 * @param E  antineutrino recon energy (MeV)
 * @param params = {H_ee_vac, a0_vac, a1_vac, Y_ee_vac, alpha}
 * @return std::vector<std::vector<double>> 
 */
std::vector<std::vector<double>> re_compute_consts(const double E, const std::vector<double>& params) {
    const double A_CC = params[4] * E; // for A_CC in [eV^2] and nuE in [MeV]
    
    // Compute new values for H_ee, Y, a0 and a1 (make sure and Y are updated after their use by others)
    const double alpha_1 = params[0] * A_CC + (1.0/3.0) * A_CC*A_CC;

    const double a0 = params[1] - params[3] * A_CC - (1.0/3.0) * params[0] * A_CC*A_CC - (2.0/27.0) * A_CC*A_CC*A_CC;
    const double a1 = params[2] - alpha_1;
    const double Y_ee = params[3] + (2.0/3.0) * alpha_1;
    const double H_ee = params[0] + (2.0/3.0) * A_CC;

    // Get eigenvalues of H, and constants X and theta
    const double arcCos = (1.0/3.0) * acos(1.5 * (a0/a1) * sqrt(- 3.0 / a1));
    const double preFact = 2.0 * sqrt(- a1 / 3.0);

    std::vector<double> eigen = {0., 0., 0.};
    std::vector<double> X = {0., 0., 0.};
    for(int i=0; i<3; ++i){
        eigen.at(i) = preFact * cos(arcCos - (2.0/3.0) * M_PI * i);
        X.at(i) = (1.0/3.0) + (eigen.at(i) * H_ee + Y_ee) / (3.0 * eigen.at(i)*eigen.at(i) + a1);
    }

    return {eigen, X};
}

/**
 * @brief Computes survival probability for an electron antineutrino, in curst with constant matter density.
 * 
 * @param E  antineutrino recon energy (MeV)
 * @param L  baseline (km)
 * @param params = {eigen[3], X[3]}
 * @return double 
 */
double survival_prob(const double E, const double L, const std::vector<std::vector<double>>& params) {

    const double scale = 1.267e3 * L / E; // for E in [MeV] and L in [km]

    const double s_10 = sin(scale * (params[0][1] - params[0][0]));
    const double s_20 = sin(scale * (params[0][2] - params[0][0]));
    const double s_21 = sin(scale * (params[0][2] - params[0][1]));

    // Compute probability
    return 4.0 * (params[1][1]*params[1][0]*s_10*s_10 + params[1][2]*params[1][0]*s_20*s_20 + params[1][2]*params[1][1]*s_21*s_21);
}

TH1D* compute_tot_reactor_spec(const std::vector<TH1D*>& hists, const std::vector<double>& L, const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13) {
    
    // Compute oscillation constants
    const std::vector<double> params = compute_oscillation_constants(fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);

    // Create histograms to sum reactor events and alphaN events, to make PDFs
    TH1D* reactor_hist = (TH1D*)(hists.at(0)->Clone());

    // Assume all the histograms have the same E binning
    unsigned int num_hists = hists.size();
    double E;
    double prob;
    std::vector<std::vector<double>> re_consts;
    for (unsigned int i = 1; i < hists.at(0)->GetXaxis()->GetNbins(); ++i) {
        E = hists.at(0)->GetXaxis()->GetBinCenter(i);
        re_consts = re_compute_consts(E, params);
        for (unsigned int j = 0; j < num_hists; ++j) {
            // Compute survival probability for particular energy at reactor baseline
            prob = survival_prob(E, L.at(j), re_consts);
            // Add value of bin from reactor core to total, scaled by survival probability
            reactor_hist->AddBinContent(i, prob * hists.at(j)->GetBinContent(i));
        }
    }

    return reactor_hist;
}

double ML_fit(const RooRealVar& E, const RooRealVar& reactor_frac, const RooRealVar& alphaN_frac, RooDataHist& dataHist, RooHistPdf& alphaN_PDF, const std::vector<TH1D*>& reactor_hists, const std::vector<double>& L, const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13) {

    // Compute total reactor IBD spectrum
    TH1D* reactor_spec = compute_tot_reactor_spec(reactor_hists, L, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);
    // Normalise histograms
    // reactor_hist->Scale(1.0 / reactor_hist->Integral(), "width");

    // Create reactor IBD PDF
    RooDataHist* tempData_react = new RooDataHist("tempData", "temporary data", E, reactor_spec);
    RooHistPdf* reactor_PDF = new RooHistPdf("PDF", "PDF", E, *tempData_react);

    // Make model
    RooAddPdf model("model", "r+a", RooArgList(*reactor_PDF, alphaN_PDF), RooArgList(reactor_frac, alphaN_frac));

    // Fit to data
    RooFitResult *result = model.fitTo(dataHist, Extended(true), PrintLevel(-1), SumW2Error(kFALSE), Save());

    // Get results
    double minll = result->minNll();

    // delete objects
    delete(tempData_react);
    delete(reactor_PDF);

    return minll;
}

TH2D* Fit_spectra(TH1D* data, const std::vector<std::vector<TH1D*>>& hists, const std::vector<double>& L, const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13) {

    //declare observables
    RooRealVar E("E", "energy", 0.9, 8);
    RooRealVar reactor_frac("reactor_frac", "reactorIBD fraction", 0.4, 0.8);
    RooRealVar alphaN_frac("alphaN_frac", "alpha-n fraction", 0.2, 0.6);

    // Make data hist (data that must be fit)
    RooDataHist dataHist("dataHist", "data hist", E, data);

    // Create alphaN PDF
    RooDataHist* tempData_alpha = new RooDataHist("tempData", "temporary data", E, hists.at(1).at(0));
    RooHistPdf* alphaN_PDF = new RooHistPdf("PDF", "PDF", E, *tempData_alpha);

    // Define varying parameters (should be evenly spaced, and contain "true" value ideally)
    std::vector<double> Dm21 = {0.2*fDmSqr21, fDmSqr21, 1.8*fDmSqr21};
    std::vector<double> Theta12 = {0.2*fSSqrTheta12, fSSqrTheta12, 1.8*fSSqrTheta12};

    // Create 2-d hist to dump data into
    double Dm21_step = Dm21.at(1) - Dm21.at(0);
    double Theta12_step = Theta12.at(1) - Theta12.at(0);

    double Dm21_lower = Dm21.at(0) - 0.5 * Dm21_step;
    double Theta12_lower = Theta12.at(0) - 0.5 * Theta12_step;
    double Dm21_upper = Dm21.at(Dm21.size()-1) + 0.5 * Dm21_step;
    double Theta12_upper = Theta12.at(Theta12.size()-1) + 0.5 * Theta12_step;
    
    TH2D* minllHist = new TH2D("minllHist", "minimised likelihood values", Theta12.size(), Theta12_lower, Theta12_upper, Dm21.size(), Dm21_lower, Dm21_upper);

    // Compute best fit Log likelihood for each set of paramters (fraction of alpha-n vs reactor IBD events is fit in each loop)
    double MLL;
    std::cout << "Looping over oscillation parameters..." << std::endl;
    for (unsigned int i = 0; i < Theta12.size(); ++i) {
        for (unsigned int j = 0; j < Dm21.size(); ++j) {
            MLL = ML_fit(E, reactor_frac, alphaN_frac, dataHist, *alphaN_PDF, hists.at(0), L, Dm21[j], fDmSqr32, Theta12[i], fSSqrTheta13);
            minllHist->SetBinContent(i+1, j+1, MLL);
        }
    }

    // delete objects
    delete(tempData_alpha);
    delete(alphaN_PDF);

    return minllHist;
}


int main(int argv, char** argc) {
    std::string PDFs_address = argc[1];
    std::string out_address = argc[2];

    // Read in file
    std::cout << "Reading in hists from file..." << std::endl;
    std::vector<std::vector<TH1D*>> hists = read_hists_from_file(PDFs_address);

    // Get baselines
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    std::cout << "Getting baselines..." << std::endl;
    std::vector<double> L = get_baselines(db, hists.at(0));

    // Get oscillation constants
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    const double fDmSqr32 = linkdb->GetD("deltamsqr32");
    const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
    const double fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    // Make fake dataset out of PDF hists (same function called to make PDFs)
    std::cout << "Creating fake dataset..." << std::endl;
    TH1D* reactor_hist = compute_tot_reactor_spec(hists.at(0), L, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);
    TH1D* data = (TH1D*)reactor_hist->Clone();
    data->Add(reactor_hist, 1.0);  // Add reactor events
    data->Add(hists.at(1).at(0), 1.0);  // Add alpha-n events

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;
    TH2D* minllHist = Fit_spectra(data, hists, L, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);

    // Write hist to file and close
    TFile *outroot = new TFile(out_address.c_str(), "RECREATE");
    minllHist->Write();
    outroot->Write();
    outroot->Close();

    return 0;
}