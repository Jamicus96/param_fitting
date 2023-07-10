#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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
#include "fitting_utils.hpp"


TH2D* Fit_spectra(reactorINFO& spectrum, TH1D& alphaN_hist, TH1D* data);
double ML_fit(reactorINFO& spectrum, RooHistPdf& alphaN_PDF, RooDataHist& dataHist, const RooRealVar& E, const RooRealVar& reactor_frac, const RooRealVar& alphaN_frac);
void read_hists_from_file(std::string file_address, std::vector<TH1D*>& reactor_hists, TH1D& alphaN_hist);
void get_baselines(RAT::DB *db, std::vector<TH1D*>& reactor_hists, std::vector<double>& baslines);
double GetReactorDistanceLLA(const double &longitude, const double&latitude, const double &altitude);
TVector3 LLAtoECEF(double longitude, double latitude, double altitude);
std::vector<std::string> SplitString(std::string str);


int main(int argv, char** argc) {
    std::string PDFs_address = argc[1];
    std::string out_address = argc[2];

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
    std::vector<double> L;
    get_baselines(db, reactor_hists, L);

    // Get oscillation constants
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    const double fDmSqr32 = linkdb->GetD("deltamsqr32");
    const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
    const double fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    // Create PDF fitting object
    reactorINFO spectrum(reactor_hists, L, fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);

    // Make fake dataset out of PDF hists (same function called to make PDFs)
    std::cout << "Creating fake dataset..." << std::endl;
    spectrum.compute_osc_reactor_spec();
    std::vector<TH1D*> osc_hist = spectrum.Get_osc_reactor_specs();
    TH1D* data = (TH1D*)osc_hist.at(0)->Clone();
    data->Add(osc_hist.at(0), 1.0);  // Add reactor events
    data->Add(&alphaN_hist, 1.0);  // Add alpha-n events

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;
    TH2D* minllHist = Fit_spectra(spectrum, alphaN_hist, data);

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
 * @param spectrum reactorINFO object that contains all the relevent reactor
 * @param alphaN_hist alphaN data
 * @param data data to fit model to
 * @return TH2D* 
 */
TH2D* Fit_spectra(reactorINFO& spectrum, TH1D& alphaN_hist, TH1D* data) {

    //declare observables
    RooRealVar E("E", "energy", 0.9, 8);
    RooRealVar reactor_frac("reactor_frac", "reactorIBD fraction", 0.4, 0.8);
    RooRealVar alphaN_frac("alphaN_frac", "alpha-n fraction", 0.2, 0.6);

    // Make data hist (data that must be fit)
    RooDataHist dataHist("dataHist", "data hist", E, data);

    // Create alphaN PDF
    RooDataHist tempData_alpha("tempData", "temporary data", E, &alphaN_hist);
    RooHistPdf alphaN_PDF("alphaN_PDF", "alphaN PDF", E, tempData_alpha);

    // Define varying parameters (evenly spaced. plot theta_12, but use s_12^2 in calculations)
    std::vector<double> Dm21;
    std::vector<double> sinTheta12;
    unsigned int N_steps = 500;
    double Dm21_min = 1E-5; double Dm21_max = 10E-5;
    double Theta12_min = 5.; double Theta12_max = 45.;  // degrees
    double Dm21_step = (Dm21_max - Dm21_min) / (double)N_steps;
    double Theta12_step = (Theta12_max - Theta12_min) / (double)(N_steps - 1);
    for (unsigned int n = 0; n < N_steps; ++n) {
        Dm21.push_back(Dm21_min + (double)n * Dm21_step);
        sinTheta12.push_back(pow(sin((Theta12_min + (double)n * Theta12_step)  * TMath::Pi() / 180.), 2));
    }

    // Create 2-d hist to dump data into
    double Dm21_lower = Dm21_min - 0.5 * Dm21_step;
    double Theta12_lower = Theta12_min - 0.5 * Theta12_step;
    double Dm21_upper = Dm21_max + 0.5 * Dm21_step;
    double Theta12_upper = Theta12_max + 0.5 * Theta12_step;
    
    TH2D* minllHist = new TH2D("minllHist", "minimised likelihood values", sinTheta12.size(), Theta12_lower, Theta12_upper, Dm21.size(), Dm21_lower, Dm21_upper);

    // Compute best fit Log likelihood for each set of paramters (fraction of alpha-n vs reactor IBD events is fit in each loop)
    double MLL;
    std::cout << "Looping over oscillation parameters..." << std::endl;
    for (unsigned int i = 0; i < sinTheta12.size(); ++i) {
        // std::cout << "i = " << i << std::endl;
        spectrum.s12_2() = sinTheta12.at(i);
        for (unsigned int j = 0; j < Dm21.size(); ++j) {
            // std::cout << "j = " << j << std::endl;
            spectrum.Dm21_2() = Dm21.at(j);
            minllHist->SetBinContent(i+1, j+1, ML_fit(spectrum, alphaN_PDF, dataHist, E, reactor_frac, alphaN_frac));
        }
    }

    // Axis titles
    minllHist->GetXaxis()->SetTitle("Theta_12");
    minllHist->GetYaxis()->SetTitle("Delta m_21^2");
    minllHist->GetYaxis()->SetTitleOffset(1); // Move x-axis label further away from the axis (0 is default)

    // Remove stats box
    minllHist->SetStats(0);

    return minllHist;
}

double ML_fit(reactorINFO& spectrum, RooHistPdf& alphaN_PDF, RooDataHist& dataHist, const RooRealVar& E, const RooRealVar& reactor_frac, const RooRealVar& alphaN_frac) {

    // Compute oscillated reactor spectrum
    spectrum.compute_osc_reactor_spec();

    // Create PDFs
    RooDataHist tempData_react("tempData", "temporary data", E, spectrum.Get_osc_reactor_specs().at(0));
    RooHistPdf reactor_PDF("reactor_PDF", "reactor PDF", E, tempData_react);

    // Add PDFs together
    RooAddPdf model = RooAddPdf("model", "r+a", RooArgList(reactor_PDF, alphaN_PDF), RooArgList(reactor_frac, alphaN_frac));

    // Fit to data
    return model.fitTo(dataHist, Extended(true), PrintLevel(-1), SumW2Error(kFALSE), Save())->minNll();
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

/**
 * @brief Split reactor name strings for getting baselines
 * 
 * @param str 
 * @return std::vector<std::string> 
 */
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
void get_baselines(RAT::DB *db, std::vector<TH1D*>& reactor_hists, std::vector<double>& baselines) {

    RAT::DBLinkPtr linkdb;
    std::vector<std::string> originReactorVect;
    std::vector<Double_t> fLatitude;
    std::vector<Double_t> fLongitute;
    std::vector<Double_t> fAltitude;
    for (unsigned int n = 0; n < reactor_hists.size(); ++n) {
        originReactorVect = SplitString(reactor_hists.at(n)->GetName());
        linkdb = db->GetLink("REACTOR",originReactorVect[0]);
        fLatitude  = linkdb->GetDArray("latitude");
        fLongitute = linkdb->GetDArray("longitude");
        fAltitude = linkdb->GetDArray("altitude");

        baselines.push_back(GetReactorDistanceLLA(fLongitute[std::stoi(originReactorVect[1])], fLatitude[std::stoi(originReactorVect[1])], fAltitude[std::stoi(originReactorVect[1])]));
    }
}
