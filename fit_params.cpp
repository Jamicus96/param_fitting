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


TH2D* Fit_spectra(PDFspec& model, RooDataHist& dataHist, const double fDmSqr32, const double fSSqrTheta13);
std::vector<std::vector<TH1D*>*> read_hists_from_file(std::string file_address);
std::vector<double>* get_baselines(RAT::DB *db, const std::vector<TH1D*>* hists);
double GetReactorDistanceLLA(const double &longitude, const double&latitude, const double &altitude);
TVector3 LLAtoECEF(double longitude, double latitude, double altitude);
std::vector<std::string> SplitString(std::string str);


int main(int argv, char** argc) {
    std::string PDFs_address = argc[1];
    std::string out_address = argc[2];

    // Read in file
    std::cout << "Reading in hists from file..." << std::endl;
    std::vector<std::vector<TH1D*>*> hists = read_hists_from_file(PDFs_address);

    // Get baselines
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    std::cout << "Getting baselines..." << std::endl;
    std::vector<double>* L = get_baselines(db, hists.at(0));

    // Get oscillation constants
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    const double fDmSqr32 = linkdb->GetD("deltamsqr32");
    const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
    const double fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    // Create PDF fitting object
    PDFspec model = PDFspec(hists.at(1)->at(0), hists.at(0), L);

    // Make fake dataset out of PDF hists (same function called to make PDFs)
    std::cout << "Creating fake dataset..." << std::endl;
    model.compute_tot_reactor_spec(fDmSqr21, fDmSqr32, fSSqrTheta12, fSSqrTheta13);
    TH1D* reactor_hist = model.Get_tot_reactor_spec();
    TH1D* data = (TH1D*)reactor_hist->Clone();
    data->Add(reactor_hist, 1.0);  // Add reactor events
    data->Add(hists.at(1)->at(0), 1.0);  // Add alpha-n events
    RooDataHist dataHist("dataHist", "data hist", model.Get_energy_var(), data);

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;
    TH2D* minllHist = Fit_spectra(model, dataHist, fDmSqr32, fSSqrTheta13);

    // Write hist to file and close
    TFile *outroot = new TFile(out_address.c_str(), "RECREATE");
    minllHist->Write();
    outroot->Write();
    outroot->Close();

    return 0;
}


/**
 * @brief Takes reactor and alpha-n data, and fits data to it for a range of Dm_21^2 and theta_12 values
 * 
 * @param model PDFspec object that contains all the relevent reactor and alpha-n information
 * @param dataHist data to fit modelt to
 * @param fDmSqr32 Dm_32^2
 * @param fSSqrTheta13 s_13^2
 * @return TH2D* 
 */
TH2D* Fit_spectra(PDFspec& model, RooDataHist& dataHist, const double fDmSqr32, const double fSSqrTheta13) {

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
        for (unsigned int j = 0; j < Dm21.size(); ++j) {
            // std::cout << "j = " << j << std::endl;
            model.compute_tot_reactor_spec(Dm21.at(j), fDmSqr32, sinTheta12.at(i), fSSqrTheta13);
            minllHist->SetBinContent(i+1, j+1, model.ML_fit(dataHist));
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


/**
 * @brief Lists all TH1D histograms from root file into a vector of vector pointers.
 * 
 * @param file_address 
 * @return std::vector<std::vector<TH1D*>*> = {&{reactor hist 1, reactor hist 2, ...}, &{alpha-n hist}}
 */
std::vector<std::vector<TH1D*>*> read_hists_from_file(std::string file_address) {

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
    std::vector<TH1D*>* reactor_hists;
    std::vector<TH1D*>* alphaN_hists;

    // Go through list of histograms and add them to temp list if they are included in name list
    while((key = (TKey*)next())){
        obj = key->ReadObj() ;
        if(obj->InheritsFrom(TH1::Class())){
            // Check which histogram in file matches name
            name = obj->GetName();
            if (name == "alphaN") {
                alphaN_hists->push_back((TH1D*)obj);
            } else {
                reactor_hists->push_back((TH1D*)obj);
            }
        }
    }

    // Error handling
    if (reactor_hists->size() == 0) {
        std::cout << "ERROR: No reactor IBD histograms!" << std::endl;
        exit(1);
    }
    if (alphaN_hists->size() == 0) {
        std::cout << "ERROR: No alpha-n histograms!" << std::endl;
        exit(1);
    }

    return {reactor_hists, alphaN_hists};
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
 * @param hists 
 * @return std::vector<double>* 
 */
std::vector<double>* get_baselines(RAT::DB *db, const std::vector<TH1D*>* hists) {

    std::vector<double>* baselines;
    RAT::DBLinkPtr linkdb;
    std::vector<std::string> originReactorVect;
    std::vector<Double_t> fLatitude;
    std::vector<Double_t> fLongitute;
    std::vector<Double_t> fAltitude;
    for (unsigned int n = 0; n < hists->size(); ++n) {
        originReactorVect = SplitString(hists->at(n)->GetName());
        linkdb = db->GetLink("REACTOR",originReactorVect[0]);
        fLatitude  = linkdb->GetDArray("latitude");
        fLongitute = linkdb->GetDArray("longitude");
        fAltitude = linkdb->GetDArray("altitude");

        baselines->push_back(GetReactorDistanceLLA(fLongitute[std::stoi(originReactorVect[1])], fLatitude[std::stoi(originReactorVect[1])], fAltitude[std::stoi(originReactorVect[1])]));
    }

    return baselines;
}
