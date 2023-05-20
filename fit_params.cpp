#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
// #include "RooDataSet.h"
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

/**
 * @brief Lists all TH1D histograms from root file into a vector.
 * 
 * @param file_address 
 * @return std::vector<TH1D*> 
 */
std::vector<TH1D*> read_hists_from_file(std::string file_address) {

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
    std::vector<TH1D*> hist_list;

    // Go through list of histograms and add them to temp list if they are included in name list
    while((key = (TKey*)next())){
        obj = key->ReadObj() ;
        if(obj->InheritsFrom(TH1::Class())){
            // Check which histogram in file matches name
            hist_list.push_back((TH1D*)obj);
        }
    }

    return hist_list;
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


TH1D* compute_tot_PDF(const std::vector<TH1D*>& hists, const double fDmSqr21, const double fDmSqr32, const double fSSqrTheta12, const double fSSqrTheta13, double IBD_frac, double alphaN_frac) {
    
    unsigned int num_hists = hists.size();
    std::string hist_name;
    for (unsigned int n = 0; n < num_hists; ++n) {
        hist_name = hists.at(n)->GetName();
    }

    // std::vector<std::string> originReactorVect = SplitString(originReactorString);
    // linkdb = db->GetLink("REACTOR",originReactorVect[0]);
    // std::vector<Double_t> fLatitude  = linkdb->GetDArray("latitude");
    // std::vector<Double_t> fLongitute = linkdb->GetDArray("longitude");
    // std::vector<Double_t> fAltitude = linkdb->GetDArray("altitude");
    // double baseline = GetReactorDistanceLLA(fLongitute[std::stoi(originReactorVect[1])], fLatitude[std::stoi(originReactorVect[1])], fAltitude[std::stoi(originReactorVect[1])]);
}

int Fit_spectra(const std::vector<TH1D*>& hists) {

    // Get oscillation paramters from ratdb
    RAT::DB *db = RAT::DB::Get();
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    // const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    const double fDmSqr32 = linkdb->GetD("deltamsqr32");
    // const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
    const double fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
}


int main(int argv, char** argc) {
    std::string PDFs_address = argc[1];

    // Read in file
    std::vector<TH1D*> hists = read_hists_from_file(PDFs_address);

    return 0;
}