#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TLine.h>
#include <TTree.h>
#include <TVector3.h>
#include <RAT/DB.hh>
#include <TRandom3.h>
#include <TMath.h>
#include <RAT/DU/Point3D.hh>
#include <RAT/DU/Utility.hh>
#include <map>
#include "cutting_utils.hpp"


// Define physical constants
double electron_mass_c2 = 0.510998910;  // [MeV]
double proton_mass_c2 = 938.272013;  // [MeV]
double neutron_mass_c2 = 939.56536;  // [MeV]

double MAX_BASELINE = 1410;  // [km] Anything beyond this averged together


void CreatePDFs_reactorIBD(TTree* EventInfo, TTree* outInfo, std::vector<TH1D*>& PDF_hists, TH2D* E_conv, const double Enu_min, const double Enu_max,
                           const double Ee_min, const double Ee_max, const unsigned int nbins, const double classifier_cut);
void CreatePDFs_alphaN(TTree* EventInfo, std::vector<TH1D*>& PDF_hists, const double Ee_min, const double Ee_max,
                       const unsigned int nbins, const double classifier_cut);
void CreatePDFs_other(TTree* EventInfo, std::vector<TH1D*>& PDF_hists, const double Ee_min, const double Ee_max,
                      const unsigned int nbins, const double classifier_cut, std::string PDF_name);
std::vector<std::string> SplitString(std::string str);
TVector3 LLAtoECEF(const double& longitude, const double& latitude, const double& altitude);
double GetReactorDistanceLLA(const double& longitude, const double& latitude, const double &altitude);

int main(int argv, char** argc) {
    std::string reactor_events_address = argc[1];
    std::string alphaN_events_address = argc[2];
    std::string geoNuTh_events_address = argc[3];
    std::string geoNuU_events_address = argc[4];
    std::string accidentals_address = argc[5];
    std::string output_file = argc[6];
    double bin_width = std::stod(argc[7]);
    double classifier_cut = std::stod(argc[8]);

    // Read in files and get their TTrees
    TFile *reactorFile = TFile::Open(reactor_events_address.c_str());
    TFile *alphaNFile = TFile::Open(alphaN_events_address.c_str());
    TFile *geoNuThFile = TFile::Open(geoNuTh_events_address.c_str());
    TFile *geoNuUFile = TFile::Open(geoNuU_events_address.c_str());
    TFile *AccUFile = TFile::Open(accidentals_address.c_str());

    TTree *reactorEventTree = (TTree *) reactorFile->Get("prompt");
    TTree *alphaNEventTree = (TTree *) alphaNFile->Get("prompt");
    TTree *geoNuThEventTree = (TTree *) geoNuThFile->Get("prompt");
    TTree *geoNuUEventTree = (TTree *) geoNuUFile->Get("prompt");
    TTree *AccEventTree = (TTree *) AccUFile->Get("promptAccidental");

    // Create recon E_e vs true E_nu 2-D hist (MeV)
    // Keep binning the same, for consistency, even though cuts may change
    double Ee_min = 0.0;
    double Ee_max = 10.0;
    unsigned int N_bins_Ee = 200;
    unsigned int Nbins = (unsigned int) ((Ee_max - Ee_min) / bin_width);

    double Enu_min = ((neutron_mass_c2 + electron_mass_c2) * (neutron_mass_c2 + electron_mass_c2) - proton_mass_c2 * proton_mass_c2) / (2.0 * proton_mass_c2);  // minimum antinu energy for IBD
    double Enu_max = Ee_max + (neutron_mass_c2 - proton_mass_c2) + (5.0 / neutron_mass_c2); // convert from zeroth order Ee, then add 5/M to include 1st order (1/M) effects
    std::cout << "E_nu_min = " << Enu_min << ", E_nu_max = " << Enu_max << std::endl;

    TH2D* E_conv = new TH2D("E_conversion", "E_conversion", Nbins, Enu_min, Enu_max, Nbins, Ee_min, Ee_max);

    // Loop through and apply tagging + cuts
    std::vector<TH1D*> PDF_hists;
    TTree *outInfo = new TTree("reactor_vals", "Reactor IBD flux fractions and baselines");
    std::cout << "Looping through reactor IBD events..." << std::endl;
    CreatePDFs_reactorIBD(reactorEventTree, outInfo, PDF_hists, E_conv, Enu_min, Enu_max, Ee_min, Ee_max, Nbins, classifier_cut);
    std::cout << "Looping through alpha-n events..." << std::endl;
    CreatePDFs_alphaN(alphaNEventTree, PDF_hists, Ee_min, Ee_max, Nbins, classifier_cut);
    std::cout << "Looping through geo-nu Thorium IBD events..." << std::endl;
    CreatePDFs_other(geoNuThEventTree, PDF_hists, Ee_min, Ee_max, Nbins, classifier_cut, "geoNu_Th");
    std::cout << "Looping through geo-nu Uranium IBD events..." << std::endl;
    CreatePDFs_other(geoNuUEventTree, PDF_hists, Ee_min, Ee_max, Nbins, classifier_cut, "geoNu_U");
    std::cout << "Looping through Accidental events..." << std::endl;
    CreatePDFs_other(AccEventTree, PDF_hists, Ee_min, Ee_max, Nbins, classifier_cut, "Accidental");

    // Normalise 2D hist along y axis (prompt_E) for each x-bin (E_nu)
    double integ;
    for (unsigned int i = 1; i <= E_conv->GetXaxis()->GetNbins(); ++i) {
        integ = E_conv->Integral(i, i, 1, Nbins);
        if (integ != 0.0) {
            for (unsigned int j = 1; j <= E_conv->GetYaxis()->GetNbins(); ++j) {
                E_conv->SetBinContent(i, j, E_conv->GetBinContent(i, j) / integ);
            }
        }
    }

    // Write un-normalised PDFs to file (scale between reactors is important. Only normalise them after adding them together)
    TFile *outroot = new TFile(output_file.c_str(), "RECREATE");
    for (unsigned int i = 0; i < PDF_hists.size(); ++i) {
        PDF_hists.at(i)->Write();
    }
    E_conv->Write();
    outInfo->Write();
    outroot->Write();
    outroot->Close();

    return 0;
}

void CreatePDFs_reactorIBD(TTree* EventInfo, TTree* outInfo, std::vector<TH1D*>& PDF_hists, TH2D* E_conv, const double Enu_min, const double Enu_max,
                           const double Ee_min, const double Ee_max, const unsigned int nbins, const double classifier_cut) {


    // Set branch addresses to unpack TTree, and set up PDF histograms
    TString *originReactor = NULL;
    Double_t parentKE1, reconEnergy, classResult;

    EventInfo->SetBranchAddress("energy", &reconEnergy);  // corrected energy
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);
    EventInfo->SetBranchAddress("parentMeta1", &originReactor);
    EventInfo->SetBranchAddress("parentKE1", &parentKE1);

    // Saved histogram for each core (in both E_e and E_nu)
    std::vector<TH1D*> promptE_hists, Enu_hists;
    std::vector<double> reactor_fluxes;  // within prompt_E_cuts
    std::map<std::string, unsigned int> reactor_map; // {"core name", std::vector index}

    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();

    /* ~~~~~~~ loop through all entries and add to histograms, sorted by reactor/core ~~~~~~~ */

    unsigned int nentries = EventInfo->GetEntries();
    std::string core_name, hist_name;
    unsigned int core_idx;
    std::vector<std::string> splt_str;
    for (unsigned int a = 0; a < nentries; ++a) {
        if (a % 1000 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);

        if (!pass_classifier(reconEnergy, classResult, classifier_cut)) continue;

        core_name = originReactor->Data();  // core name
        splt_str = SplitString(core_name);
        if (splt_str.at(0) != "BRUCE") core_name = splt_str.at(0); // All cores in reactor complexes grouped up, except at Bruce
        
        if (reactor_map.find(core_name) == reactor_map.end()) {
            // reactor not added to list yet -> add it
            reactor_map.insert({core_name, promptE_hists.size()});
            promptE_hists.push_back(new TH1D((core_name + "_promptE").c_str(), (core_name + "_promptE").c_str(), nbins, Ee_min, Ee_max));
            Enu_hists.push_back(new TH1D((core_name + "_Enu").c_str(), (core_name + "_Enu").c_str(), nbins, Enu_min, Enu_max));
            reactor_fluxes.push_back(0);
        }

        // Add event to relevant hists
        core_idx = reactor_map.at(core_name);
        promptE_hists.at(core_idx)->Fill(reconEnergy);
        Enu_hists.at(core_idx)->Fill(parentKE1);
        E_conv->Fill(parentKE1, reconEnergy); // all events go to this hist

        // If it passes prompt E cuts, add to flux too (PDFs have broader prompt E cuts)
        if (reconEnergy > IBD_MIN_PROMPT_E && reconEnergy < IBD_MAX_PROMPT_E) reactor_fluxes.at(core_idx) += 1;
    }

    /* ~~~~~~~ Loop through all reactors/cores and sort into PDFs and vectors ~~~~~~~ */

    // crated hists and numbers that will be saved
    PDF_hists.push_back(new TH1D("PWR_promptE", "PWR_promptE", nbins, Ee_min, Ee_max));
    PDF_hists.push_back(new TH1D("PWR_Enu", "PWR_Enu", nbins, Enu_min, Enu_max));
    PDF_hists.push_back(new TH1D("PHWR_promptE", "PHWR_promptE", nbins, Ee_min, Ee_max));
    PDF_hists.push_back(new TH1D("PHWR_Enu", "PHWR_Enu", nbins, Enu_min, Enu_max));

    Double_t PWR_promptE_frac;
    std::vector<Double_t> PWR_Enu_fracs, PWR_Enu_baselines;
    std::vector<Double_t> PHWR_Enu_fracs, PHWR_Enu_baselines;
    double total_flux;

    // Link values to ttree
	outInfo->Branch("PWR_promptE_frac", &PWR_promptE_frac);
    outInfo->Branch("PWR_Enu_fracs", &PWR_Enu_fracs);
    outInfo->Branch("PWR_Enu_baselines", &PWR_Enu_baselines);
    outInfo->Branch("PHWR_Enu_fracs", &PHWR_Enu_fracs);
    outInfo->Branch("PHWR_Enu_baselines", &PHWR_Enu_baselines);

    // loop through core map
    double av_baseline;
    std::vector<std::string> reactor_types;
    std::vector<Double_t> fLatitude, fLongitute, fAltitude;
    int nCores, core_num;
    std::string reactor_name;
    for (auto& x : reactor_map) {  
        core_name = x.first;
        core_idx = x.second;

        total_flux += reactor_fluxes.at(core_idx);

        if (core_name.find(" ") == std::string::npos) {
            reactor_name = core_name; // If there are no spaces
        } else {
            splt_str = SplitString(core_name);
            if (splt_str.at(0) == "BRUCE") {  // Only bruce has its core number in the name
                reactor_name = splt_str.at(0);
                core_num = std::stoi(splt_str.at(1));
            }
            else {
                reactor_name = core_name;
            }
        }

        linkdb = db->GetLink("REACTOR", reactor_name);
        fLatitude = linkdb->GetDArray("latitude");
        fLongitute = linkdb->GetDArray("longitude");
        fAltitude = linkdb->GetDArray("altitude");

        if (reactor_name == "BRUCE") {
            // Get core baseline, then add info relevant hist/vectors
            PHWR_Enu_baselines.push_back(GetReactorDistanceLLA(fLongitute.at(core_num), fLatitude.at(core_num), fAltitude.at(core_num)));
            PHWR_Enu_fracs.push_back(reactor_fluxes.at(core_idx));
            PDF_hists.at(2)->Add(promptE_hists.at(core_idx));  // PHWR_promptE
            PDF_hists.at(3)->Add(Enu_hists.at(core_idx));  // PHWR_Enu
        } else {
            // Get average reactor complex baseline
            av_baseline = 0;
            nCores = linkdb->GetI("no_cores");
            for (unsigned int iCore = 0; iCore < nCores; ++iCore) {
                av_baseline += GetReactorDistanceLLA(fLongitute.at(iCore), fLatitude.at(iCore), fAltitude.at(iCore));
            }
            av_baseline /= (double)nCores;

            // Check if it is beyond averaging threshold
            if (av_baseline > MAX_BASELINE) {
                // Add to relevant number and (ignore hist)
                PWR_promptE_frac += reactor_fluxes.at(core_idx);
            } else {
                // Get reactor type (use first core)
                reactor_types = linkdb->GetSArray("core_spectrum");
                if (reactor_types.at(0) == "PHWR") {
                    std::cout << "PHWR: " << core_name << ", av_baseline: " << av_baseline << std::endl;
                    // add info relevant hist/vectors
                    PHWR_Enu_baselines.push_back(av_baseline);
                    PHWR_Enu_fracs.push_back(reactor_fluxes.at(core_idx));
                    PDF_hists.at(2)->Add(promptE_hists.at(core_idx));  // PHWR_promptE
                    PDF_hists.at(3)->Add(Enu_hists.at(core_idx));  // PHWR_Enu
                } else {
                    // add info relevant hist/vectors
                    PWR_Enu_baselines.push_back(av_baseline);
                    PWR_Enu_fracs.push_back(reactor_fluxes.at(core_idx));
                    PDF_hists.at(0)->Add(promptE_hists.at(core_idx));  // PWR_promptE
                    PDF_hists.at(1)->Add(Enu_hists.at(core_idx));  // PWR_Enu
                }
            }
        }
    }

    // Normalise fractions, then save to TTree
    PWR_promptE_frac /= total_flux;
    for (unsigned int i = 0; i < PWR_Enu_fracs.size(); ++i) {PWR_Enu_fracs.at(i) /= total_flux;}
    for (unsigned int i = 0; i < PHWR_Enu_fracs.size(); ++i) {PHWR_Enu_fracs.at(i) /= total_flux;}

    outInfo->Fill();
}

void CreatePDFs_alphaN(TTree* EventInfo, std::vector<TH1D*>& PDF_hists, const double Ee_min, const double Ee_max,
                       const unsigned int nbins, const double classifier_cut) {

    // Set branch addresses to unpack TTree, and set up PDF histograms
    Double_t reconEnergy, classResult;

    EventInfo->SetBranchAddress("energy", &reconEnergy);  // corrected energy
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);

    PDF_hists.push_back(new TH1D("alphaN_PR", "alphaN proton recoil", nbins, Ee_min, Ee_max));
    PDF_hists.push_back(new TH1D("alphaN_C12", "alphaN 12C scatter", nbins, Ee_min, Ee_max));
    PDF_hists.push_back(new TH1D("alphaN_O16", "alphaN 16O deexcitation", nbins, Ee_min, Ee_max));

    unsigned int alphaN_O16_idx = PDF_hists.size() - 1;
    unsigned int alphaN_C12_idx = alphaN_O16_idx - 1;
    unsigned int alphaN_PR_idx = alphaN_C12_idx - 1;

    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    unsigned int nentries = EventInfo->GetEntries();
    for (unsigned int a = 0; a < nentries; ++a) {
        if (a % 1000 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);

        if (!pass_classifier(reconEnergy, classResult, classifier_cut)) continue;

        if (reconEnergy < PROTON_RECOIL_E_MAX) {
            PDF_hists.at(alphaN_PR_idx)->Fill(reconEnergy);  // proton recoil
        } else if (reconEnergy < CARBON12_SCATTER_E_MAX) {
            PDF_hists.at(alphaN_C12_idx)->Fill(reconEnergy);  // 12C scatter
        } else {
            PDF_hists.at(alphaN_O16_idx)->Fill(reconEnergy);  // 16O deexcitation
        }
    }
}

void CreatePDFs_other(TTree* EventInfo, std::vector<TH1D*>& PDF_hists, const double Ee_min, const double Ee_max,
                      const unsigned int nbins, const double classifier_cut, std::string PDF_name) {

    // Set branch addresses to unpack TTree, and set up PDF histograms
    Double_t reconEnergy, classResult;

    EventInfo->SetBranchAddress("energy", &reconEnergy);  // corrected energy
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);

    PDF_hists.push_back(new TH1D(PDF_name.c_str(), PDF_name.c_str(), nbins, Ee_min, Ee_max));

    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    unsigned int nentries = EventInfo->GetEntries();
    for (unsigned int a = 0; a < nentries; ++a) {
        if (a % 1000 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);

        if (!pass_classifier(reconEnergy, classResult, classifier_cut)) continue;

        PDF_hists.at(PDF_hists.size()-1)->Fill(reconEnergy);
    }
}


/* ~~~~~~~~~~~ Reactor Tools ~~~~~~~~~~~ */

/**
 * @brief Split reactor name strings for getting baselines
 * 
 * @param str 
 * @return std::vector<std::string> = {"reactor name", "core number"}
 */
std::vector<std::string> SplitString(std::string str) {

    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end); //each word of string now in vector
    std::string dummy = "";
    std::vector<std::string> info;

    // If there is a core number, as normal (WARNING: if ractor name has a space, but there is no core number, this will return an empty vector)
    for (int i = 0; i < tokens.size(); i++) { //combine back to reactor name, core number
        if (i == 0) {
            dummy = tokens.at(i);
            if (i != tokens.size()-2) {
                dummy += " ";
            }
        }
        else if (i != tokens.size()-1) {
            dummy += tokens.at(i);
            if (i != tokens.size()-2) {
                dummy += " ";
            }
        }
        else {
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
TVector3 LLAtoECEF(const double& longitude, const double& latitude, const double& altitude) {
  // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
  double toRad = TMath::Pi()/180.;
  double Earthradius = 6378137.0; //Radius of the Earth (in meters)
  double f = 1./298.257223563; //Flattening factor WGS84 Model
  double L, rs, x, y, z;
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
double GetReactorDistanceLLA(const double& longitude, const double& latitude, const double &altitude) {
    const TVector3 SNO_ECEF_coord_ = TVector3(672.87,-4347.18,4600.51);
    double dist = (LLAtoECEF(longitude, latitude, altitude) - SNO_ECEF_coord_).Mag();
  return dist;
}
