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
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>
#include <map>
#include "cut_utils.hpp"


// Define physical constants
double electron_mass_c2 = 0.510998910;  // MeV
double proton_mass_c2 = 938.272013;  // MeV
double neutron_mass_c2 = 939.56536;  // MeV
double CARBON12_SCATTER_E_MAX = 5.4;  // (MeV)  Could just simulate different process separately?


std::map<std::string, TH1D*> CreatePDFs(std::vector<TTree*>& EventInfos, TH2D* E_conv, const double lowenergybin, const double maxenergybin,
                                        const unsigned int nbins, const std::string data_type) {

    // Define a map of histograms {"hist name", hist}. This is to keep track of reactor IBD origin
    std::map<std::string, TH1D*> hists_map;

    // Set branch addresses to unpack TTree
    TString *originReactor = NULL;
    Double_t parentKE1, reconEnergy, reconX, reconY, reconZ, classResult;
    ULong64_t eventTime, dcApplied, dcFlagged;
    Bool_t valid;
    Int_t mcIndex;

    for (unsigned int iTree = 0; iTree < EventInfos.size(); ++iTree) {
        EventInfos.at(iTree)->SetBranchAddress("energy", &reconEnergy);
        EventInfos.at(iTree)->SetBranchAddress("posx", &reconX);
        EventInfos.at(iTree)->SetBranchAddress("posy", &reconY);
        EventInfos.at(iTree)->SetBranchAddress("posz", &reconZ);
        EventInfos.at(iTree)->SetBranchAddress("clockCount50", &eventTime);
        EventInfos.at(iTree)->SetBranchAddress("fitValid", &valid);
        EventInfos.at(iTree)->SetBranchAddress("mcIndex", &mcIndex);
        EventInfos.at(iTree)->SetBranchAddress("alphaNReactorIBD", &classResult);
        EventInfos.at(iTree)->SetBranchAddress("dcApplied", &dcApplied);
        EventInfos.at(iTree)->SetBranchAddress("dcFlagged", &dcFlagged);
        if (data_type == "reactorIBD") {
            // Want to get incoming antinu's origin core and energy
            EventInfos.at(iTree)->SetBranchAddress("parentMeta1", &originReactor);
            EventInfos.at(iTree)->SetBranchAddress("parentKE1", &parentKE1);
        } else if (data_type == "alphaN") {
            // Instead, can set all the hists for alphaN, since there aren't many
            TH1D* temp_hist_1 = new TH1D("alphaN_PR", "alphaN proton recoil", nbins, lowenergybin, maxenergybin);
            TH1D* temp_hist_2 = new TH1D("alphaN_C12", "alphaN 12C scatter", nbins, lowenergybin, maxenergybin);
            TH1D* temp_hist_3 = new TH1D("alphaN_O16", "alphaN 16O deexcitation", nbins, lowenergybin, maxenergybin);
            hists_map.insert({"alphaN_PR", temp_hist_1});
            hists_map.insert({"alphaN_C12", temp_hist_2});
            hists_map.insert({"alphaN_O16", temp_hist_3});
        } else if (data_type == "geoNu_Th") {
            // Instead, can set the hist for geo-nus, since there's only one
            TH1D* temp_hist_geoNu_1 = new TH1D("geoNu_Th", "geo-nu Thorium", nbins, lowenergybin, maxenergybin);
            hists_map.insert({"geoNu_Th", temp_hist_geoNu_1});
        } else if (data_type == "geoNu_U") {
            TH1D* temp_hist_geoNu_2 = new TH1D("geoNu_U", "geo-nu Uranium", nbins, lowenergybin, maxenergybin);
            hists_map.insert({"geoNu_U", temp_hist_geoNu_2});
        } else {
            std::cout << "ERROR: data_type is wrong. Should be `reactorIBD`, `alphaN` or `geoNu`, not `" << data_type << "`." << std::endl;
            exit(1);
        }

        RAT::DBLinkPtr linkdb;
        RAT::DB *db = RAT::DB::Get();
        db->LoadDefaults();

        /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

        unsigned int nentries = EventInfos.at(iTree)->GetEntries();
        unsigned int a = 0;

        // Loop through all events:
        // First, find prompt event that pass all cuts.
        // Then go through next 3 events to try to find a delayed event that also passes all cuts.
        while (a < nentries) {
            if (a % 1000 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

            EventInfos.at(iTree)->GetEntry(a);
            const char* hist_name;
            if (data_type == "reactorIBD") {
                hist_name = originReactor->Data();  // core name

                if (hists_map.find(hist_name) == hists_map.end()) {
                    // reactor not added to list yet -> add it
                    TH1D* temp_hist = new TH1D(hist_name, hist_name, nbins, lowenergybin, maxenergybin);
                    hists_map.insert({hist_name, temp_hist});
                }

                // Add event to E_e vs E_nu 2D hist (for reactor IBDs)
                E_conv->Fill(reconEnergy, parentKE1);

            } else if (data_type == "geoNu_Th" || data_type == "geoNu_U") {
                hist_name = data_type.c_str();
            } else {
                hist_name = "ERROR"; // To make any errors obvious, just in case (alpha-n names set later)

                if (reconEnergy < PROTON_RECOIL_E_MAX) {
                    hist_name = "alphaN_PR";  // proton recoil
                } else if (reconEnergy < CARBON12_SCATTER_E_MAX) {
                    hist_name = "alphaN_C12";  // 12C scatter
                } else {
                    hist_name = "alphaN_O16";  // 16O deexcitation
                }
            }

            // Add event to relevant histogram
            hists_map.at(hist_name)->Fill(reconEnergy);

            a++;
        }
    }

    return hists_map;
}


int main(int argv, char** argc) {
    std::string reactor_events_address = argc[1];
    std::string alphaN_events_address = argc[2];
    std::string geoNuTh_events_address = argc[3];
    std::string geoNuU_events_address = argc[4];
    std::string output_file = argc[5];
    double classifier_cut = std::stod(argc[6]);
    bool is_data = std::stoi(argc[7]);

    // Read in files and get their TTrees
    TFile *reactorFile = TFile::Open(reactor_events_address.c_str());
    TFile *alphaNFile = TFile::Open(alphaN_events_address.c_str());
    TFile *geoNuThFile = TFile::Open(geoNuTh_events_address.c_str());
    TFile *geoNuUFile = TFile::Open(geoNuU_events_address.c_str());

    TTree *reactorEventTree = (TTree *) reactorFile->Get("output");
    TTree *alphaNEventTree = (TTree *) alphaNFile->Get("output");
    TTree *geoNuThEventTree = (TTree *) geoNuThFile->Get("output");
    TTree *geoNuUEventTree = (TTree *) geoNuUFile->Get("output");

    // Create recon E_e vs true E_nu 2-D hist (MeV)
    // Keep binning the same, for consistency, even though cuts may change
    double Ee_min = 0.0;
    double Ee_max = 10.0;
    unsigned int N_bins_Ee = 100;

    double Enu_min = ((neutron_mass_c2 + electron_mass_c2) * (neutron_mass_c2 + electron_mass_c2) - proton_mass_c2 * proton_mass_c2) / (2.0 * proton_mass_c2);  // minimum antinu energy for IBD
    double Enu_max = Ee_max + (neutron_mass_c2 - proton_mass_c2) + (5.0 / neutron_mass_c2); // convert from zeroth order Ee, then add 5/M to include 1st order (1/M) effects
    unsigned int N_bins_Enu = 100;
    std::cout << "E_nu_min = " << Enu_min << ", E_nu_max = " << Enu_max << std::endl;

    TH2D* E_conv = new TH2D("E_conversion", "E_conversion", N_bins_Ee, Ee_min, Ee_max, N_bins_Enu, Enu_min, Enu_max);

    // Loop through and apply tagging + cuts
    std::vector<TH1D*> recordHists;
    std::cout << "Looping through reactor IBD events..." << std::endl;
    std::vector<TTree*> reactorEventTrees_cut;
    recordHists = Apply_tagging_and_cuts(reactorEventTree, reactorEventTrees_cut, classifier_cut, is_data);
    std::map<std::string, TH1D*> reactor_hist_map = CreatePDFs(reactorEventTrees_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "reactorIBD");

    std::cout << "Looping through alpha-n events..." << std::endl;
    std::vector<TTree*> alphaNEventTrees_cut;
    recordHists = Apply_tagging_and_cuts(alphaNEventTree, alphaNEventTrees_cut, classifier_cut, is_data);
    std::map<std::string, TH1D*> alphaN_hist_map = CreatePDFs(alphaNEventTrees_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "alphaN");

    std::cout << "Looping through geo-nu Thorium IBD events..." << std::endl;
    std::vector<TTree*> geoNuThEventTrees_cut;
    recordHists = Apply_tagging_and_cuts(geoNuThEventTree, geoNuThEventTrees_cut, classifier_cut, is_data);
    std::map<std::string, TH1D*> geoNuTh_hist_map = CreatePDFs(geoNuThEventTrees_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "geoNu_Th");

    std::cout << "Looping through geo-nu Uranium IBD events..." << std::endl;
    std::vector<TTree*> geoNuUEventTrees_cut;
    recordHists = Apply_tagging_and_cuts(geoNuUEventTree, geoNuUEventTrees_cut, classifier_cut, is_data);
    std::map<std::string, TH1D*> geoNuU_hist_map = CreatePDFs(geoNuUEventTrees_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "geoNu_U");

    // Normalise 2D hist along y axis (E_nu) for each x-bin (E_e)
    double integ;
    for (unsigned int i = 1; i <= E_conv->GetXaxis()->GetNbins(); ++i) {
        integ = E_conv->Integral(i, i, 1, N_bins_Enu);
        for (unsigned int j = 1; j <= E_conv->GetYaxis()->GetNbins(); ++j) {
            if (integ != 0.0) E_conv->SetBinContent(i, j, E_conv->GetBinContent(i, j) / integ);
        }
    }

    // Write un-normalised PDFs to file (scale between reactors is important. Only normalise them after adding them together)
    TFile *outroot = new TFile(output_file.c_str(), "RECREATE");
    for (auto& x : reactor_hist_map) {  
        x.second->Write();  // loops through all the map entries (.first = key, .second = element), and writes them to the root file
    }
    for (auto& x : alphaN_hist_map) {  
        x.second->Write();
    }
    for (auto& x : geoNuTh_hist_map) {  
        x.second->Write();
    }
    for (auto& x : geoNuU_hist_map) {  
        x.second->Write();
    }
    E_conv->Write();
    outroot->Write();
    outroot->Close();

    return 0;
}