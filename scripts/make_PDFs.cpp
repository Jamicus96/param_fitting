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
#include <map>

// Define physical constants
double electron_mass_c2 = 0.510998910;  // MeV
double proton_mass_c2 = 938.272013;  // MeV
double neutron_mass_c2 = 939.56536;  // MeV

// define max delay, since it gets used twice (for consistency)
double MAX_DELAY = 0.8E6;
double PROTON_RECOIL_E_MAX = 3.5;  // (MeV)
double CARBON12_SCATTER_E_MAX = 5.4;  // (MeV)  Could just simulate different process separately?

bool pass_prompt_cuts(const double energy, const TVector3& position) {
    if (energy < 0.5) return false;  // min energy cut (MeV)
    if (energy > 9.0) return false;  // max energy cut (MeV)
    if (position.Mag() > 5700) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(const double energy, const TVector3& position) {
    if (energy < 1.85) return false;  // min energy cut (MeV)
    if (energy > 2.4) return false;  // max energy cut (MeV)
    if (position.Mag() > 5700) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(const double delay, const TVector3& prompt_pos, const TVector3& delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < 400) return false;  // min delay cut (ns)
    if (delay > MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > 1500) return false;  // max distance cut (mm)

    return true;
}

bool pass_classifier(const double energy, const double class_result, const double class_cut) {
    // Only check classifier result for events in PDF1 (proton recoil, E < 3MeV)
    if (energy > PROTON_RECOIL_E_MAX) return true;
    if (class_result > class_cut) return true;

    return false;
}


/**
 * @brief Takes in an event tree, and compiles energy histograms.
 * Output is of the form of a map, so that for reactor IBD events, there is one
 * histogram for the events coming from each reactor core, where the key is the
 * name of the core. For other events, there is instead only one element in the
 * map, with key "alphaN". The keys are also the names of the histograms
 * themselves.
 * The point is that a PDF for each reactor can be created separately, so that
 * its baseline can be used correctly when oscillation is applied later.
 * 
 * @param EventInfo  event tree from ntuple.
 * @param classiferCut  events below this value are cut (only applied to events with E<3.5MeV so far).
 * @param is_reactorIBD  true if events are reactor IBDs, so that their provenance can be tracked.
 * @return std::map<std::string, TH1D*> 
 */
std::map<std::string, TH1D*> Apply_tagging_and_cuts(TTree* EventInfo, const double classiferCut, TH2D* E_conv, const double lowenergybin, const double maxenergybin,
                                                    const unsigned int nbins, const std::string data_type) {

    // Define a map of histograms {"hist name", hist}. This is to keep track of reactor IBD origin
    std::map<std::string, TH1D*> hists_map;

    // Set branch addresses to unpack TTree
    TString *originReactor = NULL;
    Double_t parentKE1;
    Double_t reconEnergy;
    Double_t reconX;
    Double_t reconY;
    Double_t reconZ;
    ULong64_t eventTime;
    Bool_t valid;
    Int_t mcIndex;
    Double_t classResult;

    EventInfo->SetBranchAddress("energy", &reconEnergy);
    EventInfo->SetBranchAddress("posx", &reconX);
    EventInfo->SetBranchAddress("posy", &reconY);
    EventInfo->SetBranchAddress("posz", &reconZ);
    EventInfo->SetBranchAddress("clockCount50", &eventTime);
    EventInfo->SetBranchAddress("fitValid", &valid);
    EventInfo->SetBranchAddress("mcIndex", &mcIndex);
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);
    if (data_type == "reactorIBD") {
        // Want to get incoming antinu's origin core and energy
        EventInfo->SetBranchAddress("parentMeta1", &originReactor);
        EventInfo->SetBranchAddress("parentKE1", &parentKE1);
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
    TRandom3 *rndm_dummy = new TRandom3();
    Double_t rndm;

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int nvaliddelayed = 0, nvalidpair = 0, nvalid = 0;
    double delayedTime;
    unsigned int delayedMCIndex;
    TVector3 delayedPos;
    TVector3 promptPos;
    double delay;
    unsigned int a = 0;

    // Loop through all events:
    // First, find prompt event that pass all cuts.
    // Then go through next 3 events to try to find a delayed event that also passes all cuts.
    while (a < nentries) {
        if (a % 100 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);
        const char* hist_name;
        if (data_type == "reactorIBD") {
            hist_name = originReactor->Data();  // core name
        } else if (data_type == "geoNu_Th" || data_type == "geoNu_U") {
            hist_name = data_type;
        } else {
            hist_name = "ERROR"; // To make any errors obvious, just in case (alpha-n names set later)
        }
        delayedPos = TVector3(reconX, reconY, reconZ);
        delayedTime = eventTime;
        delayedMCIndex = mcIndex;

        if (valid and pass_delayed_cuts(reconEnergy, delayedPos)) {
            nvaliddelayed++;

            // Delayed event is valid, check through the previous 100 events for event that passes prompt + classifier + tagging cuts
            for (unsigned int b = 1; b <= 100; ++b) {
                EventInfo->GetEntry(a - b);
                if (mcIndex != delayedMCIndex) continue;  // check just in case

                delay = (delayedTime - eventTime) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
                if (delay > MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

                promptPos = TVector3(reconX, reconY, reconZ);
                if (valid and pass_prompt_cuts(reconEnergy, promptPos) and pass_coincidence_cuts(delay, promptPos, delayedPos)) {
                    // Event pair survived analysis cuts
                    nvalidprompt++;
                    if (pass_classifier(reconEnergy, classResult, classiferCut)) {
                        // Event pair survived classifier cut
                        nvalid++;
                        if (data_type == "reactorIBD") {
                            if (hists_map.find(hist_name) == hists_map.end()) {
                                // reactor not added to list yet -> add it
                                TH1D* temp_hist = new TH1D(hist_name, hist_name, nbins, lowenergybin, maxenergybin);
                                hists_map.insert({hist_name, temp_hist});
                            }

                            // Add event to E_e vs E_nu 2D hist (for reactor IBDs)
                            E_conv->Fill(reconEnergy, parentKE1);

                        } else if (data_type == "alphaN") {
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
                    }
                }
            }
        }
        a++;
    }
    std::cout << "From " << nentries << " entries, number of valid delayed events: " << nvaliddelayed
              << ", number of these event pairs surviving prompt + coincidence cuts: " << nvalidpair
              << ", number of these event pairs surviving classifier cut: " << nvalid << std::endl;

    return hists_map;
}


int main(int argv, char** argc) {
    std::string reactor_events_address = argc[1];
    std::string alphaN_events_address = argc[2];
    std::string geoNuTh_events_address = argc[3];
    std::string geoNuU_events_address = argc[4];
    std::string output_file = argc[5];
    double classifier_cut = std::stod(argc[6]);

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
    std::cout << "Looping through reactor IBD events..." << std::endl; 
    std::map<std::string, TH1D*> reactor_hist_map = Apply_tagging_and_cuts(reactorEventTree, classifier_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "reactorIBD");
    std::cout << "Looping through alpha-n events..." << std::endl; 
    std::map<std::string, TH1D*> alphaN_hist_map = Apply_tagging_and_cuts(alphaNEventTree, classifier_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "alphaN");
    std::cout << "Looping through geo-nu Thorium IBD events..." << std::endl; 
    std::map<std::string, TH1D*> geoNuTh_hist_map = Apply_tagging_and_cuts(geoNuThEventTree, classifier_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "geoNu_Th");
    std::cout << "Looping through geo-nu Uranium IBD events..." << std::endl; 
    std::map<std::string, TH1D*> geoNuU_hist_map = Apply_tagging_and_cuts(geoNuUEventTree, classifier_cut, E_conv, Ee_min, Ee_max, N_bins_Ee, "geoNu_U");

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