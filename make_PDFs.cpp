#include <TFile.h>
#include <TH1.h>
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


bool pass_prompt_cuts(double energy, TVector3 position) {
    if (energy < 0.9) return false;  // min energy cut (MeV)
    if (energy > 8.0) return false;  // max energy cut (MeV)
    if (position.Mag() > 5700) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(double energy, TVector3 position) {
    if (energy < 1.7) return false;  // min energy cut (MeV)
    if (energy > 2.5) return false;  // max energy cut (MeV)
    if (position.Mag() > 5700) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(double prompt_time, double delayed_time, TVector3 prompt_pos, TVector3 delayed_pos) {
    double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < 200) return false;  // min delay cut (ns)
    if (delay >  1.1E6) return false;  // max delay cut (ns)
    if (distance > 1500) return false;  // max distance cut (mm)

    return true;
}

bool pass_classifier(double energy, double class_result, double class_cut) {
    // Only check classifier result for events in PDF1 (proton recoil, E < 3MeV)
    if (energy > 3.5) return true;
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
std::map<std::string, TH1D*> Apply_tagging_and_cuts(TTree *EventInfo, double classiferCut, bool is_reactorIBD) {

    // Define binning params
    int nbins = 100;
    double lowenergybin = 0.9;  // MeV
    double maxenergybin = 8.0;  // MeV

    // Define a map of histograms {"hist name", hist}. This is to keep track of reactor IBD origin
    std::map<std::string, TH1D*> hists_map;

    // Set branch addresses to unpack TTree
    TString *originReactor = NULL;
    // Double_t parentKE1;
    Double_t reconEnergy;
    Double_t reconX;
    Double_t reconY;
    Double_t reconZ;
    ULong64_t eventTime;
    Bool_t valid;
    Int_t mcIndex;
    Double_t classResult;

    // EventInfo->SetBranchAddress("parentKE1", &parentKE1);
    EventInfo->SetBranchAddress("energy", &reconEnergy);
    EventInfo->SetBranchAddress("posx", &reconX);
    EventInfo->SetBranchAddress("posy", &reconY);
    EventInfo->SetBranchAddress("posz", &reconZ);
    EventInfo->SetBranchAddress("clockCount50", &eventTime);
    EventInfo->SetBranchAddress("fitValid", &valid);
    EventInfo->SetBranchAddress("mcIndex", &mcIndex);
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);
    if (is_reactorIBD) {
        EventInfo->SetBranchAddress("parentMeta1", &originReactor);
    }

    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    TRandom3 *rndm_dummy = new TRandom3();
    Double_t rndm;

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int noscillatedAnalysis = 0;
    unsigned int nvalidprompt = 0;
    TVector3 promptPos;
    double promptEnergy;
    double promptTime;
    unsigned int promptMCIndex;
    TVector3 delayedPos;
    unsigned int a = 0;

    // Loop through all events:
    // First, find prompt event that pass all cuts.
    // Then go through next 3 events to try to find a delayed event that also passes all cuts.
    while (a < nentries) {
        if (a % 100 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);
        const char* core_name;
        if (is_reactorIBD) {
            core_name = originReactor->Data();
        } else {
            core_name = "alphaN";
        }
        promptEnergy = reconEnergy;
        promptPos = TVector3(reconX, reconY, reconZ);
        promptTime = eventTime;
        promptMCIndex = mcIndex;

        if (valid and pass_prompt_cuts(promptEnergy, promptPos) and pass_classifier(promptEnergy, classResult, classiferCut)) {
            nvalidprompt++;
            a++;
            EventInfo->GetEntry(a);
            if (mcIndex != promptMCIndex) continue;

            delayedPos = TVector3(reconX, reconY, reconZ);
            if (valid and pass_delayed_cuts(reconEnergy, delayedPos) and pass_coincidence_cuts(promptTime, eventTime, promptPos, delayedPos)) {
                // Event pair survived analysis cuts
                noscillatedAnalysis++;
                if (hists_map.find(core_name) == hists_map.end()) {
                    // reactor not added to list yet -> add it
                    TH1D* temp_hist = new TH1D(core_name, core_name, nbins, lowenergybin, maxenergybin);
                    hists_map.insert({core_name, temp_hist});
                }
                // Add event to relevant histogram
                hists_map.at(core_name)->Fill(promptEnergy);

            } else {
                // A different delayed event? try next event
                a++;
                EventInfo->GetEntry(a);

                delayedPos = TVector3(reconX, reconY, reconZ);
                if (valid and pass_delayed_cuts(reconEnergy, delayedPos) and pass_coincidence_cuts(promptTime, eventTime, promptPos, delayedPos)) {
                    // Event pair survived analysis cuts
                    noscillatedAnalysis++;
                    if (hists_map.find(core_name) == hists_map.end()) {
                        // reactor not added to list yet -> add it
                        TH1D* temp_hist = new TH1D(core_name, core_name, nbins, lowenergybin, maxenergybin);
                        hists_map.insert({core_name, temp_hist});
                    }
                    // Add event to relevant histogram
                    hists_map.at(core_name)->Fill(promptEnergy);

                } else {
                    // A different delayed event? try next event
                    a++;
                    EventInfo->GetEntry(a);

                    delayedPos = TVector3(reconX, reconY, reconZ);
                    if(valid and pass_delayed_cuts(reconEnergy, delayedPos) and pass_coincidence_cuts(promptTime, eventTime, promptPos, delayedPos)){
                        // Event pair survived analysis cuts
                        noscillatedAnalysis++;
                        if (hists_map.find(core_name) == hists_map.end()) {
                            // reactor not added to list yet -> add it
                            TH1D* temp_hist = new TH1D(core_name, core_name, nbins, lowenergybin, maxenergybin);
                            hists_map.insert({core_name, temp_hist});
                        }
                        // Add event to relevant histogram
                        hists_map.at(core_name)->Fill(promptEnergy);
                    }
                }
            }
        }
        a++;
    }
    std::cout  << "Number of valid prompt events: "<< nvalidprompt << " number of entries surviving analysis cuts: " << noscillatedAnalysis << std::endl;

    return hists_map;
}


int main(int argv, char** argc) {
    std::string reactor_events_address = argc[1];
    std::string alphaN_events_address = argc[2];
    // std::string geoNu_events_address = argc[3];
    std::string output_file = argc[3];
    double classifier_cut = std::stod(argc[4]);

    // Read in files and get their TTrees
    TFile *reactorFile = TFile::Open(reactor_events_address.c_str());
    TFile *alphaNFile = TFile::Open(alphaN_events_address.c_str());

    TTree *reactorEventTree = (TTree *) reactorFile->Get("output");
    TTree *alphaNEventTree = (TTree *) alphaNFile->Get("output");

    // Loop through and apply tagging + cuts
    std::cout << "Looping through reactor IBD events..." << std::endl; 
    std::map<std::string, TH1D*> reactor_hist_map = Apply_tagging_and_cuts(reactorEventTree, classifier_cut, true);
    std::cout << "Looping through alpha-n events..." << std::endl; 
    std::map<std::string, TH1D*> alphaN_hist_map = Apply_tagging_and_cuts(alphaNEventTree, classifier_cut, false);

    // Write un-normalised PDFs to file (scale between reactors is important. Only normalise them after adding them together)
    TFile *outroot = new TFile(output_file.c_str(),"RECREATE");
    for (auto& x : reactor_hist_map) {  
        x.second->Write();  // loops through all the map entries (.first = key, .second = element), and writes them to the root file
    }
    for (auto& x : alphaN_hist_map) {  
        x.second->Write();  // loops through all the map entries (.first = key, .second = element), and writes them to the root file
    }
    outroot->Write();
    outroot->Close();

    return 0;
}