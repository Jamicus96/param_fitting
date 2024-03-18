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


// #define USING_RUN_NUM
ULong64_t dcAnalysisWord = 36283883733698;  // Converted hex to decimal from 0x2100000042C2

// Define physical constants
double electron_mass_c2 = 0.510998910;  // MeV
double proton_mass_c2 = 938.272013;  // MeV
double neutron_mass_c2 = 939.56536;  // MeV

// define max delay, since it gets used twice (for consistency)
double MAX_DELAY = 0.8E6;
double PROTON_RECOIL_E_MAX = 3.5;  // (MeV)
double CARBON12_SCATTER_E_MAX = 5.4;  // (MeV)  Could just simulate different process separately?

bool pass_prompt_cuts(const double energy, const TVector3& position) {
    if (energy < 0.7) return false;  // min energy cut (MeV)
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


double EnergyCorrection(const double E, TVector3 pos, const bool is_data, RAT::DU::DetectorStateCorrection& stateCorr, RAT::DU::ReconCalibrator& e_cal) {
    // Data vs MC energy correction (Tony's)
    double Ecorr = e_cal.CalibrateEnergyRTF(is_data, E, std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()), pos.Z()); // gives the new E

    #ifdef USING_RUN_NUM
        // Correct for position coverage dependence (Logan's)
        RAT::DU::Point3D position(0, pos);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
        Ecorr /= stateCorr.GetCorrectionPos(position, 0, 0) / stateCorr.GetCorrection(9394, 0.75058); // a correction factor (divide E by it)
    #endif

    return Ecorr;
}

bool dcAppliedAndPassed(const bool is_data, ULong64_t dcApplied, ULong64_t dcFlagged) {
    if (!is_data) return true;
    if (!((dcApplied & dcAnalysisWord) == dcAnalysisWord)) return false;  // DC cut is not compatible
    if (!((dcFlagged & dcAnalysisWord) == dcFlagged)) return false;  // Did not pass DC cut
    return true;
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
void Apply_tagging_and_cuts(TTree* EventInfo, std::vector<TTree*>& CutTtrees, const double classiferCut, const bool is_data) {

    CutTtrees.push_back((TTree*)(EventInfo->CloneTree(0)));

    // Set branch addresses to unpack TTree
    Double_t reconEnergy, reconX, reconY, reconZ, classResult;
    ULong64_t eventTime, dcApplied, dcFlagged;
    Bool_t valid;
    Int_t mcIndex;

    EventInfo->SetBranchAddress("energy", &reconEnergy);
    EventInfo->SetBranchAddress("posx", &reconX);
    EventInfo->SetBranchAddress("posy", &reconY);
    EventInfo->SetBranchAddress("posz", &reconZ);
    EventInfo->SetBranchAddress("clockCount50", &eventTime);
    EventInfo->SetBranchAddress("fitValid", &valid);
    EventInfo->SetBranchAddress("mcIndex", &mcIndex);
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);
    EventInfo->SetBranchAddress("dcApplied", &dcApplied);
    EventInfo->SetBranchAddress("dcFlagged", &dcFlagged);

    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    TRandom3 *rndm_dummy = new TRandom3();
    Double_t rndm;

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    // Initialise DetectorStateCorrection (assume only one run in each file)
    RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
    RAT::DU::ReconCalibrator e_cal = RAT::DU::Utility::Get()->GetReconCalibrator();

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int nvaliddelayed = 0, nvalidpair = 0, nvalid = 0;
    double delayedTime;
    unsigned int delayedMCIndex;
    TVector3 delayedPos;
    TVector3 promptPos;
    double delay;
    bool passedDC;
    double Ecorr;
    int a = 0;
    // Loop through all events:
    // First, find prompt event that pass all cuts.
    // Then go through next 3 events to try to find a delayed event that also passes all cuts.
    while (a < nentries) {
        if (a % 1000 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);
        delayedPos = TVector3(reconX, reconY, reconZ);
        delayedTime = eventTime;
        delayedMCIndex = mcIndex;

        passedDC = dcAppliedAndPassed(is_data, dcApplied, dcFlagged);
        Ecorr = EnergyCorrection(reconEnergy, delayedPos, is_data, stateCorr, e_cal);

        if (valid and passedDC and pass_delayed_cuts(Ecorr, delayedPos)) {
            nvaliddelayed++;

            // Delayed event is valid, check through the previous 100 events for event that passes prompt + classifier + tagging cuts
            for (int b = 1; b <= 100; ++b) {
                if ((a - b) < 0) break;
                EventInfo->GetEntry(a - b);

                delay = (delayedTime - eventTime) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
                if (delay > MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

                promptPos = TVector3(reconX, reconY, reconZ);
                passedDC = dcAppliedAndPassed(is_data, dcApplied, dcFlagged);
                Ecorr = EnergyCorrection(reconEnergy, promptPos, is_data, stateCorr, e_cal);

                if (valid and passedDC and pass_prompt_cuts(Ecorr, promptPos) and pass_coincidence_cuts(delay, promptPos, delayedPos)) {
                    // Event pair survived analysis cuts
                    nvalidpair++;
                    if (pass_classifier(Ecorr, classResult, classiferCut)) {
                        // Event pair survived classifier cut
                        reconEnergy = Ecorr; // apply energy correction before re-writing it
                        std::cout << "nvalid = " << nvalid << ", CutTtrees.size() = " << CutTtrees.size() << std::endl;
                        if ((nvalid / 10000) > (CutTtrees.size()-1)) {
                            std::cout << "Adding TTree" << std::endl;
                            // Deal with ttrees getting too full to write to file
                            CutTtrees.push_back((TTree*)(EventInfo->CloneTree(0)));
                        }
                        std::cout << "CutTtrees.size() = " << CutTtrees.size() << std::endl;
                        CutTtrees.at(CutTtrees.size()-1)->Fill();
                        // EventInfo->GetEntry(a); // Also add delayed event
                        // EventInfoCut->Fill();
                        nvalid++;
                    }
                }
            }
        }
        a++;
    }
    std::cout << "From " << nentries << " entries, number of valid delayed events: " << nvaliddelayed
              << ", number of these event pairs surviving prompt + coincidence cuts: " << nvalidpair
              << ", number of these event pairs surviving classifier cut: " << nvalid << std::endl;
}


int main(int argv, char** argc) {
    std::string inputNtuple = argc[1];
    std::string outputNtuple = argc[2];
    double classifier_cut = std::stod(argc[3]);
    bool is_data = std::stoi(argc[4]);

    // Read in files and get their TTrees
    TFile *inFile = TFile::Open(inputNtuple.c_str());
    TTree *EventInfo = (TTree *) inFile->Get("output");

    // Define a new TTree to save events that pass cuts, and loop through and apply tagging + cuts
    std::cout << "Looping through events..." << std::endl;
    std::vector<TTree*> CutTtrees;
    Apply_tagging_and_cuts(EventInfo, CutTtrees, classifier_cut, is_data);

    // Write output ntuple to files, deal with ttrees getting too full to write to one file
    std::cout << "Writing TTrees..." << std::endl; 
    std::size_t SlashPos, DotPos;
    std::string filename, repoAddress, fileAddress;
    for (unsigned int i = 0; i < CutTtrees.size(); ++i) {
        // Make file name
        SlashPos = outputNtuple.find_last_of("/");
        filename = outputNtuple.substr(SlashPos+1, outputNtuple.length());
        repoAddress = outputNtuple.substr(0, SlashPos);
        DotPos = filename.find_first_of(".");
        fileAddress = repoAddress + filename.substr(0, DotPos) + "_" + std::to_string(i) + ".ntuple.root";
        TFile* outroot = new TFile(fileAddress.c_str(), "RECREATE");

        outroot->cd();
        CutTtrees.at(i)->Write();
        outroot->Close();
        delete outroot;
    }

    return 0;
}