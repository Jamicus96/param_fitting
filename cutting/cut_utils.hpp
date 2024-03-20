// header guard:
#ifndef cut_utils
#define cut_utils

// include
#include <iostream>
#include <TTree.h>
#include <TVector3.h>
#include <RAT/DB.hh>
#include <TRandom3.h>
#include <TMath.h>
#include <RAT/DU/Point3D.hh>
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>


#define VERBOSE

// #define USING_RUN_NUM
ULong64_t dcAnalysisWord = 36283883733698;  // Converted hex to decimal from 0x2100000042C2

// Define cut values
double MIN_PROMPT_E = 0.7, MAX_PROMPT_E = 9.0;  // [MeV]
double MIN_DELAYED_E = 1.85, MAX_DELAYED_E = 2.4;  // [MeV]
double FV_CUT = 5700;  // Max radius [mm]
double MIN_DELAY = 400, MAX_DELAY = 0.8E6;  // Delta t cut [ns]
double MAX_DIST = 1500;  // Delta R cut [mm]
double PROTON_RECOIL_E_MAX = 3.5;  // [MeV]

bool pass_prompt_cuts(const double energy, const TVector3& position) {
    if (energy < MIN_PROMPT_E) return false;  // min energy cut (MeV)
    if (energy > MAX_PROMPT_E) return false;  // max energy cut (MeV)
    if (position.Mag() > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(const double energy, const TVector3& position) {
    if (energy < MIN_DELAYED_E) return false;  // min energy cut (MeV)
    if (energy > MAX_DELAYED_E) return false;  // max energy cut (MeV)
    if (position.Mag() > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(const double delay, const TVector3& prompt_pos, const TVector3& delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < MIN_DELAY) return false;  // min delay cut (ns)
    if (delay > MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > MAX_DIST) return false;  // max distance cut (mm)

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


std::vector<TH1D*>& Apply_tagging_and_cuts(TTree* EventInfo, std::vector<TTree*>& CutTtrees, const double classiferCut, const bool is_data) {

    // Create some basic hists to save
    TH1D* hPromptE = new TH1D("promptE", "promptE", 200, 0, 10);
    TH1D* hDelayedE = new TH1D("delayedE", "delayedE", 200, 0, 10);
    TH1D* hPromptR = new TH1D("promptR", "promptR", 200, 0, FV_CUT);
    TH1D* hDelayedR = new TH1D("delayedR", "delayedR", 200, 0, FV_CUT);
    TH1D* hDeltaT = new TH1D("deltaT", "deltaT", 200, MIN_DELAY, MAX_DELAY);
    TH1D* hDeltaR = new TH1D("deltaR", "deltaR", 200, 0, MAX_DIST);

    // Make firsy TTree clone
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
    double promptEcorr, delayedEcorr;
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
        delayedEcorr = EnergyCorrection(reconEnergy, delayedPos, is_data, stateCorr, e_cal);

        if (valid and passedDC and pass_delayed_cuts(delayedEcorr, delayedPos)) {
            nvaliddelayed++;

            // Delayed event is valid, check through the previous 100 events for event that passes prompt + classifier + tagging cuts
            for (int b = 1; b <= 100; ++b) {
                if ((a - b) < 0) break;
                EventInfo->GetEntry(a - b);

                delay = (delayedTime - eventTime) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
                if (delay > MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

                promptPos = TVector3(reconX, reconY, reconZ);
                passedDC = dcAppliedAndPassed(is_data, dcApplied, dcFlagged);
                promptEcorr = EnergyCorrection(reconEnergy, promptPos, is_data, stateCorr, e_cal);

                if (valid and passedDC and pass_prompt_cuts(promptEcorr, promptPos) and pass_coincidence_cuts(delay, promptPos, delayedPos)) {
                    // Event pair survived analysis cuts
                    nvalidpair++;
                    if (pass_classifier(promptEcorr, classResult, classiferCut)) {
                        // Event pair survived classifier cut

                        // Add to output TTree
                        reconEnergy = promptEcorr; // apply energy correction before re-writing it
                        #ifdef VERBOSE
                            std::cout << "nvalid = " << nvalid << ", CutTtrees.size() = " << CutTtrees.size() << std::endl;
                        #endif
                        if ((nvalid / 10000) > (CutTtrees.size()-1)) {
                            #ifdef VERBOSE
                                std::cout << "Adding TTree" << std::endl;
                            #endif
                            // Deal with ttrees getting too full to write to file
                            CutTtrees.push_back((TTree*)(EventInfo->CloneTree(0)));
                        }
                        CutTtrees.at(CutTtrees.size()-1)->Fill();
                        
                        // Add info to hists
                        hPromptE->Fill(promptEcorr);
                        hDelayedE->Fill(delayedEcorr);
                        hPromptR->Fill(promptPos.Mag());
                        hDelayedR->Fill(delayedPos.Mag());
                        hDeltaT->Fill(delay);
                        hDeltaR->Fill((delayedPos - promptPos).Mag());

                        // Update
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
    
    // Package hists and return them
    std::vector<TH1D*> recordHists = {hPromptE, hDelayedE, hPromptR, hDelayedR, hDeltaT, hDeltaR};
    return recordHists;
}

//end header guard
#endif