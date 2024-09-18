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
#include "cutting_utils.hpp"

void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, const double classifier_cut, const bool is_data);


int main(int argv, char** argc) {
    std::string inputNtuple = argc[1];
    std::string previousRunNtuple = argc[2];
    std::string outputNtuple = argc[3];
    double classifier_cut = std::stod(argc[4]);
    bool is_data = std::stoi(argc[5]);

    std::cout << "Looping through events..." << std::endl;
    Apply_tagging_and_cuts(inputNtuple, previousRunNtuple, outputNtuple, classifier_cut, is_data);

    return 0;
}

void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, const double classifier_cut, const bool is_data) {

    // Read in file and create output file
    TFile* inFile = TFile::Open(inputNtuple.c_str());
    TFile* outroot = new TFile(outputNtuple.c_str(), "RECREATE");

    TTree* EventInfo = (TTree*) inFile->Get("output");
    TTree* CutPromptTree = EventInfo->CloneTree(0);
    CutPromptTree->SetObject("prompt", "prompt");
    TTree* CutDelayedTree = EventInfo->CloneTree(0);
    CutDelayedTree->SetObject("delayed", "delayed");

    // Set up db access
    RAT::DB *db = RAT::DB::Get();
    RAT::DS::Run run;

    #ifndef USING_RUN_NUM
        db->SetAirplaneModeStatus(true);
        db->LoadDefaults();
    #endif
    #ifdef USING_RUN_NUM
        db->LoadDefaults();
        db->SetServer("postgres://snoplus@pgsql.snopl.us:5400/ratdb");
        std::cout << "RAT DB tag: " << db->GetDBTag() << std::endl;
    #endif

    /* ~~~~~~~~~~~~~ If looking at end of previous ntuple (like previous run) ~~~~~~~~~~~~~ */

    Int_t lastRunNum = -999; //last run number loaded in DB
    int64_t highNhitTime = -99999999999;
    int64_t owlNhitTime = -99999999999;
    bool found_highNhit = false, found_owlNhit = false;
    
    // Read in file and find last high nhit event
    if (previousRunNtuple != "0" && previousRunNtuple != "false" && previousRunNtuple != "False") {
        TFile* previousFile = TFile::Open(previousRunNtuple.c_str());
        TTree* previousEventInfo = (TTree*) previousFile->Get("output");

        // Set branch addresses to unpack TTree
        ULong64_t prev_eventTime;
        Int_t prev_nHits, prev_owlnhits;
        Int_t prev_runNum; // run number associated with MC

        previousEventInfo->SetBranchAddress("clockCount50", &prev_eventTime);
        previousEventInfo->SetBranchAddress("runID", &prev_runNum);
        previousEventInfo->SetBranchAddress("nhits", &prev_nHits);
        previousEventInfo->SetBranchAddress("owlnhits", &prev_owlnhits);

        std::cout << "Check previous run for high nhit and owl nhit events..." << std::endl;
        int prev_nentries = previousEventInfo->GetEntries();
        for (int a = (prev_nentries-1); a >= 0; --a) {
            previousEventInfo->GetEntry(a);

            #ifdef USING_RUN_NUM
                if (lastRunNum != prev_runNum) {
                    run.SetRunID(prev_runNum);
                    db->BeginOfRun(run);
                    lastRunNum = prev_runNum;
                }
            #endif

            if (prev_nHits > 3000 && !found_highNhit) {
                highNhitTime = int64_t(prev_eventTime);
                found_highNhit = true;
                std::cout << "found!" << std::endl;
                break;
            }
            if (prev_owlnhits > 3 && !found_owlNhit) {
                owlNhitTime = int64_t(prev_eventTime);
                found_owlNhit = true;
            }
            if (found_highNhit && found_owlNhit) {
                std::cout << "found!" << std::endl;
                break;
            }
        }

    }

    /* ~~~~~~~~~~~~~ Look at main ntuple ~~~~~~~~~~~~~ */

    // Set branch addresses to unpack TTree
    Double_t reconEnergy, reconX, reconY, reconZ, classResult;
    ULong64_t eventTime, dcApplied, dcFlagged;
    Int_t nHits, GTID, owlnhits;
    Int_t runNum; //run number associated with MC
    Bool_t* valid;
    // Extra MC stuff, if available
    Double_t mcX, mcY, mcZ;
    Int_t EVindex;  // For MC: 0 = prompt, 1 = delayed
    Int_t pdg1, pdg2;

    EventInfo->SetBranchAddress("energy", &reconEnergy);
    EventInfo->SetBranchAddress("posx", &reconX);
    EventInfo->SetBranchAddress("posy", &reconY);
    EventInfo->SetBranchAddress("posz", &reconZ);
    EventInfo->SetBranchAddress("clockCount50", &eventTime);
    EventInfo->SetBranchAddress("fitValid", &valid);
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);
    EventInfo->SetBranchAddress("dcApplied", &dcApplied);
    EventInfo->SetBranchAddress("dcFlagged", &dcFlagged);
    EventInfo->SetBranchAddress("runID", &runNum);
    EventInfo->SetBranchAddress("nhits", &nHits);
    EventInfo->SetBranchAddress("owlnhits", &owlnhits);
    EventInfo->SetBranchAddress("eventID", &GTID);
    if (!is_data) {
        EventInfo->SetBranchAddress("mcPosx", &mcX);
        EventInfo->SetBranchAddress("mcPosy", &mcY);
        EventInfo->SetBranchAddress("mcPosz", &mcZ);
        EventInfo->SetBranchAddress("evIndex", &EVindex);
        EventInfo->SetBranchAddress("pdg1", &pdg1);
        EventInfo->SetBranchAddress("pdg2", &pdg2);
    }

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    // Initialise DetectorStateCorrection (assume only one run in each file)
    RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
    RAT::DU::ReconCalibrator* e_cal = RAT::DU::ReconCalibrator::Get();

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int nvaliddelayed = 0, nvalidpair = 0, nvalid = 0, nMuonCut = 0, nDCcut = 0, nValidCut = 0, negEcut = 0, nvalidMultiplicity = 0;
    unsigned int nEvtType = 0, nEvtType_insideFV_MC = 0;
    int64_t delayedTime, promptTime;
    TVector3 delayedPos, promptPos, multiplicityPos;
    double delay, highNhitDelay, owlNhitDelay;
    bool passMultiplicity;
    double promptEcorr, delayedEcorr;
    double MC_rPos;
    // Loop through all events:
    for (int a = 0; a < nentries; ++a) {
        if (a % 10000 == 0) std::cout << "Done " << ((float)a / (float)nentries) * 100.0 << "%" << std::endl;

        EventInfo->GetEntry(a);

        #ifdef USING_RUN_NUM
            if (lastRunNum != runNum) {
                run.SetRunID(runNum);
                db->BeginOfRun(run);
                lastRunNum = runNum;
            }
        #endif

        // Veto logic (accounted for in livetime)
        delayedTime = int64_t(eventTime);
        if (nHits > 3000) highNhitTime = delayedTime;
        if (owlnhits > 3) owlNhitTime = delayedTime;
        highNhitDelay = ((delayedTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
        owlNhitDelay = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
        if (highNhitDelay < highNhit_deltaT || owlNhitDelay < owlNhit_deltaT) {++nMuonCut; continue;}

        // Counting number of prompt IBD events simulated within FV, outside of veto windows
        if (!is_data) {
            if (EVindex == 0) {
                ++nEvtType;
                MC_rPos = sqrt(mcX*mcX + mcY*mcY + mcZ*mcZ);
                if (MC_rPos < FV_CUT && EVindex == 0) ++nEvtType_insideFV_MC;
            }
        }

        // Fit-valid and DC logic (accounted for in detector/cut efficiency)
        if (!valid) {++nValidCut; continue;}
        if (reconEnergy < 0) {++negEcut; continue;}
        if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

        // Compute quantities used for cuts
        delayedPos = TVector3(reconX, reconY, reconZ);
        delayedEcorr = EnergyCorrection(reconEnergy, delayedPos, is_data, stateCorr, e_cal);

        if (pass_delayed_cuts_IBD(delayedEcorr, delayedPos)) {
            nvaliddelayed++;

            // Delayed event is valid, check through the previous 1000 events for event that passes prompt + classifier + tagging cuts
            for (int b = 1; b <= 1000; ++b) {
                if ((a - b) < 0) break;
                EventInfo->GetEntry(a - b);

                promptTime = int64_t(eventTime);
                highNhitDelay = ((promptTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
                owlNhitDelay = ((promptTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
                if (highNhitDelay < highNhit_deltaT) {++nMuonCut; break;}
                if (fabs(owlNhitDelay) < owlNhit_deltaT) {++nMuonCut; continue;}
                if (!valid) {++nValidCut; continue;}
                if (reconEnergy < 0) {++negEcut; continue;}
                if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

                delay = ((delayedTime - promptTime) & 0x7FFFFFFFFFF) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
                if (delay > IBD_MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

                promptPos = TVector3(reconX, reconY, reconZ);
                promptEcorr = EnergyCorrection(reconEnergy, promptPos, is_data, stateCorr, e_cal);

                if (pass_prompt_cuts_IBD(promptEcorr, promptPos) and pass_coincidence_cuts_IBD(delay, promptPos, delayedPos)) {
                    // Event pair survived analysis cuts
                    nvalidpair++;
                    if (pass_classifier(promptEcorr, classResult, classifier_cut)) {
                        // Event pair survived classifier cut
                        nvalid++;

                        // check for multiplicity (any events around the event pairs that have E>0.4MeV and dr<2m)
                        passMultiplicity = true;
                        for (unsigned int c = 1; c < 1000; ++c) {
                            // before prompt
                            if ((a - b - c) < 0) break;
                            EventInfo->GetEntry(a - b - c);

                            delay = ((promptTime - int64_t(eventTime)) & 0x7FFFFFFFFFF) / 50E6 * 1E3; // convert number of ticks in 50MHz clock to ms
                            if (delay > 1.) break;

                            if (reconEnergy > 0.4) {
                                multiplicityPos = TVector3(reconX, reconY, reconZ);
                                if ((multiplicityPos - promptPos).Mag() < 2000. || (multiplicityPos - delayedPos).Mag() < 2000.) {
                                    passMultiplicity = false;
                                    break;
                                }
                            }
                        }
                        // after prompt is covered by before delayed
                        if (passMultiplicity) {
                            for (unsigned int c = 1; c < 1000; ++c) {
                                // before delayed
                                if (c >= b) break;
                                EventInfo->GetEntry(a - c);

                                delay = ((delayedTime - int64_t(eventTime)) & 0x7FFFFFFFFFF) / 50E6 * 1E3; // convert number of ticks in 50MHz clock to ms
                                if (delay > 1.) break;

                                if (reconEnergy > 0.4) {
                                    multiplicityPos = TVector3(reconX, reconY, reconZ);
                                    if ((multiplicityPos - promptPos).Mag() < 2000. || (multiplicityPos - delayedPos).Mag() < 2000.) {
                                        passMultiplicity = false;
                                        break;
                                    }
                                }
                            }
                        }
                        if (passMultiplicity) {
                            for (unsigned int c = 1; c < 1000; ++c) {
                                // after delayed
                                if ((a + c) > nentries) break;
                                EventInfo->GetEntry(a + c);

                                delay = ((int64_t(eventTime) - delayedTime) & 0x7FFFFFFFFFF) / 50E6 * 1E3; // convert number of ticks in 50MHz clock to ms
                                if (delay > 1.) break;

                                if (reconEnergy > 0.4) {
                                    multiplicityPos = TVector3(reconX, reconY, reconZ);
                                    if ((multiplicityPos - promptPos).Mag() < 2000. || (multiplicityPos - delayedPos).Mag() < 2000.) {
                                        passMultiplicity = false;
                                        break;
                                    }
                                }
                            }
                        }

                        if (passMultiplicity) {
                            // Add to output TTrees
                            EventInfo->GetEntry(a - b);
                            reconEnergy = promptEcorr;
                            CutPromptTree->Fill();
                            EventInfo->GetEntry(a);
                            reconEnergy = delayedEcorr;
                            CutDelayedTree->Fill();

                            ++nvalidMultiplicity;
                        }
                    }
                }
            }
        }
    }


    std::cout << "From " << nentries << " entries, number of valid delayed events: " << nvaliddelayed
              << ", number of these event pairs surviving prompt + coincidence cuts: " << nvalidpair
              << ", number of these event pairs surviving classifier cut: " << nvalid
              << ", number of that survive multiplicity cut: " << nvalidMultiplicity << std::endl;
    std::cout << "Events cut from Muon tagging: " << nMuonCut << ", invalid recon: " << nValidCut << ", negative E: " << negEcut << ", failed DC: " << nDCcut << std::endl;
    std::cout << "Last muon tag time (50MHz clock ticks): " << highNhitTime << std::endl;
    std::cout << "Number of prompt events simulated (outside veto windows): " << nEvtType << ", inside FV: " << nEvtType_insideFV_MC << std::endl;

    // Write output ntuple to files, deal with ttrees getting too full to write to one file
    std::cout << "Writing Trees..." << std::endl; 

    outroot->cd();
    CutPromptTree->Write();
    CutDelayedTree->Write();
    outroot->Close();
    delete outroot;
}
