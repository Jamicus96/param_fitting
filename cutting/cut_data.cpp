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
    if (is_data && previousRunNtuple != "0" && previousRunNtuple != "false" && previousRunNtuple != "False") {
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
    Bool_t valid;
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

    int nentries = EventInfo->GetEntries();
    int nMuonCut = 0, nValidCut = 0, negEcut = 0, nDCcut = 0;
    // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
    std::vector<std::string> passNames = {"total", "valid", "FV", "Prompt", "Delayed", "Dt", "DR", "Class", "Mult"};
    std::vector<int> dataCounter = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindexNeg = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindex0 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindex1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindex2 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindexAbove2 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindexNeg_MCFV = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindex0_MCFV = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindex1_MCFV = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindex2_MCFV = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> EVindexAbove2_MCFV = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    int64_t delayedTime, promptTime;
    TVector3 delayedPos, promptPos, multiplicityPos;
    double delay, highNhitDelay, owlNhitDelay;
    bool passMultiplicity;
    double promptEcorr, delayedEcorr;
    double delayed_MC_R, prompt_MC_R;
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
        if (is_data) {
            if (nHits > 3000) highNhitTime = delayedTime;
            if (owlnhits > 3) owlNhitTime = delayedTime;
            highNhitDelay = ((delayedTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
            owlNhitDelay = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
            if (highNhitDelay < highNhit_deltaT || owlNhitDelay < owlNhit_deltaT) {++nMuonCut; continue;}
        }

        dataCounter.at(0)++;
        if (!is_data) {
            // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
            if (EVindex < 0) EVindexNeg.at(0)++;
            if (EVindex == 0) EVindex0.at(0)++;
            if (EVindex == 1) EVindex1.at(0)++;
            if (EVindex == 2) EVindex2.at(0)++;
            if (EVindex > 2) EVindexAbove2.at(0)++;
            delayed_MC_R = sqrt(mcX*mcX + mcY*mcY + (mcZ - AV_offset)*(mcZ - AV_offset));
            if (delayed_MC_R < FV_CUT) {
                if (EVindex < 0) EVindexNeg_MCFV.at(0)++;
                if (EVindex == 0) EVindex0_MCFV.at(0)++;
                if (EVindex == 1) EVindex1_MCFV.at(0)++;
                if (EVindex == 2) EVindex2_MCFV.at(0)++;
                if (EVindex > 2) EVindexAbove2_MCFV.at(0)++;
            }
        }

        // Fit-valid and DC logic (accounted for in detector/cut efficiency)
        if (!valid) {++nValidCut; continue;}
        if (reconEnergy < 0) {++negEcut; continue;}
        if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

        dataCounter.at(1)++;
        if (!is_data) {
            // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
            if (EVindex < 0) EVindexNeg.at(1)++;
            if (EVindex == 0) EVindex0.at(1)++;
            if (EVindex == 1) EVindex1.at(1)++;
            if (EVindex == 2) EVindex2.at(1)++;
            if (EVindex > 2) EVindexAbove2.at(1)++;
            if (delayed_MC_R < FV_CUT) {
                if (EVindex < 0) EVindexNeg_MCFV.at(1)++;
                if (EVindex == 0) EVindex0_MCFV.at(1)++;
                if (EVindex == 1) EVindex1_MCFV.at(1)++;
                if (EVindex == 2) EVindex2_MCFV.at(1)++;
                if (EVindex > 2) EVindexAbove2_MCFV.at(1)++;
            }
        }

        // Compute quantities used for cuts
        delayedPos = TVector3(reconX, reconY, reconZ);
        if (!pass_FV_cut(delayedPos)) continue;

        dataCounter.at(2)++;
        if (!is_data) {
            // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
            if (EVindex < 0) EVindexNeg.at(2)++;
            if (EVindex == 0) EVindex0.at(2)++;
            if (EVindex == 1) EVindex1.at(2)++;
            if (EVindex == 2) EVindex2.at(2)++;
            if (EVindex > 2) EVindexAbove2.at(2)++;
            if (delayed_MC_R < FV_CUT) {
                if (EVindex < 0) EVindexNeg_MCFV.at(2)++;
                if (EVindex == 0) EVindex0_MCFV.at(2)++;
                if (EVindex == 1) EVindex1_MCFV.at(2)++;
                if (EVindex == 2) EVindex2_MCFV.at(2)++;
                if (EVindex > 2) EVindexAbove2_MCFV.at(2)++;
            }
        }

        delayedEcorr = EnergyCorrection(reconEnergy, delayedPos, is_data, stateCorr, e_cal);
        if (!pass_promptE_cuts_IBD(delayedEcorr)) continue;

        dataCounter.at(3)++;
        if (!is_data) {
            // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
            if (EVindex < 0) EVindexNeg.at(3)++;
            if (EVindex == 0) EVindex0.at(3)++;
            if (EVindex == 1) EVindex1.at(3)++;
            if (EVindex == 2) EVindex2.at(3)++;
            if (EVindex > 2) EVindexAbove2.at(3)++;
            if (delayed_MC_R < FV_CUT) {
                if (EVindex < 0) EVindexNeg_MCFV.at(3)++;
                if (EVindex == 0) EVindex0_MCFV.at(3)++;
                if (EVindex == 1) EVindex1_MCFV.at(3)++;
                if (EVindex == 2) EVindex2_MCFV.at(3)++;
                if (EVindex > 2) EVindexAbove2_MCFV.at(3)++;
            }
        }

        if (!pass_delayedE_cuts_IBD(delayedEcorr)) continue;

        dataCounter.at(4)++;
        if (!is_data) {
            // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
            if (EVindex < 0) EVindexNeg.at(4)++;
            if (EVindex == 0) EVindex0.at(4)++;
            if (EVindex == 1) EVindex1.at(4)++;
            if (EVindex == 2) EVindex2.at(4)++;
            if (EVindex > 2) EVindexAbove2.at(4)++;
            if (delayed_MC_R < FV_CUT) {
                if (EVindex < 0) EVindexNeg_MCFV.at(4)++;
                if (EVindex == 0) EVindex0_MCFV.at(4)++;
                if (EVindex == 1) EVindex1_MCFV.at(4)++;
                if (EVindex == 2) EVindex2_MCFV.at(4)++;
                if (EVindex > 2) EVindexAbove2_MCFV.at(4)++;
            }
        }

        // Delayed event is valid, check through the previous 1000 events for event that passes prompt + classifier + tagging cuts
        for (int b = 1; b <= 1000; ++b) {
            if ((a - b) < 0) break;
            EventInfo->GetEntry(a - b);

            promptTime = int64_t(eventTime);
            if (is_data) {
                highNhitDelay = ((promptTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
                owlNhitDelay = ((promptTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
                if (highNhitDelay < highNhit_deltaT) {++nMuonCut; break;}
                if (fabs(owlNhitDelay) < owlNhit_deltaT) {++nMuonCut; continue;}
            }
            if (!valid) {++nValidCut; continue;}
            if (reconEnergy < 0) {++negEcut; continue;}
            if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

            delay = ((delayedTime - promptTime) & 0x7FFFFFFFFFF) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
            if (delay > IBD_MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

            promptPos = TVector3(reconX, reconY, reconZ);
            promptEcorr = EnergyCorrection(reconEnergy, promptPos, is_data, stateCorr, e_cal);

            if (!pass_prompt_cuts_IBD(promptEcorr, promptPos)) continue;
            if (!pass_dT_cut_IBD(delay)) continue;

            dataCounter.at(5)++;
            if (!is_data) {
                // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
                if (EVindex < 0) EVindexNeg.at(5)++;
                if (EVindex == 0) EVindex0.at(5)++;
                if (EVindex == 1) EVindex1.at(5)++;
                if (EVindex == 2) EVindex2.at(5)++;
                if (EVindex > 2) EVindexAbove2.at(5)++;
                prompt_MC_R = sqrt(mcX*mcX + mcY*mcY + (mcZ - AV_offset)*(mcZ - AV_offset));
                if ((delayed_MC_R < FV_CUT) && (prompt_MC_R < FV_CUT)) {
                    if (EVindex < 0) EVindexNeg_MCFV.at(5)++;
                    if (EVindex == 0) EVindex0_MCFV.at(5)++;
                    if (EVindex == 1) EVindex1_MCFV.at(5)++;
                    if (EVindex == 2) EVindex2_MCFV.at(5)++;
                    if (EVindex > 2) EVindexAbove2_MCFV.at(5)++;
                }
            }

            if (!pass_dR_cut_IBD(promptPos, delayedPos)) continue;

            dataCounter.at(6)++;
            if (!is_data) {
                // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
                if (EVindex < 0) EVindexNeg.at(6)++;
                if (EVindex == 0) EVindex0.at(6)++;
                if (EVindex == 1) EVindex1.at(6)++;
                if (EVindex == 2) EVindex2.at(6)++;
                if (EVindex > 2) EVindexAbove2.at(6)++;
                if ((delayed_MC_R < FV_CUT) && (prompt_MC_R < FV_CUT)) {
                    if (EVindex < 0) EVindexNeg_MCFV.at(6)++;
                    if (EVindex == 0) EVindex0_MCFV.at(6)++;
                    if (EVindex == 1) EVindex1_MCFV.at(6)++;
                    if (EVindex == 2) EVindex2_MCFV.at(6)++;
                    if (EVindex > 2) EVindexAbove2_MCFV.at(6)++;
                }
            }

            if (!pass_classifier(promptEcorr, classResult, classifier_cut)) continue;

            dataCounter.at(7)++;
            if (!is_data) {
                // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
                if (EVindex < 0) EVindexNeg.at(7)++;
                if (EVindex == 0) EVindex0.at(7)++;
                if (EVindex == 1) EVindex1.at(7)++;
                if (EVindex == 2) EVindex2.at(7)++;
                if (EVindex > 2) EVindexAbove2.at(7)++;
                if ((delayed_MC_R < FV_CUT) && (prompt_MC_R < FV_CUT)) {
                    if (EVindex < 0) EVindexNeg_MCFV.at(7)++;
                    if (EVindex == 0) EVindex0_MCFV.at(7)++;
                    if (EVindex == 2) EVindex2_MCFV.at(7)++;
                    if (EVindex > 2) EVindexAbove2_MCFV.at(7)++;
                }
            }

            // check for multiplicity (any events around the event pairs that have E>0.4MeV and dr<2m)
            passMultiplicity = true;
            for (int c = 1; c < 1000; ++c) {
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
                for (int c = 1; c < 1000; ++c) {
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
                for (int c = 1; c < 1000; ++c) {
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

                dataCounter.at(8)++;
                if (!is_data) {
                    // total, valid, passFV, passPrompt, passDelayed, passDt, passDR, passClass, passMult
                    if (EVindex < 0) EVindexNeg.at(8)++;
                    if (EVindex == 0) EVindex0.at(8)++;
                    if (EVindex == 1) EVindex1.at(8)++;
                    if (EVindex == 2) EVindex2.at(8)++;
                    if (EVindex > 2) EVindexAbove2.at(8)++;
                    if ((delayed_MC_R < FV_CUT) && (prompt_MC_R < FV_CUT)) {
                        if (EVindex < 0) EVindexNeg_MCFV.at(8)++;
                        if (EVindex == 0) EVindex0_MCFV.at(8)++;
                        if (EVindex == 1) EVindex1_MCFV.at(8)++;
                        if (EVindex == 2) EVindex2_MCFV.at(8)++;
                        if (EVindex > 2) EVindexAbove2_MCFV.at(8)++;
                    }
                }
            }
        }
    }


    std::cout << "From " << nentries << " entries, nMuonCut: " << nMuonCut << ", nValidCut: " << nValidCut << ", negEcut: " << negEcut << ", nDCcut: " << nDCcut << std::endl;
    std::cout << "Passed All EVindex<0 EVindex==0 EVindex==1 EVindex==2 EVindex>2 EVindex<0&&rMC<5.7 EVindex==0&&rMC<5.7 EVindex==1&&rMC<5.7 EVindex==2&&rMC<5.7 EVindex>2&&rMC<5.7" << std::endl;
    for (int i = 0; i < passNames.size(); ++i) {
        std::cout << passNames.at(i) << " " << dataCounter.at(i) << " "
                  << EVindexNeg.at(i) << " " << EVindex0.at(i) << " " << EVindex1.at(i) << " " << EVindex2.at(i) << " " << EVindexAbove2.at(i) << " "
                  << EVindexNeg_MCFV.at(i) << " " << EVindex0_MCFV.at(i) << " " << EVindex1_MCFV.at(i) << " " << EVindex2_MCFV.at(i) << " " << EVindexAbove2_MCFV.at(i) << std::endl;
    }

    // Write output ntuple to files, deal with ttrees getting too full to write to one file
    std::cout << "Writing Trees..." << std::endl; 

    outroot->cd();
    CutPromptTree->Write();
    CutDelayedTree->Write();
    outroot->Close();
    delete outroot;
}
