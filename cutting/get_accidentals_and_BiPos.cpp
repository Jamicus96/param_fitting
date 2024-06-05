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


void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, std::string outputVetoTxt, const bool is_data);


int main(int argv, char** argc) {
    std::string inputNtuple = argc[1];
    std::string previousRunNtuple = argc[2];
    std::string outputNtuple = argc[3];
    std::string outputVetoTxt = argc[4];
    bool is_data = std::stoi(argc[5]);

    std::cout << "Looping through events..." << std::endl;
    Apply_tagging_and_cuts(inputNtuple, previousRunNtuple, outputNtuple, outputVetoTxt, is_data);

    return 0;
}

void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, std::string outputVetoTxt, const bool is_data) {

    // Read in file and create output file
    TFile* inFile = TFile::Open(inputNtuple.c_str());
    TFile* outroot = new TFile(outputNtuple.c_str(), "RECREATE");

    TTree* EventInfo = (TTree*) inFile->Get("output");

    TTree* AccidentalPromptTree = EventInfo->CloneTree(0);
    AccidentalPromptTree->SetObject("promptAccidental", "promptAccidental");
    TTree* AccidentalDelayedTree = EventInfo->CloneTree(0);
    AccidentalDelayedTree->SetObject("delayedAccidental", "delayedAccidental");

    TTree* BiPoPromptTree = EventInfo->CloneTree(0);
    BiPoPromptTree->SetObject("promptBiPo", "promptBiPo");
    TTree* BiPoDelayedTree = EventInfo->CloneTree(0);
    BiPoDelayedTree->SetObject("delayedBiPo", "delayedBiPo");

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
    Double_t reconEnergy, reconX, reconY, reconZ;
    ULong64_t clockCount10, clockCount50, dcApplied, dcFlagged;
    Int_t mcIndex, nHits, GTID, owlnhits, UTDays, UTSecs, UTNSecs;
    Int_t runNum; //run number associated with MC
    Bool_t* valid;

    EventInfo->SetBranchAddress("energy", &reconEnergy);
    EventInfo->SetBranchAddress("posx", &reconX);
    EventInfo->SetBranchAddress("posy", &reconY);
    EventInfo->SetBranchAddress("posz", &reconZ);
    EventInfo->SetBranchAddress("clockCount10", &clockCount10);
    EventInfo->SetBranchAddress("clockCount50", &clockCount50);
    EventInfo->SetBranchAddress("uTDays", &UTDays);
    EventInfo->SetBranchAddress("uTSecs", &UTSecs);
    EventInfo->SetBranchAddress("uTNSecs", &UTNSecs);
    EventInfo->SetBranchAddress("fitValid", &valid);
    EventInfo->SetBranchAddress("mcIndex", &mcIndex);
    EventInfo->SetBranchAddress("dcApplied", &dcApplied);
    EventInfo->SetBranchAddress("dcFlagged", &dcFlagged);
    EventInfo->SetBranchAddress("runID", &runNum);
    EventInfo->SetBranchAddress("nhits", &nHits);
    EventInfo->SetBranchAddress("owlnhits", &owlnhits);
    EventInfo->SetBranchAddress("eventID", &GTID);

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    // Open file to print results muon veto times and info to
    std::ofstream outTxtFile;
    outTxtFile.open(outputVetoTxt, std::ofstream::out | std::ofstream::trunc);
    outTxtFile << "run_number GTID Nhit clockCount50 clockCount10 UTDays UTSecs UTNSecs Veto_length" << std::endl;

    // Initialise DetectorStateCorrection (assume only one run in each file)
    RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
    RAT::DU::ReconCalibrator* e_cal = RAT::DU::ReconCalibrator::Get();

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int nvaliddelayed = 0, nvalidpair = 0, nMuonCut = 0, nDCcut = 0, nValidCut = 0, negEcut = 0;
    int64_t delayedTime, promptTime;
    TVector3 delayedPos, promptPos, multiplicityPos;
    double delay, highNhitDelay, owlNhitDelay;
    double promptEcorr, delayedEcorr;
    // Loop through all events:
    double totHighNhitTime = 0, totOwlNhitTime = 0, vetoDeltaT;
    std::vector<int> cleaned_entry, promptBiPos, delayedBiPos;
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

        delayedTime = int64_t(clockCount50);
        if (nHits > 3000) {
            // run_number GTID Nhit clockCount50 clockCount10 UTDays UTSecs UTNSecs Veto_length
            // Muon veto time converted to ns
            outTxtFile << runNum << " " << GTID << " " << nHits << " " << clockCount50 << " " << clockCount10 << " " << UTDays << " " << UTSecs << " " << (double)(UTNSecs) << " " << highNhit_deltaT * 1E9 << std::endl;
            vetoDeltaT = ((delayedTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;
            if (vetoDeltaT > highNhit_deltaT) totHighNhitTime += highNhit_deltaT;
            else totHighNhitTime += vetoDeltaT;
            highNhitTime = delayedTime;
        }
        if (owlnhits > 3) {
            vetoDeltaT = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;
            if (vetoDeltaT > owlNhit_deltaT) totOwlNhitTime += owlNhit_deltaT;
            else totOwlNhitTime += vetoDeltaT;
            owlNhitTime = delayedTime;
        }
        highNhitDelay = ((delayedTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
        owlNhitDelay = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
        if (highNhitDelay < highNhit_deltaT || owlNhitDelay < owlNhit_deltaT) {++nMuonCut; continue;}
        if (!valid) {++nValidCut; continue;}
        if (reconEnergy < 0) {++negEcut; continue;}
        if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

        // If it survives the basic fitvalid + muon veto + data cleaning cuts, add to output ntuple
        cleaned_entry.push_back(a);

        delayedPos = TVector3(reconX, reconY, reconZ);
        delayedEcorr = EnergyCorrection(reconEnergy, delayedPos, is_data, stateCorr, e_cal);
        reconEnergy = promptEcorr;  // All usable events will have corrected energies now, so don't need to re-do it

        if (pass_delayed_cuts_BiPo(delayedEcorr, delayedPos)) {
            nvaliddelayed++;

            // Delayed event is valid, check through the previous 1000 events for event that passes prompt + classifier + tagging cuts
            for (int b = 1; b <= 1000; ++b) {
                if ((a - b) < 0) break;
                EventInfo->GetEntry(a - b);

                promptTime = int64_t(clockCount50);
                highNhitDelay = ((promptTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
                owlNhitDelay = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
                if (highNhitDelay < highNhit_deltaT) {++nMuonCut; break;}
                if (fabs(owlNhitDelay) < owlNhit_deltaT) {++nMuonCut; continue;}
                if (!valid) {++nValidCut; continue;}
                if (reconEnergy < 0) {++negEcut; continue;}
                if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

                delay = ((delayedTime - promptTime) & 0x7FFFFFFFFFF) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
                if (delay > BIPO_MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

                promptPos = TVector3(reconX, reconY, reconZ);

                if (pass_prompt_cuts_BiPo(reconEnergy, promptPos) and pass_coincidence_cuts_BiPo(delay, promptPos, delayedPos)) {
                    // Event pair survived analysis cuts
                    nvalidpair++;

                    promptBiPos.push_back(a - b);
                    delayedBiPos.push_back(a);

                    // Add to output TTrees
                    EventInfo->GetEntry(a - b);
                    // reconEnergy = promptEcorr;
                    BiPoPromptTree->Fill();
                    EventInfo->GetEntry(a);
                    // reconEnergy = delayedEcorr;
                    BiPoDelayedTree->Fill();
                }
            }
        }
    }
    outTxtFile.close();

    std::cout << "[BiPos] From " << nentries << " entries, number of valid delayed events: " << nvaliddelayed
              << ", number of these event pairs surviving prompt + coincidence cuts: " << nvalidpair << std::endl;
    std::cout << "Events cut from Muon tagging: " << nMuonCut << ", invalid recon: " << nValidCut << ", negative E: " << negEcut << ", failed DC: " << nDCcut << std::endl;
    std::cout << "High nhit cut time = " << totHighNhitTime << "s, owl nhit cut time = " << totOwlNhitTime * 1E-6 << "s" << std::endl;

    // Filter out BiPos from the rest
    bool notBiPo;
    unsigned int nvalidFV = 0, nvalidprompt = 0;
    nvaliddelayed = 0;
    for (unsigned int iEntry = 0; iEntry < cleaned_entry.size(); ++iEntry) {
        notBiPo = true;
        for (unsigned int iBiPo = 0; iBiPo < promptBiPos.size(); ++iBiPo) {
            if (cleaned_entry.at(iEntry) == promptBiPos.at(iBiPo) || cleaned_entry.at(iEntry) == delayedBiPos.at(iBiPo)) {
                notBiPo = false;
                break;
            }
        }

        if (notBiPo) {
            EventInfo->GetEntry(cleaned_entry.at(iEntry));
            
            if (sqrt(reconX*reconX + reconY*reconY + (reconZ - AV_offset)*(reconZ - AV_offset)) > FV_CUT) continue;
            ++nvalidFV;

            #ifdef USING_PDF_PADDING
                if (reconEnergy < IBD_PDF_MIN_PROMPT_E || reconEnergy > IBD_PDF_MAX_PROMPT_E) continue;
            #else
                if (reconEnergy < IBD_MIN_PROMPT_E || reconEnergy > IBD_MAX_PROMPT_E) continue;
            #endif
            nvalidprompt++;
            AccidentalPromptTree->Fill();

            if (reconEnergy < IBD_MIN_DELAYED_E || reconEnergy > IBD_MAX_DELAYED_E) continue;
            nvaliddelayed++;
            AccidentalDelayedTree->Fill();
        }
    }

    std::cout << "[Accidentals] From " << cleaned_entry.size() << " cleaned entries, number of valid FV + not BiPo events: " << nvalidFV
              << ", number of valid prompt events : " << nvalidprompt << ", number of valid delayed events : " << nvaliddelayed << std::endl;

    // Write output ntuple to files, deal with ttrees getting too full to write to one file
    std::cout << "Writing Trees..." << std::endl; 

    outroot->cd();
    BiPoPromptTree->Write();
    BiPoDelayedTree->Write();
    AccidentalPromptTree->Write();
    AccidentalDelayedTree->Write();
    outroot->Close();
    inFile->Close();
    delete outroot;
}
