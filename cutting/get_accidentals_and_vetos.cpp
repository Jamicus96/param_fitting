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


void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, std::string outputVetoTxt, const double classifier_cut, const bool is_data);


int main(int argv, char** argc) {
    std::string inputNtuple = argc[1];
    std::string previousRunNtuple = argc[2];
    std::string outputNtuple = argc[3];
    std::string outputVetoTxt = argc[4];
    double classifier_cut = std::stod(argc[5]);
    bool is_data = std::stoi(argc[6]);

    std::cout << "Looping through events..." << std::endl;
    Apply_tagging_and_cuts(inputNtuple, previousRunNtuple, outputNtuple, outputVetoTxt, classifier_cut, is_data);

    return 0;
}

void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, std::string outputVetoTxt, const double classifier_cut, const bool is_data) {

    // Read in file and create output file
    TFile* inFile = TFile::Open(inputNtuple.c_str());
    TFile* outroot = new TFile(outputNtuple.c_str(), "RECREATE");

    TTree* EventInfo = (TTree*) inFile->Get("output");

    TTree* AccidentalPromptTree = EventInfo->CloneTree(0);
    AccidentalPromptTree->SetObject("promptAccidental", "promptAccidental");
    TTree* AccidentalDelayedTree = EventInfo->CloneTree(0);
    AccidentalDelayedTree->SetObject("delayedAccidental", "delayedAccidental");

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
    EventInfo->SetBranchAddress("alphaNReactorIBD", &classResult);

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    // Open file to print results muon veto times and info to
    std::ofstream outTxtFile;
    outTxtFile.open(outputVetoTxt, std::ofstream::out | std::ofstream::trunc);
    outTxtFile << "run_number GTID Nhit clockCount50 clockCount10 UTDays UTSecs UTNSecs Veto_length" << std::endl;

    // Initialise DetectorStateCorrection (assume only one run in each file)
    RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
    RAT::DU::ReconCalibrator* e_cal = RAT::DU::ReconCalibrator::Get();

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int nValid = 0, nEpos = 0, nDCpass = 0, nFVpass = 0;
    unsigned int nPromptPass = 0, nDelayedPass = 0, nClassPass = 0;
    unsigned int nMuonCut = 0;
    int64_t eventTime;
    TVector3 eventPos;
    double highNhitDelay, owlNhitDelay;
    double eventEcorr;
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

        // Veto logic (accounted for in livetime)
        eventTime = int64_t(clockCount50);
        if (nHits > 3000) {
            // run_number GTID Nhit clockCount50 clockCount10 UTDays UTSecs UTNSecs Veto_length
            // Muon veto time converted to ns
            outTxtFile << runNum << " " << GTID << " " << nHits << " " << clockCount50 << " " << clockCount10 << " " << UTDays << " " << UTSecs << " " << (double)(UTNSecs) << " " << highNhit_deltaT * 1E9 << std::endl;
            vetoDeltaT = ((eventTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;
            if (vetoDeltaT > highNhit_deltaT) totHighNhitTime += highNhit_deltaT;
            else totHighNhitTime += vetoDeltaT;
            highNhitTime = eventTime;
        }
        highNhitDelay = ((eventTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
        if (highNhitDelay < highNhit_deltaT) {++nMuonCut; continue;}
        
        if (owlnhits > 3) {
            vetoDeltaT = ((eventTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;
            if (vetoDeltaT > owlNhit_deltaT) totOwlNhitTime += owlNhit_deltaT;
            else totOwlNhitTime += vetoDeltaT;
            owlNhitTime = eventTime;
        }
        owlNhitDelay = ((eventTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
        if (owlNhitDelay < owlNhit_deltaT) {++nMuonCut; continue;}

        // Fit-valid and DC logic (accounted for in detector/cut efficiency)
        if (!valid) continue;
        nValid++;
        if (reconEnergy < 0) continue;
        nEpos++;
        if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) continue;
        nDCpass++;

        // FV cut
        eventPos = TVector3(reconX, reconY, reconZ);
        if (!pass_FV_cut(eventPos)) continue;
        nFVpass++;

        // Energy correction
        eventEcorr = EnergyCorrection(reconEnergy, eventPos, is_data, stateCorr, e_cal);

        // Prompt E cut (just for cut eff)
        if (!pass_prompt_cuts_IBD(eventEcorr)) continue;
        nPromptPass++;
        reconEnergy = eventEcorr;

        // Apply classifier to prompt events
        if (pass_classifier(eventEcorr, classResult, classifier_cut)) {
            nClassPass++;
            AccidentalPromptTree->Fill();
        }

        // delayed cut
        if (!pass_delayed_cuts_IBD(eventEcorr)) continue;
        nDelayedPass++;
        AccidentalDelayedTree->Fill();
    }
    outTxtFile.close();

    std::cout << "Cut efficiencies [nentries, nValid, nEpos, nDCpass, nFVpass, nPromptPass, nDelayedPass, nClassPass, nMuonCut]:" << std::endl;
    std::cout << nentries << ", " << nValid << ", " << nEpos << ", " << nDCpass << ", "
              << nFVpass << ", " << nPromptPass << ", " << nDelayedPass << ", "
              << nClassPass << ", " << nMuonCut << std::endl;
    std::cout << "highNhitDelay = " << highNhitDelay << ", owlNhitDelay = " << owlNhitDelay << std::endl;

    std::cout << "Last muon tag time (50MHz clock ticks): " << highNhitTime << std::endl;

    // Write output ntuple to files, deal with ttrees getting too full to write to one file
    std::cout << "Writing Trees..." << std::endl; 

    outroot->cd();
    AccidentalPromptTree->Write();
    AccidentalDelayedTree->Write();
    outroot->Close();
    inFile->Close();
    delete outroot;
}
