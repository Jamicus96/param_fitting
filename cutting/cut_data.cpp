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

void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, const double classifier_cut, const bool is_data, const std::string eventType);


int main(int argv, char** argc) {
    std::string inputNtuple = argc[1];
    std::string previousRunNtuple = argc[2];
    std::string outputNtuple = argc[3];
    double classifier_cut = std::stod(argc[4]);
    bool is_data = std::stoi(argc[5]);
    std::string eventType = argc[6];
    
    // Apply cuts and tagging!
    std::cout << "Looping through events..." << std::endl;
    Apply_tagging_and_cuts(inputNtuple, previousRunNtuple, outputNtuple, classifier_cut, is_data, eventType);

    return 0;
}

void Apply_tagging_and_cuts(std::string inputNtuple, std::string previousRunNtuple, std::string outputNtuple, const double classifier_cut, const bool is_data, const std::string eventType) {

    // Read in file and create output file
    TFile* inFile = TFile::Open(inputNtuple.c_str());
    TFile* outroot = new TFile(outputNtuple.c_str(), "RECREATE");

    TTree* EventInfo = (TTree*) inFile->Get("output");
    TTree* CutPromptTree = EventInfo->CloneTree(0);
    CutPromptTree->SetObject("prompt", "prompt");
    TTree* CutDelayedTree = EventInfo->CloneTree(0);
    CutDelayedTree->SetObject("delayed", "delayed");

    // Make histograms while at it [R_AV, delta_R, delta_T]
    unsigned int nbins = 200;
    std::vector<TH1D*> hists;
    hists.push_back(new TH1D("R_AV", "Recon R_AV [mm] (for true R_AV < 5.7 m)", nbins, 0, 6000));
    hists.push_back(new TH1D("delta_R", "Recon delta R [mm]", nbins, 0, 6000));
    hists.push_back(new TH1D("delta_T", "Recon delta T [ns]", nbins, 0, 1.5E6));

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
    std::vector<int64_t> owlNhitTimes = {-99999999999};
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
                owlNhitTimes.push_back(int64_t(prev_eventTime));
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
    Int_t pdg1, pdg2, mcIndex;

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

    EventInfo->SetBranchAddress("mcPosx", &mcX);
    EventInfo->SetBranchAddress("mcPosy", &mcY);
    EventInfo->SetBranchAddress("mcPosz", &mcZ);
    EventInfo->SetBranchAddress("evIndex", &EVindex);
    EventInfo->SetBranchAddress("pdg1", &pdg1);
    EventInfo->SetBranchAddress("pdg2", &pdg2);
    EventInfo->SetBranchAddress("mcIndex", &mcIndex);

    /* ~~~~~~~ Clean up events first ~~~~~~~ */

    unsigned int nentries = EventInfo->GetEntries();

    unsigned int numNotVetoed = 0, nEvtType = 0, numInsideFV_MC = 0, nEvtType_InsideFV_MC = 0, nValid = 0,
                 nEpos = 0, nEvtType_InsideFV_MC_Valid = 0, nFVpass = 0, nEvtType_InsideFV = 0;
    int64_t delayedTime;
    TVector3 delayedPos;
    double highNhitDelay, owlNhitDelay;

    // Loop through all events:
    std::vector<int> cleaned_entries;
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
        highNhitDelay = ((delayedTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
        if (highNhitDelay < highNhit_deltaT) continue;

        if (owlnhits > 3) owlNhitTimes.push_back(delayedTime);
        owlNhitDelay = ((delayedTime - owlNhitTimes.at(owlNhitTimes.size()-1)) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
        if (owlNhitDelay < owlNhit_deltaT) continue;

        ++numNotVetoed;
        if (check_event_type(pdg1, pdg2, eventType) && EVindex == 0) nEvtType++;

        // Counting number of events simulated within FV, outside of veto windows
        delayedPos = TVector3(mcX, mcY, mcZ);
        if (pass_FV_cut(delayedPos)) {
            numInsideFV_MC++;
            // Is it of the correct event type simulated?
            if (check_event_type(pdg1, pdg2, eventType) && EVindex == 0) {
                hists.at(0)->Fill(sqrt(reconX*reconX + reconY*reconY + (reconZ - AV_offset)*(reconZ - AV_offset))); // [R_AV, delta_R, delta_T]
                nEvtType_InsideFV_MC++;
            }
        }

        // Fit-valid and DC logic (accounted for in detector/cut efficiency)
        if (!valid) continue;
        nValid++;
        if (reconEnergy < 0) continue;
        nEpos++;

        if (pass_FV_cut(delayedPos) && check_event_type(pdg1, pdg2, eventType) && EVindex == 0) nEvtType_InsideFV_MC_Valid++;

        // FV cut
        delayedPos = TVector3(reconX, reconY, reconZ);
        if (!pass_FV_cut(delayedPos)) continue;
        nFVpass++;

        // Is it of the correct event type simulated?
        if (check_event_type(pdg1, pdg2, eventType) && EVindex == 0) nEvtType_InsideFV++;

        cleaned_entries.push_back(a);
    }

    /* ~~~~~~~ loop through events, cut loosely and tag ~~~~~~~ */

    // Initialise DetectorStateCorrection (assume only one run in each file)
    RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
    RAT::DU::ReconCalibrator* e_cal = RAT::DU::ReconCalibrator::Get();

    unsigned int nDelayedPass = 0, nPromptPass = 0, nDtPass = 0, nDrPass = 0;
    int64_t promptTime;
    TVector3 promptPos, multiplicityPos;
    double delay;
    bool passMultiplicity, passOWLnhit;
    double promptEcorr, delayedEcorr;
    // Loop through all events:
    int a, b;
    std::vector<int> prompt_entries_loose, delayed_entries_loose;
    for (int i_entry = 0; i_entry < cleaned_entries.size(); ++i_entry) {
        if (i_entry % 10000 == 0) std::cout << "Done " << ((float)i_entry / (float)cleaned_entries.size()) * 100.0 << "%" << std::endl;

        a = cleaned_entries.at(i_entry);
        EventInfo->GetEntry(a);

        #ifdef USING_RUN_NUM
            if (lastRunNum != runNum) {
                run.SetRunID(runNum);
                db->BeginOfRun(run);
                lastRunNum = runNum;
            }
        #endif

        // Energy correction
        delayedTime = int64_t(eventTime);
        delayedPos = TVector3(reconX, reconY, reconZ);
        delayedEcorr = EnergyCorrection(reconEnergy, delayedPos, is_data, stateCorr, e_cal);

        // delayed cut
        if (!(delayedEcorr > 1.3 && delayedEcorr < 3.2)) continue;
        ++nDelayedPass;

        // Delayed event is valid, check through the previous 1000 events for event that passes prompt + classifier + tagging cuts
        for (int j_entry = 1; j_entry < cleaned_entries.size(); ++j_entry) {
            if ((i_entry - j_entry) < 0) break;
            b = cleaned_entries.at(j_entry);
            EventInfo->GetEntry(a - b);

            promptTime = int64_t(eventTime);
            delay = ((delayedTime - promptTime) & 0x7FFFFFFFFFF) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
            if (delay > IBD_MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

            promptPos = TVector3(reconX, reconY, reconZ);
            promptEcorr = EnergyCorrection(reconEnergy, promptPos, is_data, stateCorr, e_cal);

            if (!pass_prompt_PDF_cuts_IBD(promptEcorr)) continue;
            nPromptPass++;
            hists.at(2)->Fill(delay); // [R_AV, delta_R, delta_T]

            // Coincidenc cuts
            if (!pass_dT_cut_IBD(delay)) continue;
            nDtPass++;
            hists.at(1)->Fill((delayedPos - promptPos).Mag()); // [R_AV, delta_R, delta_T]

            if (!pass_dR_cut_IBD(promptPos, delayedPos)) continue;
            nDrPass++;

            // Add to output TTrees
            EventInfo->GetEntry(a - b);
            reconEnergy = promptEcorr;
            prompt_entries_loose.push_back(a - b);

            EventInfo->GetEntry(a);
            reconEnergy = delayedEcorr;
            delayed_entries_loose.push_back(a);
            CutDelayedTree->Fill();  // Just record for nice plots
        }
    }

    /* ~~~~~~~ tighten delayed cuts ~~~~~~~ */

    std::vector<int> prompt_entries_mid, delayed_entries_mid;
    for (int i_entry = 0; i_entry < delayed_entries_loose.size(); ++i_entry) {
        if (i_entry % 10000 == 0) std::cout << "Done " << ((float)i_entry / (float)delayed_entries_loose.size()) * 100.0 << "%" << std::endl;

        a = delayed_entries_loose.at(i_entry);
        EventInfo->GetEntry(a);

        #ifdef USING_RUN_NUM
            if (lastRunNum != runNum) {
                run.SetRunID(runNum);
                db->BeginOfRun(run);
                lastRunNum = runNum;
            }
        #endif

        // delayed cut
        if (pass_delayed_cuts_IBD(reconEnergy)) {
            prompt_entries_mid.push_back(prompt_entries_loose.at(i_entry));
            delayed_entries_mid.push_back(a);

            // Record for PDF making
            EventInfo->GetEntry(prompt_entries_loose.at(i_entry));
            CutPromptTree->Fill();
        }
    }

    /* ~~~~~~~ tighten prompt cut, just for efficiency ~~~~~~~ */

    unsigned int num_tagged_final = 0;
    for (int i_entry = 0; i_entry < prompt_entries_mid.size(); ++i_entry) {
        if (i_entry % 10000 == 0) std::cout << "Done " << ((float)i_entry / (float)prompt_entries_mid.size()) * 100.0 << "%" << std::endl;

        a = prompt_entries_mid.at(i_entry);
        EventInfo->GetEntry(a);

        #ifdef USING_RUN_NUM
            if (lastRunNum != runNum) {
                run.SetRunID(runNum);
                db->BeginOfRun(run);
                lastRunNum = runNum;
            }
        #endif

        // delayed cut
        if (pass_prompt_cuts_IBD(reconEnergy)) ++num_tagged_final;
    }


    /* ~~~~~~~ Print info for efficiency, and record data to ntuple file ~~~~~~~ */

    std::cout << "Cut efficiencies from cleaning [nentries, numNotVetoed, nEvtType, numInsideFV_MC, nEvtType_InsideFV_MC, nValid, nEpos, nEvtType_InsideFV_MC_Valid, nFVpass, nEvtType_InsideFV]:" << std::endl;
    std::cout << nentries << ", " << numNotVetoed << ", " << nEvtType << ", " << numInsideFV_MC << ", " << nEvtType_InsideFV_MC << ", " << nValid << ", "
              << nEpos << ", " << nEvtType_InsideFV_MC_Valid << ", " << nFVpass << ", " << nEvtType_InsideFV << std::endl;

    std::cout << "Cut efficiencies from loose cuts and tagging [nDelayedPass, nPromptPass, nDtPass, nDrPass]:" << std::endl;
    std::cout << nDelayedPass << ", " << nPromptPass << ", " << nDtPass << ", " << nDrPass << std::endl;

    std::cout << "Cut efficiencies from tighter cuts [nDelayed_cut_pairs, num_tagged_final]:" << std::endl;
    std::cout << prompt_entries_mid.size() << ", " << num_tagged_final << std::endl;

    std::cout << "Last muon tag time (50MHz clock ticks): " << highNhitTime << std::endl;

    // Write output ntuple to files, deal with ttrees getting too full to write to one file
    std::cout << "Writing Trees..." << std::endl; 

    outroot->cd();
    CutPromptTree->Write();
    CutDelayedTree->Write();
    for (unsigned int i = 0; i < hists.size(); ++i) {hists.at(i)->Write();}
    outroot->Close();
    delete outroot;
}
