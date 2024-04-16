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

#define VERBOSE
#define USING_RUN_NUM
ULong64_t dcAnalysisWord = 0x2100000042C2;  // Converts hex to decimal

// Define cut values
double MIN_PROMPT_E = 0.9, MAX_PROMPT_E = 8.0;  // [MeV]
double MIN_DELAYED_E = 1.85, MAX_DELAYED_E = 2.4;  // [MeV]
double FV_CUT = 5700;  // Max radius [mm]
double MIN_DELAY = 400, MAX_DELAY = 0.8E6;  // Delta t cut [ns]
double MAX_DIST = 1500;  // Delta R cut [mm]
double PROTON_RECOIL_E_MAX = 3.5;  // [MeV]  Ideally the same as in make_PDF.cpp
double highNhit_deltaT = 20;  // [s]
double owlNhit_deltaT = 10;  // [us]
double AV_offset = 186;  // [mm] average value from https://snopl.us/monitoring/scint_level


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


bool pass_prompt_cuts(const double energy, const TVector3& position) {
    if (energy < MIN_PROMPT_E) return false;  // min energy cut (MeV)
    if (energy > MAX_PROMPT_E) return false;  // max energy cut (MeV)
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(const double energy, const TVector3& position) {
    if (energy < MIN_DELAYED_E) return false;  // min energy cut (MeV)
    if (energy > MAX_DELAYED_E) return false;  // max energy cut (MeV)
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

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


double EnergyCorrection(const double E, TVector3 pos, const bool is_data, RAT::DU::DetectorStateCorrection& stateCorr, RAT::DU::ReconCalibrator* e_cal) {
    // Data vs MC energy correction (Tony's)
    double Ecorr = e_cal->CalibrateEnergyRTF(is_data, E, std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()), pos.Z()); // gives the new E

    // #ifdef USING_RUN_NUM
    //     // Correct for position coverage dependence (Logan's)
    //     RAT::DU::Point3D position(0, pos);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
    //     Ecorr /= stateCorr.GetCorrectionPos(position, 0, 0) / stateCorr.GetCorrection(9394, 0.75058); // a correction factor (divide E by it)
    // #endif

    return Ecorr;
}

bool dcAppliedAndPassed(const bool is_data, ULong64_t dcApplied, ULong64_t dcFlagged) {
    if (!is_data) return true;
    if (!(((dcApplied & dcAnalysisWord) & dcFlagged ) == (dcApplied & dcAnalysisWord))) return false;  // failed DC cut
    return true;
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
    Int_t mcIndex, nHits, GTID, owlnhits;
    Int_t runNum; //run number associated with MC
    Bool_t* valid;

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
    EventInfo->SetBranchAddress("runID", &runNum);
    EventInfo->SetBranchAddress("nhits", &nHits);
    EventInfo->SetBranchAddress("owlnhits", &owlnhits);
    EventInfo->SetBranchAddress("eventID", &GTID);

    /* ~~~~~~~ loop through events, tag and cut ~~~~~~~ */

    // Initialise DetectorStateCorrection (assume only one run in each file)
    RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
    RAT::DU::ReconCalibrator* e_cal = RAT::DU::ReconCalibrator::Get();

    unsigned int nentries = EventInfo->GetEntries();
    unsigned int nvaliddelayed = 0, nvalidpair = 0, nvalid = 0, nMuonCut = 0, nDCcut = 0, nValidCut = 0, negEcut = 0, nvalidMultiplicity = 0;
    int64_t delayedTime, promptTime;
    TVector3 delayedPos, promptPos, multiplicityPos;
    double delay, highNhitDelay, owlNhitDelay;
    bool passMultiplicity;
    double promptEcorr, delayedEcorr;
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

        delayedTime = int64_t(eventTime);
        if (nHits > 3000) highNhitTime = delayedTime;
        if (owlnhits > 3) owlNhitTime = delayedTime;
        highNhitDelay = ((delayedTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
        owlNhitDelay = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
        if (highNhitDelay < highNhit_deltaT || owlNhitDelay < owlNhit_deltaT) {++nMuonCut; continue;}
        if (!valid) {++nValidCut; continue;}
        if (reconEnergy < 0) {++negEcut; continue;}
        if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

        delayedPos = TVector3(reconX, reconY, reconZ);
        delayedEcorr = EnergyCorrection(reconEnergy, delayedPos, is_data, stateCorr, e_cal);

        if (pass_delayed_cuts(delayedEcorr, delayedPos)) {
            nvaliddelayed++;

            // Delayed event is valid, check through the previous 1000 events for event that passes prompt + classifier + tagging cuts
            for (int b = 1; b <= 1000; ++b) {
                if ((a - b) < 0) break;
                EventInfo->GetEntry(a - b);

                promptTime = int64_t(eventTime);
                highNhitDelay = ((promptTime - highNhitTime) & 0x7FFFFFFFFFF) / 50E6;  // [s] dealing with clock rollover
                owlNhitDelay = ((delayedTime - owlNhitTime) & 0x7FFFFFFFFFF) / 50.0;  // [us] dealing with clock rollover
                if (highNhitDelay < highNhit_deltaT) {++nMuonCut; break;}
                if (fabs(owlNhitDelay) < owlNhit_deltaT) {++nMuonCut; continue;}
                if (!valid) {++nValidCut; continue;}
                if (reconEnergy < 0) {++negEcut; continue;}
                if (!dcAppliedAndPassed(is_data, dcApplied, dcFlagged)) {++nDCcut; continue;}

                delay = ((delayedTime - promptTime) & 0x7FFFFFFFFFF) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
                if (delay > MAX_DELAY) break;  // If delay becomes larger than cut, stop looking for new events

                promptPos = TVector3(reconX, reconY, reconZ);
                promptEcorr = EnergyCorrection(reconEnergy, promptPos, is_data, stateCorr, e_cal);

                if (pass_prompt_cuts(promptEcorr, promptPos) and pass_coincidence_cuts(delay, promptPos, delayedPos)) {
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

    // Write output ntuple to files, deal with ttrees getting too full to write to one file
    std::cout << "Writing Trees..." << std::endl; 

    outroot->cd();
    CutPromptTree->Write();
    CutDelayedTree->Write();
    outroot->Close();
    delete outroot;
}
