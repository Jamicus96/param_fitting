/***
Based off rat-tools/AntinuTools/scale_reactor_flux.cpp
***/

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TRandom3.h>
#include <RAT/DB.hh>
#include <TMath.h>
#include <TVector3.h>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#define USING_RUN_NUM
// #define VERBOSE

/**
 * @brief Get the livetime scaling due to muon vetos from each run
 * 
 * @param run_livetimes_address Output "livetime_per_run_*.txt" from livetime calculator
 * @param livetime_scalings map <runNum, scaling>
 */
void GetVetoLivetimeScaling(std::string run_livetimes_address, std::map<int, double>& livetime_scalings) {
    if (run_livetimes_address == "" || run_livetimes_address == "false" || run_livetimes_address == "False" || run_livetimes_address == "0") return;

    // Read in text file "<runNum> <livetime [days]> <livetime - veto time [days]"
    std::cout << "Reading in livetimes: " << run_livetimes_address << std::endl;
    std::ifstream infile(run_livetimes_address);
    if (!infile) {std::cout << "Error opening the file. Exit." << std::endl; exit(1);}

    int runNum;
    double livetime, corr_livetime;
    while (infile >> runNum >> livetime >> corr_livetime) {
        if (corr_livetime > livetime) {
            std::cout << "ERROR [corr_livetime > livetime] for run " << runNum << ": " << corr_livetime << " > " << livetime << ". Setting scaling to 1." << std::endl;
            livetime_scalings.insert({runNum, 1.0});
        }
        if (corr_livetime == 0) {
            std::cout << "ERROR [corr_livetime = 0] for run " << runNum << ": corr_livetime = " << corr_livetime << ", livetime = " << livetime << ". Setting scaling to 1." << std::endl;
            livetime_scalings.insert({runNum, 1.0});
        }
        else livetime_scalings.insert({runNum, corr_livetime / livetime});
    }

    #ifdef VERBOSE
        std::cout << "livetime_scalings:" << std::endl;
        for (auto& x : livetime_scalings) {  
            std::cout << "[" <<  x.first << "] : " << x.second << std::endl;;
        }
    #endif
}


std::vector<std::string> SplitString(std::string str, bool Osc) {
    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end); //each word of string now in vector
    std::vector<std::string> info;
    std::string dummy = "";

    for(int i=0;i<tokens.size();i++){ //combine back to reactor name, core number
        if(i==0){
            dummy = tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else if(i!=tokens.size()-1){
            dummy += tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else{
            info.push_back(dummy);
            info.push_back(tokens.at(i));
        }
    }

    return info;
}

void ScaleReactorFlux(std::string inputFile, std::string outputFile, const std::map<int, double>& livetime_scalings) {

    //load in ntuple
    TFile *reactorFile;
    try{
        reactorFile = TFile::Open(inputFile.c_str());
    }
    catch(...){
        std::cout << "Could not open input file." << std::endl;
        exit(1);
    }

    TFile *scaled_ntuple = new TFile(outputFile.c_str(),"RECREATE");

    //clone the ntuple, we want to keep the original ntuple untouched and produce a new file after scaling
    TTree *reactorEventInfo = (TTree *) reactorFile->Get("prompt"); 
    TTree *scaledReactorEventInfo = reactorEventInfo->CloneTree(0);

    //Set branch addresses for relevant parameters
    Double_t parentKE1; // antineutrino energy
    Double_t reconEnergy; //recon energy
    Int_t runNum; //run number associated with MC
    Int_t lastRunNum = -999; //last run number loaded in DB
    TString *originReactor = NULL; // string for reactor where antineutrino was produced
    Bool_t *valid;

    reactorEventInfo->SetBranchAddress("parentMeta1",&originReactor);
    reactorEventInfo->SetBranchAddress("parentKE1",&parentKE1);
    reactorEventInfo->SetBranchAddress("energy",&reconEnergy);
    reactorEventInfo->SetBranchAddress("runID",&runNum);
    reactorEventInfo->SetBranchAddress("fitValid",&valid);

    //loop over events and decide whether to keep the event or not, initialise variables and histograms first

    Int_t nentries = reactorEventInfo->GetEntries();
    TRandom3 *rndm_dummy = new TRandom3();
    Double_t rndm;
    Double_t efficiency;
    std::vector<std::string> reacts_eff_names;
    std::vector<Double_t> reacts_eff_vals;
    int num_excess_eff_evts = 0;
    int num_no_table_evts = 0;
    double av_pow_scaling = 0, av_livetime_scaling = 0;
    double livetime_scaling = 1;
    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    RAT::DS::Run run;

    db->LoadDefaults();
    db->SetServer("postgres://snoplus@pgsql.snopl.us:5400/ratdb");
    run.SetRunID(runNum);
    db->BeginOfRun(run);
    std::cout << "RAT DB tag: " << db->GetDBTag() << std::endl;

    #ifdef VERBOSE
        std::cout << "runNum = " << runNum << std::endl;
    #endif

    if (livetime_scalings.find(runNum) != livetime_scalings.end()) {
        livetime_scaling = livetime_scalings.at(runNum);
        #ifdef VERBOSE
            std::cout << "livetime_scalings.at(runNum) = " << livetime_scalings.at(runNum) << std::endl;
        #endif
    }

    //go through entries of ttree
    for(Int_t a=0; a<nentries; a++) {
        reactorEventInfo->GetEntry(a);
        const char *originReactorString = originReactor->Data();
        std::vector<std::string> originReactorVect = SplitString(originReactorString, 0);

        rndm = rndm_dummy->Rndm();

        if (lastRunNum != runNum) {
            run.SetRunID(runNum);
            db->BeginOfRun(run);
            lastRunNum = runNum;
            #ifdef VERBOSE
                std::cout << "runNum = " << runNum << std::endl;
            #endif
            if (livetime_scalings.find(runNum) != livetime_scalings.end()) {
                livetime_scaling = livetime_scalings.at(runNum);
                #ifdef VERBOSE
                    std::cout << "livetime_scalings.at(runNum) = " << livetime_scalings.at(runNum) << std::endl;
                #endif
            } else {
                livetime_scaling = 1;
            }
        }

        // Get power scaling
        linkdb = db->GetLink("REACTOR_STATUS", originReactorVect[0]);
        try{
            efficiency = linkdb->GetDArray("core_power_scale_factor")[std::stoi(originReactorVect[1])];
        }
        catch(...){
            num_no_table_evts++;
            std::cout << "Reactors status table does not exist for reactor " << originReactorString << ". Removing this event." << std::endl;
            continue;
        }

        // Rolling averages (for printing)
        av_pow_scaling += (efficiency - av_pow_scaling) / (double)(a + 1);
        av_livetime_scaling += (livetime_scaling - av_livetime_scaling) / (double)(a + 1);

        // Add veto time scaling
        efficiency *= livetime_scaling;

        // See if event survives re-scaling
        if(efficiency > 1){ // we cant add events, so for now make a note of which reactors this occurs, SOLUTION: simulate higher than design power
            num_excess_eff_evts++;
            if(std::find(reacts_eff_names.begin(), reacts_eff_names.end(),originReactorString)==reacts_eff_names.end()){
                reacts_eff_names.push_back(originReactorString);
            }
        }
        if(rndm < efficiency){ // if keeping the event, write to new ttree
            scaledReactorEventInfo->Fill();
        }
    }

    //write and close everything, give user info

    Int_t nentries_scaled = scaledReactorEventInfo->GetEntries();
    std::cout << "Raw entries: " << nentries << ", scaled entries: " << nentries_scaled << ", removed entries: " << nentries-nentries_scaled << std::endl;
    if(num_excess_eff_evts != 0){
        std::cout << "The following reactors had efficiencies > 1, events were not removed:" << std::endl;
        for(int i=0;i<reacts_eff_names.size();i++){
            std::cout << reacts_eff_names.at(i) << std::endl;
        }
        std::cout << "This was the case for " << num_excess_eff_evts << " events, equal to " << (float(num_excess_eff_evts) / nentries) * 100 << " percent of events" << std::endl;
    }
    if(num_no_table_evts != 0){
        std::cout << "Could not find a table for " << num_no_table_evts << " events, equal to " << (float(num_no_table_evts) / nentries) * 100 << " percent of events" << std::endl;
    }

    std::cout << "Average power scaling applied: " << av_pow_scaling << std::endl;
    std::cout << "Average livetime scaling applied: " << av_livetime_scaling << std::endl;

    scaled_ntuple->cd();
    scaledReactorEventInfo->Write();

    delete reactorEventInfo;
    delete scaledReactorEventInfo;
    reactorFile->Close();
    scaled_ntuple->Close();
}

//main

int main(int argv, char** argc){
    std::string inputNtuple = argc[1];
    std::string run_livetimes_address = argc[2];
    std::string outputNtuple = argc[3];

    // Get run livetime scalings from muon vetos
    std::map<int, double> livetime_scalings;
    GetVetoLivetimeScaling(run_livetimes_address, livetime_scalings);

    ScaleReactorFlux(inputNtuple, outputNtuple, livetime_scalings);
} 
