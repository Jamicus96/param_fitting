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

#define USING_RUN_NUM


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

void ScaleReactorFlux(std::string inputFile, std::string outputFile) {

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
    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    RAT::DS::Run run;

    #ifndef USING_RUN_NUM
        db->SetAirplaneModeStatus(true);
        db->LoadDefaults();
    #endif
    #ifdef USING_RUN_NUM
        db->LoadDefaults();
        db->SetServer("postgres://snoplus@pgsql.snopl.us:5400/ratdb");
        run.SetRunID(runNum);
        db->BeginOfRun(run);
        std::cout << "RAT DB tag: " << db->GetDBTag() << std::endl;
    #endif

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
        }

        linkdb = db->GetLink("REACTOR_STATUS", originReactorVect[0]);
        try{
            efficiency = linkdb->GetDArray("core_power_scale_factor")[std::stoi(originReactorVect[1])];
        }
        catch(...){
            num_no_table_evts++;
            std::cout << "Reactors status table does not exist for reactor " << originReactorString << ". Removing this event." << std::endl;
            continue;
        }
        if(efficiency > 1){ //we cant add events, so for now make a note of which reactors this occurs, SOLUTION: simulate higher than design power
            num_excess_eff_evts++;
            if(std::find(reacts_eff_names.begin(), reacts_eff_names.end(),originReactorString)==reacts_eff_names.end()){
                reacts_eff_names.push_back(originReactorString);
            }
        }
        if(rndm<efficiency){ //if keeping the event, write to new ttree
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
    std::string outputNtuple = argc[2];

    ScaleReactorFlux(inputNtuple, outputNtuple);
} 
