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
#include "cut_utils.hpp"


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
    std::vector<TH1D*> recordHists = Apply_tagging_and_cuts(EventInfo, CutTtrees, classifier_cut, is_data);

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