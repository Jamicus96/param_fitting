#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "fitVars.hpp"
#include "E_systematics.hpp"


double gaussian(const double mean, const double sigma, const double amplitude, const double x);
TH1D* GetHist(std::string root_filename, std::string hist_name);

int main(int argv, char** argc) {
    // Read in arguments (Can provide histogram from root file. Otherwise a gaussian is generated)
    double linScale = std::stod(argc[1]);
    double kB = std::stod(argc[2]);
    double kBp = std::stod(argc[3]);
    double sigPerSqrtE = std::stod(argc[4]);

    TH1D* initHist;
    if (argv >= 6) {
        initHist = GetHist(argc[5], argc[6]);
    } else {
        // Create histogram and fill with gaussian values
        initHist = new TH1D("init", "Initial Histogram", 1000, 0.0, 100.0);
        double mean = 50.0;
        double sigma = 5.0;
        double amplitude = 10.0;
        double x;
        for (unsigned int i = 1; i < initHist->GetXaxis()->GetNbins() + 1; ++i) {
            x = initHist->GetBinCenter(i);
            initHist->SetBinContent(i, gaussian(mean, sigma, amplitude, x));
        }
    }

    // Create energy systematics object, with associate variables
    FitVars* Vars = FitVars::GetInstance();
    Esys* Esysts = Esys::GetInstance();
    Vars->AddVar("linScale", linScale, 0, linScale, linScale, true);
    Vars->AddVar("kBp", kBp, 0, kBp, kBp, true);
    Vars->AddVar("sigPerSqrtE", sigPerSqrtE, 0, sigPerSqrtE, sigPerSqrtE, true);
    Esysts->AddEsys("EsysName", kB, "linScale", "kBp", "sigPerSqrtE");

    // Apply energy systematics
    TH1D* finalHist = (TH1D*)(initHist->Clone("final"));
    finalHist->Reset("ICES");
    finalHist->SetTitle("Final Histogram");
    Esysts->apply_systematics("EsysName", initHist, finalHist);

    // Save to output root file and close
    TFile *outroot = new TFile("testEsys_outRoot.root", "RECREATE");
    initHist->Write();
    finalHist->Write();
    outroot->Write();
    outroot->Close();

    return 0;
}

double gaussian(const double mean, const double sigma, const double amplitude, const double x) {
    return amplitude * exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma)) / (sqrt(2.0 * 3.14) * sigma);
}

TH1D* GetHist(std::string root_filename, std::string hist_name) {
    // Read in root file
    TFile* fin = new TFile(root_filename.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file " << root_filename << std::endl;
        exit(1);
    }

    // Get histogram
    TH1D* h = (TH1D*)fin->Get(hist_name.c_str());

    return h;
}