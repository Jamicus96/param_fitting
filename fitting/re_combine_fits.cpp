#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TTree.h>
#include <sstream>
#include <RAT/DB.hh>
#include <TRandom3.h>
#include <TKey.h>
#include <TObject.h>
#include <TList.h>
#include <TVectorD.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TText.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLatex.h>

#include "fitter.hpp"
#include "fitVars.hpp"
#include "E_systematics.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"
#include "model_Reactor.hpp"
#include "fitting_utils.hpp"


std::vector<double> combine_hists(TH2D* minllHist, const std::vector<TH2D*>& hists, const std::vector<unsigned int>& start_Dm_idx, const std::vector<unsigned int>& end_Dm_idx, const std::vector<unsigned int>& start_th_idx, const std::vector<unsigned int>& end_th_idx);
void read_hists_from_files(const std::vector<std::string>& hists_addresses, std::vector<TH2D*>& hists, std::string hist_name);
std::vector<TH2D*> likelihood_ratio_hists(TH2D* minllHist, double minimisedLikelihood);
void print_to_txt(std::string txt_fileName, TH2D* minllHist, const std::vector<TH1D*>& hists);
void GetFitSpectra(std::vector<TH1D*>& hists, std::string PDFs_address, std::string data_ntuple_address, double Dm21_2, double S_12_2, const bool use_Azimov);


int main(int argv, char** argc) {
    // file args
    std::string PDFs_address = argc[1];
    std::string data_ntuple_address = argc[2];
    std::string out_address = argc[3];
    bool use_Azimov = std::stoi(argc[4]);

    // Rest args come as: hist_address_1 start_Dm_idx_1 end_Dm_idx_1 start_th_idx_1 end_th_idx_1 hist_address_2 start_Dm_idx_2 end_Dm_idx_2 start_th_idx_2 end_th_idx_2 ...
    // These give the histogram and the index limits that it covers
    std::vector<std::string> hists_addresses;
    std::vector<unsigned int> start_Dm_idx;
    std::vector<unsigned int> end_Dm_idx;
    std::vector<unsigned int> start_th_idx;
    std::vector<unsigned int> end_th_idx;
    unsigned int option;
    for (unsigned int i = 5; i < argv; ++i) {
        option = (i - 5) % 5;
        if      (option == 0) hists_addresses.push_back(argc[i]);
        else if (option == 1) start_Dm_idx.push_back(atoi(argc[i]));
        else if (option == 2) end_Dm_idx.push_back(atoi(argc[i]));
        else if (option == 3) start_th_idx.push_back(atoi(argc[i]));
        else                  end_th_idx.push_back(atoi(argc[i]));
    }

    // Read in hists from files
    std::cout << "Reading in hists from files..." << std::endl;
    std::vector<TH2D*> hists;
    read_hists_from_files(hists_addresses, hists, "minllHist");

    // Set up new empty 2d log-likelihood histogram
    TH2D* minllHist = (TH2D*)hists.at(0)->Clone();
    minllHist->Reset("ICES");

    // Combine hist parts into one hist
    std::cout << "Fitting spectra to dataset..." << std::endl;
    std::vector<double> min_vals = combine_hists(minllHist, hists, start_Dm_idx, end_Dm_idx, start_th_idx, end_th_idx);

    std::cout << "min_ll = " << min_vals.at(0) << ", at Dm_21^2 = " << min_vals.at(1) << " and theta_12 = " << min_vals.at(2) << std::endl;

    std::vector<TH2D*> new_hists = likelihood_ratio_hists(minllHist, min_vals.at(0));

    std::vector<TH1D*> spectra;
    GetFitSpectra(spectra, PDFs_address, data_ntuple_address, min_vals.at(1), min_vals.at(2), use_Azimov);

    // Print hist to text file too
    std::string txt_fileName = out_address.substr(0, out_address.find_last_of(".")) + ".txt";
    print_to_txt(txt_fileName, new_hists.at(0), spectra);

    // Write hist to file and close
    TFile *outroot = new TFile(out_address.c_str(), "RECREATE");
    minllHist->Write();
    new_hists.at(0)->Write();
    new_hists.at(1)->Write();
    for (unsigned int i = 0; i < spectra.size(); ++i) {
        spectra.at(i)->Write();
    }
    outroot->Write();
    outroot->Close();
    delete(outroot);

    return 0;
}


std::vector<double> combine_hists(TH2D* minllHist, const std::vector<TH2D*>& hists, const std::vector<unsigned int>& start_Dm_idx, const std::vector<unsigned int>& end_Dm_idx, const std::vector<unsigned int>& start_th_idx, const std::vector<unsigned int>& end_th_idx) {

    double min_ll = 99999.;
    double min_Dm21;
    double min_Theta12;
    double content;

    unsigned int min_hist;
    bool foudZero;

    std::cout << "Looping over hists..." << std::endl;
    for (unsigned int n = 0; n < hists.size(); ++n) {
        std::cout << "Looping over bins in hist " << n << "..." << std::endl;
        foudZero = false;
        for (unsigned int i = start_th_idx.at(n); i <= end_th_idx.at(n); ++i) {
            for (unsigned int j = start_Dm_idx.at(n); j <= end_Dm_idx.at(n); ++j) {
                content = hists.at(n)->GetBinContent(i + 1, j + 1);
                minllHist->SetBinContent(i + 1, j + 1, content);

                if ((content < min_ll) && (content != 0.0)) {
                    min_ll = content;
                    min_hist = n;
                    min_Theta12 = minllHist->GetXaxis()->GetBinCenter(i + 1);
                    min_Dm21 = minllHist->GetYaxis()->GetBinCenter(j + 1);
                }
                if (content == 0.0) foudZero = true;
            }
        }
        if (foudZero) std::cout << "min_ll is zero found in hist " << n << std::endl;
    }
    std::cout << "minimum found in hist " << min_hist << std::endl;

    return {min_ll, min_Dm21, min_Theta12};
}


/**
 * @brief Read in all the 2-D histograms from list of files
 * 
 * @param hists_addresses 
 * @param hists 
 * @param hist_name 
 */
void read_hists_from_files(const std::vector<std::string>& hists_addresses, std::vector<TH2D*>& hists, std::string hist_name) {

    for (unsigned int n = 0; n < hists_addresses.size(); ++n) {
        TFile *fin = TFile::Open(hists_addresses.at(n).c_str());
        if (!fin->IsOpen()) {
            std::cout << "Cannot open input file." << std::endl;
            exit(1);
        }

        // Get histogram
        hists.push_back((TH2D*)fin->Get(hist_name.c_str()));
    }
}


std::vector<TH2D*> likelihood_ratio_hists(TH2D* minllHist, double minimisedLikelihood) {
    // now do likelihood ratio test on maximal likelihood

    std::cout << "Performing ratio test" << std::endl;

    // int deltaM21BestFitBin, theta12BestFitBin, minimisedLikelihoodBin;
    // minllHist->GetBinXYZ(minllHist->GetMinimumBin(), deltaM21BestFitBin, theta12BestFitBin, minimisedLikelihoodBin);
    // double minimisedLikelihood = minllHist->GetBinContent(minllHist->GetMinimumBin());
    double maximisedLikelihood = minllHist->GetBinContent(minllHist->GetMaximumBin());

    TH2D* sigmaMinHist = (TH2D*)minllHist->Clone();
    TH2D* sigmaMaxHist = (TH2D*)minllHist->Clone();
    sigmaMinHist->SetName("sigmaMinHist"); sigmaMinHist->SetTitle("2(l - l_min)");
    sigmaMaxHist->SetName("sigmaMaxHist"); sigmaMaxHist->SetTitle("2(l - l_max)");

    double deltaMBestFitMax, thetaBestFitMax, deltaMBestFitMin, thetaBestFitMin;

    for (unsigned int i = 1; i < minllHist->GetNbinsX() + 1; i++) {
        for (unsigned int j = 1; j < minllHist->GetNbinsY() + 1; j++) {
            sigmaMinHist->SetBinContent(i, j, 2.0 * (minllHist->GetBinContent(i, j) - minimisedLikelihood));
            sigmaMaxHist->SetBinContent(i, j, 2.0 * (minllHist->GetBinContent(i, j) - maximisedLikelihood));
            if (minimisedLikelihood == minllHist->GetBinContent(i, j)) {
                thetaBestFitMin = minllHist->GetXaxis()->GetBinCenter(i);
                deltaMBestFitMin = minllHist->GetYaxis()->GetBinCenter(j);
            }
            if (maximisedLikelihood == minllHist->GetBinContent(i, j)) {
                thetaBestFitMax = minllHist->GetXaxis()->GetBinCenter(i);
                deltaMBestFitMax = minllHist->GetYaxis()->GetBinCenter(j);
            }
        }
    }

    std::cout << "The minimum likelihood is " << minimisedLikelihood << " at deltaM21 = " << deltaMBestFitMin << " and theta12 = " << thetaBestFitMin << std::endl;
    std::cout << "The maximum likelihood is " << maximisedLikelihood << " at deltaM21 = " << deltaMBestFitMax << " and theta12 = " << thetaBestFitMax << std::endl;

    return {sigmaMinHist, sigmaMaxHist};
}

/**
 * @brief Prints hist values to text file, in format:
 * NA Dm21_0 Dm21_1 Dm21_2 ... Dm21_N
 * theta12_0 minLL_00 minLL_01 minLL_02 ... minLL_0N
 * theta12_1 minLL_10 minLL_11 minLL_12 ... minLL_1N
 * theta12_2 minLL_20 minLL_21 minLL_22 ... minLL_2N
 * ... ... ... ... ... ...
 * theta12_N minLL_N0 minLL_N1 minLL_N2 ... minLL_NN
 * 
 * @param txt_fileName 
 * @param minllHist 
 */
void print_to_txt(std::string txt_fileName, TH2D* minllHist, const std::vector<TH1D*>& hists) {
    std::ofstream datafile;
    datafile.open(txt_fileName.c_str(), std::ios::trunc);

    datafile << "# Delta log-likelihood:" << std::endl;

    double Dm21;
    double theta12;
    double minLL;
    datafile << "NA";
    for (unsigned int j = 1; j < minllHist->GetNbinsY() + 1; ++j) {
        Dm21 = minllHist->GetYaxis()->GetBinCenter(j);
        datafile << " " << Dm21;
    }
    datafile << std::endl;

    for (unsigned int i = 1; i < minllHist->GetNbinsX() + 1; ++i) {
        theta12 = minllHist->GetXaxis()->GetBinCenter(i);
        datafile << theta12;
        for (unsigned int j = 1; j < minllHist->GetNbinsY() + 1; ++j) {
            minLL = minllHist->GetBinContent(i, j);  // was: GetBinContent(i + 1, j + 1)
            datafile << " " << minLL;
        }
        datafile << std::endl;
    }

    datafile << "# Spectra:" << std::endl;

    datafile << "NA";
    for (unsigned int j = 1; j < hists.at(0)->GetNbinsX() + 1; ++j) {
        datafile << " " << hists.at(0)->GetBinCenter(j);
    }
    datafile << std::endl;
    for (unsigned int i = 0; i < hists.size(); ++i) {
        datafile << hists.at(i)->GetName();
        for (unsigned int j = 1; j < hists.at(i)->GetNbinsX() + 1; ++j) {
            datafile << " " << hists.at(i)->GetBinContent(j);
        }
        datafile << std::endl;
    }
}


void GetFitSpectra(std::vector<TH1D*>& spectra, std::string PDFs_address, std::string data_ntuple_address, double Dm21_2, double S_12_2, const bool use_Azimov) {

    // Get DB
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DB* db = RAT::DB::Get();
    db->LoadDefaults();

    // Get oscillation constants
    std::cout << "Getting oscillation parameters..." << std::endl;
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");
    const double fDmSqr32 = linkdb->GetD("deltamsqr32");
    const double fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

    // Create fitter object
    if (use_Azimov) {
        create_fitter(PDFs_address, Dm21_2, fDmSqr32, pow(sin(S_12_2  * TMath::Pi() / 180.), 2), fSSqrTheta13, db);
    } else {
        TFile* DataFile = TFile::Open(data_ntuple_address.c_str());
        TTree* DataInfo = (TTree *) DataFile->Get("prompt");
        create_fitter(PDFs_address, Dm21_2, fDmSqr32, pow(sin(S_12_2  * TMath::Pi() / 180.), 2), fSSqrTheta13, db, DataInfo);
    }
    Fitter* antinuFitter = Fitter::GetInstance();
    FitVars* Vars = FitVars::GetInstance();

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;

    double ll = antinuFitter->fit_models();
    std::cout << "ll = " << ll <<std::endl;

    // Add spectra to list
    antinuFitter->GetAllSpectra(spectra);

    // Add data hist to list
    if (use_Azimov) spectra.push_back(antinuFitter->DataHist());

    // Print some extra things
    std::cout << "Covariance matrix:" << std::endl;
    std::cout << "NA";
    for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
        if (Vars->isConstant(i)) continue;
        std::cout << "\t" << Vars->name(i);
    }
    std::cout << std::endl;
    unsigned int k = 0, l = 0;
    for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
        if (Vars->isConstant(i)) continue;
        std::cout << Vars->name(i) << "\t" << antinuFitter->GetCovarianceMatrixElement(k, 0);
        l = 0;
        for (unsigned int j = 1; j < Vars->GetNumVars(); ++j) {
            if (Vars->isConstant(j)) continue;
            std::cout << "\t" << antinuFitter->GetCovarianceMatrixElement(k, l);
            ++l;
        }
        std::cout << std::endl;
        ++k;
    }
    Double_t eplus, eminus, eparab, globcc;
    std::cout << "GetErrors:" << std::endl;
    for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
        antinuFitter->GetErrors(i, eplus, eminus, eparab, globcc);
        std::cout << Vars->name(i) << ": eplus = " << eplus << ", eminus = " << eminus << ", eparab = "
                  << eparab << ", globcc = " << globcc << std::endl;
    }

    // /* ~~~~~~~~~~~~~ Try an overall fit ~~~~~~~~~~~~~ */

    // const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    // const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");

    // double Dm212err = 0.2E-5;
    // double s122err = 0.013;
    // double min_s122 = fSSqrTheta12 - 3.*s122err;
    // double max_s122 = fSSqrTheta12 + 3.*s122err;
    // if (min_s122 < 0.0) min_s122 = 0.0;
    // if (max_s122 > 1.0) max_s122 = 1.0;
    // antinuFitter->resetVar("deltamsqr21", fDmSqr21, Dm212err, fDmSqr21 - 3.*Dm212err, fDmSqr21 + 3.*Dm212err, false);
    // antinuFitter->resetVar("sinsqrtheta12", fSSqrTheta12, s122err, min_s122, max_s122, false);

    // Reactor* ReactorMod = Reactor::GetInstance();
    // geoNu* geoNuMod = geoNu::GetInstance();
    // ReactorMod->hold_osc_params_const(false);
    // geoNuMod->hold_osc_params_const(false);

    // // Do fitting for a range of values, summarised in 2-D hist
    // std::cout << "Doing full fit..." << std::endl;

    // ll = antinuFitter->fit_models();
    // std::cout << "ll = " << ll <<std::endl;

    // // Add spectra to list
    // antinuFitter->GetAllSpectra(spectra);

    // // Add data hist to list
    // spectra.push_back(antinuFitter->DataHist());

    // // Print some extra things
    // std::cout << "Covariance matrix:" << std::endl;
    // for (unsigned int i = 0; i < Vars->GetNumVars()-4; ++i) {
    //     std::cout << antinuFitter->GetCovarianceMatrixElement(i, 0);
    //     for (unsigned int j = 1; j < Vars->GetNumVars()-4; ++j) {
    //         std::cout << "\t" << antinuFitter->GetCovarianceMatrixElement(i, j);
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "GetErrors:" << std::endl;
    // for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
    //     antinuFitter->GetErrors(i, eplus, eminus, eparab, globcc);
    //     std::cout << Vars->name(i) << ": eplus = " << eplus << ", eminus = " << eminus << ", eparab = "
    //               << eparab << ", globcc = " << globcc << std::endl;
    // }
}


TCanvas* ContourList(TH2D* minllHist){

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEFAULT SNO+ SETTINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    TStyle *snoStyle= new TStyle("snoplus","SNO+ plots style for publications");

    // use plain black on white colors
    snoStyle->SetFrameBorderMode(0);
    snoStyle->SetCanvasBorderMode(0);
    snoStyle->SetPadBorderMode(0);
    snoStyle->SetPadBorderSize(0);
    snoStyle->SetPadColor(0);
    snoStyle->SetCanvasColor(0);
    snoStyle->SetTitleColor(0);
    snoStyle->SetStatColor(0);
    snoStyle->SetFillColor(0);

    // use bold lines 
    snoStyle->SetHistLineWidth(2);
    snoStyle->SetLineWidth(2);

    // no title, stats box or fit as default
    snoStyle->SetOptTitle(0);
    //snoStyle->SetOptStat(0);
    //snoStyle->SetOptFit(0);

    // postscript dashes
    snoStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

    // text style and size
    //snoStyle->SetTextFont(132);
    //snoStyle->SetTextSize(0.24);
    snoStyle->SetLabelOffset(0.01,"x");
    snoStyle->SetTickLength(0.015,"x");
    snoStyle->SetTitleOffset(1.5,"x");
    snoStyle->SetLabelOffset(0.01,"y");
    snoStyle->SetTickLength(0.015,"y");
    snoStyle->SetTitleOffset(1.5,"y");
    snoStyle->SetLabelOffset(0.01,"z");
    snoStyle->SetTickLength(0.015,"z");
    snoStyle->SetTitleOffset(1.5,"z");
    snoStyle->SetLabelFont(132,"x");
    snoStyle->SetLabelFont(132,"y");
    snoStyle->SetLabelFont(132,"z");
    snoStyle->SetTitleFont(132,"x");
    snoStyle->SetTitleFont(132,"y");
    snoStyle->SetTitleFont(132,"z");
    snoStyle->SetLabelSize(0.04,"x");
    snoStyle->SetTitleSize(0.05,"x");
    snoStyle->SetTitleColor(1,"x");
    snoStyle->SetLabelSize(0.04,"y");
    snoStyle->SetTitleSize(0.05,"y");
    snoStyle->SetTitleColor(1,"y");
    snoStyle->SetLabelSize(0.04,"z");
    snoStyle->SetTitleSize(0.05,"z");
    snoStyle->SetTitleColor(1,"z");
    snoStyle->SetPadTickX(1);
    snoStyle->SetPadTickY(1);

    // AXIS OFFSETS
    snoStyle->SetTitleOffset(0.8,"x");
    snoStyle->SetTitleOffset(0.8,"y");
    snoStyle->SetTitleOffset(0.8,"z");

    // Legends
    snoStyle->SetLegendBorderSize(0);
    snoStyle->SetLegendFont(132);
    snoStyle->SetLegendFillColor(0);

    // graphs - set default marker to cross, rather than .
    snoStyle->SetMarkerStyle(21);  // filled square not .

    // SNO+ Preliminary label
    snoStyle->SetTextFont(132);
    snoStyle->SetTextSize(0.06);

    gROOT->SetStyle("snoplus");


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Set up canvas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
 
    TCanvas* c = new TCanvas("c","Contour List",0,0,600,600);
    
    Double_t contours[6];
    contours[0] = -0.7;
    contours[1] = -0.5;
    contours[2] = -0.1;
    contours[3] =  0.1;
    contours[4] =  0.4;
    contours[5] =  0.8;
    
    minllHist->SetContour(6, contours);
    
    // Draw contours as filled regions, and Save points
    minllHist->Draw("CONT Z LIST");
    c->Update(); // Needed to force the plotting and retrieve the contours in TGraphs
    
    // Get Contours
    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    
    if (!conts){
        printf("*** No Contours Were Extracted!\n");
        return nullptr;
    }
    
    TList* contLevel = nullptr;
    TGraph* curv     = nullptr;
    TGraph* gc       = nullptr;
    
    Int_t nGraphs    = 0;
    Int_t TotalConts = conts->GetSize();
    
    printf("TotalConts = %d\n", TotalConts);
    
    TCanvas* c1 = new TCanvas("c1","Contour List",610,0,600,600);
    c1->SetTopMargin(0.15);
    TH2F *hr = new TH2F("hr",
    "#splitline{Negative contours are returned first (highest to lowest). Positive contours are returned from}{lowest to highest. On this plot Negative contours are drawn in red and positive contours in blue.}",
    2, -2, 2, 2, 0, 6.5);
    
    hr->Draw();
    Double_t xval0, yval0, zval0;
    TLatex l;
    l.SetTextSize(0.03);
    char val[20];
    
    for(unsigned int i = 0; i < TotalConts; i++){
        contLevel = (TList*)conts->At(i);
        if (i<3) zval0 = contours[2-i];
        else     zval0 = contours[i];
        printf("Z-Level Passed in as:  Z = %f\n", zval0);
    
        // Get first graph from list on curves on this level
        curv = (TGraph*)contLevel->First();
        for(unsigned int j = 0; j < contLevel->GetSize(); j++){
            curv->GetPoint(0, xval0, yval0);
            if (zval0<0) curv->SetLineColor(kRed);
            if (zval0>0) curv->SetLineColor(kBlue);
            nGraphs ++;
            printf("\tGraph: %d  -- %d Elements\n", nGraphs,curv->GetN());
    
            // Draw clones of the graphs to avoid deletions in case the 1st
            // pad is redrawn.
            gc = (TGraph*)curv->Clone();
            gc->Draw("C");
    
            sprintf(val,"%g",zval0);
            l.DrawLatex(xval0,yval0,val);
            curv = (TGraph*)contLevel->After(curv); // Get Next graph
        }
    }
    c1->Update();
    printf("\n\n\tExtracted %d Contours and %d Graphs \n", TotalConts, nGraphs );
    gStyle->SetTitleW(0.);
    gStyle->SetTitleH(0.);
    return c1;
}

