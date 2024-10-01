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
void print_to_txt(std::string txt_fileName, TH2D* minllHist, const std::vector<TH1D*>& hists, const std::vector<double>& data);
void GetFitSpectra(std::vector<TH1D*>& hists, std::vector<double>& data, TFile *fin, TFile* DataFile, double Dm21_2, double S_12_2, const bool use_Azimov);
void overallFit(std::string txt_fileName, const bool use_Azimov);


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

    std::cout << "min_ll = " << min_vals.at(0) << ", at Dm_21^2 = " << min_vals.at(1) << " and s_12_2 = " << min_vals.at(2) << std::endl;

    std::vector<TH2D*> new_hists = likelihood_ratio_hists(minllHist, min_vals.at(0));

    std::vector<TH1D*> spectra;
    std::vector<double> data;
    TFile* DataFile = TFile::Open(data_ntuple_address.c_str());
    TFile *fin = TFile::Open(PDFs_address.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file." << std::endl;
        exit(1);
    }
    GetFitSpectra(spectra, data, fin, DataFile, min_vals.at(1), min_vals.at(2), use_Azimov);

    std::cout << "Ouside GetFitSpectra()" << std::endl;

    // Print hist to text file too
    std::string txt_fileName = out_address.substr(0, out_address.find_last_of(".")) + ".txt";
    std::cout << "txt_fileName = " << txt_fileName << std::endl;
    print_to_txt(txt_fileName, new_hists.at(0), spectra, data);

    // Write hist to file and close
    std::cout << "Writing to root file" << std::endl;
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

    // Try overall fit
    std::cout << "Trying overall fit" << std::endl;
    overallFit(txt_fileName, use_Azimov);

    // DataFile->Close();
    // fin->Close();

    return 0;
}


std::vector<double> combine_hists(TH2D* minllHist, const std::vector<TH2D*>& hists, const std::vector<unsigned int>& start_Dm_idx, const std::vector<unsigned int>& end_Dm_idx, const std::vector<unsigned int>& start_th_idx, const std::vector<unsigned int>& end_th_idx) {

    double min_ll = 99999.;
    double min_Dm21;
    double min_s_12_2;
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
                    min_s_12_2 = minllHist->GetXaxis()->GetBinCenter(i + 1);
                    min_Dm21 = minllHist->GetYaxis()->GetBinCenter(j + 1);
                }
                if (content == 0.0) foudZero = true;
            }
        }
        if (foudZero) std::cout << "min_ll is zero found in hist " << n << std::endl;
    }
    std::cout << "minimum found in hist " << min_hist << std::endl;

    return {min_ll, min_Dm21, min_s_12_2};
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
void print_to_txt(std::string txt_fileName, TH2D* minllHist, const std::vector<TH1D*>& hists, const std::vector<double>& data) {
    std::cout << "start of print_to_txt()" << std::endl;
    std::ofstream datafile;
    std::cout << "1" << std::endl;
    datafile.open(txt_fileName.c_str(), std::ios::trunc);
    std::cout << "2" << std::endl;

    datafile << "# Delta log-likelihood:" << std::endl;

    std::cout << "3" << std::endl;

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

    std::cout << "4" << std::endl;

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

    std::cout << "5" << std::endl;

    if (data.size() > 0) {
        datafile << "# Data:" << std::endl;
        datafile << data.at(0);
        for (unsigned int i = 1; i < data.size(); ++i) {
            datafile << " " << data.at(i);
        }
        datafile << std::endl;
    }

    std::cout << "6" << std::endl;
}


void GetFitSpectra(std::vector<TH1D*>& spectra, std::vector<double>& data, TFile *fin, TFile* DataFile, double Dm21_2, double S_12_2, const bool use_Azimov) {

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
        create_fitter(fin, Dm21_2, fDmSqr32, S_12_2, fSSqrTheta13, db);
    } else {
        // TFile* DataFile = TFile::Open(data_ntuple_address.c_str());
        TTree* DataInfo = (TTree *) DataFile->Get("prompt");
        create_fitter(fin, Dm21_2, fDmSqr32, S_12_2, fSSqrTheta13, db, DataInfo);
        // DataFile->Close();
    }
    Fitter* antinuFitter = Fitter::GetInstance();
    FitVars* Vars = FitVars::GetInstance();

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "Fitting spectra to dataset..." << std::endl;

    antinuFitter->fitErrors(true);
    double ll = antinuFitter->fit_models();
    std::cout << "ll = " << ll <<std::endl;

    // Add spectra to list
    antinuFitter->GetAllSpectra(spectra);

    // Add data hist to list
    if (use_Azimov) spectra.push_back(antinuFitter->DataHist());
    else antinuFitter->DataNtuple(data);

    // Print some extra things
    std::vector<std::string> VarNames;
    std::vector<std::vector<double>> CovMat;
    std::vector<std::vector<double>> CorrMat;
    unsigned int k = 0, l = 0;
    for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
        if (Vars->IsConstant(i)) continue;
        VarNames.push_back(Vars->name(i));
        CovMat.push_back({});
        CorrMat.push_back({});
        l = 0;
        for (unsigned int j = 0; j < Vars->GetNumVars(); ++j) {
            if (Vars->IsConstant(j)) continue;
            CovMat.at(k).push_back(antinuFitter->GetCovarianceMatrixElement(k, l));
            CorrMat.at(k).push_back(0);
            ++l;
        }
        ++k;
    }

    std::cout << "\nCovariance matrix:" << std::endl;
    std::cout << "NA";
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << "\t" << VarNames.at(i);
    }
    std::cout << std::endl;
    double denom;
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << VarNames.at(i);
        for (unsigned int j = 0; j < VarNames.size(); ++j) {
            std::cout << "\t" << CovMat.at(i).at(j);
            denom = CovMat.at(i).at(i) * CovMat.at(j).at(j);
            if (denom == 0) CorrMat.at(i).at(j) = 0;
            else CorrMat.at(i).at(j) = CovMat.at(i).at(j) / sqrt(fabs(CovMat.at(i).at(i) * CovMat.at(j).at(j)));
        }
        std::cout << std::endl;
    }

    // char* name;
    // Double_t value, verr, vlow, vhigh;
    // std::cout << "\nGetParameters:" << std::endl;
    // for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
    //     antinuFitter->GetParameter(i, name, value, verr, vlow, vhigh);
    //     std::cout << Vars->name(i) << ": name = " << name << ", value = " << value << ", verr = "
    //               << verr << ", vlow = " << vlow << ", vhigh = " << vhigh << std::endl;
    // }

    Double_t eplus, eminus, eparab, globcc;
    std::cout << "\nGetErrors:" << std::endl;
    for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
        antinuFitter->GetErrors(i, eplus, eminus, eparab, globcc);
        std::cout << Vars->name(i) << ": eplus = " << eplus << ", eminus = " << eminus << ", eparab = "
                  << eparab << ", globcc = " << globcc << std::endl;
    }

    std::cout << "\nCorrelation matrix:" << std::endl;
    std::cout << "NA";
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << "\t" << VarNames.at(i);
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << VarNames.at(i);
        for (unsigned int j = 0; j < VarNames.size(); ++j) {
            std::cout << "\t" << CorrMat.at(i).at(j);
        }
        std::cout << std::endl;
    }

    std::cout << "Reached end of GetFitSpectra()" << std::endl;
    return;
}


void overallFit(std::string txt_fileName, const bool use_Azimov) {
    
    Fitter* antinuFitter = Fitter::GetInstance();
    FitVars* Vars = FitVars::GetInstance();
    Reactor* ReactorMod = Reactor::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();

    // Get DB
    RAT::DB::Get()->SetAirplaneModeStatus(true);
    RAT::DB* db = RAT::DB::Get();
    db->LoadDefaults();

    // Get oscillation constants
    std::cout << "Getting oscillation parameters..." << std::endl;
    RAT::DBLinkPtr linkdb = db->GetLink("OSCILLATIONS");

    const double fDmSqr21 = linkdb->GetD("deltamsqr21");
    const double fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");

    double Dm212err = 0.18E-5;
    double s122err = 0.013;
    antinuFitter->resetVar("deltamsqr21", fDmSqr21, Dm212err, 4.E-5, 12.E-5, true, true);
    antinuFitter->resetVar("sinsqrtheta12", fSSqrTheta12, s122err, 0.05, 0.95, true, true);

    ReactorMod->hold_osc_params_const(false);
    geoNuMod->hold_osc_params_const(false);

    // Do fitting for a range of values, summarised in 2-D hist
    std::cout << "\n\nDoing full fit...\n" << std::endl;

    double ll = antinuFitter->fit_models();
    std::cout << "ll = " << ll <<std::endl;

    // Add spectra to list
    std::vector<TH1D*> spectra;
    antinuFitter->GetAllSpectra(spectra);

    // Print some extra things
    std::vector<std::string> VarNames;
    std::vector<std::vector<double>> CovMat;
    std::vector<std::vector<double>> CorrMat;
    unsigned int k = 0, l = 0;
    for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
        if (Vars->IsConstant(i)) continue;
        VarNames.push_back(Vars->name(i));
        CovMat.push_back({});
        CorrMat.push_back({});
        l = 0;
        for (unsigned int j = 0; j < Vars->GetNumVars(); ++j) {
            if (Vars->IsConstant(j)) continue;
            CovMat.at(k).push_back(antinuFitter->GetCovarianceMatrixElement(k, l));
            CorrMat.at(k).push_back(0);
            ++l;
        }
        ++k;
    }

    std::cout << "Covariance matrix:" << std::endl;
    std::cout << "NA";
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << "\t" << VarNames.at(i);
    }
    std::cout << std::endl;
    double denom;
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << VarNames.at(i);
        for (unsigned int j = 0; j < VarNames.size(); ++j) {
            std::cout << "\t" << CovMat.at(i).at(j);
            denom = CovMat.at(i).at(i) * CovMat.at(j).at(j);
            if (denom == 0) CorrMat.at(i).at(j) = 0;
            else CorrMat.at(i).at(j) = CovMat.at(i).at(j) / sqrt(fabs(CovMat.at(i).at(i) * CovMat.at(j).at(j)));
        }
        std::cout << std::endl;
    }

    Double_t eplus, eminus, eparab, globcc;
    std::cout << "GetErrors:" << std::endl;
    for (unsigned int i = 0; i < Vars->GetNumVars(); ++i) {
        antinuFitter->GetErrors(i, eplus, eminus, eparab, globcc);
        std::cout << Vars->name(i) << ": eplus = " << eplus << ", eminus = " << eminus << ", eparab = "
                  << eparab << ", globcc = " << globcc << std::endl;
    }

    std::cout << "Correlation matrix:" << std::endl;
    std::cout << "NA";
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << "\t" << VarNames.at(i);
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < VarNames.size(); ++i) {
        std::cout << VarNames.at(i);
        for (unsigned int j = 0; j < VarNames.size(); ++j) {
            std::cout << "\t" << CorrMat.at(i).at(j);
        }
        std::cout << std::endl;
    }

    // print to file
    std::ofstream datafile;
    datafile.open(txt_fileName.c_str(), std::ios::app);

    datafile << "# Global Fit Spectra:" << std::endl;

    datafile << "NA";
    for (unsigned int j = 1; j < spectra.at(0)->GetNbinsX() + 1; ++j) {
        datafile << " " << spectra.at(0)->GetBinCenter(j);
    }
    datafile << std::endl;
    for (unsigned int i = 0; i < spectra.size(); ++i) {
        std::cout << "spectra.at(" << i << ")" << std::endl;
        datafile << spectra.at(i)->GetName();
        for (unsigned int j = 1; j < spectra.at(i)->GetNbinsX() + 1; ++j) {
            datafile << " " << spectra.at(i)->GetBinContent(j);
        }
        datafile << std::endl;
    }
}