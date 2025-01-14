// header guard:
#ifndef fitter
#define fitter

// include
#include <TVirtualFitter.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TH1.h>
#include <iostream>
#include <algorithm>
#include "fitVars.hpp"
#include "model_Reactor.hpp"
#include "model_alphaN.hpp"
#include "model_geoNu.hpp"
#include "model.hpp"

class Fitter {
    private:
        static Fitter *FitterInstance_;

        static TVirtualFitter* minuit;
        static double arglist[2];

        static TH1D *data, *tot_fitModel;
        static std::vector<double> data_E;
        static double Emin, dE;
        static bool dataIsSet;
        static bool fitErrs;

        static unsigned int iBinMin, iBinMax;
        static bool setLims, binnedData, init_tot_fitModel;

        // For minuit
        static double minfuncOut, edm, errdef;
        static int nvpar, nparx;

        // class Reactor* reac;
    
    protected:
        // Constructors/desctructors
        Fitter() {
            #ifdef antinuDEBUG
                std::cout << "[Fitter::Fitter]: Creating singleton." << std::endl;
            #endif
            Initialise();
        }
        ~Fitter() {}

    public:
        // Should not be cloneable.
        Fitter(Fitter &other) = delete;

        // Should not be assignable.
        void operator=(const Fitter&) = delete;

        /**
         * This is the method that controls the access to the Fitter
         * instance. On the first run, it creates a Fitter object and places it
         * into the field. On subsequent runs, it returns the client existing
         * object stored in the field.
         */
        static Fitter *GetInstance();

        // Takes care of initialising everything for the fitter.
        // Should only be run once, unless it needs to be re-initialised
        static void Initialise();

        // Fitting functions
        static int ExecuteFit();
        static double fit_models();
        static Double_t fitFunc(Double_t* p);
        static void fitErrors(bool flag);

        // (Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag*/)
        static void fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t);
        static void SetFunc();

        double GetFitParameter(const unsigned int parIdx) {return minuit->GetParameter(parIdx);}
        double GetFitParError(const unsigned int parIdx) {return minuit->GetParError(parIdx);}
        double GetFitCovarianceMatrixElement(const unsigned int i, const unsigned int j) {return minuit->GetCovarianceMatrixElement(i, j);}
        static double ExtendedConstrainedLogLikelihood();

        // Get info functions
        static void GetAllSpectra(std::vector<TH1D*>& hists);
        static TH1D* DataHist();
        static std::vector<double> DataNtuple();
        void DataNtuple(std::vector<double>& out_data);
        static void SetData(TH1D* Data);
        static void SetData(TTree* Data);
        static void SetBinLims(const unsigned int Bin_min, const unsigned int Bin_max);

        static void GetParameter(Int_t ipar, char* name, Double_t& value, Double_t& verr, Double_t& vlow, Double_t& vhigh);
        static void GetErrors(Int_t ipar, Double_t& eplus, Double_t& eminus, Double_t& eparab, Double_t& globcc);
        static double GetCovarianceMatrixElement(int i, int j);

        static void resetVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool holdConstant = false, bool isConstrained = true);
        static void InitTotFitModHist(TH1D* exampleHist);
};

// Static methods should be defined outside the class.
Fitter* Fitter::FitterInstance_{nullptr};

TVirtualFitter* Fitter::minuit;
double Fitter::arglist[2];
TH1D *Fitter::data, *Fitter::tot_fitModel;
std::vector<double> Fitter::data_E;
double Fitter::minfuncOut, Fitter::edm, Fitter::errdef;
int Fitter::nvpar, Fitter::nparx;
bool Fitter::dataIsSet = false;
unsigned int Fitter::iBinMin, Fitter::iBinMax;
bool Fitter::setLims = false, Fitter::binnedData = false, Fitter::init_tot_fitModel = false, Fitter::fitErrs = false;
double Fitter::Emin, Fitter::dE;

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
Fitter* Fitter::GetInstance() {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::GetInstance]: Creating Instance." << std::endl;
    #endif
    if (FitterInstance_ == nullptr) {
        FitterInstance_ = new Fitter();
    }
    return FitterInstance_;
}

void Fitter::Initialise() {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::Initialise]: Initialising Minuit." << std::endl;
    #endif

    FitVars* Vars = FitVars::GetInstance();
    if (Vars->GetNumVars() == 0) {
        std::cout << "Initialising miniut with zero variables! Must set up all variables first." << std::endl;
        exit(1);
    }

    //The default minimizer is Minuit, you can also try Minuit2
    // TVirtualFitter::SetDefaultFitter("Minuit2");
    TVirtualFitter::SetDefaultFitter("Minuit");
    minuit = TVirtualFitter::Fitter(0, Vars->GetNumVars());
    
    arglist[0] = 0;
    // set print level
    minuit->ExecuteCommand("SET PRINT", arglist, 2);

    // minimiser settings
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance

    // (void (*fcn)(Int_t&, Double_t*, Double_t& f, Double_t*, Int_t))
    SetFunc();

    for (unsigned int iVar = 0; iVar < Vars->GetNumVars(); ++iVar) {
        // Sets the parameter
        if (Vars->IsConstant(iVar)) minuit->SetParameter(iVar, Vars->name(iVar).c_str(), Vars->val(iVar), 0.0, Vars->val(iVar), Vars->val(iVar));
        else minuit->SetParameter(iVar, Vars->name(iVar).c_str(), Vars->prior(iVar), Vars->err(iVar), Vars->min(iVar), Vars->max(iVar));
        // Binds the parameter's detail's addresses to fitVar's, so that changing the fitVar object changes this automatically [DOESN'T WORK]
        // minuit->GetParameter(iVar, const_cast<char*>(Vars->name(iVar).c_str()), Vars->prior(iVar), Vars->err(iVar), Vars->min(iVar), Vars->max(iVar));
    }
}

void Fitter::fitFunc_minuit_format(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t) {
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::fitFunc_minuit_format]: Running fitFunc_minuit_format()." << std::endl;
    #endif
    fval = fitFunc(p);
}

void Fitter::SetFunc() {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::SetFunc]: Bound minuit to fitFunc_minuit_format()." << std::endl;
    #endif
    minuit->SetFCN(fitFunc_minuit_format);
}

int Fitter::ExecuteFit() {
    int iErr;

    if (fitErrs) {
        arglist[0] = 2.0;  // 0, 1 (default), 2 (strategy = 2 allows Minuit to waste function calls to ensure that all values are precise)                                                                    
        iErr = minuit->ExecuteCommand("SET STR", arglist, 1);
        if (iErr != 0) return iErr;
    }

    // Adjust the error definition for the likelihood fits. 0.5 is the value for negative log-likelihood fit.
    arglist[0] = 0.5;
    iErr = minuit->ExecuteCommand("SET ERR", arglist, 1);
    if (iErr != 0) return iErr;

    // 1000 steps, Simplex method for minimization to find good starting values
    arglist[0] = 1000;
    iErr = minuit->ExecuteCommand("SIMPLEX", arglist, 1);
    if (iErr != 0) return iErr;

    // Try MIGRAD a few times
    iErr = 1;
    unsigned int num_attempts = 0;
    while (iErr) {
        ++num_attempts;
        // 1000 steps, Migrad minimization
        iErr = minuit->ExecuteCommand("MIGRAD", arglist, 1);
        if (num_attempts >= 5) break;
    }
    if (iErr != 0) return iErr;

    // Get fit errors
    if (fitErrs) {
        // Calculate uncertainties with Minos (asymmetric errors)
        iErr = minuit->ExecuteCommand("MINOS", arglist, 1);
        if (iErr != 0) return iErr;

        // Calculate Hesse Matrix (i.e. 2dn derivatives) for parameter uncertainties
        iErr = minuit->ExecuteCommand("HESSE", arglist, 1);
    }
    
    return iErr;
}

double Fitter::fit_models() {
    if (!dataIsSet) {
        std::cout << "Attempting to perform fit with no data!" << std::endl;
        exit(1);
    }

    // minimize
    #ifdef antinuDEBUG
        FitVars* Vars = FitVars::GetInstance();
        std::cout << "[Fitter::fit_models]: Running MIGRAD!" << std::endl;
        std::cout << "[Fitter::fit_models]: s_12^2 = " << Vars->val("sinsqrtheta12") << std::endl;
        std::cout << "[Fitter::fit_models]: Dm21^2 = " << Vars->val("deltamsqr21") << std::endl;
    #endif

    int iErr = 1;
    unsigned int num_attempts = 0;
    while (iErr) {
        ++num_attempts;
        // iErr = minuit->ExecuteCommand("MIGRAD", arglist, 2);
        iErr = ExecuteFit();
        if (num_attempts >= 10) break;
    }

    // Try re-initialising and running again a few times
    if (iErr != 0) {
        Initialise();
        num_attempts = 0;
        while (iErr) {
            ++num_attempts;
            // iErr = minuit->ExecuteCommand("MIGRAD", arglist, 2);
            iErr = ExecuteFit();
            if (num_attempts >= 5) break;
        }
    }

    //get result
    minuit->GetStats(minfuncOut, edm, errdef, nvpar, nparx);

    if (iErr != 0) minfuncOut = 0;
    if (std::isnan(edm)) Initialise(); // Re-initialise fitter for next use, otherwise it always fails

    std::cout << "[Fitter::fit_models]: minfuncOut = " << minfuncOut << ", edm = " << edm << ", errdef = " << errdef << ", iErr = " << iErr << ", num_attempts = " << num_attempts << std::endl;
    return minfuncOut;
}

void Fitter::fitErrors(bool flag) {fitErrs = flag;}

Double_t Fitter::fitFunc(Double_t* p) {
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::fitFunc]: Running fitFunc()." << std::endl;
    #endif
    // Reset total model spectrum to zero [[INITIALISE?]]
    if (init_tot_fitModel) tot_fitModel->Reset("ICES");

    // Compute the total spectrum from all the models, using parameters p
    FitVars* Vars = FitVars::GetInstance();
    Reactor* ReactorMod = Reactor::GetInstance();
    alphaN* alphaNMod = alphaN::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();
    Model* BasicMods = Model::GetInstance();

    if (!(ReactorMod->IsInit() || alphaNMod->IsInit() || geoNuMod->IsInit() || BasicMods->GetNumMods() > 0)) {
        std::cout << "Attempting to perform fit with no models!" << std::endl;
        exit(1);
    }

    Vars->GetVarValues(p);

    if (ReactorMod->IsInit()) {
        ReactorMod->compute_spec();
        if (!init_tot_fitModel) {
            InitTotFitModHist(ReactorMod->GetModelEsys());
            tot_fitModel->Reset("ICES");
        }
        tot_fitModel->Add(ReactorMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added ReactorMod, tot_fitModel integral = " << tot_fitModel->Integral(iBinMin, iBinMax) << std::endl;
        #endif
    }
    if (alphaNMod->IsInit()) {
        alphaNMod->compute_spec();
        if (!init_tot_fitModel) {
            InitTotFitModHist(alphaNMod->GetModelEsys());
            tot_fitModel->Reset("ICES");
        }
        tot_fitModel->Add(alphaNMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added alphaNMod, tot_fitModel integral = " << tot_fitModel->Integral(iBinMin, iBinMax) << std::endl;
        #endif
    }
    if (geoNuMod->IsInit()) {
        geoNuMod->compute_spec();
        if (!init_tot_fitModel) {
            InitTotFitModHist(geoNuMod->GetModelEsys());
            tot_fitModel->Reset("ICES");
        }
        tot_fitModel->Add(geoNuMod->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added geoNuMod, tot_fitModel integral = " << tot_fitModel->Integral(iBinMin, iBinMax) << std::endl;
        #endif
    }
    for (unsigned int iMod = 0; iMod < BasicMods->GetNumMods(); ++iMod) {
        BasicMods->compute_spec(iMod);
        if (!init_tot_fitModel) {
            InitTotFitModHist(BasicMods->GetModelEsys(iMod));
            tot_fitModel->Reset("ICES");
        }
        tot_fitModel->Add(BasicMods->GetModelEsys(iMod));
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added BasicMods[" << iMod << "], tot_fitModel integral = " << tot_fitModel->Integral(iBinMin, iBinMax) << std::endl;
        #endif
    }

    // Compute test statistic (minus, since minuit minizes, but we want max)
    return - ExtendedConstrainedLogLikelihood();
}

double Fitter::ExtendedConstrainedLogLikelihood() {
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: Running ExtendedConstrainedLogLikelihood()." << std::endl;
    #endif
    // Model PDFs have already been scaled by their respective norms, and added together
    // Assume data and models have the same binning (could generalise at some point)
    double logL = 0, binContent;

    if (binnedData) {
        // Data is a histogram, same binning as model PDFs
        for (unsigned int ibin = iBinMin; ibin < iBinMax+1; ++ibin) {
            binContent = tot_fitModel->GetBinContent(ibin);
            if (binContent > 0) logL += - binContent + data->GetBinContent(ibin) * log(binContent);
        }

    } else {
        // Data is an ntuple TTree
        for (unsigned int ibin = iBinMin; ibin < iBinMax; ++ibin) {
            binContent = tot_fitModel->GetBinContent(ibin);
            if (binContent > 0) logL += -binContent;
        }
        int E_bin;
        // Loop through all data events:
        for (unsigned int a = 0; a < data_E.size(); ++a) {
            E_bin = int((data_E.at(a) - Emin) / dE + 0.5) + 1;
            if (E_bin <= iBinMin || E_bin >= iBinMax) continue;
            binContent = tot_fitModel->GetBinContent(E_bin);
            if (binContent > 0) logL += log(binContent);
            #ifdef SUPER_DEBUG
                std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: data_E.at(" << a << ") = " << data_E.at(a) << ", E_bin = " << E_bin
                          << ", tot_fitModel->GetBinCenter(E_bin) = " << tot_fitModel->GetBinCenter(E_bin) << ", dE = " << dE << ", Emin = " << Emin
                          << ", tot_fitModel->GetBinContent(E_bin) = " << tot_fitModel->GetBinContent(E_bin) << std::endl;
            #endif
        }
    }
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: +PDFs logL = " << logL << std::endl;
    #endif

    // Add in constraints
    FitVars* Vars = FitVars::GetInstance();
    for (unsigned int iVar = 0; iVar < Vars->GetNumVars(); ++ iVar) {
        // If variable is being held constant, assume it's equal to its prior most likely value
        if (Vars->IsConstrained(iVar) && Vars->err(iVar) != 0.0) {
            // Add gaussian constraint
            logL -= 0.5 * (Vars->val(iVar) - Vars->prior(iVar))*(Vars->val(iVar) - Vars->prior(iVar)) / (Vars->err(iVar)*Vars->err(iVar));
        }
    }
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: +constraints logL = " << logL << std::endl;
    #endif

    return logL;
}

void Fitter::GetAllSpectra(std::vector<TH1D*>& hists) {

    TH1D* total_hist;

    Reactor* ReactorMod = Reactor::GetInstance();
    alphaN* alphaNMod = alphaN::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();
    Model* BasicMods = Model::GetInstance();

    if (!(ReactorMod->IsInit() || alphaNMod->IsInit() || geoNuMod->IsInit() || BasicMods->GetNumMods() > 0)) {
        std::cout << "Attempting to perform fit with no models!" << std::endl;
        exit(1);
    } else if (ReactorMod->IsInit()) {
        total_hist = (TH1D*)(ReactorMod->GetModelEsys()->Clone("Total_Model_Spectrum"));
    } else if (alphaNMod->IsInit()) {
        total_hist = (TH1D*)(alphaNMod->GetModelEsys()->Clone("Total_Model_Spectrum"));
    } else if (geoNuMod->IsInit()) {
        total_hist = (TH1D*)(geoNuMod->GetModelEsys()->Clone("Total_Model_Spectrum"));
    } else {
        total_hist = (TH1D*)(BasicMods->GetModelEsys(0)->Clone("Total_Model_Spectrum"));
    }
    total_hist->SetTitle("Total Model Spectrum");
    total_hist->Reset("ICES");

    if (ReactorMod->IsInit()) {
        ReactorMod->compute_spec();

        // Add total spectra
        hists.push_back((TH1D*)(ReactorMod->GetModelEsys()->Clone()));
        hists.at(hists.size()-1)->Reset("ICES");
        hists.at(hists.size()-1)->Add(ReactorMod->GetModelEsys());

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra (none here)
        // ReactorMod->Spectra(hists);  
    }

    if (alphaNMod->IsInit()) {
        alphaNMod->compute_spec();

        // Add total spectra
        hists.push_back((TH1D*)(alphaNMod->GetModelEsys()->Clone()));
        hists.at(hists.size()-1)->Reset("ICES");
        hists.at(hists.size()-1)->Add(alphaNMod->GetModelEsys());

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra
        alphaNMod->Spectra(hists);
    }

    if (geoNuMod->IsInit()) {
        geoNuMod->compute_spec();

        // Add total spectra
        hists.push_back((TH1D*)(geoNuMod->GetModelEsys()->Clone()));
        hists.at(hists.size()-1)->Reset("ICES");
        hists.at(hists.size()-1)->Add(geoNuMod->GetModelEsys());

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra
        geoNuMod->Spectra(hists);
    }

    for (unsigned int iMod = 0; iMod < BasicMods->GetNumMods(); ++iMod) {
        BasicMods->compute_spec(iMod);

        // Add total spectra
        hists.push_back((TH1D*)(BasicMods->GetModelEsys(iMod)->Clone()));
        hists.at(hists.size()-1)->Reset("ICES");
        hists.at(hists.size()-1)->Add(BasicMods->GetModelEsys(iMod));

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra
        BasicMods->Spectra(iMod, hists);
    }

    hists.push_back(total_hist);
}

TH1D* Fitter::DataHist() {return data;}
std::vector<double> Fitter::DataNtuple() {return data_E;}
void Fitter::DataNtuple(std::vector<double>& out_data) {
    for (unsigned int iDat = 0; iDat < data_E.size(); ++iDat) {
        out_data.push_back(data_E.at(iDat));
    }
}

void Fitter::InitTotFitModHist(TH1D* exampleHist) {
    tot_fitModel = (TH1D*)(exampleHist->Clone());
    tot_fitModel->Reset("ICES");
    tot_fitModel->SetName("tot_fit_model");
    tot_fitModel->SetTitle("Total Fit Spectrum");
    init_tot_fitModel = true;
    Emin = tot_fitModel->GetBinCenter(1);
    dE = tot_fitModel->GetBinCenter(2) - tot_fitModel->GetBinCenter(1);

    if (!setLims){
        iBinMin = 1;
        iBinMax = exampleHist->GetXaxis()->GetNbins();
    }
}

void Fitter::SetData(TH1D* Data) {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::SetData]: Setting TH1D Data." << std::endl;
    #endif
    data = Data;
    dataIsSet = true;
    binnedData = true;
}

void Fitter::SetData(TTree* Data) {
    #ifdef antinuDEBUG
        std::cout << "[Fitter::SetData]: Setting TTree Data." << std::endl;
    #endif
    data_E.clear();
    Double_t reconEnergy;
    Data->SetBranchAddress("energy", &reconEnergy);

    Int_t GTID;
    Data->SetBranchAddress("eventID", &GTID);
    std::vector<Int_t> cut_GTIDs = {7872258, 15357177, 8656140, 14585202, 13106488, 16172306, 15419458, 8780761, 8933258, 6599503, 7616837, 15373145, 13283754};
    bool is_cut;

    for (unsigned int a = 0; a < Data->GetEntries(); ++a) {
        Data->GetEntry(a);

        is_cut = false;
        for (unsigned int iGTID = 0; iGTID < cut_GTIDs.size(); ++iGTID) {
            if (GTID == cut_GTIDs.at(iGTID)) {
                is_cut = true;
                break;
            }
        }

        if (!is_cut) data_E.push_back(reconEnergy);
    }

    #ifdef antinuDEBUG
        std::cout << "[Fitter::SetData]: data_E = [" << data_E.at(0);
        for (unsigned int i = 1; i < data_E.size(); ++i) {
            std::cout << ", " << data_E.at(i);
        }
        std::cout << "]" << std::endl;
    #endif
    dataIsSet = true;
    binnedData = false;
}

void Fitter::SetBinLims(const unsigned int Bin_min, const unsigned int Bin_max) {
    iBinMin = Bin_min; iBinMax = Bin_max;
    setLims = true;
    #ifdef antinuDEBUG
        std::cout << "[Fitter::SetBinLims]: iBinMin = " << iBinMin << ", iBinMax = " << iBinMax << std::endl;
    #endif

    // Set models to use the same bin limits
    Reactor* ReactorMod = Reactor::GetInstance();
    alphaN* alphaNMod = alphaN::GetInstance();
    geoNu* geoNuMod = geoNu::GetInstance();
    Model* BasicMods = Model::GetInstance();
    
    if (ReactorMod->IsInit()) {
        ReactorMod->SetBinLims(iBinMin, iBinMax);
    } else {
        std::cout << "[Fitter::SetBinLims]: WARNING! Setting fitter bin limits before initialising Reactor model. Bin limits will not be set there too." << std::endl;
    }
    if (alphaNMod->IsInit()) {
        alphaNMod->SetBinLims(iBinMin, iBinMax);
    } else {
        std::cout << "[Fitter::SetBinLims]: WARNING! Setting fitter bin limits before initialising alphaN model. Bin limits will not be set there too." << std::endl;
    }
    if (geoNuMod->IsInit()) {
        geoNuMod->SetBinLims(iBinMin, iBinMax);
    } else {
        std::cout << "[Fitter::SetBinLims]: WARNING! Setting fitter bin limits before initialising geoNu model. Bin limits will not be set there too." << std::endl;
    }
    if (BasicMods->GetNumMods() > 0) {
        for (unsigned int iMod = 0; iMod < BasicMods->GetNumMods(); ++iMod) {
            BasicMods->SetBinLims(iMod, iBinMin, iBinMax);
        }
    } else {
        std::cout << "[Fitter::SetBinLims]: WARNING! Setting fitter bin limits before initialising basic models. Bin limits will not be set there too." << std::endl;
    }
}

void Fitter::GetParameter(Int_t ipar, char* name, Double_t& value, Double_t& verr, Double_t& vlow, Double_t& vhigh) {minuit->GetParameter(ipar, name, value, verr, vlow, vhigh);}
void Fitter::GetErrors(Int_t ipar, Double_t& eplus, Double_t& eminus, Double_t& eparab, Double_t& globcc) {minuit->GetErrors(ipar, eplus, eminus, eparab, globcc);}
double Fitter::GetCovarianceMatrixElement(int i, int j) {return minuit->GetCovarianceMatrixElement(i, j);}

void Fitter::resetVar(std::string Parname, double Value, double Verr, double Vlow, double Vhigh, bool holdConstant, bool isConstrained) {
    FitVars* Vars = FitVars::GetInstance();
    Vars->resetVar(Parname, Value, Verr, Vlow, Vhigh, holdConstant, isConstrained);

    minuit->SetParameter(Vars->findIdx(Parname), Parname.c_str(), Value, Verr, Vlow, Vhigh);
}

//end header guard
#endif


/* ~~~~~ Saving other potentially useful commands here ~~~~~ */

// unsigned int NtotParams = minuit->GetNumberTotalParameters();
// unsigned int NfreeParams = minuit->GetNumberFreeParameters();

// double chi2, edm, errdef;
// int nvpar, nparx;
// minuit->GetStats(chi2, edm, errdef, nvpar, nparx);

// void TFitter::GetConfidenceIntervals(Int_t n, Int_t ndim, const Double_t *x, Double_t *ci, Double_t cl)