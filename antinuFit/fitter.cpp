#include "fitter.hpp"



TVirtualFitter* Fitter::minuit;
double Fitter::arglist[2];


TH1D* Fitter::data, *Fitter::tot_fitModel;
unsigned int Fitter::numBins;

double Fitter::minfuncOut, Fitter::edm, Fitter::errdef;
int Fitter::nvpar, Fitter::nparx;

FitVars Fitter::Vars;
Esys Fitter::Esysts;
std::vector<Model*> Fitter::Mods;
std::vector<TH1D*> Fitter::hists;

unsigned int Fitter::numMods;
unsigned int Fitter::numHists;


Fitter::Fitter(TH1D& Data) {

    data = Data;
    numBins = data.GetXaxis()->GetNbins();
    tot_fitModel = (TH1D*)(data.Clone());
    tot_fitModel->SetName("tot_fit_model");
    tot_fitModel->SetTitle("Total Fit Spectrum");

    //The default minimizer is Minuit, you can also try Minuit2
    // TVirtualFitter::SetDefaultFitter("Minuit2");
    TVirtualFitter::SetDefaultFitter("Minuit");
    minuit = TVirtualFitter::Fitter(0, Vars.GetNumVars());
    
    arglist[0] = 0;
    // set print level
    minuit->ExecuteCommand("SET PRINT", arglist, 2);

    // minimiser settings
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance

    // (void (*fcn)(Int_t&, Double_t*, Double_t& f, Double_t*, Int_t))
    SetFunc();
}

Fitter::Fitter() {}

Double_t Fitter::fitFunc(Double_t* p) {
    // Reset total model spectrum to zero
    tot_fitModel->Reset("ICES");

    // Compute the total spectrum from all the models, using parameters p
    Vars.GetVarValues(p);
    for (unsigned int iModel = 0; iModel < numMods; ++iModel) {
        Mods.at(iModel)->compute_spec();
        tot_fitModel->Add(Mods.at(iModel)->GetModelEsys());
        #ifdef SUPER_DEBUG
            std::cout << "[Fitter::fitFunc]: Added model " << iModel << ", tot_fitModel integral = " << tot_fitModel->Integral() << std::endl;
        #endif
    }

    // Compute test statistic (minus, since minuit minizes, but we want max)
    return - ExtendedConstrainedLogLikelihood();
}

double Fitter::ExtendedConstrainedLogLikelihood() {
    // Model PDFs have already been scaled by their respective norms, and added together
    // Assume data and models have the same binning (could generalise at some point)
    double logL = - tot_fitModel->Integral();
    for (unsigned int ibin = 1; ibin < numBins+1; ++ibin) {
        if (tot_fitModel->GetBinContent(ibin) > 0) logL += data->GetBinContent(ibin) * log(tot_fitModel->GetBinContent(ibin));
    }
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: +PDFs logL = " << logL << std::endl;
    #endif

    // Add in constraints
    for (unsigned int iVar = 0; iVar < numVars; ++ iVar) {
        // If variable is being held constant, assume it's equal to its prior most likely value
        if (!Vars.isConstant(iVar)) {
            // Add gaussian constraint
            logL -= 0.5 * (Vars.val(iVar) - Vars.val.prior(iVar))*(Vars.val(iVar) - Vars.val.prior(iVar)) / (Vars.err(iVar)*Vars.err(iVar));
        }
    }
    #ifdef SUPER_DEBUG
        std::cout << "[Fitter::ExtendedConstrainedLogLikelihood]: +constraints logL = " << logL << std::endl;
    #endif

    return logL;
}

double Fitter::fit_models() {
    // minimize
    #ifdef DEBUG
        std::cout << "[Fitter::fit_models]: Running MIGRAD!" << std::endl;
        std::cout << "[Fitter::fit_models]: s_12^2 = " << variables.at(variables.size() - 2)->val() << std::endl;
        std::cout << "[Fitter::fit_models]: Dm21^2 = " << variables.at(variables.size() - 4)->val() << std::endl;
    #endif

    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    //get result
    minuit->GetStats(minfuncOut, edm, errdef, nvpar, nparx);
    std::cout << "[Fitter::fit_models]: minfuncOut = " << minfuncOut << std::endl;
    return minfuncOut;
}

void Fitter::GetAllSpectra(std::vector<TH1D*>& hists) {
    TH1D* total_hist = (TH1D*)(Mods.at(0)->Spectrum()->Clone("Total_Model_Spectrum"));
    total_hist->SetTitle("Total Model Spectrum");
    for (unsigned int iModel = 0; iModel < numMods; ++iModel) {
        // Add total spectra
        hists.push_back((TH1D*)(Mods.at(iModel)->Spectrum()->Clone()));
        hists.at(hists.size()-1)->Add(Mods.at(iModel)->Spectrum());

        // Add to total overall specrtrum
        total_hist->Add(hists.at(hists.size()-1));

        // Add constituent spectra
        Mods.at(iModel)->Spectra(hists);
    }
    hists.push_back(total_hist);
}
