#include "fitter.hpp"



Fitter::Fitter(TH1D* Data, std::vector<FitVar*>& Variables, std::vector<Model*>& Models,
               const unsigned int MaxNparams) : Fitter::Fitter(Data, MaxNparams) {

    variables = Variables;
    numVars = variables.size();
    var_bestFit_errs.resize(numVars);

    models = Models;
    numModels = models.size();

    for (unsigned int i = 0; i < numVars; ++i) {
        // Sets the parameter
        minuit->SetParameter(i, variables.at(i)->name(), variables.at(i)->prior(), variables.at(i)->err(), variables.at(i)->min(), variables.at(i)->max());
        // Binds the parameter's detail's addresses to fitVar's, so that changing the fitVar object changes this automatically
        minuit->GetParameter(i, variables.at(i)->name(), variables.at(i)->prior(), variables.at(i)->err(), variables.at(i)->min(), variables.at(i)->max())	
        variables.at(i)->ParIdx() = i;
    }
}

Fitter::Fitter(TH1D* Data, const unsigned int MaxNparams) : Fitter::Fitter(MaxNparams) {

    data = Data;
    numBins = data->GetXaxis()->GetNbins();
    tot_fitModel = (TH1D*)(data->Clone());
    tot_fitModel->SetName("tot_fit_model");
    tot_fitModel->SetTitle("Total Fit Spectrum");
}

Fitter::Fitter(const unsigned int MaxNparams) {

    maxNparams = MaxNparams;
    if (numVars > maxNparams) maxNparams = numVars;
    var_bestFits = new Double_t[maxNparams];

    //The default minimizer is Minuit, you can also try Minuit2
    TVirtualFitter::SetDefaultFitter("Minuit2");
    // TVirtualFitter::SetDefaultFitter("Minuit");
    minuit = TVirtualFitter::Fitter(0, maxNparams);
    
    arglist[0] = 0;
    // set print level
    minuit->ExecuteCommand("SET PRINT", arglist, 2);

    // minimiser settings
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance

    // (void (*fcn)(Int_t&, Double_t*, Double_t& f, Double_t*, Int_t))
    minuit->SetFCN(fitFunc);
}

Double_t Fitter::fitFunc(Double_t* p) {
    // Reset total model spectrum to zero
    std::fill(tot_fitModel.begin(), tot_fitModel.end(), 0.0);

    // Compute the total spectrum from all the models, using parameters p
    for (unsigned int iModel = 0; i < numModels; ++iModel) {
        models.at(iModel)->compute_spec(Double_t* p);
        tot_fitModel->Add(models.at(iModel)->Spectrum());
    }

    // Compute test statistic
}

void Fitter::fit_models() {
    // minimize
    minuit->ExecuteCommand("MIGRAD", arglist, 2);
}

double  Fitter::fit_models() {
    // minimize
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    //get result
    return this->fitFunc(var_bestFits);
}