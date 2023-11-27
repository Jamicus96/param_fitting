#include "fitter.hpp"



Fitter::Fitter(TH1D* Data, std::vector<FitVar*>& Variables, std::vector<Model*>& Models) : Fitter::Fitter(Data) {

    numVars = Variables.size();
    variables.resize(numVars);
    var_bestFit_errs = new double[numVars];

    for (unsigned int i = 0; i < numVars; ++i) {
        variables.at(i) = Variables.at(i);
        // Sets the parameter
        minuit->SetParameter(i, variables.at(i)->name().c_str(), variables.at(i)->prior(), variables.at(i)->err(), variables.at(i)->min(), variables.at(i)->max());
        // Binds the parameter's detail's addresses to fitVar's, so that changing the fitVar object changes this automatically
        minuit->GetParameter(i, const_cast<char*>(variables.at(i)->name().c_str()), variables.at(i)->prior(), variables.at(i)->err(), variables.at(i)->min(), variables.at(i)->max());
        variables.at(i)->ParIdx() = i;
    }

    numModels = Models.size();
    models.resize(numModels);
    for (unsigned int i = 0; i < numModels; ++i) {
        models.at(i) = Models.at(i);
    }
}

Fitter::Fitter(TH1D* Data) : Fitter::Fitter() {

    data = Data;
    numBins = data->GetXaxis()->GetNbins();
    tot_fitModel = (TH1D*)(data->Clone());
    tot_fitModel->SetName("tot_fit_model");
    tot_fitModel->SetTitle("Total Fit Spectrum");
}

Fitter::Fitter() {

    var_bestFits = new Double_t[numVars];

    //The default minimizer is Minuit, you can also try Minuit2
    TVirtualFitter::SetDefaultFitter("Minuit2");
    // TVirtualFitter::SetDefaultFitter("Minuit");
    minuit = TVirtualFitter::Fitter(0, numVars);
    
    arglist[0] = 0;
    // set print level
    minuit->ExecuteCommand("SET PRINT", arglist, 2);

    // minimiser settings
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance

    // (void (*fcn)(Int_t&, Double_t*, Double_t& f, Double_t*, Int_t))
    SetFunc();
}

Double_t Fitter::fitFunc(Double_t* p) {
    // Reset total model spectrum to zero
    tot_fitModel->Reset("ICES");

    // Compute the total spectrum from all the models, using parameters p
    for (unsigned int iModel = 0; iModel < numModels; ++iModel) {
        models.at(iModel)->compute_spec(p);
        tot_fitModel->Add(models.at(iModel)->Spectrum());
    }

    // Compute test statistic
    return ExtendedConstrainedLogLikelihood(p);
}

double Fitter::ExtendedConstrainedLogLikelihood(Double_t* p) {
    // Model PDFs have already been scaled by their respective norms, and added together
    double logL = - tot_fitModel->Integral();
    for (unsigned int ibin = 1; ibin < numBins+1; ++ibin) {
        logL += data->GetBinContent(ibin) * log(tot_fitModel->GetBinContent(ibin));
    }
    // Add in constraints
    for (unsigned int iVar = 0; iVar < numVars; ++ iVar) {
        // If variable is being held constant, assume it's equal to its prior most likely value
        if (!(variables.at(iVar)->isConstant())) {
            // Add gaussian constraint
            logL += 0.5 * (p[iVar] - variables.at(iVar)->prior())*(p[iVar] - variables.at(iVar)->prior()) / (variables.at(iVar)->err()*variables.at(iVar)->err());
        }
    }
    return logL;
}

double Fitter::fit_models() {
    // minimize
    minuit->ExecuteCommand("MIGRAD", arglist, 2);

    //get result
    return fitFunc(var_bestFits);
}