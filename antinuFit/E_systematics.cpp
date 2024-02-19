#include "E_systematics.hpp"


/**
 * @brief Construct a new Esys:: Esys object.
 * The basic idea of this class is to apply linear, non-linear and smearing energy corrections
 * to Model::model_spec, outputting the modified spectrum to Model::model_spec_sys.
 * 
 * @param dc  scaling: E' = (1 + fDc) * E
 * @param kB  non-linearity: E' = E * (1 + fkB * E) / (1 + (fkB + fDkB) * E)
 * @param dkB  Gaussian std = fSigPerRootE * sqrt(E)
 * @param sigPerRootE 
 */
Esys::Esys(double kB, FitVar* dc, FitVar* dkB, FitVar* sigPerRootE) {
    // Save variables
    vDc = dc;
    vDkB = dkB;
    vSigPerRootE = sigPerRootE;

    fDc = vDc->val();
    fDkB = vDkB->val();
    fSigPerRootE = vSigPerRootE->val();

    bIsInit = false;
}

void Esys::operator = (const Esys& systematic) {
    fkB = systematic.fkB;
    vDc = systematic.vDc;
    vDkB = systematic.vDkB;
    vSigPerRootE = systematic.vSigPerRootE;
    fDc = systematic.fDc;
    fDkB = systematic.fDkB;
    fSigPerRootE = systematic.fSigPerRootE;

    fEmin = systematic.fEmin;
    fDE = systematic.fDE;
    fEratio = systematic.fEratio;
    iNumBins = systematic.iNumBins;
    model_spec_scaled = systematic.model_spec_scaled;
    bIsInit = systematic.bIsInit;
}

/**
 * @brief Initialises histogram-dependent quantities.
 * Should be run in model daughter classes, after Model::model_spec has been set up.
 * Assumes histograms have constant bin sizes.
 * 
 */
void Esys::initialise(TH1D* example_hist) {
    // Record binning information
    fEmin = example_hist->GetBinCenter(1);
    fDE = example_hist->GetBinCenter(2) - example_hist->GetBinCenter(1);
    fEratio = fEmin / fDE;
    iNumBins = example_hist->GetXaxis()->GetNbins();

    // Set up empty histogram
    model_spec_scaled =  (TH1D*)(example_hist->Clone());  // temporary hist for internal use

    bIsInit = true;
}

/**
 * @brief Applies all three systematic corrections in a row: linear scaling -> non-linear scaling -> smearing
 * 
 */
void Esys::apply_systematics(Double_t* p, TH1D* INhist, TH1D* OUThist) {
    this->GetVarValues(p);
    model_spec_scaled->Reset("ICES");

    this->Esys::apply_scaling(INhist);
    this->Esys::apply_smearing(OUThist);
}

/**
 * @brief Applied the linear and non linear scaling (in that order)
 * 
 */
void Esys::apply_scaling(TH1D* INhist) {
    unsigned int j_min;
    unsigned int j_max;
    double weight;
    for (unsigned int i = 1; i < iNumBins+1; ++i) {
        j_min = (unsigned int)(this->inv_scaling(fEmin + fDE * (i - 0.5)) / fDE + 0.5 - fEratio);
        j_max = (unsigned int)(this->inv_scaling(fEmin + fDE * (i + 0.5)) / fDE + 0.5 - fEratio);

        if (j_min == j_max) {
            // Bin j maps over the whole of bin i (and possibly some of its neighbouring bins)
            // If scaling is a trivial transformation, it produces weight = 1
            weight = (this->inv_scaling(fEmin + fDE * (j_min + 0.5)) - this->inv_scaling(fEmin + fDE * (j_min - 0.5))) / fDE;
            model_spec_scaled->AddBinContent(INhist->GetBinContent(j_min), weight);
        } else {
            weight = this->inv_scaling(fEmin + fDE * (j_min + 0.5)) / fDE - fEratio + (i - 0.5);
            model_spec_scaled->AddBinContent(INhist->GetBinContent(j_min), weight);
            for (unsigned int j = j_min+1; j < j_max; ++j) {
                model_spec_scaled->AddBinContent(INhist->GetBinContent(j), 1.0);
            }
            weight = fEratio + (i + 0.5) - this->inv_scaling(fEmin + fDE * (j_min - 0.5)) / fDE;
            model_spec_scaled->AddBinContent(INhist->GetBinContent(j_max), weight);
        }
    }
}

/**
 * @brief Applies smearing
 * 
 */
void Esys::apply_smearing(TH1D* OUThist) {
    double sigma, weight;
    unsigned int num_bins_5sigmas;
    int j_min, j_max;
    for (unsigned int i = 1; i < iNumBins+1; ++i) {
        sigma = vars.at(iSigPerRootE) * std::sqrt(model_spec_scaled->GetBinContent(i));
        num_bins_5sigmas = (unsigned int)(sigma / fDE);

        j_min = i - num_bins_5sigmas;
        j_max = i + num_bins_5sigmas;
        if (j_min < 1) j_min = 1;
        if (j_max > iNumBins) j_max = iNumBins;
        for (unsigned int j = j_min; j < j_max+1; ++j) {
            weight = this->integ_normal(fDE * (j - i - 0.5) / sigma, fDE * (j - i + 0.5) / sigma);
            OUThist->AddBinContent(model_spec_scaled->GetBinContent(j), weight);
        }
    }
}


// double scaling(double E) {
//     double lin_scaled = (1.0 + vars.at(iDc)) * E;
//     return ((1.0 + fkB * lin_scaled) / (1.0 + (fkB + vars.at(vDkB)) * lin_scaled)) * lin_scaled;
// }

double Esys::inv_scaling(double E) {
    double kBp = fkB + vars.at(iDkB);
    return (kBp * E - 1.0 + std::sqrt((1.0 - kBp * E) * (1.0 - kBp * E) + 4.0 * fkB * E)) / (2.0 * fkB * (1.0 + vars.at(iDc)))
}

double Esys::integ_normal(double x1, double x2) {
    return 0.5 * (std::erfc(x1 / std::sqrt(2.0)) - std::erfc(x2 / std::sqrt(2.0)));
}
