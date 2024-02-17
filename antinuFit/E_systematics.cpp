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
Esys::Esys(double dc, double kB, double dkB, double sigPerRootE) {
    // Assign variables
    fDc = dc;
    fkB = kB;
    fDkB = dkB;
    fSigPerRootE = sigPerRootE;
}

void Esys::operator = (const Esys& systematic) {
    fDc = systematic.fDc;
    fkB = systematic.fkB;
    fDkB = systematic.fDkB;
    fSigPerRootE = systematic.fSigPerRootE;
}

/**
 * @brief Initialises histogram-dependent quantities.
 * Should be run in model daughter classes, after Model::model_spec has been set up.
 * Assumes histograms have constant bin sizes.
 * 
 */
void Esys::initialise() {
    // Record binning information
    fEmin = model_spec->GetBinCenter(1);
    fDE = model_spec->GetBinCenter(2) - model_spec->GetBinCenter(1);
    fEratio = fEmin / fDE;
    iNumBins = model_spec->GetXaxis()->GetNbins();

    // Set up empty histograms
    model_spec_scaled =  (TH1D*)(model_spec->Clone());  // temporary hist for internal use
    model_spec_sys = (TH1D*)(model_spec->Clone());  // output histogram, used by Model
}

/**
 * @brief Applies all three systematic corrections in a row: linear scaling -> non-linear scaling -> smearing
 * 
 */
void Esys::apply_systematics() {
    Esys::apply_scaling();
    Esys::apply_smearing();
}

double scaling(double E) {
    double lin_scaled = (1.0 + fDc) * E;
    return ((1.0 + fkB * lin_scaled) / (1.0 + (fkB + fDkB) * lin_scaled)) * lin_scaled;
}

double inv_scaling(double E) {
    double kBp = fkB + fDkB;
    return (kBp * E - 1.0 + std::sqrt((1.0 - kBp * E) * (1.0 - kBp * E) + 4.0 * fkB * E)) / (2.0 * fkB * (1.0 + fDc))
}

/**
 * @brief Applied the linear and non linear scaling (in that order)
 * 
 */
void Esys::apply_scaling() {
    model_spec_scaled->Reset("ICES");
    unsigned int j_min;
    unsigned int j_max;
    double weight;
    for (unsigned int i = 1; i < iNumBins+1; ++i) {
        j_min = (unsigned int)(inv_scaling(fEmin + fDE * (i - 0.5)) / fDE + 0.5 - fEratio);
        j_max = (unsigned int)(inv_scaling(fEmin + fDE * (i + 0.5)) / fDE + 0.5 - fEratio);

        if (j_min == j_max) {
            // Bin j maps over the whole of bin i (and possibly some of its neighbouring bins)
            // If scaling is a trivial transformation, it produces weight = 1
            weight = (inv_scaling(fEmin + fDE * (j_min + 0.5)) - inv_scaling(fEmin + fDE * (j_min - 0.5))) / fDE;
            model_spec_scaled->AddBinContent(model_spec->GetBinContent(j_min), weight);
        } else {
            weight = inv_scaling(fEmin + fDE * (j_min + 0.5)) / fDE - fEratio + (i - 0.5);
            model_spec_scaled->AddBinContent(model_spec->GetBinContent(j_min), weight);
            for (unsigned int j = j_min+1; j < j_max; ++j) {
                model_spec_scaled->AddBinContent(model_spec->GetBinContent(j), 1.0);
            }
            weight = fEratio + (i + 0.5) - inv_scaling(fEmin + fDE * (j_min - 0.5)) / fDE;
            model_spec_scaled->AddBinContent(model_spec->GetBinContent(j_max), weight);
        }
    }
}

double integ_normal(double x1, double x2) {
    return 0.5 * (std::erfc(x1 / std::sqrt(2.0)) - std::erfc(x2 / std::sqrt(2.0)));
}

/**
 * @brief Applies smearing
 * 
 */
void Esys::apply_smearing() {
    model_spec_sys->Reset("ICES");

    double sigma, weight;
    unsigned int num_bins_5sigmas;
    int j_min, j_max;
    for (unsigned int i = 1; i < iNumBins+1; ++i) {
        sigma = fSigPerRootE * std::sqrt(model_spec_scaled->GetBinContent(i));
        num_bins_5sigmas = (unsigned int)(sigma / fDE);

        j_min = i - num_bins_5sigmas;
        j_max = i + num_bins_5sigmas;
        if (j_min < 1) j_min = 1;
        if (j_max > iNumBins) j_max = iNumBins;
        for (unsigned int j = j_min; j < j_max+1; ++j) {
            weight = integ_normal(fDE * (j - i - 0.5) / sigma, fDE * (j - i + 0.5) / sigma);
            model_spec_sys->AddBinContent(model_spec_scaled->GetBinContent(j), weight);
        }
    }
}
