#include "E_systematics.hpp"



Esys::Esys(const Esys& systematic) {
    fKB = systematic.fKB;
    iC = systematic.iC;
    iKBp = systematic.iKBp;
    iSigPerRootE = systematic.iSigPerRootE;

    fEmin = systematic.fEmin;
    fDE = systematic.fDE;
    fEratio = systematic.fEratio;
    iNumBins = systematic.iNumBins;
    model_spec_scaled = systematic.model_spec_scaled;
    bIsInit = systematic.bIsInit;
}

void Esys::operator = (const Esys& systematic) {
    fKB = systematic.fKB;
    iC = systematic.iC;
    iKBp = systematic.iKBp;
    iSigPerRootE = systematic.iSigPerRootE;

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
    hists.push_back((TH1D*)(example_hist->Clone()));  // temporary hist for internal use
    tempHist_idx = hists.size() - 1;

    bIsInit = true;
}


/**
 * @brief Applies all three systematic corrections in a row: linear scaling -> non-linear scaling -> smearing
 * 
 */
void Esys::apply_systematics(TH1D* INhist, TH1D* OUThist) {
    if (!bIsInit) this->initialise(INhist);
    hists.at(tempHist_idx)->Reset("ICES");

    this->Esys::apply_scaling(INhist);
    this->Esys::apply_smearing(OUThist);
}

/**
 * @brief Applied the linear and non linear scaling (in that order)
 * 
 */
void Esys::apply_scaling(TH1D* INhist) {
    if (fKB == Vars.val(iKBp) && Vars.val(iC) == 1.0) {
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::apply_scaling]: No scaling." << std::endl;
        #endif
        // No scaling
        model_spec_scaled->Add(INhist);
    } else {
        unsigned int j_min;
        unsigned int j_max;
        double weight;
        for (unsigned int i = 1; i < iNumBins+1; ++i) {
            j_min = (unsigned int)(this->inv_scaling(fEmin + fDE * (i - 0.5)) / fDE + 0.5 - fEratio);
            j_max = (unsigned int)(this->inv_scaling(fEmin + fDE * (i + 0.5)) / fDE + 0.5 - fEratio);

            #ifdef SUPER_DEBUG
                std::cout << "[Esys::apply_scaling]: i = " << i << ", j € {" << j_min << ", " << j_max << "}." << std::endl;
            #endif

            if (j_min == j_max) {
                // Bin j maps over the whole of bin i (and possibly some of its neighbouring bins)
                // If scaling is a trivial transformation, it produces weight = 1
                weight = (this->inv_scaling(fEmin + fDE * (j_min + 0.5)) - this->inv_scaling(fEmin + fDE * (j_min - 0.5))) / fDE;
                hists.at(tempHist_idx)->AddBinContent(INhist->GetBinContent(j_min), weight);
            } else {
                weight = this->inv_scaling(fEmin + fDE * (j_min + 0.5)) / fDE - fEratio + (i - 0.5);
                hists.at(tempHist_idx)->AddBinContent(INhist->GetBinContent(j_min), weight);
                for (unsigned int j = j_min+1; j < j_max; ++j) {
                    hists.at(tempHist_idx)->AddBinContent(INhist->GetBinContent(j), 1.0);
                }
                weight = fEratio + (i + 0.5) - this->inv_scaling(fEmin + fDE * (j_min - 0.5)) / fDE;
                hists.at(tempHist_idx)->AddBinContent(INhist->GetBinContent(j_max), weight);
            }
        }
    }
}

/**
 * @brief Applies smearing
 * 
 */
void Esys::apply_smearing(TH1D* OUThist) {
    if (5.0 * Vars.val(iSigPerRootE) * std::sqrt(hists.at(tempHist_idx)->GetBinContent(iNumBins)) <= fDE) {
        // 5 sigmas in bin with highest energy is still smaller than bin width -> no smearing (includes sigma=0 case)
        OUThist->Add(hists.at(tempHist_idx));
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::apply_smearing]: No smearing." << std::endl;
        #endif
    } else {
        double sigma, weight;
        unsigned int num_bins_5sigmas;
        int j_min, j_max;
        for (unsigned int i = 1; i < iNumBins+1; ++i) {
            sigma = Vars.val(iSigPerRootE) * std::sqrt(hists.at(tempHist_idx)->GetBinContent(i));
            num_bins_5sigmas = (unsigned int)(5.0 * sigma / fDE);

            j_min = i - num_bins_5sigmas;
            j_max = i + num_bins_5sigmas;
            if (j_min < 1) j_min = 1;
            if (j_max > iNumBins) j_max = iNumBins;

            #ifdef SUPER_DEBUG
                std::cout << "[Esys::apply_smearing]: i = " << i << ", j € {" << j_min << ", " << j_max << "}." << std::endl;
            #endif

            for (unsigned int j = j_min; j < j_max+1; ++j) {
                weight = this->integ_normal(fDE * (j - i - 0.5) / sigma, fDE * (j - i + 0.5) / sigma);
                OUThist->AddBinContent(hists.at(tempHist_idx)->GetBinContent(j), weight);
            }
        }
    }
}


// double scaling(double E) {
//     double lin_scaled = (1.0 + vars.at(iDc)) * E;
//     return ((1.0 + fkB * lin_scaled) / (1.0 + (fkB + vars.at(vDkB)) * lin_scaled)) * lin_scaled;
// }

double Esys::inv_scaling(double E) {
    if (fKB == Vars.val(iKBp)) {
        // Only linear scaling
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::inv_scaling]: Only linear scaling. E = " << E << ", fC = " << fC << std::endl;
        #endif
        return E / fC;
    } else {
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::inv_scaling]: E = " << E << ", fC = " << fC << ", fKBp = " << fKBp << ", fKB = " << fKB << std::endl;
        #endif
        return (Vars.val(iKBp) * E - 1.0 + std::sqrt((1.0 - Vars.val(iKBp) * E) * (1.0 - Vars.val(iKBp) * E) + 4.0 * fKB * E)) / (2.0 * fKB * Vars.val(iC));
    }
}

double Esys::integ_normal(double x1, double x2) {
    return 0.5 * (std::erfc(x1 / std::sqrt(2.0)) - std::erfc(x2 / std::sqrt(2.0)));
}
