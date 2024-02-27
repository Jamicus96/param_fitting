#include "E_systematics.hpp"



void Esys::AddEsys(const std::string name, const double kB, const unsigned int linScale_idx, const unsigned int kBp_idx, const unsigned int sigPerRootE_idx) {
    ++numEsysts;

    names.push_back(name);
    fKB.push_back(kB); iC.push_back(linScale_idx);
    iKBp.push_back(kBp_idx); iSigPerRootE.push_back(sigPerRootE_idx);

    fEmin.push_back(0.0); fDE.push_back(0.0);
    fEratio.push_back(0.0); iNumBins.push_back(0);
    tempHist_idx.push_back(0); bIsInit.push_back(false);
}

/**
 * @brief Initialises histogram-dependent quantities.
 * Should be run in model daughter classes, after Model::model_spec has been set up.
 * Assumes histograms have constant bin sizes.
 * 
 */
void Esys::initialise(const unsigned int idx, TH1D* example_hist) {
    // Record binning information
    fEmin.at(idx) = example_hist->GetBinCenter(1);
    fDE.at(idx) = example_hist->GetBinCenter(2) - example_hist->GetBinCenter(1);
    fEratio.at(idx) = fEmin.at(idx) / fDE.at(idx);
    iNumBins.at(idx) = example_hist->GetXaxis()->GetNbins();

    // Set up empty histogram
    tempHist = (TH1D*)(example_hist->Clone());  // temporary hist for internal use
    bIsInit.at(idx) = true;
}


/**
 * @brief Applies all three systematic corrections in a row: linear scaling -> non-linear scaling -> smearing
 * 
 */
void Esys::apply_systematics(const unsigned int idx, TH1D* INhist, TH1D* OUThist) {
    if (!bIsInit.at(idx)) this->initialise(idx, INhist);
    tempHist->Reset("ICES");

    this->Esys::apply_scaling(idx, INhist);
    this->Esys::apply_smearing(idx, OUThist);
}

/**
 * @brief Applied the linear and non linear scaling (in that order)
 * 
 */
void Esys::apply_scaling(const unsigned int idx, TH1D* INhist) {
    if (fKB.at(idx) == Vars.val(iKBp.at(idx)) && Vars.val(iC.at(idx)) == 1.0) {
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::apply_scaling]: No scaling." << std::endl;
        #endif
        // No scaling
        tempHist->Add(INhist);
    } else {
        unsigned int j_min;
        unsigned int j_max;
        double weight;
        for (unsigned int i = 1; i < iNumBins.at(idx)+1; ++i) {
            j_min = (unsigned int)(this->inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (i - 0.5)) / fDE.at(idx) + 0.5 - fEratio.at(idx));
            j_max = (unsigned int)(this->inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (i + 0.5)) / fDE.at(idx) + 0.5 - fEratio.at(idx));

            #ifdef SUPER_DEBUG
                std::cout << "[Esys::apply_scaling]: i = " << i << ", j € {" << j_min << ", " << j_max << "}." << std::endl;
            #endif

            if (j_min == j_max) {
                // Bin j maps over the whole of bin i (and possibly some of its neighbouring bins)
                // If scaling is a trivial transformation, it produces weight = 1
                weight = (this->inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (j_min + 0.5)) - this->inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (j_min - 0.5))) / fDE.at(idx);
                tempHist->AddBinContent(INhist->GetBinContent(j_min), weight);
            } else {
                weight = this->inv_scaling(fEmin.at(idx) + fDE.at(idx) * (j_min + 0.5)) / fDE.at(idx) - fEratio.at(idx) + (i - 0.5);
                tempHist->AddBinContent(INhist->GetBinContent(j_min), weight);
                for (unsigned int j = j_min+1; j < j_max; ++j) {
                    tempHist->AddBinContent(INhist->GetBinContent(j), 1.0);
                }
                weight = fEratio.at(idx) + (i + 0.5) - this->inv_scaling(fEmin.at(idx) + fDE.at(idx) * (j_min - 0.5)) / fDE.at(idx);
                tempHist->AddBinContent(INhist->GetBinContent(j_max), weight);
            }
        }
    }
}

/**
 * @brief Applies smearing
 * 
 */
void Esys::apply_smearing(const unsigned int idx, TH1D* OUThist) {
    if (5.0 * Vars.val(iSigPerRootE.at(idx)) * std::sqrt(tempHist->GetBinContent(iNumBins.at(idx))) <= fDE.at(idx)) {
        // 5 sigmas in bin with highest energy is still smaller than bin width -> no smearing (includes sigma=0 case)
        OUThist->Add(tempHist);
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::apply_smearing]: No smearing." << std::endl;
        #endif
    } else {
        double sigma, weight;
        unsigned int num_bins_5sigmas;
        int j_min, j_max;
        for (unsigned int i = 1; i < iNumBins.at(idx)+1; ++i) {
            sigma = Vars.val(iSigPerRootE.at(idx)) * std::sqrt(tempHist->GetBinContent(i));
            num_bins_5sigmas = (unsigned int)(5.0 * sigma / fDE.at(idx));

            j_min = i - num_bins_5sigmas;
            j_max = i + num_bins_5sigmas;
            if (j_min < 1) j_min = 1;
            if (j_max > iNumBins.at(idx)) j_max = iNumBins.at(idx);

            #ifdef SUPER_DEBUG
                std::cout << "[Esys::apply_smearing]: i = " << i << ", j € {" << j_min << ", " << j_max << "}." << std::endl;
            #endif

            for (unsigned int j = j_min; j < j_max+1; ++j) {
                weight = this->integ_normal(fDE.at(idx) * (j - i - 0.5) / sigma, fDE.at(idx) * (j - i + 0.5) / sigma);
                OUThist->AddBinContent(tempHist->GetBinContent(j), weight);
            }
        }
    }
}


// double scaling(double E) {
//     double lin_scaled = (1.0 + vars.at(iDc)) * E;
//     return ((1.0 + fkB * lin_scaled) / (1.0 + (fkB + vars.at(vDkB)) * lin_scaled)) * lin_scaled;
// }

double Esys::inv_scaling(const unsigned int idx, const double E) {
    if (fKB.at(idx) == Vars.val(iKBp.at(idx))) {
        // Only linear scaling
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::inv_scaling]: Only linear scaling. E = " << E << ", fC = " << fC << std::endl;
        #endif
        return E / Vars.val(iKBp.at(iC));
    } else {
        #ifdef SUPER_DEBUG
            std::cout << "[Esys::inv_scaling]: E = " << E << ", fC = " << fC << ", fKBp = " << fKBp << ", fKB = " << fKB << std::endl;
        #endif
        return (Vars.val(iKBp.at(idx)) * E - 1.0 + std::sqrt((1.0 - Vars.val(iKBp.at(idx)) * E) * (1.0 - Vars.val(iKBp.at(idx)) * E) + 4.0 * fKB.at(idx) * E)) / (2.0 * fKB.at(idx) * Vars.val(iC.at(idx)));
    }
}

double Esys::integ_normal(const unsigned int idx, const double x1, const double x2) {
    return 0.5 * (std::erfc(x1 / std::sqrt(2.0)) - std::erfc(x2 / std::sqrt(2.0)));
}
