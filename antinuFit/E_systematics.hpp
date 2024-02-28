// header guard:
#ifndef E_systematics
#define E_systematics

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "fitVars.hpp"


class Esys {
    private:
        static Esys *EsysInstance_;

        unsigned int numEsysts;
        std::vector<std::string> names;
        // Scaling: E' = (1 + fDc) * E
        // Non-linearity: E' = E * (1 + fkB * E) / (1 + (fkB + fDkB) * E)
        // Smearing: Gaussian std = fSigPerRootE * sqrt(E)
        std::vector<unsigned int> iC, iKBp, iSigPerRootE;
        std::vector<double> fKB;

        // Middle of lowest energy bin, bin width, and the ratio fEmin/fDE
        std::vector<double> fEmin, fDE, fEratio;
        std::vector<unsigned int> iNumBins;
        TH1D* tempHist;
        std::vector<bool> bIsInit;
    
    protected:
        Esys(): numEsysts(0) {}
        ~Esys() {}

    public:
        // Should not be cloneable.
        Esys(Esys &other) = delete;
        // Should not be assignable.
        void operator=(const Esys &) = delete;
        /**
         * This is the static method that controls the access to the Esys
         * instance. On the first run, it creates a Esys object and places it
         * into the static field. On subsequent runs, it returns the client existing
         * object stored in the static field.
         */
        static Esys *GetInstance();

        void AddEsys(const std::string name, const double kB, const unsigned int linScale_idx, const unsigned int kBp_idx, const unsigned int sigPerRootE_idx) {
            ++numEsysts;

            names.push_back(name);
            fKB.push_back(kB); iC.push_back(linScale_idx);
            iKBp.push_back(kBp_idx); iSigPerRootE.push_back(sigPerRootE_idx);

            fEmin.push_back(0.0); fDE.push_back(0.0);
            fEratio.push_back(0.0); iNumBins.push_back(0);
            bIsInit.push_back(false);
        }
        void AddEsys(const std::string name, const double kB, const std::string linScale_name, const std::string kBp_name, const std::string sigPerRootE_name) {
            FitVars* Vars = FitVars::GetInstance();
            AddEsys(name, kB, Vars->findIdx(linScale_name), Vars->findIdx(kBp_name), Vars->findIdx(sigPerRootE_name));
        }

        // Member functions
        unsigned int findIdx(const std::string EsysName) {
            for (unsigned int EsysIdx = 0; EsysIdx < numEsysts; ++EsysIdx) {
                if (EsysName == names.at(EsysIdx)) return EsysIdx;
            }
            std::cout << "[Esys::findIdx]: E-systematic '" << EsysName << "' not found!" << std::endl;
            exit(1);
            return 0;
        }

        void initialise(const unsigned int idx, TH1D* example_hist) {
            // Record binning information
            fEmin.at(idx) = example_hist->GetBinCenter(1);
            fDE.at(idx) = example_hist->GetBinCenter(2) - example_hist->GetBinCenter(1);
            fEratio.at(idx) = fEmin.at(idx) / fDE.at(idx);
            iNumBins.at(idx) = example_hist->GetXaxis()->GetNbins();

            // Set up empty histogram
            tempHist = (TH1D*)(example_hist->Clone());  // temporary hist for internal use
            bIsInit.at(idx) = true;
        }
        
        void initialise(std::string Parname, TH1D* example_hist) {
            FitVars* Vars = FitVars::GetInstance();
            initialise(Vars->findIdx(Parname), example_hist);
            }

        unsigned int GetNumEsysts() {return numEsysts;}

        void apply_systematics(const unsigned int idx, TH1D* INhist, TH1D* OUThist) {
            if (!bIsInit.at(idx)) initialise(idx, INhist);
            tempHist->Reset("ICES");

            Esys::apply_scaling(idx, INhist);
            Esys::apply_smearing(idx, OUThist);
        }
        
        void apply_scaling(const unsigned int idx, TH1D* INhist) {
            FitVars* Vars = FitVars::GetInstance();
            if (fKB.at(idx) == Vars->val(iKBp.at(idx)) && Vars->val(iC.at(idx)) == 1.0) {
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
                    j_min = (unsigned int)(inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (i - 0.5)) / fDE.at(idx) + 0.5 - fEratio.at(idx));
                    j_max = (unsigned int)(inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (i + 0.5)) / fDE.at(idx) + 0.5 - fEratio.at(idx));

                    #ifdef SUPER_DEBUG
                        std::cout << "[Esys::apply_scaling]: i = " << i << ", j € {" << j_min << ", " << j_max << "}." << std::endl;
                    #endif

                    if (j_min == j_max) {
                        // Bin j maps over the whole of bin i (and possibly some of its neighbouring bins)
                        // If scaling is a trivial transformation, it produces weight = 1
                        weight = (inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (j_min + 0.5)) - inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (j_min - 0.5))) / fDE.at(idx);
                        tempHist->AddBinContent(INhist->GetBinContent(j_min), weight);
                    } else {
                        weight = inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (j_min + 0.5)) / fDE.at(idx) - fEratio.at(idx) + (i - 0.5);
                        tempHist->AddBinContent(INhist->GetBinContent(j_min), weight);
                        for (unsigned int j = j_min+1; j < j_max; ++j) {
                            tempHist->AddBinContent(INhist->GetBinContent(j), 1.0);
                        }
                        weight = fEratio.at(idx) + (i + 0.5) - inv_scaling(idx, fEmin.at(idx) + fDE.at(idx) * (j_min - 0.5)) / fDE.at(idx);
                        tempHist->AddBinContent(INhist->GetBinContent(j_max), weight);
                    }
                }
            }
        }
        
        void apply_smearing(const unsigned int idx, TH1D* OUThist) {
            FitVars* Vars = FitVars::GetInstance();
            if (5.0 * Vars->val(iSigPerRootE.at(idx)) * std::sqrt(tempHist->GetBinContent(iNumBins.at(idx))) <= fDE.at(idx)) {
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
                    sigma = Vars->val(iSigPerRootE.at(idx)) * std::sqrt(tempHist->GetBinContent(i));
                    num_bins_5sigmas = (unsigned int)(5.0 * sigma / fDE.at(idx));

                    j_min = i - num_bins_5sigmas;
                    j_max = i + num_bins_5sigmas;
                    if (j_min < 1) j_min = 1;
                    if (j_max > iNumBins.at(idx)) j_max = iNumBins.at(idx);

                    #ifdef SUPER_DEBUG
                        std::cout << "[Esys::apply_smearing]: i = " << i << ", j € {" << j_min << ", " << j_max << "}." << std::endl;
                    #endif

                    for (unsigned int j = j_min; j < j_max+1; ++j) {
                        weight = integ_normal(idx, fDE.at(idx) * (j - i - 0.5) / sigma, fDE.at(idx) * (j - i + 0.5) / sigma);
                        OUThist->AddBinContent(tempHist->GetBinContent(j), weight);
                    }
                }
            }
        }

        void apply_systematics(std::string Parname, TH1D* INhist, TH1D* OUThist) {
            FitVars* Vars = FitVars::GetInstance();
            apply_systematics(Vars->findIdx(Parname), INhist, OUThist);
        }

        void apply_scaling(std::string Parname, TH1D* INhist) {
            FitVars* Vars = FitVars::GetInstance();
            apply_scaling(Vars->findIdx(Parname), INhist);
        }

        void apply_smearing(std::string Parname, TH1D* OUThist) {
            FitVars* Vars = FitVars::GetInstance();
            apply_smearing(Vars->findIdx(Parname), OUThist);
        }

        double inv_scaling(const unsigned int idx, const double E) {
            FitVars* Vars = FitVars::GetInstance();
            if (fKB.at(idx) == Vars->val(iKBp.at(idx))) {
                // Only linear scaling
                #ifdef SUPER_DEBUG
                    std::cout << "[Esys::inv_scaling]: Only linear scaling. E = " << E << ", fC = " << fC << std::endl;
                #endif
                return E / Vars->val(iC.at(idx));
            } else {
                #ifdef SUPER_DEBUG
                    std::cout << "[Esys::inv_scaling]: E = " << E << ", fC = " << fC << ", fKBp = " << fKBp << ", fKB = " << fKB << std::endl;
                #endif
                return (Vars->val(iKBp.at(idx)) * E - 1.0 + std::sqrt((1.0 - Vars->val(iKBp.at(idx)) * E) * (1.0 - Vars->val(iKBp.at(idx)) * E) + 4.0 * fKB.at(idx) * E)) / (2.0 * fKB.at(idx) * Vars->val(iC.at(idx)));
            }
        }

        double integ_normal(const unsigned int idx, const double x1, const double x2) {
            return 0.5 * (std::erfc(x1 / std::sqrt(2.0)) - std::erfc(x2 / std::sqrt(2.0)));
        }

        double inv_scaling(std::string Parname, const double E) {
            FitVars* Vars = FitVars::GetInstance();
            return inv_scaling(Vars->findIdx(Parname), E);
        }

        double integ_normal(std::string Parname, const double x1, const double x2) {
            FitVars* Vars = FitVars::GetInstance();
            return integ_normal(Vars->findIdx(Parname), x1, x2);
        }
};

// Static methods should be defined outside the class.
Esys* Esys::EsysInstance_{nullptr};

/**
 * The first time we call GetInstance we will lock the storage location
 *      and then we make sure again that the variable is null and then we
 *      set the value. RU:
 */
Esys *Esys::GetInstance() {
    if (EsysInstance_ == nullptr) {
        EsysInstance_ = new Esys();
    }
    return EsysInstance_;
}

//end header guard
#endif