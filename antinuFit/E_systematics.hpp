// header guard:
#ifndef E_systematics
#define E_systematics

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "fitVar.hpp"
#include "fitter.hpp"


class Esys : Fitter {
    private:
        // Scaling: E' = (1 + fDc) * E
        // Non-linearity: E' = E * (1 + fkB * E) / (1 + (fkB + fDkB) * E)
        // Smearing: Gaussian std = fSigPerRootE * sqrt(E)
        unsigned int iC, iKBp, iSigPerRootE;
        double fKB;
        unsigned int Esys_idx;

        // Middle of lowest energy bin, bin width, and the ratio fEmin/fDE
        double fEmin, fDE, fEratio;
        unsigned int iNumBins;

        unsigned int tempHist_idx;
        bool bIsInit;

    public:
        // Constructors
        Esys(const Esys& systematic);
        Esys(const double kB, const unsigned int linScale_idx, const unsigned int kBp_idx, const unsigned int sigPerRootE_idx) {fKB = kB; iC = linScale_idx; iKBp = kBp_idx; iSigPerRootE = sigPerRootE_idx;};
        Esys(unsigned int EsysIdx, const double kB, const unsigned int linScale_idx, const unsigned int kBp_idx, const unsigned int sigPerRootE_idx) : Esys(kB, linScale_idx, kBp_idx, sigPerRootE_idx) {Esys_idx = EsysIdx;};
        void operator = (const Esys& systematic);

        // Member functions
        void SetEsysIdx(unsigned int EsysIdx) {Esys_idx = EsysIdx;};
        void initialise(TH1D* example_hist);

        void apply_systematics(TH1D* INhist, TH1D* OUThist);
        void apply_scaling(TH1D* INhist);
        void apply_smearing(TH1D* OUThist);

        double inv_scaling(double E);
        double integ_normal(double x1, double x2);

        // Destructor
        ~Esys() {};
};

//end header guard
#endif