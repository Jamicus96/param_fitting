// header guard:
#ifndef E_systematics
#define E_systematics

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "fitVar.hpp"


class Esys {
    private:
        // Scaling: E' = (1 + fDc) * E
        // Non-linearity: E' = E * (1 + fkB * E) / (1 + (fkB + fDkB) * E)
        // Smearing: Gaussian std = fSigPerRootE * sqrt(E)
        FitVar *vC, *vKBp, *vSigPerRootE;
        double fC, fKBp, fSigPerRootE;
        double fKB;

        // Middle of lowest energy bin, bin width, and the ratio fEmin/fDE
        double fEmin, fDE, fEratio;
        unsigned int iNumBins;

        TH1D* model_spec_scaled;
        bool bIsInit;

    public:
        // Constructors
        Esys(double kB, FitVar* linScale, FitVar* kBp, FitVar* sigPerRootE);
        void operator = (const Esys& systematic);

        // Member functions
        void initialise(TH1D* example_hist);

        void GetVarValues(Double_t* p);

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