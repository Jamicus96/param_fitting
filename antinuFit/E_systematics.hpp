// header guard:
#ifndef E_systematics
#define E_systematics

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "model.hpp"
#include "fitVar.hpp"


class Esys: public Model {
    private:
        // Scaling: E' = (1 + fDc) * E
        // Non-linearity: E' = E * (1 + fkB * E) / (1 + (fkB + fDkB) * E)
        // Smearing: Gaussian std = fSigPerRootE * sqrt(E)
        double fDc, fkB, fDkB, fSigPerRootE;

        // Middle of lowest energy bin, bin width, and the ratio fEmin/fDE
        double fEmin, fDE, fEratio;
        unsigned int iNumBins;

        TH1D* model_spec_scaled;

    public:
        // Constructors
        Esys(double dc, double kB, double dkB, double sigPerRootE);
        void operator = (const Esys& systematic);

        // Member functions
        void initialise();
        void apply_systematics();
        void apply_scaling();
        void apply_smearing();

        // Destructor
        ~Esys() {};
};

//end header guard
#endif