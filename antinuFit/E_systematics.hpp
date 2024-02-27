// header guard:
#ifndef E_systematics
#define E_systematics

// include
#include <iostream>
#include <sstream>
#include <TH1.h>
#include "fitVars.hpp"
#include "fitter.hpp"


class Esys : Fitter {
    private:
        static std::vector<std::string> names;
        // Scaling: E' = (1 + fDc) * E
        // Non-linearity: E' = E * (1 + fkB * E) / (1 + (fkB + fDkB) * E)
        // Smearing: Gaussian std = fSigPerRootE * sqrt(E)
        static std::vector<unsigned int> iC, iKBp, iSigPerRootE;
        static std::vector<double> fKB;

        // Middle of lowest energy bin, bin width, and the ratio fEmin/fDE
        static std::vector<double> fEmin, fDE, fEratio;
        static std::vector<unsigned int> iNumBins, tempHist_idx;
        static std::vector<bool> bIsInit;

    public:
        // Constructors
        Esys() {numEsysts = 0;};
        static void AddEsys(const std::string name, const double kB, const unsigned int linScale_idx, const unsigned int kBp_idx, const unsigned int sigPerRootE_idx);
        static void AddEsys(const std::string name, const double kB, const std::string linScale_name, const std::string kBp_name, const std::string sigPerRootE_name) {
            AddEsys(name, kB, Vars.findIdx(linScale_name), Vars.findIdx(kBp_name), Vars.findIdx(sigPerRootE_name));
        };

        // Member functions
        static unsigned int findIdx(const std::string EsysName) {
            for (unsigned int EsysIdx = 0; EsysIdx < numEsysts; ++EsysIdx) {
                if (EsysName == names.at(EsysIdx)) return EsysIdx;
            }
            std::cout << "[Esys::findIdx]: E-systematic '" << EsysName << "' not found!" << std::endl;
            exit(1);
            return 0;
        };

        static void initialise(const unsigned int idx, TH1D* example_hist);
        static void initialise(std::string Parname, TH1D* example_hist) {initialise(Vars.findIdx(Parname), example_hist)};

        static unsigned int GetNumEsysts() {return numEsysts;};

        static void apply_systematics(const unsigned int idx, TH1D* INhist, TH1D* OUThist);
        static void apply_scaling(const unsigned int idx, TH1D* INhist);
        static void apply_smearing(const unsigned int idx, TH1D* OUThist);

        static void apply_systematics(std::string Parname, TH1D* INhist, TH1D* OUThist) {apply_systematics(Vars.findIdx(Parname), INhist, OUThist)};
        static void apply_scaling(std::string Parname, TH1D* INhist) {apply_scaling(Vars.findIdx(Parname), INhist)};
        static void apply_smearing(std::string Parname, TH1D* OUThist) {apply_smearing(Vars.findIdx(Parname), OUThist)};

        static double inv_scaling(const unsigned int idx, const double E);
        static double integ_normal(const unsigned int idx, const double x1, const double x2);

        static double inv_scaling(std::string Parname, const double E) {return inv_scaling(Vars.findIdx(Parname), E)};
        static double integ_normal(std::string Parname, const double x1, const double x2) {return integ_normal(Vars.findIdx(Parname), x1, x2)};

        // Destructor
        ~Esys() {};
};

//end header guard
#endif