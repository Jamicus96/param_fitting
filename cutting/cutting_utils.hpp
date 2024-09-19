#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>

#define VERBOSE
#define USING_RUN_NUM
// #define USING_PDF_PADDING  // use broader prompt E cuts for IBDs, to make PDFs suitable for energy systematics

ULong64_t dcAnalysisWord = 0x2100000042C2;  // Converts hex to decimal

// Global DC and FV cuts + AV offset
double highNhit_deltaT = 20;  // [s]
double owlNhit_deltaT = 10;  // [us]
double AV_offset = 186;  // [mm] average value from https://snopl.us/monitoring/scint_level
double FV_CUT = 5700;  // Max radius [mm]

// Define BiPo cut values (taken from Ziping docdb-7812 slide 5, except the Delta_t > 2us cut since don't care about BiPo212 contamination)
double BIPO_MIN_PROMPT_E = 0.55, BIPO_MAX_PROMPT_E = 3.5;  // [MeV]
double BIPO_MIN_DELAYED_E = 0.6, BIPO_MAX_DELAYED_E = 1.1;  // [MeV]
double BIPO_MIN_DELAY = 0, BIPO_MAX_DELAY = 1E6;  // Delta t cut [ns]
double BIPO_MAX_DIST = 800;  // Delta R cut [mm]

// IBD cuts
double IBD_MIN_PROMPT_E = 0.9, IBD_MAX_PROMPT_E = 8.0;  // [MeV]
double IBD_PDF_MIN_PROMPT_E = 0.5, IBD_PDF_MAX_PROMPT_E = 9.0;  // [MeV], extra padding for PDFs
double IBD_MIN_DELAYED_E = 1.85, IBD_MAX_DELAYED_E = 2.4;  // [MeV]
double IBD_MIN_DELAY = 400, IBD_MAX_DELAY = 0.8E6;  // Delta t cut [ns]
double IBD_MAX_DIST = 1500;  // Delta R cut [mm]

// alpha-n process split energies
double CARBON12_SCATTER_E_MAX = 5.4;  // [MeV]  Could just simulate different process separately?
double PROTON_RECOIL_E_MAX = 3.5;  // [MeV]  Ideally the same as in cut_data.cpp


/* ~~~~~~~~~~~~ General functions (E correction, DC, classifier) ~~~~~~~~~~~~ */

double EnergyCorrection(const double E, TVector3 pos, const bool is_data, RAT::DU::DetectorStateCorrection& stateCorr, RAT::DU::ReconCalibrator* e_cal) {
    // Data vs MC energy correction (Tony's)
    double Ecorr = e_cal->CalibrateEnergyRTF(is_data, E, std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()), pos.Z()); // gives the new E

    // #ifdef USING_RUN_NUM
    //     // Correct for position coverage dependence (Logan's)
    //     RAT::DU::Point3D position(0, pos);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
    //     Ecorr /= stateCorr.GetCorrectionPos(position, 0, 0) / stateCorr.GetCorrection(9394, 0.75058); // a correction factor (divide E by it)
    // #endif

    return Ecorr;
}

bool dcAppliedAndPassed(const bool is_data, ULong64_t dcApplied, ULong64_t dcFlagged) {
    if (!is_data) return true;
    if (!(((dcApplied & dcAnalysisWord) & dcFlagged ) == (dcApplied & dcAnalysisWord))) return false;  // failed DC cut
    return true;
}

bool pass_classifier(const double energy, const double class_result, const double class_cut) {
    // Only check classifier result for events in PDF1 (proton recoil, E < 3MeV)
    if (energy > PROTON_RECOIL_E_MAX) return true;
    if (class_result > class_cut) return true;

    return false;
}

/* ~~~~~~~~~~~~ IBD cuts ~~~~~~~~~~~~ */

bool pass_prompt_cuts_IBD(const double energy, const TVector3& position) {
    #ifdef USING_PDF_PADDING
        if (energy < IBD_PDF_MIN_PROMPT_E) return false;  // min energy cut (MeV)
        if (energy > IBD_PDF_MAX_PROMPT_E) return false;  // max energy cut (MeV)
    #else
        if (energy < IBD_MIN_PROMPT_E) return false;  // min energy cut (MeV)
        if (energy > IBD_MAX_PROMPT_E) return false;  // max energy cut (MeV)
    #endif
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts_IBD(const double energy, const TVector3& position) {
    if (energy < IBD_MIN_DELAYED_E) return false;  // min energy cut (MeV)
    if (energy > IBD_MAX_DELAYED_E) return false;  // max energy cut (MeV)
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts_IBD(const double delay, const TVector3& prompt_pos, const TVector3& delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < IBD_MIN_DELAY) return false;  // min delay cut (ns)
    if (delay > IBD_MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > IBD_MAX_DIST) return false;  // max distance cut (mm)

    return true;
}

/* ~~~~~~~~~~~~ BiPo cuts ~~~~~~~~~~~~ */

bool pass_prompt_cuts_BiPo(const double energy, const TVector3& position) {
    if (energy < BIPO_MIN_PROMPT_E) return false;  // min energy cut (MeV)
    if (energy > BIPO_MAX_PROMPT_E) return false;  // max energy cut (MeV)
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts_BiPo(const double energy, const TVector3& position) {
    if (energy < BIPO_MIN_DELAYED_E) return false;  // min energy cut (MeV)
    if (energy > BIPO_MAX_DELAYED_E) return false;  // max energy cut (MeV)
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts_BiPo(const double delay, const TVector3& prompt_pos, const TVector3& delayed_pos) {
    // double delay = (delayed_time - prompt_time) / 50E6 * 1E9; // convert number of ticks in 50MHz clock to ns
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < BIPO_MIN_DELAY) return false;  // min delay cut (ns)
    if (delay > BIPO_MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > BIPO_MAX_DIST) return false;  // max distance cut (mm)

    return true;
}