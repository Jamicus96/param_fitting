////////////////////////////////////////////////////////////////////
/// \file
///
/// \brief Re-write...
///
/// \author James Page <j.page@sussex.ac.uk>
///
/// REVISION HISTORY:\n
///
/// \details EV Calibrated hit times are plotted minus transit times
/// based on the MC position or the fitted position.
/// Multiple PDFs can be made at the same time, and cuts can be performed
/// for all events, or specific to each PDF.
///
/// To compile: g++ -g -std=c++1y getProjVec.cpp -o getProjVec.exe `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
///
////////////////////////////////////////////////////////////////////

#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include "RAT/DU/Point3D.hh"
#include <RAT/DS/FitResult.hh>

#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include <RAT/DB.hh>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TVectorD.h>

#include <string>
#include <fstream>
#include "cutting_utils.hpp"

#define LPC_uses_Point3D  // Backwards compatibility: Comment out if RAT version < 7.0.11 (Light Path Calculation using TVector3 vs Point3D)


void Apply_tagging_and_cuts(std::string input_root_address, double classifier_cut, int runNum, int GTID, bool is_data, bool verbose, std::string fitName = "");

/* ~~~~~~~~~~~~~~~~~~~~~~ MAIN FUNCTION ~~~~~~~~~~~~~~~~~~~~~ */

int main(int argc, char** argv) {
    std::string input_root_address = argv[1];
    double classifier_cut = std::stod(argv[2]);
    int runNum = std::stoi(argv[3]);
    int GTID = std::stoi(argv[4]);  // only look for specific GTID
    bool is_data = std::stoi(argv[5]);
    bool verbose = std::stoi(argv[6]);

    // Loop through files to get info from every event (including t_res), and print all to text file
    if (verbose) {std::cout << "Getting info..." << std::endl;}
    Apply_tagging_and_cuts(input_root_address, classifier_cut, runNum, GTID, is_data, verbose);

    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~ CUT FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~ */


bool pass_prompt_cuts(double energy, TVector3 position) {
    if (energy < IBD_MIN_PROMPT_E) return false;  // min energy cut (MeV)
    if (energy > IBD_MAX_PROMPT_E) return false;  // max energy cut (MeV)
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_delayed_cuts(double energy, TVector3 position) {
    if (energy < IBD_MIN_DELAYED_E) return false;  // min energy cut (MeV)
    if (energy > IBD_MAX_DELAYED_E) return false;  // max energy cut (MeV)
    if (sqrt(position.X()*position.X() + position.Y()*position.Y() + (position.Z() - AV_offset)*(position.Z() - AV_offset)) > FV_CUT) return false;  // FV cut (mm)

    return true;
}

bool pass_coincidence_cuts(double delay, TVector3 prompt_pos, TVector3 delayed_pos) {
    double distance = (delayed_pos - prompt_pos).Mag();

    if (delay < IBD_MIN_DELAY) return false;  // min delay cut (ns)
    if (delay > IBD_MAX_DELAY) return false;  // max delay cut (ns)
    if (distance > IBD_MAX_DIST) return false;  // max distance cut (mm)

    return true;
}



/* ~~~~~~~~~~~~~~~~~~~~~~ PRIMARY FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~ */

void Apply_tagging_and_cuts(std::string input_root_address, double classifier_cut, int runNum, int GTID, bool is_data, bool verbose, std::string fitName) {
    if (verbose) {std::cout << "Running Apply_tagging_and_cuts()" << std::endl;}

    std::cout << "Reading in file: " << input_root_address << std::endl;
    RAT::DU::DSReader dsReader(input_root_address);

    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();
    RAT::DU::DetectorStateCorrection stateCorr = RAT::DU::Utility::Get()->GetDetectorStateCorrection();
    RAT::DU::ReconCalibrator e_cal = RAT::DU::Utility::Get()->GetReconCalibrator();

    /*********** Prepare classifier info ***********/

    // Set up db access
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    RAT::DS::Run run;

    #ifndef USING_RUN_NUM
        db->SetAirplaneModeStatus(true);
        db->LoadDefaults();
    #endif
    #ifdef USING_RUN_NUM
        db->LoadDefaults();
        db->SetServer("postgres://snoplus@pgsql.snopl.us:5400/ratdb");
        std::cout << "RAT DB tag: " << db->GetDBTag() << std::endl;

        if (runNum > 0) {
            run.SetRunID(runNum);
            db->BeginOfRun(run);
        }
    #endif


    RAT::DBLinkPtr linkdb = db->GetLink("ALPHAN_IBD_PROMPT_FISHER", "labppo_2p2_scintillator");

    std::vector<double> fTimes = linkdb->GetDArray("times");
    std::vector<double> fA_vec = linkdb->GetDArray("a_vec");

    /*********** Loop through all files, entries, events, and PMTs ***********/
    
    
    std::cout << "Looping through entries..." << std::endl;
    // loops through entries
    bool GTID_found = false;
    for (unsigned int iEntry = 0; iEntry < dsReader.GetEntryCount(); ++iEntry) {
        if (iEntry % 10000 == 0) std::cout << "Done " << ((float)iEntry / (float)dsReader.GetEntryCount()) * 100.0 << "%" << std::endl;
        RAT::DS::Entry rDS = dsReader.GetEntry(iEntry);

        // loops through events
        for (unsigned int iEV = 0; iEV < rDS.GetEVCount(); ++iEV) {
            if (verbose) std::cout << "iEV = " << iEV << std::endl;
            RAT::DS::EV rEV = rDS.GetEV(iEV);

            // only look out for specific event
            if (rEV.GetGTID() == GTID) {
                GTID_found = true;
                std::cout << "Found GTID " << GTID << "!" << std::endl;

                if (fitName == "") fitName = rEV.GetDefaultFitName();
                try {
                    // Get recon info
                    RAT::DS::FitResult fitResult = rEV.GetFitResult(fitName);
                    if (!fitResult.GetValid()) std::cout << "Event " << GTID << " did not reconstruct!" << std::endl;; // fit invalid
                    RAT::DS::FitVertex& rVertex = fitResult.GetVertex(0);
                    if (!(rVertex.ValidPosition() && rVertex.ValidTime() && rVertex.ValidEnergy())) std::cout << "Event " << GTID << " did not reconstruct!" << std::endl; // fit invalid

                    double Ecorr = rVertex.GetEnergy();
                    double vertex_time = rVertex.GetTime();
                    TVector3 pos = rVertex.GetPosition();

                    // Data vs MC energy correction (Tony's)
                    Ecorr = e_cal.CalibrateEnergyRTF(is_data, Ecorr, std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()), pos.Z()); // gives the new E

                    #ifdef LPC_uses_Point3D
                        // Correct for position coverage dependence (Logan's)
                        RAT::DU::Point3D position(0, pos);  // position of event [mm] (as Point3D in PSUP coordinates, see system_id in POINT3D_SHIFTS tables)
                        // Ecorr /= stateCorr.GetCorrectionPos(position, 0, 0) / stateCorr.GetCorrection(9394, 0.75058); // a correction factor (divide E by it)
                    #else
                        TVector3 position = pos;
                    #endif

                    std::cout << "Got reconstructed info!" << std::endl;

                    // Adapted from RAT:
                    double Result = 0.0;
                    double fTimeStep = fTimes[1] - fTimes[0];
                    unsigned int index;
                    double corTime;
                    RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
                    for (size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); ++iPMT) {
                        RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT(iPMT);
                        // Use new time residual calculator
                        corTime = fTRCalc.CalcTimeResidual(pmtCal, position, vertex_time);
                        if (corTime > fTimes.back() || corTime < fTimes[0]) continue; // ignore t_res outside specified range
                        index = (corTime - fTimes[0]) / fTimeStep;
                        Result += fA_vec[index];
                    }

                    // Compute R component of Fisher discriminant statistic
                    Result += pos.Mag() * fA_vec[fA_vec.size() - 1];

                    std::cout << "Event " << GTID << " classifier result: " << Result << std::endl;
                    std::cout << "Event " << GTID << " classifier pass: " << (Result > -8.81) << std::endl;
                }
                catch (RAT::DS::DataNotFound&) {std::cout << "Event " << GTID << " did not reconstruct!" << std::endl;;}  // no fit data
                catch (RAT::DS::FitCollection::NoResultError&) {std::cout << "Event " << GTID << " did not reconstruct!" << std::endl;;} // no fit result by the name of fitName
                catch (RAT::DS::FitResult::NoVertexError&) {std::cout << "Event " << GTID << " did not reconstruct!" << std::endl;;} // no fit vertex
                catch (RAT::DS::FitVertex::NoValueError&) {std::cout << "Event " << GTID << " did not reconstruct!" << std::endl;;} // position or time missing
            }
        }
    }
    dsReader.Delete();

    std::cout << "GTID_found: " << GTID_found << std::endl;
}
