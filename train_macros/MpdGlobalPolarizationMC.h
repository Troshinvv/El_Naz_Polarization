#ifndef MPDGLOBALPOLARIZATIONMC_H
#define MPDGLOBALPOLARIZATIONMC_H

#include <deque>

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdHelix.h"
#include "FairMCEventHeader.h"
#include "MpdAnalysisTask2.h"

class MpdGlobalPolarizationMC : public MpdAnalysisTask2 {
   ClassDef(MpdGlobalPolarizationMC, 1);

public:
   MpdGlobalPolarizationMC() {}
   MpdGlobalPolarizationMC(const char *name, const char *outputName = "taskName");
   ~MpdGlobalPolarizationMC() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

private:
   
   int NITER_CENT;       // number of centrality bins (defined for 4, 7, 10)
   int NITER;            // number of angular bins
   int cent_cut_choice;  // choice of centrality cut (0 - no cent cut, 1 - with cent cut)
   double cent_cut;      // value for centrality cut
   int particle_choice;  // particle choice (pdg of hyperon) --- currently Lambda (3122) or anti-Lambda (-3122))
  
   // Event properties
   std::string mOutFile = "histos.root";

   TClonesArray         *mMCTracks        = nullptr;

   // Histograms
   TList mHistoList;

   // General QA
   TH1F *mhCentrality         = nullptr;
   TH1D *NCentr               = nullptr;
   TH1D *Resolution_EP1_true  = nullptr;
   TH1D *Resolution_EP1_exp   = nullptr;

   TH1D **Lpolar_y;            // P_{y} component of polarization vector for full hyperons (primary+secondary)  
   TH1D **Lpolar_y_prim;       // P_{y} component of polarization vector for primary hyperons
   TH1D **PstarRP_hist;        // Angular distributions (delta phi* of proton) for full hyperons from MCTracks (w.r.t. reaction plane)
   TH1D **PstarRP_hist_prim;   // Angular distributions (delta phi* of proton) for primary hyperons from MCTracks (w.r.t. reaction plane)
   TH1D **PstarEP_hist;        // Angular distributions (delta phi* of proton) for full hyperons from MCTracks (w.r.t. event plane)
   TH1D **PstarEP_hist_prim;   // Angular distributions (delta phi* of proton) for primary hyperons from MCTracks (w.r.t. event plane)

   //general parameters:
   float pi             = TMath::Pi();
   int pdgCodePr        = 2212; // pdg of proton
   int pdgCodeAPr       = -2212;// pdg of antiproton
   int pdgCodeNeg       = -211; // pdg of pi-
   int pdgCodePos       = 211;  // pdg of pi+
   int pdgCodeL0        = 3122; // pdg of lambda (1.11568)
   int pdgCodeAL0       = -3122;// pdg of antilambda (1.11568)

   double massPr      = 0.938272;  // mass of proton
   double massL0      = 1.11568;   // mass of lambda (1.11568)

   int pdgCodeHyperon;       // pdg of analyzed hyperon (should be set in the config file as particle_choice)
   int pdgCodeDaughter;      // pdg of daughter particle (proton in case of Lambda, antiproton in case of antiLambda)
   double massDaughter;      // mass of daughter particle
   double massHyperon;       // mass of analyzed hyperon

   double phiRP_mc;        // RP angle 
   double phiEP_mc;        // 1st-order EP angle (from FHCal)
   double ResEP_mc;        // cos(PhiEP_FHCal_Full - PhiRP) for calculation of true 1st-order EP resolution
   double ResEPSub_mc;     // cos(PhiEP_FHCal_S - PhiEP_FHCal_N) for calculation of reconstructed 1st-order EP resolution 
   double Centrality_tpc;  // Centrality from TPC

   
   double *SubEvRes1;      // cos(PhiEP_FHCal_S - PhiEP_FHCal_N) for calculation of reconstructed 1st-order EP resolution 
   double *ResEP1_true; 	// cos(PhiEP_FHCal_Full - PhiRP) for calculation of true 1st-order EP resolution

   int *centrality_min;
   int *centrality_max;
   double *_CentrBins;

   /**
    * @brief Select or reject event (implement different cuts)
    * 
    * @param event    event
    */
   bool selectEvent(MpdAnalysisEvent &event);

   /**
    * @brief Fill the necessary histograms for the output
    * 
    * @param event    event
    */
   void fillHistograms(MpdAnalysisEvent &event);
   
   /**
    * @brief Find angle of proton in the lambda frame
    * 
    * @param vPr                Proton vector
    * @param vLamb              Lambda vector
    * @param cos_prot           cos theta of proton
    * @param phi_prot_star      phi of proton
    */
   void FindPolarAngle(TVector3 &vPr, TVector3 &vLamb, double &cos_prot, double &phi_prot_star);
   double* init_double_array (const int n, const double fmt...);
   int* init_int_array (const int n, const int fmt...);

};
#endif
