#ifndef MPDGLOBALPOLARIZATION_H
#define MPDGLOBALPOLARIZATION_H

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

class MpdGlobalPolarization : public MpdAnalysisTask2 {
   ClassDef(MpdGlobalPolarization, 1);

public:
   MpdGlobalPolarization() {}
   MpdGlobalPolarization(const char *name, const char *outputName = "taskName");
   ~MpdGlobalPolarization() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

private:
   // Event selection cuts
   float mZvtxCut;    //(V) event selection cut (cm)

   // Track selection cuts (corresponding to centrality determination)
   int   mNofHitsCut;   //(V) minimal number of hits to accept track
   float mEtaCut;  //(V) maximal pseudorapidity accepted
   float mPtminCut;  //(V) minimal pt used in analysis
   float mDcaCut;  //(V) maximal DCA accepted

   int NITER_CENT;    // number of centrality bins (defined for 4, 7, 10)
   int NITER;   // number of angular bins
   int cent_cut_choice;    // choice of centrality cut (0 - no cent cut, 1 - with cent cut)
   double cent_cut;   // value for centrality cut

   std::string particle_choice; // particle choice (Lambda or anti-Lambda)

   // Event properties
   TVector3 mPrimaryVertex;
   std::string mOutFile = "histos.root";

   TClonesArray         *mMCTracks        = nullptr;
   TClonesArray         *ZDCHits          = nullptr;

   // Histograms
   TList mHistoList;

   // General QA
   TH1D *mhEvents             = nullptr;
   TH1F *mhVertex             = nullptr;
   TH1F *mhCentrality         = nullptr;
   TH1D *NCentr               = nullptr;
   TH1D *Resolution_EP1_true  = nullptr;
   TH1D *Resolution_EP1_exp   = nullptr;

   TProfile **Lpolar_y;
   TProfile **Lpolar_y_prim;
   TH1D **PstarRP_hist;
   TH1D **PstarRP_hist_prim;
   TH1D **PstarEP_hist;
   TH1D **PstarEP_hist_prim;

   //general parameters:
   float pi             = TMath::Pi();
   int _N_MODULES_TOTAL = 90;
   int _N_ARM           = 2;
   int pdgCodePr        = 2212; // pdg of proton
   int pdgCodeAPr       = -2212;// pdg of antiproton
   int pdgCodeNeg       = -211; // pdg of pi-
   int pdgCodePos       = 211;  // pdg of pi+
   int pdgCodeL0        = 3122; // pdg of lambda (1.11568)
   int pdgCodeAL0       = -3122;// pdg of antilambda (1.11568)

   static constexpr short nMixEventCent = 10;  //(V) number of bins in centrality for mixing
   int mCenBin = 0;

   double b0, phiRP_mc, phiEP_mc, ResEP_mc, ResEPSub_mc, Centrality_tpc;
   double *SubEvRes1, *ResEP1_true; 	

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
