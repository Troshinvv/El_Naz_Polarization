#ifndef MPDGLOBALPOLARIZATIONRECO_H
#define MPDGLOBALPOLARIZATIONRECO_H

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
#include "MpdTpcKalmanFilter.h"
#include "MpdKalmanTrack.h"
#include "TpcSectorGeoAZ.h"
#include "BaseTpcSectorGeo.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdParticle.h"

#include "MpdAnalysisTask2.h"
#include "MpdLambdaPol.h"
#include "hyperLinkDef.h"

class MpdGlobalPolarizationRECO : public MpdAnalysisTask2 {
   ClassDef(MpdGlobalPolarizationRECO, 1);

public:
   MpdGlobalPolarizationRECO() {}
   MpdGlobalPolarizationRECO(const char *name, const char *outputName, const char *analysis_choice, const char *selection_choice);
   ~MpdGlobalPolarizationRECO() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

private:
   std::string analysis_choice;  // analysis choice (selection or analysis)
   std::string selection_choice; // selection choice (selection using dca/chi/omega)

   // Event selection cuts
   float mZvtxCut;      //(V) event selection cut (cm)
   int NITER_CENT;      // number of centrality bins (defined for 4, 7, 10)
   int NITER;           // number of angular bins
   int cent_cut_choice; // choice of centrality cut (0 - no cent cut, 1 - with cent cut)
   double cent_cut;     // value for centrality cut
   int nMix;            // number of events to mix

   // Track selection cuts (corresponding to centrality determination)
   int   mNofHitsCut;   //(V) minimal number of hits to accept track
   float mEtaCut;       //(V) maximal pseudorapidity accepted
   float mPtminCut;     //(V) minimal pt used in analysis
   float mDcaCut;       //(V) maximal DCA accepted

   // PID parameters
   double sigM;         // sigma for M
   double sigE;         // sigma for E
   double energy;       // energy
   double coef;         // coefficient   

   std::string particle_choice; // particle choice (Lambda or anti-Lambda)
   std::string generator;       // generator name
   std::string tracking;        // tracking name
   std::string MCFile;          // MC file with geometry

   int NITER_Selections;      // number of bins for selection parameters
   double omega_start;        // starting value for omega_2
   double omega_step;         // step value for chi_p

   double chi_pi_start;       // starting value for chi_pi (> this)
   double chi_p_start;        // starting value for chi_p (> this)
   double chi_V0_start;       // starting value for chi_V0 (< this)
   double lambda_path_start;  // starting value for lambda_path (> this)
   double lambda_angle_start; // starting value for lambda_angle (< this)

   double chi_pi_step;        // step value for chi_pi
   double chi_p_step;         // step value for chi_p
   double chi_V0_step;        // step value for chi_V0
   double lambda_path_step;   // step value for lambda_path
   double lambda_angle_step;  // step value for lambda_angle

   std::string selections_values;  // file with selection values

   // Event properties
   TVector3 mPrimaryVertex;

   TClonesArray         *mMCTracks        = nullptr;
   TClonesArray         *mKalmanTracks    = nullptr;
   TClonesArray         *mMpdGlobalTracks = nullptr;
   TClonesArray         *ZDCHits          = nullptr;
   TClonesArray         *tpcPoints        = nullptr;
   TClonesArray         *fTofMatches      = nullptr;

   MpdTpcKalmanFilter   *recoTpc          = nullptr;
   BaseTpcSectorGeo *secGeo;
   TChain *simMC;
   TBranch *tpcSimB;

   // Histograms
   TList mHistoList;

   // General QA
   TH1D *mhEvents             = nullptr;
   TH1F *mhVertex             = nullptr;
   TH1F *mhCentrality         = nullptr;
   TH1D *NCentr               = nullptr;
   TH1D *Resolution_EP1_true  = nullptr;
   TH1D *Resolution_EP1_exp   = nullptr;
   TH1D *hMassL               = nullptr;
   TH1D *hMassLsig            = nullptr;
   TH1D *hMassLbkg            = nullptr;
   TH1D *Lpolar_full          = nullptr;
   TH1D *hPIDflag             = nullptr;  
   TH1D *hLambFlag            = nullptr;  
   TH1D *hXiFlag              = nullptr;  
   TH1D *hPtProt              = nullptr;  
   TH1D *hPtProtT             = nullptr; 
   TH1D *hPtProtF             = nullptr; 

   TTree *results;

   TH1D **hm0_full;
   TH1D *hm0_Full;
   TH1D *hm0_before_full;
   TH1D ***hm0;
   TH1D **hm0_before;
   TH1D **hm0_after;

   TH1D **Lpolar;
   TH1D **Lpolar_prim;
   TH1D **PstarEP_hist;
   TH1D **PstarEP_hist_prim;
   TH1D **PstarRP_hist;
   TH1D **PstarRP_hist_prim;
   TH1D **PstarRP_hist_MC;
   TH1D **PstarRP_hist_MC_prim;
   TH1D **Dca_pion;
   TH1D **Dca_proton;
   TH1D **Chi_pion;
   TH1D **Chi_proton;
   TH1D **Dca_lambda;
   TH1D **Chi_lambda;
   TH1D **Dca_v0;
   TH1D **Chi_v0;
   TH1D **Path_hist;
   TH1D **Angle_hist; 
	

   //general parameters:
   double pi             = TMath::Pi();
   int _N_MODULES_TOTAL = 90;
   int _N_ARM           = 2;
   int pdgCodePr        = 2212; // pdg of proton
   int pdgCodeAPr       = -2212;// pdg of antiproton
   int pdgCodeNeg       = -211; // pdg of pi-
   int pdgCodePos       = 211;  // pdg of pi+
   int pdgCodeL0        = 3122; // pdg of lambda (1.11568)
   int pdgCodeAL0       = -3122;// pdg of antilambda (1.11568)
   int pdgCodeXi        = 3312; // pdg of Xi- (1.3213 GeV/c)
   int pdgCodeAXi       = -3312; // pdg of antiXi- (1.3213 GeV/c)
   int pdgCodeKm        = -321; // pdg of K-

   //defining the parameters for the analysis (should check this at some point):
   const Double_t gC2p      = 3.; //4.; //9.;           //chi2 of p to PV
   const Double_t gC2pi     = 5.; //5.; //11.;          //chi2 of pion to PV

   const Double_t gDCAp     = 0.;           //cm - DCA of p to PV
   const Double_t gDCApi    = 0.;           //cm - DCA of pion to PV
   const Double_t gPathL    = 0.0;          //cm - path to Lambda decay
   const Double_t gC2L      = 25.; //9999.;  //chi2 between pion & p in V0

   double xmin_anglemin = 0.;
	double xmax_anglemax = 2.*pi;

   int idMax;
   double b0, phiRP_mc, phiEP_mc, ResEP_mc, ResEPSub_mc, Centrality_tpc;
   double omega_value_full;
   double *SubEvRes1, *ResEP1_true; 	
   double *omega_value;
   double *chi_pi_value;
	double *chi_p_value;
	double *chi_V0_value;
	double *lambda_path_value;
	double *lambda_angle_value;
   double *angle_min;
   double *angle_max;

   //define the variables for the tree:
   float massh, pth, ph, etah, yh, chi2h, disth, path, c2pv;
   float etas[2], mcps[2], ps[2], pts[2], chi2s[2], dcas[2], c2s[2], probs[2], chi2sL[2], dcasL[2],  angle;
   float mcthetas[2], thetas[2], mcphis[2], phis[2];
   float dca, omega1, omega2, omegaL[3], cosA, cosAmc, polarhx, polarhy, polarhz, phi_star, phi_star_MC, phi_Lam;
   int fEvNo, origs[2], qs[2], dstNo[2], layMx[2], ntr, nLamb, nLamb_MC;
   
   vector<vector<double> > vecL1;
   vector<pair<double,double> > fVecL1, fVecL2;
   std::vector<MpdLambdaPol> vLambdas;
   std::vector<tuple<int,float,float,float> > fvLambMpdgPtEtaY;
   std::vector<tuple<int,float,float,float> > *vvvLpt = &fvLambMpdgPtEtaY;
   std::vector<tuple<int,float,float,float> > fvXiMpdgPtEtaY;
   std::vector<tuple<int,float,float,float> > *fvvvXipt = &fvXiMpdgPtEtaY;
   multimap<int, MpdTpcKalmanTrack> fMapPiEvent; // for event mixing
   map<int,MpdVertex> fMapVertexEvent; // for event mixing
   map<int,int> ids, moths, pdgs, fLays;
   map<int,double> pots, ths, rads;
   MpdVertex *fMpdVert;
   MpdPid *fPid;

   int moth, pdg;
   int *centrality_min;
   int *centrality_max;
   double *_CentrBins;
   
   //functions to initialize arrays:
   double* init_double_array (const int n, const double fmt...);
   int* init_int_array (const int n, const int fmt...);

   /**
    * @brief Select or reject event (implement different cuts)
    * 
    * @param event    event
    */
   bool selectEvent(MpdAnalysisEvent &event);

   /**
    * @brief Calculate some particle properties for MC Lambda 
    * 
    * @param event    event
    */
   void ParticleMCProperties(MpdAnalysisEvent &event);

   /**
    * @brief Collect tracks
    * 
    * @param event    event
    */
   void CollectTracks(MpdAnalysisEvent &event);

   /**
    * @brief Fill the necessary histograms for the output
    * 
    * @param event    event
    */
   void fillHistograms(MpdAnalysisEvent &event);

   /**
    * @brief Calculate Lambda Acceptance
    * 
    * @param event    event
    * @param idMax    maxTPC track ID
    */
   void CalculateLambdaAcceptance(MpdAnalysisEvent &event, int idMax);

   /**
    * @brief Construct a helix
    * 
    * @param tr   track
    * @return MpdHelix  returns helix
    */
   MpdHelix MakeHelix(const MpdTpcKalmanTrack *tr);

   /**
    * @brief Construct a helix
    * 
    * @param part   particle
    * @return MpdHelix  returns helix
    */
   MpdHelix MakeHelix(const MpdParticle *part); 

   /**
    * @brief Reconstruction Efficiency
    * 
    * @param vecP       Proton vector
    * @param vecPi      Pion vector
    * @param use_pid    switch to toggle if particle identification is provided. Does some other stuff if it isnt
    */
   void RecoEff(vector<int> &vecP, vector<int> &vecPi, bool use_pid = false);

   /**
    * @brief Apply PID
    * 
    * @param pid       PID
    * @param vecP      Proton vector
    * @param vecPi     Pion vector 
    */
   void ApplyPid(MpdPid *pid, vector<int> &vecP, vector<int> &vecPi);

   /**
    * @brief construct Lambda candidates (proton/pion pairs)
    * 
    * @param vecP    Proton vector
    * @param vecPi   Pion vector
    * @param vecL    Lambda vector
    * @param phiEP   Event plane angle
    */
   void BuildLambda(vector<int> &vecP, vector<int> &vecPi, vector<MpdParticle*> &vecL, double &phiEP);
   
   /**
    * @brief Collect particles (pions, protons, kaons)
    * 
    * @param vecPi      Pion vector
    * @param vecK       Kaon vector
    * @param vecP       Proton vector
    */
   void CollectParticles(vector<int> &vecPi, vector<int> &vecK, vector<int> &vecP);

   /**
    * @brief Find angle of proton in the lambda frame
    * 
    * @param lamb       Lambda vector
    * @param vPart      Proton vector
    */
   void FindPolarAngle(MpdParticle &lamb, vector<MpdParticle*> &vPart); 
};
#endif
