#ifndef MPDGLOBALPOLARIZATIONRECO_H
#define MPDGLOBALPOLARIZATIONRECO_H

#include <deque>

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "MpdEvent.h"
#include "MpdVertex.h"
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
/**
 * @brief Wagon for analysis of RECO global polarization
 *
 */
class MpdGlobalPolarizationRECO : public MpdAnalysisTask2 {
   ClassDef(MpdGlobalPolarizationRECO, 1);

public:
   MpdGlobalPolarizationRECO() {}
   MpdGlobalPolarizationRECO(const char *name, const char *outputName, const char *analysis_choice,
                             const char *selection_choice);
   ~MpdGlobalPolarizationRECO() {} // TODO: normal descructor with cleaning off histos

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

private:
   std::string analysis_choice;  // analysis choice (selection or analysis)
   std::string selection_choice; // selection choice (selection using dca/chi/omega)

   // Event selection cuts
   float  mZvtxCut;        //(V) event selection cut (cm)
   int    NITER_CENT;      // number of centrality bins (defined for 4, 7, 10)
   int    NITER;           // number of angular bins
   int    NITER_ETA;       // number of eta bins
   int    NITER_PT;        // number of pt bins
   int    cent_cut_choice; // choice of centrality cut (0 - no cent cut, 1 - with cent cut)
   double cent_cut;        // value for centrality cut
   int    particle_choice; // particle choice (pdg of hyperon) --- currently Lambda (3122) or anti-Lambda (-3122))
   int    nMix;            // number of events to mix

   // Track selection cuts
   int mNofHitsCut; //(V) minimal number of hits to accept track

   // PID parameters
   double sigM;   // sigma for M
   double sigE;   // sigma for E
   double energy; // energy
   double coef;   // coefficient

   std::string generator; // generator name
   std::string tracking;  // tracking name
   std::string MCFile;    // MC file with geometry

   // Parameters for topological selection
   int         NITER_Selections;   // number of bins for selection parameters
   double      omega_start;        // starting value for omega_2
   double      omega_step;         // step value for chi_p
   double      chi_pi_start;       // starting value for chi_pi (> this)
   double      chi_p_start;        // starting value for chi_p (> this)
   double      chi_V0_start;       // starting value for chi_V0 (< this)
   double      lambda_path_start;  // starting value for lambda_path (> this)
   double      lambda_angle_start; // starting value for lambda_angle (< this)
   double      chi_pi_step;        // step value for chi_pi
   double      chi_p_step;         // step value for chi_p
   double      chi_V0_step;        // step value for chi_V0
   double      lambda_path_step;   // step value for lambda_path
   double      lambda_angle_step;  // step value for lambda_angle
   std::string selections_values;  // file with selection values
   double     *omega_value;

   // Parameters for selection from topology selection file (MB)
   double omega_value_full;
   double chi_pi_value_full;
   double chi_p_value_full;
   double chi_V0_value_full;
   double lambda_path_value_full;
   double lambda_angle_value_full;

   // Event properties
   TVector3 mPrimaryVertex;

   TClonesArray *mMCTracks     = nullptr;
   TClonesArray *mKalmanTracks = nullptr;
   TClonesArray *ZDCHits       = nullptr;
   TClonesArray *tpcPoints     = nullptr;
   TClonesArray *fTofMatches   = nullptr;

   MpdTpcKalmanFilter *recoTpc = nullptr;
   BaseTpcSectorGeo   *secGeo;
   TChain             *simMC;
   TBranch            *tpcSimB;
   MpdVertex          *fMpdVert;
   MpdPid             *fPid;

   // Histograms
   TList mHistoList;

   // General QA (analysis)
   TH1D *hEvents              = nullptr;
   TH1F *hVertex              = nullptr;
   TH1F *hCentrality          = nullptr;
   TH1D *hNevCentr            = nullptr;
   TH1D *hResolution_EP1_true = nullptr;
   TH1D *hResolution_EP1_reco = nullptr;
   TH1D *hMassL               = nullptr;
   TH1D *hMassLsig            = nullptr;
   TH1D *hMassLbkg            = nullptr;
   TH1D *hPIDflag             = nullptr;
   TH1D *hLambFlag            = nullptr;
   TH1D *hXiFlag              = nullptr;
   TH1D *hPtProt              = nullptr;
   TH1D *hPtProtT             = nullptr;
   TH1D *hPtProtF             = nullptr;

   TH1D *hHypRecopT_Prim = nullptr;
   TH1D *hHypRecopT_Full = nullptr;
   TH1D *hHypRecoEta_Prim = nullptr;
   TH1D *hHypRecoEta_Full = nullptr;
   
   TH1D *hHypMCpT_Prim = nullptr;
   TH1D *hHypMCEta_Prim = nullptr;
   
   TH2F *hHypRecopTEta_Full        = nullptr;
   TH2F *hHypRecopTEta_Prim        = nullptr;
   TH2F *hHypMCpTEta        = nullptr;

   // Tree for chi selection
   TTree *results_tree;

   // Invariant mass histograms
   TH1D **hm0_full;       // Omega2 selection: Full MB invariant mass histograms, after selection
   TH1D **hm0_full_mixed; // Omega2 selection: Full MB invariant mass histograms (mixed background), after selection

   TH1D *hm0_Full;            // Omega2/Chi analysis: Full MB invariant mass histograms, after selection
   TH1D *hm0_before_full;     // Omega2/Chi selection/analysis: Full MB invariant mass histograms, before selection
   TH1D *hm0_before_full_mix; // Omega2/Chi selection/analysis: Full MB invariant mass histograms, before selection
                              // (mixed background)

   TH1D ***hm0; // Omega2/Chi analysis: Invariant mass histograms in bins of delta phi* of proton for each centrality
                // bin, after selection
   TH1D ***hm0_mixed; // Omega2/Chi analysis: Invariant mass histograms in bins of delta phi* of proton for each
                      // centrality bin, after selection (mixed background)
   TH1D **hm0_before; // Omega2/Chi analysis: Invariant mass histograms in centrality bins, before selection
   TH1D **hm0_after;  // Omega2/Chi analysis: Invariant mass histograms in centrality bins, after selection

   // Histograms for polarization determination
   TH1D **hPolarY_Full;     // P_{y} component of polarization vector for full hyperons (primary+secondary)
   TH1D **hPolarY_Prim;     // P_{y} component of polarization vector for primary hyperons
   TH1D **hDeltaPhiEP_Full; // Angular distributions (delta phi* of proton, where phi* is reconstructed angle) for full
                            // hyperons from associated MCTracks (w.r.t. event plane)
   TH1D **hDeltaPhiEP_Prim; // Angular distributions (delta phi* of proton, where phi* is reconstructed angle) for
                            // primary hyperons from associated MCTracks (w.r.t. event plane)
   TH1D **hDeltaPhiRP_Full; // Angular distributions (delta phi* of proton, where phi* is reconstructed angle) for full
                            // hyperons from associated MCTracks (w.r.t. reaction plane)
   TH1D **hDeltaPhiRP_Prim; // Angular distributions (delta phi* of proton, where phi* is reconstructed angle) for
                            // primary hyperons from associated MCTracks (w.r.t. reaction plane)
   TH1D **hDeltaPhiRP_MC_Full; // Angular distributions (delta phi* of proton, where phi* is MC angle) for full hyperons
                               // from associated MCTracks (w.r.t. reaction plane)
   TH1D **hDeltaPhiRP_MC_Prim; // Angular distributions (delta phi* of proton, where phi* is MC angle) for primary
                               // hyperons from associated MCTracks (w.r.t. reaction plane)

   // Histograms to check pt and eta dependence of polarization
   TProfile **hPolvsPt;  // P_{y} component of polarization vector for full hyperons (primary+secondary) vs transverse
                         // momentum p_{T} of hyperon
   TProfile **hPolvsEta; // P_{y} component of polarization vector for full hyperons (primary+secondary) vs
                         // pseudorapidity eta of hyperon
   TH1D ***hm0_ptbin;  // Omega2/Chi analysis: Invariant mass histograms in bins of delta phi* of proton for each p_{T}
                       // bin for centrality interval 20-50%, after selection
   TH1D ***hm0_etabin; // Omega2/Chi analysis: Invariant mass histograms in bins of delta phi* of proton for each eta
                       // bin for centrality interval 20-50%, after selection
   TH1D ***hm0_ptbin_mixed;  // Omega2/Chi analysis: Invariant mass histograms in bins of delta phi* of proton for each
                             // p_{T} bin for centrality interval 20-50%, after selection (mixed background)
   TH1D ***hm0_etabin_mixed; // Omega2/Chi analysis: Invariant mass histograms in bins of delta phi* of proton for each
                             // eta bin for centrality interval 20-50%, after selection (mixed background)

   TH1D ***hm0_pt_eta_bin;
   TH1D ***hPolarY_Full_pt_eta_bin;
   TH1D ***hPolarY_Prim_pt_eta_bin;
   TH1D ***hDeltaPhiEP_Full_pt_eta_bin;
   TH1D ***hDeltaPhiRP_Full_pt_eta_bin;
   TH1D ***hDeltaPhiEP_Prim_pt_eta_bin;
   TH1D ***hDeltaPhiRP_Prim_pt_eta_bin;
   TProfile ***hv1EP_Full_pt_eta_bin;
   TProfile ***hv1RP_Full_pt_eta_bin;
   TProfile ***hv1EP_Prim_pt_eta_bin;
   TProfile ***hv1RP_Prim_pt_eta_bin;
   // general parameters:
   double pi         = TMath::Pi();
   int    pdgCodePr  = 2212;  // pdg of proton
   int    pdgCodeAPr = -2212; // pdg of antiproton
   int    pdgCodeNeg = -211;  // pdg of pi-
   int    pdgCodePos = 211;   // pdg of pi+
   int    pdgCodeL0  = 3122;  // pdg of lambda (1.11568)
   int    pdgCodeAL0 = -3122; // pdg of antilambda (1.11568)
   int    pdgCodeXi  = 3312;  // pdg of Xi- (1.3213 GeV/c)
   int    pdgCodeAXi = -3312; // pdg of antiXi- (1.3213 GeV/c)
   int    pdgCodeKm  = -321;  // pdg of K-

   double massPr      = 0.938272; // mass of proton
   double massPi      = 0.139570; // mass of pion
   double massL0      = 1.11568;  // mass of lambda (1.11568)
   double ptbin_min   = 0.0;      // ptbin_min
   double ptbin_max   = 3.5;      // ptbin_max
   double ptbin_step  = 0.5;      // ptbin_step
   double etabin_min  = -1.5;     // etabin_min
   double etabin_max  = 1.5;      // etabin_max
   double etabin_step = 0.5;      // etabin_step

   int pdgCodeHyperon;        // pdg of analyzed hyperon (should be set in the config file as particle_choice)
   int pdgCodeDaughterBar;    // pdg of daughter particle baryon (proton in case of Lambda, antiproton in case of
                              // antiLambda)
   int    pdgCodeDaughterMes; // pdg of daughter particle meson (pi- in case of Lambda, pi+ in case of antiLambda)
   int    pdgCodeMotherHyp;   // pdg of mother of hyperon (Xi)
   double massHyperon;        // mass of analyzed hyperon
   double massDaughterBar;    // mass of daughter particle
   double massDaughterMes;    // mass of daughter particle

   // defining the parameters for the analysis:
   const double gC2p  = 3.; // 4.; //9.;           //chi2 of p to PV
   const double gC2pi = 5.; // 5.; //11.;          //chi2 of pion to PV

   const double gDCAp  = 0.;  // cm - DCA of p to PV
   const double gDCApi = 0.;  // cm - DCA of pion to PV
   const double gPathL = 0.0; // cm - path to Lambda decay
   const double gC2L   = 25.; // 9999.;             //chi2 between pion & p in V0

   double xmin_anglemin = 0.;
   double xmax_anglemax = 2. * pi;

   double phiRP;          // RP angle
   double phiEP;          // 1st-order EP angle (from FHCal)
   double ResEP;          // cos(PhiEP_FHCal_Full - PhiRP) for calculation of true 1st-order EP resolution
   double ResEPSub;       // cos(PhiEP_FHCal_S - PhiEP_FHCal_N) for calculation of reconstructed 1st-order EP resolution
   double Centrality_tpc; // Centrality from TPC
   double b0;             // impact parameter

   double *SubEvRes1;   // cos(PhiEP_FHCal_S - PhiEP_FHCal_N) for calculation of reconstructed 1st-order EP resolution
   double *ResEP1_true; // cos(PhiEP_FHCal_Full - PhiRP) for calculation of true 1st-order EP resolution

   int idMax;

   double *angle_min; // minimum values for the range of angular distributions for polarization determination
   double *angle_max; // maximum values for the range of angular distributions for polarization determination

   // define the variables for the lambda object
   float massh, pth, ph, etah, yh, chi2h, disth, path, c2pv;
   float etas[2], mcps[2], ps[2], pts[2], chi2s[2], dcas[2], c2s[2], probs[2], chi2sL[2], dcasL[2], angle;
   float mcthetas[2], thetas[2], mcphis[2], phis[2];
   float dca, omega1, omega2, omegaL[3], cosA, cosAmc, polarhx, polarhy, polarhz, phi_star, phi_star_MC, phi_Lam;
   int   fEvNo, origs[2], qs[2], dstNo[2], layMx[2], ntr, nLamb, nLamb_MC;

   vector<vector<double>>                        vecL1;
   vector<pair<double, double>>                  fVecL1, fVecL2;
   std::vector<MpdLambdaPol>                     vLambdas;
   std::vector<MpdLambdaPol>                    *fvvvL;
   std::vector<tuple<int, float, float, float>>  fvLambMpdgPtEtaY;
   std::vector<tuple<int, float, float, float>> *fvvvLpt;
   std::vector<tuple<int, float, float, float>>  fvXiMpdgPtEtaY;
   std::vector<tuple<int, float, float, float>> *fvvvXipt;
   multimap<int, MpdTpcKalmanTrack>              fMapPiEvent;     // for event mixing
   map<int, MpdVertex>                           fMapVertexEvent; // for event mixing
   map<int, int>                                 ids, moths, pdgs, fLays;
   map<int, double>                              pots, ths, rads;

   int     moth, pdg;
   int    *centrality_min;
   int    *centrality_max;
   double *_CentrBins;

   double *pt_min;
   double *pt_max;
   double *eta_min;
   double *eta_max;

   // functions to initialize arrays:
   double *init_double_array(const int n, const double fmt...);
   int    *init_int_array(const int n, const int fmt...);

   /**
    * @brief Select or reject event (implement different cuts)
    *
    * @param event    event
    */
   bool selectEvent(MpdAnalysisEvent &event);

   /**
    * @brief Calculate MC information for Lambda and Xi hyperons and save it as ntuple (pdg,pt,eta,y)
    *
    * @param event    event
    */
   void ParticleMCProperties(MpdAnalysisEvent &event);

   /**
    * @brief Check tracks: find maximal reached layer No., exclude clones, find maximal track id
    *
    * @param event    event
    */
   void CheckTracks(MpdAnalysisEvent &event);

   /**
    * @brief Fill the necessary histograms for the output
    *
    * @param event    event
    */
   void fillHistograms(MpdAnalysisEvent &event);

   /**
    * @brief Calculate Lambda and Xi acceptance from MC data and save in maps, fill Flags for Lambda and Xi
    *
    * @param event    event
    */
   void CalculateLambdaAcceptance(MpdAnalysisEvent &event);

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
    * @brief Reconstruction Efficiency, if no PID. Otherwise fills Flags to check PID influence
    *
    * @param vecP       Proton vector
    * @param vecPi      Pion vector
    * @param use_pid    switch to toggle if particle identification is provided. Does some other stuff if it isnt
    */
   void RecoEff(vector<int> &vecP, vector<int> &vecPi, bool use_pid = false);

   /**
    * @brief Apply PID to get protons and pions (vecP and vecPi)
    *
    * @param vecP      Proton vector
    * @param vecPi     Pion vector
    */
   void ApplyPid(vector<int> &vecP, vector<int> &vecPi);

   /**
    * @brief Construct Lambda candidates (proton/pion pairs) as Lambda hyperon object
    *
    * @param vecP    Proton vector
    * @param vecPi   Pion vector
    * @param vecL    Lambda vector
    * @param phiEP   Event plane angle
    */
   void BuildLambda(vector<int> &vecP, vector<int> &vecPi, vector<MpdParticle *> &vecL, double &phiEP);

   /**
    * @brief Collect particles (pions, protons) from MpdTpcKalmanTrack checking the information from associated MC track
    *
    * @param vecPi      Pion vector
    * @param vecK       Kaon vector
    * @param vecP       Proton vector
    */
   void CollectParticles(vector<int> &vecPi, vector<int> &vecK, vector<int> &vecP);

   /**
    * @brief Find azimuthal angle of proton in the lambda frame
    *
    * @param lamb       Lambda vector
    * @param vPart      Proton vector
    */
   void FindPolarAngle(MpdParticle &lamb, vector<MpdParticle *> &vPart);
};
#endif
