#include <iostream>
#include <fstream> // std::ifstream

#include "MpdVertex.h"
#include "MpdEvent.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "MpdMCTrack.h"
#include "MpdTofMatchingData.h"
#include "MpdGlobalPolarizationRECO.h"
#include "MpdLambdaPol.h"

ClassImp(MpdGlobalPolarizationRECO);

MpdGlobalPolarizationRECO::MpdGlobalPolarizationRECO(const char *name, const char *outputName,
                                                     const char *analysis_choice, const char *selection_choice)
   : MpdAnalysisTask2(name, outputName)
{
   readParameters(name);
   param("mZvtxCut", mZvtxCut, 130.0);
   param("mNofHitsCut", mNofHitsCut, 10);
   param("NITER_CENT", NITER_CENT, 4);
   param("NITER", NITER, 20);
   param("NITER_ETA", NITER_ETA, 6);
   param("NITER_PT", NITER_PT, 5);
   param("cent_cut_choice", cent_cut_choice, 0);
   param("cent_cut", cent_cut, 70.0);
   param("particle_choice", particle_choice, 3122);
   param("nMix", nMix, 5);
   param("MCFile", MCFile, "");
   param("sigM", sigM, 3.0);
   param("sigE", sigE, 3.0);
   param("energy", energy, 9.2);
   param("coef", coef, 1.0);
   param("generator", generator, "PHSD");
   param("tracking", tracking, "CFHM");
   param("NITER_Selections", NITER_Selections, 30);
   param("omega_start", omega_start, 1.4);
   param("omega_step", omega_step, 0.1);
   param("chi_pi_start", chi_pi_start, 7.6);
   param("chi_p_start", chi_p_start, 4.2);
   param("chi_V0_start", chi_V0_start, 5.6);
   param("lambda_path_start", lambda_path_start, 1.6);
   param("lambda_angle_start", lambda_angle_start, 0.06);
   param("chi_pi_step", chi_pi_step, 0.2);
   param("chi_p_step", chi_p_step, 0.2);
   param("chi_V0_step", chi_V0_step, 0.2);
   param("lambda_path_step", lambda_path_step, 0.2);
   param("lambda_angle_step", lambda_angle_step, 0.02);
   param("selections_values", selections_values, "");

   this->analysis_choice  = analysis_choice;
   this->selection_choice = selection_choice;
}

void MpdGlobalPolarizationRECO::UserInit()
{
   cout << analysis_choice << " using " << selection_choice << " selection" << endl;
   // Initializing list of output histograms
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);
   fEvNo = -1;

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // Choice of analyzed particle
   pdgCodeHyperon = particle_choice;
   if (pdgCodeHyperon == pdgCodeL0) {
      cout << "You have chosen to analyze Lambda hyperons: "
           << " pdg: " << pdgCodeHyperon << endl;
      pdgCodeDaughterBar = pdgCodePr;
      pdgCodeDaughterMes = pdgCodeNeg;
      pdgCodeMotherHyp   = pdgCodeXi;
      massHyperon        = massL0;
      massDaughterBar    = massPr;
      massDaughterMes    = massPi;
   } else if (pdgCodeHyperon == pdgCodeAL0) {
      cout << "You have chosen to analyze anti-Lambda hyperons: "
           << " pdg: " << pdgCodeHyperon << endl;
      pdgCodeDaughterBar = pdgCodeAPr;
      pdgCodeDaughterMes = pdgCodePos;
      pdgCodeMotherHyp   = pdgCodeAXi;
      massHyperon        = massL0;
      massDaughterBar    = massPr;
      massDaughterMes    = massPi;
   } else {
      cout << "This pdg code for particle_choice is not defined! Please provide the definition in the code." << endl;
      exit(1);
   }
   cout << "massHyperon: " << massHyperon << "; massDaughterBar: " << massDaughterBar
        << "; massDaughterMes: " << massDaughterMes << endl;

   // Initializing centrality bins, dependent on the number of bins
   if (NITER_CENT == 4) {
      centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
      centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
      _CentrBins     = init_double_array(5, 0, 0., 10., 20., 50., 100.);
   } else if (NITER_CENT == 7) {
      centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
      centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
      _CentrBins     = init_double_array(8, 0, 0., 10., 20., 30., 40., 50., 60., 70.);
   } else if (NITER_CENT == 10) {
      centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
      centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
      _CentrBins     = init_double_array(11, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.);
   } else {
      cout << "This values of centrality bins is not defined! Please provide the definition in the code." << endl;
      exit(1);
   }

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      cout << "iter_cent = " << iter_cent << "; centrality_min = " << centrality_min[iter_cent]
           << "; _CentrBins = " << _CentrBins[iter_cent] << endl;
   }

   // Initializing pt bins
   if (NITER_PT == 5) {
      pt_min = init_double_array(5, 0, 0., 0.5, 1.0, 1.5, 2.0);
      pt_max = init_double_array(5, 0, 0.5, 1.0, 1.5, 2.0, 3.0);
   } else {
      cout << "This values of NITER_PT is not defined! Please put NITER_PT = 5." << endl;
      exit(1);
   }

   for (int iter_pt = 0; iter_pt < NITER_PT; iter_pt++) {
      cout << "iter_pt = " << iter_pt << "; pt_min = " << pt_min[iter_pt] << "; pt_max = " << pt_max[iter_pt] << endl;
   }

   // Initializing eta bins
   eta_min = new double[NITER_ETA];
   eta_max = new double[NITER_ETA];
   for (int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++) {
      eta_min[iter_eta] = etabin_min + etabin_step * iter_eta;
      eta_max[iter_eta] = etabin_min + etabin_step * (iter_eta + 1);
      cout << "iter_eta = " << iter_eta << "; eta_min = " << eta_min[iter_eta] << "; eta_max = " << eta_max[iter_eta]
           << endl;
   }

   // Initializing PID
   fPid = new MpdPid(sigM, sigE, energy, coef, generator, tracking, "pikaprdetrhe3he4");

   // Creating the necessary histograms
   // Only the ones required for "selection" or "analysis" will be saved

   hEvents     = new TH1D("hEvents", "Number of events", 10, 0., 10.);
   hVertex     = new TH1F("hVertex", "Event vertex distribution", 100, -200., 200.);
   hCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);

   hNevCentr            = new TH1D("hNevCentr", "Events in centrality bins", NITER_CENT, _CentrBins);
   hResolution_EP1_true = new TH1D("hResolution_EP1_true", "True EP1 resolution", NITER_CENT, _CentrBins);
   hResolution_EP1_reco = new TH1D("hResolution_EP1_reco", "Reco EP1 resolution", NITER_CENT, _CentrBins);

   hMassL    = new TH1D("hMassL", "Lambda mass", 50, 1.070, 1.170);
   hMassLsig = new TH1D("hMassLsig", "Lambda mass (signal)", 50, 1.070, 1.170);
   hMassLbkg = new TH1D("hMassLbkg", "Lambda mass (bckg.)", 50, 1.070, 1.170);
   hPIDflag  = new TH1D("hPIDflag", "PID flags", 12, 0, 12);
   hLambFlag = new TH1D("hLambFlag", "Flags for Lambda", 14, 0, 14);
   hXiFlag   = new TH1D("hXiFlag", "Flags for Xi", 14, 0, 14);
   hPtProt   = new TH1D("hPtProt", "Proton Pt", 20, 0, 5);
   hPtProtT  = new TH1D("hPtProtT", "True Proton Pt", 20, 0, 5);
   hPtProtF  = new TH1D("hPtProtF", "False Proton Pt", 20, 0, 5);

   hHypRecopT_Prim = new TH1D("hHypRecopT_Prim","Reco pT Primary", 300,0,3);
   hHypRecopT_Full = new TH1D("hHypRecopT_Full","Reco pT Full", 300,0,3);
   hHypRecoEta_Prim = new TH1D("hHypRecoEta_Prim","Reco Rapidity Primary", 300,-1.5,1.5);
   hHypRecoEta_Full = new TH1D("hHypRecoEta_Full","Reco Rapidity Full", 300,-1.5,1.5);
   hHypMCpT_Prim = new TH1D("hHypMCpT_Prim","MC pT Primary", 300,0,3);
   hHypMCEta_Prim = new TH1D("hHypMCEta_Prim","MC Rapidity Primary", 300,-1.5,1.5);

  
hHypRecopTEta_Full =new TH2F("hHypRecopTEta_Full","Reco pT vs Eta Full",600,-3,3,400,0,4); 
hHypRecopTEta_Prim =new TH2F("hHypRecopTEta_Prim","Reco pT vs Eta Prim",600,-3,3,400,0,4);
hHypMCpTEta =new TH2F("hHypMCpTEta","MC pT vs Eta",600,-3,3,400,0,4);

   fvvvL   = &vLambdas;
   fvvvLpt = &fvLambMpdgPtEtaY;

   results_tree = new TTree("event", "Event");
   results_tree->Branch("b0", &b0, "b0/D");                                                  // impact parameter
   results_tree->Branch("Centrality_tpc", &Centrality_tpc, "Centrality_tpc/D");              // event centrality
   TBranch *br = results_tree->Branch("l0", "std::vector<MpdLambdaPol>", &fvvvL);            // lambda candidates
   results_tree->Branch("ptetayl0", "std::vector<tuple<int,float,float,float> >", &fvvvLpt); // lambda phase space (MC)
   results_tree->Branch("nLamb", &nLamb, "nLamb/I");          // number of Lambda (in collection)
   results_tree->Branch("nLamb_MC", &nLamb_MC, "nLamb_MC/I"); // number of Lambda (MC)

   if (analysis_choice == "analysis") {
      fOutputList->Add(hEvents);
      fOutputList->Add(hVertex);
      fOutputList->Add(hCentrality);
      fOutputList->Add(hNevCentr);
      fOutputList->Add(hResolution_EP1_true);
      fOutputList->Add(hResolution_EP1_reco);
      fOutputList->Add(hMassL);
      fOutputList->Add(hMassLsig);
      fOutputList->Add(hMassLbkg);
      fOutputList->Add(hPIDflag);
      fOutputList->Add(hLambFlag);
      fOutputList->Add(hXiFlag);

      fOutputList->Add(hHypRecopT_Prim);
      fOutputList->Add(hHypRecopT_Full);
      fOutputList->Add(hHypRecoEta_Prim);
      fOutputList->Add(hHypRecoEta_Full);
      fOutputList->Add(hHypMCpT_Prim);
       fOutputList->Add(hHypMCEta_Prim);


      fOutputList->Add(hHypRecopTEta_Full);
      fOutputList->Add(hHypRecopTEta_Prim);
      fOutputList->Add(hHypMCpTEta);

      hm0_Full = new TH1D("hm0_Full", "Full (MB) Invariant Mass of Hyperon", 100, 1.07, 1.17);
      fOutputList->Add(hm0_Full);
      hm0_before_full =
         new TH1D("hm0_before_full", "Full (MB) Invariant Mass of Hyperon (before selection)", 100, 1.07, 1.17);
      fOutputList->Add(hm0_before_full);
      hm0_before          = new TH1D *[NITER_CENT];
      hm0_after           = new TH1D *[NITER_CENT];
      hPolarY_Full        = new TH1D *[NITER_CENT];
      hPolarY_Prim        = new TH1D *[NITER_CENT];
      hDeltaPhiEP_Full    = new TH1D *[NITER_CENT];
      hDeltaPhiEP_Prim    = new TH1D *[NITER_CENT];
      hDeltaPhiRP_Full    = new TH1D *[NITER_CENT];
      hDeltaPhiRP_Prim    = new TH1D *[NITER_CENT];
      hDeltaPhiRP_MC_Full = new TH1D *[NITER_CENT];
      hDeltaPhiRP_MC_Prim = new TH1D *[NITER_CENT];
      hm0                 = new TH1D **[NITER_CENT];
      hm0_mixed           = new TH1D **[NITER_CENT];
      angle_min           = new double[NITER];
      angle_max           = new double[NITER];

      hm0_ptbin        = new TH1D **[NITER_PT];
      hm0_etabin       = new TH1D **[NITER_ETA];
      hm0_ptbin_mixed  = new TH1D **[NITER_PT];
      hm0_etabin_mixed = new TH1D **[NITER_ETA];

      hm0_pt_eta_bin = new TH1D **[NITER_PT];
      hPolarY_Full_pt_eta_bin = new TH1D **[NITER_PT];
      hPolarY_Prim_pt_eta_bin = new TH1D **[NITER_PT];
      hDeltaPhiEP_Full_pt_eta_bin = new TH1D **[NITER_PT];
   hDeltaPhiRP_Full_pt_eta_bin = new TH1D **[NITER_PT];
   hDeltaPhiEP_Prim_pt_eta_bin = new TH1D **[NITER_PT];
   hDeltaPhiRP_Prim_pt_eta_bin = new TH1D **[NITER_PT];
   hv1EP_Full_pt_eta_bin = new TProfile **[NITER_PT];
   hv1RP_Full_pt_eta_bin = new TProfile **[NITER_PT];
   hv1EP_Prim_pt_eta_bin = new TProfile **[NITER_PT];
   hv1RP_Prim_pt_eta_bin = new TProfile **[NITER_PT];

      hPolvsPt  = new TProfile *[NITER_CENT];
      hPolvsEta = new TProfile *[NITER_CENT];

      double step_angle = (xmax_anglemax - xmin_anglemin) / NITER;
      cout << "xmin_anglemin = " << xmin_anglemin << "; xmax_anglemax = " << xmax_anglemax << endl;
      for (int iter = 0; iter < NITER; iter++) {
         angle_min[iter] = xmin_anglemin + step_angle * iter;
         angle_max[iter] = xmin_anglemin + step_angle * (iter + 1);
         cout << "iter = " << iter << "; angle_min = " << angle_min[iter] << "; angle_max = " << angle_max[iter]
              << endl;
      }

      // reading the optimal topology selection parameters from the file
      if (selection_choice == "omega2") {
         // reading the omega_2 values from the file:
         cout << "Topology selection using omega_2 parameter" << endl;
         ifstream selections_file;
         selections_file.open(selections_values);
         if (selections_file.fail()) {
            cout << "File with selection values does not exist! Please run the 'selection' choice first! Exiting... "
                 << endl;
            exit(1);
         }
         double test = -1.;
         selections_file >> omega_value_full >> test;
         selections_file.close();
         if (test != -1.) {
            cout << "File must have only one value, did you load the correct one? Exiting... " << endl;
            exit(1);
         }
         cout << "omega_value_full =  " << omega_value_full << endl;
      } else if (selection_choice == "chi") {
         // reading the chi selection values from the file:
         cout << "Topology selection using chi selection parameters" << endl;
         ifstream selections_file;
         selections_file.open(selections_values);
         if (selections_file.fail()) {
            cout << "File with selection values does not exist! Please run the 'selection' choice first! Exiting... "
                 << endl;
            exit(1);
         }
         chi_pi_value_full = chi_p_value_full = chi_V0_value_full = lambda_path_value_full = lambda_angle_value_full =
            -1;
         selections_file >> chi_pi_value_full >> chi_p_value_full >> chi_V0_value_full >> lambda_path_value_full >>
            lambda_angle_value_full;
         selections_file.close();
         if (lambda_angle_value_full == -1.) {
            cout << "File must have five values, did you load the correct one? Exiting... " << endl;
            exit(1);
         }
         cout << "chi_pi_value_full =  " << chi_pi_value_full << "; chi_p_value_full =  " << chi_p_value_full
              << "; chi_V0_value_full =  " << chi_V0_value_full
              << "; lambda_path_value_full =  " << lambda_path_value_full
              << "; lambda_angle_value_full =  " << lambda_angle_value_full << endl;
      } else {
         cout << "No such selection_choice defined yet! Please choose either 'omega2' or 'chi'! Exiting... " << endl;
         exit(1);
      }

      for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
         hm0[iter_cent]       = new TH1D *[NITER];
         hm0_mixed[iter_cent] = new TH1D *[NITER];
         hm0_before[iter_cent] =
            new TH1D(Form("hm0_before_%d", iter_cent), Form("hm0_before_%d", iter_cent), 100, 1.07, 1.17);
         fOutputList->Add(hm0_before[iter_cent]);
         hm0_after[iter_cent] =
            new TH1D(Form("hm0_after_%d", iter_cent), Form("hm0_after_%d", iter_cent), 100, 1.07, 1.17);
         fOutputList->Add(hm0_after[iter_cent]);
         hPolarY_Full[iter_cent] =
            new TH1D(Form("hPolarY_Full_%d", iter_cent), Form("hPolarY_Full_%d", iter_cent), 100, -1., 1.);
         fOutputList->Add(hPolarY_Full[iter_cent]);
         hPolarY_Prim[iter_cent] =
            new TH1D(Form("hPolarY_Prim_%d", iter_cent), Form("hPolarY_Prim_%d", iter_cent), 100, -1., 1.);
         fOutputList->Add(hPolarY_Prim[iter_cent]);
         hDeltaPhiEP_Full[iter_cent] = new TH1D(Form("hDeltaPhiEP_Full_%d", iter_cent),
                                                Form("hDeltaPhiEP_Full_%d", iter_cent), NITER, 0., 2. * pi);
         fOutputList->Add(hDeltaPhiEP_Full[iter_cent]);
         hDeltaPhiEP_Prim[iter_cent] = new TH1D(Form("hDeltaPhiEP_Prim_%d", iter_cent),
                                                Form("hDeltaPhiEP_Prim_%d", iter_cent), NITER, 0., 2. * pi);
         fOutputList->Add(hDeltaPhiEP_Prim[iter_cent]);
         hDeltaPhiRP_Full[iter_cent] = new TH1D(Form("hDeltaPhiRP_Full_%d", iter_cent),
                                                Form("hDeltaPhiRP_Full_%d", iter_cent), NITER, 0., 2. * pi);
         fOutputList->Add(hDeltaPhiRP_Full[iter_cent]);
         hDeltaPhiRP_Prim[iter_cent] = new TH1D(Form("hDeltaPhiRP_Prim_%d", iter_cent),
                                                Form("hDeltaPhiRP_Prim_%d", iter_cent), NITER, 0., 2. * pi);
         fOutputList->Add(hDeltaPhiRP_Prim[iter_cent]);
         hDeltaPhiRP_MC_Full[iter_cent] = new TH1D(Form("hDeltaPhiRP_MC_Full_%d", iter_cent),
                                                   Form("hDeltaPhiRP_MC_Full_%d", iter_cent), NITER, 0., 2. * pi);
         fOutputList->Add(hDeltaPhiRP_MC_Full[iter_cent]);
         hDeltaPhiRP_MC_Prim[iter_cent] = new TH1D(Form("hDeltaPhiRP_MC_Prim_%d", iter_cent),
                                                   Form("hDeltaPhiRP_MC_Prim_%d", iter_cent), NITER, 0., 2. * pi);
         fOutputList->Add(hDeltaPhiRP_MC_Prim[iter_cent]);
         hPolvsPt[iter_cent] =
            new TProfile(Form("hPolvsPt_%d", iter_cent), Form("hPolvsPt_%d", iter_cent), NITER, 0., 3.0);
         fOutputList->Add(hPolvsPt[iter_cent]);
         hPolvsEta[iter_cent] =
            new TProfile(Form("hPolvsEta_%d", iter_cent), Form("hPolvsEta_%d", iter_cent), NITER, -1.5, 1.5);
         fOutputList->Add(hPolvsEta[iter_cent]);

         for (int iter = 0; iter < NITER; iter++) {
            hm0[iter_cent][iter] =
               new TH1D(Form("hm0_%d_%d", iter_cent, iter), Form("hm0_%d_%d", iter_cent, iter), 100, 1.07, 1.17);
            fOutputList->Add(hm0[iter_cent][iter]);
            hm0_mixed[iter_cent][iter] = new TH1D(Form("hm0_mixed_%d_%d", iter_cent, iter),
                                                  Form("hm0_mixed%d_%d", iter_cent, iter), 100, 1.07, 1.17);
            fOutputList->Add(hm0_mixed[iter_cent][iter]);
         }
      }
      for (int iter_pt = 0; iter_pt < NITER_PT; iter_pt++) {
         hm0_ptbin[iter_pt]       = new TH1D *[NITER];
         hm0_ptbin_mixed[iter_pt] = new TH1D *[NITER];
         for (int iter = 0; iter < NITER; iter++) {
            hm0_ptbin[iter_pt][iter] = new TH1D(Form("hm0_ptbin_%d_%d", iter_pt, iter),
                                                Form("hm0_ptbin_%d_%d", iter_pt, iter), 100, 1.07, 1.17);
            fOutputList->Add(hm0_ptbin[iter_pt][iter]);
            hm0_ptbin_mixed[iter_pt][iter] = new TH1D(Form("hm0_ptbin_mixed_%d_%d", iter_pt, iter),
                                                      Form("hm0_ptbin_mixed_%d_%d", iter_pt, iter), 100, 1.07, 1.17);
            fOutputList->Add(hm0_ptbin_mixed[iter_pt][iter]);
         }
      }
      for (int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++) {
         hm0_etabin[iter_eta]       = new TH1D *[NITER];
         hm0_etabin_mixed[iter_eta] = new TH1D *[NITER];
         for (int iter = 0; iter < NITER; iter++) {
            hm0_etabin[iter_eta][iter] = new TH1D(Form("hm0_etabin_%d_%d", iter_eta, iter),
                                                  Form("hm0_etabin_%d_%d", iter_eta, iter), 100, 1.07, 1.17);
            fOutputList->Add(hm0_etabin[iter_eta][iter]);
            hm0_etabin_mixed[iter_eta][iter] =
               new TH1D(Form("hm0_etabin_mixed_%d_%d", iter_eta, iter), Form("hm0_etabin_mixed_%d_%d", iter_eta, iter),
                        100, 1.07, 1.17);
            fOutputList->Add(hm0_etabin_mixed[iter_eta][iter]);
         }
      }
      cout<< "my cycle"<<endl;
      for (int iter_pt = 0; iter_pt < NITER_PT; iter_pt++) {
           hm0_pt_eta_bin[iter_pt] = new TH1D *[NITER_ETA];
           hPolarY_Full_pt_eta_bin[iter_pt] = new TH1D *[NITER_ETA];
           hPolarY_Prim_pt_eta_bin[iter_pt] = new TH1D *[NITER_ETA];
           hDeltaPhiEP_Full_pt_eta_bin[iter_pt] = new TH1D *[NITER_ETA];
           hDeltaPhiRP_Full_pt_eta_bin[iter_pt] = new TH1D *[NITER_ETA];
   hDeltaPhiEP_Prim_pt_eta_bin[iter_pt] = new TH1D *[NITER_ETA];
   hDeltaPhiRP_Prim_pt_eta_bin[iter_pt] = new TH1D *[NITER_ETA];
   hv1EP_Full_pt_eta_bin[iter_pt] = new TProfile *[NITER_ETA];
   hv1RP_Full_pt_eta_bin[iter_pt] = new TProfile *[NITER_ETA];
   hv1EP_Prim_pt_eta_bin[iter_pt] = new TProfile *[NITER_ETA];
   hv1RP_Prim_pt_eta_bin[iter_pt] = new TProfile *[NITER_ETA];
         for (int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++) {
            hm0_pt_eta_bin[iter_pt][iter_eta] = new TH1D (Form("hm0_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hm0_pt_eta_bin%d_%d", iter_pt, iter_eta), 100, 1.07, 1.17);
            fOutputList->Add(hm0_pt_eta_bin[iter_pt][iter_eta]);
            hPolarY_Full_pt_eta_bin[iter_pt][iter_eta] = new TH1D (Form("hPolarY_Full_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hPolarY_Full_pt_eta_bin%d_%d", iter_pt, iter_eta), 100, -1., 1.);
            fOutputList->Add(hPolarY_Full_pt_eta_bin[iter_pt][iter_eta]);
            hPolarY_Prim_pt_eta_bin[iter_pt][iter_eta] = new TH1D (Form("hPolarY_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hPolarY_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta), 100, -1., 1.);
            fOutputList->Add(hPolarY_Prim_pt_eta_bin[iter_pt][iter_eta]);


            hDeltaPhiEP_Full_pt_eta_bin[iter_pt][iter_eta] = new TH1D (Form("hDeltaPhiEP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hDeltaPhiEP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta), NITER, 0., 2. * pi);
            fOutputList->Add(hDeltaPhiEP_Full_pt_eta_bin[iter_pt][iter_eta]);
            hDeltaPhiRP_Full_pt_eta_bin[iter_pt][iter_eta] = new TH1D (Form("hDeltaPhiRP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hDeltaPhiRP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta), NITER, 0., 2. * pi);
            fOutputList->Add(hDeltaPhiRP_Full_pt_eta_bin[iter_pt][iter_eta]);
            hDeltaPhiRP_Prim_pt_eta_bin[iter_pt][iter_eta] = new TH1D (Form("hDeltaPhiRP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hDeltaPhiRP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta), NITER, 0., 2. * pi);
            fOutputList->Add(hDeltaPhiRP_Prim_pt_eta_bin[iter_pt][iter_eta]);
            hDeltaPhiEP_Prim_pt_eta_bin[iter_pt][iter_eta] = new TH1D (Form("hDeltaPhiEP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hDeltaPhiEP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta), NITER, 0., 2. * pi);
            fOutputList->Add(hDeltaPhiEP_Prim_pt_eta_bin[iter_pt][iter_eta]);


            hv1EP_Prim_pt_eta_bin[iter_pt][iter_eta] = new TProfile (Form("hv1EP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hv1EP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta), 100, 1.07, 1.17);
            fOutputList->Add(hv1EP_Prim_pt_eta_bin[iter_pt][iter_eta]);
            hv1RP_Prim_pt_eta_bin[iter_pt][iter_eta] = new TProfile (Form("hv1RP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hv1RP_Prim_pt_eta_bin%d_%d", iter_pt, iter_eta), 100, 1.07, 1.17);
            fOutputList->Add(hv1RP_Prim_pt_eta_bin[iter_pt][iter_eta]);
            hv1EP_Full_pt_eta_bin[iter_pt][iter_eta] = new TProfile (Form("hv1EP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hv1EP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta), 100, 1.07, 1.17);
            fOutputList->Add(hv1EP_Full_pt_eta_bin[iter_pt][iter_eta]);
            hv1RP_Full_pt_eta_bin[iter_pt][iter_eta] = new TProfile (Form("hv1RP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta),
                                                  Form("hv1RP_Full_pt_eta_bin%d_%d", iter_pt, iter_eta), 100, 1.07, 1.17);
            fOutputList->Add(hv1RP_Full_pt_eta_bin[iter_pt][iter_eta]);
         } 
      }
cout<< "my cycle ends here"<<endl;

   } else if (analysis_choice == "selection") {
      if (selection_choice == "omega2") {
         cout << "Topology selection using omega_2 parameter" << endl;
         hm0_full        = new TH1D *[NITER_Selections];
         hm0_full_mixed  = new TH1D *[NITER_Selections];
         hm0_before_full = new TH1D("hm0_before_full", "hm0_before_full", 100, 1.07, 1.17);
         fOutputList->Add(hm0_before_full);
         hm0_before_full_mix = new TH1D("hm0_before_full_mix", "hm0_before_full_mix", 100, 1.07, 1.17);
         fOutputList->Add(hm0_before_full_mix);
         omega_value = new double[NITER_Selections];
         for (int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++) {
            omega_value[iter_sel] = omega_start + omega_step * iter_sel;
            cout << "iter_sel = " << iter_sel << "; omega_value = " << omega_value[iter_sel] << endl;
         }
         for (int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++) {
            hm0_full[iter_sel] =
               new TH1D(Form("hm0_full_%d", iter_sel), Form("hm0_full_%d", iter_sel), 100, 1.07, 1.17);
            fOutputList->Add(hm0_full[iter_sel]);
            hm0_full_mixed[iter_sel] =
               new TH1D(Form("hm0_full_mixed_%d", iter_sel), Form("hm0_full_mixed_%d", iter_sel), 100, 1.07, 1.17);
            fOutputList->Add(hm0_full_mixed[iter_sel]);
         }
      } else if (selection_choice == "chi") {
         cout << "Topology selection using chi parameters" << endl;
         // Just save the tree
         fOutputList->Add(results_tree);
      }
   } else {
      cout << "No such analysis_choice defined yet! Please choose either 'selection' or 'analysis'! Exiting... "
           << endl;
      exit(1);
   }

   ResEP1_true = new double[NITER_CENT];
   SubEvRes1   = new double[NITER_CENT];

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      ResEP1_true[iter_cent] = 0.;
      SubEvRes1[iter_cent]   = 0.;
   }
   phiRP    = 0.;
   phiEP    = 0.;
   ResEP    = 0.;
   ResEPSub = 0.;
}

void MpdGlobalPolarizationRECO::ProcessEvent(MpdAnalysisEvent &event)
{
   if (!selectEvent(event)) {
      return;
   }

   // Calculate the reaction/event plane and its resolution
   b0       = event.fMCEventHeader->GetB();
   phiRP    = event.fMCEventHeader->GetRotZ();
   phiRP    = TMath::ATan2(TMath::Sin(phiRP), TMath::Cos(phiRP));
   phiEP    = event.fMpdEP.GetPhiEP_FHCal_F_all();
   ResEP    = TMath::Cos(phiEP - phiRP);
   ResEPSub = TMath::Cos(event.fMpdEP.GetPhiEP_FHCal_S_all() - event.fMpdEP.GetPhiEP_FHCal_N_all());

   // The following reconstruction code is largely based on the Hyperons code by A. Zinchenko (Hyperons wagon)
   // Remove first stored event in case of event mixing
   if (nMix && fMapVertexEvent.size() > nMix) {
      int iev0 = fMapVertexEvent.begin()->first;
      fMapVertexEvent.erase(iev0);
      fMapPiEvent.erase(iev0);
   }
   fMapVertexEvent[fEvNo] = *fMpdVert;

   ParticleMCProperties(event);

   TArrayI *indxs = fMpdVert->GetIndices();
   int      nPrim = indxs->GetSize();
   set<int> indxVert;
   for (int k = 0; k < nPrim; k++) indxVert.insert((*indxs)[k]);
   if (fEvNo % 100 == 0) {
      cout << " *** Event No: " << fEvNo << ", reco tracks in TPC: "
           << " " << mKalmanTracks->GetEntriesFast() << ", vertices: " << event.fVertex->GetEntriesFast() << endl;
      cout << " Number of primary (used for vertex reco) tracks: " << indxVert.size() << endl;
   }

   CheckTracks(event);
   CalculateLambdaAcceptance(event);

   // Check Track quality
   for (int j = 0; j < mKalmanTracks->GetEntriesFast(); j++) {
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      int    id     = tr->GetTrackID();
      double etaRec = tr->Momentum3().Eta();
      if (TMath::Abs(fLays[id]) < -41 ||
          TMath::Abs(etaRec) > 13) // here's a placeholder for possible cuts on lays and eta
         tr->SetChi2(-9.);
      if (tr->GetNofHits() < mNofHitsCut && abs(tr->Momentum3().Eta())>1.5) tr->SetChi2(-9.);
   }

   // Collect "good" pions and protons
   vector<int> vecPi, vecK, vecP;
   CollectParticles(vecPi, vecK, vecP);
   if (fEvNo % 100 == 0) cout << " Number of protons(anti-protons), pi: " << vecP.size() << " " << vecPi.size() << endl;
   RecoEff(vecP, vecPi, 1);
   if (fEvNo % 100 == 0) cout << " Number of protons(anti-protons), pi: " << vecP.size() << " " << vecPi.size() << endl;

   // Apply PID to find reconstructed tracks
   ApplyPid(vecP, vecPi);
   vector<MpdParticle *> vecL;
   vecL.clear();
   vLambdas.clear();

   // Build Lambda hyperon object (construct Lambda candidates from proton-pion pairs)
   BuildLambda(vecP, vecPi, vecL, phiEP);
   nLamb = vecL.size();

   // Fill tree for chi selection
   results_tree->Fill();

   // Fill the necessary histograms
   fillHistograms(event);

   for (int ipart = 0; ipart < nLamb; ipart++) delete vecL[ipart];

   fvLambMpdgPtEtaY.clear();
   fvXiMpdgPtEtaY.clear();
}

void MpdGlobalPolarizationRECO::Finish()
{
   cout << "Finish() ..." << endl;
}

bool MpdGlobalPolarizationRECO::selectEvent(MpdAnalysisEvent &event)
{
   static bool isInitialized = false;
   fEvNo++;
   if (!isInitialized) {

      if (MCFile == "") {
         cout << "Parameter MCFile must be set in the config file!" << endl;
         exit(1);
      }
      cout << "Reading Geo from simulation file for track refit ... " << endl;
      simMC = new TChain("mpdsim");
      simMC->AddFile(TString(MCFile)); // using for now the one from Request 25 (same geometry)
      simMC->SetName("mpdsim1");
      TFile fileMC(simMC->GetListOfFiles()->First()->GetTitle());
      fileMC.Get("FairGeoParSet");
      TClonesArray *tpcPoints = (TClonesArray *)fileMC.FindObjectAny("TpcPoint");
      simMC->SetBranchAddress("TpcPoint", &tpcPoints);
      TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
      secGeo           = new TpcSectorGeoAZ();
      recoTpc          = new MpdTpcKalmanFilter(*secGeo, "TPC Kalman filter");
      recoTpc->FillGeoScheme();
      isInitialized = true;
   }

   hEvents->Fill(0.5); // Number of full events

   mMCTracks     = event.fMCTrack;
   mKalmanTracks = event.fTPCKalmanTrack;
   fMpdVert      = (MpdVertex *)event.fVertex->First();
   fMpdVert->Position(mPrimaryVertex);
   fTofMatches = event.fTOFMatching;

   int nMC   = mMCTracks->GetEntriesFast();
   int nTrMc = 0;
   for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
      MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
      if (pr->GetMotherId() == -1) {
         nTrMc++;
      }
   }
   if (nTrMc <= 2 * 209) // Just nucleons of Bi+Bi --> Request30-PHSD
   {
      cout << "nTrMc clause - empty event!" << endl;
      cout << " nTrMc " << nTrMc << "; nMC " << nMC << endl;
      return false;
   }

   hEvents->Fill(1.5); // Number of events with filled vertex

   if (!event.fVertex) // if event vertex not filled, skip event
   {
      return false;
   }

   if (mPrimaryVertex.Z() == 0) // not reconstructed (==0)
   {
      return false;
   }

   if (fabs(mPrimaryVertex.Z()) > mZvtxCut) // beyond the limits
   {
      return false;
   }

   hEvents->Fill(2.5); // Number of events after vertex checks (where the vertex was filled, reconstructed and is within
                       // the limits (less than the vertex cut))

   Centrality_tpc = event.getCentrTPC();

   if (Centrality_tpc < 0 || Centrality_tpc >= 100) // TPC centrality not defined
   {
      return false;
   }

   hEvents->Fill(3.5); // Number of events, satisfying all the vertex criteria and having defined centrality

   return true;
}

void MpdGlobalPolarizationRECO::ParticleMCProperties(MpdAnalysisEvent &event)
{
   TVector3 genVert;
   event.fMCEventHeader->GetVertex(genVert);
   nLamb_MC = 0;
   for (int j = 0; j < mMCTracks->GetEntriesFast(); j++) {
      TVector3    mom;
      MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(j);
      mcTr->GetMomentum(mom);
      TVector3 pos;
      double   r = 0.0;
      if (mcTr->GetPdgCode() == pdgCodeMotherHyp) {
         // Check production vertex
         int mpdg = -1;
         if (mcTr->GetMotherId() >= 0) {
            MpdMCTrack *moth = (MpdMCTrack *)mMCTracks->UncheckedAt(mcTr->GetMotherId());
            mpdg             = moth->GetPdgCode();
         }
         mcTr->GetStartVertex(pos);
         pos -= genVert;
         r = pos.Mag();
         if (r < 50.0) {
            hXiFlag->Fill(0);
            double pt = mom.Pt();
            if (mcTr->GetMotherId() == -1) pt *= -1; // negative pT for primaries
            double eta = (TMath::Abs(pt) > 0.001) ? mom.Eta() : TMath::Sign(100., mom.Z());
            fvXiMpdgPtEtaY.push_back(make_tuple(mpdg, pt, eta, mcTr->GetRapidity()));
         }
      } else if (mcTr->GetPdgCode() == pdgCodeHyperon) {
         // Check production vertex

         int mpdg = -1;
         if (mcTr->GetMotherId() >= 0) {
            MpdMCTrack *moth = (MpdMCTrack *)mMCTracks->UncheckedAt(mcTr->GetMotherId());
            mpdg             = moth->GetPdgCode();
         }
         mcTr->GetStartVertex(pos);
         pos -= genVert;
         r = pos.Mag();
         if (r < 50.0) {
            // Production vertex constraint 50 cm
            hLambFlag->Fill(0);
            double pt = mom.Pt();
            if (mcTr->GetMotherId() < 0) pt *= -1; // negative pT for primaries
            double eta = (TMath::Abs(pt) > 0.001) ? mom.Eta() : TMath::Sign(100., mom.Z());
            fvLambMpdgPtEtaY.push_back(make_tuple(mpdg, pt, eta, mcTr->GetRapidity()));
         }
         nLamb_MC++;
         if(mcTr->GetMotherId() < 0){
         hHypMCpTEta->Fill(mcTr->GetRapidity(),mom.Pt());
         if(TMath::Abs(mcTr->GetRapidity())<1. )
         hHypMCpT_Prim->Fill(mom.Pt());
         if(mom.Pt()<3 && mom.Pt()>0.5)
         hHypMCEta_Prim->Fill(mcTr->GetRapidity());
}
      }
   }
}

void MpdGlobalPolarizationRECO::CheckTracks(MpdAnalysisEvent &event)
{
   ids.clear();
   moths.clear();
   pdgs.clear();
   pots.clear();
   ths.clear();
   rads.clear();
   fLays.clear();
   map<int, FairMCPoint *>       points;
   map<int, MpdTpcKalmanTrack *> tracks;

   // Get max. reached layer No.
   for (int j = 0; j < mKalmanTracks->GetEntriesFast(); j++) {
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
      int                id = tr->GetTrackID();
      ids[id]++;
      MpdKalmanHit *hit = (MpdKalmanHit *)tr->GetTrHits()->First();
      if (fLays.find(id) == fLays.end())
         fLays[id] = hit->GetLayer();
      else
         fLays[id] = TMath::Max(hit->GetLayer(), fLays[id]);
   }

   // Exclude "clones" (multiple loops)
   for (int j = 0; j < mKalmanTracks->GetEntriesFast(); j++) {
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
      int                id = tr->GetTrackID();
      if (tracks.find(id) == tracks.end()) tracks[id] = tr;
      // Get track info
      TClonesArray *hits  = tr->GetTrHits();
      int           nHits = hits->GetEntriesFast();

      for (int ih = nHits - 1; ih >= 0; ih--) {
         MpdKalmanHit *hit = (MpdKalmanHit *)hits->UncheckedAt(ih);
         if (hit->GetUniqueID()) continue;
         // No MC points
         if (ids[id] > 1 && points[id]) {
            // More than 1 reco track with the same ID - take the one closer to z = 0
            if (TMath::Abs(tr->GetParam(1)) < TMath::Abs(tracks[id]->GetParam(1))) {
               // Exclude previous track from further consideration
               tracks[id]->SetChi2(-9.);
               tracks[id] = tr;
            } else {
               tr->SetChi2(-9.); // exclude this track from further consideration
               break;
            }
         }
         points[id] = (FairMCPoint *)0x1;
      }      // for (Int_t ih = nHits-1; ih >= 0;
      break; // unclear what this break does here
   }         // for (Int_t j = 0; j < nITS;

   idMax = ids.rbegin()->first;
}
void MpdGlobalPolarizationRECO::CalculateLambdaAcceptance(MpdAnalysisEvent &event)
{
   TVector3 genVert;
   event.fMCEventHeader->GetVertex(genVert);
   multimap<int, int> mapLamb, mapXi;
   for (int j = 0; j <= idMax; j++) {
      TVector3    mom;
      MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(j);
      mcTr->GetMomentum(mom);
      int mothID = mcTr->GetMotherId();
      if (mothID == -1 && fLays.find(j) != fLays.end() && fLays[j] != 0) {
         fLays[j] = -fLays[j]; // flag primary tracks
      }
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      moths[j] = -1;
      pdgs[j]  = mcTr->GetPdgCode();
      if (mothID >= 0) {
         // Check lambda production vertex ( < 50 cm)
         MpdMCTrack *moth = (MpdMCTrack *)mMCTracks->UncheckedAt(mothID);
         moth->GetStartVertex(pos);
         pos -= genVert;
         if (pos.Mag() < 50.0) {
            moths[j] = moth->GetPdgCode();
            if (moths[j] == pdgCodeHyperon && (pdgs[j] == pdgCodeDaughterBar || pdgs[j] == pdgCodeDaughterMes))
               mapLamb.insert(pair<int, int>(mothID, j));
         }
         if (moths[j] == pdgCodeMotherHyp && (pdgs[j] == pdgCodeHyperon || pdgs[j] == pdgCodeDaughterMes))
            mapXi.insert(pair<int, int>(mothID, j));
      }
   }
   multimap<int, int>::iterator                                     mit, mit1;
   pair<multimap<int, int>::iterator, multimap<int, int>::iterator> ret;

   mit = mapLamb.begin();
   while (mit != mapLamb.end()) {
      int mothID = mit->first;
      if (mapLamb.count(mothID) != 2) // only one decay particle
      {
         mit = mapLamb.upper_bound(mothID);
         continue;
      }
      ret         = mapLamb.equal_range(mothID);
      int nppi[2] = {0}, nok = 0;
      int nok1 = 0, nok2 = 0, nok3 = 0;

      for (mit1 = ret.first; mit1 != ret.second; mit1++) {
         TVector3    mom;
         MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(mit1->second);
         if (mcTr->GetPdgCode() == pdgCodeDaughterBar)
            nppi[0] = 1;
         else if (mcTr->GetPdgCode() == pdgCodeDaughterMes)
            nppi[1] = 1;
         mcTr->GetMomentum(mom);
         if (mom.Pt() < 0.001) continue;
         if (TMath::Abs(mom.Eta()) < 1.3) ++nok;
         if ((TMath::Abs(mom.Eta()) < 1.3) && mom.Pt() > 0.05) ++nok1;
         if ((TMath::Abs(mom.Eta()) < 1.3) && mom.Pt() > 0.1) ++nok2;
         if ((TMath::Abs(mom.Eta()) < 1.3) && mom.Pt() > 0.2) ++nok3;
      }
      if (nppi[0] != 1 || nppi[1] != 1) {
         // not p - p- decay
         cout << " Wrong decay mode !!! " << endl;
         mit = mapLamb.upper_bound(mothID);
         continue;
      }
      if (nppi[0] == 1 && nppi[1] == 1) hLambFlag->Fill(1);
      if (nok == 2) hLambFlag->Fill(2);
      if (nok1 == 2) hLambFlag->Fill(4);
      if (nok2 == 2) hLambFlag->Fill(6);
      if (nok3 == 2) hLambFlag->Fill(8);

      // Check Xi-
      MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(mothID);
      int         gmID = mcTr->GetMotherId();

      if (mapXi.find(gmID) != mapXi.end()) {
         ret = mapXi.equal_range(gmID);
         for (mit1 = ret.first; mit1 != ret.second; mit1++) {
            TVector3    mom;
            MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(mit1->second);
            if (mcTr->GetPdgCode() != pdgCodeDaughterMes) continue;
            hXiFlag->Fill(1);
            mcTr->GetMomentum(mom);
            if (mom.Pt() < 0.001) continue;
            if (TMath::Abs(mom.Eta()) < 1.3 && nok == 2) hXiFlag->Fill(2);
            if (TMath::Abs(mom.Eta()) < 1.3 && mom.Pt() > 0.05 && nok1 == 2) hXiFlag->Fill(4);
            if (TMath::Abs(mom.Eta()) < 1.3 && mom.Pt() > 0.1 && nok2 == 2) hXiFlag->Fill(6);
            if (TMath::Abs(mom.Eta()) < 1.3 && mom.Pt() > 0.2 && nok3 == 2) hXiFlag->Fill(8);
         }
      }
      mit = mapLamb.upper_bound(mothID);
   } // while (mit != mapLamb.end())
}

void MpdGlobalPolarizationRECO::fillHistograms(MpdAnalysisEvent &event)
{
   hVertex->Fill(mPrimaryVertex.Z());
   hCentrality->Fill(Centrality_tpc);
   hNevCentr->Fill(Centrality_tpc);

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      if (Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue;

      ResEP1_true[iter_cent] += ResEP;
      SubEvRes1[iter_cent] += ResEPSub;
   }

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      hResolution_EP1_true->SetBinContent(iter_cent + 1, ResEP1_true[iter_cent]);
      hResolution_EP1_reco->SetBinContent(iter_cent + 1, SubEvRes1[iter_cent]);
   }
   if (analysis_choice == "selection") {
      for (int i = 0; i < vLambdas.size(); i++) // cycle for reco Lambda
      {
         if (selection_choice == "omega2") {
            MpdLambdaPol *lamb = (MpdLambdaPol *)&vLambdas.at(i);
            if (lamb->origs[0] > -8) hm0_before_full->Fill(lamb->massh);
            for (int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++) {
               if (lamb->omega2 > omega_value[iter_sel]) {
                  if (lamb->origs[0] < -8) hm0_full_mixed[iter_sel]->Fill(lamb->massh);
                  if (lamb->origs[0] > -8) hm0_full[iter_sel]->Fill(lamb->massh);
               }
            }

         } else if (selection_choice == "chi") {
            // so far nothing, as we are using the tree and post-analysis to find selection values
         }
      }

   } else if (analysis_choice == "analysis") {
      for (int i = 0; i < vLambdas.size(); i++) // cycle for reco Lambda
      {
         MpdLambdaPol *lamb = (MpdLambdaPol *)&vLambdas.at(i);
         if (lamb->origs[0] > -8) {
            hm0_before_full->Fill(lamb->massh);
            if (selection_choice == "omega2") {
               if (lamb->omega2 > omega_value_full){ hm0_Full->Fill(lamb->massh); 
if(lamb->origs[0] ==1){
if(TMath::Abs(lamb->etah)<1.)
hHypRecopT_Prim->Fill(lamb->pth);
if(lamb->pth>0.5 && lamb->pth<3)
hHypRecoEta_Prim->Fill(lamb->etah);
hHypRecopTEta_Prim->Fill(lamb->etah,lamb->pth);}
if(lamb->origs[0]>0){
hHypRecopTEta_Full->Fill(lamb->etah,lamb->pth);
if(TMath::Abs(lamb->etah)<1.)
hHypRecopT_Full->Fill(lamb->pth);
if(lamb->pth>0.5 && lamb->pth<3)
hHypRecoEta_Full->Fill(lamb->etah);}
           /* if (Centrality_tpc >= centrality_max[2] || Centrality_tpc < centrality_min[2]) continue;
            for (int iter_pt = 0; iter_pt < NITER_PT; iter_pt++) {
               for (int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++) {
                  if (lamb->etah >= eta_max[iter_eta] || lamb->etah < eta_min[iter_eta]) continue;
                  if (lamb->pth >= pt_max[iter_pt] || lamb->pth < pt_min[iter_pt]) continue;
                   double phi_diff_hist   = phiEP - lamb->phi_star;
               double phi_diff_histRP = phiRP - lamb->phi_star;
               if (phi_diff_hist < 0) phi_diff_hist = phi_diff_hist + 2. * pi;
               if (phi_diff_histRP < 0) phi_diff_histRP = phi_diff_histRP + 2. * pi;
                  
                   hm0_pt_eta_bin[iter_pt][iter_eta]->Fill(lamb->massh);
                   hPolarY_Full_pt_eta_bin[iter_pt][iter_eta]->Fill(lamb->polarhy);
                   hDeltaPhiEP_Full_pt_eta_bin[iter_pt][iter_eta]->Fill(phi_diff_hist);
                   hDeltaPhiRP_Full_pt_eta_bin[iter_pt][iter_eta]->Fill(phi_diff_histRP);
                   hv1EP_Full_pt_eta_bin[iter_pt][iter_eta]->Fill(lamb->massh,TMath::Cos(phi_Lam - phiEP));
                   hv1RP_Full_pt_eta_bin[iter_pt][iter_eta]->Fill(lamb->massh,TMath::Cos(phi_Lam - phiRP));
                   if(lamb->origs[0] == 1) // true lambda
                   {hPolarY_Prim_pt_eta_bin[iter_pt][iter_eta]->Fill(lamb->polarhy);
                   hDeltaPhiEP_Prim_pt_eta_bin[iter_pt][iter_eta]->Fill(phi_diff_hist);
                   hDeltaPhiRP_Prim_pt_eta_bin[iter_pt][iter_eta]->Fill(phi_diff_histRP);
                   hv1EP_Prim_pt_eta_bin[iter_pt][iter_eta]->Fill(lamb->massh,TMath::Cos(phi_Lam - phiEP));
                   hv1RP_Prim_pt_eta_bin[iter_pt][iter_eta]->Fill(lamb->massh,TMath::Cos(phi_Lam - phiRP));}
                  }
                }*/
               
}            


            } else if (selection_choice == "chi") {
               if ((lamb->chi2s[0] > chi_pi_value_full) && (lamb->chi2s[1] > chi_p_value_full) &&
                   (lamb->chi2h < chi_V0_value_full) && (lamb->path > lambda_path_value_full) &&
                   (lamb->angle < lambda_angle_value_full))
                  hm0_Full->Fill(lamb->massh);
            }
         }
         for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
            if (Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue;

            if (lamb->origs[0] > 0) {
               hPolarY_Full[iter_cent]->Fill(lamb->polarhy);

               double phi_diff_hist   = phiEP - lamb->phi_star;
               double phi_diff_histRP = phiRP - lamb->phi_star;
               double phi_diff_MC     = phiRP - lamb->phi_star_MC;
               if (phi_diff_hist < 0) phi_diff_hist = phi_diff_hist + 2. * pi;
               if (phi_diff_histRP < 0) phi_diff_histRP = phi_diff_histRP + 2. * pi;
               if (phi_diff_MC < 0) phi_diff_MC = phi_diff_MC + 2. * pi;

               hDeltaPhiEP_Full[iter_cent]->Fill(phi_diff_hist);
               hDeltaPhiRP_Full[iter_cent]->Fill(phi_diff_histRP);
               hDeltaPhiRP_MC_Full[iter_cent]->Fill(phi_diff_MC);

               hPolvsPt[iter_cent]->Fill(lamb->pth, lamb->polarhy);
               hPolvsEta[iter_cent]->Fill(lamb->etah, lamb->polarhy);
            }
            if (lamb->origs[0] == 1) // true lambda
            {
               hPolarY_Prim[iter_cent]->Fill(lamb->polarhy);

               float phi_diff_hist   = phiEP - lamb->phi_star;
               float phi_diff_histRP = phiRP - lamb->phi_star;
               float phi_diff_MC     = phiRP - lamb->phi_star_MC;
               if (phi_diff_hist < 0) phi_diff_hist = phi_diff_hist + 2. * pi;
               if (phi_diff_histRP < 0) phi_diff_histRP = phi_diff_histRP + 2. * pi;
               if (phi_diff_MC < 0) phi_diff_MC = phi_diff_MC + 2. * pi;

               hDeltaPhiEP_Prim[iter_cent]->Fill(phi_diff_hist);
               hDeltaPhiRP_Prim[iter_cent]->Fill(phi_diff_histRP);
               hDeltaPhiRP_MC_Prim[iter_cent]->Fill(phi_diff_MC);
            }
            if (lamb->origs[0] > -8) hm0_before[iter_cent]->Fill(lamb->massh);

            if (selection_choice == "omega2") {
               if (lamb->omega2 <= omega_value_full)
                  continue; // now using only the value for full dataset (MinBias), not for each centrality bin
            } else if (selection_choice == "chi") {
               if (lamb->chi2s[0] <= chi_pi_value_full) continue;
               if (lamb->chi2s[1] <= chi_p_value_full) continue;
               if (lamb->chi2h >= chi_V0_value_full) continue;
               if (lamb->path <= lambda_path_value_full) continue;
               if (lamb->angle >= lambda_angle_value_full) continue;
            }

            if (lamb->origs[0] > -8) hm0_after[iter_cent]->Fill(lamb->massh);
            double phi_diff = phiEP - lamb->phi_star;
            if (phi_diff < 0) phi_diff = phi_diff + 2. * pi;

            for (int iter = 0; iter < NITER; iter++) {
               if (phi_diff < angle_max[iter] && phi_diff >= angle_min[iter]) {
                  if (lamb->origs[0] > -8) hm0[iter_cent][iter]->Fill(lamb->massh);
                  if (lamb->origs[0] < -8) hm0_mixed[iter_cent][iter]->Fill(lamb->massh);
               }
            }

            // here we collect the histograms for pT and eta dependence, for centrality region 20-50%:
            if (NITER_CENT != 4) continue;
            if (iter_cent != 2) continue;
            for (int iter_pt = 0; iter_pt < NITER_PT; iter_pt++) {
               if (lamb->pth >= pt_max[iter_pt] || lamb->pth < pt_min[iter_pt]) continue;
               for (int iter = 0; iter < NITER; iter++) {
                  if (phi_diff < angle_max[iter] && phi_diff >= angle_min[iter]) {
                     if (lamb->origs[0] > -8) hm0_ptbin[iter_pt][iter]->Fill(lamb->massh);
                     if (lamb->origs[0] < -8) hm0_ptbin_mixed[iter_pt][iter]->Fill(lamb->massh);
                  }
               }
            }
            for (int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++) {
               if (lamb->etah >= eta_max[iter_eta] || lamb->etah < eta_min[iter_eta]) continue;
               for (int iter = 0; iter < NITER; iter++) {
                  if (phi_diff < angle_max[iter] && phi_diff >= angle_min[iter]) {
                     if (lamb->origs[0] > -8) hm0_etabin[iter_eta][iter]->Fill(lamb->massh);
                     if (lamb->origs[0] < -8) hm0_etabin_mixed[iter_eta][iter]->Fill(lamb->massh);
                  }
               }
            }
         }
      }
   }
}

void MpdGlobalPolarizationRECO::CollectParticles(vector<int> &vecPi, vector<int> &vecK, vector<int> &vecP)
{
   for (Int_t j = 0; j < mKalmanTracks->GetEntriesFast(); j++) {
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      int         id   = tr->GetTrackID();
      MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(id);
      // !!!
      if (mcTr->GetMotherId() == 0 &&
          ((MpdMCTrack *)mMCTracks->UncheckedAt(0))->GetPdgCode() == 1010010030) // !!! decay product of artificial H3L
         continue;
      // !!!
      if (mcTr->GetPdgCode() == pdgCodeDaughterBar &&
          tr->Charge() == TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeDaughterBar)->Charge() / 3))
         vecP.push_back(j);
      else if (mcTr->GetPdgCode() == pdgCodeDaughterMes &&
               tr->Charge() == TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeDaughterMes)->Charge() / 3))
         vecPi.push_back(j);
   }
}
MpdHelix MpdGlobalPolarizationRECO::MakeHelix(const MpdTpcKalmanTrack *tr)
{
   double r   = tr->GetPosNew();
   double phi = tr->GetParam(0) / r;
   double x   = r * TMath::Cos(phi);
   double y   = r * TMath::Sin(phi);
   double dip = tr->GetParam(3);
   double cur = 0.3 * 0.01 * 5 / 10; // 5 kG
   cur *= TMath::Abs(tr->GetParam(4));
   TVector3 o(x, y, tr->GetParam(1));
   int      h = (int)TMath::Sign(1.1, tr->GetParam(4));
   MpdHelix helix(cur, dip, tr->GetParam(2) - TMath::PiOver2() * h, o, h);
   return helix;
}

MpdHelix MpdGlobalPolarizationRECO::MakeHelix(const MpdParticle *part)
{
   double dip = TMath::PiOver2() - part->Theta();
   double cur = TMath::Abs(part->GetMeas(4));
   if (part->GetCharge() == 0) cur = numeric_limits<double>::epsilon();
   int      h     = (int)TMath::Sign(1.1, part->GetMeas(4));
   double   phase = part->GetMeas(2) - TMath::PiOver2() * h;
   double   x     = part->GetXY(0);
   double   y     = part->GetXY(1);
   TVector3 o(x, y, part->GetMeas(1));
   MpdHelix helix(cur, dip, phase, o, h);
   return helix;
}

void MpdGlobalPolarizationRECO::ApplyPid(vector<int> &vecP, vector<int> &vecPi)
{
   // Fill the maps for the Lambda and Xi flags:
   map<int, set<int>> mapL, mapXi;
   map<int, set<int>> mapL13, mapXi13;

   int nP  = vecP.size();
   int nPi = vecPi.size();

   for (int ip = 0; ip < nP; ip++) {
      MpdTpcKalmanTrack *trP = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecP[ip]);
      if (trP->GetUniqueID() > 0) {
         // Lambda decay product
         int mid = trP->GetUniqueID();
         if (mapL.find(mid) == mapL.end()) {
            set<int> aaa;
            mapL[mid] = aaa;
         }
         mapL[mid].insert(vecP[ip]);
         if (TMath::Abs(trP->Momentum3().Eta()) < 1.3) {
            if (mapL13.find(mid) == mapL13.end()) {
               set<int> aaa;
               mapL13[mid] = aaa;
            }
            mapL13[mid].insert(vecP[ip]);
         }
         trP->SetUniqueID(0); // reset
      }
      if (trP->GetVertex().GetUniqueID() > 0) {
         // Xi- decay product
         int mid = trP->GetVertex().GetUniqueID();
         if (mapXi.find(mid) == mapXi.end()) {
            set<int> aaa;
            mapXi[mid] = aaa;
         }
         mapXi[mid].insert(vecP[ip]);
      }
   }

   for (int ip = 0; ip < nPi; ip++) {
      MpdTpcKalmanTrack *trP = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecPi[ip]);
      if (trP->GetUniqueID() > 0) {
         // Lambda decay product
         int mid = trP->GetUniqueID();
         if (mapL.find(mid) == mapL.end()) {
            set<int> aaa;
            mapL[mid] = aaa;
         }
         mapL[mid].insert(vecPi[ip]);
         if (TMath::Abs(trP->Momentum3().Eta()) < 1.3) {
            if (mapL13.find(mid) == mapL13.end()) {
               set<int> aaa;
               mapL13[mid] = aaa;
            }
            mapL13[mid].insert(vecPi[ip]);
         }
         trP->SetUniqueID(0); // reset
      }
      if (trP->GetVertex().GetUniqueID() > 0) {
         // Xi- decay product
         int mid = trP->GetVertex().GetUniqueID();
         if (mapXi.find(mid) == mapXi.end()) {
            set<int> aaa;
            mapXi[mid] = aaa;
         }
         mapXi[mid].insert(vecPi[ip]);
      }
   }

   for (map<int, set<int>>::iterator mit = mapL.begin(); mit != mapL.end(); mit++)
      if (mit->second.size() != 2) mit->second.insert(-999); // not 2 decay products reconstructed - add fake indx

   for (map<int, set<int>>::iterator mit = mapL13.begin(); mit != mapL13.end(); mit++)
      if (mit->second.size() != 2) mit->second.insert(-999); // not 2 decay products reconstructed - add fake indx

   // Get TOF matches
   int           nTofMatch = fTofMatches->GetEntriesFast();
   map<int, int> mapTof;

   for (int itof = 0; itof < nTofMatch; itof++) {
      MpdTofMatchingData *match        = (MpdTofMatchingData *)fTofMatches->UncheckedAt(itof);
      mapTof[match->GetKFTrackIndex()] = itof;
   }

   vecP.clear();
   vecPi.clear();

   // Refill the veP and vecPi vectors using the information from PID:
   int nTracks = mKalmanTracks->GetEntriesFast();
   for (int j = 0; j < nTracks; j++) {
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
      // Cut out bad tracks
      if (tr->GetChi2() < -8) continue;
      int         id     = tr->GetTrackID();
      MpdMCTrack *mcTr   = (MpdMCTrack *)mMCTracks->UncheckedAt(id);
      int         mothId = mcTr->GetMotherId();
      int         uid    = tr->GetVertex().GetUniqueID();
      int         ret = 0, eta13 = 0, charge = tr->Charge();
      bool        m2Flag = kFALSE; // flag to check whether there is m2 from TOF
      double      trEta  = TMath::Abs(tr->Momentum3().Eta());
      double      dedx   = tr->GetDedx();
      double      m2     = -1;

      // Check whether there is m2 from TOF
      if (mapTof.count(j) > 0) {
         MpdTofMatchingData *match = (MpdTofMatchingData *)fTofMatches->UncheckedAt(mapTof[j]);
         m2                        = match->GetMass2();
         m2Flag                    = kTRUE;
      }

      // Change flag if eta of the track is less than 1.3
      if (TMath::Abs(tr->Momentum3().Eta()) < 1.3) eta13 = 1;
      // dE/dx+TOF (pid if we have both dedx and m2)
      if (m2Flag) ret = fPid->FillProbs(tr->Momentum(), dedx, m2, charge);
      // only dE/dx (pid if we have only dedx)
      if (ret == 0) ret = fPid->FillProbs(tr->Momentum(), dedx, charge);
      //!!! No PID for (anti)protons above 2.5 GeV/c !!!
      if (tr->Momentum() > 2.5 && charge * pdgCodeDaughterBar > 0) ret = 1;
      // if there is no PID, but the MC track exists, fill Flags
      if (ret == 0 && eta13) {
         if (mcTr->GetPdgCode() == pdgCodeDaughterMes) hPIDflag->Fill(2.1); // lost pion
         if (mcTr->GetPdgCode() == pdgCodeDaughterBar) hPIDflag->Fill(6.1); // lost proton
      }

      // negative threshold values for extra check (the value of prob should not be less than 0)
      double piThr   = -0.75;
      double probThr = -0.60;

      // looking for pi- (in case of Lambda) or pi+ (in case of ALambda)
      if (pdgCodeHyperon * charge < 0) {
         double prob = fPid->GetProbPi(); // for nsig method can be either 1 or 0 -> if 1 we found a particle
         if (prob > piThr && prob > fPid->GetProbKa() && prob > fPid->GetProbPr() && prob > fPid->GetProbDe() &&
             prob > fPid->GetProbTr() && prob > fPid->GetProbHe3() && prob > fPid->GetProbHe4()) {
            // "pion"
            if (mcTr->GetPdgCode() == pdgCodeDaughterMes && eta13)
               hPIDflag->Fill(0.1); // correct pion
            else if (mcTr->GetPdgCode() != pdgCodeDaughterMes && eta13)
               hPIDflag->Fill(1.1); // false pion
            //
            if (mapL.find(mothId + 1) != mapL.end() && mapL[mothId + 1].find(j) != mapL[mothId + 1].end())
               mapL[mothId + 1].erase(j);
            if (mapL13.find(mothId + 1) != mapL13.end() && mapL13[mothId + 1].find(j) != mapL13[mothId + 1].end())
               mapL13[mothId + 1].erase(j);
            if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end()) mapXi[uid].erase(j);
            //
            double chi2 = TMath::Min(tr->GetChi2Vertex(), 999.);
            if (chi2 < gC2pi) // less than 5
               continue;
            vecPi.push_back(j);
         } else if (mcTr->GetPdgCode() == pdgCodeDaughterMes && eta13)
            hPIDflag->Fill(2.1); // lost pion
      } else                     // now looking for proton (in case of Lambda) or anti-proton (in case of ALambda)
      {
         if (trEta < 0.5 && mcTr->GetPdgCode() == pdgCodeDaughterBar) hPtProt->Fill(tr->Pt());
         double prob = fPid->GetProbPr();
         if (tr->Momentum() > 2.5 && charge * pdgCodeDaughterBar > 0) prob = 9.9; //!!! force (anti)proton !!! AZ-130423
         if (prob > probThr && prob > fPid->GetProbKa() && prob > fPid->GetProbPi() && prob > fPid->GetProbDe()) {
            // "proton"
            if (mcTr->GetPdgCode() == pdgCodeDaughterBar) {
               if (eta13) hPIDflag->Fill(4.1); // correct proton
               if (trEta < 0.5) hPtProtT->Fill(tr->Pt());
            } else if (mcTr->GetPdgCode() != pdgCodeDaughterBar) {
               if (eta13) hPIDflag->Fill(5.1); // false proton
               if (trEta < 0.5) hPtProtF->Fill(tr->Pt());
            }
            //
            if (mapL.find(mothId + 1) != mapL.end() && mapL[mothId + 1].find(j) != mapL[mothId + 1].end())
               mapL[mothId + 1].erase(j);
            if (mapL13.find(mothId + 1) != mapL13.end() && mapL13[mothId + 1].find(j) != mapL13[mothId + 1].end())
               mapL13[mothId + 1].erase(j);
            if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end()) mapXi[uid].erase(j);
            // Refit for proton track:
            MpdTpcKalmanTrack trCor = *tr;
            trCor.SetDirection(MpdKalmanTrack::kInward);
            int ok = 0;
            ok     = recoTpc->Refit(&trCor, massPr, 1); // refit
            if (!ok) continue;
            MpdParticle prot(trCor, 0);
            prot.SetPdg(pdgCodeDaughterBar);
            prot.SetMass(); // this is not necessary, SetPdg should call SetMass
            double chi2 = TMath::Min(prot.Chi2Vertex(fMpdVert), 999.);
            if (chi2 < gC2p) // more than 3
               continue;
            vecP.push_back(j);
         } else if (mcTr->GetPdgCode() == pdgCodeDaughterBar && eta13)
            hPIDflag->Fill(6.1); // lost proton
      }
   }

   // Fill flags if there are still particles in maps, which were not identified by PID
   int nLok    = 0;
   int nXiok   = 0;
   int nLok13  = 0;
   int nXiok13 = 0;
   for (map<int, set<int>>::iterator mit = mapL.begin(); mit != mapL.end(); mit++) {
      if (mit->second.size() == 0) nLok++;
   }
   for (map<int, set<int>>::iterator mit = mapL13.begin(); mit != mapL13.end(); mit++) {
      if (mit->second.size() == 0) nLok13++;
   }
   for (map<int, set<int>>::iterator mit = mapXi.begin(); mit != mapXi.end(); mit++) {
      if (mit->second.size() == 0) nXiok++;
   }

   hLambFlag->Fill(13, nLok);
   hLambFlag->Fill(11, nLok13);
   hXiFlag->Fill(11, nXiok);
}

void MpdGlobalPolarizationRECO::RecoEff(vector<int> &vecP, vector<int> &vecPi, bool use_pid)
{
   int nPi = vecPi.size();
   int nP  = vecP.size();

   for (int ip = nP - 1; ip >= 0; ip--) // Proton
   {
      MpdTpcKalmanTrack *trP    = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecP[ip]);
      MpdMCTrack        *mcTr   = (MpdMCTrack *)mMCTracks->UncheckedAt(trP->GetTrackID());
      int                mothId = mcTr->GetMotherId();
      if (mothId < 0) continue;
      MpdMCTrack *moth = (MpdMCTrack *)mMCTracks->UncheckedAt(mothId);
      if (moth->GetPdgCode() == pdgCodeHyperon) // Proton from Lambda
      {
         Int_t mp = mothId;
         for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) // Pion
         {
            MpdTpcKalmanTrack *trPi   = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecPi[jpi]);
            MpdMCTrack        *mcTr   = (MpdMCTrack *)mMCTracks->UncheckedAt(trPi->GetTrackID());
            int                mothId = mcTr->GetMotherId();
            if (mothId < 0) continue;
            MpdMCTrack *moth = (MpdMCTrack *)mMCTracks->UncheckedAt(mothId);
            if (moth->GetPdgCode() == pdgCodeHyperon && mp == mothId) {
               hLambFlag->Fill(12);
               if (TMath::Abs(trP->Momentum3().Eta()) < 1.3 && TMath::Abs(trPi->Momentum3().Eta()) < 1.3)
                  hLambFlag->Fill(10);
               // AZ - flag decay tracks to check PID influence later
               trP->SetUniqueID(mothId + 1);
               trPi->SetUniqueID(mothId + 1);
               int gmId = moth->GetMotherId();
               if (gmId >= 0) {
                  MpdMCTrack *gmoth = (MpdMCTrack *)mMCTracks->UncheckedAt(gmId);
                  if (gmoth->GetPdgCode() == pdgCodeMotherHyp) {
                     for (int kpi = nPi - 1; kpi >= 0; kpi--) // Pion
                     {
                        MpdTpcKalmanTrack *trK    = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecPi[kpi]);
                        MpdMCTrack        *mcTr   = (MpdMCTrack *)mMCTracks->UncheckedAt(trK->GetTrackID());
                        int                mothId = mcTr->GetMotherId();
                        if (mothId < 0) continue;
                        MpdMCTrack *moth = (MpdMCTrack *)mMCTracks->UncheckedAt(mothId);
                        if (moth->GetPdgCode() == pdgCodeMotherHyp && gmId == mothId) {
                           hXiFlag->Fill(10);
                           // AZ - flag decay tracks to check PID influence later
                           trK->GetVertex().SetUniqueID(mothId + 1);
                           trP->GetVertex().SetUniqueID(mothId + 1);
                           trPi->GetVertex().SetUniqueID(mothId + 1);
                           break;
                        }
                     }
                  } // if (gmoth->GetPdgCode() == pdgCodeMotherHyp)
               }
               break;
            }
         }
      }
   }

   if (use_pid) return; // skip the rest if PID is used

   for (int ip = nP - 1; ip >= 0; ip--) // Proton
   {
      MpdTpcKalmanTrack *trP   = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecP[ip]);
      MpdTpcKalmanTrack  trCor = *trP;
      trCor.SetDirection(MpdKalmanTrack::kInward);
      recoTpc->Refit(&trCor, massPr, 1); // refit
      MpdParticle prot(trCor, vecP[ip]);
      prot.SetPdg(pdgCodeDaughterBar);
      prot.SetMass(); // this is not necessary, SetPdg should call SetMass

      double chi2 = TMath::Min(prot.Chi2Vertex(fMpdVert), 999.);
      if (chi2 < gC2p) vecP.erase(vecP.begin() + ip);
   }

   if (nP) {
      for (Int_t jpi = nPi - 1; jpi >= 0; jpi--) // Pion
      {
         MpdTpcKalmanTrack *trPi  = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecPi[jpi]);
         MpdTpcKalmanTrack  trCor = *trPi;

         double chi2 = TMath::Min(trPi->GetChi2Vertex(), 999.);
         if (chi2 < gC2pi) vecPi.erase(vecPi.begin() + jpi);
      }
   }
}

void MpdGlobalPolarizationRECO::BuildLambda(vector<int> &vecP, vector<int> &vecPi, vector<MpdParticle *> &vecL,
                                            double &phiEP)
{
   qs[0] = 0, qs[1] = 0, origs[0] = 0, origs[1] = 0, layMx[0] = 0, layMx[0] = 0, dstNo[0] = 0, dstNo[1] = 0;
   etas[0] = 0, etas[1] = 0, ps[0] = 0, ps[1] = 0, pts[0] = 0, pts[1] = 0, chi2s[0] = 0, chi2s[1] = 0, dcas[0] = 0,
   dcas[1] = 0, c2s[0] = 0, c2s[1] = 0;
   mcps[0] = 0, mcps[1] = 0, mcthetas[0] = 0, mcthetas[1] = 0, thetas[0] = 0, thetas[1] = 0, mcphis[0] = 0,
   mcphis[0] = 0, phis[0] = 0, phis[1] = 0;
   path = 0, massh = 0, chi2h = 0, angle = 0, pth = 0, ph = 0, etah = 0, disth = 0, yh = 0;
   dca = 0, omega1 = 0, omega2 = 0, cosA = 0, cosAmc = 0, polarhx = 0, polarhy = 0, polarhz = 0;
   phi_star = 0, phi_star_MC = 0, phi_Lam = 0;

   int                   nPi     = vecPi.size();
   int                   nP      = vecP.size();
   int                   saveMix = 0;
   int                   mpdg    = 0;
   vector<MpdParticle *> vPart;
   vecL1.clear();
   fVecL1.clear();
   fVecL2.clear();

   // Cycle through identified protons
   for (int ip = 0; ip < nP; ip++) {
      MpdTpcKalmanTrack *trP  = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecP[ip]);
      MpdMCTrack        *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(trP->GetTrackID());
      TVector3           mcmom1;
      mcTr->GetMomentum(mcmom1);
      mcps[1]                  = mcmom1.Mag();
      mcphis[1]                = mcmom1.Phi();
      mcthetas[1]              = mcmom1.Theta();
      int               mothId = mcTr->GetMotherId();
      MpdTpcKalmanTrack trCor  = *trP;
      trCor.SetDirection(MpdKalmanTrack::kInward);
      recoTpc->Refit(&trCor, massPr, 1); // refit
      MpdParticle prot(trCor, vecP[ip]);
      prot.SetPdg(pdgCodeDaughterBar);
      prot.SetMass(); // this is not necessary, SetPdg should call SetMass
      qs[1]          = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeDaughterBar)->Charge() / 3);
      etas[1]        = trP->Momentum3().Eta();
      ps[1]          = trP->Momentum();
      pts[1]         = trP->Pt();
      chi2s[1]       = TMath::Min(prot.Chi2Vertex(fMpdVert), 9999.);
      c2s[1]         = trP->GetChi2() / (trP->GetNofTrHits() * 2 - 5);
      layMx[1]       = TMath::Abs(fLays[trP->GetTrackID()]);
      MpdHelix helix = MakeHelix(trP);
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      double   s = helix.pathLength(mPrimaryVertex);
      pca        = helix.at(s);
      pca -= mPrimaryVertex;
      dcas[1]  = pca.Mag();
      origs[1] = 0;
      if (mothId >= 0 && ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetPdgCode() == pdgCodeHyperon)
         origs[1] = -1; // from lambda
      saveMix++;

      for (int jpi = 0; jpi < nPi; jpi++) // Pion
      {
         MpdTpcKalmanTrack *trPi  = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(vecPi[jpi]);
         MpdMCTrack        *mcTr1 = (MpdMCTrack *)mMCTracks->UncheckedAt(trPi->GetTrackID());
         TVector3           mcmom2;
         mcTr1->GetMomentum(mcmom2);
         mcps[0]     = mcmom2.Mag();
         mcphis[0]   = mcmom2.Phi();
         mcthetas[0] = mcmom2.Theta();
         int mothId1 = mcTr1->GetMotherId();
         origs[0]    = 0;
         if (mothId1 >= 0 && ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId1))->GetPdgCode() == pdgCodeHyperon)
            origs[0] = -1; // from lambda
         MpdParticle *pion = new MpdParticle(*trPi, vecPi[jpi]);
         pion->SetPdg(pdgCodeDaughterMes);
         pion->SetMass(); // this is not necessary, SetPdg should call SetMass
         if (nMix > 0 && saveMix == 1) fMapPiEvent.insert(pair<int, MpdTpcKalmanTrack>(fEvNo, *trPi));
         vPart.clear();
         vPart.push_back(new MpdParticle(prot));
         vPart.push_back(pion);
         MpdParticle lambPart;
         double      chi2 = lambPart.BuildMother(vPart);
         TVector3    v0(lambPart.Getx()(0, 0), lambPart.Getx()(1, 0), lambPart.Getx()(2, 0));
         v0 -= mPrimaryVertex;
         double decay = v0.Mag();
         path         = TMath::Sign(decay, v0 * lambPart.Momentum3());
         massh        = lambPart.GetMass();

         if (chi2 >= 0 && chi2 < gC2L && path > gPathL && massh < 1.2) {
            if (origs[1] > 0) origs[1] = -1;
            MpdMCTrack *moth = NULL;
            hMassL->Fill(lambPart.GetMass());
            if (mothId != mothId1 || mothId < 0) {
               hMassLbkg->Fill(lambPart.GetMass());
            } else {
               if (origs[0] == -1) {
                  hMassLsig->Fill(lambPart.GetMass());
                  origs[0] = origs[1] = 1;
                  moth                = (MpdMCTrack *)mMCTracks->UncheckedAt(mothId);
               } else
                  hMassLbkg->Fill(lambPart.GetMass());
            }

            qs[0]    = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeDaughterMes)->Charge() / 3); // pion
            etas[0]  = trPi->Momentum3().Eta();
            ps[0]    = trPi->Momentum();
            pts[0]   = trPi->Pt();
            chi2s[0] = TMath::Min(pion->Chi2Vertex(fMpdVert), 9999.);
            c2s[0]   = trPi->GetChi2() / (trPi->GetNofTrHits() * 2 - 5);
            layMx[0] = TMath::Abs(fLays[trPi->GetTrackID()]);
            MpdHelix helix1 = MakeHelix(trPi);
            // Get 3-D DCA to primary vertex
            s   = helix1.pathLength(mPrimaryVertex);
            pca = helix1.at(s);
            pca -= mPrimaryVertex;
            dcas[0] = pca.Mag();

            chi2h   = chi2;
            angle   = v0.Angle(lambPart.Momentum3());
            pth     = lambPart.Pt();       // reconstructed
            ph      = lambPart.Momentum(); // reconstructed
            phi_Lam = lambPart.Phi();      // reconstructed
            if (pth > 0.001)
               etah = lambPart.Momentum3().Eta();
            else
               etah = TMath::Sign(100., lambPart.Momentum3().Z());
            pair<double, double> paths = helix.pathLengths(helix1);
            TVector3             p1    = helix.at(paths.first);
            TVector3             p2    = helix1.at(paths.second);
            p1 -= p2;
            disth = p1.Mag(); // closest distance between daughters

            // Get 3-D DCA of lambda to primary vertex
            MpdHelix helix2 = MakeHelix(&lambPart);
            s               = helix2.pathLength(mPrimaryVertex);
            pca             = helix2.at(s);
            pca -= mPrimaryVertex;
            dca    = pca.Mag();
            c2pv   = TMath::Min(lambPart.Chi2Vertex(fMpdVert), 9999.);
            omega1 = dcas[0] * dcas[1] / (dca * dca + disth * disth);
            omega2 = TMath::Sqrt(chi2s[0] * chi2s[1]) / (c2pv + chi2h);

            dstNo[0] = vecPi[jpi]; // pion index
            dstNo[1] = vecP[ip];   // proton index

            if (lambPart.GetMass() >= 1.10518 && lambPart.GetMass() <= 1.12668) // lambda mass +- 5*2.15 MeV
            {
               vecL.push_back(new MpdParticle(lambPart));
               vector<double> lambPars(6);
               lambPars[0] = disth;
               lambPars[1] = angle;
               for (Int_t jl = 0; jl < 2; jl++) {
                  lambPars[jl + 2] = chi2s[jl];
                  lambPars[jl + 4] = dcas[jl];
               }
               vecL1.push_back(lambPars);
               fVecL1.push_back(pair<double, double>(disth, angle));
               fVecL2.push_back(pair<double, double>(chi2s[0], chi2s[1]));
            }

            lambPart.SetMass(massHyperon); // set true mass
            yh   = lambPart.Rapidity();
            mpdg = -1;
            // Check mother of lambda
            if (origs[0] == 1) {
               int gMothId = moth->GetMotherId();
               if (gMothId >= 0) {
                  origs[0] = origs[1] = 2; // secondary lambda
                  mpdg                = ((MpdMCTrack *)mMCTracks->UncheckedAt(gMothId))->GetPdgCode();
               }
            }

            polarhx = 0.0;
            polarhy = 0.0;
            polarhz = 0.0;
            if (origs[0] > 0) {
               double weight_pol = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetWeight();
               polarhx           = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetPolar(0);
               polarhy           = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetPolar(1);
               polarhz           = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetPolar(2);
               TVector3 polar_changed(polarhx, polarhy, polarhz);
               if (phiEP != 0.) polar_changed.RotateZ(-phiEP);
               polarhx = weight_pol * polar_changed.X();
               polarhy = weight_pol * polar_changed.Y();
               polarhz = weight_pol * polar_changed.Z();
            }

            FindPolarAngle(lambPart, vPart);
            MpdLambdaPol l0(massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas, mcthetas, thetas, mcphis, phis,
                            mcps, ps, pts, chi2s, dcas, dca, c2pv, omega1, omega2, cosA, cosAmc, polarhx, polarhy,
                            polarhz, phi_star, phi_star_MC, phi_Lam, c2s, origs, qs, layMx, mpdg);
            vLambdas.push_back(l0);
         } // if (chi2 >= 0 && chi2 < gC2L...

         int nPart = vPart.size();
         for (int ipart = 0; ipart < nPart; ++ipart) delete vPart[ipart];
      } // for (int jpi = 0; jpi < nPi;
      if (nMix == 0) continue;
      // Event mixing
      for (map<int, MpdVertex>::iterator mit = fMapVertexEvent.begin(); mit != fMapVertexEvent.end(); mit++) {
         if (mit->first == fEvNo) break;
         double dz = mit->second.GetZ() - fMapVertexEvent.rbegin()->second.GetZ();
         pair<multimap<int, MpdTpcKalmanTrack>::iterator, multimap<int, MpdTpcKalmanTrack>::iterator> piter =
            fMapPiEvent.equal_range(mit->first);

         for (multimap<int, MpdTpcKalmanTrack>::iterator mit1 = piter.first; mit1 != piter.second; mit1++) {
            MpdTpcKalmanTrack piTr = mit1->second;
            double            z0tr = piTr.GetParam(1);
            piTr.SetParam(1, z0tr - dz);
            origs[0]          = -9;
            MpdParticle *pion = new MpdParticle(piTr);
            pion->SetPdg(pdgCodeDaughterMes);
            pion->SetMass(); // this is not necessary, SetPdg should call SetMass

            vPart.clear();
            vPart.push_back(new MpdParticle(prot));
            vPart.push_back(pion);

            MpdParticle lambPart;
            double      chi2 = lambPart.BuildMother(vPart);
            TVector3    v0(lambPart.Getx()(0, 0), lambPart.Getx()(1, 0), lambPart.Getx()(2, 0));
            v0 -= mPrimaryVertex;
            double decay = v0.Mag();
            path         = TMath::Sign(decay, v0 * lambPart.Momentum3());
            massh        = lambPart.GetMass();

            if (chi2 >= 0 && chi2 < gC2L && path > gPathL && massh < 1.2) {
               if (origs[1] > 0) origs[1] = -1;
               qs[0]    = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeDaughterMes)->Charge() / 3); // pion
               etas[0]  = piTr.Momentum3().Eta();
               ps[0]    = piTr.Momentum();
               pts[0]   = piTr.Pt();
               chi2s[0] = TMath::Min(pion->Chi2Vertex(fMpdVert), 999.);
               c2s[0]   = piTr.GetChi2() / (piTr.GetNofTrHits() * 2 - 5);
               MpdHelix helix1 = MakeHelix(&piTr);
               // Get 3-D DCA to primary vertex
               s   = helix1.pathLength(mPrimaryVertex);
               pca = helix1.at(s);
               pca -= mPrimaryVertex;
               dcas[0] = pca.Mag();

               chi2h   = chi2;
               angle   = v0.Angle(lambPart.Momentum3());
               pth     = lambPart.Pt();       // reconstructed
               ph      = lambPart.Momentum(); // reconstructed
               phi_Lam = lambPart.Phi();      // reconstructed
               if (pth > 0.001)
                  etah = lambPart.Momentum3().Eta();
               else
                  etah = TMath::Sign(100., lambPart.Momentum3().Z());
               lambPart.SetMass(massHyperon); // set true mass
               yh                         = lambPart.Rapidity();
               pair<double, double> paths = helix.pathLengths(helix1);
               TVector3             p1    = helix.at(paths.first);
               TVector3             p2    = helix1.at(paths.second);
               p1 -= p2;
               disth = p1.Mag(); // closest distance between daughters
               mpdg  = -9;       // mixing

               // Get 3-D DCA of lambda to primary vertex
               MpdHelix helix2 = MakeHelix(&lambPart);
               s               = helix2.pathLength(mPrimaryVertex);
               pca             = helix2.at(s);
               pca -= mPrimaryVertex;
               dca    = pca.Mag();
               c2pv   = TMath::Min(lambPart.Chi2Vertex(fMpdVert), 9999.);
               omega1 = dcas[0] * dcas[1] / (dca * dca + disth * disth);
               omega2 = TMath::Sqrt(chi2s[0] * chi2s[1]) / (c2pv + chi2h);

               polarhx = 0.0;
               polarhy = 0.0;
               polarhz = 0.0;
               if (origs[0] > 0) {
                  double weight_pol = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetWeight();
                  polarhx           = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetPolar(0);
                  polarhy           = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetPolar(1);
                  polarhz           = ((MpdMCTrack *)mMCTracks->UncheckedAt(mothId))->GetPolar(2);
                  TVector3 polar_changed(polarhx, polarhy, polarhz);
                  if (phiEP != 0.) polar_changed.RotateZ(-phiEP);
                  polarhx = weight_pol * polar_changed.X();
                  polarhy = weight_pol * polar_changed.Y();
                  polarhz = weight_pol * polar_changed.Z();
               }

               FindPolarAngle(lambPart, vPart);

               MpdLambdaPol l0(massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas, mcthetas, thetas, mcphis,
                               phis, mcps, ps, pts, chi2s, dcas, dca, c2pv, omega1, omega2, cosA, cosAmc, polarhx,
                               polarhy, polarhz, phi_star, phi_star_MC, phi_Lam, c2s, origs, qs, layMx, mpdg);
               vLambdas.push_back(l0);
            } // if (chi2 >= 0 && chi2 < gC2L

            int nPart = vPart.size();
            for (int ipart = 0; ipart < nPart; ipart++) delete vPart[ipart];
         } // for (multimap<int,MpdTpcKalmanTrack>::iterator mit1 = piter.first;
      }
   } // for (int ip = 0; ip < nP;
}

void MpdGlobalPolarizationRECO::FindPolarAngle(MpdParticle &lamb, vector<MpdParticle *> &vPart)
{
   TVector3       vPr, vPi, vLamb;
   TLorentzVector prLor, piLor, lambLor;

   // True (exact) parameters
   vPr.SetMagThetaPhi(mcps[1], mcthetas[1], mcphis[1]);
   vPi.SetMagThetaPhi(mcps[0], mcthetas[0], mcphis[0]);

   prLor.SetVectM(vPr, massDaughterBar);
   piLor.SetVectM(vPi, massDaughterMes);

   lambLor = prLor + piLor;

   TVector3 boostV_MC;
   boostV_MC = lambLor.BoostVector();
   boostV_MC *= -1;

   prLor.Boost(boostV_MC);
   vPr = prLor.Vect();

   cosAmc      = vPr.CosTheta();
   phi_star_MC = vPr.Phi();

   // Reco parameters
   vLamb.SetMagThetaPhi(lamb.Momentum(), lamb.Theta(), lamb.Phi());
   lambLor.SetVectM(vLamb, massHyperon);

   MpdParticle *prot = vPart[0];
   vPr.SetMagThetaPhi(prot->Momentum(), prot->Theta(), prot->Phi());

   prLor.SetVectM(vPr, massDaughterBar);

   TVector3 boostV;
   boostV = lambLor.BoostVector();
   boostV *= -1;

   prLor.Boost(boostV);
   vPr = prLor.Vect();

   cosA     = vPr.CosTheta();
   phi_star = vPr.Phi(); // azimuthal angle of proton in the lambda frame (phi*)
}
double *MpdGlobalPolarizationRECO::init_double_array(const int n, const double fmt...)
{
   va_list args;
   va_start(args, fmt);

   double *ret = new double[n];

   for (int i = 0; i < n; i++) {
      ret[i] = va_arg(args, double);
   }

   return ret;
}

int *MpdGlobalPolarizationRECO::init_int_array(const int n, const int fmt...)
{
   va_list args;
   va_start(args, fmt);

   int *ret = new int[n];

   for (int i = 0; i < n; i++) {
      ret[i] = va_arg(args, int);
   }

   return ret;
}
