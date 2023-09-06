#include <iostream>
#include <fstream> // std::ifstream

#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdMCTrack.h"
#include "MpdGlobalPolarizationMC.h"
#include "TFile.h"

ClassImp(MpdGlobalPolarizationMC);

MpdGlobalPolarizationMC::MpdGlobalPolarizationMC(const char *name, const char *outputName)
   : MpdAnalysisTask2(name, outputName)
{
   readParameters(name);
   param("NITER_CENT", NITER_CENT, 4);
   param("NITER", NITER, 20);
   param("cent_cut_choice", cent_cut_choice, 0);
   param("cent_cut", cent_cut, 100.0);
   param("particle_choice", particle_choice, 3122);
}

void MpdGlobalPolarizationMC::UserInit()
{
   // Initializing list of output histograms
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // Choice of analyzed particle
   pdgCodeHyperon = particle_choice;
   if (pdgCodeHyperon == pdgCodeL0) {
      cout << "You have chosen to analyze Lambda hyperons: "
           << " pdg: " << pdgCodeHyperon << endl;
      pdgCodeDaughter = pdgCodePr;
      massHyperon     = massL0;
      massDaughter    = massPr;
   } else if (pdgCodeHyperon == pdgCodeAL0) {
      cout << "You have chosen to analyze anti-Lambda hyperons: "
           << " pdg: " << pdgCodeHyperon << endl;
      pdgCodeDaughter = pdgCodeAPr;
      massHyperon     = massL0;
      massDaughter    = massPr;
   } else {
      cout << "This pdg code for particle_choice is not defined! Please provide the definition in the code." << endl;
      exit(0);
   }
   cout << "massHyperon: " << massHyperon << endl;
   cout << "massDaughter: " << massDaughter << endl;

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
      exit(0);
   }

   // Initializing output histograms
   hCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);
   fOutputList->Add(hCentrality);

   hNevCentr = new TH1D("hNevCentr", "Events in centrality bins", NITER_CENT, _CentrBins);
   fOutputList->Add(hNevCentr);

   hResolution_EP1_true = new TH1D("hResolution_EP1_true", "True EP1 resolution", NITER_CENT, _CentrBins);
   fOutputList->Add(hResolution_EP1_true);

   hResolution_EP1_reco = new TH1D("hResolution_EP1_reco", "Reco EP1 resolution", NITER_CENT, _CentrBins);
   fOutputList->Add(hResolution_EP1_reco);

   hPolarY_Full     = new TH1D *[NITER_CENT];
   hPolarY_Prim     = new TH1D *[NITER_CENT];
   hDeltaPhiRP_Full = new TH1D *[NITER_CENT];
   hDeltaPhiRP_Prim = new TH1D *[NITER_CENT];
   hDeltaPhiEP_Full = new TH1D *[NITER_CENT];
   hDeltaPhiEP_Prim = new TH1D *[NITER_CENT];

   hPolarY_Full_pt_eta_bin = new TH1D ***[NITER_CENT];
   hPolarY_Prim_pt_eta_bin = new TH1D ***[NITER_CENT];

   hPolarY2_Prim_pt_eta_bin = new TH1D ***[NITER_CENT];


   hinvmass = new TH1D *[NITER_CENT];
   hpT_Full = new TH1D *[NITER_CENT];
   hpT_Prim = new TH1D *[NITER_CENT];
   heta_Full = new TH1D *[NITER_CENT];
   heta_Prim = new TH1D *[NITER_CENT];
   hv1EtaPhi_Full = new TProfile **[NITER_CENT];
   hv1EtaPhi_Prim = new TProfile **[NITER_CENT];
   hv1EtaPsi_Full = new TProfile **[NITER_CENT];
   hv1EtaPsi_Prim = new TProfile **[NITER_CENT];

   hv1pTPhi_Full = new TProfile **[NITER_CENT];
   hv1pTPhi_Prim = new TProfile **[NITER_CENT];
   hv1pTPsi_Full = new TProfile **[NITER_CENT];
   hv1pTPsi_Prim = new TProfile **[NITER_CENT];

   hPv1pTPsi_Prim = new TProfile **[NITER_CENT];
   hPv1EtaPsi_Prim = new TProfile **[NITER_CENT];
   hNv1pTPsi_Prim = new TProfile **[NITER_CENT];
   hNv1EtaPsi_Prim = new TProfile **[NITER_CENT];
   hNPpTPsi_Prim = new TProfile **[NITER_CENT];
   hNPEtaPsi_Prim = new TProfile **[NITER_CENT];

   hPpTPsi_Prim = new TProfile **[NITER_CENT];
   hPEtaPsi_Prim = new TProfile **[NITER_CENT];
   hP2pTPsi_Prim = new TProfile **[NITER_CENT];
   hP2EtaPsi_Prim = new TProfile **[NITER_CENT];

   hv12pTPsi_Prim = new TProfile **[NITER_CENT];
   hv12EtaPsi_Prim = new TProfile **[NITER_CENT];


   hP4pTPsi_Prim = new TProfile **[NITER_CENT];
   hP4EtaPsi_Prim = new TProfile **[NITER_CENT];
   hv14pTPsi_Prim = new TProfile **[NITER_CENT];
   hv14EtaPsi_Prim = new TProfile **[NITER_CENT];

   hP2v1pTPsi_Prim = new TProfile **[NITER_CENT];
   hP2v1EtaPsi_Prim = new TProfile **[NITER_CENT];
   hPv12pTPsi_Prim = new TProfile **[NITER_CENT];
   hPv12EtaPsi_Prim = new TProfile **[NITER_CENT];
   hPv1Psi_Prim = new TProfile *[NITER_CENT];
   hPv1Psi_Prim_cut = new TProfile *[NITER_CENT];
   hPv1Psi_Prim_Neg = new TProfile *[NITER_CENT];
   hPv1Psi_Prim_cut_Neg = new TProfile *[NITER_CENT];

   hP_invmass_Psi = new TProfile *[NITER_CENT];
   hv1_invmass_Psi = new TProfile *[NITER_CENT];
   const double BinEdges[] = {0.,10.,20.,50.,100.};
   hN = new TProfile("hN","hN;centrality;N",4,BinEdges);
   hN2 = new TProfile("hN2","hN2;centrality;N2",4,BinEdges);
   fOutputList->Add(hN);
   fOutputList->Add(hN2);


   

   const int NITER_PT = 5;
   const int NITER_ETA = 6;
   const double pt_edges[6] = {0.,0.5,1.,1.5,2.,3.};
   const double eta_edges[7] = {-1.5,-1.,-0.5,0.,0.5,1.,1.5};
   hv1CentPhi_Full = new TProfile ("hv1CentPhi_Full","v1(cent);cent;v1",4,BinEdges);
   fOutputList->Add(hv1CentPhi_Full);
   hv1CentPhi_Prim = new TProfile ("hv1CentPhi_Prim","v1(cent);cent;v1",4,BinEdges);
   fOutputList->Add(hv1CentPhi_Prim);
   hv1CentPsi_Full = new TProfile ("hv1CentPsi_Full","v1(cent);cent;v1",4,BinEdges);
   fOutputList->Add(hv1CentPsi_Full);
   hv1CentPsi_Prim = new TProfile ("hv1CentPsi_Prim","v1(cent);cent;v1",4,BinEdges);
   fOutputList->Add(hv1CentPsi_Prim);

   ResEP1_true = new double[NITER_CENT];
   SubEvRes1   = new double[NITER_CENT];

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      hP_invmass_Psi[iter_cent] = new TProfile(Form("hP_invmass_Psi_cent%d",iter_cent),Form("hP_invmass_Psi_cent%d;inv_mass, GeV/c;sin(phiRP-phiP)",iter_cent),100,1.07,1.17);
      fOutputList->Add(hP_invmass_Psi[iter_cent]);
      hv1_invmass_Psi[iter_cent] = new TProfile(Form("hv1_invmass_Psi_cent%d",iter_cent),Form("hv1_invmass_Psi_cent%d;inv_mass, GeV/c;v1",iter_cent),100,1.07,1.17);
      fOutputList->Add(hv1_invmass_Psi[iter_cent]);

      hPv1Psi_Prim[iter_cent] = new TProfile(Form("hPv1Psi_Prim_cent%d",iter_cent),Form("hPv1Psi_Prim_cent%d;v1;Py,%",iter_cent),5,-1,1);
      fOutputList->Add(hPv1Psi_Prim[iter_cent]);
      hPv1Psi_Prim_cut[iter_cent] = new TProfile(Form("hPv1Psi_Prim_cut_cent%d",iter_cent),Form("hPv1Psi_Prim_cut_cent%d;v1;Py,%",iter_cent),5,-1,1);
      fOutputList->Add(hPv1Psi_Prim_cut[iter_cent]);
      hPv1Psi_Prim_Neg[iter_cent] = new TProfile(Form("hPv1Psi_Prim_Neg_cent%d",iter_cent),Form("hPv1Psi_Prim_Neg_cent%d;v1;Py,%",iter_cent),5,-1,1);
      fOutputList->Add(hPv1Psi_Prim_Neg[iter_cent]);
      hPv1Psi_Prim_cut_Neg[iter_cent] = new TProfile(Form("hPv1Psi_Prim_cut_Neg_cent%d",iter_cent),Form("hPv1Psi_Prim_cut_Neg_cent%d;v1;Py,%",iter_cent),5,-1,1);
      fOutputList->Add(hPv1Psi_Prim_cut_Neg[iter_cent]);

      hPolarY_Full[iter_cent] =
         new TH1D(Form("hPolarY_Full_cent%d", iter_cent), Form("hPolarY_Full_cent%d", iter_cent), 100, -1., 1.);
      fOutputList->Add(hPolarY_Full[iter_cent]);
      hPolarY_Prim[iter_cent] =
         new TH1D(Form("hPolarY_Prim_cent%d", iter_cent), Form("hPolarY_Prim_cent%d", iter_cent), 100, -1., 1.);
      fOutputList->Add(hPolarY_Prim[iter_cent]);
      hDeltaPhiRP_Full[iter_cent] = new TH1D(Form("hDeltaPhiRP_Full_cent%d", iter_cent),
                                             Form("hDeltaPhiRP_Full_cent%d", iter_cent), NITER, 0., 2. * pi);
      fOutputList->Add(hDeltaPhiRP_Full[iter_cent]);
      hDeltaPhiRP_Prim[iter_cent] = new TH1D(Form("hDeltaPhiRP_Prim_cent%d", iter_cent),
                                             Form("hDeltaPhiRP_Prim_cent%d", iter_cent), NITER, 0., 2. * pi);
      fOutputList->Add(hDeltaPhiRP_Prim[iter_cent]);
      hDeltaPhiEP_Full[iter_cent] = new TH1D(Form("hDeltaPhiEP_Full_cent%d", iter_cent),
                                             Form("hDeltaPhiEP_Full_cent%d", iter_cent), NITER, 0., 2. * pi);
      fOutputList->Add(hDeltaPhiEP_Full[iter_cent]);
      hDeltaPhiEP_Prim[iter_cent] = new TH1D(Form("hDeltaPhiEP_Prim_cent%d", iter_cent),
                                             Form("hDeltaPhiEP_Prim_cent%d", iter_cent), NITER, 0., 2. * pi);
      fOutputList->Add(hDeltaPhiEP_Prim[iter_cent]);

      hpT_Full[iter_cent] =
         new TH1D(Form("hpT_Full_cent%d", iter_cent), Form("hpT_Full_cent%d", iter_cent), 300, 0., 3.);
      fOutputList->Add(hpT_Full[iter_cent]);
      hpT_Prim[iter_cent] =
         new TH1D(Form("hpT_Prim_cent%d", iter_cent), Form("hpT_Prim_cent%d", iter_cent), 300, 0., 3.);
      fOutputList->Add(hpT_Prim[iter_cent]);
      heta_Full[iter_cent] =
         new TH1D(Form("heta_Full_cent%d", iter_cent), Form("heta_Full_cent%d", iter_cent), 300, -1.5, 1.5);
      fOutputList->Add(heta_Full[iter_cent]);
      heta_Prim[iter_cent] =
         new TH1D(Form("heta_Prim_cent%d", iter_cent), Form("heta_Prim_cent%d", iter_cent), 300, -1.5, 1.5);
      fOutputList->Add(heta_Prim[iter_cent]);
      hinvmass[iter_cent] = 
          new TH1D(Form("hinvmass_cent%d", iter_cent), Form("hinvmass_cent%d", iter_cent), 100, 1.070, 1.170);
      fOutputList->Add(hinvmass[iter_cent]);
     
      

      hv1EtaPhi_Full[iter_cent] = new TProfile *[NITER_PT];
   hv1EtaPhi_Prim[iter_cent] = new TProfile *[NITER_PT];
   hv1EtaPsi_Full[iter_cent] = new TProfile *[NITER_PT];
   hv1EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];

   hv1pTPhi_Full[iter_cent] = new TProfile *[NITER_ETA];
   hv1pTPhi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hv1pTPsi_Full[iter_cent] = new TProfile *[NITER_ETA];
   hv1pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];  


   hPv1pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hPv1EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];
   hNv1pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hNv1EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];
   hNPpTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hNPEtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];

   hPpTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hPEtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];
   hP2pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hP2EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];
   hv12pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hv12EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];

   hP4pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hP4EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];
   hv14pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hv14EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];

   hP2v1pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hP2v1EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];
   hPv12pTPsi_Prim[iter_cent] = new TProfile *[NITER_ETA];
   hPv12EtaPsi_Prim[iter_cent] = new TProfile *[NITER_PT];

   for(int iter_pt =0;iter_pt<NITER_PT;iter_pt++){
      hv1EtaPhi_Full[iter_cent][iter_pt] = new TProfile (Form("hv1EtaPhi_Full%d_%d", iter_cent,iter_pt), Form("hv1EtaPhi_Full%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hv1EtaPhi_Full[iter_cent][iter_pt]);
            hv1EtaPhi_Prim[iter_cent][iter_pt] = new TProfile (Form("hv1EtaPhi_Prim%d_%d", iter_cent,iter_pt), Form("hv1EtaPhi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hv1EtaPhi_Prim[iter_cent][iter_pt]);
            hv1EtaPsi_Full[iter_cent][iter_pt] = new TProfile (Form("hv1EtaPsi_Full%d_%d", iter_cent,iter_pt), Form("hv1EtaPsi_Full%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hv1EtaPsi_Full[iter_cent][iter_pt]);
            hv1EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hv1EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hv1EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hv1EtaPsi_Prim[iter_cent][iter_pt]);

            hPv1EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hPv1EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hPv1EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hPv1EtaPsi_Prim[iter_cent][iter_pt]);
            hNPEtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hNPEtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hNPEtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hNPEtaPsi_Prim[iter_cent][iter_pt]);
            hPEtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hPEtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hPEtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hPEtaPsi_Prim[iter_cent][iter_pt]);
            hP2EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hP2EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hP2EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hP2EtaPsi_Prim[iter_cent][iter_pt]);
            hNv1EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hNv1EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hNv1EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hNv1EtaPsi_Prim[iter_cent][iter_pt]);
            hv12EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hv12EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hv12EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hv12EtaPsi_Prim[iter_cent][iter_pt]);


            hPv12EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hPv12EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hPv12EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hPv12EtaPsi_Prim[iter_cent][iter_pt]);
            hP2v1EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hP2v1EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hP2v1EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hP2v1EtaPsi_Prim[iter_cent][iter_pt]);
            hv14EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hv14EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hv14EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hv14EtaPsi_Prim[iter_cent][iter_pt]);
            hP4EtaPsi_Prim[iter_cent][iter_pt] = new TProfile (Form("hP4EtaPsi_Prim%d_%d", iter_cent,iter_pt), Form("hP4EtaPsi_Prim%d_%d", iter_cent,iter_pt), 15, -1.5, 1.5);
            fOutputList->Add(hP4EtaPsi_Prim[iter_cent][iter_pt]);
    }
    for(int iter_eta =0;iter_eta<NITER_ETA;iter_eta++){
      hv1pTPhi_Full[iter_cent][iter_eta] = new TProfile (Form("hv1pTPhi_Full%d_%d", iter_cent,iter_eta), Form("hv1pTPhi_Full%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hv1pTPhi_Full[iter_cent][iter_eta]);
            hv1pTPhi_Prim[iter_cent][iter_eta] = new TProfile (Form("hv1pTPhi_Prim%d_%d", iter_cent,iter_eta), Form("hv1pTPhi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hv1pTPhi_Prim[iter_cent][iter_eta]);
            hv1pTPsi_Full[iter_cent][iter_eta] = new TProfile (Form("hv1pTPsi_Full%d_%d", iter_cent,iter_eta), Form("hv1pTPsi_Full%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hv1pTPsi_Full[iter_cent][iter_eta]);
            hv1pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hv1pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hv1pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hv1pTPsi_Prim[iter_cent][iter_eta]);

            hPv1pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hPv1pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hPv1pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hPv1pTPsi_Prim[iter_cent][iter_eta]);
            hNPpTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hNPpTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hNPpTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hNPpTPsi_Prim[iter_cent][iter_eta]);
            hPpTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hPpTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hPpTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hPpTPsi_Prim[iter_cent][iter_eta]);
            hP2pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hP2pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hP2pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hP2pTPsi_Prim[iter_cent][iter_eta]);
            hNv1pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hNv1pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hNv1pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hNv1pTPsi_Prim[iter_cent][iter_eta]);
            hv12pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hv12pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hv12pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hv12pTPsi_Prim[iter_cent][iter_eta]);

            hv14pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hv14pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hv14pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hv14pTPsi_Prim[iter_cent][iter_eta]);
            hP4pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hP4pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hP4pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hP4pTPsi_Prim[iter_cent][iter_eta]);

            hP2v1pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hP2v1pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hP2v1pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hP2v1pTPsi_Prim[iter_cent][iter_eta]);
            hPv12pTPsi_Prim[iter_cent][iter_eta] = new TProfile (Form("hPv12pTPsi_Prim%d_%d", iter_cent,iter_eta), Form("hPv12pTPsi_Prim%d_%d", iter_cent,iter_eta), 15, 0., 3.);
            fOutputList->Add(hPv12pTPsi_Prim[iter_cent][iter_eta]);
    }

      ResEP1_true[iter_cent] = 0.;
      SubEvRes1[iter_cent]   = 0.;
   }
   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++){
      hPolarY_Full_pt_eta_bin[iter_cent] = new TH1D **[NITER_PT];
      hPolarY_Prim_pt_eta_bin[iter_cent] = new TH1D **[NITER_PT];
      hPolarY2_Prim_pt_eta_bin[iter_cent] = new TH1D **[NITER_PT];
      for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++){
         hPolarY_Full_pt_eta_bin[iter_cent][iter_pt] = new TH1D *[NITER_ETA];
         hPolarY_Prim_pt_eta_bin[iter_cent][iter_pt] = new TH1D *[NITER_ETA];
         hPolarY2_Prim_pt_eta_bin[iter_cent][iter_pt] = new TH1D *[NITER_ETA];
         for(int iter_eta =0; iter_eta< NITER_ETA;iter_eta++){
            hPolarY_Full_pt_eta_bin[iter_cent][iter_pt][iter_eta] =
         new TH1D(Form("hPolarY_Full_pt_eta_bin%d_%d_%d", iter_cent,iter_pt,iter_eta), Form("hPolarY_Full_pt_eta_bin%d_%d_%d", iter_cent,iter_pt,iter_eta), 100, -1., 1.);
            fOutputList->Add(hPolarY_Full_pt_eta_bin[iter_cent][iter_pt][iter_eta]);
            hPolarY_Prim_pt_eta_bin[iter_cent][iter_pt][iter_eta] =
         new TH1D(Form("hPolarY_Prim_pt_eta_bin%d_%d_%d", iter_cent,iter_pt,iter_eta), Form("hPolarY_Prim_pt_eta_bin%d_%d_%d", iter_cent,iter_pt,iter_eta), 100, -1., 1.);
            fOutputList->Add(hPolarY_Prim_pt_eta_bin[iter_cent][iter_pt][iter_eta]);

          hPolarY2_Prim_pt_eta_bin[iter_cent][iter_pt][iter_eta] =
         new TH1D(Form("hPolarY2_Prim_pt_eta_bin%d_%d_%d", iter_cent,iter_pt,iter_eta), Form("hPolarY2_Prim_pt_eta_bin%d_%d_%d", iter_cent,iter_pt,iter_eta), 100, -1., 1.);
            fOutputList->Add(hPolarY2_Prim_pt_eta_bin[iter_cent][iter_pt][iter_eta]);
         }
      }   
   }
   
   phiRP    = 0.;
   phiEP    = 0.;
   ResEP    = 0.;
   ResEPSub = 0.;
}

void MpdGlobalPolarizationMC::ProcessEvent(MpdAnalysisEvent &event)
{
   if (!selectEvent(event)) return;

   // Calculate the reaction/event plane and its resolution
   phiRP    = event.fMCEventHeader->GetRotZ();
   phiRP    = TMath::ATan2(TMath::Sin(phiRP), TMath::Cos(phiRP));
   phiEP    = event.fMpdEP.GetPhiEP_FHCal_F_all();
   ResEP    = TMath::Cos(phiEP - phiRP);
   ResEPSub = TMath::Cos(event.fMpdEP.GetPhiEP_FHCal_S_all() - event.fMpdEP.GetPhiEP_FHCal_N_all());
   

   fillHistograms(event);
}

void MpdGlobalPolarizationMC::Finish()
{
   cout << "Finish() ..." << endl;
}

//--------------------------------------
bool MpdGlobalPolarizationMC::selectEvent(MpdAnalysisEvent &event)
{
   mMCTracks = event.fMCTrack;

   Centrality_tpc = event.getCentrTPC();

   if (Centrality_tpc < 0 || Centrality_tpc >= 100) // TPC centrality not defined
      return false;

   if (cent_cut_choice == 1) {
      cout << "Cutting out events with more than " << cent_cut << "% centrality!" << endl;
      if (Centrality_tpc > cent_cut) return false;
   }

   return true;
}

void MpdGlobalPolarizationMC::fillHistograms(MpdAnalysisEvent &event)
{  
   const int NITER_PT = 5;
   const int NITER_ETA = 6;
   const double pt_edges[6] = {0.,0.5,1.,1.5,2.,3.};
   const double eta_edges[7] = {-1.5,-1.,-0.5,0.,0.5,1.,1.5};


   hCentrality->Fill(Centrality_tpc);
   hNevCentr->Fill(Centrality_tpc);
   int nMC = mMCTracks->GetEntriesFast();

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      if (Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue;
      vector<int> vecP, vecPi;

      ResEP1_true[iter_cent] += ResEP;
      SubEvRes1[iter_cent] += ResEPSub;
      int N =0;

for (int j = 0; j < nMC; ++j) // iterating over MCTracks
      {
         MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(j);
         if (mcTr->GetPdgCode() == pdgCodeDaughter) {
            int mcTr_MotherID = mcTr->GetMotherId();
            if (mcTr_MotherID < 0) continue;
            MpdMCTrack *mcTr_Mother = (MpdMCTrack *)mMCTracks->UncheckedAt(mcTr_MotherID);
            if (mcTr_Mother->GetPdgCode() ==
                pdgCodeHyperon)
            {
             int mcTr_Lam_MotherID = mcTr_Mother->GetMotherId();
             if(mcTr_Lam_MotherID<0)
             N++;
            }
           }

if(mcTr->GetMotherId()!=-1 && TMath::Abs(mcTr->GetRapidity())<1.5){
   int mcTr_MotherID = mcTr->GetMotherId();
   MpdMCTrack *mcTr_Mother = (MpdMCTrack *)mMCTracks->UncheckedAt(mcTr_MotherID);
   if (mcTr_Mother->GetPdgCode() == pdgCodeHyperon){  
if(mcTr->GetPdgCode() == 2212)
vecP.push_back(j);
else if(mcTr->GetPdgCode() == -211)
vecPi.push_back(j);
  }
}
}
int nP = vecP.size();
int nPi = vecPi.size();
for(int ip = 0; ip<nP;ip++)
{
 MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(vecP[ip]);
 TVector3 momP;
 mcTr->GetMomentum(momP);
 double p_E = mcTr->GetEnergy();
double p_phi = momP.Phi();
for (int ipi =0; ipi<nPi;ipi++){
MpdMCTrack *mcTr2 = (MpdMCTrack *)mMCTracks->UncheckedAt(vecPi[ipi]);
 TVector3 momPi;
 mcTr2->GetMomentum(momPi);
 double pi_E = mcTr2->GetEnergy();
 double lam_pt;
 double lam_y;
 double phi_pair;
lam_pt = (momPi+momP).Pt();
lam_y = 0.5*TMath::Log(((p_E+pi_E)+(momPi+momP).Z())/((p_E+pi_E)-(momPi+momP).Z()));
phi_pair = (momPi+momP).Phi();
double inv_mass;
inv_mass = TMath::Sqrt((p_E+pi_E)*(p_E+pi_E) - (momP+momPi)*(momP+momPi));
hP_invmass_Psi[iter_cent]->Fill(inv_mass,TMath::Sin(phiRP-p_phi));
hinvmass[iter_cent]->Fill(inv_mass);
if(lam_pt>0.5)
if(lam_pt<2)
if(lam_y>0.5)
if(lam_y<1)
hv1_invmass_Psi[iter_cent]->Fill(inv_mass, TMath::Cos(phi_pair-phiRP));
}
 
  
}
vecP.clear();
vecPi.clear();
hN->Fill(Centrality_tpc,N);
hN2->Fill(Centrality_tpc,N*N);

      for (int j = 0; j < nMC; ++j) // iterating over MCTracks
      {
         TVector3    mom;
         TVector3    mom_moth;
         MpdMCTrack *mcTr = (MpdMCTrack *)mMCTracks->UncheckedAt(j);
         mcTr->GetMomentum(mom);
         if (mcTr->GetPdgCode() == pdgCodeDaughter) {
            int mcTr_MotherID = mcTr->GetMotherId();
            if (mcTr_MotherID < 0) continue; // if it's a primary particle (e.g. proton), we don't need it

            // Get the track, corresponding to the mother ID
            MpdMCTrack *mcTr_Mother = (MpdMCTrack *)mMCTracks->UncheckedAt(mcTr_MotherID);
            if (mcTr_Mother->GetPdgCode() ==
                pdgCodeHyperon) // choose only the case, when it's the mother we need (e.g. Lambda)
            {
               mcTr_Mother->GetMomentum(mom_moth);
               int mcTr_Lam_MotherID = mcTr_Mother->GetMotherId(); // ID of the mother of hyperon (e.g. Lambda)

               // Determine the polarization vector from the hyperon MC information
               double polar_x    = 0.;
               double polar_y    = 0.;
               double polar_z    = 0.;
               double weight_pol = 0.;

               weight_pol = mcTr_Mother->GetWeight();
               polar_x    = mcTr_Mother->GetPolar(0);
               polar_y    = mcTr_Mother->GetPolar(1);
               polar_z    = mcTr_Mother->GetPolar(2);

               TVector3 polar_vector(polar_x, polar_y, polar_z); // polarization vector for distributions

               // Rotate the polarization vector back w.r.t. the generated reaction plane
               if (phiRP != 0.) polar_vector.RotateZ(-phiRP);

               // Restore the values of the original polarization vector (from the model)
               polar_x = weight_pol * polar_vector.X();
               polar_y = weight_pol * polar_vector.Y();
               polar_z = weight_pol * polar_vector.Z();

               // Fill model global polarization distribution for primary hyperons
               if (mcTr_Lam_MotherID < 0) {
                  hPolarY_Prim[iter_cent]->Fill(polar_y);
               }

               // Fill model global polarization distribution for full hyperons
               hPolarY_Full[iter_cent]->Fill(polar_y);

               // daughter (e.g. proton) values (MagThetaPhi)
               double p_prot   = mom.Mag();
               double eta_prot = mom.Theta();
               double phi_prot = mom.Phi();
               // hyperon (e.g. Lambda) values (MagThetaPhi)
               double   p_lam   = mom_moth.Mag();
               double   eta_lam = mom_moth.Theta();
               double   phi_lam = mom_moth.Phi();
               TVector3 vPr, vLamb;
               vPr.SetMagThetaPhi(p_prot, eta_prot, phi_prot);
               vLamb.SetMagThetaPhi(p_lam, eta_lam, phi_lam);

               // calculate the costheta and phi of daughter particle in the rest frame of the hyperon
               double cos_prot      = 0.0;
               double phi_prot_star = 0.0;
               FindPolarAngle(vPr, vLamb, cos_prot, phi_prot_star);

               // Calculate the difference between the RP(EP) angle and the azimuthal angle of proton
               double phi_diff_RP = 0.;
               double phi_diff_EP = 0.;
               phi_diff_RP        = phiRP - phi_prot_star;
               phi_diff_EP        = phiEP - phi_prot_star;

               if (phi_diff_RP < 0) phi_diff_RP = phi_diff_RP + 2. * pi;
               if (phi_diff_EP < 0) phi_diff_EP = phi_diff_EP + 2. * pi;

               // Fill angular distribution for all hyperons
               hDeltaPhiRP_Full[iter_cent]->Fill(phi_diff_RP);
               hDeltaPhiEP_Full[iter_cent]->Fill(phi_diff_EP);
                
               // Fill QA and flow distribution for all hyperons
               hpT_Full[iter_cent]->Fill(mcTr_Mother->GetPt());
               heta_Full[iter_cent]->Fill(mcTr_Mother->GetRapidity());
               for(int iter_pt=0;iter_pt<NITER_PT;iter_pt++){
                  if((mcTr_Mother->GetPt())>pt_edges[iter_pt] && (mcTr_Mother->GetPt())<pt_edges[iter_pt+1]){
                  hv1EtaPhi_Full[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),TMath::Cos(phi_lam - phiEP));
                   hv1EtaPsi_Full[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),TMath::Cos(phi_lam - phiRP));
                      }
               }
               for(int iter_eta=0;iter_eta<NITER_ETA;iter_eta++){
                  if((mcTr_Mother->GetRapidity())>eta_edges[iter_eta] && (mcTr_Mother->GetRapidity())<eta_edges[iter_eta+1]){
                   hv1pTPhi_Full[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),TMath::Cos(phi_lam - phiEP));
               hv1pTPsi_Full[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),TMath::Cos(phi_lam - phiRP));}}
               if(mcTr_Mother->GetRapidity()>0 && mcTr_Mother->GetRapidity()<1.5 && mcTr_Mother->GetPt()>0.5 && mcTr_Mother->GetPt()<3.){
               hv1CentPhi_Full->Fill(Centrality_tpc,TMath::Cos(phi_lam - phiEP));
               hv1CentPsi_Full->Fill(Centrality_tpc,TMath::Cos(phi_lam - phiRP));}

               for(int iter_pt=0;iter_pt<NITER_PT;iter_pt++){
                 if((mcTr_Mother->GetPt())>pt_edges[iter_pt] && (mcTr_Mother->GetPt())<pt_edges[iter_pt+1]){
                   for(int iter_eta=0;iter_eta<NITER_ETA;iter_eta++){
                          if((mcTr_Mother->GetRapidity())>eta_edges[iter_eta] && (mcTr_Mother->GetRapidity())<eta_edges[iter_eta+1]){
                             hPolarY_Full_pt_eta_bin[iter_cent][iter_pt][iter_eta]->Fill(polar_y);
                           }
                        }
                     }
                }


               // Fill angular distribution for primary hyperons
               if (mcTr_Lam_MotherID < 0) {
                  hDeltaPhiRP_Prim[iter_cent]->Fill(phi_diff_RP);
                  hDeltaPhiEP_Prim[iter_cent]->Fill(phi_diff_EP);

               // Fill QA and flow distribution for primary hyperons
               hpT_Prim[iter_cent]->Fill(mcTr_Mother->GetPt());
               heta_Prim[iter_cent]->Fill(mcTr_Mother->GetRapidity());
               if(TMath::Abs(mcTr_Mother->GetRapidity())<1 && polar_z>0)
               hPv1Psi_Prim[iter_cent]->Fill(TMath::Cos(phi_lam - phiRP),polar_y);
               if(TMath::Abs(mcTr_Mother->GetRapidity())<1 && polar_z<0)
               hPv1Psi_Prim_Neg[iter_cent]->Fill(TMath::Cos(phi_lam - phiRP),polar_y);
               if(TMath::Abs(mcTr_Mother->GetRapidity())<1 && TMath::Abs(mcTr_Mother->GetRapidity())<0.5 && mcTr_Mother->GetPt()>0.5 && mcTr_Mother->GetPt()<3 && polar_z>0)
               hPv1Psi_Prim_cut[iter_cent]->Fill(TMath::Cos(phi_lam - phiRP),polar_y);
               if(TMath::Abs(mcTr_Mother->GetRapidity())<1 && TMath::Abs(mcTr_Mother->GetRapidity())<0.5 && mcTr_Mother->GetPt()>0.5 && mcTr_Mother->GetPt()<3 && polar_z<0)
               hPv1Psi_Prim_cut_Neg[iter_cent]->Fill(TMath::Cos(phi_lam - phiRP),polar_y);
               
               for(int iter_pt=0;iter_pt<NITER_PT;iter_pt++){
                  if((mcTr_Mother->GetPt())>pt_edges[iter_pt] && (mcTr_Mother->GetPt())<pt_edges[iter_pt+1]){
                  hv1EtaPhi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),TMath::Cos(phi_lam - phiEP));
                   hv1EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),TMath::Cos(phi_lam - phiRP));
                    hv12EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP)));
                   hPv1EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),polar_y*TMath::Cos(phi_lam - phiRP));
                   hNPEtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),polar_y*N);
                   hPEtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),polar_y);
                   hP2EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),polar_y*polar_y);
                   hNv1EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),N*TMath::Cos(phi_lam - phiRP));

                   hv14EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP)));
                   hP4EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),polar_y*polar_y*polar_y*polar_y);
                   hPv12EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),polar_y*TMath::Cos(phi_lam - phiRP)*TMath::Cos(phi_lam - phiRP));
                   hP2v1EtaPsi_Prim[iter_cent][iter_pt]->Fill(mcTr_Mother->GetRapidity(),polar_y*polar_y*TMath::Cos(phi_lam - phiRP));
                      }
               }
               for(int iter_eta=0;iter_eta<NITER_ETA;iter_eta++){
                  if((mcTr_Mother->GetRapidity())>eta_edges[iter_eta] && (mcTr_Mother->GetRapidity())<eta_edges[iter_eta+1]){
                   hv1pTPhi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),TMath::Cos(phi_lam - phiRP));

               hv1pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),TMath::Cos(phi_lam - phiRP));
               hv12pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP)));  
               hPv1pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),polar_y*TMath::Cos(phi_lam - phiRP));
               hNPpTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),polar_y*N);
               hPpTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),polar_y);
               hP2pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),polar_y*polar_y);
               hNv1pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),N*TMath::Cos(phi_lam - phiRP));


               hv14pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP))*(TMath::Cos(phi_lam - phiRP)));
               hP4pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),polar_y*polar_y*polar_y*polar_y);
               hP2v1pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),polar_y*TMath::Cos(phi_lam - phiRP)*polar_y);
               hPv12pTPsi_Prim[iter_cent][iter_eta]->Fill(mcTr_Mother->GetPt(),polar_y*TMath::Cos(phi_lam - phiRP)*TMath::Cos(phi_lam - phiRP));
                }}
               if(mcTr_Mother->GetRapidity()>0 && mcTr_Mother->GetRapidity()<1.5 && mcTr_Mother->GetPt()>0.5 && mcTr_Mother->GetPt()<3.){
               hv1CentPhi_Prim->Fill(Centrality_tpc,TMath::Cos(phi_lam - phiEP));
               hv1CentPsi_Prim->Fill(Centrality_tpc,TMath::Cos(phi_lam - phiRP));}


               for(int iter_pt=0;iter_pt<NITER_PT;iter_pt++){
                 if((mcTr_Mother->GetPt())>pt_edges[iter_pt] && (mcTr_Mother->GetPt())<pt_edges[iter_pt+1]){
                   for(int iter_eta=0;iter_eta<NITER_ETA;iter_eta++){
                          if((mcTr_Mother->GetRapidity())>eta_edges[iter_eta] && (mcTr_Mother->GetRapidity())<eta_edges[iter_eta+1]){
                             hPolarY_Prim_pt_eta_bin[iter_cent][iter_pt][iter_eta]->Fill(polar_y);
                             hPolarY2_Prim_pt_eta_bin[iter_cent][iter_pt][iter_eta]->Fill(polar_y*polar_y);
                           }
                        }
                     }
                }
               }
            } // mcTr_Mother->GetPdgCode() == pdgCodeHyperon
         }    // mcTr->GetPdgCode() == pdgCodeDaughter
      }       // cycle over MC tracks
   }

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      hResolution_EP1_true->SetBinContent(iter_cent + 1, ResEP1_true[iter_cent]);
      hResolution_EP1_reco->SetBinContent(iter_cent + 1, SubEvRes1[iter_cent]);
   }
}
double *MpdGlobalPolarizationMC::init_double_array(const int n, const double fmt...)
{
   va_list args;
   va_start(args, fmt);

   double *ret = new double[n];

   for (int i = 0; i < n; i++) {
      ret[i] = va_arg(args, double);
   }

   return ret;
}

int *MpdGlobalPolarizationMC::init_int_array(const int n, const int fmt...)
{
   va_list args;
   va_start(args, fmt);

   int *ret = new int[n];

   for (int i = 0; i < n; i++) {
      ret[i] = va_arg(args, int);
   }

   return ret;
}

void MpdGlobalPolarizationMC::FindPolarAngle(TVector3 &vPr, TVector3 &vLamb, double &cos_prot, double &phi_prot_star)
{
   cos_prot      = 0.;
   phi_prot_star = 0.;

   TLorentzVector prLor, lambLor;
   prLor.SetVectM(vPr, massDaughter);
   lambLor.SetVectM(vLamb, massHyperon);
   TVector3 boostV;
   boostV = lambLor.BoostVector();
   boostV *= -1;
   prLor.Boost(boostV);
   vPr = prLor.Vect();

   cos_prot      = vPr.CosTheta(); // cos theta
   phi_prot_star = vPr.Phi();      // phi
}
