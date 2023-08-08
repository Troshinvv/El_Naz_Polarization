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
   param("cent_cut", cent_cut, 70.0);
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

   ResEP1_true = new double[NITER_CENT];
   SubEvRes1   = new double[NITER_CENT];

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
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

      ResEP1_true[iter_cent] = 0.;
      SubEvRes1[iter_cent]   = 0.;
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
   hCentrality->Fill(Centrality_tpc);
   hNevCentr->Fill(Centrality_tpc);
   int nMC = mMCTracks->GetEntriesFast();

   for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++) {
      if (Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue;

      ResEP1_true[iter_cent] += ResEP;
      SubEvRes1[iter_cent] += ResEPSub;

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

               // Fill angular distribution for primary hyperons
               if (mcTr_Lam_MotherID < 0) {
                  hDeltaPhiRP_Prim[iter_cent]->Fill(phi_diff_RP);
                  hDeltaPhiEP_Prim[iter_cent]->Fill(phi_diff_EP);
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
