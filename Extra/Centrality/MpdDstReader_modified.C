///The MpdDstReader from the centrality framework with a couple of modifications:
//Extra histograms: hRefMult and hRefMult_DCA (to have both versions - with and without DCA cut in the same file), hRefMult_MC (using MCTracks for comparison), hPtReco and hPtMC (to compare p_T distributions)
//Absolute value of p_T used for the cut - necessary for GlobalTracks, which store negative values (see histogram hPtReco)

#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>

#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TRandom3.h>

// Define whether new or old mpdroot version is used
#define _NEW_MCSTACK_
//#define _OLD_MCSTACK_

#include <FairMCEventHeader.h>
#ifdef _OLD_MCSTACK_
#include <FairMCTrack.h>
#endif
#ifdef _NEW_MCSTACK_
#include <MpdMCTrack.h>
#endif
#include <MpdEvent.h>
#include <MpdZdcDigi.h>
#include <MpdPid.h>
#include <MpdTrack.h>
#include <MpdKalmanTrack.h>
#include <MpdVertex.h>

#include "commonFunctions.C"
void MpdDstReader_modified()
{
  TStopwatch timer;
  timer.Start();

  //Request 30: official production for the polarization analysis
  TString inname = "/eos/nica/mpd/sim/data/qa/req30v2/data/*.root", outname = "/scratch2/nazarova/CentralityFramework/Multiplicity_histos/MultHisto_Dataset2_PHSD_9GeV_abspt.root"; 
		
  TChain inChain("mpdsim");
  inChain.Add(inname);

  const Double_t cut_pt = 0.15; // default: 0.15 GeV/c
  const Double_t cut_eta = 0.5; // default: 0.5
  const Int_t cut_nhits = 16;   // default: 16
  const Double_t dca_cut = 0.5; // default: 0.5 cm

  TH1F *hRefMult = new TH1F("hRefMultSTAR","hRefMultSTAR",2500,0,2500);
  TH1F *hRefMult_DCA = new TH1F("hRefMultSTAR_DCA","hRefMultSTAR_DCA",2500,0,2500);
  TH2F *hBvsRefMult = new TH2F("hBvsRefMult","hBvsRefMult",2500,0,2500,200,0.,20.);
	TH1F *hRefMult_MC = new TH1F("hRefMultSTAR_MC","hRefMultSTAR_MC",2500,0,2500);
	
	TH1F *hPtReco = new TH1F("hPtReco","hPtReco",100,-5.,5.);
	TH1F *hPtMC = new TH1F("hPtMC","hPtMC",100,-5.,5.);
	
  FairMCEventHeader *MCHeader;
  TClonesArray *MCTracks;
  MpdEvent *MPDEvent;
  TClonesArray *MpdGlobalTracks;

  MCHeader = nullptr;
  MCTracks = nullptr;
  MPDEvent = nullptr;

  inChain.SetBranchAddress("MCEventHeader.", &MCHeader);
  inChain.SetBranchAddress("MCTrack", &MCTracks);
  inChain.SetBranchAddress("MPDEvent.", &MPDEvent);
  
  std::vector<Long64_t> vEvents;
  std::random_device rd;
  std::mt19937 g(rd());

  Long64_t Nentries = 5e5; // or (Long64_t) inChain->GetEntries();
  //Long64_t Nentries = (Long64_t) inChain.GetEntries();

  // Starting event loop
  Long64_t Nevents;
  Nevents = (Nentries < 5e5) ? Nentries : 5e5;

  Int_t refMult, refMult_DCA, refMult_MC, nhits;
  Double_t pt, eta;
  Long64_t iEvent;
  TRandom3 *rnd = new TRandom3();
  Double_t prob_skip;
  
  for (Long64_t jentry=0; jentry<Nevents;jentry++)
  {
    iEvent = jentry;

    inChain.GetEntry(iEvent);
    if (jentry%1000 == 0) std::cout << "Event ["
      << jentry << "/" << Nevents << "]" << std::endl;
    refMult = 0;
    refMult_DCA = 0;
    MpdGlobalTracks = (TClonesArray*) MPDEvent->GetGlobalTracks();
    Int_t ntracks = MpdGlobalTracks->GetEntriesFast();
    for (int iTr=0; iTr<ntracks; iTr++)
    {
      auto mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(iTr);
#ifdef _OLD_MCSTACK_
      auto mctrack = (FairMCTrack*) MCTracks->UncheckedAt(mpdtrack->GetID());
#endif
#ifdef _NEW_MCSTACK_
      auto mctrack = (MpdMCTrack*) MCTracks->UncheckedAt(mpdtrack->GetID());
#endif
    
      pt = mpdtrack->GetPt();
      eta = mpdtrack->GetEta();
      nhits = mpdtrack->GetNofHits();

	hPtReco->Fill(pt);

      if (TMath::Abs(pt) < cut_pt) continue; //absolute value of pt here (necessary for GlobalTracks)
      if (TMath::Abs(eta) > cut_eta) continue;
      if (nhits < cut_nhits) continue;
      
	refMult++; // Reco multiplicity without DCA cut
	
      // Primary track selection
      if (TMath::Sqrt(TMath::Power(mpdtrack->GetDCAX(),2) + TMath::Power(mpdtrack->GetDCAY(),2) + TMath::Power(mpdtrack->GetDCAZ(),2)) > dca_cut) continue;

      refMult_DCA++; // Reco multiplicity with DCA cut
    }
    Int_t nMC = MCTracks->GetEntriesFast();
    refMult_MC = 0;
    TVector3 mom; 
    for (Int_t j = 0; j < nMC; ++j) 
	{
		MpdMCTrack* mcTr = (MpdMCTrack*) MCTracks->UncheckedAt(j);
		if ((TMath::Abs(mcTr->GetPdgCode()) != 211) && (TMath::Abs(mcTr->GetPdgCode()) != 2212) && (TMath::Abs(mcTr->GetPdgCode()) != 321)) continue;
		Int_t mothId = mcTr->GetMotherId();
		if (mcTr->GetMotherId() != -1) continue; //only primary
		
		mcTr->GetMomentum(mom);
		Double_t pt_MC = mom.Pt();
		Double_t eta_MC = mom.Eta();
		
		hPtMC->Fill(pt_MC);
		
		if (TMath::Abs(pt_MC) < cut_pt) continue;
		if (TMath::Abs(eta_MC) > cut_eta) continue;
		refMult_MC++; // MC multiplicity (for primary tracks of protons, pions and kaons)
	}
    hRefMult->Fill(refMult);
    hRefMult_DCA->Fill(refMult_DCA);
    hRefMult_MC->Fill(refMult_MC);
    hBvsRefMult->Fill(refMult,MCHeader->GetB());
  }
  TFile *fo = new TFile(outname,"recreate");
  fo->cd();
  hRefMult->Write();
  hRefMult_DCA->Write();
  hRefMult_MC->Write();
  hBvsRefMult->Write();
  hPtReco->Write();
  hPtMC->Write();
  fo->Close();

  timer.Stop();
  timer.Print();
}
