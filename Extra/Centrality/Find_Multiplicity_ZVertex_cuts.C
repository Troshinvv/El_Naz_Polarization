//Creates different plots for the completed centrality calibration
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <Riostream.h>
#include <TLatex.h>

#include <iostream>
using namespace std;

void Find_Multiplicity_ZVertex_cuts(TString inname = "/scratch2/nazarova/polarization-analysis-framework/output_vertexstudy_default/*.root",TString outname = "/scratch2/nazarova/CentralityFramework/Multiplicity_histos/MultHisto_Vertexcuts.root")
{
	TChain chain("event");
	chain.Add(inname);
	
	//variables:
	Float_t b0, Z_vertex, Z_vertex_Reco;
	Int_t n_tracks_mpd, n_tracks_mpd_DCA_default, n_tracks_mpd_DCA_liza, n_tracks_mpd_DCA_victor;
	
	//set adress to the branches of the input tree:
	chain.SetBranchAddress("b0", &b0);
	chain.SetBranchAddress("n_tracks_mpd", &n_tracks_mpd);
	chain.SetBranchAddress("n_tracks_mpd_DCA_default", &n_tracks_mpd_DCA_default);
	chain.SetBranchAddress("n_tracks_mpd_DCA_liza", &n_tracks_mpd_DCA_liza);
	chain.SetBranchAddress("n_tracks_mpd_DCA_victor", &n_tracks_mpd_DCA_victor);
	chain.SetBranchAddress("Z_vertex", &Z_vertex);
	chain.SetBranchAddress("Z_vertex_Reco", &Z_vertex_Reco);
	
	TH1F *hRefMult_DCA_default = new TH1F("hRefMult_DCA_default","hRefMult_DCA_default",2500,0,2500);
	TH1F *hRefMult_DCA_default_vertex_50 = new TH1F("hRefMult_DCA_default_vertex_50","hRefMult_DCA_default_vertex_50",2500,0,2500);
	TH1F *hRefMult_DCA_default_vertex_10 = new TH1F("hRefMult_DCA_default_vertex_10","hRefMult_DCA_default_vertex_10",2500,0,2500);
	
	TH1F *hRefMult_DCA_liza = new TH1F("hRefMult_DCA_liza","hRefMult_DCA_liza",2500,0,2500);
	TH1F *hRefMult_DCA_liza_vertex_50 = new TH1F("hRefMult_DCA_liza_vertex_50","hRefMult_DCA_liza_vertex_50",2500,0,2500);
	TH1F *hRefMult_DCA_liza_vertex_10 = new TH1F("hRefMult_DCA_liza_vertex_10","hRefMult_DCA_liza_vertex_10",2500,0,2500);
	
	TH1F *hRefMult_DCA_victor = new TH1F("hRefMult_DCA_victor","hRefMult_DCA_victor",2500,0,2500);
	TH1F *hRefMult_DCA_victor_vertex_50 = new TH1F("hRefMult_DCA_victor_vertex_50","hRefMult_DCA_victor_vertex_50",2500,0,2500);
	TH1F *hRefMult_DCA_victor_vertex_10 = new TH1F("hRefMult_DCA_victor_vertex_10","hRefMult_DCA_victor_vertex_10",2500,0,2500);
	
	Int_t full_events = chain.GetEntries();
	cout << " Number of events = " << full_events << endl;
	Int_t fill_counter = 0;
	for(Int_t iev = 0; iev < full_events; ++iev) //cycle for entries
	{   	
		chain.GetEntry(iev); // read event
		//chain.GetEntry(3*iev); // read event
		if (iev%1000 == 0) std::cout << "Event [" << iev << "/" << full_events << "]" << std::endl;
		/*
		if(fill_counter < 500000)
		{
			hRefMult_DCA->Fill(n_tracks_mpd_DCA);	
		}
		fill_counter++;
		*/
		
		hRefMult_DCA_default->Fill(n_tracks_mpd_DCA_default);
		hRefMult_DCA_liza->Fill(n_tracks_mpd_DCA_liza);
		hRefMult_DCA_victor->Fill(n_tracks_mpd_DCA_victor);
		
		if(Z_vertex == 0) continue;
		
		if(TMath::Abs(Z_vertex) < 50) 
		{
			hRefMult_DCA_default_vertex_50->Fill(n_tracks_mpd_DCA_default);
			hRefMult_DCA_liza_vertex_50->Fill(n_tracks_mpd_DCA_liza);
			hRefMult_DCA_victor_vertex_50->Fill(n_tracks_mpd_DCA_victor);
		}
		
		if(TMath::Abs(Z_vertex) < 10) 
		{
			hRefMult_DCA_default_vertex_10->Fill(n_tracks_mpd_DCA_default);
			hRefMult_DCA_liza_vertex_10->Fill(n_tracks_mpd_DCA_liza);
			hRefMult_DCA_victor_vertex_10->Fill(n_tracks_mpd_DCA_victor);
		}
	}
	
	TFile *fo = new TFile(outname,"recreate");
	fo->cd();
	hRefMult_DCA_default->Write();
	hRefMult_DCA_liza->Write();
	hRefMult_DCA_victor->Write();
	
	hRefMult_DCA_default_vertex_50->Write();
	hRefMult_DCA_liza_vertex_50->Write();
	hRefMult_DCA_victor_vertex_50->Write();
	
	hRefMult_DCA_default_vertex_10->Write();
	hRefMult_DCA_liza_vertex_10->Write();
	hRefMult_DCA_victor_vertex_10->Write();
}
