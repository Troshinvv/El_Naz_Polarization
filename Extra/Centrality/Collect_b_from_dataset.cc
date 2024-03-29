///Constructs a tree with impact parameter, multiplicity with/without using DCA cut
///Use the corresponding script to run the macro on the cluster
//Provide the 

// MPD includes: 
#include "TpcPoint.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdHelix.h"
#include "MpdVertex.h"
#include "MpdParticle.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdMotherFitterPart.h"
#include "MpdTrackFinderIts5spd.h"
#include "MpdMCEventHeader.h"
#include "MpdZdcDigi.h"
#include "MpdPid.h"
#include "MpdMCTrack.h"

// CBM includes:
#include "FairMCPoint.h"
#include "FairRunAna.h"

// ROOT includes:
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <Riostream.h>

#include <set>
#include <map>
#include <tuple>
#include <vector>

#include <iostream>

using namespace std;

#include "TStopwatch.h"

MpdVertex *mpdVert;
TVector3 vtxN, momN, primVert;
TClonesArray *itsTracks, *mcTracks, *mpdTracks;

//#define ITS
#ifdef ITS
typedef MpdItsKalmanTrack AzTrack;
#else
typedef MpdTpcKalmanTrack AzTrack;
#endif

//define the variables for the tree:
Float_t b0, Z_vertex, Z_vertex_Reco;
Int_t evNo, n_tracks_mpd, n_tracks_mpd_DCA;

//function to check input file:
bool checkfile (char* filename) 
{
	TFile *file_for_check = TFile::Open(filename);
	// if exist
    if (!file_for_check)
    {
        cout << "no specified file: " << filename << endl;
        return false;
    } 
    // if zombie
	if (file_for_check->IsZombie())
	{
        cout << "file: " << filename << " is a Zombie!" << endl;
        return false;
    }
	file_for_check->Close();
	return true;
}
//starting the main Analyzer:
int main(int argc, char** argv)
{
	TStopwatch timer;
    timer.Start();
    
	Int_t n1 = 0;
	Int_t n2 = 0;
	Int_t skipFiles = 0;
	Int_t iset = 1;
	
	Int_t _N_Hits = 16; // default value
	Double_t cut_pt = 0.15; // default value (GeV/c)
	Double_t cut_eta = 0.5; // default value
	Double_t dca_cut = 0.5; // default: 0.5 cm
	
	TChain chain("mpdsim");
	TString outname = ""; 
	
	for (int i = 0 ; i < argc ; i++) 
	{
		string tmp(argv[i]);
		if (tmp.compare("-events") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				n2 = number;
			}
		}
		if (tmp.compare("-cut_pt") == 0 && i+1 < argc) 
		{
			float number = atof(argv[i+1]);
			if (number > 0) 
			{
				cut_pt = number;
			}
		}
		if (tmp.compare("-cut_eta") == 0 && i+1 < argc) 
		{
			float number = atof(argv[i+1]);
			if (number > 0) 
			{
				cut_eta = number;
			}
		}
		if (tmp.compare("-dca_cut") == 0 && i+1 < argc) 
		{
			float number = atof(argv[i+1]);
			if (number > 0) 
			{
				dca_cut = number;
			}
		}
		if (tmp.compare("-_N_Hits") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				_N_Hits = number;
			}
		}		
		if (tmp.compare("-inname") == 0 && i+1 < argc) 
		{
			int count = 0;
			while(i+count+1 < argc && argv[i+count+1][0] != '-')
			{
				//add the file to the chain after checking:
				if (checkfile(argv[i+count+1]))
				{
					cout << "adding to the chain: " << argv[i+count+1] << endl;
					chain.Add(argv[i+count+1]);
				}
				count++;
			}
			cout << "number of files = " << count << endl;
		}
		if (tmp.compare("-outname") == 0 && i+1 < argc) 
		{
			outname = argv[i+1];
		}
	}
	
	if(outname.Length() == 0)
	{
		cerr << "!! Please provide the path for the output file !!" << endl;
		return 1;
	}
	
	cout << "events chosen = " << n2 << "; N_Hits = " << _N_Hits << "; cut_pt = " << cut_pt << "; cut_eta = " << cut_eta << "; dca_cut = " << dca_cut << endl;
	
	cout << "entries = " << chain.GetEntries() << endl;
	
	//return(0); //for debugging
	
	TChain *simITS = (TChain*) gROOT->FindObject("mpdsim");
	TFile fileITS(simITS->GetListOfFiles()->First()->GetTitle());
	TClonesArray *vtxs = (TClonesArray*) fileITS.FindObjectAny("Vertex");
	simITS->SetBranchAddress("Vertex",&vtxs);
	TBranch *vtxB = simITS->GetBranch("Vertex");
	
	itsTracks = (TClonesArray*) fileITS.FindObjectAny("TpcKalmanTrack");
	simITS->SetBranchAddress("TpcKalmanTrack",&itsTracks);
	TBranch *itsRecoB = simITS->GetBranch("TpcKalmanTrack");
	
	MpdEvent *event = 0x0;
	simITS->SetBranchAddress("MPDEvent.", &event);
	MpdMCEventHeader *mcHeader = 0x0;
	simITS->SetBranchAddress("MCEventHeader.", &mcHeader); 
	

	//output file:
	TFile out(outname,"recreate");
	
	TTree *tree = new TTree("event","Event");
	tree->Branch("b0",&b0,"b0/F"); //impact parameter
	tree->Branch("n_tracks_mpd",&n_tracks_mpd,"n_tracks_mpd/I"); // number of tracks (multiplicity) in TPC for centrality determination
	tree->Branch("n_tracks_mpd_DCA",&n_tracks_mpd_DCA,"n_tracks_mpd_DCA/I"); // number of tracks (multiplicity) in TPC for centrality determination, using DCA cut
	tree->Branch("Z_vertex",&Z_vertex,"Z_vertex/F"); //Z_vertex
	tree->Branch("Z_vertex_Reco",&Z_vertex_Reco,"Z_vertex_Reco/F"); //Z_vertex_Reco
	
	Int_t events = simITS->GetEntries();
	cout << " Number of events = " << events << " Number n2 = " << n2 << endl;
	if (n2 != 0) events = TMath::Min (events, n2);
	cout << " Number of events = " << events << endl;
	
	for (Int_t i = 0; i < events; ++i) 
	{
		if (i < n1) continue;
		simITS->GetEntry(i);
		evNo = i + 1;
		cout << " Event " << evNo << " out of " << events << endl;
		TVector3 genVert;
		mcHeader->GetVertex(genVert); 
    	
		Int_t nMpdTr = 0;
		if (event) mpdTracks = event->GetGlobalTracks();
		if (mpdTracks) nMpdTr = mpdTracks->GetEntriesFast();
		
		MpdVertex *vtx = (MpdVertex*) vtxs->First();
		mpdVert = vtx;
		Z_vertex = vtx->GetZ();
		
		//MpdVertex *vtx_Reco = (MpdVertex*) event->First();
		//mpdVert = vtx;
		Z_vertex_Reco = event->GetPrimaryVerticesZ();
		
		//cout << "Z_vertex = " << Z_vertex << "; Z_vertex_Reco = " << Z_vertex_Reco << endl;
	
		// Loop over DST tracks
		//calculating centrality in TPC (new calibration)
		Int_t k_check = 0;
		Int_t k_check_DCA = 0;
		for (Int_t j = 0; j < nMpdTr; ++j) 
		{
			MpdTrack *mpdTr = (MpdTrack*) mpdTracks->UncheckedAt(j);
			
			if (TMath::Abs(mpdTr->GetEta())<cut_eta && TMath::Abs(mpdTr->GetPt())>cut_pt && mpdTr->GetNofHits()>_N_Hits) 
			{
				k_check++;
				if (TMath::Sqrt(TMath::Power(mpdTr->GetDCAX(),2) + TMath::Power(mpdTr->GetDCAY(),2) + TMath::Power(mpdTr->GetDCAZ(),2)) < dca_cut) {k_check_DCA++;}
			}
			
		}
		n_tracks_mpd = k_check;
		n_tracks_mpd_DCA = k_check_DCA;
		
		b0 = mcHeader->GetB();	
		
		tree->Fill();
		
	} // for (Int_t i = 0; i < events;
	
	out.Write();
	out.Close();
	
	timer.Stop();
    Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
	cout << "End of the program" << endl;
	
	return(0);
}
