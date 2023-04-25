///Analyzer for Event plane test on the dataset
//Calculate for them the necessary distributions
//compile with: make MCTest_Lambda, start with the corresponding sge script
//e.g. ./MCTest_Lambda -events 10 -NITER_CENT 4 -NITER 20 -inname $INFILE -outname $OUTFILE -centname $CENTFILE
// MPD includes: 
#include "TpcPoint.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdHelix.h"
#include "MpdMCTrack.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdMCEventHeader.h"
#include "MpdZdcDigi.h"

// CBM includes:
#include "FairMCPoint.h"

// ROOT includes:
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TRandom2.h>
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

// Global parameters:
#define pi TMath::Pi()
#define _N_MODULES_TOTAL 90
#define _N_ARM 2 
const int CENT = 10; //amount of centrality intervals (for the centrality determination technique)
//pdg numbers 
Int_t pdgCodePr= 2212; // proton
Int_t pdgCodeAPr= -2212; // antiproton
Int_t pdgCodeNeg= -211; // pi-
Int_t pdgCodePos= 211; // pi+
Int_t pdgCodeL0 = 3122; // lambda (1.11568)
Int_t pdgCodeAL0 = -3122; // antilambda (1.11568)

Int_t pdgCodeNeutr= 2112; // neutron
Int_t pdgCodeANeutr= -2112; // antineutron

TClonesArray *itsTracks, *mcTracks, *mpdTracks;
TClonesArray *ZDCHits;
MpdZdcDigi* ZDCHit;
//define the functions:

Float_t GetCentrality(Int_t multiplicity);
Double_t GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID);
Double_t GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy);
Double_t GetPsiFullZdc(Float_t* zdc_energy, Int_t n);
void GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy);
Double_t* GetAngles();

Float_t b0;
Int_t evNo, n_tracks_mpd, n_tracks_mpd_DCA;
Float_t phiEP_1_ZDC, ResEP_1_ZDC, ResEPSub_1_ZDC;
Float_t ZDC_energy_mpd[_N_MODULES_TOTAL];	

TRandom *RNG = new TRandom();
Int_t Border_max[CENT];
Int_t Border_min[CENT];
Float_t Centrality_bin[CENT];

Float_t phiEP_mc;

double* init_double_array (const int n, const double fmt...)
{
	va_list args;
    va_start(args, fmt);
    
	double* ret = new double[n];
     
    for (int i=0 ; i<n ; i++) 
    {
		ret[i] = va_arg(args, double);
    }
 
	return ret;
}

int* init_int_array (const int n, const int fmt...)
{
	va_list args;
    va_start(args, fmt);
    
	int* ret = new int[n];
     
    for (int i=0 ; i<n ; i++) 
    {
		ret[i] = va_arg(args, int);
    }
 
	return ret;
}
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
	Int_t n1 = 0;
	Int_t n2 = 0;
	Int_t skipFiles = 0;
	Int_t iset = 1;
	
	Int_t _N_Hits = 16; // default value
	Double_t cut_pt = 0.15; // default value (GeV/c)
	Double_t cut_eta = 0.5; // default value
	Double_t dca_cut = 0.5; // default: 0.5 cm
	Int_t dca_choice = 0; // choice of DCA (0 - no DCA cut, 1 - with DCA cut)
	
	int NITER = 20; //amount of delta(phi) and cos(theta) cuts 
	int NITER_CENT = 4; //amount of centrality intervals (for the analysis)	
	
	TChain chain("mpdsim");
	TString outname = ""; 
	TString centname = "";
	
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
		if (tmp.compare("-dca_choice") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				dca_choice = number;
			}
		}	
		if (tmp.compare("-NITER") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				NITER = number;
			}
		}
		if (tmp.compare("-NITER_CENT") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				NITER_CENT = number;
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
		if (tmp.compare("-centname") == 0 && i+1 < argc) 
		{
			if (checkfile(argv[i+1]))
				centname = argv[i+1];
		}
	}
	
	if(outname.Length() == 0)
	{
		cerr << "!! Please provide the path for the output file !!" << endl;
		return 1;
	}
	if(centname.Length() == 0)
	{
		cerr << "!! Please provide the path for the centrality file !!" << endl;
		return 1;
	}
	cout << "events chosen = " << n2 << "; N_Hits = " << _N_Hits << "; cut_pt = " << cut_pt << "; cut_eta = " << cut_eta << "; dca_cut = " << dca_cut << "; dca_choice = " << dca_choice << endl;
	
	cout << "entries = " << chain.GetEntries() << endl;
	double SubEvRes1[NITER_CENT], ResEP1_true[NITER_CENT]; 
	
	int *centrality_min;
	int *centrality_max;
	double *_CentrBins;
	if (NITER_CENT == 4)
	{		
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		_CentrBins = init_double_array(5, 0, 0.,10.,20.,50.,100.);
	}else if (NITER_CENT == 7)
	{
		centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
		_CentrBins = init_double_array(8, 0, 0., 10., 20., 30., 40., 50., 60., 70.);
	}else if (NITER_CENT == 10)
	{
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		_CentrBins = init_double_array(11, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.);
	}else 
	{
		cout << "This values of centrality bins is not defined! Please provide the definition in the code." << endl;
		return 1;
	}
	
	Int_t Ncc, MinBorder, MaxBorder;
	Float_t MinPercent, MaxPercent;
	
	TFile *file_centrality = new TFile(centname);
	//reading the tree
	TTree *tree_centrality = (TTree *)file_centrality->Get("Result");    
	tree_centrality->SetBranchAddress("Ncc", &Ncc);
	tree_centrality->SetBranchAddress("MinPercent", &MinPercent);
	tree_centrality->SetBranchAddress("MaxPercent", &MaxPercent);
	tree_centrality->SetBranchAddress("MinBorder", &MinBorder);
	tree_centrality->SetBranchAddress("MaxBorder", &MaxBorder);
	cout << "Amount of events in centrality tree = " << tree_centrality->GetEntries() << endl;
	for (Int_t i = 0; i < tree_centrality->GetEntries(); ++i) 
	{
		tree_centrality->GetEntry(i);
		cout << "i = " << i << "; Ncc = " << Ncc << "; MinPercent = " << MinPercent << "; MaxPercent = " << MaxPercent << "; MinBorder = " << MinBorder << "; MaxBorder = " << MaxBorder << endl;
		Border_max[i] = MaxBorder;
		Border_min[i] = MinBorder;
		Centrality_bin[i] = MinPercent + 0.5*(MaxPercent - MinPercent);
	}
	
	for (Int_t i = 0; i < tree_centrality->GetEntries(); ++i)
	{	
		cout << "Border_min = " << Border_min[i] << "; Border_max = " << Border_max[i] << "; Centrality_bin = " << Centrality_bin[i] << endl;
	}
	
	TChain *simITS = (TChain*) gROOT->FindObject("mpdsim");
	TFile fileITS(simITS->GetListOfFiles()->First()->GetTitle());
	
	itsTracks = (TClonesArray*) fileITS.FindObjectAny("TpcKalmanTrack");
	simITS->SetBranchAddress("TpcKalmanTrack",&itsTracks);
	TBranch *itsRecoB = simITS->GetBranch("TpcKalmanTrack");
	
	MpdEvent *event = 0x0;
	simITS->SetBranchAddress("MPDEvent.", &event);
	MpdMCEventHeader *mcHeader = 0x0;
	simITS->SetBranchAddress("MCEventHeader.", &mcHeader); 
	
	mcTracks = (TClonesArray*) fileITS.FindObjectAny("MCTrack");
	simITS->SetBranchAddress("MCTrack",&mcTracks);
	TBranch *mcBranch = simITS->GetBranch("MCTrack");	
	
	ZDCHits = 0;
	simITS->SetBranchAddress("ZdcDigi",&ZDCHits);
	
	//output file:
	TFile out(outname,"recreate");
	
	TH1F *Psi_EP_diff_cent[NITER_CENT];
	char *Psi_EP_diff_cent_int = new char[NITER_CENT];
	
	TH1F *hEta_full = new TH1F("hEta_full","hEta_full",100,-7.,7.);
	TH1F *hEta_full_lowb = new TH1F("hEta_full_lowb","hEta_full_lowb",100,-7.,7.);
	TH1F *hEta_full_highb = new TH1F("hEta_full_highb","hEta_full_highb",100,-7.,7.);
	
	TH1F *hEta_pions_full = new TH1F("hEta_pions_full","hEta_pions_full",100,-7.,7.);
	TH1F *hEta_protons_full = new TH1F("hEta_protons_full","hEta_protons_full",100,-7.,7.);
	TH1F *hEta_neutrons_full = new TH1F("hEta_neutrons_full","hEta_neutrons_full",100,-7.,7.);
	
	TH1F *hEta_pions_lowb = new TH1F("hEta_pions_lowb","hEta_pions_lowb",100,-7.,7.);
	TH1F *hEta_protons_lowb = new TH1F("hEta_protons_lowb","hEta_protons_lowb",100,-7.,7.);
	TH1F *hEta_neutrons_lowb = new TH1F("hEta_neutrons_lowb","hEta_neutrons_lowb",100,-7.,7.);
	
	TH1F *hEta_pions_highb = new TH1F("hEta_pions_highb","hEta_pions_highb",100,-7.,7.);
	TH1F *hEta_protons_highb = new TH1F("hEta_protons_highb","hEta_protons_highb",100,-7.,7.);
	TH1F *hEta_neutrons_highb = new TH1F("hEta_neutrons_highb","hEta_neutrons_highb",100,-7.,7.);
	
	TH1F *hEta_pions_full_prim = new TH1F("hEta_pions_full_prim","hEta_pions_full_prim",100,-7.,7.);
	TH1F *hEta_protons_full_prim = new TH1F("hEta_protons_full_prim","hEta_protons_full_prim",100,-7.,7.);
	TH1F *hEta_neutrons_full_prim = new TH1F("hEta_neutrons_full_prim","hEta_neutrons_full_prim",100,-7.,7.);
	
	TH1F *hEta_pions_lowb_prim = new TH1F("hEta_pions_lowb_prim","hEta_pions_lowb_prim",100,-7.,7.);
	TH1F *hEta_protons_lowb_prim = new TH1F("hEta_protons_lowb_prim","hEta_protons_lowb_prim",100,-7.,7.);
	TH1F *hEta_neutrons_lowb_prim = new TH1F("hEta_neutrons_lowb_prim","hEta_neutrons_lowb_prim",100,-7.,7.);
	
	TH1F *hEta_pions_highb_prim = new TH1F("hEta_pions_highb_prim","hEta_pions_highb_prim",100,-7.,7.);
	TH1F *hEta_protons_highb_prim = new TH1F("hEta_protons_highb_prim","hEta_protons_highb_prim",100,-7.,7.);
	TH1F *hEta_neutrons_highb_prim = new TH1F("hEta_neutrons_highb_prim","hEta_neutrons_highb_prim",100,-7.,7.);
	
	TH1F *hEvent_counter = new TH1F("hEvent_counter","hEvent_counter",2,0.,2.);
	
	TH1F *NCentr = new TH1F("NCentr","NCentr",NITER_CENT,_CentrBins);
	NCentr->SetYTitle("Entries");
	NCentr->SetXTitle("Centrality TPC");
	NCentr->SetLineColor(kBlack);
	
	TH1D *Resolution_EP1_true = new TH1D("Resolution_EP1_true","Resolution_EP1_true",NITER_CENT,_CentrBins);
	Resolution_EP1_true->SetYTitle("R_{EP}^{1} = < cos(#Psi^{1}_{EP} - #Psi_{RP}) >");
	Resolution_EP1_true->SetXTitle("Centrality, [%]");
	Resolution_EP1_true->SetMarkerStyle(20);
	Resolution_EP1_true->SetLineColor(kBlack);
	Resolution_EP1_true->SetMarkerColor(kBlack);
	Resolution_EP1_true->SetMarkerSize(2);
	Resolution_EP1_true->SetLineWidth(2);
	
	TH1D *Resolution_EP1_exp = new TH1D("Resolution_EP1_exp","Resolution_EP1_exp",NITER_CENT,_CentrBins);
	Resolution_EP1_exp->SetYTitle("R_{EP}^{1} = < cos(#Psi^{1}_{EP} - #Psi_{RP}) >");
	Resolution_EP1_exp->SetXTitle("Centrality, [%]");
	Resolution_EP1_exp->SetMarkerStyle(20);
	Resolution_EP1_exp->SetLineColor(kBlue);
	Resolution_EP1_exp->SetMarkerColor(kBlue);
	Resolution_EP1_exp->SetMarkerSize(2);
	Resolution_EP1_exp->SetLineWidth(2);
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sprintf(Psi_EP_diff_cent_int,"Psi_EP_diff_cent_%d",iter_cent); 
		Psi_EP_diff_cent[iter_cent] = new TH1F(Psi_EP_diff_cent_int,Psi_EP_diff_cent_int,100,-600.,600.);
		Psi_EP_diff_cent[iter_cent]->SetYTitle("Entries");
		Psi_EP_diff_cent[iter_cent]->SetXTitle("#Psi^{1}_{EP} - #Psi_{RP} [deg]");
		Psi_EP_diff_cent[iter_cent]->SetLineColor(kBlack);
		Psi_EP_diff_cent[iter_cent]->SetLineWidth(2);
	}
	Int_t Centrality_tpc;
	
	Int_t events = simITS->GetEntries();
	cout << " Number of events = " << events << " Number n2 = " << n2 << endl;
	if (n2 != 0) events = TMath::Min (events, n2);
	cout << " Number of events = " << events << endl;
	
	Int_t ev_counter_lowb = 0;
	Int_t ev_counter_highb = 0;
	
	for (Int_t i = 0; i < events; ++i) 
	{
		if (i < n1) continue;
		simITS->GetEntry(i);
		evNo = i + 1;
		cout << " Event " << evNo << " out of " << events << endl;
		
		TVector3 genVert;
		mcHeader->GetVertex(genVert); 
		Int_t nMC = mcTracks->GetEntriesFast();
		cout << " nMC " << nMC << endl;
		
		Int_t nMpdTr = 0;
		if (event) mpdTracks = event->GetGlobalTracks();
		if (mpdTracks) nMpdTr = mpdTracks->GetEntriesFast();
		
		//Loop over DST tracks
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
		if(n_tracks_mpd == 0) continue;
		
		//Extracting impact parameter and RP angle:
		b0 = mcHeader->GetB();
			
		if(b0 < 6) ev_counter_lowb++;
		if(b0 > 6) ev_counter_highb++;
		
		for(Int_t j = 0; j < nMC; j++)
		{
			MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(j);
			Int_t pdg = mcTr->GetPdgCode();
			
			Float_t thetatrack = TMath::ATan2(mcTr->GetPt(),mcTr->GetPz());
			Float_t eta = -TMath::Log(TMath::Tan(thetatrack/2.));
			
			// Check production vertex
			TVector3 pos;
			Double_t r = 0.0;
			if (mcTr->GetMotherId() >= 0) 
			{
				mcTr->GetStartVertex(pos);
				pos -= genVert;
				r = pos.Mag();
				if (r >= 50.0) continue; 
			}
			
			hEta_full->Fill(eta);
			if(b0 < 6) hEta_full_lowb->Fill(eta);
			if(b0 > 6) hEta_full_highb->Fill(eta);
			switch (abs(pdg))
			{
				case 211: // pions
				hEta_pions_full->Fill(eta);
				if(b0 < 6) hEta_pions_lowb->Fill(eta);
				if(b0 > 6) hEta_pions_highb->Fill(eta);
				
				if((mcTr->GetMotherId() == -1)) hEta_pions_full_prim->Fill(eta);
				if((mcTr->GetMotherId() == -1) && (b0 < 6)) hEta_pions_lowb_prim->Fill(eta);
				if((mcTr->GetMotherId() == -1) && (b0 > 6)) hEta_pions_highb_prim->Fill(eta);
				
				break;
				case 2112: // neutron
				hEta_neutrons_full->Fill(eta);
				if(b0 < 6) hEta_neutrons_lowb->Fill(eta);
				if(b0 > 6) hEta_neutrons_highb->Fill(eta);
				
				if((mcTr->GetMotherId() == -1)) hEta_neutrons_full_prim->Fill(eta);
				if((mcTr->GetMotherId() == -1) && (b0 < 6)) hEta_neutrons_lowb_prim->Fill(eta);
				if((mcTr->GetMotherId() == -1) && (b0 > 6)) hEta_neutrons_highb_prim->Fill(eta);
				
				break;
				case 2212: // protons
				hEta_protons_full->Fill(eta);
				if(b0 < 6) hEta_protons_lowb->Fill(eta);
				if(b0 > 6) hEta_protons_highb->Fill(eta);
				
				if((mcTr->GetMotherId() == -1)) hEta_protons_full_prim->Fill(eta);
				if((mcTr->GetMotherId() == -1) && (b0 < 6)) hEta_protons_lowb_prim->Fill(eta);
				if((mcTr->GetMotherId() == -1) && (b0 > 6)) hEta_protons_highb_prim->Fill(eta);
				
				break;
				default: break;
			}
		}
		
		// calculating EP angles:
		phiEP_mc = mcHeader->GetRotZ();
		phiEP_mc = TMath::ATan2(TMath::Sin(phiEP_mc),TMath::Cos(phiEP_mc));
		
		// calculating ZDC energy loss:
		for (long int j = 0; j < 90; ++j)
		{
			ZDC_energy_mpd[j] = 0;
		}
		Int_t number_of_zdchits = ZDCHits->GetEntries();
		for (Int_t zdchit_number = 0; zdchit_number < number_of_zdchits; ++zdchit_number)
		{
			ZDCHit = (MpdZdcDigi*) ZDCHits->At(zdchit_number);
			Int_t detector_ID = ZDCHit->GetDetectorID();//1,2
			Int_t module_ID = ZDCHit->GetModuleID();//1-45
			Double_t energy_deposit_per_hit = ZDCHit->GetELoss();
				
			ZDC_energy_mpd[(detector_ID - 1)*45 + module_ID - 1] += energy_deposit_per_hit;
		}
		//if both ZDCs have some signal then calculate EP and its resolutions using ZDC
		if ((GetTotalEnergy(ZDC_energy_mpd, 0) != 0)&&(GetTotalEnergy(ZDC_energy_mpd, 1) != 0))
		{
			phiEP_1_ZDC = GetPsiFullZdc(ZDC_energy_mpd, 1);
			ResEP_1_ZDC = TMath::Cos(phiEP_1_ZDC - phiEP_mc);
				
			Double_t qx,qy;
			Double_t psi_N_R_1 = GetPsiHalfZdc(ZDC_energy_mpd, 0, 1, qx, qy);
			Double_t psi_N_L_1 = GetPsiHalfZdc(ZDC_energy_mpd, 1, 1, qx, qy);
			ResEPSub_1_ZDC = TMath::Cos(psi_N_R_1 - psi_N_L_1);
		}
		
		//Estimate centrality from the TPC multiplicity
		if (dca_choice == 0)
		{
			Centrality_tpc = GetCentrality(n_tracks_mpd);
		}else if(dca_choice == 1)
		{
			Centrality_tpc = GetCentrality(n_tracks_mpd_DCA);
		}else
		{
			cout << "This value of DCA choice is not defined!" << endl;
			return 1;
		}
		NCentr->Fill(Centrality_tpc);	
		
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			if(Centrality_tpc < centrality_max[iter_cent] && Centrality_tpc >= centrality_min[iter_cent])
			{
				ResEP1_true[iter_cent] += ResEP_1_ZDC;
				SubEvRes1[iter_cent] += ResEPSub_1_ZDC;
				Psi_EP_diff_cent[iter_cent]->Fill(180.*(phiEP_1_ZDC - phiEP_mc)/3.14);
			}
		}
	} // for (Int_t i = 0; i < events;
	
	cout << "ev_counter_lowb = " << ev_counter_lowb << endl;
	cout << "ev_counter_highb = " << ev_counter_highb << endl;
	hEvent_counter->SetBinContent(1,ev_counter_lowb);
	hEvent_counter->SetBinContent(2,ev_counter_highb);
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);
		Resolution_EP1_exp->SetBinContent(iter_cent+1,SubEvRes1[iter_cent]);
	}

	out.Write();
	out.Close();
	
	return(0);
}
//get centrality from TPC (through track multiplicity)
Float_t GetCentrality(Int_t multiplicity)
{
	Float_t centrality = -1.;
	
	for (int i = 0; i < CENT; ++i)
	{
		if(i == 0)
		{
			if(multiplicity >= Border_min[i])
			{
				centrality = Centrality_bin[i];
				Double_t random_number_new = RNG->Uniform(Centrality_bin[i]-5,Centrality_bin[i]+5);
				centrality = random_number_new;
			}
		}else if(i == CENT-1)
		{
			if(multiplicity >= Border_min[i] && multiplicity <= Border_max[i])
			{
				centrality = Centrality_bin[i];
				Double_t random_number_new = RNG->Uniform(Centrality_bin[i]-5,Centrality_bin[i]+5);
				centrality = random_number_new;
			}
		}else 
		if(multiplicity >= Border_min[i] && multiplicity < Border_max[i])
		{
			centrality = Centrality_bin[i];
			Double_t random_number_new = RNG->Uniform(Centrality_bin[i]-5,Centrality_bin[i]+5);
			centrality = random_number_new;
		}
	}
	return centrality;
}
Double_t GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy)
{
	Double_t Qcos=0., Qsin=0.;
	GetQsZdc(zdc_energy,zdc_ID,n,Qcos,Qsin);

	qx = Qcos;
	qy = Qsin;
	Double_t PsiEP = (1/(Double_t)n) * TMath::ATan2(Qsin,Qcos);
	return PsiEP;
}

//
Double_t GetPsiFullZdc(Float_t* zdc_energy, Int_t n)
{
	Double_t QcosR=0., QsinR=0.;
	GetQsZdc(zdc_energy,0,n,QcosR, QsinR);

	Double_t QcosL=0., QsinL=0.;
	GetQsZdc(zdc_energy,1,n,QcosL,QsinL);

	Double_t psiEP = TMath::ATan2(QsinR + QsinL,QcosR + QcosL)/(Double_t) n; // (-pi/n,pi/n]
    return psiEP;
}
//
void GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy)
{
	Double_t *phi_angle_of_modules = GetAngles();
	Double_t Qcos = 0 , Qsin = 0;
	Double_t w_sign = 0;
	Double_t total_energy = GetTotalEnergy(zdc_energy,zdc_ID);

	if (harm == 2) w_sign = 1;
	else if (harm == 1)
	{
		if (zdc_ID == 0) w_sign = 1;
		else if (zdc_ID == 1) w_sign = -1;
	}


	for (int module = _N_MODULES_TOTAL/2 * zdc_ID; module < _N_MODULES_TOTAL/2 * (zdc_ID + 1); ++module)
	{
		if ((module==22) || (module==67)) continue;
//		if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
//				(i==66)||(i==68)||(i==74)) continue;
		Qcos += w_sign*zdc_energy[module] / total_energy * TMath::Cos(harm*phi_angle_of_modules[module]);
        Qsin += w_sign*zdc_energy[module] / total_energy * TMath::Sin(harm*phi_angle_of_modules[module]);
	}

	Qx = Qcos; Qy = Qsin;
}
//
Double_t GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID)
{
	Double_t total_energy = 0.;
	for (int i = _N_MODULES_TOTAL/2 * zdc_ID; i < _N_MODULES_TOTAL/2 * (zdc_ID + 1); ++i)
	{
		if ((i==22) || (i==67)) continue;
//		if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
//				(i==66)||(i==68)||(i==74)) continue;
		total_energy += zdc_energy[i];
	}
	return total_energy;
}
//
Double_t* GetAngles()
{
	Double_t *phi_angle_of_module = new Double_t[_N_MODULES_TOTAL];

	for (int i = 0; i < _N_ARM; ++i)
	{
		Int_t x_axis_switch;
		if (i == 0) x_axis_switch = 1;
		else if (i == 1) x_axis_switch = -1;

		for (Int_t j = 0; j < _N_MODULES_TOTAL/2; ++j)
		{
			Double_t y = 0, x = 0;

			if ((j>=0) && (j<=4))
			{
				y = 45., x = (j-2)*15.;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = TMath::ATan2(y,x_axis_switch*x);
			}
			else if ((j>=5) && (j<=39))
			{
				y = (3-(j+2)/7)*15, x = (3-(j+2)%7)*15;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = TMath::ATan2(y,x_axis_switch*x);
			}
			else if ((j>=40) && (j<=44))
			{
				y = -45. , x = (j-42)*15.;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = TMath::ATan2(y,x_axis_switch*x);
			}
		}
	}

	return phi_angle_of_module;
}
