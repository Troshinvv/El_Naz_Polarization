///Analyzer for Lambda MCtest (via MCTracks info) on the dataset
//Choose protons -> then only the ones who come from Lambda decays 
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

TClonesArray *itsTracks, *mcTracks, *mpdTracks;
TClonesArray *ZDCHits;
MpdZdcDigi* ZDCHit;
//define the functions:

Float_t GetCentrality(Int_t multiplicity);
void FindPolarAngle(TVector3 &vPr, TVector3 &vLamb, Double_t &cos_prot, Double_t &phi_prot_star);
void Calculate_transverse_cos(TVector3 &vPr, TVector3 &vLamb, Double_t &cos_trans);

Double_t GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID);
Double_t GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy);
Double_t GetPsiFullZdc(Float_t* zdc_energy, Int_t n);
void GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy);
Double_t* GetAngles();

//define the variables we need:
Double_t costh, phi, cos_prot, phi_prot_star, phi_Lam, Weight_pol, cos_trans;
Float_t b0;
Int_t evNo, n_tracks_mpd, n_tracks_mpd_DCA;
Float_t phiEP_1_ZDC, ResEP_1_ZDC, ResEPSub_1_ZDC;
Float_t ZDC_energy_mpd[_N_MODULES_TOTAL];	
	
TRandom *RNG = new TRandom();
Int_t Border_max[CENT];
Int_t Border_min[CENT];
Float_t Centrality_bin[CENT];

Float_t phiEP_mc;
TVector3 polar3_new;

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
	Double_t cut_eta = 1.0; // default value
	Double_t dca_cut = 1.0; // default: 0.5 cm
	Int_t dca_choice = 0; // choice of DCA (0 - no DCA cut, 1 - with DCA cut)
	double cent_cut = 70.; // value for centrality cut
	Int_t cent_cut_choice = 0; // choice of centrality cut (0 - no cent cut, 1 - with cent cut)
	
	int NITER = 20; //amount of delta(phi) and cos(theta) cuts 
	int NITER_CENT = 4; //amount of centrality intervals (for the analysis)	
	
	TChain chain("mpdsim");
	TString outname = ""; 
	TString centname = "";
	TString particle_choice = "Lambda"; 
	TString angle_choice = "RP"; 
	
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
		if (tmp.compare("-angle_choice") == 0 && i+1 < argc) 
		{
			angle_choice = argv[i+1];
			cout << "angle_choice: " << argv[i+1] << endl;
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
		if (tmp.compare("-particle_choice") == 0 && i+1 < argc) 
		{
			particle_choice = argv[i+1];
			cout << "particle: " << argv[i+1] << endl;
		}
		if (tmp.compare("-CENT_CUT") == 0 && i+1 < argc) 
		{
			double number = atof(argv[i+1]);
			if (number > 0) 
			{
				cent_cut = number;
			}
		}
		if (tmp.compare("-cent_cut_choice") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				cent_cut_choice = number;
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
	cout << "events chosen = " << n2 << "; N_Hits = " << _N_Hits << "; cut_pt = " << cut_pt << "; cut_eta = " << cut_eta << "; dca_cut = " << dca_cut << "; dca_choice = " << dca_choice << "; particle_choice = " << particle_choice << "; angle_choice = " << angle_choice << endl;
	
	cout << "entries full = " << chain.GetEntries() << endl;
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
	
	//histograms for the output:
	
	TH1D *Lpolar_x[NITER_CENT], *Lpolar_y[NITER_CENT], *Lpolar_z[NITER_CENT], *Lpolar_norm[NITER_CENT], *Lpolar_y_prim[NITER_CENT];
	char *lpolar_x_int = new char[NITER_CENT];
	char *lpolar_y_int = new char[NITER_CENT];
	char *lpolar_z_int = new char[NITER_CENT];
	char *lpolar_norm_int = new char[NITER_CENT];	
	char *lpolar_y_prim_int = new char[NITER_CENT];
	
	TH1D *CosTheta_hist[NITER_CENT], *Pstar_hist[NITER_CENT], *PstarRP_hist[NITER_CENT];
	char *CosTheta_hist_int = new char[NITER_CENT];
	char *Pstar_hist_int = new char[NITER_CENT];
	char *PstarRP_hist_int = new char[NITER_CENT];
	
	TH1D *CosTheta_hist_prim[NITER_CENT], *Pstar_hist_prim[NITER_CENT], *PstarRP_hist_prim[NITER_CENT];
	char *CosTheta_hist_prim_int = new char[NITER_CENT];
	char *Pstar_hist_prim_int = new char[NITER_CENT];
	char *PstarRP_hist_prim_int = new char[NITER_CENT];
	
	TH1D *CosTheta_trans[NITER_CENT];
	char *CosTheta_trans_int = new char[NITER_CENT];
	
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
	
	TH1F *NPositive = new TH1F("Positive Polarization","Positive Polarization",NITER_CENT,_CentrBins);
	NPositive->SetYTitle("N_{H_pos}");
	NPositive->SetXTitle("Centrality TPC");
	NPositive->SetLineColor(kBlack);
	TH1F *NNegative = new TH1F("Negative Polarization","Negative Polarization",NITER_CENT,_CentrBins);
	NNegative->SetYTitle("N_{H_neg}");
	NNegative->SetXTitle("Centrality TPC");
	NNegative->SetLineColor(kBlack);
	TH1F *NHyperon = new TH1F("Amount of particles","Amount of particles",NITER_CENT,_CentrBins);
	NHyperon->SetYTitle("N_{H}");
	NHyperon->SetXTitle("Centrality TPC");
	NHyperon->SetLineColor(kBlack);
	Float_t Number_Positive[NITER_CENT];
	Float_t Number_Negative[NITER_CENT];
	Float_t Number_Lambda[NITER_CENT];
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Number_Positive[iter_cent] = 0;
		Number_Negative[iter_cent] = 0;
		Number_Lambda[iter_cent] = 0;
		//Distributions of components of polarization vector for full hyperons (primary+secondary) and P_{y} component for primary hyperons
		sprintf(lpolar_x_int,"LPolar_x_%d",iter_cent);      
		Lpolar_x[iter_cent] = new TH1D(lpolar_x_int,lpolar_x_int,100,-1.,1.);
		sprintf(lpolar_y_int,"LPolar_y_%d",iter_cent);      
		Lpolar_y[iter_cent] = new TH1D(lpolar_y_int,lpolar_y_int,100,-1.,1.);
		sprintf(lpolar_z_int,"LPolar_z_%d",iter_cent);      
		Lpolar_z[iter_cent] = new TH1D(lpolar_z_int,lpolar_z_int,100,-1.,1.);
		sprintf(lpolar_norm_int,"LPolar_norm_%d",iter_cent);      
		Lpolar_norm[iter_cent] = new TH1D(lpolar_norm_int,lpolar_norm_int,100,-1.,1.);
		
		sprintf(lpolar_y_prim_int,"LPolar_y_prim_%d",iter_cent);      
		Lpolar_y_prim[iter_cent] = new TH1D(lpolar_y_prim_int,lpolar_y_prim_int,100,-1.,1.);
		
		//Angular distributions for full hyperons from MCTracks
		sprintf(CosTheta_hist_int,"CosTheta_hist_%d",iter_cent);      
		CosTheta_hist[iter_cent] = new TH1D(CosTheta_hist_int,CosTheta_hist_int,NITER,-1.0,1.0);
		sprintf(Pstar_hist_int,"Pstar_hist_%d",iter_cent);      
		Pstar_hist[iter_cent] = new TH1D(Pstar_hist_int,Pstar_hist_int,NITER,0.,2.*pi);
		sprintf(PstarRP_hist_int,"PstarRP_hist_%d",iter_cent);      
		PstarRP_hist[iter_cent] = new TH1D(PstarRP_hist_int,PstarRP_hist_int,NITER,0.,2.*pi);
		
		//Angular distributions for primary hyperons from MCTracks
		sprintf(CosTheta_hist_prim_int,"CosTheta_hist_prim_%d",iter_cent);      
		CosTheta_hist_prim[iter_cent] = new TH1D(CosTheta_hist_prim_int,CosTheta_hist_prim_int,NITER,-1.0,1.0);
		sprintf(Pstar_hist_prim_int,"Pstar_hist_prim_%d",iter_cent);      
		Pstar_hist_prim[iter_cent] = new TH1D(Pstar_hist_prim_int,Pstar_hist_prim_int,NITER,0.,2.*pi);
		sprintf(PstarRP_hist_prim_int,"PstarRP_hist_prim_%d",iter_cent);      
		PstarRP_hist_prim[iter_cent] = new TH1D(PstarRP_hist_prim_int,PstarRP_hist_prim_int,NITER,0.,2.*pi);
		
		//Angular distribution for transverse polarization
		sprintf(CosTheta_trans_int,"CosTheta_trans_%d",iter_cent);      
		CosTheta_trans[iter_cent] = new TH1D(CosTheta_trans_int,CosTheta_trans_int,NITER,-1.0,1.0);
		
	}
	
	Int_t Centrality_tpc;
	
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
		
		//Extracting impact parameter and RP angle:
		b0 = mcHeader->GetB();
		if (angle_choice == "RP")
		{
			phiEP_mc = mcHeader->GetRotZ();
			phiEP_mc = TMath::ATan2(TMath::Sin(phiEP_mc),TMath::Cos(phiEP_mc));
		}else if (angle_choice == "EP")
		{
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
		}else
		{
			cout << "This angle choice is not defined! Please choose either RP or EP" << endl;
			return 1;
		}
		//Estimate centrality from the TPC multiplicity
		if (dca_choice == 0)
		{
			if(n_tracks_mpd == 0) continue;
			Centrality_tpc = GetCentrality(n_tracks_mpd);
		}else if(dca_choice == 1)
		{
			if(n_tracks_mpd_DCA == 0) continue;
			Centrality_tpc = GetCentrality(n_tracks_mpd_DCA);
		}else
		{
			cout << "This value of DCA choice is not defined!" << endl;
			return 1;
		}
		if (cent_cut_choice == 0)
		{
		}else if (cent_cut_choice == 1)
		{
			if(Centrality_tpc > cent_cut) continue;
		}else
		{
			cout << "This value of cent_cut_choice is not defined!" << endl;
			return 1;
		}
		NCentr->Fill(Centrality_tpc);
		
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			if(Centrality_tpc < centrality_max[iter_cent] && Centrality_tpc >= centrality_min[iter_cent])
			{
				ResEP1_true[iter_cent] += ResEP_1_ZDC;
				SubEvRes1[iter_cent] += ResEPSub_1_ZDC;
				for (Int_t j = 0; j < nMC; ++j) 
				{
					TVector3 mom; 
					TVector3 mom_moth; 
					MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(j);
					mcTr->GetMomentum(mom);		
					
					//introducing cut on pt and eta of proton to check:
					//if (TMath::Abs(mom.Eta())>1.3) continue;	
					
					if (particle_choice == "Lambda")
					{	
						if (mcTr->GetPdgCode() == pdgCodePr) 
						{
							int mcTr_MotherID = mcTr->GetMotherId();
							if (mcTr_MotherID < 0) continue; // if it's a primary proton, we don't need it													
							MpdMCTrack* mcTr_Mother = (MpdMCTrack*) mcTracks->UncheckedAt(mcTr_MotherID); // get the track, corresponding to the ID of protons' mother
							if(mcTr_Mother->GetPdgCode() == pdgCodeL0) // choose only the case, when Lambda is the mother
							{
								mcTr_Mother->GetMomentum(mom_moth);
								
								//introducing cuts on pt and eta of lambda to check:
								//cout << "mom.Eta() = " << mom_moth.Eta() << endl;
								//cout << "mom_moth.Eta() = " << mom_moth.Eta() << "; " << mom_moth.Pt() << endl;
								
								//if (TMath::Abs(mom_moth.Eta())>=1. && mom_moth.Pt()>=3 && mom_moth.Pt()<=0.4) continue;
								
								int mcTr_Lam_MotherID = mcTr_Mother->GetMotherId(); // ID of the mother of Lambda
								Float_t polar_x = 0.;
								Float_t polar_y = 0.;
								Float_t polar_z = 0.;
								Float_t pol_length = 0.;
								Int_t sign = 1;
								
								Float_t weight_pol = mcTr_Mother->GetWeight();
								polar_x = mcTr_Mother->GetPolar(0);
								polar_y = mcTr_Mother->GetPolar(1);
								polar_z = mcTr_Mother->GetPolar(2);
								sign = TMath::Sign(1, mcTr_Mother->GetPdgCode());
									
								TVector3 polar_changed(polar_x, polar_y, polar_z); //polarization vector for distributions		
								if (phiEP_mc != 0.) polar_changed.RotateZ(-phiEP_mc);
								polar_x = weight_pol*polar_changed.X();
								polar_y = weight_pol*polar_changed.Y();
								polar_z = weight_pol*polar_changed.Z();
								pol_length = TMath::Sqrt(polar_x*polar_x + polar_y*polar_y + polar_z*polar_z);
								//for primary Lambda
								if (mcTr_Lam_MotherID < 0) 
								{
									Lpolar_y_prim[iter_cent]->Fill(polar_y);
									Lpolar_y_prim[iter_cent]->SetYTitle("Entries");
									Lpolar_y_prim[iter_cent]->SetXTitle("P_{y}");
									Lpolar_y_prim[iter_cent]->SetLineColor(kBlack);	
								}
								Lpolar_x[iter_cent]->Fill(polar_x);
								Lpolar_x[iter_cent]->SetYTitle("Entries");
								Lpolar_x[iter_cent]->SetXTitle("P_{x}");
								Lpolar_x[iter_cent]->SetLineColor(kBlack);	
									
								Lpolar_y[iter_cent]->Fill(polar_y);
								Lpolar_y[iter_cent]->SetYTitle("Entries");
								Lpolar_y[iter_cent]->SetXTitle("P_{y}");
								Lpolar_y[iter_cent]->SetLineColor(kBlack);	
									
								Lpolar_z[iter_cent]->Fill(polar_z);
								Lpolar_z[iter_cent]->SetYTitle("Entries");
								Lpolar_z[iter_cent]->SetXTitle("P_{z}");
								Lpolar_z[iter_cent]->SetLineColor(kBlack);	
									
								Lpolar_norm[iter_cent]->Fill(pol_length);
								Lpolar_norm[iter_cent]->SetYTitle("Entries");
								Lpolar_norm[iter_cent]->SetXTitle("|P|");
								Lpolar_norm[iter_cent]->SetLineColor(kBlack);
									
								//proton values (MagThetaPhi)
								Float_t p_prot = mom.Mag();
								Float_t eta_prot = mom.Theta();
								Float_t phi_prot = mom.Phi();
								//lambda values (MagThetaPhi)	
								Float_t p_lam = mom_moth.Mag();
								Float_t eta_lam = mom_moth.Theta();
								Float_t phi_lam = mom_moth.Phi();
								TVector3 vPr, vLamb;
								vPr.SetMagThetaPhi(p_prot, eta_prot, phi_prot);
								vLamb.SetMagThetaPhi(p_lam, eta_lam, phi_lam);
									
								//calculate the costheta and phi from dataset:
								cos_prot = 0.0; 
								phi_prot_star = 0.0; 
								FindPolarAngle (vPr, vLamb, cos_prot, phi_prot_star);
								float phi_star = phi_prot_star;
								float phi_diff = 0.;
								if (angle_choice == "RP")
								{
									phi_diff = phiEP_mc - phi_prot_star;
								}else if (angle_choice == "EP")
								{
									phi_diff = phiEP_1_ZDC - phi_prot_star;		
								}else
								{
									cout << "This angle choice is not defined! Please choose either RP or EP" << endl;
									return 1;
								}
								
								if (phi_star < 0) phi_star = phi_star + 2.*pi;
								if (phi_diff < 0) phi_diff = phi_diff + 2.*pi;
								
								//fill distribution for all Lambda	
								CosTheta_hist[iter_cent]->Fill(cos_prot);
								Pstar_hist[iter_cent]->Fill(phi_star);
								PstarRP_hist[iter_cent]->Fill(phi_diff);
								
								//fill distribution for primary Lambda
								if (mcTr_Lam_MotherID < 0) 
								{
									CosTheta_hist_prim[iter_cent]->Fill(cos_prot);
									Pstar_hist_prim[iter_cent]->Fill(phi_star);
									PstarRP_hist_prim[iter_cent]->Fill(phi_diff);
								}
								//calculate the costheta for transverse polarization:
								cos_trans = 0.0; 
								Calculate_transverse_cos(vPr, vLamb, cos_trans);
								CosTheta_trans[iter_cent]->Fill(cos_trans);
									
								//number of positively/negatively polarized Lambdas:
								if(polar_changed.Y() > 0) Number_Positive[iter_cent]++;
								if(polar_changed.Y() < 0) Number_Negative[iter_cent]++;							
								//number of found Lambdas:
								Number_Lambda[iter_cent]++;									
							} // mcTr_Mother->GetPdgCode() == pdgCodeL0
								
						} // mcTr->GetPdgCode() == pdgCodePr	
					}else if (particle_choice == "ALambda")
					{
						if (mcTr->GetPdgCode() == pdgCodeAPr) 
						{
							int mcTr_MotherID = mcTr->GetMotherId();
							if (mcTr_MotherID < 0) continue; // if it's a primary proton, we don't need it													
							MpdMCTrack* mcTr_Mother = (MpdMCTrack*) mcTracks->UncheckedAt(mcTr_MotherID); // get the track, corresponding to the ID of protons' mother
							if(mcTr_Mother->GetPdgCode() == pdgCodeAL0) // choose only the case, when Lambda is the mother
							{
								mcTr_Mother->GetMomentum(mom_moth);
								
								//introducing cuts on pt and eta of lambda to check:
								//cout << "mom.Eta() = " << mom_moth.Eta() << endl;
								//cout << "mom_moth.Eta() = " << mom_moth.Eta() << "; " << mom_moth.Pt() << endl;
								
								//if (TMath::Abs(mom_moth.Eta())>=1. && mom_moth.Pt()>=3 && mom_moth.Pt()<=0.4) continue;
								
								int mcTr_Lam_MotherID = mcTr_Mother->GetMotherId(); // ID of the mother of Lambda
								Float_t polar_x = 0.;
								Float_t polar_y = 0.;
								Float_t polar_z = 0.;
								Float_t pol_length = 0.;
								Int_t sign = 1;
								
								Float_t weight_pol = mcTr_Mother->GetWeight();
								polar_x = mcTr_Mother->GetPolar(0);
								polar_y = mcTr_Mother->GetPolar(1);
								polar_z = mcTr_Mother->GetPolar(2);
								sign = TMath::Sign(1, mcTr_Mother->GetPdgCode());
									
								TVector3 polar_changed(polar_x, polar_y, polar_z); //polarization vector for distributions		
								if (phiEP_mc != 0.) polar_changed.RotateZ(-phiEP_mc);
								polar_x = weight_pol*polar_changed.X();
								polar_y = weight_pol*polar_changed.Y();
								polar_z = weight_pol*polar_changed.Z();
								pol_length = TMath::Sqrt(polar_x*polar_x + polar_y*polar_y + polar_z*polar_z);
								//for primary Lambda
								if (mcTr_Lam_MotherID < 0) 
								{
									Lpolar_y_prim[iter_cent]->Fill(polar_y);
									Lpolar_y_prim[iter_cent]->SetYTitle("Entries");
									Lpolar_y_prim[iter_cent]->SetXTitle("P_{y}");
									Lpolar_y_prim[iter_cent]->SetLineColor(kBlack);	
								}
								Lpolar_x[iter_cent]->Fill(polar_x);
								Lpolar_x[iter_cent]->SetYTitle("Entries");
								Lpolar_x[iter_cent]->SetXTitle("P_{x}");
								Lpolar_x[iter_cent]->SetLineColor(kBlack);	
									
								Lpolar_y[iter_cent]->Fill(polar_y);
								Lpolar_y[iter_cent]->SetYTitle("Entries");
								Lpolar_y[iter_cent]->SetXTitle("P_{y}");
								Lpolar_y[iter_cent]->SetLineColor(kBlack);	
									
								Lpolar_z[iter_cent]->Fill(polar_z);
								Lpolar_z[iter_cent]->SetYTitle("Entries");
								Lpolar_z[iter_cent]->SetXTitle("P_{z}");
								Lpolar_z[iter_cent]->SetLineColor(kBlack);	
									
								Lpolar_norm[iter_cent]->Fill(pol_length);
								Lpolar_norm[iter_cent]->SetYTitle("Entries");
								Lpolar_norm[iter_cent]->SetXTitle("|P|");
								Lpolar_norm[iter_cent]->SetLineColor(kBlack);
									
								//proton values (MagThetaPhi)
								Float_t p_prot = mom.Mag();
								Float_t eta_prot = mom.Theta();
								Float_t phi_prot = mom.Phi();
								//lambda values (MagThetaPhi)	
								Float_t p_lam = mom_moth.Mag();
								Float_t eta_lam = mom_moth.Theta();
								Float_t phi_lam = mom_moth.Phi();
								TVector3 vPr, vLamb;
								vPr.SetMagThetaPhi(p_prot, eta_prot, phi_prot);
								vLamb.SetMagThetaPhi(p_lam, eta_lam, phi_lam);
									
								//calculate the costheta and phi from dataset:
								cos_prot = 0.0; 
								phi_prot_star = 0.0; 
								FindPolarAngle (vPr, vLamb, cos_prot, phi_prot_star);
								float phi_star = phi_prot_star;
								float phi_diff = 0.;
								
								if (angle_choice == "RP")
								{
									phi_diff = phiEP_mc - phi_prot_star;
								}else if (angle_choice == "EP")
								{
									phi_diff = phiEP_1_ZDC - phi_prot_star;		
								}else
								{
									cout << "This angle choice is not defined! Please choose either RP or EP" << endl;
									return 1;
								}
								
								if (phi_star < 0) phi_star = phi_star + 2.*pi;
								if (phi_diff < 0) phi_diff = phi_diff + 2.*pi;
								
								//fill distribution for all Lambda	
								CosTheta_hist[iter_cent]->Fill(cos_prot);
								Pstar_hist[iter_cent]->Fill(phi_star);
								PstarRP_hist[iter_cent]->Fill(phi_diff);
								
								//fill distribution for primary Lambda
								if (mcTr_Lam_MotherID < 0) 
								{
									CosTheta_hist_prim[iter_cent]->Fill(cos_prot);
									Pstar_hist_prim[iter_cent]->Fill(phi_star);
									PstarRP_hist_prim[iter_cent]->Fill(phi_diff);
								}
								//calculate the costheta for transverse polarization:
								cos_trans = 0.0; 
								Calculate_transverse_cos(vPr, vLamb, cos_trans);
								CosTheta_trans[iter_cent]->Fill(cos_trans);
									
								//number of positively/negatively polarized Lambdas:
								if(polar_changed.Y() > 0) Number_Positive[iter_cent]++;
								if(polar_changed.Y() < 0) Number_Negative[iter_cent]++;							
								//number of found Lambdas:
								Number_Lambda[iter_cent]++;									
							} // mcTr_Mother->GetPdgCode() == pdgCodeL0
								
						} // mcTr->GetPdgCode() == pdgCodeAPr	
					}else
					{
						cout << "This particle choice is not defined! Please provide the definition in the code." << endl;
						return 1;
					}
				} // j < nMC
			} // centrality bins 
		} // iter_cent < NITER_CENT
		
	} // i < events
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		NPositive->SetBinContent(iter_cent+1,Number_Positive[iter_cent]);
		NNegative->SetBinContent(iter_cent+1,Number_Negative[iter_cent]);
		NHyperon->SetBinContent(iter_cent+1,Number_Lambda[iter_cent]);	
		if (angle_choice == "EP")
		{
			Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);
			Resolution_EP1_exp->SetBinContent(iter_cent+1,SubEvRes1[iter_cent]);
		}
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
//Get costheta and phi of the daughter particle
void FindPolarAngle(TVector3 &vPr, TVector3 &vLamb, Double_t &cos_prot, Double_t &phi_prot_star)
{
	//Compute azimuthal angle of proton in the lambda frame
	cos_prot = 0.;
	phi_prot_star = 0.;
	
	TLorentzVector prLor, lambLor;
      
	prLor.SetVectM(vPr, 0.938272);
	lambLor.SetVectM(vLamb, 1.11568);

	//standard Lorentz boost
	TVector3 boostV;
	boostV = lambLor.BoostVector();
	boostV *= -1;
  	
  	prLor.Boost(boostV);
	vPr = prLor.Vect();	
	
	//calculating the azimuthal angle of proton in the lambda frame (phi*):
	
	cos_prot = vPr.CosTheta(); // cos theta
	phi_prot_star = vPr.Phi();
}
//Get costheta for transverse polarization
void Calculate_transverse_cos(TVector3 &vPr, TVector3 &vLamb, Double_t &cos_trans)
{
	//Compute decay proton angle w.r.t. lambda decay plane
	cos_trans = 0.;

	TLorentzVector prLor, lambLor;
      
	prLor.SetVectM(vPr, 0.938272);
	lambLor.SetVectM(vLamb, 1.11568);

	//standard Lorentz boost
	TVector3 boostV;
	boostV = lambLor.BoostVector();
	boostV *= -1;
	
	TVector3 vBeam(0,0,1);
  	TVector3 vPlane = vBeam.Cross(vLamb);
 
  	prLor.Boost(boostV);
	vPr = prLor.Vect();	
	
	//calculating the transverse cos theta:
	cos_trans = vPr * vPlane / (vPr.Mag() * vPlane.Mag());

}
//
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
