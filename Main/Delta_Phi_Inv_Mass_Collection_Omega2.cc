//Constructs necessary histograms for analysis

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdarg>

#include "TCanvas.h"
#include "TFile.h"
#include <TCut.h>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include <TString.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TChain.h>
#include <TGraphErrors.h>
using namespace std;

#include "Struct_L0_v1.h"
#define pi TMath::Pi()
const int CENT = 10; //amount of centrality intervals (for the centrality determination technique)

Float_t GetCentrality(Int_t multiplicity);

TRandom *RNG = new TRandom();
Int_t Border_max[CENT];
Int_t Border_min[CENT];
Float_t Centrality_bin[CENT];

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

int main(int argc, char** argv)
{
	
//default values for all the dimensions
	int NITER = 20; //amount of cos theta cuts (delta(phi))
	int NITER_CENT = 4; //amount of centrality intervals (for the analysis)	
	int n_ev = 0; // number of events to analyze (default 0 - means all)
	float cut_ptmin = 0.15; // minimal pt for analysis
	float cut_ptmax = 10.; // maximal pt for analysis
	float cut_eta = 1.; // eta cut for analysis
//default values for centrality cuts	
//	Int_t _N_Hits = 16; // default value
//	Double_t cut_pt = 0.15; // default value (GeV/c)
//	Double_t cut_eta = 0.5; // default value
//	Double_t dca_cut = 0.5; // default: 0.5 cm
	Int_t dca_choice = 0; // choice of DCA (0 - no DCA cut, 1 - with DCA cut)
//create the chain of input files:
	TChain chain("event");
	
//strings for the output filename and centrality filename:
	TString outname = ""; 
	TString centname = "";
	TString selections_values = "";

//checking command line arguments, in case I want different optional parameters (NITER, NITER_CORR, NITER_CENT, events (0 means all), inname, outname, centname)
	for (int i = 0 ; i < argc ; i++) 
	{
		string tmp(argv[i]);
		if (tmp.compare("-help") == 0) 
		{
			cout << "Choose the path for the input and output files, and optional parameters, in case they differ from the default ones." << endl;
			
			cout << "Usage: ./Reader_phsd_pol_output -NITER 6 -NITER_CORR 10 -NITER_CENT 9 -events 100 -inname /home/elizaveta/data/PHSD/pol/PHSD_new_polarization_transfer_full_pol_part*.root -outname ../data_files/PHSD_readout_final_cent_tpc_4bins_pol.root -centname /home/elizaveta/data/Centrality_calibration/PHSD_FINAL.root" << endl;
			//for laptop: -inname: /home/liza/data/phsd_full/Calibrated_centrality/pol/PHSD_new_polarization_transfer_full_pol_part*.root -centname: /home/liza/data/centrality/PHSD_FINAL.root
			return 0;
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
		if (tmp.compare("-events") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				n_ev = number;
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
		if (tmp.compare("-selectionsfile") == 0 && i+1 < argc) 
		{
			selections_values = argv[i+1];
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
	
	if(selections_values.Length() == 0)
	{
		cerr << "!! Please provide the path for the selections file !!" << endl;
		return 1;
	}

	//cout << "events chosen = " << n_ev << "; NITER = " << NITER << "; NITER_CENT = " << NITER_CENT << "; N_Hits = " << _N_Hits << "; cut_pt = " << cut_pt << "; cut_eta = " << cut_eta << "; dca_cut = " << dca_cut << "; dca_choice = " << dca_choice << endl;
	cout << "events chosen = " << n_ev << "; NITER = " << NITER << "; NITER_CENT = " << NITER_CENT <<  "; dca_choice = " << dca_choice << endl;
	
	cout << "entries in chain: = " << chain.GetEntries() << endl;
	
//define the parameters, histograms and other stuff:
	double y_nLamb[NITER_CENT], y_nLamb_MC[NITER_CENT], Lcount[NITER_CENT], cos_Lamb[NITER_CENT];
	
	double ResEP1_exp[NITER_CENT], SubEvRes1[NITER_CENT], ResEP1_true[NITER_CENT], ResEP1_corr[NITER_CENT]; 

	TH2D *Lambda_PhaseSpace[NITER_CENT], *Lambda_PhaseSpace_true[NITER_CENT];
	char *int_Lambda_PhaseSpace = new char[NITER_CENT];
	char *int_Lambda_PhaseSpace_true = new char[NITER_CENT];

	TH1D *hm0[NITER_CENT];
	TH1D *hm0_before_selection[NITER_CENT];
	char *int_hm0 = new char[NITER_CENT];
	char *int_hm0_before_selection = new char[NITER_CENT];

	TH1D *Lambda_InvMass[NITER_CENT][NITER];
	char *int_Lambda_InvMass = new char[NITER_CENT];

	TH1D *Lpolar[NITER_CENT];
	TH1D *Lpolar_prim[NITER_CENT]; 
	char *polar_int = new char[NITER_CENT];
	char *polar_prim_int = new char[NITER_CENT];

	TH1D *cosTheta_hist[NITER_CENT], *cosTheta_MC_hist[NITER_CENT], *DirFlow_Lambda[NITER_CENT];
	char *cosTheta_hist_int = new char[NITER_CENT];
	char *cosTheta_MC_hist_int = new char[NITER_CENT];
	char *DirFlow_Lambda_int = new char[NITER_CENT];
	
	TH1D *PstarEP_hist[NITER_CENT], *PstarEP_hist_prim[NITER_CENT], *PstarRP_hist[NITER_CENT], *PstarRP_hist_prim[NITER_CENT], *PstarRP_hist_MC[NITER_CENT], *PstarRP_hist_MC_prim[NITER_CENT];
	char *PstarEP_hist_int = new char[NITER_CENT];
	char *PstarEP_hist_prim_int = new char[NITER_CENT];
	char *PstarRP_hist_int = new char[NITER_CENT];
	char *PstarRP_hist_prim_int = new char[NITER_CENT];
	char *PstarRP_hist_MC_int = new char[NITER_CENT];
	char *PstarRP_hist_MC_prim_int = new char[NITER_CENT];
	
	TLatex latex;
	latex.SetNDC();

//calculated values for the omega2 cut, and the cut values for centrality:	
	double *omega2_cutvalues;
	int *centrality_min;
	int *centrality_max;
	double *_CentrBins;
	double *noErr;
	
	if (NITER_CENT == 4)
	{				
		omega2_cutvalues = init_double_array(4, 0, 0.0, 0.0, 0.0, 0.0); //assigning zero values for omega_2 cut selection
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		_CentrBins = init_double_array(5, 0, 0.,10.,20.,50.,100.);
		noErr = init_double_array(4, 0, 0., 0., 0., 0.);		
	}
	if (NITER_CENT == 7) //need to recheck - the last value is fishy
	{		
		omega2_cutvalues = init_double_array(7, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //assigning zero values for omega_2 cut selection
		centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
		_CentrBins = init_double_array(8, 0, 0., 10., 20., 30., 40., 50., 60., 70.);
		noErr = init_double_array(7, 0, 0., 0., 0., 0., 0., 0., 0.);		
	}
	if (NITER_CENT == 8) //need to recheck - the last value is fishy
	{
		omega2_cutvalues = init_double_array(8, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //assigning zero values for omega_2 cut selection
		centrality_min = init_int_array(8, 0, 0, 10, 20, 30, 40, 50, 60, 70);
		centrality_max = init_int_array(8, 0, 10, 20, 30, 40, 50, 60, 70, 80);
		_CentrBins = init_double_array(9, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80.);
		noErr = init_double_array(8, 0, 0., 0., 0., 0., 0., 0., 0., 0.);		
	}
	if (NITER_CENT == 10)
	{
		omega2_cutvalues = init_double_array(10, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //assigning zero values for omega_2 cut selection
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		_CentrBins = init_double_array(11, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.);
		noErr = init_double_array(10, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);		
	}
	
	//reading the omega_2 values from the file:
	ifstream selections_file;
	selections_file.open(selections_values);
	Int_t nlines = 0;
	while (1) 
	{
		selections_file >> omega2_cutvalues[nlines];
		if (!selections_file.good()) break;
		nlines++;
	}
	printf("found %d points\n",nlines);
	if(nlines != NITER_CENT)
	{
		cout << "Amount of omega_2 selection cuts is not equal to the centrality bins: " << "; nlines = " << nlines << "; NITER_CENT = " << NITER_CENT << endl;
		return 1;
	}
	selections_file.close();
	
	for (int i = 0; i < NITER_CENT; i++)
	{
		cout << "omega2_cutvalues[i] = " <<  omega2_cutvalues[i] << "; centrality_min[i] = " <<  centrality_min[i] <<  "; centrality_max[i] = " <<  centrality_max[i] <<  "; _CentrBins[i] = " <<  _CentrBins[i] << endl;
	}
	//return 1; //debugging
	
	double *cos_min;
	double *cos_max;
	double *angle_min;
	double *angle_max;
	
//range (-1.0;1.0) for 6 or 20 cos(theta) bins, angle range (0.;6.28) for 6 or 20 bins
	if (NITER == 6)
	{		
		cos_min = init_double_array(6, 0, -1., -0.666, -0.333, 0., 0.333, 0.666);
		cos_max = init_double_array(6, 0, -0.666, -0.333, 0., 0.333, 0.666, 1.);
		angle_min = init_double_array(6, 0, 0., 1.046, 2.092, 3.138, 4.184, 5.23);
		angle_max = init_double_array(6, 0, 1.046, 2.092, 3.138, 4.184, 5.23, 6.28);
	}
	if (NITER == 12)
	{	
		/*	
		double xmin_anglemin = 0.;
		double xmax_anglemax = 2.*pi;
		double xmin_cosmin = -1.;
		double xmax_cosmax = 1.;
		
		double step_angle = (xmax_anglemax - xmin_anglemin)/NITER;
		double step_cos = (xmax_cosmax - xmin_cosmin)/NITER;
		
		cos_min = init_double_array(12, xmin_cosmin, xmin_cosmin + 1.*step_cos, xmin_cosmin + 2.*step_cos, xmin_cosmin + 3.*step_cos, xmin_cosmin + 4.*step_cos, xmin_cosmin + 5.*step_cos, xmin_cosmin + 6.*step_cos, xmin_cosmin + 7.*step_cos, xmin_cosmin + 8.*step_cos, xmin_cosmin + 9.*step_cos, xmin_cosmin + 10.*step_cos, xmin_cosmin + 11.*step_cos);
		cos_max = init_double_array(12, xmin_cosmin + 1.*step_cos, xmin_cosmin + 2.*step_cos, xmin_cosmin + 3.*step_cos, xmin_cosmin + 4.*step_cos, xmin_cosmin + 5.*step_cos, xmin_cosmin + 6.*step_cos, xmin_cosmin + 7.*step_cos, xmin_cosmin + 8.*step_cos, xmin_cosmin + 9.*step_cos, xmin_cosmin + 10.*step_cos, xmin_cosmin + 11.*step_cos, xmax_cosmax);
		angle_min = init_double_array(12, xmin_anglemin, xmin_anglemin + 1.*step_angle, xmin_anglemin + 2.*step_angle, xmin_anglemin + 3.*step_angle, xmin_anglemin + 4.*step_angle, xmin_anglemin + 5.*step_angle, xmin_anglemin + 6.*step_angle, xmin_anglemin + 7.*step_angle, xmin_anglemin + 8.*step_angle, xmin_anglemin + 9.*step_angle, xmin_anglemin + 10.*step_angle, xmin_anglemin + 11.*step_angle);
		angle_max = init_double_array(12, xmin_anglemin + 1.*step_angle, xmin_anglemin + 2.*step_angle, xmin_anglemin + 3.*step_angle, xmin_anglemin + 4.*step_angle, xmin_anglemin + 5.*step_angle, xmin_anglemin + 6.*step_angle, xmin_anglemin + 7.*step_angle, xmin_anglemin + 8.*step_angle, xmin_anglemin + 9.*step_angle, xmin_anglemin + 10.*step_angle, xmin_anglemin + 11.*step_angle, xmax_anglemax);
		*/
		cos_min = init_double_array(12, 0, -1., -0.833, -0.666, -0.499, -0.332, -0.165, 0.002, 0.169, 0.336, 0.503, 0.67, 0.837);
		cos_max = init_double_array(12, 0, -0.833, -0.666, -0.499, -0.332, -0.165, 0.002, 0.169, 0.336, 0.503, 0.67, 0.837, 1.);
		
		angle_min = init_double_array(12, 0, 0., 0.523, 1.046, 1.569, 2.092, 2.615, 3.138, 3.661, 4.184, 4.707, 5.23, 5.753);
		angle_max = init_double_array(12, 0, 0.523, 1.046, 1.569, 2.092, 2.615, 3.138, 3.661, 4.184, 4.707, 5.23, 5.753, 6.28);
	}
	if (NITER == 20)
	{		
		cos_min = init_double_array(20, 0, -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9);
		cos_max = init_double_array(20, 0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.);
		
		angle_min = init_double_array(20, 0, 0., 0.314, 0.628, 0.942, 1.256, 1.57, 1.884, 2.198, 2.512, 2.826, 3.14, 3.454, 3.768, 4.082, 4.396, 4.71, 5.024, 5.338, 5.652, 5.966);
		angle_max = init_double_array(20, 0, 0.314, 0.628, 0.942, 1.256, 1.57, 1.884, 2.198, 2.512, 2.826, 3.14, 3.454, 3.768, 4.082, 4.396, 4.71, 5.024, 5.338, 5.652, 5.966, 6.28);
	}
	for (int i = 0; i < NITER; i++)
	{
		cout << "cos_min[i] = " <<  cos_min[i] << "; cos_max[i] = " <<  cos_max[i] <<  "; angle_min[i] = " <<  angle_min[i] <<  "; angle_max[i] = " <<  angle_max[i] << endl;
	}
	
	double centrality_bin[NITER_CENT]; 	
	double polar_par_value[NITER_CENT];
	double polar_par_value_err[NITER_CENT];
	double polar_prim_par_value[NITER_CENT];
	double polar_prim_par_value_err[NITER_CENT];
	int lost_lambda[NITER_CENT];
	
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
	
//define variables for the input tree:
	Float_t b0, phiEP_mc, phiEP_mc_or, phiEP_1_ZDC, phiEP_2_ZDC, ResEP_1_ZDC, ResEP_2_ZDC, ResEPSub_1_ZDC, ResEPSub_2_ZDC;
	Int_t Centrality_tpc, n_tracks_mpd, n_tracks_mpd_DCA, nLamb, nLamb_MC;
	Float_t phiEP_1_ZDC_corr, phiEP_2_ZDC_corr;
	
//set adresses to the branches:
	vector<L0>  *l0 = 0;
	chain.SetBranchAddress("l0", &l0);  
	chain.SetBranchAddress("b0", &b0);
	chain.SetBranchAddress("phiEP_mc", &phiEP_mc);
	chain.SetBranchAddress("phiEP_mc_or", &phiEP_mc_or);
	chain.SetBranchAddress("n_tracks_mpd", &n_tracks_mpd);
	chain.SetBranchAddress("n_tracks_mpd_DCA", &n_tracks_mpd_DCA);
	chain.SetBranchAddress("nLamb", &nLamb);
	chain.SetBranchAddress("nLamb_MC", &nLamb_MC);
	chain.SetBranchAddress("phiEP_1_ZDC", &phiEP_1_ZDC);
	chain.SetBranchAddress("phiEP_2_ZDC", &phiEP_2_ZDC);
	chain.SetBranchAddress("ResEP_1_ZDC", &ResEP_1_ZDC);
	chain.SetBranchAddress("ResEP_2_ZDC", &ResEP_2_ZDC);
	chain.SetBranchAddress("ResEPSub_1_ZDC", &ResEPSub_1_ZDC);
	chain.SetBranchAddress("ResEPSub_2_ZDC", &ResEPSub_2_ZDC);
	
//output file:		
	TFile out(outname,"recreate");
	
//the necessary histograms for analysis/tests
	TH1F *NCentr = new TH1F("NCentr","NCentr",NITER_CENT,_CentrBins);
	NCentr->SetYTitle("Entries");
	NCentr->SetXTitle("Centrality TPC");
	NCentr->SetLineColor(kBlack);
	
	TH2D *Lambda_PhaseSpace_MB = new TH2D("Lambda_PhaseSpace_MB","Lambda_PhaseSpace_MB",100,-1.5,1.5,100.,0.,3.);
	TH2D *Lambda_PhaseSpace_MB_true = new TH2D("Lambda_PhaseSpace_MB_true","Lambda_PhaseSpace_MB_true",100,-1.5,1.5,100.,0.,3.);
	
	TH1D *dNLambda = new TH1D("dNLambda","dNLambda",NITER_CENT,_CentrBins);
	dNLambda->SetYTitle("N_{#Lambda}/N_{events}");
	dNLambda->SetXTitle("Centrality, [%]");
	dNLambda->SetMarkerStyle(20);
	dNLambda->SetLineColor(kBlack);
	dNLambda->SetMarkerColor(kBlack);
	dNLambda->SetMarkerSize(2);
	dNLambda->SetLineWidth(2);
	
	TH1D *dNLambda_MC = new TH1D("dNLambda_MC","dNLambda_MC",NITER_CENT,_CentrBins);
	dNLambda_MC->SetYTitle("N_{#Lambda}/N_{events}");
	dNLambda_MC->SetXTitle("Centrality, [%]");
	dNLambda_MC->SetMarkerStyle(22);
	dNLambda_MC->SetLineColor(kRed);
	dNLambda_MC->SetMarkerColor(kRed);
	dNLambda_MC->SetMarkerSize(2);
	dNLambda_MC->SetLineWidth(2);
	
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
	
	TH1D *Counts_hist = new TH1D("Counts_hist","Counts_hist",NITER_CENT,_CentrBins);
	Counts_hist->SetYTitle("Counts");
	Counts_hist->SetXTitle("Centrality, [%]");
	Counts_hist->SetMarkerStyle(20);
	Counts_hist->SetLineColor(kBlue);
	Counts_hist->SetMarkerColor(kBlue);
	Counts_hist->SetMarkerSize(2);
	Counts_hist->SetLineWidth(2);
	
	TH1D *Lost_Lambdas = new TH1D("Lost_Lambdas","Lost_Lambdas",NITER_CENT,_CentrBins);
	Lost_Lambdas->SetYTitle("Counts");
	Lost_Lambdas->SetXTitle("Centrality, [%]");
	Lost_Lambdas->SetMarkerStyle(20);
	Lost_Lambdas->SetLineColor(kBlue);
	Lost_Lambdas->SetMarkerColor(kBlue);
	Lost_Lambdas->SetMarkerSize(2);
	Lost_Lambdas->SetLineWidth(2);
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Lcount[iter_cent] = 0;
		y_nLamb[iter_cent] = 0;
		y_nLamb_MC[iter_cent] = 0;
		cos_Lamb[iter_cent] = 0;
		lost_lambda[iter_cent] = 0;
				
		sprintf(int_Lambda_PhaseSpace,"Lambda_PhaseSpace_%d",iter_cent);
		Lambda_PhaseSpace[iter_cent] = new TH2D(int_Lambda_PhaseSpace,int_Lambda_PhaseSpace,100,-1.5,1.5,100.,0.,3.);
		sprintf(int_Lambda_PhaseSpace_true,"Lambda_PhaseSpace_true_%d",iter_cent);
		Lambda_PhaseSpace_true[iter_cent] = new TH2D(int_Lambda_PhaseSpace_true,int_Lambda_PhaseSpace_true,100,-1.5,1.5,100.,0.,3.);
		
		sprintf(int_hm0,"hm0_%d",iter_cent);
		hm0[iter_cent] = new TH1D(int_hm0,int_hm0,100, 1.07, 1.17);
		sprintf(int_hm0_before_selection,"hm0_before_selection_%d",iter_cent);
		hm0_before_selection[iter_cent] = new TH1D(int_hm0_before_selection,int_hm0_before_selection,100, 1.07, 1.17);
		
		sprintf(polar_int,"LPolar_%d",iter_cent);      
		Lpolar[iter_cent] = new TH1D(polar_int,polar_int,100,-1.,1.);
		
		sprintf(polar_prim_int,"LPolar_prim_%d",iter_cent);      
		Lpolar_prim[iter_cent] = new TH1D(polar_prim_int,polar_prim_int,100,-1.,1.);
		
		sprintf(cosTheta_hist_int,"cosTheta_hist_%d",iter_cent);      
		cosTheta_hist[iter_cent] = new TH1D(cosTheta_hist_int,cosTheta_hist_int,NITER,-1.0,1.0);
		
		sprintf(cosTheta_MC_hist_int,"cosTheta_MC_hist_%d",iter_cent);      
		cosTheta_MC_hist[iter_cent] = new TH1D(cosTheta_MC_hist_int,cosTheta_MC_hist_int,NITER,-1.0,1.0);
		
		//redid as I had in MC tests:
		sprintf(PstarRP_hist_int,"PstarRP_hist_%d",iter_cent);      
		PstarRP_hist[iter_cent] = new TH1D(PstarRP_hist_int,PstarRP_hist_int,NITER,0.,2.*pi);
		
		sprintf(PstarRP_hist_prim_int,"PstarRP_hist_prim_%d",iter_cent);      
		PstarRP_hist_prim[iter_cent] = new TH1D(PstarRP_hist_prim_int,PstarRP_hist_prim_int,NITER,0.,2.*pi);
		
		sprintf(PstarEP_hist_int,"PstarEP_hist_%d",iter_cent);      
		PstarEP_hist[iter_cent] = new TH1D(PstarEP_hist_int,PstarEP_hist_int,NITER,0.,2.*pi);
		
		sprintf(PstarEP_hist_prim_int,"PstarEP_hist_prim_%d",iter_cent);      
		PstarEP_hist_prim[iter_cent] = new TH1D(PstarEP_hist_prim_int,PstarEP_hist_prim_int,NITER,0.,2.*pi);
		
		sprintf(PstarRP_hist_MC_int,"PstarRP_hist_MC_%d",iter_cent);      
		PstarRP_hist_MC[iter_cent] = new TH1D(PstarRP_hist_MC_int,PstarRP_hist_MC_int,NITER,0.,2.*pi);
		
		sprintf(PstarRP_hist_MC_prim_int,"PstarRP_hist_MC_prim_prim_%d",iter_cent);      
		PstarRP_hist_MC_prim[iter_cent] = new TH1D(PstarRP_hist_MC_prim_int,PstarRP_hist_MC_prim_int,NITER,0.,2.*pi);
		
		sprintf(DirFlow_Lambda_int,"DirFlow_Lambda_%d",iter_cent);      
		DirFlow_Lambda[iter_cent] = new TH1D(DirFlow_Lambda_int,DirFlow_Lambda_int,6,-1.3,1.3);	
		
		for(int iter_cos = 0; iter_cos < NITER; iter_cos++)
		{	
			sprintf(int_Lambda_InvMass,"Lambda_mass_%d_%d",iter_cent,iter_cos);
			Lambda_InvMass[iter_cent][iter_cos] = new TH1D(int_Lambda_InvMass,int_Lambda_InvMass,100, 1.07, 1.17);
		}
		
	}
	
	Int_t full_events = chain.GetEntries();
	cout << " Full number of Events = " << full_events << "; Chosen number of events = " << n_ev << endl;
	if (n_ev != 0) full_events = TMath::Min (full_events, n_ev);
	cout << " Number of events = " << full_events << endl;
	
	Int_t count_ev = 0; //for check
	Int_t zero_events = 0; //for check
	for(Int_t iev = 0; iev < full_events; ++iev) //cycle for entries
	{   	
		chain.GetEntry(iev); // read current event in chain
		cout << " Event: " << iev << "\t\t\t\t\r"<<flush;
		
		if(dca_choice == 0) 
		{
			Centrality_tpc = GetCentrality(n_tracks_mpd);
			if(n_tracks_mpd == 0) 
			{ 
				zero_events++; 
				continue;
			}
		}else 
		if(dca_choice == 1)
		{
			Centrality_tpc = GetCentrality(n_tracks_mpd_DCA);
			if(n_tracks_mpd_DCA == 0) 
			{ 
				zero_events++; 
				continue;
			}
		}else 
		{
			cout << " Please choose a correct version of DCA choice " << endl;
			return 1;
		}	
		
		NCentr->Fill(Centrality_tpc);
		
		L0 *lamb = nullptr;
		
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			if(Centrality_tpc < centrality_max[iter_cent] && Centrality_tpc >= centrality_min[iter_cent])
			{
				Lcount[iter_cent]++;
				y_nLamb[iter_cent] += nLamb;
				y_nLamb_MC[iter_cent] += nLamb_MC; 
				SubEvRes1[iter_cent] += ResEPSub_1_ZDC;
				ResEP1_true[iter_cent] += ResEP_1_ZDC;
				
				//cout << "l0->size() = " << l0->size() << endl;
				//cycle for the reconstructed lambda:
				for(UInt_t i = 0; i < l0->size(); ++i) 
				{
					lamb = (L0*) &l0->at(i);
					hm0_before_selection[iter_cent]->Fill(lamb->massh);	
					//check omega2:				
					//cout << "omega2 = " << lamb->omega2 << endl;
					//cout << "omega2 = " << omega2_cutvalues[iter_cent] << endl;
					if(lamb->omega2>omega2_cutvalues[iter_cent]) 
					{							
						Lambda_PhaseSpace_MB->GetYaxis()->SetTitle("p_{T}, GeV/c");
						Lambda_PhaseSpace_MB->GetXaxis()->SetTitle("y");
						Lambda_PhaseSpace_MB->Fill(lamb->yh,lamb->pth);
						Lambda_PhaseSpace_MB->SetTitle("Lambda Phase Space MB");
						
						Lambda_PhaseSpace[iter_cent]->GetYaxis()->SetTitle("p_{T}, GeV/c");
						Lambda_PhaseSpace[iter_cent]->GetXaxis()->SetTitle("y");
						Lambda_PhaseSpace[iter_cent]->Fill(lamb->yh,lamb->pth);
						Lambda_PhaseSpace[iter_cent]->SetTitle("Lambda Phase Space");
						
						hm0[iter_cent]->Fill(lamb->massh); 
							
						float phi_diff = phiEP_1_ZDC - lamb->phi_star;
						float phi_diff_MC = phiEP_mc - lamb->phi_star_MC;
						if (phi_diff < 0) phi_diff = phi_diff + 2.*pi;
						if (phi_diff_MC < 0) phi_diff_MC = phi_diff_MC + 2.*pi;
												
						for(int iter_cos = 0; iter_cos < NITER; iter_cos++)
						{
							if(phi_diff < angle_max[iter_cos] && phi_diff >= angle_min[iter_cos])
							{
								Lambda_InvMass[iter_cent][iter_cos]->Fill(lamb->massh);
							}
						}
							
						
					}
					if (lamb->origs[0] > 0)
					{	
						Lambda_PhaseSpace_MB_true->GetYaxis()->SetTitle("p_{T}, GeV/c");
						Lambda_PhaseSpace_MB_true->GetXaxis()->SetTitle("y");
						Lambda_PhaseSpace_MB_true->Fill(lamb->yh,lamb->pth);
						Lambda_PhaseSpace_MB_true->SetTitle("Lambda (true) Phase Space MB");
						
						Lambda_PhaseSpace_true[iter_cent]->GetYaxis()->SetTitle("p_{T}, GeV/c");
						Lambda_PhaseSpace_true[iter_cent]->GetXaxis()->SetTitle("y");
						Lambda_PhaseSpace_true[iter_cent]->Fill(lamb->yh,lamb->pth);
						Lambda_PhaseSpace_true[iter_cent]->SetTitle("Lambda (true) Phase Space");
						
						Lpolar[iter_cent]->Fill(lamb->polarhy);
						Lpolar[iter_cent]->SetYTitle("Entries");
						Lpolar[iter_cent]->SetXTitle("P_{y}");
						Lpolar[iter_cent]->SetLineColor(kBlack);
							
						cosTheta_hist[iter_cent]->Fill(lamb->cosA);
						cosTheta_MC_hist[iter_cent]->Fill(lamb->cosAmc);
						
						float phi_diff_hist = phiEP_1_ZDC - lamb->phi_star;
						float phi_diff_histRP = phiEP_mc - lamb->phi_star;
						float phi_diff_MC = phiEP_mc - lamb->phi_star_MC;
						if (phi_diff_hist < 0) phi_diff_hist = phi_diff_hist + 2.*pi;
						if (phi_diff_histRP < 0) phi_diff_histRP = phi_diff_histRP + 2.*pi;	
						if (phi_diff_MC < 0) phi_diff_MC = phi_diff_MC + 2.*pi;	
							
						PstarEP_hist[iter_cent]->Fill(phi_diff_hist); 
						PstarRP_hist[iter_cent]->Fill(phi_diff_histRP);
						PstarRP_hist_MC[iter_cent]->Fill(phi_diff_MC);
							
						DirFlow_Lambda[iter_cent]->Fill(lamb->yh); 
						
						if(lamb->omega2<=omega2_cutvalues[iter_cent]) 
						{
							cout << "Lost true Lambda!" << endl;
							lost_lambda[iter_cent]++;
						}							
					}
					if (lamb->origs[0] == 1) //true lambda (primary?)
					{	
						Lpolar_prim[iter_cent]->Fill(lamb->polarhy);
						Lpolar_prim[iter_cent]->SetYTitle("Entries");
						Lpolar_prim[iter_cent]->SetXTitle("P_{y}");
						Lpolar_prim[iter_cent]->SetLineColor(kRed);
						
						float phi_diff_hist = phiEP_1_ZDC - lamb->phi_star;
						float phi_diff_histRP = phiEP_mc - lamb->phi_star;
						float phi_diff_MC = phiEP_mc - lamb->phi_star_MC;
						if (phi_diff_hist < 0) phi_diff_hist = phi_diff_hist + 2.*pi;
						if (phi_diff_histRP < 0) phi_diff_histRP = phi_diff_histRP + 2.*pi;	
						if (phi_diff_MC < 0) phi_diff_MC = phi_diff_MC + 2.*pi;	
							
						PstarEP_hist_prim[iter_cent]->Fill(phi_diff_hist); 
						PstarRP_hist_prim[iter_cent]->Fill(phi_diff_histRP);
						PstarRP_hist_MC_prim[iter_cent]->Fill(phi_diff_MC);
					}
				}
			}
		}
		count_ev++;	
	}
	cout << "Actual number of analyzed events = " << count_ev << endl;
	cout << "Tossed out events = " << zero_events << endl;
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		cout << "iter = " << iter_cent << "; y_nLamb = " << y_nLamb[iter_cent] << "; y_nLamb_MC = " << y_nLamb_MC[iter_cent] << "; Lcount = " << Lcount[iter_cent] << endl;
		y_nLamb[iter_cent] = y_nLamb[iter_cent]/Lcount[iter_cent];
		y_nLamb_MC[iter_cent] = y_nLamb_MC[iter_cent]/Lcount[iter_cent];
		cout << "iter = " << iter_cent << "; y_nLamb (per event) = " << y_nLamb[iter_cent] <<"; y_nLamb_MC (per event) = " << y_nLamb_MC[iter_cent] << endl; 
		dNLambda->SetBinContent(iter_cent+1,y_nLamb[iter_cent]);
		dNLambda_MC->SetBinContent(iter_cent+1,y_nLamb_MC[iter_cent]);
		
		Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);
		Resolution_EP1_exp->SetBinContent(iter_cent+1,SubEvRes1[iter_cent]);
		Counts_hist->SetBinContent(iter_cent+1,Lcount[iter_cent]);
		Lost_Lambdas->SetBinContent(iter_cent+1,lost_lambda[iter_cent]);
		
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		cout << "iter_cent = " << iter_cent << "; Lost Lambdas: " << lost_lambda[iter_cent] << endl;
	}
	out.Write();
	out.Close();
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
