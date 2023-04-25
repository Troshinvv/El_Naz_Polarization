//Macro to find optimal value of omega_2 parameter

#include <iostream>
#include <vector>

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
#include <TRandom.h>
#include <TChain.h>
#include <TGraph.h>
using namespace std;

#include "Struct_L0_v1.h"

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
	int NITER = 20; // amount of test values of omega_2
	int NITER_CENT = 4; // amount of centrality intervals (for the analysis)	
	int n_ev = 0; // number of events to analyze (default 0 - means all)
	double cent_cut = 70.; // value for centrality cut
	Int_t cent_cut_choice = 0; // choice of centrality cut (0 - no cent cut, 1 - with cent cut)
	
	
//create the chain of input files:
	TChain chain("event");
	
//strings for the output filename and centrality filename:
	TString outname = ""; 
	TString centname = "";
	
	//checking command line arguments, in case I want different optional parameters (NITER, NITER_CORR, NITER_CENT, events (0 means all), inname, outname, centname)
	for (int i = 0 ; i < argc ; i++) 
	{
		string tmp(argv[i]);
		if (tmp.compare("-help") == 0) 
		{
			cout << "Choose the path for the input and output files, and optional parameters, in case they differ from the default ones." << endl;
			
			cout << "Usage: ./Omega2_test_CENT -NITER 10 -NITER_CENT 10 -events 100 -inname /scratch2/nazarova/phsd_output/Anal*.root -outname PHSD_omegatest.root -centname /scratch2/nazarova/CentralityFramework/Framework/results/PHSD_9GeV_b20_noDCA/FINAL.root" << endl;
			return 0;
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
	
	cout << "events chosen = " << n_ev << "; NITER = " << NITER << "; NITER_CENT = " << NITER_CENT << "; CENT_CUT = " << cent_cut << "; cent_cut_choice = " << cent_cut_choice << endl;
	
	//histrograms of selection parameters for full set:
	TH1D *Dca_pion[NITER_CENT], *Dca_proton[NITER_CENT], *Chi_pion[NITER_CENT], *Chi_proton[NITER_CENT], *Dca_lambda[NITER_CENT], *Chi_lambda[NITER_CENT], *Dca_v0[NITER_CENT], *Chi_v0[NITER_CENT], *Path_hist[NITER_CENT], *Angle_hist[NITER_CENT];
	char *int_Dca_pion = new char[NITER_CENT];
	char *int_Dca_proton = new char[NITER_CENT];
	char *int_Chi_pion = new char[NITER_CENT];
	char *int_Chi_proton = new char[NITER_CENT];
	char *int_Dca_lambda = new char[NITER_CENT];
	char *int_Chi_lambda = new char[NITER_CENT];
	char *int_Dca_v0 = new char[NITER_CENT];
	char *int_Chi_v0 = new char[NITER_CENT];
	char *int_Path_hist = new char[NITER_CENT];
	char *int_Angle_hist = new char[NITER_CENT];
	
	//histrograms of selection parameters for true lambda:
	TH1D *Dca_pion_true[NITER_CENT], *Dca_proton_true[NITER_CENT], *Chi_pion_true[NITER_CENT], *Chi_proton_true[NITER_CENT], *Dca_lambda_true[NITER_CENT], *Chi_lambda_true[NITER_CENT], *Dca_v0_true[NITER_CENT], *Chi_v0_true[NITER_CENT], *Path_hist_true[NITER_CENT], *Angle_hist_true[NITER_CENT];
	char *int_Dca_pion_true = new char[NITER_CENT];
	char *int_Dca_proton_true = new char[NITER_CENT];
	char *int_Chi_pion_true = new char[NITER_CENT];
	char *int_Chi_proton_true = new char[NITER_CENT];
	char *int_Dca_lambda_true = new char[NITER_CENT];
	char *int_Chi_lambda_true = new char[NITER_CENT];
	char *int_Dca_v0_true = new char[NITER_CENT];
	char *int_Chi_v0_true = new char[NITER_CENT];
	char *int_Path_hist_true = new char[NITER_CENT];
	char *int_Angle_hist_true = new char[NITER_CENT];

	int *centrality_min;
	int *centrality_max;
	double *_CentrBins;
	if (NITER_CENT == 4)
	{		
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		_CentrBins = init_double_array(5, 0, 0.,10.,20.,50.,100.);	
	}
	if (NITER_CENT == 7)
	{		
		centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
		_CentrBins = init_double_array(8, 0, 0., 10., 20., 30., 40., 50., 60., 70.);	
	}
	if (NITER_CENT == 8)
	{		
		centrality_min = init_int_array(8, 0, 0, 10, 20, 30, 40, 50, 60, 70);
		centrality_max = init_int_array(8, 0, 10, 20, 30, 40, 50, 60, 70, 80);
		_CentrBins = init_double_array(9, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80.);	
	}
	if (NITER_CENT == 10)
	{
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		_CentrBins = init_double_array(11, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.);
	}
	
	for (int i = 0; i < NITER_CENT; i++)
	{
		cout << "centrality_min[i] = " <<  centrality_min[i] <<  "; centrality_max[i] = " <<  centrality_max[i] <<  "; _CentrBins[i] = " <<  _CentrBins[i] << endl;
	}
	
	TLatex latex;
	latex.SetNDC();
	
	double centrality_bin[NITER_CENT]; 	
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
	
	Int_t Centrality_tpc, n_tracks_mpd, n_tracks_mpd_DCA;
	Float_t b0;
	
	//set adress to l0:
	vector<L0>  *l0 = 0;
	chain.SetBranchAddress("l0", &l0);  
	chain.SetBranchAddress("b0", &b0);
	chain.SetBranchAddress("n_tracks_mpd", &n_tracks_mpd);
	chain.SetBranchAddress("n_tracks_mpd_DCA", &n_tracks_mpd_DCA);
	
//create output file:
	TFile out(outname,"recreate");
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sprintf(int_Dca_pion,"Dca_pion_%d",iter_cent);
		Dca_pion[iter_cent] = new TH1D(int_Dca_pion,int_Dca_pion,100, 0.0, 160.);
		Dca_pion[iter_cent]->SetYTitle("Entries");
		Dca_pion[iter_cent]->SetXTitle("DCA_{#pi}");
		sprintf(int_Dca_proton,"Dca_proton_%d",iter_cent);
		Dca_proton[iter_cent] = new TH1D(int_Dca_proton,int_Dca_proton,100, 0.0, 160.);
		Dca_proton[iter_cent]->SetYTitle("Entries");
		Dca_proton[iter_cent]->SetXTitle("DCA_{p}");
		sprintf(int_Chi_pion,"Chi_pion_%d",iter_cent);
		Chi_pion[iter_cent] = new TH1D(int_Chi_pion,int_Chi_pion,100, 0.0, 10000.);
		Chi_pion[iter_cent]->SetYTitle("Entries");
		Chi_pion[iter_cent]->SetXTitle("#chi^{2}_{#pi}");
		sprintf(int_Chi_proton,"Chi_proton_%d",iter_cent);
		Chi_proton[iter_cent] = new TH1D(int_Chi_proton,int_Chi_proton,100, 0.0, 10000.);
		Chi_proton[iter_cent]->SetYTitle("Entries");
		Chi_proton[iter_cent]->SetXTitle("#chi^{2}_{p}");
		sprintf(int_Dca_lambda,"Dca_lambda_%d",iter_cent);
		Dca_lambda[iter_cent] = new TH1D(int_Dca_lambda,int_Dca_lambda,100, 0.0, 200.);
		Dca_lambda[iter_cent]->SetYTitle("Entries");
		Dca_lambda[iter_cent]->SetXTitle("DCA_{#Lambda}");
		sprintf(int_Chi_lambda,"Chi_lambda_%d",iter_cent);
		Chi_lambda[iter_cent] = new TH1D(int_Chi_lambda,int_Chi_lambda,100, 0.0, 10000.);
		Chi_lambda[iter_cent]->SetYTitle("Entries");
		Chi_lambda[iter_cent]->SetXTitle("#chi^{2}_{#Lambda}");
		sprintf(int_Dca_v0,"Dca_v0_%d",iter_cent);
		Dca_v0[iter_cent] = new TH1D(int_Dca_v0,int_Dca_v0,100, 0.0, 140.);
		Dca_v0[iter_cent]->SetYTitle("Entries");
		Dca_v0[iter_cent]->SetXTitle("DCA_{V_{0}}");
		sprintf(int_Chi_v0,"Chi_v0_%d",iter_cent);
		Chi_v0[iter_cent] = new TH1D(int_Chi_v0,int_Chi_v0,100, 0.0, 25.);
		Chi_v0[iter_cent]->SetYTitle("Entries");
		Chi_v0[iter_cent]->SetXTitle("#chi^{2}_{V_{0}}");
		sprintf(int_Path_hist,"Path_hist_%d",iter_cent);
		Path_hist[iter_cent] = new TH1D(int_Path_hist,int_Path_hist,100, 0.0, 300.);
		Path_hist[iter_cent]->SetYTitle("Entries");
		Path_hist[iter_cent]->SetXTitle("Path");
		sprintf(int_Angle_hist,"Angle_hist_%d",iter_cent);
		Angle_hist[iter_cent] = new TH1D(int_Angle_hist,int_Angle_hist,100, 0.0, 1.6);
		Angle_hist[iter_cent]->SetYTitle("Entries");
		Angle_hist[iter_cent]->SetXTitle("Angle");
		
		
		sprintf(int_Dca_pion_true,"Dca_pion_true_%d",iter_cent);
		Dca_pion_true[iter_cent] = new TH1D(int_Dca_pion_true,int_Dca_pion_true,100, 0.0, 160.);
		Dca_pion_true[iter_cent]->SetYTitle("Entries");
		Dca_pion_true[iter_cent]->SetXTitle("DCA_{#pi}");
		sprintf(int_Dca_proton_true,"Dca_proton_true_%d",iter_cent);
		Dca_proton_true[iter_cent] = new TH1D(int_Dca_proton_true,int_Dca_proton_true,100, 0.0, 160.);
		Dca_proton_true[iter_cent]->SetYTitle("Entries");
		Dca_proton_true[iter_cent]->SetXTitle("DCA_{p}");
		sprintf(int_Chi_pion_true,"Chi_pion_true_%d",iter_cent);
		Chi_pion_true[iter_cent] = new TH1D(int_Chi_pion_true,int_Chi_pion_true,100, 0.0, 10000.);
		Chi_pion_true[iter_cent]->SetYTitle("Entries");
		Chi_pion_true[iter_cent]->SetXTitle("#chi^{2}_{#pi}");
		sprintf(int_Chi_proton_true,"Chi_proton_true_%d",iter_cent);
		Chi_proton_true[iter_cent] = new TH1D(int_Chi_proton_true,int_Chi_proton_true,100, 0.0, 10000.);
		Chi_proton_true[iter_cent]->SetYTitle("Entries");
		Chi_proton_true[iter_cent]->SetXTitle("#chi^{2}_{p}");
		sprintf(int_Dca_lambda_true,"Dca_lambda_true_%d",iter_cent);
		Dca_lambda_true[iter_cent] = new TH1D(int_Dca_lambda_true,int_Dca_lambda_true,100, 0.0, 200.);
		Dca_lambda_true[iter_cent]->SetYTitle("Entries");
		Dca_lambda_true[iter_cent]->SetXTitle("DCA_{#Lambda}");
		sprintf(int_Chi_lambda_true,"Chi_lambda_true_%d",iter_cent);
		Chi_lambda_true[iter_cent] = new TH1D(int_Chi_lambda_true,int_Chi_lambda_true,100, 0.0, 10000.);
		Chi_lambda_true[iter_cent]->SetYTitle("Entries");
		Chi_lambda_true[iter_cent]->SetXTitle("#chi^{2}_{#Lambda}");
		sprintf(int_Dca_v0_true,"Dca_v0_true_%d",iter_cent);
		Dca_v0_true[iter_cent] = new TH1D(int_Dca_v0_true,int_Dca_v0_true,100, 0.0, 140.);
		Dca_v0_true[iter_cent]->SetYTitle("Entries");
		Dca_v0_true[iter_cent]->SetXTitle("DCA_{V_{0}}");
		sprintf(int_Chi_v0_true,"Chi_v0_true_%d",iter_cent);
		Chi_v0_true[iter_cent] = new TH1D(int_Chi_v0_true,int_Chi_v0_true,100, 0.0, 25.);
		Chi_v0_true[iter_cent]->SetYTitle("Entries");
		Chi_v0_true[iter_cent]->SetXTitle("DCA_{V_{0}}");
		sprintf(int_Path_hist_true,"Path_hist_true_%d",iter_cent);
		Path_hist_true[iter_cent] = new TH1D(int_Path_hist_true,int_Path_hist_true,100, 0.0, 300.);
		Path_hist_true[iter_cent]->SetYTitle("Entries");
		Path_hist_true[iter_cent]->SetXTitle("Path");
		sprintf(int_Angle_hist_true,"Angle_hist_true_%d",iter_cent);
		Angle_hist_true[iter_cent] = new TH1D(int_Angle_hist_true,int_Angle_hist_true,100, 0.0, 1.6);
		Angle_hist_true[iter_cent]->SetYTitle("Entries");
		Angle_hist_true[iter_cent]->SetXTitle("Path");
	}
	
	Int_t full_events = chain.GetEntries();
	cout << " Full number of Events = " << full_events << "; Chosen number of events = " << n_ev << endl;
	if (n_ev != 0) full_events = TMath::Min (full_events, n_ev);
	cout << " Number of events = " << full_events << endl;
	
	Int_t zero_events = 0; //for check
	for(Int_t iev = 0; iev < full_events; ++iev) //cycle for entries
	{  
		chain.GetEntry(iev); // read event
		//cout << " Event: " << iev << "\t\t\t\t\r"<<flush;
		if (iev%1000 == 0) std::cout << "Event [" << iev << "/" << full_events << "]" << std::endl;
		
		Centrality_tpc = GetCentrality(n_tracks_mpd);
		if(n_tracks_mpd == 0) 
		{ 
			zero_events++; 
			continue;
		}
//this one is just in case, the previous one takes care of all the "empty events"		
		if(Centrality_tpc < 0 || Centrality_tpc > 100) 
		{
			zero_events++; 
			continue;
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
		
		L0 *lamb = nullptr;
		
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			if(Centrality_tpc < centrality_max[iter_cent] && Centrality_tpc >= centrality_min[iter_cent])
			{	
				for(UInt_t i = 0; i < l0->size(); ++i) //cycle for reco l0
				{
					lamb = (L0*) &l0->at(i);
					Dca_pion[iter_cent]->Fill(lamb->dcas[0]);
					Dca_proton[iter_cent]->Fill(lamb->dcas[1]);
					Chi_pion[iter_cent]->Fill(lamb->chi2s[0]);
					Chi_proton[iter_cent]->Fill(lamb->chi2s[1]);
					Dca_lambda[iter_cent]->Fill(lamb->dca);
					Chi_lambda[iter_cent]->Fill(lamb->c2pv);
					Dca_v0[iter_cent]->Fill(lamb->disth);
					Chi_v0[iter_cent]->Fill(lamb->chi2h);
					Path_hist[iter_cent]->Fill(lamb->path);
					Angle_hist[iter_cent]->Fill(lamb->angle);	
						
					if (lamb->origs[0] > 0)
					{
						Dca_pion_true[iter_cent]->Fill(lamb->dcas[0]);
						Dca_proton_true[iter_cent]->Fill(lamb->dcas[1]);
						Chi_pion_true[iter_cent]->Fill(lamb->chi2s[0]);
						Chi_proton_true[iter_cent]->Fill(lamb->chi2s[1]);
						Dca_lambda_true[iter_cent]->Fill(lamb->dca);
						Chi_lambda_true[iter_cent]->Fill(lamb->c2pv);
						Dca_v0_true[iter_cent]->Fill(lamb->disth);
						Chi_v0_true[iter_cent]->Fill(lamb->chi2h);
						Path_hist_true[iter_cent]->Fill(lamb->path);
						Angle_hist_true[iter_cent]->Fill(lamb->angle);
					}
    
				}
			}
		}
    
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Dca_pion[iter_cent]->Write();
		Dca_proton[iter_cent]->Write();
		Chi_pion[iter_cent]->Write();
		Chi_proton[iter_cent]->Write();
		Dca_lambda[iter_cent]->Write();
		Chi_lambda[iter_cent]->Write();
		Dca_v0[iter_cent]->Write();
		Chi_v0[iter_cent]->Write();
		Path_hist[iter_cent]->Write();
		Angle_hist[iter_cent]->Write();
		
		Dca_pion_true[iter_cent]->Write();
		Dca_proton_true[iter_cent]->Write();
		Chi_pion_true[iter_cent]->Write();
		Chi_proton_true[iter_cent]->Write();
		Dca_lambda_true[iter_cent]->Write();
		Chi_lambda_true[iter_cent]->Write();
		Dca_v0_true[iter_cent]->Write();
		Chi_v0_true[iter_cent]->Write();
		Path_hist_true[iter_cent]->Write();
		Angle_hist_true[iter_cent]->Write();
	}
	
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
