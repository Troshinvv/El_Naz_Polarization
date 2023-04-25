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
	int NITER = 10; // amount of test for variables for selection: - check!!!
	int NITER_CENT = 4; // amount of centrality intervals (for the analysis)	
	int n_ev = 0; // number of events to analyze (default 0 - means all)
	double cent_cut = 70.; // value for centrality cut
	Int_t cent_cut_choice = 0; // choice of centrality cut (0 - no cent cut, 1 - with cent cut)

	//starting values for variables for selection: - check!!!
	double chi_pi_start = 0.; // starting value for chi_pi
	double chi_p_start = 0.; // starting value for chi_p
	double chi_V0_start = 0.; // starting value for chi_V0
	double lambda_path_start = 0.; // starting value for lambda_path
	double lambda_angle_start = 0.; // starting value for lambda_angle
	
	//step values for variables for selection: - check!!!
	double chi_pi_step = 1.; // step value for chi_pi
	double chi_p_step = 1.; // step value for chi_p
	double chi_V0_step = 1.; // step value for chi_V0
	double lambda_path_step = 1.; // step value for lambda_path
	double lambda_angle_step = 0.1; // step value for lambda_angle
	
	
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
	
	double chi_pi_value[NITER];
	double chi_p_value[NITER];
	double chi_V0_value[NITER];
	double lambda_path_value[NITER];
	double lambda_angle_value[NITER];
	
	for(int iter = 0; iter < NITER; iter++)
	{
		chi_pi_value[iter] = chi_pi_start + chi_pi_step*iter;
		chi_p_value[iter] = chi_p_start + chi_p_step*iter;
		chi_V0_value[iter] = chi_V0_start + chi_V0_step*iter;
		lambda_path_value[iter] = lambda_path_start + lambda_path_step*iter;
		lambda_angle_value[iter] = lambda_angle_start + lambda_angle_step*iter;
	}
	
	for(int iter = 0; iter < NITER; iter++)
	{
		cout << "iter = " << iter << "; chi_pi_value = " << chi_pi_value[iter] << "; chi_p_value = " << chi_p_value[iter] << "; chi_V0_value = " << chi_V0_value[iter] << "; lambda_path_value = " << lambda_path_value[iter] << "; lambda_angle_value = " << lambda_angle_value[iter] << endl;
	}
	
	//return 1; //debug
	
	TH1D *hm0_before[NITER_CENT]; // Invariant mass before selection
	TH1D *hm0[NITER_CENT][NITER][NITER][NITER][NITER][NITER]; // Invariant mass after selection
	char *int_hm0 = new char[NITER];
	char *int_hm0_before = new char[NITER_CENT];

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
		sprintf(int_hm0_before,"hm0_before_%d",iter_cent);
		hm0_before[iter_cent] = new TH1D(int_hm0_before,int_hm0_before,100, 1.07, 1.17);
		hm0_before[iter_cent]->SetYTitle("Entries");
		hm0_before[iter_cent]->SetXTitle("M_{inv}, GeV/c^{2}");
		
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							sprintf(int_hm0,"hm0_cent%d_pi%d_p%d_V%d_path%d_angle%d",iter_cent,iter_pi,iter_p,iter_V0,iter_path,iter_angle);
							
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = new TH1D(int_hm0,int_hm0,100, 1.07, 1.17);
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetTitle("Lmass after selection");
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetYTitle("Entries");
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetXTitle("M_{inv}, GeV/c^{2}");
						}
					}
				}
			}
		}
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
					hm0_before[iter_cent]->Fill(lamb->massh);
					
					for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
					{
						for(int iter_p = 0; iter_p < NITER; iter_p++)
						{
							for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
							{
								for(int iter_path = 0; iter_path < NITER; iter_path++)
								{
									for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
									{
										if(lamb->chi2s[0] > chi_pi_value[iter_pi] && lamb->chi2s[1] > chi_p_value[iter_p] && lamb->chi2h < chi_V0_value[iter_V0] && lamb->path > lambda_path_value[iter_path] && lamb->angle < lambda_angle_value[iter_angle])
										{
											hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Fill(lamb->massh);
										}
									}
								}
							}
						}
					}    
				}
			}
		}
    
	}

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		hm0_before[iter_cent]->Write();
		
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Write();
						}
					}
				}
			}
		} 
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
