///Macro to collect invariant mass histograms with different selection parameters values (chi_pi, chi_p, chi_V0, path, angle) for each centrality bin

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

#include "MpdLambdaPol.h"

void CreateFillWrite(int NITER, int n_ev, double cent_cut, Int_t cent_cut_choice, int iter_cent, int iter_pi, int iter_p, TChain *chain, int *centrality_min, int *centrality_max, double *chi_pi_value, double *chi_p_value, double *chi_V0_value, double *lambda_path_value, double *lambda_angle_value, TString fullpath);

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
	double chi_pi_start = 6.6; // starting value for chi_pi (> this)
	double chi_p_start = 3.6; // starting value for chi_p (> this)
	double chi_V0_start = 5.6; // starting value for chi_V0 (< this)
	double lambda_path_start = 1.6; // starting value for lambda_path (> this)
	double lambda_angle_start = 0.06; // starting value for lambda_angle (< this)
	
	//step values for variables for selection: - check!!!
	double chi_pi_step = 0.2; // step value for chi_pi
	double chi_p_step = 0.2; // step value for chi_p
	double chi_V0_step = 0.2; // step value for chi_V0
	double lambda_path_step = 0.2; // step value for lambda_path
	double lambda_angle_step = 0.02; // step value for lambda_angle
	
	//create the chain of input files:
	TChain *chain = new TChain("event");
	
	//strings for the output filename and centrality filename:
	TString outnamedir = ""; 
	TString outname = ""; 
	
	//checking command line arguments, in case I want different optional parameters (NITER, NITER_CENT, events (0 means all), inname, outname)
	for (int i = 0 ; i < argc ; i++) 
	{
		string tmp(argv[i]);
		if (tmp.compare("-help") == 0) 
		{
			cout << "Choose the path for the input and output files, and optional parameters, in case they differ from the default ones." << endl;
			
			cout << "Usage: ./Selections_histogram_collection_train -NITER 10 -NITER_CENT 4 -events 100 -inname $INFILE -outnamedir $OUTFILEDIR -outname $OUTFILE" << endl;
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
					chain->Add(argv[i+count+1]);
				}
				count++;
			}
			cout << "number of files = " << count << endl;
		}
		if (tmp.compare("-outnamedir") == 0 && i+1 < argc) 
		{
			outnamedir = argv[i+1];
		}
		if (tmp.compare("-outname") == 0 && i+1 < argc) 
		{
			outname = argv[i+1];
		}
	}
	
	if(outnamedir.Length() == 0)
	{
		cerr << "!! Please provide the path for the outnamedir file !!" << endl;
		return 1;
	}
	
	if(outname.Length() == 0)
	{
		cerr << "!! Please provide the path for the output file !!" << endl;
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
	
//create the histograms and save then in the output file:
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				TString outdir = Form("cent%d_chipi%d_chip%d",iter_cent,iter_pi,iter_p); 
				TString fullpath;
				fullpath = outnamedir;
				fullpath += '/';
				fullpath += outdir;
				fullpath += '/';
				fullpath += outname;
				CreateFillWrite(NITER,n_ev,cent_cut,cent_cut_choice,iter_cent,iter_pi,iter_p,chain,centrality_min, centrality_max, chi_pi_value, chi_p_value, chi_V0_value,lambda_path_value,lambda_angle_value,fullpath);	
			}
		}
	}
	
}
//Function which creates output file, creates and fills histograms, writes them into the file
void CreateFillWrite(int NITER, int n_ev, double cent_cut, int cent_cut_choice,int iter_cent, int iter_pi, int iter_p, TChain *chain, int *centrality_min, int *centrality_max, double *chi_pi_value, double *chi_p_value, double *chi_V0_value, double *lambda_path_value, double *lambda_angle_value, TString fullpath)
{
	double Centrality_tpc;
	double b0;
	
	//set adress to l0:
	vector<MpdLambdaPol>  *l0 = 0;
	chain->SetBranchAddress("l0", &l0);  
	chain->SetBranchAddress("b0", &b0);  
	chain->SetBranchAddress("Centrality_tpc", &Centrality_tpc);
	
	Int_t full_events = chain->GetEntries();
	cout << " Full number of Events = " << full_events << "; Chosen number of events = " << n_ev << endl;
	if (n_ev != 0) full_events = TMath::Min (full_events, n_ev);
	cout << " Number of events = " << full_events << endl;
	
	TFile outfile(fullpath,"recreate");
	
	TH1D *hm0[NITER][NITER][NITER]; // Invariant mass after selection
	char *int_hm0 = new char[100];
		
	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				sprintf(int_hm0,"hm0_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle);
							
				hm0[iter_V0][iter_path][iter_angle] = new TH1D(int_hm0,int_hm0,100, 1.07, 1.17);
				hm0[iter_V0][iter_path][iter_angle]->SetTitle("Lmass after selection");
				hm0[iter_V0][iter_path][iter_angle]->SetYTitle("Entries");
				hm0[iter_V0][iter_path][iter_angle]->SetXTitle("M_{inv}, GeV/c^{2}");
			}
		}
	}
	
	//now fill the necessary histograms for this output file:
	for(Int_t iev = 0; iev < full_events; ++iev) //cycle for entries
	{  
		chain->GetEntry(iev); // read event
		if (iev%1000 == 0) std::cout << "File " << fullpath << ", Event [" << iev << "/" << full_events << "]" << std::endl;
		
		//cout << "Centrality_tpc = " << Centrality_tpc << endl;
		
		if (cent_cut_choice == 0)
		{
		}else if (cent_cut_choice == 1)
		{
			if(Centrality_tpc > cent_cut) continue;
		}else
		{
			cout << "This value of cent_cut_choice is not defined!" << endl;
			exit(1);
		}
		
		MpdLambdaPol *lamb = nullptr;
	
		for(int i = 0; i < l0->size(); ++i) //cycle for reco l0
		{
			if(Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue; 
			lamb = (MpdLambdaPol*) &l0->at(i);
			if(lamb->chi2s[0] > chi_pi_value[iter_pi] && lamb->chi2s[1] > chi_p_value[iter_p])	
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							if(lamb->chi2h < chi_V0_value[iter_V0] && lamb->path > lambda_path_value[iter_path] && lamb->angle < lambda_angle_value[iter_angle])
							{
								hm0[iter_V0][iter_path][iter_angle]->Fill(lamb->massh);
							}
						}
					}
				}
			}
		}
	}

	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				hm0[iter_V0][iter_path][iter_angle]->Write();
			}
		}
	}
	outfile.Close();
}
