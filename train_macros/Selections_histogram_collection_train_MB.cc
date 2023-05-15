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

void CreateFillWrite(int NITER, int n_ev, int iter_pi, int iter_p, TChain *chains, double *chi_pi_value, double *chi_p_value, double *chi_V0_value, double *lambda_path_value, double *lambda_angle_value, TString fullpath);

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
	int n_ev = 0; // number of events to analyze (default 0 - means all)
	
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
	
	cout << "events chosen = " << n_ev << "; NITER = " << NITER << endl;
	
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
	
	TLatex latex;
	latex.SetNDC();
	
//create the histograms and save then in the output file:
	
	for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
	{
		for(int iter_p = 0; iter_p < NITER; iter_p++)
		{
			TString outdir = Form("chipi%d_chip%d",iter_pi,iter_p); 
			TString fullpath;
			fullpath = outnamedir;
			fullpath += '/';
			fullpath += outdir;
			fullpath += '/';
			fullpath += outname;
			CreateFillWrite(NITER,n_ev,iter_pi,iter_p,chain, chi_pi_value, chi_p_value, chi_V0_value,lambda_path_value,lambda_angle_value,fullpath);	
		}
	}
	
}
//Function which creates output file, creates and fills histograms, writes them into the file
void CreateFillWrite(int NITER, int n_ev, int iter_pi, int iter_p, TChain *chain, double *chi_pi_value, double *chi_p_value, double *chi_V0_value, double *lambda_path_value, double *lambda_angle_value, TString fullpath)
{
	
	//set adress to l0:
	vector<MpdLambdaPol>  *l0 = 0;
	chain->SetBranchAddress("l0", &l0);  
	
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
		
		MpdLambdaPol *lamb = nullptr;
	
		for(int i = 0; i < l0->size(); ++i) //cycle for reco l0
		{
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
