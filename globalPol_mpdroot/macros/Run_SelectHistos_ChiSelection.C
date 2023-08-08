/**
 * @file Run_SelectHistos_ChiSelection.C
 * @author Elizaveta Nazarova
 * @brief Macro to collect invariant mass histograms with different chi-selection parameters values
 * @version 0.1
 * @date 2023-07-15
 * Make sure that you create folders of type chipi%d_chip%d in the output directory (e.g. for NITER=10: mkdir chipi{0..9}_chip{0..9})
 */
#include <iostream> 
#include <fstream>  
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
/**
 * @brief Main function - Macro to collect invariant mass histograms with different chi-selection parameters values
 * 
 * @param inname        name of input file with tree for chi selection mode of RECO wagon (e.g. OutputTree_GlopalPolRECO.root)
 * @param configname    name of the config file (e.g. pGlobalPolRECO.txt)
 * @param outnamedir    absolute path of the directory for the output file (e.g. /outdir)
 * @param outname       name of output file (e.g. OutputChiSelection_GlopalPolRECO.root)
 * @param NITER         number of different values of parameters for selection
 * @param n_ev          number of events to analyze (default 0 - means all)
 */
void Run_SelectHistos_ChiSelection(TString inname = "OutputTree_GlopalPolRECO.root", TString configname = "pGlobalPolRECO.txt", TString outnamedir = "", TString outname = "", const int NITER = 10, const int n_ev = 0)
{
	gROOT->LoadMacro("mpdloadlibs.C");
	gROOT->ProcessLine("mpdloadlibs()");

	cout << "Reading input text file: " << configname << endl;
	double chi_pi_start, chi_p_start, chi_V0_start, lambda_path_start, lambda_angle_start; // starting values 
	double chi_pi_step, chi_p_step, chi_V0_step, lambda_path_step, lambda_angle_step; // step values

	ifstream ifs(configname);
	if (!ifs) 
	{
		cout << "File " << configname << " can not be opened" << endl;
		return 1;
	}

	string a;
	double b;
	while (ifs.good()) 
	{
		ifs >> a;
		// if comment, skip to the end of line
		if (a.find_first_of('#') == 0 || a.find_first_of("//") == 0) 
		{
			ifs.ignore(999, '\n');
			continue;
		} 
		else if (a == "chi_pi_start") 
		{
			ifs >> b;
			chi_pi_start = b;
		}
		else if (a == "chi_p_start") 
		{
			ifs >> b;
			chi_p_start = b;
		}
		else if (a == "chi_V0_start") 
		{
			ifs >> b;
			chi_V0_start = b;
		}
		else if (a == "lambda_path_start") 
		{
			ifs >> b;
			lambda_path_start = b;
		}
		else if (a == "lambda_angle_start") 
		{
			ifs >> b;
			lambda_angle_start = b;
		}
		else if (a == "chi_pi_step") 
		{
			ifs >> b;
			chi_pi_step = b;
		}
		else if (a == "chi_p_step") 
		{
			ifs >> b;
			chi_p_step = b;
		}
		else if (a == "chi_V0_step") 
		{
			ifs >> b;
			chi_V0_step = b;
		}
		else if (a == "lambda_path_step") 
		{
			ifs >> b;
			lambda_path_step = b;
		}
		else if (a == "lambda_angle_step") 
		{
			ifs >> b;
			lambda_angle_step = b;
		}
		cout << a << "\t" << b << endl;
	}
	ifs.close();
	
	// create the chain of input files:
	TChain *chain = new TChain("event");
	if (checkfile(inname))
	{
		cout << "adding to the chain: " << inname << endl;
		chain->Add(inname);
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
	
	//set adress to Lambda object:
	vector<MpdLambdaPol>  *l0 = 0;
	chain->SetBranchAddress("l0", &l0);  
	
	int full_events = chain->GetEntries();
	cout << " Full number of Events = " << full_events << "; Chosen number of events = " << n_ev << endl;
	if (n_ev != 0) full_events = TMath::Min (full_events, n_ev);
	cout << " Number of events = " << full_events << endl;
	
	TFile outfile(fullpath,"recreate");
	
	TH1D *hm0[NITER][NITER][NITER]; // Invariant mass after selection
	
	TH1D *hm0_mixed[NITER][NITER][NITER]; // Invariant mass after selection (mixed)
		
	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				hm0[iter_V0][iter_path][iter_angle] = new TH1D(Form("hm0_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),Form("hm0_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),100, 1.07, 1.17);
				hm0[iter_V0][iter_path][iter_angle]->SetTitle("Lmass after selection");
				hm0[iter_V0][iter_path][iter_angle]->SetYTitle("Entries");
				hm0[iter_V0][iter_path][iter_angle]->SetXTitle("M_{inv}, GeV/c^{2}");
							
				hm0_mixed[iter_V0][iter_path][iter_angle] = new TH1D(Form("hm0_mixed_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),Form("hm0_mixed_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),100, 1.07, 1.17);
				hm0_mixed[iter_V0][iter_path][iter_angle]->SetTitle("Lmass mixed after selection");
				hm0_mixed[iter_V0][iter_path][iter_angle]->SetYTitle("Entries");
				hm0_mixed[iter_V0][iter_path][iter_angle]->SetXTitle("M_{inv}, GeV/c^{2}");
			}
		}
	}
	
	//now fill the necessary histograms for this output file:
	for(int iev = 0; iev < full_events; ++iev) //cycle for entries
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
								if (lamb->origs[0] > -8)
									hm0[iter_V0][iter_path][iter_angle]->Fill(lamb->massh);
								if (lamb->origs[0] < -8)
									hm0_mixed[iter_V0][iter_path][iter_angle]->Fill(lamb->massh);
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
				hm0_mixed[iter_V0][iter_path][iter_angle]->Write();
			}
		}
	}
	outfile.Close();
}
