///Using the calculated tree with impact parameter and multiplicity from the dataset, construct the distributions of model's average impact parameter in each centrality region
///Provide the choice of DCA (with/without) and the paths for input files, centrality file and output file
//start with (e.g.): make Impact_par_from_tpc_centrality && ./Impact_par_from_tpc_centrality -dca 0 -inname /scratch2/nazarova/phsd_output/*.root -centname /scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_noDCA/FINAL.root -outname /scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_noDCA/B_average_model.root

#include <iostream>
#include <vector>

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
#include <TChain.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TGraphErrors.h>
using namespace std;

const int NITER_CENT = 10; //amount of centrality intervals

Float_t GetCentrality(Int_t multiplicity);

TRandom *RNG = new TRandom();
Int_t Border_max[NITER_CENT];
Int_t Border_min[NITER_CENT];
Float_t Centrality_bin[NITER_CENT];

TH1D *Centrality[NITER_CENT], *Vertex_Z[NITER_CENT];
char *cent_int = new char[NITER_CENT];
char *Vertex_Z_int = new char[NITER_CENT];
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
int main(int argc, char** argv)
{
	TLatex latex;
	latex.SetNDC();	
	
	const int centrality_min[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
	const int centrality_max[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	const float _CentrBins[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
	const double noErr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		
	double centrality_bin[NITER_CENT]; 
	
	TString inname = "/scratch2/nazarova/phsd_output/*.root";
	TChain chain("event");
	TString outname = ""; 
	TString centname = "";
	Int_t n2 = 0; 
	Int_t dca_choice = 0; // default (no DCA)
	
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
		if (tmp.compare("-dca") == 0 && i+1 < argc) 
		{
			float number = atof(argv[i+1]);
			if (number > 0) 
			{
				dca_choice = number;
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
		if (tmp.compare("-centname") == 0 && i+1 < argc) 
		{
			centname = argv[i+1];
		}
		if (tmp.compare("-outname") == 0 && i+1 < argc) 
		{
			outname = argv[i+1];
		}
	}
	
	if(centname.Length() == 0)
	{
		cerr << "!! Please provide the path for the centrality file !!" << endl;
		return 1;
	}
	
	if(outname.Length() == 0)
	{
		cerr << "!! Please provide the path for the output file !!" << endl;
		return 1;
	}
	
	cout << "events chosen = " << n2 << endl;
	
	if(dca_choice == 0)  {cout << " Without using DCA cut" << endl;}
	if(dca_choice == 1)  {cout << " Using DCA cut" << endl;}
	
	Int_t Ncc, MinBorder, MaxBorder;
	Float_t MinPercent, MaxPercent;
	TFile *file_centrality = new TFile(centname);
	//reading the centrality tree
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
	//if(Border_min[NITER_CENT-1] > 100)
	
	//Border_min[NITER_CENT-1] = 1;
	//Border_max[NITER_CENT-1] = 2;
	
	//Border_min[NITER_CENT-2] = 2;
	
	//Border_max[NITER_CENT-1] = 4;
	for (Int_t i = 0; i < NITER_CENT; ++i)
	{	
		cout << "Border_min = " << Border_min[i] << "; Border_max = " << Border_max[i] << "; Centrality_bin = " << Centrality_bin[i] << endl;
	}
	//return(0);
	//variables:
	Float_t b0, Z_vertex, Centrality_tpc;
	Int_t n_tracks_mpd, n_tracks_mpd_DCA;
	
	//set adress to the branches of the input tree:
	chain.SetBranchAddress("b0", &b0);
	chain.SetBranchAddress("n_tracks_mpd", &n_tracks_mpd);
	chain.SetBranchAddress("n_tracks_mpd_DCA", &n_tracks_mpd_DCA);
	chain.SetBranchAddress("Z_vertex", &Z_vertex);
	
	//the necessary histograms for analysis/tests
	TH1F *NCentr = new TH1F("NCentr","NCentr",NITER_CENT,_CentrBins);
	NCentr->SetYTitle("Entries");
	NCentr->SetXTitle("Centrality TPC");
	NCentr->SetLineColor(kBlack);
	
	TH1D* B_average_VS_Centrality = new TH1D("B_average_VS_Centrality", "B_average_VS_Centrality", NITER_CENT,0,100);
	TH2D *hVertex_Z_vs_Cent = new TH2D("hVertex_Z_vs_Cent","hVertex_Z_vs_Cent",100,-150,150,NITER_CENT,0,100);
	TH2D *hRecoVertex_Z_vs_Cent = new TH2D("hRecoVertex_Z_vs_Cent","hRecoVertex_Z_vs_Cent",100,-200,200,NITER_CENT,0,1600);
	
	TH1D* B_plot = new TH1D("B_plot", "B_plot", 200,0,20);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 468,0.,468);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 481,0.,481);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 471,0.,471.9);
	
	//for b(0,12), DCA<1.0, cut 30
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 380,0.,380.9);
	
	//for b(0,16), DCA<1.0, cut 25, cut 30
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 396,0.,396.5);
	
	//for b(0,16), DCA<0.5, cut 25, cut 30
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 315,0.,315.9);	
	
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 380,0.,380.9);
	
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 396,0.,396.5);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 674,0.,674.7);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 442,0.,442);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 748,0.,748.8);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 670,0.,670.0);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 421,0.,421.2);
	//TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 739,0.,739.7);
	TH1D *Multiplicity_rescaled = new TH1D("Multiplicity_rescaled", "Multiplicity_rescaled", 664,0.,664.3);
	
	TFile out(outname,"recreate");
	
	for(int iter = 0; iter < NITER_CENT; iter++)
	{
		centrality_bin[iter] = centrality_min[iter] + (centrality_max[iter] - centrality_min[iter])/2;
		
		sprintf(cent_int,"Centrality_%d",iter);                               
		Centrality[iter] = new TH1D(cent_int,cent_int,100, 0., 20.);
		Centrality[iter]->SetYTitle("Entries");
		Centrality[iter]->SetXTitle("b, fm");
		Centrality[iter]->SetLineWidth(2);
		
		sprintf(Vertex_Z_int,"Vertex_Z_%d",iter);                               
		Vertex_Z[iter] = new TH1D(Vertex_Z_int,Vertex_Z_int,100,-150,150);
		Vertex_Z[iter]->SetYTitle("Entries");
		Vertex_Z[iter]->SetXTitle("b, fm");
		Vertex_Z[iter]->SetLineWidth(2);	
	}
	
	Int_t full_events = chain.GetEntries();
	cout << " Number of events = " << full_events << " Number n2 = " << n2 << endl;
	if (n2 != 0) full_events = TMath::Min (full_events, n2);
	cout << " Number of events = " << full_events << endl;

	Int_t fill_counter = 0;
	Int_t event_counter = 0;

	Int_t count_ev = 0; //for check
	Int_t zero_events = 0; //for check
	for(Int_t iev = 0; iev < full_events; ++iev) //cycle for entries
	{   	
		chain.GetEntry(iev); // read event
		if (iev%1000 == 0) std::cout << "Event [" << iev << "/" << full_events << "]" << std::endl;
		if(dca_choice == 0) 
		{
			Centrality_tpc = GetCentrality(n_tracks_mpd);
			//filling a plot for multiplicity comparison
			if(fill_counter < 500000)
			{
				Multiplicity_rescaled->Fill(n_tracks_mpd);		
			}
			fill_counter++;
			if(n_tracks_mpd == 0) 
			{ 
				zero_events++; 
				continue;
			}
		}else 
		if(dca_choice == 1)
		{
			Centrality_tpc = GetCentrality(n_tracks_mpd_DCA);
			//filling a plot for multiplicity comparison
			if(fill_counter < 500000)
			{
				Multiplicity_rescaled->Fill(n_tracks_mpd_DCA);		
			}
			fill_counter++;
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
		
//		cout << "fill_counter = " << fill_counter << endl;
		B_plot->Fill(b0);	
		
		NCentr->Fill(Centrality_tpc);
		hVertex_Z_vs_Cent->Fill(Z_vertex,Centrality_tpc);
		for(int iter = 0; iter < NITER_CENT; iter++)
		{		
			if(iter == NITER_CENT-1)
			{
				if(Centrality_tpc >= centrality_min[iter])
				{
					Centrality[iter]->Fill(b0); 
					Vertex_Z[iter]->Fill(Z_vertex); 
				}
			}else if(Centrality_tpc < centrality_max[iter] && Centrality_tpc >= centrality_min[iter])
			{
				Centrality[iter]->Fill(b0); 
				Vertex_Z[iter]->Fill(Z_vertex); 
			}
		}
		
		count_ev++;				
		
	}
	cout << " Full number of Events = " << count_ev << endl;
	
	for(int iter = 0; iter < NITER_CENT; iter++)
	{
		double B_av=Centrality[iter]->GetMean();
		double RMS_B_av=Centrality[iter]->GetRMS();
		cout << "B_av = " << B_av << " ; RMS = " << RMS_B_av << endl;
		/*
		if(iter == NITER_CENT -1)
		{
			B_average_VS_Centrality->SetBinContent(iter+1,0.);
			B_average_VS_Centrality->SetBinError(iter+1,0.);
		}else*/
		B_average_VS_Centrality->SetBinContent(iter+1,B_av);
		B_average_VS_Centrality->SetBinError(iter+1,RMS_B_av);
		Vertex_Z[iter]->Write();
		
	}
	B_average_VS_Centrality->SetLineWidth(2);
	B_average_VS_Centrality->SetMarkerSize(2);
	B_average_VS_Centrality->SetLineColor(kGreen);
	B_average_VS_Centrality->SetMarkerColor(kGreen);
	B_average_VS_Centrality->SetMarkerStyle(22);
	B_average_VS_Centrality->Write();
	NCentr->Write();
	hVertex_Z_vs_Cent->Write();
	B_plot->Write();
	Multiplicity_rescaled->Write();
	
	cout << "zero_events = " << zero_events << endl;
	
	out.Close();
	
}
//get centrality from TPC (through track multiplicity)
Float_t GetCentrality(Int_t multiplicity)
{
	Float_t centrality = -1.;
	
	for (int i = 0; i < NITER_CENT; ++i)
	{
		if(i == 0)
		{
			if(multiplicity >= Border_min[i])
			{
				centrality = Centrality_bin[i];
				Double_t random_number_new = RNG->Uniform(Centrality_bin[i]-4.9,Centrality_bin[i]+4.9);
				centrality = random_number_new;
			}
		/*}else if(i == NITER_CENT-1)
		{
			if(multiplicity >= Border_min[i] && multiplicity < Border_max[i])
			{			
				centrality = Centrality_bin[i];
				Double_t random_number_new = RNG->Uniform(Centrality_bin[i]-4.9,Centrality_bin[i]+4.9);
				centrality = random_number_new;
				//cout << "Cent_int = " << i << " ; multiplicity = " << multiplicity << " ; centrality = " << centrality << endl;
			}*/
		}else 
		if(multiplicity >= Border_min[i] && multiplicity < Border_max[i])
		{
			centrality = Centrality_bin[i];
			Double_t random_number_new = RNG->Uniform(Centrality_bin[i]-4.9,Centrality_bin[i]+4.9);
			centrality = random_number_new;
		}
	}
	return centrality;
}
/*
Float_t GetCentrality(Int_t multiplicity)
{
	Float_t centrality = -1.;
	
	for (int i = 0; i < NITER_CENT; ++i)
	{
		if(multiplicity > Border_min[i] && multiplicity <= Border_max[i])
		{
			centrality = Centrality_bin[i];
			Double_t random_number = RNG->Uniform(Centrality_bin[i]-5,Centrality_bin[i]+5);
			centrality = random_number;
		}
	}
	return centrality;
}
*/
