//creates distribution of invariant mass and performs the fitting procedure (Legendre polinoms + Gaus for Signal), then separate background fit in the sidebands, then cutting off the background, for the full phase space, with different selection criteria parameters to find the optimal set of values
//writes the values to the txt file (for each centrality bin)

#include <TBranch.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TTree.h>
#include <Riostream.h>
#include <TLatex.h>

#include <iostream>
using namespace std;

void FitFile(int NITER, TString fullpath, Double_t xmin, Double_t nsig_bckg, Double_t nsig, double &max_value_file, int &max_value_iter_V0_file, int &max_value_iter_path_file, int &max_value_iter_angle_file);

//fit function for background (Legendre polinoms (L0 - L4))
Double_t background(Double_t *x, Double_t *par) {
	Double_t vxmix = 1.07;
	Double_t vxmax = 1.17;
	Double_t e = 2.*(x[0] - vxmax)/(vxmax - vxmix) - 1.;
	return par[0]*(1. + par[1]*e + par[2]*0.5*(3.*TMath::Power(e,2) - 1.) + par[3]*0.5*e*(5.*TMath::Power(e,2) - 3.) + par[4]*0.125*(35.*TMath::Power(e,4) - 30.*TMath::Power(e,2) + 3.));
}

//fit function for signal (gaus)
Double_t gaussian(Double_t *x, Double_t *par) {
   return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2)));
}

// Function for the fit of background and signal
Double_t fitting_function(Double_t *x, Double_t *par) {
	Double_t vxmix = 1.07;
	Double_t vxmax = 1.17;
	Double_t e = 2.*(x[0] - vxmax)/(vxmax - vxmix) - 1.;
	return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*e + par[5]*0.5*(3.*TMath::Power(e,2) - 1.) + par[6]*0.5*e*(5.*TMath::Power(e,2) - 3.) + par[7]*0.125*(35.*TMath::Power(e,4) - 30.*TMath::Power(e,2) + 3.));
}

void LegendreFit_selections_train_MB(TString inFileDir = "Directory_Path", TString outname = "Selections_Full.root", int NITER = 10, double xmin = 1.086, double nsig_bckg = 7.0, double nsig = 4.0)
{
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	
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
	
	double max_value_final = 0.;
	int max_value_iter_pi_final = 0;
	int max_value_iter_p_final = 0;
	int max_value_iter_V0_final = 0;
	int max_value_iter_path_final = 0;
	int max_value_iter_angle_final = 0;
	
	//finding the directory and the files in it, then running on each file the fitting routine to find the best value of significance:
	
	for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
	{
		for(int iter_p = 0; iter_p < NITER; iter_p++)
		{
			TString outdir = Form("chipi%d_chip%d",iter_pi,iter_p); 
			TString fullpath;
			fullpath = inFileDir;
			fullpath += '/';
			fullpath += outdir;
			fullpath += '/';
			fullpath += outname;
			
			double max_value_file = 0.;
			int max_value_iter_V0_file = 0;
			int max_value_iter_path_file = 0;
			int max_value_iter_angle_file = 0;

			FitFile(NITER,fullpath,xmin,nsig_bckg,nsig,max_value_file,max_value_iter_V0_file,max_value_iter_path_file,max_value_iter_angle_file);	
			
			//now need to compare all the obtained values (best values from each directory/file), to find optimal values for each centrality bin:
			
			if(max_value_file > max_value_final)
			{
				max_value_final = max_value_file;
				max_value_iter_pi_final = iter_pi;
				max_value_iter_p_final = iter_p;
				max_value_iter_V0_final = max_value_iter_V0_file;
				max_value_iter_path_final = max_value_iter_path_file;
				max_value_iter_angle_final = max_value_iter_angle_file;
			}
			
			cout << "File = " << fullpath << "; Highest SSB ratio = " << max_value_file << "; max_value_iter_V0 = " << max_value_iter_V0_file << "; max_value_iter_path = " << max_value_iter_path_file << "; max_value_iter_angle = " << max_value_iter_angle_file <<  endl;			
		}
	}
	
	ofstream output_file_selections;
	output_file_selections.open ("ChiSelection_values_MB.txt");
	
	ofstream output_iterations_selections;
	output_iterations_selections.open ("ChiSelection_iterations_MB.txt");
	
	//Final values
	cout << "Highest SSB ratio = " << max_value_final << "; max_value_iter_pi = " << max_value_iter_pi_final << "; max_value_iter_p = " << max_value_iter_p_final << "; max_value_iter_V0 = " << max_value_iter_V0_final << "; max_value_iter_path = " << max_value_iter_path_final << "; max_value_iter_angle = " << max_value_iter_angle_final <<  endl;
	cout << "Values: " << "; chi_pi_value = " << chi_pi_value[max_value_iter_pi_final] << "; chi_p_value = " << chi_p_value[max_value_iter_p_final] << "; chi_V0_value = " << chi_V0_value[max_value_iter_V0_final] << "; lambda_path_value = " << lambda_path_value[max_value_iter_path_final] << "; lambda_angle_value = " << lambda_angle_value[max_value_iter_angle_final] <<  endl;
	
	output_file_selections << chi_pi_value[max_value_iter_pi_final] << "  " << chi_p_value[max_value_iter_p_final] << "  " << chi_V0_value[max_value_iter_V0_final] << "  " << lambda_path_value[max_value_iter_path_final] << "  " << lambda_angle_value[max_value_iter_angle_final] << "\n"; //write to file
	
	output_iterations_selections << max_value_iter_pi_final << "  " << max_value_iter_p_final << "  " << max_value_iter_V0_final << "  " << max_value_iter_path_final << "  " << max_value_iter_angle_final << "\n"; //write to file
	if(max_value_iter_pi_final == NITER-1 || max_value_iter_pi_final == 0) 
	{
		cout << "chi_pi value at the end of range! " << "; iter = " << max_value_iter_pi_final << "; value = " << chi_pi_value[max_value_iter_pi_final] << endl;
	}
	if(max_value_iter_p_final == NITER-1 || max_value_iter_p_final == 0) 
	{
		cout << "chi_p value at the end of range! " << "; iter = " << max_value_iter_p_final << "; value = " << chi_p_value[max_value_iter_p_final] << endl;
	}
	if(max_value_iter_V0_final == NITER-1 || max_value_iter_V0_final == 0) 
	{
		cout << "chi_V0 value at the end of range! " << "; iter = " << max_value_iter_V0_final << "; value = " << chi_V0_value[max_value_iter_V0_final] << endl;
	}
	if(max_value_iter_path_final == NITER-1 || max_value_iter_path_final == 0) 
	{
		cout << "path value at the end of range! " << "; iter = " << max_value_iter_path_final << "; value = " << lambda_path_value[max_value_iter_path_final] << endl;
	}
	if(max_value_iter_angle_final == NITER-1 || max_value_iter_angle_final == 0) 
	{
		cout << "angle value at the end of range! " << "; iter = " << max_value_iter_angle_final << "; value = " << lambda_angle_value[max_value_iter_angle_final] << endl;
	}
	
	output_file_selections.close();
	output_iterations_selections.close();
	
}
//Function which opens a file, performs the fitting of all the histograms inside, then outputs the corresponding best values for significance, and the values of iterations (which file it was)
void FitFile(int NITER, TString fullpath, Double_t xmin, Double_t nsig_bckg, Double_t nsig, double &max_value_file, int &max_value_iter_V0_file, int &max_value_iter_path_file, int &max_value_iter_angle_file)
{			
	TH1D *hm0[NITER][NITER][NITER], *hm0_signal[NITER][NITER][NITER], *hm0_bckg[NITER][NITER][NITER];
	TF1 *fitting_fnc[NITER][NITER][NITER], *backFcn[NITER][NITER][NITER], *signalFcn[NITER][NITER][NITER];
	char *int_fitting_fnc = new char[100];
	char *int_backFcn = new char[100];
	char *int_signalFcn = new char[100];

	double sum_sig[NITER][NITER][NITER];
	double sum_full[NITER][NITER][NITER];
	double entries_new[NITER][NITER][NITER];
	double entries_old[NITER][NITER][NITER];
	double efficiency[NITER][NITER][NITER];
	double ratio_SB[NITER][NITER][NITER];
	double ratio_SSB[NITER][NITER][NITER];

	double xmax[NITER][NITER][NITER];
	int bin_left[NITER][NITER][NITER];
	int bin_right[NITER][NITER][NITER];

	//TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);

	TFile *myFile_data = new TFile(fullpath);
	
	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				hm0[iter_V0][iter_path][iter_angle] = (TH1D*) myFile_data->Get(Form("hm0_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle));
				hm0_signal[iter_V0][iter_path][iter_angle] = (TH1D*)hm0[iter_V0][iter_path][iter_angle]->Clone("hm0sig_V%d_path%d_angle%d");
				hm0_bckg[iter_V0][iter_path][iter_angle] = (TH1D*)hm0[iter_V0][iter_path][iter_angle]->Clone("hm0bckg_V%d_path%d_angle%d");
				xmax[iter_V0][iter_path][iter_angle] = hm0[iter_V0][iter_path][iter_angle]->GetXaxis()->GetXmax();
			}
		}
	}
	
	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				sprintf(int_fitting_fnc,"int_fitting_fnc_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle);
				fitting_fnc[iter_V0][iter_path][iter_angle] = new TF1(int_fitting_fnc,fitting_function,xmin,xmax[iter_V0][iter_path][iter_angle],8);
			
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetNpx(500);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetLineWidth(4);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetLineColor(2); //red color for fitting line
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
				Int_t npar = fitting_fnc[iter_V0][iter_path][iter_angle]->GetNumberFreeParameters();
				for (Int_t ip = 0; ip < npar; ++ip) fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(ip,0);	

				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(0,10000);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(1,1.116);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(2,0.002);
							
				//change to the canvas:
				//c1->cd();	
				//c1->Clear();
				hm0[iter_V0][iter_path][iter_angle]->SetMinimum(1.0E-20);
				hm0[iter_V0][iter_path][iter_angle]->GetXaxis()->SetRangeUser(1.07,1.17);
				hm0[iter_V0][iter_path][iter_angle]->Fit(int_fitting_fnc,"wq0","",1.086,1.18);
				//hm0[iter_V0][iter_path][iter_angle]->Fit(int_fitting_fnc,"qn","",1.086,1.18);
							
				Float_t mass = fitting_fnc[iter_V0][iter_path][iter_angle]->GetParameter(npar-7);
				Float_t sigma = TMath::Abs(fitting_fnc[iter_V0][iter_path][iter_angle]->GetParameter(npar-6));
						
				Double_t xmin_cut = mass - nsig_bckg * sigma;
				Double_t xmax_cut = mass + nsig_bckg * sigma;
						
				Int_t imin = hm0[iter_V0][iter_path][iter_angle]->FindBin(xmin_cut);
				Int_t imax = hm0[iter_V0][iter_path][iter_angle]->FindBin(xmax_cut);
							
				//improve the picture:
				sprintf(int_backFcn,"backFcn_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle);
				backFcn[iter_V0][iter_path][iter_angle] = new TF1(int_backFcn,background,xmin,xmax[iter_V0][iter_path][iter_angle],5);
				backFcn[iter_V0][iter_path][iter_angle]->SetLineColor(kMagenta);
							
				Double_t errs[200] = {0};
				Int_t nbins = hm0[iter_V0][iter_path][iter_angle]->GetNbinsX();
				for (Int_t ib = 1; ib <= nbins; ++ib) 
				{
					if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg[iter_V0][iter_path][iter_angle]->GetBinContent(ib));
				}
				hm0_bckg[iter_V0][iter_path][iter_angle]->SetError(errs);
				hm0_bckg[iter_V0][iter_path][iter_angle]->SetLineColor(kMagenta);
				hm0_bckg[iter_V0][iter_path][iter_angle]->Fit(int_backFcn,"q0","same",1.086,1.17);
				//hm0_bckg[iter_V0][iter_path][iter_angle]->Fit(int_backFcn,"qn","same",1.086,1.17);
							
				sprintf(int_signalFcn,"signalFcn_cut1_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle);
				signalFcn[iter_V0][iter_path][iter_angle] = new TF1(int_signalFcn,gaussian,xmin,xmax[iter_V0][iter_path][iter_angle],3);
				signalFcn[iter_V0][iter_path][iter_angle]->SetLineColor(kBlue);
				signalFcn[iter_V0][iter_path][iter_angle]->SetNpx(500);
							
				hm0_signal[iter_V0][iter_path][iter_angle]->Sumw2();	
				hm0_signal[iter_V0][iter_path][iter_angle]->Add(backFcn[iter_V0][iter_path][iter_angle], -1);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetLineColor(kBlack);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetLineWidth(2);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetMarkerStyle(2);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetMarkerSize(1);
				signalFcn[iter_V0][iter_path][iter_angle]->SetParameters(hm0[iter_V0][iter_path][iter_angle]->GetMaximum(),1.1157,0.003);
				hm0_signal[iter_V0][iter_path][iter_angle]->Fit(int_signalFcn,"q0","same",1.105,1.125);
				//hm0_signal[iter_V0][iter_path][iter_angle]->Fit(int_signalFcn,"qn","same",1.105,1.125);
							
				mass = signalFcn[iter_V0][iter_path][iter_angle]->GetParameter(1);
				sigma = signalFcn[iter_V0][iter_path][iter_angle]->GetParameter(2);
							
				Double_t left = mass - nsig * sigma;
				Double_t right = mass + nsig * sigma;
				left = TMath::Max(left,1.105);
				right = TMath::Min(right,1.13);
							
				bin_left[iter_V0][iter_path][iter_angle] = hm0_signal[iter_V0][iter_path][iter_angle]->FindBin(left);
				bin_right[iter_V0][iter_path][iter_angle] = hm0_signal[iter_V0][iter_path][iter_angle]->FindBin(right);
							
				entries_new[iter_V0][iter_path][iter_angle] = hm0_signal[iter_V0][iter_path][iter_angle]->Integral(bin_left[iter_V0][iter_path][iter_angle],bin_right[iter_V0][iter_path][iter_angle]);
				entries_old[iter_V0][iter_path][iter_angle] = hm0[iter_V0][iter_path][iter_angle]->GetEntries();			
			}
		}
	}
			
	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				sum_sig[iter_V0][iter_path][iter_angle] = 0.;
				sum_full[iter_V0][iter_path][iter_angle] = 0.;
			}
		}
	}
	
	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				efficiency[iter_V0][iter_path][iter_angle] = 100.*entries_new[iter_V0][iter_path][iter_angle]/entries_old[iter_V0][iter_path][iter_angle];
			
				for (Int_t ib = bin_left[iter_V0][iter_path][iter_angle]; ib <= bin_right[iter_V0][iter_path][iter_angle]; ++ib) 
				{
					sum_sig[iter_V0][iter_path][iter_angle] += hm0_signal[iter_V0][iter_path][iter_angle]->GetBinContent(ib); 
					sum_full[iter_V0][iter_path][iter_angle] += hm0[iter_V0][iter_path][iter_angle]->GetBinContent(ib);
				}
				ratio_SB[iter_V0][iter_path][iter_angle] = sum_sig[iter_V0][iter_path][iter_angle]/(sum_full[iter_V0][iter_path][iter_angle] - sum_sig[iter_V0][iter_path][iter_angle]);
				ratio_SSB[iter_V0][iter_path][iter_angle] = sum_sig[iter_V0][iter_path][iter_angle]/TMath::Sqrt(sum_full[iter_V0][iter_path][iter_angle]);
			}
		}
	}
	
	//find the highest values for each centrality interval and the corresponding values of selections parameters:
	double max_value = ratio_SSB[0][0][0];
	int max_value_iter_V0 = 0;
	int max_value_iter_path = 0;
	int max_value_iter_angle = 0;
	
	for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
	{
		for(int iter_path = 0; iter_path < NITER; iter_path++)
		{
			for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
			{
				if(ratio_SSB[iter_V0][iter_path][iter_angle] > max_value)
				{
					max_value = ratio_SSB[iter_V0][iter_path][iter_angle];
					max_value_iter_V0 = iter_V0;
					max_value_iter_path = iter_path;
					max_value_iter_angle = iter_angle;
				}
			}
		}
	}
	
	max_value_file = max_value;
	max_value_iter_V0_file = max_value_iter_V0;
	max_value_iter_path_file = max_value_iter_path;
	max_value_iter_angle_file = max_value_iter_angle;

	cout << "File = " << fullpath << "; Highest SSB ratio = " << max_value << endl;
	/*
	if (c1) 
	{ 
		c1->Close(); 
		delete c1;
		c1 = 0; 
	}*/
	myFile_data->Close();
}
