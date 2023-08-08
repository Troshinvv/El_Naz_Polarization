/**
 * @file mChi_Selection_SidebandsFit_MB.C
 * @author Elizaveta Nazarova
 * @brief Analysis macro for RECO GlobalPol wagon, selection mode, chi selection
 * @version 0.1
 * @date 2023-07-15
 * Using obtained invariant mass distributions from the "selection" mode of GlobalPolarizationRECO wagon, performs background subtraction based on sidebands method to extract the signal
 * Make sure that the input values correspond to the ones in the config file used to run the wagon!
 * Finds the optimal values of 5 selection parameters for chi selection for full MB, corresponding to the maximum significance value
 * Writes the obtained values to the .txt file
 */
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

#include <iostream> // std::cout
#include <fstream>  // std::ifstream
using namespace std;

void FitFile(int NITER, TString fullpath, double xmin, double nsig_bckg, double nsig, double &max_value_file, int &max_value_iter_V0_file, int &max_value_iter_path_file, int &max_value_iter_angle_file);

//fit function for background (Legendre polinoms (L0 - L4))
double background(double *x, double *par) {
	double vxmix = 1.07;
	double vxmax = 1.17;
	double e = 2.*(x[0] - vxmax)/(vxmax - vxmix) - 1.;
	return par[0]*(1. + par[1]*e + par[2]*0.5*(3.*TMath::Power(e,2) - 1.) + par[3]*0.5*e*(5.*TMath::Power(e,2) - 3.) + par[4]*0.125*(35.*TMath::Power(e,4) - 30.*TMath::Power(e,2) + 3.));
}

//fit function for signal (gaus)
double gaussian(double *x, double *par) {
   return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2)));
}

// Function for the fit of background and signal
double fitting_function(double *x, double *par) {
	double vxmix = 1.07;
	double vxmax = 1.17;
	double e = 2.*(x[0] - vxmax)/(vxmax - vxmix) - 1.;
	return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*e + par[5]*0.5*(3.*TMath::Power(e,2) - 1.) + par[6]*0.5*e*(5.*TMath::Power(e,2) - 3.) + par[7]*0.125*(35.*TMath::Power(e,4) - 30.*TMath::Power(e,2) + 3.));
}
/**
 * @brief Analysis macro for RECO GlobalPol wagon to obtain optimal chi selection values using sidebands method for background subtraction
 * 
 * @param inFileDir            absolute path of the directory with the input files (e.g. /outdir)
 * @param configname           name of the config file (e.g. pGlobalPolRECO.txt)
 * @param inname               name of merged input file (e.g. Selections_Full.root)
 * @param outFile_values       name of output file for saving the obtained selection values
 * @param outFile_iterations   name of output file for saving the iterations of obtained selection values
 * @param NITER                number of different values of parameters for selection (default: 10)
 * @param xmin                 minimal value of invariant mass (X-axis) for fitting (default: 1.086)
 * @param nsig_bckg            value to obtain the range for sidebands for background fit around the peak (+/- nsig_bckg*sigma) (default: 7.0)
 * @param nsig                 value to obtain the range for signal extraction around the peak (+/- nsig*sigma) (default: 3.0)
 */
void mChi_Selection_SidebandsFit_MB(TString inFileDir = "/outdir", TString configname = "pGlobalPolRECO.txt", TString inname = "Selections_Full.root", TString outFile_values = "ChiSelection_values_MB_3sigma.txt", TString outFile_iterations = "ChiSelection_iterations_MB_3sigma.txt", int NITER = 10, double xmin = 1.086, double nsig_bckg = 7.0, double nsig = 3.0)
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
	
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);

	//finding the directory and the files in it, then running on each file the fitting routine to find the highest value of significance:
	
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
			fullpath += inname;
			
			double max_value_file = 0.;
			int max_value_iter_V0_file = 0;
			int max_value_iter_path_file = 0;
			int max_value_iter_angle_file = 0;

			FitFile(NITER,fullpath,xmin,nsig_bckg,nsig,max_value_file,max_value_iter_V0_file,max_value_iter_path_file,max_value_iter_angle_file);	
			
			// now need to compare all the obtained values (best values from each directory/file), to find optimal values for each centrality bin:
			
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
	output_file_selections.open(outFile_values);
	
	ofstream output_iterations_selections;
	output_iterations_selections.open(outFile_iterations);
	
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
void FitFile(int NITER, TString fullpath, double xmin, double nsig_bckg, double nsig, double &max_value_file, int &max_value_iter_V0_file, int &max_value_iter_path_file, int &max_value_iter_angle_file)
{			
	TH1D *hm0[NITER][NITER][NITER], *hm0_signal[NITER][NITER][NITER], *hm0_bckg[NITER][NITER][NITER];
	TF1 *fitting_fnc[NITER][NITER][NITER], *backFcn[NITER][NITER][NITER], *signalFcn[NITER][NITER][NITER];

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
	
	double int_signal[NITER][NITER][NITER];
	double int_background[NITER][NITER][NITER];

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
				fitting_fnc[iter_V0][iter_path][iter_angle] = new TF1(Form("int_fitting_fnc_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),fitting_function,xmin,xmax[iter_V0][iter_path][iter_angle],8);
			
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetNpx(500);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetLineWidth(4);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetLineColor(2);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
				int npar = fitting_fnc[iter_V0][iter_path][iter_angle]->GetNumberFreeParameters();
				for (int ip = 0; ip < npar; ++ip) fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(ip,0);	

				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(0,10000);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(1,1.116);
				fitting_fnc[iter_V0][iter_path][iter_angle]->SetParameter(2,0.002);
							
				hm0[iter_V0][iter_path][iter_angle]->SetMinimum(1.0E-20);
				hm0[iter_V0][iter_path][iter_angle]->GetXaxis()->SetRangeUser(1.07,1.17);
				hm0[iter_V0][iter_path][iter_angle]->Fit(Form("int_fitting_fnc_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),"wq0","",xmin,xmax[iter_V0][iter_path][iter_angle]);
							
				double mass = fitting_fnc[iter_V0][iter_path][iter_angle]->GetParameter(npar-7);
				double sigma = TMath::Abs(fitting_fnc[iter_V0][iter_path][iter_angle]->GetParameter(npar-6));
						
				double xmin_cut = mass - nsig_bckg * sigma;
				double xmax_cut = mass + nsig_bckg * sigma;
						
				int imin = hm0[iter_V0][iter_path][iter_angle]->FindBin(xmin_cut);
				int imax = hm0[iter_V0][iter_path][iter_angle]->FindBin(xmax_cut);
							
				//improve the picture:
				backFcn[iter_V0][iter_path][iter_angle] = new TF1(Form("backFcn_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),background,xmin,xmax[iter_V0][iter_path][iter_angle],5);
				backFcn[iter_V0][iter_path][iter_angle]->SetLineColor(kMagenta);
							
				double errs[200] = {0};
				int nbins = hm0[iter_V0][iter_path][iter_angle]->GetNbinsX();
				for (int ib = 1; ib <= nbins; ++ib) 
				{
					if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg[iter_V0][iter_path][iter_angle]->GetBinContent(ib));
				}
				hm0_bckg[iter_V0][iter_path][iter_angle]->SetError(errs);
				hm0_bckg[iter_V0][iter_path][iter_angle]->SetLineColor(kMagenta);
				hm0_bckg[iter_V0][iter_path][iter_angle]->Fit(Form("backFcn_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),"q0","same",xmin,xmax[iter_V0][iter_path][iter_angle]);
							
				signalFcn[iter_V0][iter_path][iter_angle] = new TF1(Form("signalFcn_cut1_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),gaussian,xmin,xmax[iter_V0][iter_path][iter_angle],3);
				signalFcn[iter_V0][iter_path][iter_angle]->SetLineColor(kBlue);
				signalFcn[iter_V0][iter_path][iter_angle]->SetNpx(500);
							
				hm0_signal[iter_V0][iter_path][iter_angle]->Sumw2();	
				hm0_signal[iter_V0][iter_path][iter_angle]->Add(backFcn[iter_V0][iter_path][iter_angle], -1);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetLineColor(kBlack);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetLineWidth(2);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetMarkerStyle(2);
				hm0_signal[iter_V0][iter_path][iter_angle]->SetMarkerSize(1);
				signalFcn[iter_V0][iter_path][iter_angle]->SetParameters(hm0[iter_V0][iter_path][iter_angle]->GetMaximum(),1.1157,0.003);
				hm0_signal[iter_V0][iter_path][iter_angle]->Fit(Form("signalFcn_cut1_V%d_path%d_angle%d",iter_V0,iter_path,iter_angle),"q0","same",xmin,xmax[iter_V0][iter_path][iter_angle]);
							
				mass = signalFcn[iter_V0][iter_path][iter_angle]->GetParameter(1);
				sigma = signalFcn[iter_V0][iter_path][iter_angle]->GetParameter(2);
							
				double left = mass - nsig * sigma;
				double right = mass + nsig * sigma;
				
				int_signal[iter_V0][iter_path][iter_angle] = signalFcn[iter_V0][iter_path][iter_angle]->Integral(left,right);
				int_background[iter_V0][iter_path][iter_angle] = backFcn[iter_V0][iter_path][iter_angle]->Integral(left,right);
							
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
			
				for (int ib = bin_left[iter_V0][iter_path][iter_angle]; ib <= bin_right[iter_V0][iter_path][iter_angle]; ++ib) 
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
	myFile_data->Close();
}
