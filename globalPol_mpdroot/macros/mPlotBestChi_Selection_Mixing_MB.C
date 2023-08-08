/**
 * @file mPlotBestChi_Selection_SidebandsFit_MB.C
 * @author Elizaveta Nazarova
 * @brief Plotting macro for RECO GlobalPol wagon, selection mode, chi selection
 * @version 0.1
 * @date 2023-07-15
 * Plot the performance of the event mixing method for signal extraction for the obtained optimal values of selection parameters
 * Make sure that the input values correspond to the ones in the config file used 
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
/**
 * @brief Plotting macro for RECO GlobalPol wagon using optimal chi selection values from event mixing method for background subtraction
 * 
 * @param inFileDir            absolute path of the directory with the input files (e.g. /outdir)
 * @param configname           name of the config file (e.g. pGlobalPolRECO.txt)
 * @param inname               name of merged input file (e.g. Selections_Full.root)
 * @param selections_file_name name of input file with the iterations of obtained optimalselection values
 * @param NITER                number of different values of parameters for selection (default: 10)
 * @param xmin                 minimal value of invariant mass (X-axis) for fitting (default: 1.086)
 * @param nsig_bckg            value to obtain the range for sidebands for background fit around the peak (+/- nsig_bckg*sigma) (default: 7.0)
 * @param nsig                 value to obtain the range for signal extraction around the peak (+/- nsig*sigma) (default: 3.0)
 */
void mPlotBestChi_Selection_Mixing_MB(TString inFileDir = "/outdir", TString configname = "pGlobalPolRECO.txt", TString inname = "Selections_Full.root", TString selections_file_name = "Selections_iterations.txt", int NITER = 10, double xmin = 1.09, double nsig = 3.0)
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
	
	cout << "Reading selection iterations from file:" << endl;
	ifstream selections_file;
	selections_file.open(selections_file_name);
	if(selections_file.fail())
	{
		cout << "File with selection values does not exist! Exiting... " << endl;
		exit(0);
	}
	selections_file >> max_value_iter_pi_final >> max_value_iter_p_final >> max_value_iter_V0_final >> max_value_iter_path_final >> max_value_iter_angle_final;
	selections_file.close();
	cout << "max_value_iter_pi_final =  " << max_value_iter_pi_final << "; max_value_iter_p_final =  " << max_value_iter_p_final << "; max_value_iter_V0_final =  " << max_value_iter_V0_final << "; max_value_iter_path_final =  " << max_value_iter_path_final << "; max_value_iter_angle_final =  " << max_value_iter_angle_final << endl;
		
	//now lets plot the final fits:
	TH1D *hm0, *hm0_signal, *hm0_bckg;
	TF1 *backFcn, *signalFcn;

	double sum_sig = 0.;
	double sum_full = 0.;
	double entries_new;
	double entries_old;
	double efficiency;
	double ratio_SB;
	double ratio_SSB;

	double xmax;
	int bin_left;
	int bin_right;
	
	double int_signal;
	double int_background;

	TCanvas *c1 = new TCanvas("canvas_c1","canvas_c1",0,0,600,600);

	TString outdir = Form("chipi%d_chip%d",max_value_iter_pi_final,max_value_iter_p_final); 
	TString fullpath;
	fullpath = inFileDir;
	fullpath += '/';
	fullpath += outdir;
	fullpath += '/';
	fullpath += outname;
	
	TFile *myFile_data = new TFile(fullpath);
	hm0 = (TH1D*) myFile_data->Get(Form("hm0_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final));
	hm0_signal = (TH1D*)hm0->Clone("hm0sig_V%d_path%d_angle%d");
	hm0_bckg = (TH1D*) myFile_data->Get(Form("hm0_mixed_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final));
	hm0_bckg->Scale(1./5.);
	xmax = hm0->GetXaxis()->GetXmax();

	int bin_left_scale = hm0_bckg->FindBin(1.14);
	int bin_right_scale = hm0_bckg->FindBin(1.16);
	double sum_hm0_scale = 0.;
	double sum_hm0_mixed_scale = 0.;
	for (int ib = bin_left_scale; ib <= bin_right_scale; ++ib) 
	{
		sum_hm0_scale += hm0->GetBinContent(ib); 
		sum_hm0_mixed_scale += hm0_bckg->GetBinContent(ib);
	}
	hm0_bckg->Scale(sum_hm0_scale/sum_hm0_mixed_scale);	
	backFcn = new TF1(Form("backFcn_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),background,xmin,xmax,5);
	backFcn->SetLineColor(kRed);
	backFcn->SetLineWidth(4);
							
	//change to the canvas:
	c1->cd();	
	hm0->SetMinimum(1.0E-20);
	hm0->GetXaxis()->SetRangeUser(1.07,1.17);
	hm0->Draw();
	hm0_bckg->Draw("same");
	hm0_bckg->Fit(Form("backFcn_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),"q","same",xmin,xmax);
	
	signalFcn = new TF1(Form("signalFcn_cut1_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),gaussian,xmin,xmax,3);
	signalFcn->SetLineColor(kBlue);
	signalFcn->SetNpx(500);
	hm0_signal->Sumw2();	
	hm0_signal->Add(backFcn, -1);
	hm0_signal->SetLineColor(kBlack);
	hm0_signal->SetLineWidth(2);
	hm0_signal->SetMarkerStyle(2);
	hm0_signal->SetMarkerSize(1);
	signalFcn->SetParameters(hm0->GetMaximum(),1.1157,0.003);
	hm0_signal->Fit(Form("signalFcn_cut1_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),"q0","same",xmin,xmax);
				
	double mass = signalFcn->GetParameter(1);
	double sigma = signalFcn->GetParameter(2);
	double left = mass - nsig * sigma;
	double right = mass + nsig * sigma;
	
	int_signal = signalFcn->Integral(left,right);
	int_background = backFcn->Integral(left,right);
				
	bin_left = hm0_signal->FindBin(left);
	bin_right = hm0_signal->FindBin(right);
				
	entries_new = hm0_signal->Integral(bin_left,bin_right);
	entries_old = hm0->GetEntries();	

	efficiency = 100.*entries_new/entries_old;	
	
	for (int ib = bin_left; ib <= bin_right; ++ib) 
	{
		sum_sig += hm0_signal->GetBinContent(ib); 
		sum_full += hm0->GetBinContent(ib);
	}
	ratio_SB = sum_sig/(sum_full - sum_sig);
	ratio_SSB = sum_sig/TMath::Sqrt(sum_full);
	
	TLine *line_left = new TLine(left,0,left,0.6*hm0->GetMaximum());
	line_left->SetLineColor(kBlack);
	line_left->SetLineWidth(3);
	line_left->SetLineStyle(2);
	line_left->Draw("same");
	TLine *line_right = new TLine(right,0,right,0.6*hm0->GetMaximum());
	line_right->SetLineColor(kBlack);
	line_right->SetLineWidth(3);
	line_right->SetLineStyle(2);
	line_right->Draw("same");
	TLegend *legend=new TLegend(0.6,0.65,0.88,0.9);
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->SetBorderSize(0);
	legend->AddEntry(hm0,"Data","lp");
	legend->AddEntry(backFcn,"Bckg fit","l");	
	legend->AddEntry(line_left,"Cut-off (signal)","l");
	legend->Draw("same");
		
	TLatex *title_SSB = new TLatex(1.075,0.92*hm0->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB)); 
	title_SSB->Draw("same");
	TLatex *title_SB = new TLatex(1.075,0.76*hm0->GetMaximum(), Form("S/B = %.2f",ratio_SB)); 
	title_SB->Draw("same");
	TLatex *title_efficiency = new TLatex(1.075,0.66*hm0->GetMaximum(), Form("eff. = %.2f [%%]",efficiency)); 
	title_efficiency->Draw("same");
	
	cout << "File = " << fullpath << "; Highest SSB ratio = " << ratio_SSB << endl;
}
