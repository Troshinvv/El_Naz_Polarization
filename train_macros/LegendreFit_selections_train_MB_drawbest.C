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

void LegendreFit_selections_train_MB_drawbest(TString inFileDir = "Directory_Path", TString outname = "Selections_Full.root", TString selections_file_name = "Selections_iterations.txt", int NITER = 10, double xmin = 1.086, double nsig_bckg = 7.0, double nsig = 4.0)
{
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	
	//starting values for variables for selection: - check!!!
	double chi_pi_start = 8.2; // starting value for chi_pi (> this)
	double chi_p_start = 5.2; // starting value for chi_p (> this)
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
	
	int max_value_iter_pi_final = 0;
	int max_value_iter_p_final = 0;
	int max_value_iter_V0_final = 0;
	int max_value_iter_path_final = 0;
	int max_value_iter_angle_final = 0;
	
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
	TF1 *fitting_fnc, *backFcn, *signalFcn;

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
	hm0_bckg = (TH1D*)hm0->Clone("hm0bckg_V%d_path%d_angle%d");
	xmax = hm0->GetXaxis()->GetXmax();

	fitting_fnc = new TF1(Form("int_fitting_fnc_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),fitting_function,xmin,xmax,8);

	fitting_fnc->SetNpx(500);
	fitting_fnc->SetLineWidth(4);
	fitting_fnc->SetLineColor(2); //red color for fitting line
	fitting_fnc->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
	Int_t npar = fitting_fnc->GetNumberFreeParameters();
	for (Int_t ip = 0; ip < npar; ++ip) fitting_fnc->SetParameter(ip,0);	

	fitting_fnc->SetParameter(0,10000);
	fitting_fnc->SetParameter(1,1.116);
	fitting_fnc->SetParameter(2,0.002);
				
	//change to the canvas:
	c1->cd();	
	hm0->SetMinimum(1.0E-20);
	hm0->GetXaxis()->SetRangeUser(1.07,1.17);
	hm0->Draw();
	hm0->Fit(Form("int_fitting_fnc_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),"wq0","",1.086,1.18);
				
	Float_t mass = fitting_fnc->GetParameter(npar-7);
	Float_t sigma = TMath::Abs(fitting_fnc->GetParameter(npar-6));
			
	Double_t xmin_cut = mass - nsig_bckg * sigma;
	Double_t xmax_cut = mass + nsig_bckg * sigma;
			
	Int_t imin = hm0->FindBin(xmin_cut);
	Int_t imax = hm0->FindBin(xmax_cut);
				
	//improve the picture:
	backFcn = new TF1(Form("backFcn_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),background,xmin,xmax,5);
	backFcn->SetLineColor(kRed);
	backFcn->SetLineWidth(4);
				
	Double_t errs[200] = {0};
	Int_t nbins = hm0->GetNbinsX();
	for (Int_t ib = 1; ib <= nbins; ++ib) 
	{
		if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg->GetBinContent(ib));
	}
	hm0_bckg->SetError(errs);
	hm0_bckg->SetLineColor(kMagenta);
	hm0_bckg->Fit(Form("backFcn_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),"q","same",1.086,1.17);
				
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
	hm0_signal->Fit(Form("signalFcn_cut1_V%d_path%d_angle%d",max_value_iter_V0_final,max_value_iter_path_final,max_value_iter_angle_final),"q0","same",1.105,1.125);
				
	mass = signalFcn->GetParameter(1);
	sigma = signalFcn->GetParameter(2);
				
	Double_t left = mass - nsig * sigma;
	Double_t right = mass + nsig * sigma;
	left = TMath::Max(left,1.105);
	right = TMath::Min(right,1.13);
				
	bin_left = hm0_signal->FindBin(left);
	bin_right = hm0_signal->FindBin(right);
				
	entries_new = hm0_signal->Integral(bin_left,bin_right);
	entries_old = hm0->GetEntries();			

	efficiency = 100.*entries_new/entries_old;	
	
	for (Int_t ib = bin_left; ib <= bin_right; ++ib) 
	{
		sum_sig += hm0_signal->GetBinContent(ib); 
		sum_full += hm0->GetBinContent(ib);
	}
	ratio_SB = sum_sig/(sum_full - sum_sig);
	ratio_SSB = sum_sig/TMath::Sqrt(sum_full);
	
	int muh1 = 0.6*hm0_signal->GetMaximum();
	TLine *line_left = new TLine(left,0,left,muh1);
	line_left->SetLineColor(kBlack);
	line_left->SetLineWidth(2);
	line_left->SetLineStyle(2);
	line_left->Draw("same");
	
	TLine *line_right = new TLine(right,0,right,muh1);
	line_right->SetLineColor(kBlack);
	line_right->SetLineWidth(2);
	line_right->SetLineStyle(2);
	line_right->Draw("same");
	
	TLine *line_backleft = new TLine(xmin_cut,0,xmin_cut,0.6*hm0->GetMaximum());
	line_backleft->SetLineColor(kBlue+1);
	line_backleft->SetLineWidth(2);
	line_backleft->SetLineStyle(9);
	line_backleft->Draw("same");
	
	TLine *line_backright = new TLine(xmax_cut,0,xmax_cut,0.6*hm0->GetMaximum());
	line_backright->SetLineColor(kBlue+1);
	line_backright->SetLineWidth(2);
	line_backright->SetLineStyle(9);
	line_backright->Draw("same");
	
	TLegend *legend=new TLegend(0.6,0.65,0.88,0.9);
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->SetBorderSize(0);
	legend->AddEntry(hm0,"Data","lp");
	legend->AddEntry(backFcn,"Bckg fit","l");	
	legend->AddEntry(line_left,"Cut-off (signal)","l");
	legend->AddEntry(line_backleft,"Cut-off (bckg)","l");
	legend->Draw("same");
		
	TLatex *title_SSB = new TLatex(1.075,0.86*hm0->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB)); 
	title_SSB->Draw("same");
	TLatex *title_efficiency = new TLatex(1.075,0.72*hm0->GetMaximum(), Form("eff. = %.2f",efficiency)); 
	title_efficiency->Draw("same");
	
	cout << "File = " << fullpath << "; Highest SSB ratio = " << ratio_SSB << endl;
		
	
}
