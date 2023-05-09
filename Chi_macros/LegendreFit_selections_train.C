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

void LegendreFit_selections_v3(TString inFileDir = "Directory_Path", TString outname = "Selections_Full.root", int NITER_CENT = 4, int NITER = 10, double xmin = 1.086, double nsig_bckg = 7.0, double nsig = 3.0)
{
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	
	const char *cent_values[4] = {"0 - 10 %", "10 - 20 %", "20 - 50 %", "50 - 100 %"};
	
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
	
	double max_value_final[NITER_CENT];
	int max_value_iter_pi_final[NITER_CENT];
	int max_value_iter_p_final[NITER_CENT];
	int max_value_iter_V0_final[NITER_CENT];
	int max_value_iter_path_final[NITER_CENT];
	int max_value_iter_angle_final[NITER_CENT];
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		max_value_final[iter_cent] = 0.;
		max_value_iter_pi_final[iter_cent] = 0;
		max_value_iter_p_final[iter_cent] = 0;
		max_value_iter_V0_final[iter_cent] = 0;
		max_value_iter_path_final[iter_cent] = 0;
		max_value_iter_angle_final[iter_cent] = 0;
	}
	//exit(0);
	
	//finding the directory and the files in it, then running on each file the fitting routine to find the best value of significance:
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				TString outdir = Form("cent%d_chipi%d_chip%d",iter_cent,iter_pi,iter_p); 
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
				
				if(max_value_file > max_value_final[iter_cent])
				{
					max_value_final[iter_cent] = max_value_file;
					max_value_iter_pi_final[iter_cent] = iter_pi;
					max_value_iter_p_final[iter_cent] = iter_p;
					max_value_iter_V0_final[iter_cent] = max_value_iter_V0_file;
					max_value_iter_path_final[iter_cent] = max_value_iter_path_file;
					max_value_iter_angle_final[iter_cent] = max_value_iter_angle_file;
				}
				
				cout << "File = " << fullpath << "; Highest SSB ratio = " << max_value_file << "; max_value_iter_V0 = " << max_value_iter_V0_file << "; max_value_iter_path = " << max_value_iter_path_file << "; max_value_iter_angle = " << max_value_iter_angle_file <<  endl;			
			}
		}
	}
	
	ofstream output_file_selections;
	output_file_selections.open (Form("ChiSelection_values_Cent_%d.txt",NITER_CENT));
	
	//Final values (for each centrality interval):
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		cout << "iter_cent = " << iter_cent << "; Highest SSB ratio = " << max_value_final[iter_cent] << "; max_value_iter_pi = " << max_value_iter_pi_final[iter_cent] << "; max_value_iter_p = " << max_value_iter_p_final[iter_cent] << "; max_value_iter_V0 = " << max_value_iter_V0_final[iter_cent] << "; max_value_iter_path = " << max_value_iter_path_final[iter_cent] << "; max_value_iter_angle = " << max_value_iter_angle_final[iter_cent] <<  endl;
		
		output_file_selections << chi_pi_value[max_value_iter_pi_final[iter_cent]] << "  " << chi_p_value[max_value_iter_p_final[iter_cent]] << "  " << chi_V0_value[max_value_iter_V0_final[iter_cent]] << "  " << lambda_path_value[max_value_iter_path_final[iter_cent]] << "  " << lambda_angle_value[max_value_iter_angle_final[iter_cent]] << "\n"; //write to file
		
		if(max_value_iter_pi_final[iter_cent] == NITER-1 || max_value_iter_pi_final[iter_cent] == 0) 
		{
			cout << "chi_pi value at the end of range! " << "; iter = " << max_value_iter_pi_final[iter_cent] << "; value = " << chi_pi_value[max_value_iter_pi_final[iter_cent]] << endl;
		}
		if(max_value_iter_p_final[iter_cent] == NITER-1 || max_value_iter_p_final[iter_cent] == 0) 
		{
			cout << "chi_p value at the end of range! " << "; iter = " << max_value_iter_p_final[iter_cent] << "; value = " << chi_p_value[max_value_iter_p_final[iter_cent]] << endl;
		}
		if(max_value_iter_V0_final[iter_cent] == NITER-1 || max_value_iter_V0_final[iter_cent] == 0) 
		{
			cout << "chi_V0 value at the end of range! " << "; iter = " << max_value_iter_V0_final[iter_cent] << "; value = " << chi_V0_value[max_value_iter_V0_final[iter_cent]] << endl;
		}
		if(max_value_iter_path_final[iter_cent] == NITER-1 || max_value_iter_path_final[iter_cent] == 0) 
		{
			cout << "path value at the end of range! " << "; iter = " << max_value_iter_path_final[iter_cent] << "; value = " << lambda_path_value[max_value_iter_path_final[iter_cent]] << endl;
		}
		if(max_value_iter_angle_final[iter_cent] == NITER-1 || max_value_iter_angle_final[iter_cent] == 0) 
		{
			cout << "angle value at the end of range! " << "; iter = " << max_value_iter_angle_final[iter_cent] << "; value = " << lambda_angle_value[max_value_iter_angle_final[iter_cent]] << endl;
		}
	}
	output_file_selections.close();
	
	//now lets plot the final fits:
	TH1D *hm0[NITER_CENT], *hm0_signal[NITER_CENT], *hm0_bckg[NITER_CENT];
	TF1 *fitting_fnc[NITER_CENT], *backFcn[NITER_CENT], *signalFcn[NITER_CENT];
	char *int_fitting_fnc = new char[100];
	char *int_backFcn = new char[100];
	char *int_signalFcn = new char[100];

	double sum_sig[NITER_CENT];
	double sum_full[NITER_CENT];
	double entries_new[NITER_CENT];
	double entries_old[NITER_CENT];
	double efficiency[NITER_CENT];
	double ratio_SB[NITER_CENT];
	double ratio_SSB[NITER_CENT];

	double xmax[NITER_CENT];
	int bin_left[NITER_CENT];
	int bin_right[NITER_CENT];
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sum_sig[iter_cent] = 0.;
		sum_full[iter_cent] = 0.;
	}

	TCanvas *c1[NITER_CENT];
	char *int_canvas_c1 = new char[100];

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sprintf(int_canvas_c1,"canvas_c1_%d",iter_cent);
		c1[iter_cent] = new TCanvas(int_canvas_c1,int_canvas_c1,0,0,600,600);
	
		TString outdir = Form("cent%d_chipi%d_chip%d",iter_cent,max_value_iter_pi_final[iter_cent],max_value_iter_p_final[iter_cent]); 
		TString fullpath;
		fullpath = inFileDir;
		fullpath += '/';
		fullpath += outdir;
		fullpath += '/';
		fullpath += outname;
		
		TFile *myFile_data = new TFile(fullpath);
		hm0[iter_cent] = (TH1D*) myFile_data->Get(Form("hm0_V%d_path%d_angle%d",max_value_iter_V0_final[iter_cent],max_value_iter_path_final[iter_cent],max_value_iter_angle_final[iter_cent]));
		hm0_signal[iter_cent] = (TH1D*)hm0[iter_cent]->Clone("hm0sig_V%d_path%d_angle%d");
		hm0_bckg[iter_cent] = (TH1D*)hm0[iter_cent]->Clone("hm0bckg_V%d_path%d_angle%d");
		xmax[iter_cent] = hm0[iter_cent]->GetXaxis()->GetXmax();
	
		sprintf(int_fitting_fnc,"int_fitting_fnc_V%d_path%d_angle%d",max_value_iter_V0_final[iter_cent],max_value_iter_path_final[iter_cent],max_value_iter_angle_final[iter_cent]);
		fitting_fnc[iter_cent] = new TF1(int_fitting_fnc,fitting_function,xmin,xmax[iter_cent],8);
	
		fitting_fnc[iter_cent]->SetNpx(500);
		fitting_fnc[iter_cent]->SetLineWidth(4);
		fitting_fnc[iter_cent]->SetLineColor(2); //red color for fitting line
		fitting_fnc[iter_cent]->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
		Int_t npar = fitting_fnc[iter_cent]->GetNumberFreeParameters();
		for (Int_t ip = 0; ip < npar; ++ip) fitting_fnc[iter_cent]->SetParameter(ip,0);	

		fitting_fnc[iter_cent]->SetParameter(0,10000);
		fitting_fnc[iter_cent]->SetParameter(1,1.116);
		fitting_fnc[iter_cent]->SetParameter(2,0.002);
					
		//change to the canvas:
		c1[iter_cent]->cd();	
		hm0[iter_cent]->SetMinimum(1.0E-20);
		hm0[iter_cent]->GetXaxis()->SetRangeUser(1.07,1.17);
		hm0[iter_cent]->Fit(int_fitting_fnc,"wq","",1.086,1.18);
					
		Float_t mass = fitting_fnc[iter_cent]->GetParameter(npar-7);
		Float_t sigma = TMath::Abs(fitting_fnc[iter_cent]->GetParameter(npar-6));
				
		Double_t xmin_cut = mass - nsig_bckg * sigma;
		Double_t xmax_cut = mass + nsig_bckg * sigma;
				
		Int_t imin = hm0[iter_cent]->FindBin(xmin_cut);
		Int_t imax = hm0[iter_cent]->FindBin(xmax_cut);
					
		//improve the picture:
		sprintf(int_backFcn,"backFcn_V%d_path%d_angle%d",max_value_iter_V0_final[iter_cent],max_value_iter_path_final[iter_cent],max_value_iter_angle_final[iter_cent]);
		backFcn[iter_cent] = new TF1(int_backFcn,background,xmin,xmax[iter_cent],5);
		backFcn[iter_cent]->SetLineColor(kMagenta);
					
		Double_t errs[200] = {0};
		Int_t nbins = hm0[iter_cent]->GetNbinsX();
		for (Int_t ib = 1; ib <= nbins; ++ib) 
		{
			if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg[iter_cent]->GetBinContent(ib));
		}
		hm0_bckg[iter_cent]->SetError(errs);
		hm0_bckg[iter_cent]->SetLineColor(kMagenta);
		hm0_bckg[iter_cent]->Fit(int_backFcn,"q","same",1.086,1.17);
					
		sprintf(int_signalFcn,"signalFcn_cut1_V%d_path%d_angle%d",max_value_iter_V0_final[iter_cent],max_value_iter_path_final[iter_cent],max_value_iter_angle_final[iter_cent]);
		signalFcn[iter_cent] = new TF1(int_signalFcn,gaussian,xmin,xmax[iter_cent],3);
		signalFcn[iter_cent]->SetLineColor(kBlue);
		signalFcn[iter_cent]->SetNpx(500);
					
		hm0_signal[iter_cent]->Sumw2();	
		hm0_signal[iter_cent]->Add(backFcn[iter_cent], -1);
		hm0_signal[iter_cent]->SetLineColor(kBlack);
		hm0_signal[iter_cent]->SetLineWidth(2);
		hm0_signal[iter_cent]->SetMarkerStyle(2);
		hm0_signal[iter_cent]->SetMarkerSize(1);
		signalFcn[iter_cent]->SetParameters(hm0[iter_cent]->GetMaximum(),1.1157,0.003);
		hm0_signal[iter_cent]->Fit(int_signalFcn,"q","same",1.105,1.125);
					
		mass = signalFcn[iter_cent]->GetParameter(1);
		sigma = signalFcn[iter_cent]->GetParameter(2);
					
		Double_t left = mass - nsig * sigma;
		Double_t right = mass + nsig * sigma;
		left = TMath::Max(left,1.105);
		right = TMath::Min(right,1.13);
					
		bin_left[iter_cent] = hm0_signal[iter_cent]->FindBin(left);
		bin_right[iter_cent] = hm0_signal[iter_cent]->FindBin(right);
					
		entries_new[iter_cent] = hm0_signal[iter_cent]->Integral(bin_left[iter_cent],bin_right[iter_cent]);
		entries_old[iter_cent] = hm0[iter_cent]->GetEntries();			

		efficiency[iter_cent] = 100.*entries_new[iter_cent]/entries_old[iter_cent];	
		
		for (Int_t ib = bin_left[iter_cent]; ib <= bin_right[iter_cent]; ++ib) 
		{
			sum_sig[iter_cent] += hm0_signal[iter_cent]->GetBinContent(ib); 
			sum_full[iter_cent] += hm0[iter_cent]->GetBinContent(ib);
		}
		ratio_SB[iter_cent] = sum_sig[iter_cent]/(sum_full[iter_cent] - sum_sig[iter_cent]);
		ratio_SSB[iter_cent] = sum_sig[iter_cent]/TMath::Sqrt(sum_full[iter_cent]);
		
		int muh1 = hm0_signal[iter_cent]->GetMaximum();
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
		
		TLine *line_backleft = new TLine(xmin_cut,0,xmin_cut,0.6*hm0[iter_cent]->GetMaximum());
		line_backleft->SetLineColor(kBlue+1);
		line_backleft->SetLineWidth(2);
		line_backleft->SetLineStyle(9);
		line_backleft->Draw("same");
		
		TLine *line_backright = new TLine(xmax_cut,0,xmax_cut,0.6*hm0[iter_cent]->GetMaximum());
		line_backright->SetLineColor(kBlue+1);
		line_backright->SetLineWidth(2);
		line_backright->SetLineStyle(9);
		line_backright->Draw("same");
		
		TLegend *legend=new TLegend(0.55,0.6,0.85,0.85);
		legend->SetTextFont(72);
		legend->SetTextSize(0.04);
		legend->SetBorderSize(0);
		legend->AddEntry(hm0[iter_cent],"Data","lp");
		legend->AddEntry(fitting_fnc[iter_cent],"Global fit","l");
		legend->AddEntry(backFcn[iter_cent],"Bckg fit","l");
		legend->AddEntry(signalFcn[iter_cent],"Signal","l");		
		legend->AddEntry(line_left,"Cut-off (signal)","l");
		legend->AddEntry(line_backleft,"Cut-off (bckg)","l");
		legend->Draw("same");
		
		TLatex *title1 = new TLatex(1.08,0.96*hm0[iter_cent]->GetMaximum(), cent_values[iter_cent]); 
		title1->Draw("same");
			
		TLatex *title_SSB = new TLatex(1.075,0.86*hm0[iter_cent]->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB[iter_cent])); 
		title_SSB->Draw("same");
		TLatex *title_efficiency = new TLatex(1.075,0.72*hm0[iter_cent]->GetMaximum(), Form("eff. = %.2f",efficiency[iter_cent])); 
		title_efficiency->Draw("same");
		
		cout << "File = " << fullpath << "; Highest SSB ratio = " << ratio_SSB[iter_cent] << endl;
		
	}
	
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

	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);

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
				c1->cd();	
				c1->Clear();
				hm0[iter_V0][iter_path][iter_angle]->SetMinimum(1.0E-20);
				hm0[iter_V0][iter_path][iter_angle]->GetXaxis()->SetRangeUser(1.07,1.17);
				hm0[iter_V0][iter_path][iter_angle]->Fit(int_fitting_fnc,"wqo","",1.086,1.18);
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
				hm0_bckg[iter_V0][iter_path][iter_angle]->Fit(int_backFcn,"q","same",1.086,1.17);
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
				hm0_signal[iter_V0][iter_path][iter_angle]->Fit(int_signalFcn,"q","same",1.105,1.125);
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
	if (c1) 
	{ 
		c1->Close(); 
		delete c1;
		c1 = 0; 
	}
	myFile_data->Close();
}
