//creates distribution of invariant mass and performs the fitting procedure (Legendre polinoms + Gaus for Signal), then separate background fit in the sidebands,
// then cutting off the background, for the full phase space, with different selection criteria omega_2
//finds the optimal value of omega_2 for each centrality bin, corresponding to the maximum significance value
//writes the values to the txt file

#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <Riostream.h>
#include <TLatex.h>

#include <iostream>
using namespace std;
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

void LegendreFit_omega2(TString inFile = "PHSD_omegatest.root",double omega_value_start = 1.4, double omega_step = 0.1, int NITER_CENT = 4, int NITER = 30, double xmin = 1.086, double nsig_bckg = 7.0, double nsig = 4.0)
{
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	
	TH1D *hm0[NITER_CENT][NITER], *hm0_signal[NITER_CENT][NITER], *hm0_bckg[NITER_CENT][NITER];
	TF1 *fitting_fnc[NITER_CENT][NITER], *backFcn[NITER_CENT][NITER], *signalFcn[NITER_CENT][NITER];
	char *int_fitting_fnc = new char[100];
	char *int_backFcn = new char[100];
	char *int_signalFcn = new char[100];

	TH1D *hm0_full[NITER], *hm0_signal_full[NITER], *hm0_bckg_full[NITER];
	TF1 *fitting_fnc_full[NITER], *backFcn_full[NITER], *signalFcn_full[NITER];
	char *int_fitting_fnc_full = new char[100];
	char *int_backFcn_full = new char[100];
	char *int_signalFcn_full = new char[100];

	double sum_sig[NITER_CENT][NITER];
	double sum_full[NITER_CENT][NITER];
	double entries_new[NITER_CENT][NITER];
	double entries_old[NITER_CENT][NITER];
	double efficiency[NITER_CENT][NITER];
	double ratio_SB[NITER_CENT][NITER];
	double ratio_SSB[NITER_CENT][NITER];

	double sum_sig_full[NITER];
	double sum_full_full[NITER];
	double entries_new_full[NITER];
	double entries_old_full[NITER];
	double efficiency_full[NITER];
	double ratio_SB_full[NITER];
	double ratio_SSB_full[NITER];

	double xmax[NITER_CENT][NITER];
	int bin_left[NITER_CENT][NITER];
	int bin_right[NITER_CENT][NITER];

	double xmax_full[NITER];
	int bin_left_full[NITER];
	int bin_right_full[NITER];

	TCanvas *c1[NITER_CENT];
	char *int_canvas_c1 = new char[100];
	
	TCanvas *c2 = new TCanvas("c2","c2",0,0,600,600);
	if(NITER == 10) 
	{
		c2->Divide(5,2);
	}else
	if(NITER == 20)
	{
		c2->Divide(5,4);
	}else
	if(NITER == 30)
	{
		c2->Divide(5,6);
	}else {cout << "not set up value of omegas: " << NITER << endl;}

	double omega_value[NITER];
	
	
	for(int iter = 0; iter < NITER; iter++)
	{
		omega_value[iter] = omega_value_start + omega_step*iter;
	}
	const char *omega_values[30] = {"#omega_{2} > 1.4", "#omega_{2} > 1.5", "#omega_{2} > 1.6", "#omega_{2} > 1.7", "#omega_{2} > 1.8", "#omega_{2} > 1.9", "#omega_{2} > 2.0", "#omega_{2} > 2.1", "#omega_{2} > 2.2", "#omega_{2} > 2.3", "#omega_{2} > 2.4", "#omega_{2} > 2.5", "#omega_{2} > 2.6", "#omega_{2} > 2.7", "#omega_{2} > 2.8", "#omega_{2} > 2.9", "#omega_{2} > 3.0", "#omega_{2} > 3.1", "#omega_{2} > 3.2", "#omega_{2} > 3.3", "#omega_{2} > 3.4", "#omega_{2} > 3.5", "#omega_{2} > 3.6", "#omega_{2} > 3.7", "#omega_{2} > 3.8", "#omega_{2} > 3.9", "#omega_{2} > 4.0", "#omega_{2} > 4.1", "#omega_{2} > 4.2", "#omega_{2} > 4.3"}; 
	for (int i = 0; i < NITER; i++)
	{
		cout << "omega2_cutvalues[i] = " <<  omega_value[i] << "; omega_values[i] = " <<  omega_values[i] << endl;
	} 

	TFile *myFile_data = new TFile(inFile);

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sprintf(int_canvas_c1,"canvas_c1_%d",iter_cent);
		c1[iter_cent] = new TCanvas(int_canvas_c1,int_canvas_c1,0,0,600,600);
		if(NITER == 10) 
		{
			c1[iter_cent]->Divide(5,2);
		}else
		if(NITER == 20)
		{
			c1[iter_cent]->Divide(5,4);
		}else
		if(NITER == 30)
		{
			c1[iter_cent]->Divide(5,6);
		}else {cout << "not set up value of omegas: " << NITER << endl; continue;}
	
		for(int iter = 0; iter < NITER; iter++)
		{
			hm0[iter_cent][iter] = (TH1D*) myFile_data->Get(Form("hm0_%d_%d",iter_cent,iter));
			hm0_signal[iter_cent][iter] = (TH1D*)hm0[iter_cent][iter]->Clone("hm0sig_%d");
			hm0_bckg[iter_cent][iter] = (TH1D*)hm0[iter_cent][iter]->Clone("hm0bckg_%d");
			xmax[iter_cent][iter] = hm0[iter_cent][iter]->GetXaxis()->GetXmax();
			cout << "iter = " << iter << "; xmax = " << xmax[iter_cent][iter] << endl;
		}
	}

	for(int iter = 0; iter < NITER; iter++)
	{
		hm0_full[iter] = (TH1D*) myFile_data->Get(Form("hm0_full_%d",iter));
		hm0_signal_full[iter] = (TH1D*)hm0_full[iter]->Clone("hm0_fullsig_%d");
		hm0_bckg_full[iter] = (TH1D*)hm0_full[iter]->Clone("hm0_fullbckg_%d");
		xmax_full[iter] = hm0_full[iter]->GetXaxis()->GetXmax();
		cout << "iter = " << iter << "; xmax = " << xmax_full[iter] << endl;
	}

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter = 0; iter < NITER; iter++)
		{
			sprintf(int_fitting_fnc,"int_fitting_fnc_%d_%d",iter_cent,iter);
			fitting_fnc[iter_cent][iter] = new TF1(int_fitting_fnc,fitting_function,xmin,xmax[iter_cent][iter],8);
			
			fitting_fnc[iter_cent][iter]->SetNpx(500);
			fitting_fnc[iter_cent][iter]->SetLineWidth(4);
			fitting_fnc[iter_cent][iter]->SetLineColor(2); //red color for fitting line
			fitting_fnc[iter_cent][iter]->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
			Int_t npar = fitting_fnc[iter_cent][iter]->GetNumberFreeParameters();
			for (Int_t ip = 0; ip < npar; ++ip) fitting_fnc[iter_cent][iter]->SetParameter(ip,0);	

			fitting_fnc[iter_cent][iter]->SetParameter(0,10000);
			fitting_fnc[iter_cent][iter]->SetParameter(1,1.116);
			fitting_fnc[iter_cent][iter]->SetParameter(2,0.002);
			
			//change to the canvas:
			c1[iter_cent]->cd(iter+1);	
			hm0[iter_cent][iter]->SetMinimum(1.0E-20);
			hm0[iter_cent][iter]->GetXaxis()->SetRangeUser(1.07,1.17);
			hm0[iter_cent][iter]->Fit(int_fitting_fnc,"w","",1.086,1.18);
			
			Float_t mass = fitting_fnc[iter_cent][iter]->GetParameter(npar-7);
			Float_t sigma = TMath::Abs(fitting_fnc[iter_cent][iter]->GetParameter(npar-6));
		
			Double_t xmin_cut = mass - nsig_bckg * sigma;
			Double_t xmax_cut = mass + nsig_bckg * sigma;
		
			Int_t imin = hm0[iter_cent][iter]->FindBin(xmin_cut);
			Int_t imax = hm0[iter_cent][iter]->FindBin(xmax_cut);
			
			//improve the picture:
			sprintf(int_backFcn,"backFcn_%d_%d",iter_cent,iter);
			backFcn[iter_cent][iter] = new TF1(int_backFcn,background,xmin,xmax[iter_cent][iter],5);
			backFcn[iter_cent][iter]->SetLineColor(kMagenta);
			
			Double_t errs[200] = {0};
			Int_t nbins = hm0[iter_cent][iter]->GetNbinsX();
			for (Int_t ib = 1; ib <= nbins; ++ib) 
			{
				if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg[iter_cent][iter]->GetBinContent(ib));
			}
			hm0_bckg[iter_cent][iter]->SetError(errs);
			hm0_bckg[iter_cent][iter]->SetLineColor(kMagenta);
			hm0_bckg[iter_cent][iter]->Fit(int_backFcn,"","same",1.086,1.17);
			
			sprintf(int_signalFcn,"signalFcn_cut1_%d_%d",iter_cent,iter);
			signalFcn[iter_cent][iter] = new TF1(int_signalFcn,gaussian,xmin,xmax[iter_cent][iter],3);
			signalFcn[iter_cent][iter]->SetLineColor(kBlue);
			signalFcn[iter_cent][iter]->SetNpx(500);
			
			hm0_signal[iter_cent][iter]->Sumw2();	
			hm0_signal[iter_cent][iter]->Add(backFcn[iter_cent][iter], -1);
			hm0_signal[iter_cent][iter]->SetLineColor(kBlack);
			hm0_signal[iter_cent][iter]->SetLineWidth(2);
			hm0_signal[iter_cent][iter]->SetMarkerStyle(2);
			hm0_signal[iter_cent][iter]->SetMarkerSize(1);
			signalFcn[iter_cent][iter]->SetParameters(hm0[iter_cent][iter]->GetMaximum(),1.1157,0.003);
			hm0_signal[iter_cent][iter]->Fit(int_signalFcn,"","same",1.105,1.125);
			
			mass = signalFcn[iter_cent][iter]->GetParameter(1);
			sigma = signalFcn[iter_cent][iter]->GetParameter(2);
			
			Double_t left = mass - nsig * sigma;
			Double_t right = mass + nsig * sigma;
			left = TMath::Max(left,1.105);
			right = TMath::Min(right,1.13);
			
			bin_left[iter_cent][iter] = hm0_signal[iter_cent][iter]->FindBin(left);
			bin_right[iter_cent][iter] = hm0_signal[iter_cent][iter]->FindBin(right);
			
			entries_new[iter_cent][iter] = hm0_signal[iter_cent][iter]->Integral(bin_left[iter_cent][iter],bin_right[iter_cent][iter]);
			entries_old[iter_cent][iter] = hm0[iter_cent][iter]->GetEntries();
			
			int muh1 = hm0_signal[iter_cent][iter]->GetMaximum();
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
			
			TLegend *legend=new TLegend(0.15,0.65,0.35,0.85);
			legend->SetTextFont(72);
			legend->SetTextSize(0.04);
			legend->SetBorderSize(0);
			legend->AddEntry(hm0[iter_cent][iter],"Data","lp");
			legend->AddEntry(backFcn[iter_cent][iter],"Background","l");
			legend->AddEntry(signalFcn[iter_cent][iter],"Signal","l");
			legend->AddEntry(fitting_fnc[iter_cent][iter],"Global Fit","l");
			legend->AddEntry(line_left,"Cut-off range","l");
			legend->Draw("same");
			
			TLatex *title1 = new TLatex(1.13,0.96*hm0[iter_cent][iter]->GetMaximum(), omega_values[iter]); 
			title1->Draw("same");
		}
	}
	for(int iter = 0; iter < NITER; iter++)
	{
		sprintf(int_fitting_fnc_full,"int_fitting_fnc_full_%d",iter);
		fitting_fnc_full[iter] = new TF1(int_fitting_fnc_full,fitting_function,xmin,xmax_full[iter],8);
		
		fitting_fnc_full[iter]->SetNpx(500);
		fitting_fnc_full[iter]->SetLineWidth(4);
		fitting_fnc_full[iter]->SetLineColor(2); //red color for fitting line
		fitting_fnc_full[iter]->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
		Int_t npar = fitting_fnc_full[iter]->GetNumberFreeParameters();
		for (Int_t ip = 0; ip < npar; ++ip) fitting_fnc_full[iter]->SetParameter(ip,0);	

		fitting_fnc_full[iter]->SetParameter(0,10000);
		fitting_fnc_full[iter]->SetParameter(1,1.116);
		fitting_fnc_full[iter]->SetParameter(2,0.002);
		
		//change to the canvas:
		c2->cd(iter+1);	
		hm0_full[iter]->SetMinimum(1.0E-20);
		hm0_full[iter]->GetXaxis()->SetRangeUser(1.07,1.17);
		hm0_full[iter]->Fit(int_fitting_fnc_full,"w","",1.086,1.18);
		
		Float_t mass = fitting_fnc_full[iter]->GetParameter(npar-7);
		Float_t sigma = TMath::Abs(fitting_fnc_full[iter]->GetParameter(npar-6));
	
		Double_t xmin_cut = mass - nsig_bckg * sigma;
		Double_t xmax_cut = mass + nsig_bckg * sigma;
	
		Int_t imin = hm0_full[iter]->FindBin(xmin_cut);
		Int_t imax = hm0_full[iter]->FindBin(xmax_cut);
		
		//improve the picture:
		sprintf(int_backFcn_full,"backFcn_%d",iter);
		backFcn_full[iter] = new TF1(int_backFcn,background,xmin,xmax_full[iter],5);
		backFcn_full[iter]->SetLineColor(kMagenta);
		
		Double_t errs[200] = {0};
		Int_t nbins = hm0_full[iter]->GetNbinsX();
		for (Int_t ib = 1; ib <= nbins; ++ib) 
		{
			if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg_full[iter]->GetBinContent(ib));
		}
		hm0_bckg_full[iter]->SetError(errs);
		hm0_bckg_full[iter]->SetLineColor(kMagenta);
		hm0_bckg_full[iter]->Fit(int_backFcn_full,"","same",1.086,1.17);
		
		sprintf(int_signalFcn_full,"signalFcn_cut1_%d",iter);
		signalFcn_full[iter] = new TF1(int_signalFcn_full,gaussian,xmin,xmax_full[iter],3);
		signalFcn_full[iter]->SetLineColor(kBlue);
		signalFcn_full[iter]->SetNpx(500);
		
		hm0_signal_full[iter]->Sumw2();	
		hm0_signal_full[iter]->Add(backFcn_full[iter], -1);
		hm0_signal_full[iter]->SetLineColor(kBlack);
		hm0_signal_full[iter]->SetLineWidth(2);
		hm0_signal_full[iter]->SetMarkerStyle(2);
		hm0_signal_full[iter]->SetMarkerSize(1);
		signalFcn_full[iter]->SetParameters(hm0_full[iter]->GetMaximum(),1.1157,0.003);
		hm0_signal_full[iter]->Fit(int_signalFcn_full,"","same",1.105,1.125);
		
		mass = signalFcn_full[iter]->GetParameter(1);
		sigma = signalFcn_full[iter]->GetParameter(2);
		
		Double_t left = mass - nsig * sigma;
		Double_t right = mass + nsig * sigma;
		left = TMath::Max(left,1.105);
		right = TMath::Min(right,1.13);
		
		bin_left_full[iter] = hm0_signal_full[iter]->FindBin(left);
		bin_right_full[iter] = hm0_signal_full[iter]->FindBin(right);
		
		entries_new_full[iter] = hm0_signal_full[iter]->Integral(bin_left_full[iter],bin_right_full[iter]);
		entries_old_full[iter] = hm0_full[iter]->GetEntries();
		
		int muh1 = hm0_signal_full[iter]->GetMaximum();
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
		
		TLegend *legend=new TLegend(0.15,0.65,0.35,0.85);
		legend->SetTextFont(72);
		legend->SetTextSize(0.04);
		legend->SetBorderSize(0);
		legend->AddEntry(hm0_full[iter],"Data","lp");
		legend->AddEntry(backFcn_full[iter],"Background","l");
		legend->AddEntry(signalFcn_full[iter],"Signal","l");
		legend->AddEntry(fitting_fnc_full[iter],"Global Fit","l");
		legend->AddEntry(line_left,"Cut-off range","l");
		legend->Draw("same");
		
		TLatex *title1 = new TLatex(1.13,0.96*hm0_full[iter]->GetMaximum(), omega_values[iter]); 
		title1->Draw("same");
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter = 0; iter < NITER; iter++)
		{
			sum_sig[iter_cent][iter] = 0.;
			sum_full[iter_cent][iter] = 0.;
		}
	}
	for(int iter = 0; iter < NITER; iter++)
	{
		sum_sig_full[iter] = 0.;
		sum_full_full[iter] = 0.;
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter = 0; iter < NITER; iter++)
		{
			efficiency[iter_cent][iter] = 100.*entries_new[iter_cent][iter]/entries_old[iter_cent][iter];
			
			for (Int_t ib = bin_left[iter_cent][iter]; ib <= bin_right[iter_cent][iter]; ++ib) {
				sum_sig[iter_cent][iter] += hm0_signal[iter_cent][iter]->GetBinContent(ib);
				sum_full[iter_cent][iter] += hm0[iter_cent][iter]->GetBinContent(ib);
			}
			cout << "Sum of bins; Signal: " << sum_sig[iter_cent][iter] << "; Full: " << sum_full[iter_cent][iter] << endl;
			ratio_SB[iter_cent][iter] = sum_sig[iter_cent][iter]/(sum_full[iter_cent][iter] - sum_sig[iter_cent][iter]);
			ratio_SSB[iter_cent][iter] = sum_sig[iter_cent][iter]/TMath::Sqrt(sum_full[iter_cent][iter]);
			cout << "efficiency = " << efficiency[iter_cent][iter] << "; S/B = " << ratio_SB[iter_cent][iter] << "; S/(#sqrt{S+B}) = " << ratio_SSB[iter_cent][iter] << endl;
		}
	}
	for(int iter = 0; iter < NITER; iter++)
	{
		efficiency_full[iter] = 100.*entries_new_full[iter]/entries_old_full[iter];
		
		for (Int_t ib = bin_left_full[iter]; ib <= bin_right_full[iter]; ++ib) {
			sum_sig_full[iter] += hm0_signal_full[iter]->GetBinContent(ib);
			sum_full_full[iter] += hm0_full[iter]->GetBinContent(ib);
		}
		cout << "Sum of bins; Signal: " << sum_sig_full[iter] << "; Full: " << sum_full_full[iter] << endl;
		ratio_SB_full[iter] = sum_sig_full[iter]/(sum_full_full[iter] - sum_sig_full[iter]);
		ratio_SSB_full[iter] = sum_sig_full[iter]/TMath::Sqrt(sum_full_full[iter]);
		cout << "efficiency = " << efficiency_full[iter] << "; S/B = " << ratio_SB_full[iter] << "; S/(#sqrt{S+B}) = " << ratio_SSB_full[iter] << endl;
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter = 0; iter < NITER; iter++)
		{
			c1[iter_cent]->cd(iter+1);	
			//TLatex *title_SB = new TLatex(1.13,0.86*hm0[iter]->GetMaximum(), Form("S/B = %f",ratio_SB[iter])); 
			//title_SB->Draw("same");
			TLatex *title_SSB = new TLatex(1.13,0.86*hm0[iter_cent][iter]->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB[iter_cent][iter])); 
			title_SSB->Draw("same");
			TLatex *title_efficiency = new TLatex(1.13,0.76*hm0[iter_cent][iter]->GetMaximum(), Form("eff. = %.2f",efficiency[iter_cent][iter])); 
			title_efficiency->Draw("same");
		
		}
		//c1[iter_cent]->SaveAs(Form("%d.pdf",iter_cent));
	}
	
	ofstream output_file_omega;
	output_file_omega.open (Form("Omega2_values_Cent_%d.txt",NITER_CENT));

//output the highest values for each centrality interval and the corresponding omega_2 value;
	double max_value_full = ratio_SSB_full[0];
	int max_value_iter_full = 0;
	for(int iter = 0; iter < NITER; iter++)
	{
		if(ratio_SSB_full[iter] > max_value_full)
		{
			max_value_full = ratio_SSB_full[iter];
			max_value_iter_full = iter;
		}
	}
	cout << "; Highest SSB ratio = " << max_value_full << "; iter = " << max_value_iter_full << "; Omega_2 = " << omega_value[max_value_iter_full] << " --- " << omega_values[max_value_iter_full] << endl;
	
	output_file_omega << omega_value[max_value_iter_full] << "\n"; //write to file
	
	double max_value[NITER_CENT];
	int max_value_iter[NITER_CENT];
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		max_value[iter_cent] = ratio_SSB[iter_cent][0];
		max_value_iter[iter_cent] = 0;
		for(int iter = 0; iter < NITER; iter++)
		{
			if(ratio_SSB[iter_cent][iter] > max_value[iter_cent])
			{
				max_value[iter_cent] = ratio_SSB[iter_cent][iter];
				max_value_iter[iter_cent] = iter;
			}
		}
		cout << "centrality interval = " << iter_cent << "; Highest SSB ratio = " << max_value[iter_cent] << "; iter = " << max_value_iter[iter_cent] << "; Omega_2 = " << omega_value[max_value_iter[iter_cent]] << " --- " << omega_values[max_value_iter[iter_cent]] << endl;
		
		output_file_omega << omega_value[max_value_iter[iter_cent]] << "\n"; //write to file
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		if(max_value_iter[iter_cent] == NITER-1 ) 
		{
			cout << "Omega_2 value at the last bin!!! " << "; iter = " << max_value_iter[iter_cent] << endl;
		}
		if(max_value_iter[iter_cent] == 0 ) 
		{
			cout << "Omega_2 value at the first bin!!! " << "; iter_cent = " << iter_cent << "; iter = " << max_value_iter[iter_cent] << endl;
		}
	}
	output_file_omega.close();
	
}
