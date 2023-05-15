//creates distribution of invariant mass in different ranges of cos(theta), for 6 bins of theta
//optimized version 
//fitting procedure (Legendre polinoms + Gaus for Signal), then separate background fit in the sidebands, then cutting off the background
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

Double_t ResEventPlane(Double_t chi);
Double_t Chi(Double_t res);
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

//fit function for background (Legendre polinoms (L0 - L4))
Double_t background(Double_t *x, Double_t *par) {
	Double_t vxmix = 1.08;
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
	//Double_t vxmix = 1.07;
	Double_t vxmix = 1.08;
	Double_t vxmax = 1.17;
	Double_t e = 2.*(x[0] - vxmax)/(vxmax - vxmix) - 1.;
	return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*e + par[5]*0.5*(3.*TMath::Power(e,2) - 1.) + par[6]*0.5*e*(5.*TMath::Power(e,2) - 3.) + par[7]*0.125*(35.*TMath::Power(e,4) - 30.*TMath::Power(e,2) + 3.));
}

// Function for the fit the angular distribution
Double_t fitting_function_angular(Double_t *x, Double_t *par) {
	
	return par[0]*(1. + 2.*par[1]*TMath::Sin(x[0]) + 2.*par[2]*TMath::Cos(x[0]) + 2.*par[3]*TMath::Sin(2.*x[0]) + 2.*par[4]*TMath::Cos(2.*x[0]));
}
void Legendre_Fit_Delta_Phi_Inv_Mass(TString infile = "/scratch2/nazarova/tests/PHSD_output.root", TString outfile = "invmass_output.root", const int NITER_CENT = 4, const int NITER = 20, const double nsig = 4.0, Double_t xmin = 1.086, Double_t nsig_bckg = 7.0)
{
	TH1D *hm0[NITER_CENT][NITER], *hm0_signal[NITER_CENT][NITER], *hm0_bckg[NITER_CENT][NITER];
	TF1 *fitting_fnc[NITER_CENT][NITER], *backFcn[NITER_CENT][NITER], *signalFcn[NITER_CENT][NITER], *fitting_fnc_angular[NITER_CENT];
	char *int_fitting_fnc = new char[100];
	char *int_backFcn = new char[100];
	char *int_signalFcn = new char[100];
	char *int_fitting_fnc_angular = new char[100];

	TH1D *PstarRP_hist[NITER_CENT], *PstarRP_hist_MC[NITER_CENT], *Lpolar[NITER_CENT], *Lpolar_prim[NITER_CENT];
	TH1D *dNLambda, *dNLambda_MC, *dNLambda_Reco;

	TH1D *angle_diff_inv_mass[NITER_CENT];
	char *int_angle_diff_inv_mass = new char[100];

	double sum_sig[NITER_CENT][NITER];
	double sum_full[NITER_CENT][NITER];
	double entries_new[NITER_CENT][NITER];
	double entries_old[NITER_CENT][NITER];
	double efficiency[NITER_CENT][NITER];
	double ratio_SB[NITER_CENT][NITER];
	double ratio_SSB[NITER_CENT][NITER];

	double sin_mean[NITER_CENT];
	double sin_mean_hist[NITER_CENT];
	double pol_mean[NITER_CENT];
	double pol_mean_hist[NITER_CENT];
	double pol_mean_exp[NITER_CENT];

	double xmax[NITER_CENT][NITER];
	int bin_left[NITER_CENT][NITER];
	int bin_right[NITER_CENT][NITER];

	TCanvas *c1[NITER_CENT];
	char *int_canvas_c1 = new char[100];
	TCanvas *c2[NITER_CENT];
	char *int_canvas_c2 = new char[100];
	
	double par1_fix[NITER_CENT], par2_fix[NITER_CENT], par4_fix[NITER_CENT], par5_fix[NITER_CENT], par6_fix[NITER_CENT], par7_fix[NITER_CENT];
	double err1_fix[NITER_CENT], err2_fix[NITER_CENT], err4_fix[NITER_CENT], err5_fix[NITER_CENT], err6_fix[NITER_CENT], err7_fix[NITER_CENT];
	
	double ResEP1_exp[NITER_CENT], ResEP1_true[NITER_CENT], SubEvRes1[NITER_CENT];

	//gStyle->SetOptFit();
	gStyle->SetOptStat(0000);
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gStyle->SetOptTitle(0);
	
	double *_CentrBins;
	if (NITER_CENT == 4)
	{		
		_CentrBins = init_double_array(5, 0, 0.,10.,20.,50.,100.);
	}
	if (NITER_CENT == 8) //need to recheck - the last value is fishy
	{
		_CentrBins = init_double_array(9, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80.);	
	}
	if (NITER_CENT == 10)
	{
		_CentrBins = init_double_array(11, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.);		
	}
	
	Double_t nbx_PstarRP_hist[NITER_CENT];
	Double_t binmin_PstarRP_hist[NITER_CENT];
	Double_t binmax_PstarRP_hist[NITER_CENT];
	
	Double_t NEv_cent[NITER_CENT];
	
	const char *sin_interval[20] = {"bin1","bin2","bin3","bin4","bin5","bin6","bin7","bin8","bin9","bin10","bin11","bin12","bin13","bin14","bin15","bin16","bin17","bin18","bin19","bin20"};

	TFile *myFile_data = new TFile(infile);
	
	//dNLambda_MC = (TH1D*)myFile_data->Get("dNLambda_MC");
	//dNLambda = (TH1D*)myFile_data->Get("dNLambda");
	TH1F *NCentr = (TH1F*) myFile_data->Get("hNevCentr");
	TH1D *Resolution_EP1_exp_sub = (TH1D*) myFile_data->Get("hResolution_EP1_reco");
	TH1D *Resolution_EP1_true = (TH1D*) myFile_data->Get("hResolution_EP1_true");
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sprintf(int_canvas_c1,"canvas_c1_%d",iter_cent);
		c1[iter_cent] = new TCanvas(int_canvas_c1,int_canvas_c1,0,0,600,600);
		if (NITER == 6)
		{
			c1[iter_cent]->Divide(3,2);
		}
		if (NITER == 12)
		{
			c1[iter_cent]->Divide(6,2);
		}
		if (NITER == 20)
		{
			c1[iter_cent]->Divide(5,4);
		}
		if (NITER == 24)
		{
			c1[iter_cent]->Divide(6,4);
		}
		
		sprintf(int_canvas_c2,"canvas_c2_%d",iter_cent);
		c2[iter_cent] = new TCanvas(int_canvas_c2,int_canvas_c2,0,0,600,600);
		
		PstarRP_hist[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_%d",iter_cent));
		PstarRP_hist[iter_cent]->SetMarkerStyle(2);
		PstarRP_hist[iter_cent]->SetLineColor(kRed);
		PstarRP_hist[iter_cent]->SetMarkerSize(4);
		PstarRP_hist[iter_cent]->SetLineWidth(2);
		PstarRP_hist[iter_cent]->SetMinimum(1.0E-20);
		nbx_PstarRP_hist[iter_cent] = PstarRP_hist[iter_cent]->GetNbinsX();
		binmin_PstarRP_hist[iter_cent] = PstarRP_hist[iter_cent]->GetXaxis()->GetXmin();
		binmax_PstarRP_hist[iter_cent] = PstarRP_hist[iter_cent]->GetXaxis()->GetXmax();
		
		PstarRP_hist_MC[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_MC_%d",iter_cent));
		//Lpolar[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_%d",iter_cent));
		//Lpolar_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_prim_%d",iter_cent));
		Lpolar[iter_cent] = (TH1D*) myFile_data->Get(Form("Lpolar_%d",iter_cent));
		Lpolar_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("Lpolar_prim_%d",iter_cent));
		
		for(int iter = 0; iter < NITER; iter++)
		{
			hm0[iter_cent][iter] = (TH1D*) myFile_data->Get(Form("hm0_%d_%d",iter_cent,iter));
			hm0[iter_cent][iter]->SetYTitle("Entries");
			hm0[iter_cent][iter]->SetXTitle("M_{inv}, GeV/c^{2}");	
			hm0_signal[iter_cent][iter] = (TH1D*)hm0[iter_cent][iter]->Clone("Lambda_sig_%d");
			hm0_bckg[iter_cent][iter] = (TH1D*)hm0[iter_cent][iter]->Clone("Lambda_bckg_%d");
			xmax[iter_cent][iter] = hm0[iter_cent][iter]->GetXaxis()->GetXmax();
		}
		
		sprintf(int_angle_diff_inv_mass,"Angle_diff_inv_mass_%d",iter_cent);
		angle_diff_inv_mass[iter_cent] = new TH1D(int_angle_diff_inv_mass,int_angle_diff_inv_mass,nbx_PstarRP_hist[iter_cent],0.,6.28);
		angle_diff_inv_mass[iter_cent]->SetYTitle("Entries");
		angle_diff_inv_mass[iter_cent]->SetXTitle("#Psi^{1}_{EP} - #phi");
		angle_diff_inv_mass[iter_cent]->SetMarkerStyle(2);
		angle_diff_inv_mass[iter_cent]->SetLineColor(kBlack);
		angle_diff_inv_mass[iter_cent]->SetMarkerSize(4);
		angle_diff_inv_mass[iter_cent]->SetLineWidth(2);
		
		cout << iter_cent << ": centrality = " << NCentr->GetBinCenter(iter_cent+1) << endl;
		NEv_cent[iter_cent] = NCentr->GetBinContent(iter_cent+1);
		cout << "  number of events = " << NEv_cent[iter_cent] << endl;
		
		par1_fix[iter_cent] = 0.;
		par2_fix[iter_cent] = 0.;
		par4_fix[iter_cent] = 0.;
		par5_fix[iter_cent] = 0.;
		par6_fix[iter_cent] = 0.;
		par7_fix[iter_cent] = 0.;
		
		err1_fix[iter_cent] = 0.;
		err2_fix[iter_cent] = 0.;
		err4_fix[iter_cent] = 0.;
		err5_fix[iter_cent] = 0.;
		err6_fix[iter_cent] = 0.;
		err7_fix[iter_cent] = 0.;
	}
	
	dNLambda_Reco = new TH1D("dNLambda_Reco","dNLambda_Reco",NITER_CENT,_CentrBins);
	dNLambda_Reco->SetYTitle("N_{#Lambda}/N_{events}");
	dNLambda_Reco->SetXTitle("Centrality, [%]");
	dNLambda_Reco->SetMarkerStyle(20);
	dNLambda_Reco->SetLineColor(kBlue);
	dNLambda_Reco->SetMarkerColor(kBlue);
	dNLambda_Reco->SetMarkerSize(2);
	dNLambda_Reco->SetLineWidth(2);
	
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
			
			if(iter == 0) 
			{
				fitting_fnc[iter_cent][iter]->SetParameter(0,1000);
				fitting_fnc[iter_cent][iter]->SetParameter(1,1.116);
				fitting_fnc[iter_cent][iter]->SetParameter(2,0.002);
			}
			if(iter != 0) 
			{
				fitting_fnc[iter_cent][iter]->SetParameter(1,par1_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(2,par2_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(4,par4_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(5,par5_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(6,par6_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(7,par7_fix[iter_cent]);
				
			}
			//change to the canvas:
			c1[iter_cent]->cd(iter+1);	
			hm0[iter_cent][iter]->SetMinimum(1.0E-20);
			hm0[iter_cent][iter]->GetXaxis()->SetRangeUser(1.07,1.17);
			hm0[iter_cent][iter]->Fit(int_fitting_fnc,"w","",1.086,1.18);			
			
			if(iter == 0) 
			{
				par1_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParameter(1);
				par2_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParameter(2);
				par4_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParameter(4);
				par5_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParameter(5);
				par6_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParameter(6);
				par7_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParameter(7);
				
				//cout << "iter_cent = " << iter_cent << "; par1_fix = " << par1_fix[iter_cent] << endl;
				
				err1_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParError(1);
				err2_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParError(2);
				err4_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParError(4);
				err5_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParError(5);
				err6_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParError(6);
				err7_fix[iter_cent] = fitting_fnc[iter_cent][iter]->GetParError(7);
				
				fitting_fnc[iter_cent][iter]->SetParameter(1,par1_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(2,par2_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(4,par4_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(5,par5_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(6,par6_fix[iter_cent]);
				fitting_fnc[iter_cent][iter]->SetParameter(7,par7_fix[iter_cent]);
				
				hm0[iter_cent][iter]->Fit(int_fitting_fnc,"w","",1.086,1.18);	
			}
			
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
			
			TLatex *title1 = new TLatex(1.16,0.96*hm0[iter_cent][iter]->GetMaximum(), sin_interval[iter]); 
			title1->Draw("same");
		}
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter = 0; iter < NITER; iter++)
		{
			sum_sig[iter_cent][iter] = 0.;
			sum_full[iter_cent][iter] = 0.;
		}
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
	double sum_lambda_rec[NITER_CENT];
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sum_lambda_rec[iter_cent] = 0.;
		for(int iter = 0; iter < NITER; iter++)
		{
			c1[iter_cent]->cd(iter+1);	
			TLatex *title_SSB = new TLatex(1.13,0.86*hm0[iter_cent][iter]->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB[iter_cent][iter])); 
			title_SSB->Draw("same");
			TLatex *title_efficiency = new TLatex(1.13,0.76*hm0[iter_cent][iter]->GetMaximum(), Form("eff. = %.2f",efficiency[iter_cent][iter])); 
			title_efficiency->Draw("same");
		
		}
		Int_t i_cut1 = PstarRP_hist[iter_cent]->FindBin(binmin_PstarRP_hist[iter_cent]);
		
		for(Int_t ibin = 0; ibin < NITER; ibin++)
		{
			
			angle_diff_inv_mass[iter_cent]->SetBinContent(i_cut1,entries_new[iter_cent][ibin]);
			angle_diff_inv_mass[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(sum_full[iter_cent][ibin]));
			sum_lambda_rec[iter_cent] += entries_new[iter_cent][ibin];
			/*
			angle_diff_inv_mass[iter_cent]->SetBinContent(i_cut1,10.*entries_new[iter_cent][ibin]);
			angle_diff_inv_mass[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(10.*sum_full[iter_cent][ibin]));
			sum_lambda_rec[iter_cent] += 10.*entries_new[iter_cent][ibin];*/
			i_cut1++;
		}	
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		cout << "iter = " << iter_cent << "; Number of Lambda = " << sum_lambda_rec[iter_cent] << endl;
	}
	
	TFile out(outfile,"recreate");
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		//angle_diff_inv_mass[iter_cent]->Scale(1./NEv_cent[iter_cent],"width");
		angle_diff_inv_mass[iter_cent]->Write();
		PstarRP_hist[iter_cent]->Write();
		PstarRP_hist_MC[iter_cent]->Write();
		Lpolar[iter_cent]->Write();
		Lpolar_prim[iter_cent]->Write();
		
		cout << "par1 = " << par1_fix[iter_cent] << " +/- " << err1_fix[iter_cent] << endl;
		cout << "par2 = " << par2_fix[iter_cent] << " +/- " << err2_fix[iter_cent] << endl;
		cout << "par4 = " << par4_fix[iter_cent] << " +/- " << err4_fix[iter_cent] << endl;
		cout << "par5 = " << par5_fix[iter_cent] << " +/- " << err5_fix[iter_cent] << endl;
		cout << "par6 = " << par6_fix[iter_cent] << " +/- " << err6_fix[iter_cent] << endl;
		cout << "par7 = " << par1_fix[iter_cent] << " +/- " << err7_fix[iter_cent] << endl;
		
		//dNLambda_MC->SetBinContent(iter_cent+1,dNLambda_MC->GetBinContent(iter_cent+1)/NEv_cent[iter_cent]);
		//dNLambda->SetBinContent(iter_cent+1,dNLambda->GetBinContent(iter_cent+1)/NEv_cent[iter_cent]);
		ResEP1_true[iter_cent] = Resolution_EP1_true->GetBinContent(iter_cent+1)/NEv_cent[iter_cent];
		SubEvRes1[iter_cent] = Resolution_EP1_exp_sub->GetBinContent(iter_cent+1)/NEv_cent[iter_cent];
		
		//cout << "ResEP1_true[iter_cent] = " << ResEP1_true[iter_cent] << endl;
		//cout << "SubEvRes1[iter_cent] = " << SubEvRes1[iter_cent] << endl;
		//cout << "NEv_cent[iter_cent] = " << NEv_cent[iter_cent] << endl;
		
		SubEvRes1[iter_cent] = TMath::Sqrt(SubEvRes1[iter_cent]);
		ResEP1_exp[iter_cent] = ResEventPlane(TMath::Sqrt(2.)*Chi(SubEvRes1[iter_cent]));
		Resolution_EP1_exp_sub->SetBinContent(iter_cent+1,ResEP1_exp[iter_cent]);
		Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);
		
		cout << "ResEP1_exp[iter_cent] = " << ResEP1_exp[iter_cent] << endl;
	}
	
	//dNLambda->Write();
	//dNLambda_MC->Write();
	dNLambda_Reco->Write();
	Resolution_EP1_exp_sub->Write();
	Resolution_EP1_true->Write();
	
	out.Close();
}
Double_t ResEventPlane(Double_t chi)
{
  // plane resolution as function of chi
  Double_t con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
  Double_t arg = chi * chi / 4.;
  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res ;
}

Double_t Chi(Double_t res)
{
  // chi from the event plane resolution
 
  Double_t chi   = 2.0;
  Double_t delta = 1.0;
  for(int i = 0; i < 15; i++)// for(int i = 0; i < 50; i++) ---- check!!!
  {
   if(ResEventPlane(chi) < res) { chi = chi + delta ; }
   else                         { chi = chi - delta ; }
   delta = delta / 2.;
  }
	cout << "chi = " << chi << endl;
  return chi ;
}
