//creates distributions of ratios pol/non-pol over cos(theta)
//top plots - cos distributions (reco and true) for polarized/non-polar cases
//bottom plots - their ratios
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

//fit function for the angle difference (from the HADES presentation at SQM 2019)
Double_t fitting_fnc_2orders(Double_t *x, Double_t *par) 
{
   return par[0]*(1. + 2.*par[1]*TMath::Sin(x[0]) + 2.*par[2]*TMath::Cos(x[0]) + 2.*par[3]*TMath::Sin(2.*x[0]) + 2.*par[4]*TMath::Cos(2.*x[0]));
}
Double_t fitting_fnc_3orders(Double_t *x, Double_t *par) 
{
   return par[0]*(1. + 2.*par[1]*TMath::Sin(x[0]) + 2.*par[2]*TMath::Cos(x[0]) + 2.*par[3]*TMath::Sin(2.*x[0]) + 2.*par[4]*TMath::Cos(2.*x[0]) + 2.*par[5]*TMath::Sin(3.*x[0]) + 2.*par[6]*TMath::Cos(3.*x[0]));
}

Double_t fitting_fnc_wo_constants(Double_t *x, Double_t *par) 
{
   return par[0]*(1. + 2.*3.14*0.732*par[1]*TMath::Sin(x[0])/8. + 2.*par[2]*TMath::Cos(x[0]) + 2.*par[3]*TMath::Sin(2.*x[0]) + 2.*par[4]*TMath::Cos(2.*x[0]));
}

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
void Global_Polarization_Fit_Delta_Phi(TString inFile = "invmass_output.root", const int NITER_CENT = 4)
{
	TH1D *angle_diff_inv_mass[NITER_CENT], *angle_diff_inv_mass_clone[NITER_CENT], *PstarRP_hist[NITER_CENT], *PstarRP_hist_MC[NITER_CENT], *Lpolar[NITER_CENT], *Lpolar_prim[NITER_CENT];
	char *int_angle_diff_inv_mass_clone = new char[100];

	TH1D *h_ratio_inv;

	TGraphErrors *Polar_hist_MC, *Polar_prim_hist_MC, *Polar_hist_RECO, *Polar_hist_RECO_wo_constants;

	double polar_par_MC[NITER_CENT], polar_par_MC_err[NITER_CENT], polar_prim_par_MC[NITER_CENT], polar_prim_par_MC_err[NITER_CENT], polar_par_RECO[NITER_CENT], polar_par_RECO_err[NITER_CENT], polar_par_MC_RECO[NITER_CENT], polar_par_MC_RECO_err[NITER_CENT], polar_par_RECO_wo_constants[NITER_CENT], polar_par_RECO_wo_constants_err[NITER_CENT];
	double centrality_bin[NITER_CENT]; 

	double resolution_reco[NITER_CENT], resolution_true[NITER_CENT], xmin_angle_diff_reco[NITER_CENT], xmax_angle_diff_reco[NITER_CENT], xmin_angle_diff_MC[NITER_CENT], xmax_angle_diff_MC[NITER_CENT];

	TF1 *fitting_fnc_angle_diff_reco[NITER_CENT], *fitting_fnc_angle_diff_MC[NITER_CENT], *fitting_fnc_angle_diff_reco_2[NITER_CENT];
	char *int_fitting_fnc_angle_diff_reco = new char[100];
	char *int_fitting_fnc_angle_diff_MC = new char[100];
	char *int_fitting_fnc_angle_diff_reco_2 = new char[100];

	double mean_pol_reco[NITER_CENT], mean_pol_MC[NITER_CENT], mean_pol_reco_err[NITER_CENT], mean_pol_MC_err[NITER_CENT];
	
	double par0_fix[NITER_CENT], par2_fix[NITER_CENT], par3_fix[NITER_CENT], par4_fix[NITER_CENT], par0_fix_2[NITER_CENT], par2_fix_2[NITER_CENT], par3_fix_2[NITER_CENT], par4_fix_2[NITER_CENT];

	TCanvas *c1[NITER_CENT];
	char *int_canvas_c1 = new char[100];
	
	TCanvas *c2 = new TCanvas("c2","c2",0,0,600,600);
	c2->Divide(2,2);
	//c2->Divide(4,2);
	TCanvas *c3 = new TCanvas("c3","c3",0,0,600,600);
	
	TCanvas *c4 = new TCanvas("c4","c4",0,0,600,600);
	if(NITER_CENT == 4)
	{
		c4->Divide(2,2);
	}else if(NITER_CENT == 7)
	{
		c4->Divide(4,2);
	}else if(NITER_CENT == 10)
	{
		c4->Divide(5,2);
	}else {cout << "weird stuff" << endl; return(0);}
	
	gStyle->SetOptFit();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	gROOT->ForceStyle();
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	int *centrality_min;
	int *centrality_max;
	double *noErr;
	
	if (NITER_CENT == 4)
	{
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		noErr = init_double_array(4, 0, 0., 0., 0., 0.);		
	}
	if (NITER_CENT == 7)
	{
		centrality_min = init_int_array(8, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(8, 0, 10, 20, 30, 40, 50, 60, 70);
		noErr = init_double_array(8, 0, 0., 0., 0., 0., 0., 0., 0.);		
	}
	if (NITER_CENT == 8)
	{
		centrality_min = init_int_array(8, 0, 0, 10, 20, 30, 40, 50, 60, 70);
		centrality_max = init_int_array(8, 0, 10, 20, 30, 40, 50, 60, 70, 80);
		noErr = init_double_array(8, 0, 0., 0., 0., 0., 0., 0., 0., 0.);		
	}
	if (NITER_CENT == 10)
	{
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		noErr = init_double_array(10, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
	}
	
//input file:
	TFile *myFile_data = new TFile(inFile);

	TH1D *Resolution_EP1_exp_sub = (TH1D*) myFile_data->Get("Resolution_EP1_exp");
	TH1D *Resolution_EP1_true = (TH1D*) myFile_data->Get("Resolution_EP1_true");
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		resolution_reco[iter_cent] = Resolution_EP1_exp_sub->GetBinContent(iter_cent+1);
		//resolution_reco[iter_cent] = 1.;
		resolution_true[iter_cent] = Resolution_EP1_true->GetBinContent(iter_cent+1);
		cout << "iter_cent = " << iter_cent << "; res_reco = " << resolution_reco[iter_cent] << "; res_true = " << resolution_true[iter_cent] << endl;
		//resolution_reco[iter_cent] = resolution_true[iter_cent];
		
		par0_fix[iter_cent] = 0.;
		par2_fix[iter_cent] = 0.;
		par3_fix[iter_cent] = 0.;
		par4_fix[iter_cent] = 0.;
		
		par0_fix_2[iter_cent] = 0.;
		par2_fix_2[iter_cent] = 0.;
		par3_fix_2[iter_cent] = 0.;
		par4_fix_2[iter_cent] = 0.;
		
	}
	
	//resolution_true[2] = 0.828737;
	//resolution_true[3] = 0.6283145;
//	return(0); //debugging
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sprintf(int_canvas_c1,"canvas_c1_%d",iter_cent);
		c1[iter_cent] = new TCanvas(int_canvas_c1,int_canvas_c1,0,0,600,600);
		c1[iter_cent]->Divide(2,1);
		
		centrality_bin[iter_cent] = centrality_min[iter_cent] + (centrality_max[iter_cent] - centrality_min[iter_cent])/2;
		
		angle_diff_inv_mass[iter_cent] = (TH1D*) myFile_data->Get(Form("Angle_diff_inv_mass_%d",iter_cent));
		angle_diff_inv_mass[iter_cent]->SetYTitle("dN/d(#Delta#phi_{p}^{*})");
		angle_diff_inv_mass[iter_cent]->SetXTitle("#Psi^{1}_{EP} - #phi_{p}^{*}");
		angle_diff_inv_mass[iter_cent]->SetMarkerStyle(20);
		angle_diff_inv_mass[iter_cent]->SetLineColor(kBlack);
		angle_diff_inv_mass[iter_cent]->SetMarkerSize(2);
		angle_diff_inv_mass[iter_cent]->SetLineWidth(2);
		sprintf(int_angle_diff_inv_mass_clone,"Angle_diff_inv_mass_clone_%d",iter_cent);
		angle_diff_inv_mass_clone[iter_cent] = (TH1D*)angle_diff_inv_mass[iter_cent]->Clone(int_angle_diff_inv_mass_clone);
		
		PstarRP_hist[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_%d",iter_cent));
		PstarRP_hist_MC[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_MC_%d",iter_cent));
		//Lpolar[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_%d",iter_cent));
		//Lpolar_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_prim_%d",iter_cent));
		Lpolar[iter_cent] = (TH1D*) myFile_data->Get(Form("Lpolar_%d",iter_cent));
		Lpolar_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("Lpolar_prim_%d",iter_cent));
		
		polar_par_MC[iter_cent] = -100.*Lpolar[iter_cent]->GetMean();
		polar_par_MC_err[iter_cent] = 100.*Lpolar[iter_cent]->GetMeanError();
		polar_prim_par_MC[iter_cent] = -100.*Lpolar_prim[iter_cent]->GetMean();
		polar_prim_par_MC_err[iter_cent] = 100.*Lpolar_prim[iter_cent]->GetMeanError();
		
		xmin_angle_diff_reco[iter_cent] = angle_diff_inv_mass[iter_cent]->GetXaxis()->GetXmin();
		xmax_angle_diff_reco[iter_cent] = angle_diff_inv_mass[iter_cent]->GetXaxis()->GetXmax();
		xmin_angle_diff_MC[iter_cent] = PstarRP_hist[iter_cent]->GetXaxis()->GetXmin();
		xmax_angle_diff_MC[iter_cent] = PstarRP_hist[iter_cent]->GetXaxis()->GetXmax();
		
		sprintf(int_fitting_fnc_angle_diff_reco,"fitting_fnc_angle_diff_reco_%d",iter_cent);
		fitting_fnc_angle_diff_reco[iter_cent] = new TF1(int_fitting_fnc_angle_diff_reco,fitting_fnc_2orders,xmin_angle_diff_reco[iter_cent],xmax_angle_diff_reco[iter_cent],5);
		fitting_fnc_angle_diff_reco[iter_cent]->SetLineWidth(4);
		fitting_fnc_angle_diff_reco[iter_cent]->SetLineColor(2); //red color for fitting line		
		//fitting_fnc_angle_diff_reco[iter_cent]->SetParLimits(1,0.,1.);
		
		
		sprintf(int_fitting_fnc_angle_diff_reco_2,"fitting_fnc_angle_diff_reco_2_%d",iter_cent);
		fitting_fnc_angle_diff_reco_2[iter_cent] = new TF1(int_fitting_fnc_angle_diff_reco_2,fitting_fnc_wo_constants,xmin_angle_diff_reco[iter_cent],xmax_angle_diff_reco[iter_cent],5);
		fitting_fnc_angle_diff_reco_2[iter_cent]->SetLineWidth(4);
		fitting_fnc_angle_diff_reco_2[iter_cent]->SetLineColor(2); //red color for fitting line		
	//	fitting_fnc_angle_diff_reco_2[iter_cent]->SetParLimits(1,0.,1.);
		
		sprintf(int_fitting_fnc_angle_diff_MC,"fitting_fnc_angle_diff_MC_%d",iter_cent);
		fitting_fnc_angle_diff_MC[iter_cent] = new TF1(int_fitting_fnc_angle_diff_MC,fitting_fnc_2orders,xmin_angle_diff_MC[iter_cent],xmax_angle_diff_MC[iter_cent],5);
		fitting_fnc_angle_diff_MC[iter_cent]->SetLineWidth(4);
		fitting_fnc_angle_diff_MC[iter_cent]->SetLineColor(2); //red color for fitting line
		fitting_fnc_angle_diff_MC[iter_cent]->SetParameters(10000., 0.01, 0.0003, 0.0002, 0.002);
		
		c1[iter_cent]->cd(1);	
		angle_diff_inv_mass[iter_cent]->SetYTitle("dN_{#Lambda}/d#Delta #phi");
		angle_diff_inv_mass[iter_cent]->Draw("p9");
		angle_diff_inv_mass[iter_cent]->Fit(int_fitting_fnc_angle_diff_reco);
		//angle_diff_inv_mass[iter_cent]->Fit(int_fitting_fnc_angle_diff_reco,"wl","",xmin_angle_diff_reco[iter_cent],xmax_angle_diff_reco[iter_cent]);
		
		par0_fix[iter_cent] = fitting_fnc_angle_diff_reco[iter_cent]->GetParameter(0);
		par2_fix[iter_cent] = fitting_fnc_angle_diff_reco[iter_cent]->GetParameter(2);
		par3_fix[iter_cent] = fitting_fnc_angle_diff_reco[iter_cent]->GetParameter(3);
		par4_fix[iter_cent] = fitting_fnc_angle_diff_reco[iter_cent]->GetParameter(4);
		
		fitting_fnc_angle_diff_reco[iter_cent]->SetParameter(0,par0_fix[iter_cent]);
		fitting_fnc_angle_diff_reco[iter_cent]->SetParameter(2,par2_fix[iter_cent]);
		fitting_fnc_angle_diff_reco[iter_cent]->SetParameter(3,par3_fix[iter_cent]);
		fitting_fnc_angle_diff_reco[iter_cent]->SetParameter(4,par4_fix[iter_cent]);
		
		//angle_diff_inv_mass[iter_cent]->Fit(int_fitting_fnc_angle_diff_reco,"wl","",xmin_angle_diff_reco[iter_cent],xmax_angle_diff_reco[iter_cent]);		
		
		cout << " ----------- Reco ----------- " << endl;
		cout << "True polarization (full) = " << -100.*Lpolar[iter_cent]->GetMean() << " +/- " << 100.*Lpolar[iter_cent]->GetMeanError() << endl;
		cout << "True polarization (primary) = " << -100.*Lpolar_prim[iter_cent]->GetMean() << " +/- " << 100.*Lpolar_prim[iter_cent]->GetMeanError() << endl;
		cout << "Parameter p1 = " << fitting_fnc_angle_diff_reco[iter_cent]->GetParameter(1) << " +/- " << fitting_fnc_angle_diff_reco[iter_cent]->GetParError(1) << endl;
		
		mean_pol_reco[iter_cent] = (8./(3.14*0.732))*(TMath::Abs(fitting_fnc_angle_diff_reco[iter_cent]->GetParameter(1))/resolution_reco[iter_cent]);
		mean_pol_reco_err[iter_cent] = (8./(3.14*0.732))*(TMath::Abs(fitting_fnc_angle_diff_reco[iter_cent]->GetParError(1))/resolution_reco[iter_cent]);
		
		cout << "Reco pol = " << 100.*mean_pol_reco[iter_cent] << " +/- " << 100.*mean_pol_reco_err[iter_cent] << endl;
		cout << " ----------- Reco ----------- " << endl;
		
		polar_par_RECO[iter_cent] = 100.*mean_pol_reco[iter_cent];
		polar_par_RECO_err[iter_cent] = 100.*mean_pol_reco_err[iter_cent];		
		
		TPaveText *t1 = new TPaveText(.2,.2,.65,.3);
		t1->SetTextAlign(22);
		t1->SetTextColor(kRed+2);
		t1->SetTextFont(72);
		t1->SetTextSize(0.04);
		t1->Paint("NDC");
		TText *t1_1 = t1->AddText(Form("p_{1} = %.4f +/- %.4f",fitting_fnc_angle_diff_reco[iter_cent]->GetParameter(1),fitting_fnc_angle_diff_reco[iter_cent]->GetParError(1)));
		TText *t1_2 = t1->AddText(Form("P_{H} = %.4f +/- %.4f",polar_par_RECO[iter_cent],polar_par_RECO_err[iter_cent]));
		t1->Draw("same");
		
		c1[iter_cent]->cd(2);	
		angle_diff_inv_mass_clone[iter_cent]->SetYTitle("dN_{#Lambda}/d#Delta #phi");
		angle_diff_inv_mass_clone[iter_cent]->Draw("p9");
		angle_diff_inv_mass_clone[iter_cent]->Fit(int_fitting_fnc_angle_diff_reco);
		//angle_diff_inv_mass_clone[iter_cent]->Fit(int_fitting_fnc_angle_diff_reco_2,"wl","",xmin_angle_diff_reco[iter_cent],xmax_angle_diff_reco[iter_cent]);
		
		par0_fix_2[iter_cent] = fitting_fnc_angle_diff_reco_2[iter_cent]->GetParameter(0);
		par2_fix_2[iter_cent] = fitting_fnc_angle_diff_reco_2[iter_cent]->GetParameter(2);
		par3_fix_2[iter_cent] = fitting_fnc_angle_diff_reco_2[iter_cent]->GetParameter(3);
		par4_fix_2[iter_cent] = fitting_fnc_angle_diff_reco_2[iter_cent]->GetParameter(4);
		
		fitting_fnc_angle_diff_reco_2[iter_cent]->SetParameter(0,par0_fix_2[iter_cent]);
		fitting_fnc_angle_diff_reco_2[iter_cent]->SetParameter(2,par2_fix_2[iter_cent]);
		fitting_fnc_angle_diff_reco_2[iter_cent]->SetParameter(3,par3_fix_2[iter_cent]);
		fitting_fnc_angle_diff_reco_2[iter_cent]->SetParameter(4,par4_fix_2[iter_cent]);
		
		//angle_diff_inv_mass_clone[iter_cent]->Fit(int_fitting_fnc_angle_diff_reco_2,"wl","",xmin_angle_diff_reco[iter_cent],xmax_angle_diff_reco[iter_cent]);
		
		polar_par_RECO_wo_constants[iter_cent] = 100.*(8./(3.14*0.732))*(TMath::Abs(fitting_fnc_angle_diff_reco_2[iter_cent]->GetParameter(1))/resolution_reco[iter_cent]);
		polar_par_RECO_wo_constants_err[iter_cent] = 100.*(8./(3.14*0.732))*(TMath::Abs(fitting_fnc_angle_diff_reco_2[iter_cent]->GetParError(1))/resolution_reco[iter_cent]);		
		
		TPaveText *t2 = new TPaveText(.2,.2,.65,.3);
		t2->SetTextAlign(22);
		t2->SetTextColor(kRed+2);
		t2->SetTextFont(72);
		t2->SetTextSize(0.04);
		t2->Paint("NDC");
		TText *t2_1 = t2->AddText(Form("p_{1} = %.4f +/- %.4f",fitting_fnc_angle_diff_reco_2[iter_cent]->GetParameter(1),fitting_fnc_angle_diff_reco_2[iter_cent]->GetParError(1)));
		TText *t2_2 = t2->AddText(Form("P_{H} = %.4f +/- %.4f",polar_par_RECO_wo_constants[iter_cent],polar_par_RECO_wo_constants_err[iter_cent]));
		t2->Draw("same");
		
		c4->cd(iter_cent+1);	
		
		angle_diff_inv_mass[iter_cent]->Draw();
		
		t1->Draw("same");
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		cout << "iter_cent = " << iter_cent << "; centrality_bin = " << centrality_bin[iter_cent] << "; polar_par_MC = " << polar_par_RECO[iter_cent] << "; polar_par_MC_err = " << polar_par_RECO_err[iter_cent] << endl;
	}
	
	Polar_hist_MC = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_MC, noErr, polar_par_MC_err);
	Polar_hist_MC->SetName("Polar_MC");
	Polar_hist_MC->SetTitle("Polar_MC");
	Polar_hist_MC->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_hist_MC->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_hist_MC->SetLineColor(kBlack);
	Polar_hist_MC->SetMarkerColor(kBlack);
	Polar_hist_MC->SetMarkerSize(2);
	Polar_hist_MC->SetMarkerStyle(24);	
	
	Polar_prim_hist_MC = new TGraphErrors(NITER_CENT, centrality_bin, polar_prim_par_MC, noErr, polar_prim_par_MC_err);
	Polar_prim_hist_MC->SetName("Polar_MC_prim");
	Polar_prim_hist_MC->SetTitle("Polar_MC_prim");
	Polar_prim_hist_MC->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_prim_hist_MC->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_prim_hist_MC->SetLineColor(kRed+1);
	Polar_prim_hist_MC->SetMarkerColor(kRed+1);
	Polar_prim_hist_MC->SetMarkerSize(2);
	Polar_prim_hist_MC->SetMarkerStyle(26);
	
	Polar_hist_RECO = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_RECO, noErr, polar_par_RECO_err);
	Polar_hist_RECO->SetName("Polar_RECO");
	Polar_hist_RECO->SetTitle("Polar_RECO");
	Polar_hist_RECO->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_hist_RECO->SetLineColor(kBlack);
	Polar_hist_RECO->SetMarkerColor(kBlack);
	Polar_hist_RECO->SetMarkerSize(2);
	Polar_hist_RECO->SetMarkerStyle(20);	
	
	Polar_hist_RECO_wo_constants = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_RECO_wo_constants, noErr, polar_par_RECO_wo_constants_err);
	Polar_hist_RECO_wo_constants->SetName("Polar_RECO");
	Polar_hist_RECO_wo_constants->SetTitle("Polar_RECO");
	Polar_hist_RECO_wo_constants->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_hist_RECO_wo_constants->SetLineColor(kBlack);
	Polar_hist_RECO_wo_constants->SetMarkerColor(kBlack);
	Polar_hist_RECO_wo_constants->SetMarkerSize(2);
	Polar_hist_RECO_wo_constants->SetMarkerStyle(20);
	
	double star_cent_bin[] = {centrality_bin[2]+1.5};
	double star_cent_bin_err[] = {0.};
//For Lambda:
	double star_pol_par[] = {1.179}; // corrected value for 11.5GeV
	double star_pol_err[] = {0.347}; // corrected error
	
	double star_pol_par_nores[] = {0.693}; // corrected value for 11.5GeV (nores correction)
	double star_pol_err_nores[] = {0.202}; // corrected error (nores correction)
	
//For ALambda:
/*	double star_pol_par[] = {1.580}; // corrected value for 11.5GeV
	double star_pol_err[] = {1.106}; // corrected error
	
	double star_pol_par_nores[] = {0.921}; // corrected value for 11.5GeV (nores correction)
	double star_pol_err_nores[] = {0.649}; // corrected error (nores correction)*/
	
	TGraphErrors *Polar_hist_STAR = new TGraphErrors(1, star_cent_bin, star_pol_par, star_cent_bin_err, star_pol_err);
	Polar_hist_STAR->SetName("Polar_STAR");
	Polar_hist_STAR->SetTitle("Polar_STAR");
	Polar_hist_STAR->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_hist_STAR->SetLineColor(kBlue+1);
	Polar_hist_STAR->SetMarkerColor(kBlue+1);
	Polar_hist_STAR->SetMarkerSize(2);
	Polar_hist_STAR->SetMarkerStyle(29);	
	
	TGraphErrors *Polar_hist_STAR_nores = new TGraphErrors(1, star_cent_bin, star_pol_par_nores, star_cent_bin_err, star_pol_err_nores);
	Polar_hist_STAR_nores->SetName("Polar_STAR");
	Polar_hist_STAR_nores->SetTitle("Polar_STAR");
	Polar_hist_STAR_nores->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_hist_STAR_nores->SetLineColor(kRed+1);
	Polar_hist_STAR_nores->SetMarkerColor(kRed+1);
	Polar_hist_STAR_nores->SetMarkerSize(2);
	Polar_hist_STAR_nores->SetMarkerStyle(29);	
	
	c2->cd();
	Polar_prim_hist_MC->GetYaxis()->SetRangeUser(-1.0,5.);
	Polar_prim_hist_MC->Draw("ap");
	Polar_hist_MC->Draw("psame");
	Polar_hist_RECO->Draw("psame");
	Polar_hist_STAR->Draw("psame");
	Polar_hist_STAR_nores->Draw("psame");
	
	TLegend *legend2_1=new TLegend(0.15,0.65,0.35,0.85);
	legend2_1->SetTextFont(72);
	legend2_1->SetTextSize(0.04);
	legend2_1->SetBorderSize(0);
	legend2_1->AddEntry(Polar_hist_MC,"MC (full)","p");
	legend2_1->AddEntry(Polar_prim_hist_MC,"MC (primary)","p");
	legend2_1->AddEntry(Polar_hist_RECO,"Reco","p");
//	legend2_1->AddEntry(Polar_hist_STAR,"STAR, Nature548 (2017) 62","p");	
	legend2_1->AddEntry(Polar_hist_STAR,"STAR (11.5GeV)","p");
	legend2_1->AddEntry(Polar_hist_STAR_nores,"STAR (w/o res.)","p");
	legend2_1->Draw("same");
	
	TLine *line = new TLine(0.,0.,80.,0.);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw("same");	
	
	c3->cd();
	Polar_prim_hist_MC->GetYaxis()->SetRangeUser(-1.0,5.);
	Polar_prim_hist_MC->Draw("ap");
	Polar_hist_MC->Draw("psame");
	Polar_hist_RECO_wo_constants->Draw("psame");
	
	TLegend *legend3_1=new TLegend(0.15,0.7,0.35,0.9);
	legend3_1->SetTextFont(72);
	legend3_1->SetTextSize(0.04);
	legend3_1->SetBorderSize(0);
	legend3_1->AddEntry(Polar_hist_MC,"MC (full)","p");
	legend3_1->AddEntry(Polar_prim_hist_MC,"MC (primary)","p");
	legend3_1->AddEntry(Polar_hist_RECO,"Reco (pol)","p");
//	legend3_1->AddEntry(Polar_hist_STAR,"STAR, Nature548 (2017) 62","p");	
	legend3_1->Draw("same");
	
	line->Draw("same");	

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		double n_L = 0.;
		for(int iter_bin = 0; iter_bin < 20; iter_bin++)
		{
			n_L += angle_diff_inv_mass[iter_cent]->GetBinContent(iter_bin+1);
		}
		cout << "N_L = " << n_L << endl;
	}
	
	TFile out("RecoPol_dots.root","recreate");
	Polar_prim_hist_MC->Write();
	Polar_hist_MC->GetYaxis()->SetRangeUser(-1.0,5.);
	Polar_hist_RECO->GetYaxis()->SetRangeUser(-1.0,5.);
	Polar_hist_MC->Write();
	Polar_hist_RECO->Write();
	
	out.Close();
}

