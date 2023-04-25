///Plotting distributions of global polarization from the MCTest for a single particle (Lambda/ALambda)
//Distributions of model polarization (-P_{y}) in different centrality bins (for primary and full particles)
//Fitting of the angular distributions (\Psi_{RP} - \phi or \Psi_{EP} - \phi) for primary and full particles
//Comparison of mean values of global polarization: MC (full/primary), MCTracks (full/primary), STAR result for 11.5GeV with/without resolution correction

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
#define pi TMath::Pi()
Double_t ResEventPlane(Double_t chi);
Double_t Chi(Double_t res);
//fit function for the angle difference (delta(phi))
Double_t fitting_fnc_2orders(Double_t *x, Double_t *par) 
{
   return par[0]*(1. + 2.*par[1]*TMath::Sin(x[0]) + 2.*par[2]*TMath::Cos(x[0]) + 2.*par[3]*TMath::Sin(2.*x[0]) + 2.*par[4]*TMath::Cos(2.*x[0]));
}
Double_t fitting_fnc_3orders(Double_t *x, Double_t *par) 
{
   return par[0]*(1. + 2.*par[1]*TMath::Sin(x[0]) + 2.*par[2]*TMath::Cos(x[0]) + 2.*par[3]*TMath::Sin(2.*x[0]) + 2.*par[4]*TMath::Cos(2.*x[0]) + 2.*par[5]*TMath::Sin(3.*x[0]) + 2.*par[6]*TMath::Cos(3.*x[0]));
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
void Plot_single_particle(TString inFile = "Output_MCTest_Lambda.root", TString fitting_choice = "2orders", TString angle_choice = "RP", TString particle_choice = "Lambda", const int NITER_CENT = 4)
{
	TH1D *Lpolar_y[NITER_CENT], *Lpolar_y_prim[NITER_CENT];
	TH1D *CosTheta_hist[NITER_CENT], *Pstar_hist[NITER_CENT], *PstarRP_hist[NITER_CENT], *CosTheta_hist_prim[NITER_CENT], *Pstar_hist_prim[NITER_CENT], *PstarRP_hist_prim[NITER_CENT];
	
	TF1 *fnc_PstarRP[NITER_CENT], *fnc_PstarRP_prim[NITER_CENT];
	char *int_fnc_PstarRP = new char[NITER_CENT];
	char *int_fnc_PstarRP_prim = new char[NITER_CENT];
	
	Double_t xmin_CosTheta[NITER_CENT], xmax_CosTheta[NITER_CENT];
	Double_t xmin_Pstar[NITER_CENT], xmax_Pstar[NITER_CENT];
	
	double polar_par_hist[NITER_CENT], polar_par_hist_err[NITER_CENT], polar_par_hist_prim[NITER_CENT], polar_par_hist_prim_err[NITER_CENT], polar_y_mean[NITER_CENT], polar_y_mean_err[NITER_CENT], polar_y_mean_prim[NITER_CENT], polar_y_mean_prim_err[NITER_CENT];
	
	double binmin_PstarRP_hist[NITER_CENT];
	
//canvases for costheta distributions
	TCanvas *c1 = new TCanvas("Fitting Full","Fitting Full",0,0,600,600);
	TCanvas *c2 = new TCanvas("Fitting Primary","Fitting Primary",0,0,600,600);
	TCanvas *c3 = new TCanvas("P_{y} distributions (full)","P_{y} distributions (full)",0,0,600,600);
	TCanvas *c4 = new TCanvas("P_{y} distributions (primary)","P_{y} distributions (primary)",0,0,600,600);
	TCanvas *c5 = new TCanvas("Mean polarization","Mean polarization",0,0,600,600);
	
	char **cent_interval;
	if(NITER_CENT == 4)
	{
		c1->Divide(NITER_CENT/2,2);
		c2->Divide(NITER_CENT/2,2);
		c3->Divide(NITER_CENT/2,2);
		c4->Divide(NITER_CENT/2,2);
		cent_interval = (char *[]){"0 - 10%","10 - 20%","20 - 50%","50 - 100%"};
	}else if(NITER_CENT == 7)
	{
		c1->Divide(4,2);
		c2->Divide(4,2);
		c3->Divide(4,2);
		c4->Divide(4,2);
		cent_interval = (char *[]){"0 - 10%","10 - 20%","20 - 30%","30 - 40%","40 - 50%","50 - 60%","60 - 70%"};
	}else if(NITER_CENT == 10)
	{
		c1->Divide(NITER_CENT/2,2);
		c2->Divide(NITER_CENT/2,2);
		c3->Divide(NITER_CENT/2,2);
		c4->Divide(NITER_CENT/2,2);
		cent_interval = (char *[]){"0 - 10%","10 - 20%","20 - 30%","30 - 40%","40 - 50%","50 - 60%","60 - 70%","70 - 80%","80 - 90%","90 - 100%"};
	}else
	{
		cout << "This centrality binning is not defined!" << endl;
		return 1;
	}
	
	gStyle->SetOptFit();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	gROOT->ForceStyle();
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
		
	
		
	//const char *cent_interval[] = {"0 - 10%","10 - 20%","20 - 50%","50 - 100%"};
		
//input file:
	TFile *myFile_data = new TFile(inFile);

	TH1F *NCentr = (TH1F*) myFile_data->Get("NCentr");
	TH1F *NPositive = (TH1F*) myFile_data->Get("Positive Polarization");
	TH1F *NNegative = (TH1F*) myFile_data->Get("Negative Polarization");
	TH1D *Resolution_EP1_true, *Resolution_EP1_exp;
	double NEv_cent[NITER_CENT], ResEP1_true[NITER_CENT],ResEP1_exp[NITER_CENT],SubEvRes1[NITER_CENT];
	if (angle_choice == "RP")
	{
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			ResEP1_exp[iter_cent] = 1.;
		}
	}else if (angle_choice == "EP")
	{
		Resolution_EP1_true = (TH1D*) myFile_data->Get("Resolution_EP1_true");
		Resolution_EP1_exp = (TH1D*) myFile_data->Get("Resolution_EP1_exp");
		
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			NEv_cent[iter_cent] = NCentr->GetBinContent(iter_cent+1);
			ResEP1_true[iter_cent] = Resolution_EP1_true->GetBinContent(iter_cent+1)/NEv_cent[iter_cent];
			SubEvRes1[iter_cent] = Resolution_EP1_exp->GetBinContent(iter_cent+1)/NEv_cent[iter_cent];
			SubEvRes1[iter_cent] = TMath::Sqrt(SubEvRes1[iter_cent]);
			ResEP1_exp[iter_cent] = ResEventPlane(TMath::Sqrt(2.)*Chi(SubEvRes1[iter_cent]));
			cout << "iter_cent = " << iter_cent << "; ResEP1_true = " << ResEP1_true[iter_cent] << "; ResEP1_exp = " << ResEP1_exp[iter_cent] << endl;
		}
	}else
	{
		cout << "This angle choice is not defined! Please choose either RP or EP" << endl;
		return 1;
	}
			
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		Lpolar_y[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_%d",iter_cent));
		Lpolar_y_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_%d",iter_cent));		
		
		CosTheta_hist[iter_cent] = (TH1D*) myFile_data->Get(Form("CosTheta_hist_%d",iter_cent));
		Pstar_hist[iter_cent] = (TH1D*) myFile_data->Get(Form("Pstar_hist_%d",iter_cent));
		PstarRP_hist[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_%d",iter_cent));	
		
		CosTheta_hist_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("CosTheta_hist_prim_%d",iter_cent));
		Pstar_hist_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("Pstar_hist_prim_%d",iter_cent));
		PstarRP_hist_prim[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_prim_%d",iter_cent));	
		binmin_PstarRP_hist[iter_cent] = PstarRP_hist[iter_cent]->GetXaxis()->GetXmin();
	}
	
	//try to set errors:
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Int_t i_cut1 = PstarRP_hist[iter_cent]->FindBin(binmin_PstarRP_hist[iter_cent]);
		for(int iter_bin = 0; iter_bin < NITER_CENT; iter_bin++)
		{
			PstarRP_hist[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(PstarRP_hist[iter_cent]->GetBinContent(i_cut1)));
			i_cut1++;
		}
	}	
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		xmin_CosTheta[iter_cent] = CosTheta_hist[iter_cent]->GetXaxis()->GetXmin();
		xmax_CosTheta[iter_cent] = CosTheta_hist[iter_cent]->GetXaxis()->GetXmax();
		xmin_Pstar[iter_cent] = PstarRP_hist[iter_cent]->GetXaxis()->GetXmin();
		xmax_Pstar[iter_cent] = PstarRP_hist[iter_cent]->GetXaxis()->GetXmax();
		
		//create fitting functions for each histogram (delta(phi)_{RP} distributions):		
		
		sprintf(int_fnc_PstarRP,"fitting_fnc_PstarRP_%d",iter_cent);
		sprintf(int_fnc_PstarRP_prim,"fitting_fnc_PstarRP_prim_%d",iter_cent);
		if (fitting_choice == "2orders")
		{
			fnc_PstarRP[iter_cent] = new TF1(int_fnc_PstarRP,fitting_fnc_2orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
			fnc_PstarRP_prim[iter_cent] = new TF1(int_fnc_PstarRP_prim,fitting_fnc_2orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
		}else if (fitting_choice == "3orders")
		{
			fnc_PstarRP[iter_cent] = new TF1(int_fnc_PstarRP,fitting_fnc_3orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
			fnc_PstarRP_prim[iter_cent] = new TF1(int_fnc_PstarRP_prim,fitting_fnc_3orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
		}else
		{
			cout << "This fitting choice is not defined! Please choose either 2orders or 3orders" << endl;
			return 1;
		}
		
		fnc_PstarRP[iter_cent]->SetLineWidth(4);
		fnc_PstarRP[iter_cent]->SetLineColor(2); //red color for fitting line
		fnc_PstarRP[iter_cent]->SetParameters(10000., 0.01, 0.0003, 0.0002, 0.002);	
		fnc_PstarRP_prim[iter_cent]->SetLineWidth(4);
		fnc_PstarRP_prim[iter_cent]->SetLineColor(2); //red color for fitting line
		fnc_PstarRP_prim[iter_cent]->SetParameters(10000., 0.01, 0.0003, 0.0002, 0.002);
		
		//plotting:
		c1->cd(iter_cent+1);	
		
		PstarRP_hist[iter_cent]->SetYTitle("dN/d(#Delta#phi_{p}^{*})");
		PstarRP_hist[iter_cent]->SetXTitle("#Delta#phi_{p}^{*}");
		PstarRP_hist[iter_cent]->Draw("p9");
		PstarRP_hist[iter_cent]->Fit(int_fnc_PstarRP);
		
		PstarRP_hist[iter_cent]->Fit(int_fnc_PstarRP,"L");	
		//PstarRP_hist[iter_cent]->Fit(int_fnc_PstarRP,"w");	
		//PstarRP_hist[iter_cent]->Fit(int_fnc_PstarRP,"w","",xmin_Pstar[iter_cent],xmax_Pstar[iter_cent]);
		
		polar_par_hist[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		polar_par_hist_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
		TPaveText *t1_1 = new TPaveText(.13,.75,.5,.94);
		t1_1->SetTextAlign(22);
		t1_1->SetTextColor(kRed+2);
		t1_1->SetTextFont(72);
		t1_1->SetTextSize(0.04);
		t1_1->Paint("NDC");
		TText *t1_1_1 = t1_1->AddText(cent_interval[iter_cent]);
		TText *t1_1_2 = t1_1->AddText(Form("<P_{#Lambda}> (MCTracks) = %.4f +/- %.4f",polar_par_hist[iter_cent],polar_par_hist_err[iter_cent]));
		TText *t1_1_3 = t1_1->AddText(Form("<P_{#Lambda}> (full MC) = %.4f",-100.*Lpolar_y[iter_cent]->GetMean()));
		t1_1->Draw("same");

		
		c2->cd(iter_cent+1);	
				
		PstarRP_hist_prim[iter_cent]->SetYTitle("dN/d(#Delta#phi_{p}^{*})");
		PstarRP_hist_prim[iter_cent]->SetXTitle("#Delta#phi_{p}^{*}");
		PstarRP_hist_prim[iter_cent]->Draw("p9");
		PstarRP_hist_prim[iter_cent]->Fit(int_fnc_PstarRP_prim);
		polar_par_hist_prim[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_prim[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		polar_par_hist_prim_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_prim[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
		
		polar_y_mean[iter_cent] = -100.*Lpolar_y[iter_cent]->GetMean();
		polar_y_mean_err[iter_cent] = -100.*Lpolar_y[iter_cent]->GetMeanError();
		polar_y_mean_prim[iter_cent] = -100.*Lpolar_y_prim[iter_cent]->GetMean();
		polar_y_mean_prim_err[iter_cent] = -100.*Lpolar_y_prim[iter_cent]->GetMeanError();
		
		TPaveText *t2_2 = new TPaveText(.13,.75,.5,.94);
		t2_2->SetTextAlign(22);
		t2_2->SetTextColor(kRed+2);
		t2_2->SetTextFont(72);
		t2_2->SetTextSize(0.04);
		t2_2->Paint("NDC");
		TText *t2_2_1 = t2_2->AddText(cent_interval[iter_cent]);
		TText *t2_2_2 = t2_2->AddText(Form("<P_{#Lambda}> (MCTracks) = %.4f +/- %.4f",polar_par_hist_prim[iter_cent],polar_par_hist_prim_err[iter_cent]));
		TText *t2_2_3 = t2_2->AddText(Form("<P_{#Lambda}> (prim MC) = %.4f",-100.*Lpolar_y_prim[iter_cent]->GetMean()));
		t2_2->Draw("same");
		
		c3->cd(iter_cent+1);	
		
		Lpolar_y[iter_cent]->Draw();
		t1_1->Draw("same");
		
		TPaveText *t3_1 = new TPaveText(.7,.75,.94,.94);
		t3_1->SetTextAlign(22);
		t3_1->SetTextColor(kRed+2);
		t3_1->SetTextFont(72);
		t3_1->SetTextSize(0.04);
		t3_1->Paint("NDC");
		TText *t3_1_1 = t3_1->AddText(Form("N_{#Lambda} = %1.2e",Lpolar_y[iter_cent]->GetEntries()));
		t3_1->Draw("same");
		
		c4->cd(iter_cent+1);	
		
		Lpolar_y_prim[iter_cent]->Draw();
		t2_2->Draw("same");		
		
		TPaveText *t4_1 = new TPaveText(.7,.75,.94,.94);
		t4_1->SetTextAlign(22);
		t4_1->SetTextColor(kRed+2);
		t4_1->SetTextFont(72);
		t4_1->SetTextSize(0.04);
		t4_1->Paint("NDC");
		TText *t4_1_1 = t4_1->AddText(Form("N_{#Lambda} = %1.2e",Lpolar_y_prim[iter_cent]->GetEntries()));
		t4_1->Draw("same");
	}
	
	int *centrality_min;
	int *centrality_max;
	double *noErr;
	double centrality_bin[NITER_CENT]; 
	if (NITER_CENT == 4)
	{
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		noErr = init_double_array(4, 0, 0., 0., 0., 0.);		
	}else if (NITER_CENT == 7)
	{
		centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
		noErr = init_double_array(7, 0, 0., 0., 0., 0., 0., 0., 0.);		
	}else if (NITER_CENT == 10)
	{
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		noErr = init_double_array(10, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
	}else
	{
		cout << "This centrality binning is not defined!" << endl;
		return 1;
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		centrality_bin[iter_cent] = centrality_min[iter_cent] + (centrality_max[iter_cent] - centrality_min[iter_cent])/2;
	}
	double star_cent_bin[] = {centrality_bin[2]+1.5};
	double star_cent_bin_err[] = {0.};

	double star_pol_par[1], star_pol_err[1], star_pol_par_nores[1], star_pol_err_nores[1];
//For Lambda:
	if (particle_choice == "Lambda")
	{
		star_pol_par[0] = 1.179; // corrected value for 11.5GeV
		star_pol_err[0] = 0.347; // corrected error
		
		star_pol_par_nores[0] = 0.693; // corrected value for 11.5GeV (nores correction)
		star_pol_err_nores[0] = 0.202; // corrected error (nores correction)
	}else if (particle_choice == "ALambda")
	{
		star_pol_par[0] = 1.580; // corrected value for 11.5GeV
		star_pol_err[0] = 1.106; // corrected error
		
		star_pol_par_nores[0] = 0.921; // corrected value for 11.5GeV (nores correction)
		star_pol_err_nores[0] = 0.649; // corrected error (nores correction)
	}else
	{
		cout << "This particle choice is not defined! Please provide the definition in the code." << endl;
		return 1;
	}
						
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
	
	TGraphErrors *Polar_mean_MC = new TGraphErrors(NITER_CENT, centrality_bin, polar_y_mean, noErr, polar_y_mean_err);
	Polar_mean_MC->SetName("Polar_mean_MC");
	Polar_mean_MC->SetTitle("Polar_mean_MC");
	Polar_mean_MC->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MC->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MC->SetLineColor(kRed+1);
	Polar_mean_MC->SetMarkerColor(kRed+1);
	Polar_mean_MC->SetMarkerSize(2);
	Polar_mean_MC->SetMarkerStyle(24);	

	TGraphErrors *Polar_mean_MC_prim = new TGraphErrors(NITER_CENT, centrality_bin, polar_y_mean_prim, noErr, polar_y_mean_prim_err);
	Polar_mean_MC_prim->SetName("Polar_mean_MC");
	Polar_mean_MC_prim->SetTitle("Polar_mean_MC");
	Polar_mean_MC_prim->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MC_prim->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MC_prim->SetLineColor(kRed+1);
	Polar_mean_MC_prim->SetMarkerColor(kRed+1);
	Polar_mean_MC_prim->SetMarkerSize(2);
	Polar_mean_MC_prim->SetMarkerStyle(26);
	
	TGraphErrors *Polar_mean_MCtracks = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_hist, noErr, polar_par_hist_err);
	Polar_mean_MCtracks->SetName("Polar_mean_MCtracks");
	Polar_mean_MCtracks->SetTitle("Polar_mean_MCtracks");
	Polar_mean_MCtracks->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MCtracks->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MCtracks->SetLineColor(kBlack);
	Polar_mean_MCtracks->SetMarkerColor(kBlack);
	Polar_mean_MCtracks->SetMarkerSize(2);
	Polar_mean_MCtracks->SetMarkerStyle(20);
	
	TGraphErrors *Polar_mean_MCtracks_prim = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_hist_prim, noErr, polar_par_hist_prim_err);
	Polar_mean_MCtracks_prim->SetName("Polar_mean_MCtracks");
	Polar_mean_MCtracks_prim->SetTitle("Polar_mean_MCtracks");
	Polar_mean_MCtracks_prim->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MCtracks_prim->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MCtracks_prim->SetLineColor(kBlack);
	Polar_mean_MCtracks_prim->SetMarkerColor(kBlack);
	Polar_mean_MCtracks_prim->SetMarkerSize(2);
	Polar_mean_MCtracks_prim->SetMarkerStyle(22);  //for primary
	
	TLine *line = new TLine(0.,0.,80.,0.);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw("same");	
	
	c5->cd();	
	Polar_mean_MC->GetYaxis()->SetRangeUser(-1.0,6.);
	Polar_mean_MC->Draw("ap");
	Polar_mean_MCtracks->Draw("psame");
	Polar_mean_MCtracks_prim->Draw("psame");
	Polar_mean_MC_prim->Draw("psame");
	Polar_hist_STAR->Draw("psame");
	Polar_hist_STAR_nores->Draw("psame");
	
	TLegend *legend1_1=new TLegend(0.15,0.65,0.35,0.9);
	legend1_1->SetTextFont(72);
	legend1_1->SetTextSize(0.04);
	legend1_1->SetBorderSize(0);
	legend1_1->AddEntry(Polar_mean_MC,"MC (full)","p");
	legend1_1->AddEntry(Polar_mean_MCtracks,"MCTracks (full)","p");	
	legend1_1->AddEntry(Polar_mean_MC_prim,"MC (primary)","p");
	legend1_1->AddEntry(Polar_mean_MCtracks_prim,"MCTracks (primary)","p");
	legend1_1->AddEntry(Polar_hist_STAR,"STAR (11.5GeV)","p");
	legend1_1->AddEntry(Polar_hist_STAR_nores,"STAR (w/o res.)","p");
	legend1_1->Draw("same");
	
	line->Draw("same");	
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		float Npol_pos = NPositive->GetBinContent(iter_cent+1);
		float Npol_neg = NNegative->GetBinContent(iter_cent+1);
		float Ncentr_number = NCentr->GetBinContent(iter_cent+1);
		
		cout << "cent = " << cent_interval[iter_cent] << "; Ncentr_number = " << Ncentr_number << "; number of positive L = " << Npol_pos << "; umber of negative L = " << Npol_neg << endl;
		cout << "Difference = " << (TMath::Abs(Npol_pos - Npol_neg))/2. << endl;
		float polar_number = 100.*TMath::Abs((Npol_pos - Npol_neg)/(Npol_pos + Npol_neg));
		//float polar_number = Npol_pos/(Npol_pos + Npol_neg);
		cout << "Polarization = " << polar_number << endl;
	}
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
