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
void Plot_Full(TString inFile = "Output_MCTest_Lambda.root", TString fitting_choice = "2orders", TString angle_choice = "RP", const int NITER_CENT = 4)
{
	TH1D *Lpolar_y_Lambda[NITER_CENT], *Lpolar_y_prim_Lambda[NITER_CENT], *Lpolar_y_ALambda[NITER_CENT], *Lpolar_y_prim_ALambda[NITER_CENT];
	TH1D *Lpolar_y_Xi[NITER_CENT], *Lpolar_y_prim_Xi[NITER_CENT], *Lpolar_y_AXi[NITER_CENT], *Lpolar_y_prim_AXi[NITER_CENT], *Lpolar_y_Xi0[NITER_CENT], *Lpolar_y_prim_Xi0[NITER_CENT], *Lpolar_y_AXi0[NITER_CENT], *Lpolar_y_prim_AXi0[NITER_CENT], *Lpolar_y_Sig0[NITER_CENT], *Lpolar_y_prim_Sig0[NITER_CENT], *Lpolar_y_ASig0[NITER_CENT], *Lpolar_y_prim_ASig0[NITER_CENT];
	
	TH1D *CosTheta_hist_Lambda[NITER_CENT], *Pstar_hist_Lambda[NITER_CENT], *PstarRP_hist_Lambda[NITER_CENT], *CosTheta_hist_prim_Lambda[NITER_CENT], *Pstar_hist_prim_Lambda[NITER_CENT], *PstarRP_hist_prim_Lambda[NITER_CENT];
	
	TH1D *CosTheta_hist_ALambda[NITER_CENT], *Pstar_hist_ALambda[NITER_CENT], *PstarRP_hist_ALambda[NITER_CENT], *CosTheta_hist_prim_ALambda[NITER_CENT], *Pstar_hist_prim_ALambda[NITER_CENT], *PstarRP_hist_prim_ALambda[NITER_CENT];
	
	TH1D *Lpolar_y_Lambda_from_Xi[NITER_CENT], *Lpolar_y_ALambda_from_AXi[NITER_CENT], *Lpolar_y_Lambda_from_Xi0[NITER_CENT], *Lpolar_y_ALambda_from_AXi0[NITER_CENT], *Lpolar_y_Lambda_from_Sig0[NITER_CENT], *Lpolar_y_ALambda_from_ASig0[NITER_CENT];
	
	TF1 *fnc_PstarRP_Lambda[NITER_CENT], *fnc_PstarRP_prim_Lambda[NITER_CENT];
	char *int_fnc_PstarRP_Lambda = new char[NITER_CENT];
	char *int_fnc_PstarRP_prim_Lambda = new char[NITER_CENT];
	
	TF1 *fnc_PstarRP_ALambda[NITER_CENT], *fnc_PstarRP_prim_ALambda[NITER_CENT];
	char *int_fnc_PstarRP_ALambda = new char[NITER_CENT];
	char *int_fnc_PstarRP_prim_ALambda = new char[NITER_CENT];
	
	Double_t xmin_CosTheta[NITER_CENT], xmax_CosTheta[NITER_CENT];
	Double_t xmin_Pstar[NITER_CENT], xmax_Pstar[NITER_CENT];
	
	double polar_par_hist_Lambda[NITER_CENT], polar_par_hist_err_Lambda[NITER_CENT], polar_par_hist_prim_Lambda[NITER_CENT], polar_par_hist_prim_err_Lambda[NITER_CENT], polar_y_mean_Lambda[NITER_CENT], polar_y_mean_err_Lambda[NITER_CENT], polar_y_mean_prim_Lambda[NITER_CENT], polar_y_mean_prim_err_Lambda[NITER_CENT];
	
	double polar_par_hist_ALambda[NITER_CENT], polar_par_hist_err_ALambda[NITER_CENT], polar_par_hist_prim_ALambda[NITER_CENT], polar_par_hist_prim_err_ALambda[NITER_CENT], polar_y_mean_ALambda[NITER_CENT], polar_y_mean_err_ALambda[NITER_CENT], polar_y_mean_prim_ALambda[NITER_CENT], polar_y_mean_prim_err_ALambda[NITER_CENT];
	
//canvases for costheta distributions
	TCanvas *c1 = new TCanvas("(Lambda) Fitting Full","(Lambda) Fitting Full",0,0,600,600);
	TCanvas *c2 = new TCanvas("(Lambda) Fitting Primary","(Lambda) Fitting Primary",0,0,600,600);
	TCanvas *c3 = new TCanvas("(Lambda) P_{y} (full)","(Lambda) P_{y} (full)",0,0,600,600);
	TCanvas *c4 = new TCanvas("(Lambda) P_{y} (primary)","(Lambda) P_{y} (primary)",0,0,600,600);
	TCanvas *c5 = new TCanvas("(Lambda) Mean polarization","(Lambda) Mean polarization",0,0,600,600);
	
	TCanvas *c6 = new TCanvas("(ALambda) Fitting Full","(ALambda) Fitting Full",0,0,600,600);
	TCanvas *c7 = new TCanvas("(ALambda) Fitting Primary","(ALambda) Fitting Primary",0,0,600,600);
	TCanvas *c8 = new TCanvas("(ALambda) P_{y} (full)","(ALambda) P_{y} (full)",0,0,600,600);
	TCanvas *c9 = new TCanvas("(ALambda) P_{y} (primary)","(ALambda) P_{y} (primary)",0,0,600,600);
	TCanvas *c10 = new TCanvas("(ALambda) Mean polarization","(ALambda) Mean polarization",0,0,600,600);
	
	TCanvas *c11 = new TCanvas("(Xi) P_{y} (full)","(Xi) P_{y} (full)",0,0,600,600);
	TCanvas *c12 = new TCanvas("(Xi) P_{y} (primary)","(Xi) P_{y} (primary)",0,0,600,600);
	TCanvas *c13 = new TCanvas("(Xi0) P_{y} (full)","(Xi0) P_{y} (full)",0,0,600,600);
	TCanvas *c14 = new TCanvas("(Xi0) P_{y} (primary)","(Xi0) P_{y} (primary)",0,0,600,600);
	TCanvas *c15 = new TCanvas("(Sig0) P_{y} (full)","(Sig0) P_{y} (full)",0,0,600,600);
	TCanvas *c16 = new TCanvas("(Sig0) P_{y} (primary)","(Sig0) P_{y} (primary)",0,0,600,600);
	
	TCanvas *c17 = new TCanvas("(AXi) P_{y} (full)","(AXi) P_{y} (full)",0,0,600,600);
	TCanvas *c18 = new TCanvas("(AXi) P_{y} (primary)","(AXi) P_{y} (primary)",0,0,600,600);
	TCanvas *c19 = new TCanvas("(AXi0) P_{y} (full)","(AXi0) P_{y} (full)",0,0,600,600);
	TCanvas *c20 = new TCanvas("(AXi0) P_{y} (primary)","(AXi0) P_{y} (primary)",0,0,600,600);
	TCanvas *c21 = new TCanvas("(ASig0) P_{y} (full)","(ASig0) P_{y} (full)",0,0,600,600);
	TCanvas *c22 = new TCanvas("(ASig0) P_{y} (primary)","(ASig0) P_{y} (primary)",0,0,600,600);
	
	TCanvas *c23 = new TCanvas("(Lambda and ALambda) Mean polarization","(Lambda and ALambda) Mean polarization)",0,0,600,600);
	
	char *cent_4bins[] = {"0 - 10%","10 - 20%","20 - 50%","50 - 100%"};
	char *cent_7bins[] = {"0 - 10%","10 - 20%","20 - 30%","30 - 40%","40 - 50%","50 - 60%","60 - 70%"};
	char *cent_10bins[] = {"0 - 10%","10 - 20%","20 - 30%","30 - 40%","40 - 50%","50 - 60%","60 - 70%","70 - 80%","80 - 90%","90 - 100%"};
	
	char *cent_interval[NITER_CENT];
	if(NITER_CENT == 4)
	{
		c1->Divide(NITER_CENT/2,2);
		c2->Divide(NITER_CENT/2,2);
		c3->Divide(NITER_CENT/2,2);
		c4->Divide(NITER_CENT/2,2);
		c6->Divide(NITER_CENT/2,2);
		c7->Divide(NITER_CENT/2,2);
		c8->Divide(NITER_CENT/2,2);
		c9->Divide(NITER_CENT/2,2);
		c11->Divide(NITER_CENT/2,2);
		c12->Divide(NITER_CENT/2,2);
		c13->Divide(NITER_CENT/2,2);
		c14->Divide(NITER_CENT/2,2);
		c15->Divide(NITER_CENT/2,2);
		c16->Divide(NITER_CENT/2,2);
		c17->Divide(NITER_CENT/2,2);
		c18->Divide(NITER_CENT/2,2);
		c19->Divide(NITER_CENT/2,2);
		c20->Divide(NITER_CENT/2,2);
		c21->Divide(NITER_CENT/2,2);
		c22->Divide(NITER_CENT/2,2);
		
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			cent_interval[iter_cent] = cent_4bins[iter_cent];
		}
	}else if(NITER_CENT == 7)
	{
		c1->Divide(4,2);
		c2->Divide(4,2);
		c3->Divide(4,2);
		c4->Divide(4,2);
		c6->Divide(4,2);
		c7->Divide(4,2);
		c8->Divide(4,2);
		c9->Divide(4,2);
		c11->Divide(4,2);
		c12->Divide(4,2);
		c13->Divide(4,2);
		c14->Divide(4,2);
		c15->Divide(4,2);
		c16->Divide(4,2);
		c17->Divide(4,2);
		c18->Divide(4,2);
		c19->Divide(4,2);
		c20->Divide(4,2);
		c21->Divide(4,2);
		c22->Divide(4,2);
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			cent_interval[iter_cent] = cent_7bins[iter_cent];
		}
	}else if(NITER_CENT == 10)
	{
		c1->Divide(NITER_CENT/2,2);
		c2->Divide(NITER_CENT/2,2);
		c3->Divide(NITER_CENT/2,2);
		c4->Divide(NITER_CENT/2,2);
		c6->Divide(NITER_CENT/2,2);
		c7->Divide(NITER_CENT/2,2);
		c8->Divide(NITER_CENT/2,2);
		c9->Divide(NITER_CENT/2,2);
		c11->Divide(NITER_CENT/2,2);
		c12->Divide(NITER_CENT/2,2);
		c13->Divide(NITER_CENT/2,2);
		c14->Divide(NITER_CENT/2,2);
		c15->Divide(NITER_CENT/2,2);
		c16->Divide(NITER_CENT/2,2);
		c17->Divide(NITER_CENT/2,2);
		c18->Divide(NITER_CENT/2,2);
		c19->Divide(NITER_CENT/2,2);
		c20->Divide(NITER_CENT/2,2);
		c21->Divide(NITER_CENT/2,2);
		c22->Divide(NITER_CENT/2,2);
		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			cent_interval[iter_cent] = cent_10bins[iter_cent];
		}
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
	TH1F *NPositive_Lambda = (TH1F*) myFile_data->Get("Lambda Positive Polarization");
	TH1F *NNegative_Lambda = (TH1F*) myFile_data->Get("Lambda Negative Polarization");
	TH1F *NPositive_ALambda = (TH1F*) myFile_data->Get("ALambda Positive Polarization");
	TH1F *NNegative_ALambda = (TH1F*) myFile_data->Get("ALambda Negative Polarization");
	
	TH1F *NLambda = (TH1F*) myFile_data->Get("Amount of Lambda");
	TH1F *NALambda = (TH1F*) myFile_data->Get("Amount of ALambda");
	TH1F *NXi = (TH1F*) myFile_data->Get("Amount of Xi");
	TH1F *NAXi = (TH1F*) myFile_data->Get("Amount of AXi");
	TH1F *NXi0 = (TH1F*) myFile_data->Get("Amount of Xi0");
	TH1F *NAXi0 = (TH1F*) myFile_data->Get("Amount of AXi0");
	TH1F *NSig0 = (TH1F*) myFile_data->Get("Amount of Sig0");
	TH1F *NASig0 = (TH1F*) myFile_data->Get("Amount of ASig0");
	
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
		Lpolar_y_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_Lambda_%d",iter_cent));
		Lpolar_y_prim_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_Lambda_%d",iter_cent));	
		
		Lpolar_y_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_ALambda_%d",iter_cent));
		Lpolar_y_prim_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_ALambda_%d",iter_cent));	
		
		Lpolar_y_Xi[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_Xi_%d",iter_cent));
		Lpolar_y_prim_Xi[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_Xi_%d",iter_cent));
		
		Lpolar_y_AXi[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_AXi_%d",iter_cent));
		Lpolar_y_prim_AXi[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_AXi_%d",iter_cent));
		
		Lpolar_y_Xi0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_Xi0_%d",iter_cent));
		Lpolar_y_prim_Xi0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_Xi0_%d",iter_cent));
		
		Lpolar_y_AXi0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_AXi0_%d",iter_cent));
		Lpolar_y_prim_AXi0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_AXi0_%d",iter_cent));
		
		Lpolar_y_Sig0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_Sig0_%d",iter_cent));
		Lpolar_y_prim_Sig0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_Sig0_%d",iter_cent));
		
		Lpolar_y_ASig0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_ASig0_%d",iter_cent));
		Lpolar_y_prim_ASig0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_prim_ASig0_%d",iter_cent));
		
		Lpolar_y_Lambda_from_Xi[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_Lambda_from_Xi_%d",iter_cent));
		Lpolar_y_ALambda_from_AXi[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_ALambda_from_AXi_%d",iter_cent));
		Lpolar_y_Lambda_from_Xi0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_Lambda_from_Xi0_%d",iter_cent));
		Lpolar_y_ALambda_from_AXi0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_ALambda_from_AXi0_%d",iter_cent));
		Lpolar_y_Lambda_from_Sig0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_Lambda_from_Sig0_%d",iter_cent));
		Lpolar_y_ALambda_from_ASig0[iter_cent] = (TH1D*) myFile_data->Get(Form("LPolar_y_ALambda_from_ASig0_%d",iter_cent));
		
		CosTheta_hist_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("CosTheta_hist_Lambda_%d",iter_cent));
		Pstar_hist_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("Pstar_hist_Lambda_%d",iter_cent));
		PstarRP_hist_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_Lambda_%d",iter_cent));
		
		CosTheta_hist_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("CosTheta_hist_ALambda_%d",iter_cent));
		Pstar_hist_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("Pstar_hist_ALambda_%d",iter_cent));
		PstarRP_hist_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_ALambda_%d",iter_cent));		
		
		CosTheta_hist_prim_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("CosTheta_hist_prim_Lambda_%d",iter_cent));
		Pstar_hist_prim_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("Pstar_hist_prim_Lambda_%d",iter_cent));
		PstarRP_hist_prim_Lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_prim_Lambda_%d",iter_cent));	
		
		CosTheta_hist_prim_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("CosTheta_hist_prim_ALambda_%d",iter_cent));
		Pstar_hist_prim_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("Pstar_hist_prim_ALambda_%d",iter_cent));
		PstarRP_hist_prim_ALambda[iter_cent] = (TH1D*) myFile_data->Get(Form("PstarRP_hist_prim_ALambda_%d",iter_cent));	
		
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		xmin_CosTheta[iter_cent] = CosTheta_hist_Lambda[iter_cent]->GetXaxis()->GetXmin();
		xmax_CosTheta[iter_cent] = CosTheta_hist_Lambda[iter_cent]->GetXaxis()->GetXmax();
		xmin_Pstar[iter_cent] = PstarRP_hist_Lambda[iter_cent]->GetXaxis()->GetXmin();
		xmax_Pstar[iter_cent] = PstarRP_hist_Lambda[iter_cent]->GetXaxis()->GetXmax();
		
		//create fitting functions for each histogram (delta(phi)_{RP} distributions):		
		
		//For Lambda:
		sprintf(int_fnc_PstarRP_Lambda,"fitting_fnc_PstarRP_Lambda_%d",iter_cent);
		sprintf(int_fnc_PstarRP_prim_Lambda,"fitting_fnc_PstarRP_prim_Lambda_%d",iter_cent);
		if (fitting_choice == "2orders")
		{
			fnc_PstarRP_Lambda[iter_cent] = new TF1(int_fnc_PstarRP_Lambda,fitting_fnc_2orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
			fnc_PstarRP_prim_Lambda[iter_cent] = new TF1(int_fnc_PstarRP_prim_Lambda,fitting_fnc_2orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
		}else if (fitting_choice == "3orders")
		{
			fnc_PstarRP_Lambda[iter_cent] = new TF1(int_fnc_PstarRP_Lambda,fitting_fnc_3orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
			fnc_PstarRP_prim_Lambda[iter_cent] = new TF1(int_fnc_PstarRP_prim_Lambda,fitting_fnc_3orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
		}else
		{
			cout << "This fitting choice is not defined! Please choose either 2orders or 3orders" << endl;
			return 1;
		}
		
		fnc_PstarRP_Lambda[iter_cent]->SetLineWidth(4);
		fnc_PstarRP_Lambda[iter_cent]->SetLineColor(2); //red color for fitting line
		fnc_PstarRP_Lambda[iter_cent]->SetParameters(10000., 0.01, 0.0003, 0.0002, 0.002);	
		fnc_PstarRP_prim_Lambda[iter_cent]->SetLineWidth(4);
		fnc_PstarRP_prim_Lambda[iter_cent]->SetLineColor(2); //red color for fitting line
		fnc_PstarRP_prim_Lambda[iter_cent]->SetParameters(10000., 0.01, 0.0003, 0.0002, 0.002);
		
		//plotting:
		c1->cd(iter_cent+1);	
		
		PstarRP_hist_Lambda[iter_cent]->SetYTitle("dN/d(#Delta#phi_{p}^{*})");
		PstarRP_hist_Lambda[iter_cent]->SetXTitle("#Delta#phi_{p}^{*}");
		PstarRP_hist_Lambda[iter_cent]->Draw("p9");
		PstarRP_hist_Lambda[iter_cent]->Fit(int_fnc_PstarRP_Lambda);
		
		//PstarRP_hist_Lambda[iter_cent]->Fit(int_fnc_PstarRP_Lambda,"L");	//PstarRP_hist[iter_cent]->Fit(int_fnc_PstarRP_dataset,"w","",xmin_Pstar[iter_cent],xmax_Pstar[iter_cent]);
		
		polar_par_hist_Lambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_Lambda[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		polar_par_hist_err_Lambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_Lambda[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
		TPaveText *t1_1 = new TPaveText(.13,.75,.5,.94);
		t1_1->SetTextAlign(22);
		t1_1->SetTextColor(kRed+2);
		t1_1->SetTextFont(72);
		t1_1->SetTextSize(0.04);
		t1_1->Paint("NDC");
		TText *t1_1_1 = t1_1->AddText(cent_interval[iter_cent]);
		TText *t1_1_2 = t1_1->AddText(Form("<P_{#Lambda}> (MCTracks) = %.4f +/- %.4f",polar_par_hist_Lambda[iter_cent],polar_par_hist_err_Lambda[iter_cent]));
		TText *t1_1_3 = t1_1->AddText(Form("<P_{#Lambda}> (full MC) = %.4f",-100.*Lpolar_y_Lambda[iter_cent]->GetMean()));
		t1_1->Draw("same");

		
		c2->cd(iter_cent+1);	
				
		PstarRP_hist_prim_Lambda[iter_cent]->SetYTitle("dN/d(#Delta#phi_{p}^{*})");
		PstarRP_hist_prim_Lambda[iter_cent]->SetXTitle("#Delta#phi_{p}^{*}");
		PstarRP_hist_prim_Lambda[iter_cent]->Draw("p9");
		PstarRP_hist_prim_Lambda[iter_cent]->Fit(int_fnc_PstarRP_prim_Lambda);
		polar_par_hist_prim_Lambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_prim_Lambda[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		polar_par_hist_prim_err_Lambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_prim_Lambda[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
		
		polar_y_mean_Lambda[iter_cent] = -100.*Lpolar_y_Lambda[iter_cent]->GetMean();
		polar_y_mean_err_Lambda[iter_cent] = -100.*Lpolar_y_Lambda[iter_cent]->GetMeanError();
		polar_y_mean_prim_Lambda[iter_cent] = -100.*Lpolar_y_prim_Lambda[iter_cent]->GetMean();
		polar_y_mean_prim_err_Lambda[iter_cent] = -100.*Lpolar_y_prim_Lambda[iter_cent]->GetMeanError();
		
		TPaveText *t2_2 = new TPaveText(.13,.75,.5,.94);
		t2_2->SetTextAlign(22);
		t2_2->SetTextColor(kRed+2);
		t2_2->SetTextFont(72);
		t2_2->SetTextSize(0.04);
		t2_2->Paint("NDC");
		TText *t2_2_1 = t2_2->AddText(cent_interval[iter_cent]);
		TText *t2_2_2 = t2_2->AddText(Form("<P_{#Lambda}> (MCTracks) = %.4f +/- %.4f",polar_par_hist_prim_Lambda[iter_cent],polar_par_hist_prim_err_Lambda[iter_cent]));
		TText *t2_2_3 = t2_2->AddText(Form("<P_{#Lambda}> (prim MC) = %.4f",-100.*Lpolar_y_prim_Lambda[iter_cent]->GetMean()));
		t2_2->Draw("same");
		
		c3->cd(iter_cent+1);	
		
		Lpolar_y_Lambda[iter_cent]->Draw();
		t1_1->Draw("same");
		
		TPaveText *t3_1 = new TPaveText(.7,.75,.94,.94);
		t3_1->SetTextAlign(22);
		t3_1->SetTextColor(kRed+2);
		t3_1->SetTextFont(72);
		t3_1->SetTextSize(0.04);
		t3_1->Paint("NDC");
		TText *t3_1_1 = t3_1->AddText(Form("N_{#Lambda} = %1.2e",Lpolar_y_Lambda[iter_cent]->GetEntries()));
		t3_1->Draw("same");
		
		c4->cd(iter_cent+1);	
		
		Lpolar_y_prim_Lambda[iter_cent]->Draw();
		t2_2->Draw("same");		
		
		TPaveText *t4_1 = new TPaveText(.7,.75,.94,.94);
		t4_1->SetTextAlign(22);
		t4_1->SetTextColor(kRed+2);
		t4_1->SetTextFont(72);
		t4_1->SetTextSize(0.04);
		t4_1->Paint("NDC");
		TText *t4_1_1 = t4_1->AddText(Form("N_{#Lambda} = %1.2e",Lpolar_y_prim_Lambda[iter_cent]->GetEntries()));
		t4_1->Draw("same");
		
		//Now for ALambda:
		sprintf(int_fnc_PstarRP_ALambda,"fitting_fnc_PstarRP_ALambda_%d",iter_cent);
		sprintf(int_fnc_PstarRP_prim_ALambda,"fitting_fnc_PstarRP_prim_ALambda_%d",iter_cent);
		if (fitting_choice == "2orders")
		{
			fnc_PstarRP_ALambda[iter_cent] = new TF1(int_fnc_PstarRP_ALambda,fitting_fnc_2orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
			fnc_PstarRP_prim_ALambda[iter_cent] = new TF1(int_fnc_PstarRP_prim_ALambda,fitting_fnc_2orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
		}else if (fitting_choice == "3orders")
		{
			fnc_PstarRP_ALambda[iter_cent] = new TF1(int_fnc_PstarRP_ALambda,fitting_fnc_3orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
			fnc_PstarRP_prim_ALambda[iter_cent] = new TF1(int_fnc_PstarRP_prim_ALambda,fitting_fnc_3orders,xmin_Pstar[iter_cent],xmax_Pstar[iter_cent],5);
		}else
		{
			cout << "This fitting choice is not defined! Please choose either 2orders or 3orders" << endl;
			return 1;
		}
		
		fnc_PstarRP_ALambda[iter_cent]->SetLineWidth(4);
		fnc_PstarRP_ALambda[iter_cent]->SetLineColor(2); //red color for fitting line
		fnc_PstarRP_ALambda[iter_cent]->SetParameters(10000., 0.01, 0.0003, 0.0002, 0.002);	
		fnc_PstarRP_prim_ALambda[iter_cent]->SetLineWidth(4);
		fnc_PstarRP_prim_ALambda[iter_cent]->SetLineColor(2); //red color for fitting line
		fnc_PstarRP_prim_ALambda[iter_cent]->SetParameters(10000., 0.01, 0.0003, 0.0002, 0.002);
		
		//plotting:
		c6->cd(iter_cent+1);	
		
		PstarRP_hist_ALambda[iter_cent]->SetYTitle("dN/d(#Delta#phi_{p}^{*})");
		PstarRP_hist_ALambda[iter_cent]->SetXTitle("#Delta#phi_{p}^{*}");
		PstarRP_hist_ALambda[iter_cent]->Draw("p9");
		PstarRP_hist_ALambda[iter_cent]->Fit(int_fnc_PstarRP_ALambda);
		
		//PstarRP_hist_Lambda[iter_cent]->Fit(int_fnc_PstarRP_Lambda,"L");	//PstarRP_hist[iter_cent]->Fit(int_fnc_PstarRP_dataset,"w","",xmin_Pstar[iter_cent],xmax_Pstar[iter_cent]);
		
		polar_par_hist_ALambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_ALambda[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		polar_par_hist_err_ALambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_ALambda[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
		TPaveText *t6_1 = new TPaveText(.13,.75,.5,.94);
		t6_1->SetTextAlign(22);
		t6_1->SetTextColor(kRed+2);
		t6_1->SetTextFont(72);
		t6_1->SetTextSize(0.04);
		t6_1->Paint("NDC");
		TText *t6_1_1 = t6_1->AddText(cent_interval[iter_cent]);
		TText *t6_1_2 = t6_1->AddText(Form("<P_{#Lambda}> (MCTracks) = %.4f +/- %.4f",polar_par_hist_ALambda[iter_cent],polar_par_hist_err_ALambda[iter_cent]));
		TText *t6_1_3 = t6_1->AddText(Form("<P_{#Lambda}> (full MC) = %.4f",-100.*Lpolar_y_ALambda[iter_cent]->GetMean()));
		t6_1->Draw("same");

		
		c7->cd(iter_cent+1);	
				
		PstarRP_hist_prim_ALambda[iter_cent]->SetYTitle("dN/d(#Delta#phi_{p}^{*})");
		PstarRP_hist_prim_ALambda[iter_cent]->SetXTitle("#Delta#phi_{p}^{*}");
		PstarRP_hist_prim_ALambda[iter_cent]->Draw("p9");
		PstarRP_hist_prim_ALambda[iter_cent]->Fit(int_fnc_PstarRP_prim_ALambda);
		polar_par_hist_prim_ALambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_prim_ALambda[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		polar_par_hist_prim_err_ALambda[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_PstarRP_prim_ALambda[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
		
		polar_y_mean_ALambda[iter_cent] = -100.*Lpolar_y_ALambda[iter_cent]->GetMean();
		polar_y_mean_err_ALambda[iter_cent] = -100.*Lpolar_y_ALambda[iter_cent]->GetMeanError();
		polar_y_mean_prim_ALambda[iter_cent] = -100.*Lpolar_y_prim_ALambda[iter_cent]->GetMean();
		polar_y_mean_prim_err_ALambda[iter_cent] = -100.*Lpolar_y_prim_ALambda[iter_cent]->GetMeanError();
		
		TPaveText *t7_2 = new TPaveText(.13,.75,.5,.94);
		t7_2->SetTextAlign(22);
		t7_2->SetTextColor(kRed+2);
		t7_2->SetTextFont(72);
		t7_2->SetTextSize(0.04);
		t7_2->Paint("NDC");
		TText *t7_2_1 = t7_2->AddText(cent_interval[iter_cent]);
		TText *t7_2_2 = t7_2->AddText(Form("<P_{#bar#Lambda}> (MCTracks) = %.4f +/- %.4f",polar_par_hist_prim_ALambda[iter_cent],polar_par_hist_prim_err_Lambda[iter_cent]));
		TText *t7_2_3 = t7_2->AddText(Form("<P_{#bar#Lambda}> (prim MC) = %.4f",-100.*Lpolar_y_prim_ALambda[iter_cent]->GetMean()));
		t7_2->Draw("same");
		
		c8->cd(iter_cent+1);	
		
		Lpolar_y_ALambda[iter_cent]->Draw();
		t6_1->Draw("same");
		
		TPaveText *t8_1 = new TPaveText(.7,.75,.94,.94);
		t8_1->SetTextAlign(22);
		t8_1->SetTextColor(kRed+2);
		t8_1->SetTextFont(72);
		t8_1->SetTextSize(0.04);
		t8_1->Paint("NDC");
		TText *t8_1_1 = t8_1->AddText(Form("N_{#bar#Lambda} = %1.2e",Lpolar_y_ALambda[iter_cent]->GetEntries()));
		t8_1->Draw("same");
		
		c9->cd(iter_cent+1);	
		
		Lpolar_y_prim_ALambda[iter_cent]->Draw();
		t7_2->Draw("same");		
		
		TPaveText *t9_1 = new TPaveText(.7,.75,.94,.94);
		t9_1->SetTextAlign(22);
		t9_1->SetTextColor(kRed+2);
		t9_1->SetTextFont(72);
		t9_1->SetTextSize(0.04);
		t9_1->Paint("NDC");
		TText *t9_1_1 = t9_1->AddText(Form("N_{#bar#Lambda} = %1.2e",Lpolar_y_prim_ALambda[iter_cent]->GetEntries()));
		t9_1->Draw("same");
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

	double star_pol_par_Lambda[] = {1.179}; // corrected value for 11.5GeV
	double star_pol_err_Lambda[] = {0.347}; // corrected error
	
	double star_pol_par_nores_Lambda[] = {0.693}; // corrected value for 11.5GeV (nores correction)
	double star_pol_err_nores_Lambda[] = {0.202}; // corrected error (nores correction)
	
	double star_pol_par_ALambda[] = {1.580}; // corrected value for 11.5GeV
	double star_pol_err_ALambda[] = {1.106}; // corrected error
	
	double star_pol_par_nores_ALambda[] = {0.921}; // corrected value for 11.5GeV (nores correction)
	double star_pol_err_nores_ALambda[] = {0.649}; // corrected error (nores correction)
						
	TGraphErrors *Polar_hist_STAR_Lambda = new TGraphErrors(1, star_cent_bin, star_pol_par_Lambda, star_cent_bin_err, star_pol_err_Lambda);
	Polar_hist_STAR_Lambda->SetName("Polar_STAR Lambda");
	Polar_hist_STAR_Lambda->SetTitle("Polar_STAR Lambda");
	Polar_hist_STAR_Lambda->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_hist_STAR_Lambda->SetLineColor(kBlue+1);
	Polar_hist_STAR_Lambda->SetMarkerColor(kBlue+1);
	Polar_hist_STAR_Lambda->SetMarkerSize(2);
	Polar_hist_STAR_Lambda->SetMarkerStyle(29);	
	
	TGraphErrors *Polar_hist_STAR_nores_Lambda = new TGraphErrors(1, star_cent_bin, star_pol_par_nores_Lambda, star_cent_bin_err, star_pol_err_nores_Lambda);
	Polar_hist_STAR_nores_Lambda->SetName("Polar_STAR Lambda");
	Polar_hist_STAR_nores_Lambda->SetTitle("Polar_STAR Lambda");
	Polar_hist_STAR_nores_Lambda->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_hist_STAR_nores_Lambda->SetLineColor(kRed+1);
	Polar_hist_STAR_nores_Lambda->SetMarkerColor(kRed+1);
	Polar_hist_STAR_nores_Lambda->SetMarkerSize(2);
	Polar_hist_STAR_nores_Lambda->SetMarkerStyle(29);	
	
	TGraphErrors *Polar_hist_STAR_ALambda = new TGraphErrors(1, star_cent_bin, star_pol_par_ALambda, star_cent_bin_err, star_pol_err_ALambda);
	Polar_hist_STAR_ALambda->SetName("Polar_STAR ALambda");
	Polar_hist_STAR_ALambda->SetTitle("Polar_STAR ALambda");
	Polar_hist_STAR_ALambda->GetYaxis()->SetTitle("P_{#bar#Lambda}, [%]");
	Polar_hist_STAR_ALambda->SetLineColor(kBlue+1);
	Polar_hist_STAR_ALambda->SetMarkerColor(kBlue+1);
	Polar_hist_STAR_ALambda->SetMarkerSize(2);
	Polar_hist_STAR_ALambda->SetMarkerStyle(29);	
	
	TGraphErrors *Polar_hist_STAR_nores_ALambda = new TGraphErrors(1, star_cent_bin, star_pol_par_nores_ALambda, star_cent_bin_err, star_pol_err_nores_ALambda);
	Polar_hist_STAR_nores_ALambda->SetName("Polar_STAR ALambda");
	Polar_hist_STAR_nores_ALambda->SetTitle("Polar_STAR ALambda");
	Polar_hist_STAR_nores_ALambda->GetYaxis()->SetTitle("P_{#bar#Lambda}, [%]");
	Polar_hist_STAR_nores_ALambda->SetLineColor(kRed+1);
	Polar_hist_STAR_nores_ALambda->SetMarkerColor(kRed+1);
	Polar_hist_STAR_nores_ALambda->SetMarkerSize(2);
	Polar_hist_STAR_nores_ALambda->SetMarkerStyle(29);	
	
	TGraphErrors *Polar_mean_MC_Lambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_y_mean_Lambda, noErr, polar_y_mean_err_Lambda);
	Polar_mean_MC_Lambda->SetName("Polar_mean_MC");
	Polar_mean_MC_Lambda->SetTitle("Polar_mean_MC");
	Polar_mean_MC_Lambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MC_Lambda->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MC_Lambda->SetLineColor(kRed+1);
	Polar_mean_MC_Lambda->SetMarkerColor(kRed+1);
	Polar_mean_MC_Lambda->SetMarkerSize(2);
	Polar_mean_MC_Lambda->SetMarkerStyle(24);	

	TGraphErrors *Polar_mean_MC_prim_Lambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_y_mean_prim_Lambda, noErr, polar_y_mean_prim_err_Lambda);
	Polar_mean_MC_prim_Lambda->SetName("Polar_mean_MC");
	Polar_mean_MC_prim_Lambda->SetTitle("Polar_mean_MC");
	Polar_mean_MC_prim_Lambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MC_prim_Lambda->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MC_prim_Lambda->SetLineColor(kRed+1);
	Polar_mean_MC_prim_Lambda->SetMarkerColor(kRed+1);
	Polar_mean_MC_prim_Lambda->SetMarkerSize(2);
	Polar_mean_MC_prim_Lambda->SetMarkerStyle(26);
	
	TGraphErrors *Polar_mean_MCtracks_Lambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_hist_Lambda, noErr, polar_par_hist_err_Lambda);
	Polar_mean_MCtracks_Lambda->SetName("Polar_mean_MCtracks");
	Polar_mean_MCtracks_Lambda->SetTitle("Polar_mean_MCtracks");
	Polar_mean_MCtracks_Lambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MCtracks_Lambda->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MCtracks_Lambda->SetLineColor(kBlack);
	Polar_mean_MCtracks_Lambda->SetMarkerColor(kBlack);
	Polar_mean_MCtracks_Lambda->SetMarkerSize(2);
	Polar_mean_MCtracks_Lambda->SetMarkerStyle(20);
	
	TGraphErrors *Polar_mean_MCtracks_prim_Lambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_hist_prim_Lambda, noErr, polar_par_hist_prim_err_Lambda);
	Polar_mean_MCtracks_prim_Lambda->SetName("Polar_mean_MCtracks");
	Polar_mean_MCtracks_prim_Lambda->SetTitle("Polar_mean_MCtracks");
	Polar_mean_MCtracks_prim_Lambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MCtracks_prim_Lambda->GetYaxis()->SetTitle("P_{#Lambda}, [%]");
	Polar_mean_MCtracks_prim_Lambda->SetLineColor(kBlack);
	Polar_mean_MCtracks_prim_Lambda->SetMarkerColor(kBlack);
	Polar_mean_MCtracks_prim_Lambda->SetMarkerSize(2);
	Polar_mean_MCtracks_prim_Lambda->SetMarkerStyle(22);  //for primary
	
	TGraphErrors *Polar_mean_MC_ALambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_y_mean_ALambda, noErr, polar_y_mean_err_ALambda);
	Polar_mean_MC_ALambda->SetName("Polar_mean_MC");
	Polar_mean_MC_ALambda->SetTitle("Polar_mean_MC");
	Polar_mean_MC_ALambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MC_ALambda->GetYaxis()->SetTitle("P_{#bar#Lambda}, [%]");
	Polar_mean_MC_ALambda->SetLineColor(kRed+1);
	Polar_mean_MC_ALambda->SetMarkerColor(kRed+1);
	Polar_mean_MC_ALambda->SetMarkerSize(2);
	Polar_mean_MC_ALambda->SetMarkerStyle(24);	

	TGraphErrors *Polar_mean_MC_prim_ALambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_y_mean_prim_ALambda, noErr, polar_y_mean_prim_err_ALambda);
	Polar_mean_MC_prim_ALambda->SetName("Polar_mean_MC");
	Polar_mean_MC_prim_ALambda->SetTitle("Polar_mean_MC");
	Polar_mean_MC_prim_ALambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MC_prim_ALambda->GetYaxis()->SetTitle("P_{#bar#Lambda}, [%]");
	Polar_mean_MC_prim_ALambda->SetLineColor(kRed+1);
	Polar_mean_MC_prim_ALambda->SetMarkerColor(kRed+1);
	Polar_mean_MC_prim_ALambda->SetMarkerSize(2);
	Polar_mean_MC_prim_ALambda->SetMarkerStyle(26);
	
	TGraphErrors *Polar_mean_MCtracks_ALambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_hist_ALambda, noErr, polar_par_hist_err_ALambda);
	Polar_mean_MCtracks_ALambda->SetName("Polar_mean_MCtracks");
	Polar_mean_MCtracks_ALambda->SetTitle("Polar_mean_MCtracks");
	Polar_mean_MCtracks_ALambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MCtracks_ALambda->GetYaxis()->SetTitle("P_{#bar#Lambda}, [%]");
	Polar_mean_MCtracks_ALambda->SetLineColor(kBlack);
	Polar_mean_MCtracks_ALambda->SetMarkerColor(kBlack);
	Polar_mean_MCtracks_ALambda->SetMarkerSize(2);
	Polar_mean_MCtracks_ALambda->SetMarkerStyle(20);
	
	TGraphErrors *Polar_mean_MCtracks_prim_ALambda = new TGraphErrors(NITER_CENT, centrality_bin, polar_par_hist_prim_ALambda, noErr, polar_par_hist_prim_err_ALambda);
	Polar_mean_MCtracks_prim_ALambda->SetName("Polar_mean_MCtracks");
	Polar_mean_MCtracks_prim_ALambda->SetTitle("Polar_mean_MCtracks");
	Polar_mean_MCtracks_prim_ALambda->GetXaxis()->SetTitle("Centrality, [%]");
	Polar_mean_MCtracks_prim_ALambda->GetYaxis()->SetTitle("P_{#bar#Lambda}, [%]");
	Polar_mean_MCtracks_prim_ALambda->SetLineColor(kBlack);
	Polar_mean_MCtracks_prim_ALambda->SetMarkerColor(kBlack);
	Polar_mean_MCtracks_prim_ALambda->SetMarkerSize(2);
	Polar_mean_MCtracks_prim_ALambda->SetMarkerStyle(22);  //for primary
	
	TLine *line = new TLine(0.,0.,80.,0.);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw("same");	
	
	//Mean Polarization for Lambda:
	
	c5->cd();	
	Polar_mean_MC_Lambda->GetYaxis()->SetRangeUser(-1.0,6.);
	Polar_mean_MC_Lambda->Draw("ap");
	Polar_mean_MCtracks_Lambda->Draw("psame");
	Polar_mean_MCtracks_prim_Lambda->Draw("psame");
	Polar_mean_MC_prim_Lambda->Draw("psame");
	Polar_hist_STAR_Lambda->Draw("psame");
	Polar_hist_STAR_nores_Lambda->Draw("psame");
	
	TLegend *legend1_1=new TLegend(0.15,0.65,0.35,0.9);
	legend1_1->SetTextFont(72);
	legend1_1->SetTextSize(0.04);
	legend1_1->SetBorderSize(0);
	legend1_1->AddEntry(Polar_mean_MC_Lambda,"MC (full)","p");
	legend1_1->AddEntry(Polar_mean_MCtracks_Lambda,"MCTracks (full)","p");	
	legend1_1->AddEntry(Polar_mean_MC_prim_Lambda,"MC (primary)","p");
	legend1_1->AddEntry(Polar_mean_MCtracks_prim_Lambda,"MCTracks (primary)","p");
	legend1_1->AddEntry(Polar_hist_STAR_Lambda,"STAR (11.5GeV)","p");
	legend1_1->AddEntry(Polar_hist_STAR_nores_Lambda,"STAR (w/o res.)","p");
	legend1_1->Draw("same");
	
	line->Draw("same");	
	
	//Mean Polarization for ALambda:
	
	//need to change, that it doesn't affect the other plots!
	c10->cd();	
	Polar_mean_MC_ALambda->GetYaxis()->SetRangeUser(-1.0,6.);
	Polar_mean_MC_ALambda->GetYaxis()->SetTitle("P_{H}, [%]");
	Polar_mean_MC_ALambda->Draw("ap");
	Polar_mean_MCtracks_ALambda->Draw("psame");
	Polar_mean_MCtracks_prim_ALambda->Draw("psame");
	Polar_mean_MC_prim_ALambda->Draw("psame");
	Polar_hist_STAR_ALambda->Draw("psame");
	Polar_hist_STAR_nores_ALambda->Draw("psame");
	
	TLegend *legend10_1=new TLegend(0.15,0.65,0.35,0.9);
	legend10_1->SetTextFont(72);
	legend10_1->SetTextSize(0.04);
	legend10_1->SetBorderSize(0);
	legend10_1->AddEntry(Polar_mean_MC_ALambda,"MC (full)","p");
	legend10_1->AddEntry(Polar_mean_MCtracks_ALambda,"MCTracks (full)","p");	
	legend10_1->AddEntry(Polar_mean_MC_prim_ALambda,"MC (primary)","p");
	legend10_1->AddEntry(Polar_mean_MCtracks_prim_ALambda,"MCTracks (primary)","p");
	legend10_1->AddEntry(Polar_hist_STAR_ALambda,"STAR (11.5GeV)","p");
	legend10_1->AddEntry(Polar_hist_STAR_nores_ALambda,"STAR (w/o res.)","p");
	legend10_1->Draw("same");
	
	line->Draw("same");	
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		c11->cd(iter_cent+1);	
		
		Lpolar_y_Xi[iter_cent]->Draw();
		Lpolar_y_Lambda_from_Xi[iter_cent]->SetLineColor(kRed);	
		Lpolar_y_Lambda_from_Xi[iter_cent]->Draw("same");
		
		TPaveText *t11_1 = new TPaveText(.13,.75,.5,.94);
		t11_1->SetTextAlign(22);
		t11_1->SetTextColor(kRed+2);
		t11_1->SetTextFont(72);
		t11_1->SetTextSize(0.04);
		t11_1->Paint("NDC");
		TText *t11_1_1 = t11_1->AddText(cent_interval[iter_cent]);
		TText *t11_1_2 = t11_1->AddText(Form("<P_{#Lambda}> (#Xi) = %.4f",-100.*Lpolar_y_Xi[iter_cent]->GetMean()));
		TText *t11_1_3 = t11_1->AddText(Form("<P_{#Lambda}> (#Lambda from #Xi) = %.4f",-100.*Lpolar_y_Lambda_from_Xi[iter_cent]->GetMean()));
		t11_1->Draw("same");
		
		c12->cd(iter_cent+1);	
		
		Lpolar_y_prim_Xi0[iter_cent]->Draw();
		Lpolar_y_Lambda_from_Xi[iter_cent]->SetLineColor(kRed);	
		Lpolar_y_Lambda_from_Xi[iter_cent]->Draw("same");
		
		TPaveText *t12_1 = new TPaveText(.13,.75,.5,.94);
		t12_1->SetTextAlign(22);
		t12_1->SetTextColor(kRed+2);
		t12_1->SetTextFont(72);
		t12_1->SetTextSize(0.04);
		t12_1->Paint("NDC");
		TText *t12_1_1 = t12_1->AddText(cent_interval[iter_cent]);
		TText *t12_1_2 = t12_1->AddText(Form("<P_{#Lambda}> (#Xi) = %.4f",-100.*Lpolar_y_prim_Xi0[iter_cent]->GetMean()));
		TText *t12_1_3 = t12_1->AddText(Form("<P_{#Lambda}> (#Lambda from #Xi) = %.4f",-100.*Lpolar_y_Lambda_from_Xi[iter_cent]->GetMean()));
		t12_1->Draw("same");
	}
	
	c23->cd();	
	Polar_mean_MC_ALambda->GetYaxis()->SetRangeUser(-1.0,6.);
	Polar_mean_MC_ALambda->SetMarkerStyle(24);
	Polar_mean_MC_ALambda->SetMarkerColor(kBlack);
	Polar_mean_MC_ALambda->Draw("ap");
	Polar_mean_MC_prim_ALambda->SetMarkerStyle(26);
	Polar_mean_MC_prim_ALambda->Draw("psame");
	Polar_hist_STAR_ALambda->SetMarkerStyle(30);
	Polar_hist_STAR_ALambda->Draw("psame");
	
	Polar_mean_MC_Lambda->SetMarkerStyle(20);
	Polar_mean_MC_Lambda->SetMarkerColor(kBlack);
	Polar_mean_MC_Lambda->Draw("psame");
	Polar_mean_MC_prim_Lambda->SetMarkerStyle(22);
	Polar_mean_MC_prim_Lambda->Draw("psame");
	Polar_hist_STAR_Lambda->Draw("psame");
	
	TLegend *legend23_1=new TLegend(0.15,0.65,0.35,0.9);
	legend23_1->SetTextFont(72);
	legend23_1->SetTextSize(0.04);
	legend23_1->SetBorderSize(0);
	legend23_1->AddEntry(Polar_mean_MC_Lambda,"#Lambda MC (full)","p");
	legend23_1->AddEntry(Polar_mean_MC_ALambda,"#bar{#Lambda} MC (full)","p");	
	legend23_1->AddEntry(Polar_mean_MC_prim_Lambda,"#Lambda MC (primary)","p");
	legend23_1->AddEntry(Polar_mean_MC_prim_ALambda,"#bar{#Lambda} MC (primary)","p");
	legend23_1->AddEntry(Polar_hist_STAR_Lambda,"#Lambda STAR (11.5GeV)","p");
	legend23_1->AddEntry(Polar_hist_STAR_ALambda,"#bar{#Lambda} STAR (11.5GeV)","p");
	legend23_1->Draw("same");
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		float Npol_pos_Lambda = NPositive_Lambda->GetBinContent(iter_cent+1);
		float Npol_neg_Lambda = NNegative_Lambda->GetBinContent(iter_cent+1);
		float Ncentr_number = NCentr->GetBinContent(iter_cent+1);
		
		cout << "cent = " << cent_interval[iter_cent] << "; Ncentr_number = " << Ncentr_number << "; number of positive L = " << Npol_pos_Lambda << "; umber of negative L = " << Npol_neg_Lambda << endl;
		cout << "Difference = " << (TMath::Abs(Npol_pos_Lambda - Npol_neg_Lambda))/2. << endl;
		float polar_number_Lambda = 100.*TMath::Abs((Npol_pos_Lambda - Npol_neg_Lambda)/(Npol_pos_Lambda + Npol_neg_Lambda));
		//float polar_number = Npol_pos/(Npol_pos + Npol_neg);
		cout << "Polarization = " << polar_number_Lambda << endl;
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
