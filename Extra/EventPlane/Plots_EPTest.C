///Plotting distributions of EPTest

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
void Plots_EPTest(TString inFile = "Output_EPTest.root", TString generator_name = "PHSD", const int NITER_CENT = 4)
{
	TH1F *Psi_EP_diff_cent_phsd[NITER_CENT];
	TF1 *fitting_fnc_phsd[NITER_CENT];
	char *fitting_fnc_int_phsd = new char[NITER_CENT];
	double sigma_deg_phsd[NITER_CENT];
	
	//canvases for costheta distributions
	TCanvas *c1 = new TCanvas("Fitting Difference","Fitting Difference",0,0,600,600);
	TCanvas *c2 = new TCanvas("Plotting differece","Plotting differece",0,0,600,600);
	TCanvas *c3 = new TCanvas("Plotting resolution","Plotting resolution",0,0,600,600);
	
	char **cent_interval;
	if(NITER_CENT == 4)
	{
		c1->Divide(NITER_CENT/2,2);
		cent_interval = (char *[]){"0 - 10%","10 - 20%","20 - 50%","50 - 100%"};
	}else if(NITER_CENT == 7)
	{
		c1->Divide(4,2);
		cent_interval = (char *[]){"0 - 10%","10 - 20%","20 - 30%","30 - 40%","40 - 50%","50 - 60%","60 - 70%"};
	}else if(NITER_CENT == 10)
	{
		c1->Divide(NITER_CENT/2,2);
		cent_interval = (char *[]){"0 - 10%","10 - 20%","20 - 30%","30 - 40%","40 - 50%","50 - 60%","60 - 70%","70 - 80%","80 - 90%","90 - 100%"};
	}else
	{
		cout << "This centrality binning is not defined!" << endl;
		return 1;
	}
	
	int *centrality_min;
	int *centrality_max;
	double *noErr;
	double *_CentrBins;
	if (NITER_CENT == 4)
	{
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		noErr = init_double_array(4, 0, 0., 0., 0., 0.);		
		_CentrBins = init_double_array(5, 0, 0.,10.,20.,50.,100.);
	}else if (NITER_CENT == 7)
	{
		centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
		noErr = init_double_array(7, 0, 0., 0., 0., 0., 0., 0., 0.);	
		_CentrBins = init_double_array(8, 0, 0., 10., 20., 30., 40., 50., 60., 70.);	
	}else if (NITER_CENT == 10)
	{
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		noErr = init_double_array(10, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
		_CentrBins = init_double_array(11, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.);
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
	
//input file:
	TFile *myFile_data = new TFile(inFile);

	TH1F *NCentr = (TH1F*) myFile_data->Get("NCentr");	
	double NEv_cent[NITER_CENT], ResEP1_true[NITER_CENT],ResEP1_exp[NITER_CENT],SubEvRes1[NITER_CENT];
	TH1D *Resolution_EP1_true_old = (TH1D*) myFile_data->Get("Resolution_EP1_true");
	TH1D *Resolution_EP1_exp_old = (TH1D*) myFile_data->Get("Resolution_EP1_exp");
		
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Psi_EP_diff_cent_phsd[iter_cent] = (TH1F*)myFile_data->Get(Form("Psi_EP_diff_cent_%d",iter_cent));
		NEv_cent[iter_cent] = NCentr->GetBinContent(iter_cent+1);
		ResEP1_true[iter_cent] = Resolution_EP1_true_old->GetBinContent(iter_cent+1)/NEv_cent[iter_cent];
		SubEvRes1[iter_cent] = Resolution_EP1_exp_old->GetBinContent(iter_cent+1)/NEv_cent[iter_cent];
		SubEvRes1[iter_cent] = TMath::Sqrt(SubEvRes1[iter_cent]);
		ResEP1_exp[iter_cent] = ResEventPlane(TMath::Sqrt(2.)*Chi(SubEvRes1[iter_cent]));
		cout << "iter_cent = " << iter_cent << "; ResEP1_true = " << ResEP1_true[iter_cent] << "; ResEP1_exp = " << ResEP1_exp[iter_cent] << endl;
	}
		
	TH1D *Resolution_phsd_deg = new TH1D("Resolution_phsd_deg","Resolution_phsd_deg",NITER_CENT,_CentrBins);
	Resolution_phsd_deg->SetYTitle("#sigma(#Psi^{1}_{EP} - #Psi_{RP}) [deg]");
	Resolution_phsd_deg->SetXTitle("Centrality, [%]");
	Resolution_phsd_deg->SetMarkerStyle(20);
	Resolution_phsd_deg->SetLineColor(kRed+1);
	Resolution_phsd_deg->SetMarkerColor(kRed+1);
	Resolution_phsd_deg->SetMarkerSize(2);
	Resolution_phsd_deg->SetLineWidth(2);
	
	TH1D *Resolution_EP1_true = new TH1D("Resolution_EP_true","Resolution_EP_true",NITER_CENT,_CentrBins);
	Resolution_EP1_true->SetYTitle("R_{EP}^{1} = < cos(#Psi^{1}_{EP} - #Psi_{RP}) >");
	Resolution_EP1_true->SetXTitle("Centrality, [%]");
	Resolution_EP1_true->SetMarkerStyle(20);
	Resolution_EP1_true->SetLineColor(kBlack);
	Resolution_EP1_true->SetMarkerColor(kBlack);
	Resolution_EP1_true->SetMarkerSize(2);
	Resolution_EP1_true->SetLineWidth(2);
	
	TH1D *Resolution_EP1_exp = new TH1D("Resolution_EP_exp","Resolution_EP_exp",NITER_CENT,_CentrBins);
	Resolution_EP1_exp->SetYTitle("R_{EP}^{1} = < cos(#Psi^{1}_{EP} - #Psi_{RP}) >");
	Resolution_EP1_exp->SetXTitle("Centrality, [%]");
	Resolution_EP1_exp->SetMarkerStyle(23);
	Resolution_EP1_exp->SetLineColor(kRed+1);
	Resolution_EP1_exp->SetMarkerColor(kRed+1);
	Resolution_EP1_exp->SetMarkerSize(2);
	Resolution_EP1_exp->SetLineWidth(2);
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);	
		Resolution_EP1_exp->SetBinContent(iter_cent+1,ResEP1_exp[iter_cent]);	
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{	
		
		sprintf(fitting_fnc_int_phsd,"fitting_fnc_phsd_%d",iter_cent);
		fitting_fnc_phsd[iter_cent] = new TF1(fitting_fnc_int_phsd,"gaus",-50., 50.);
		fitting_fnc_phsd[iter_cent]->SetLineStyle(2);
		fitting_fnc_phsd[iter_cent]->SetLineWidth(3);
		fitting_fnc_phsd[iter_cent]->SetLineColor(kRed); //red color for fitting line
		
		gStyle->SetOptFit(1);
		c1->cd(iter_cent+1);	

		Psi_EP_diff_cent_phsd[iter_cent]->Fit(fitting_fnc_int_phsd,"","",-100.,100.);
		sigma_deg_phsd[iter_cent] = fitting_fnc_phsd[iter_cent]->GetParameter(2);		
	}	
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Resolution_phsd_deg->SetBinContent(iter_cent+1,sigma_deg_phsd[iter_cent]);
	}
	
	c2->cd();
	Resolution_phsd_deg->GetYaxis()->SetRangeUser(0.,100.);
	Resolution_phsd_deg->Draw("lp");
	
	TLegend *legend2=new TLegend(0.7,0.75,0.8,0.85);
	legend2->SetTextFont(72);
	legend2->SetTextSize(0.04);
	legend2->SetBorderSize(0);
	legend2->AddEntry(Resolution_phsd_deg,generator_name,"p");
	legend2->Draw("same");
	
	c3->cd();
	TLegend *legend3=new TLegend(0.35,0.15,0.55,0.35);
	legend3->SetTextFont(72);
	legend3->SetTextSize(0.04);
	legend3->SetBorderSize(0);
	
	Resolution_EP1_true->GetYaxis()->SetRangeUser(0.,1.);
	Resolution_EP1_true->SetYTitle("R_{EP}^{1} = < cos(#Psi^{1}_{EP} - #Psi_{RP}) >");
	Resolution_EP1_true->Draw("psame");
	Resolution_EP1_exp->Draw("psame");
	
	legend3->AddEntry(Resolution_EP1_true,"True","lp");
	legend3->AddEntry(Resolution_EP1_exp,"Reco","lp");
	legend3->Draw("same");
	legend2->Draw("same");
	
	double resolution_4intervals[4];	
	for(int iter = 0; iter < 4; iter++)
	{	
		resolution_4intervals[iter] = 0.;
	}
	if(NITER_CENT == 4)
	{	
		for(int iter = 0; iter < 4; iter++)
		{	
			resolution_4intervals[iter] = Resolution_EP1_exp->GetBinContent(iter+1);
		}
					
	}else if(NITER_CENT == 7)
	{
		resolution_4intervals[0] = Resolution_EP1_exp->GetBinContent(1);
		resolution_4intervals[1] = Resolution_EP1_exp->GetBinContent(2);
		resolution_4intervals[2] = (Resolution_EP1_exp->GetBinContent(3) + Resolution_EP1_exp->GetBinContent(4) + Resolution_EP1_exp->GetBinContent(5))/3.;
		resolution_4intervals[3] = (Resolution_EP1_exp->GetBinContent(6) + Resolution_EP1_exp->GetBinContent(7))/2.;
	}else if(NITER_CENT == 10)
	{
		resolution_4intervals[0] = Resolution_EP1_exp->GetBinContent(1);
		resolution_4intervals[1] = Resolution_EP1_exp->GetBinContent(2);
		resolution_4intervals[2] = (Resolution_EP1_exp->GetBinContent(3) + Resolution_EP1_exp->GetBinContent(4) + Resolution_EP1_exp->GetBinContent(5))/3.;
		resolution_4intervals[3] = (Resolution_EP1_exp->GetBinContent(6) + Resolution_EP1_exp->GetBinContent(7) + Resolution_EP1_exp->GetBinContent(8) + Resolution_EP1_exp->GetBinContent(9) + Resolution_EP1_exp->GetBinContent(10))/4.;
	}else
	{
		cout << "This centrality binning is not defined!" << endl;
		return 1;
	}
	
	for(int iter = 0; iter < 4; iter++)
	{
		cout << "iter = " << iter << "; resolution = " << resolution_4intervals[iter] << endl;
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
