/**
 * @file mPlotMC.C
 * @author Elizaveta Nazarova
 * @brief Plotting macro for MC GlobalPol wagon
 * @version 0.1
 * @date 2023-07-15
 * Plotting the final results of the GlobalPolarizationMC wagon, obtained using macro mAnalyzeMC.C
 * Values from STAR experiment at 20-50% centrality are shown for comparison only with the choice of NITER_CENT = 4 
 */
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
/**
 * @brief Main function for plotting MC results
 * 
 * @param inFile 			input file
 * @param NITER_CENT 		number of centrality bins
 * @param particle_pdg 		pdg of analyzed hyperon
 * @param angle_choice 		choice of angle (RP or EP)
 */
void mPlotMC(TString inFile = "Output_GlopalPolMC_Final.root", const int NITER_CENT = 4, int particle_pdg = 3122, TString angle_choice = "RP")
{
	const char *cent_interval_4bins[] = {"0 - 10 %","10 - 20 %","20 - 50 %","50 - 100 %"};
	const char *cent_interval_7bins[] = {"0 - 10 %","10 - 20 %","20 - 30 %","30 - 40 %","40 - 50 %","50 - 60 %","60 - 70 %"};
	const char *cent_interval_10bins[] = {"0 - 10 %","10 - 20 %","20 - 30 %","30 - 40 %","40 - 50 %","50 - 60 %","60 - 70 %","70 - 80 %","80 - 90 %","90 - 100 %"};
	
	TString title_Yaxis_fit = "";
	TString title_Yaxis_meanpol = "";
	if (particle_pdg == 3122)
	{
		title_Yaxis_fit = "dN_{#Lambda}/d#Delta #phi";
		title_Yaxis_meanpol = "P_{#Lambda} [%]";
	}
	else if (particle_pdg == -3122)
	{
		title_Yaxis_fit = "dN_{#bar#Lambda}/d#Delta #phi";
		title_Yaxis_meanpol = "P_{#bar#Lambda} [%]";
	}
	else
	{
		cout << "This particle_pdg is not defined!" << endl;
		return;
	}
	
	// create canvas and divide them depending on the number of centrality bins
	TCanvas *c1 = new TCanvas("Fitting Full","Fitting Full",0,0,600,600);
	TCanvas *c2 = new TCanvas("Fitting Primary","Fitting Primary",0,0,600,600);
	TCanvas *c3 = new TCanvas("P_{y} distributions (full)","P_{y} distributions (full)",0,0,600,600);
	TCanvas *c4 = new TCanvas("P_{y} distributions (primary)","P_{y} distributions (primary)",0,0,600,600);
	TCanvas *c5 = new TCanvas("Mean polarization (full)","Mean polarization (full)",0,0,600,600);
	TCanvas *c6 = new TCanvas("Mean polarization (all)","Mean polarization (all)",0,0,600,600);
	
	if(NITER_CENT == 4 || NITER_CENT == 10)
	{
		c1->Divide(NITER_CENT/2,2);
		c2->Divide(NITER_CENT/2,2);
		c3->Divide(NITER_CENT/2,2);
		c4->Divide(NITER_CENT/2,2);
	}
	else if(NITER_CENT == 7)
	{
		c1->Divide(4,2);
		c2->Divide(4,2);
		c3->Divide(4,2);
		c4->Divide(4,2);
	}
	else
	{
		cout << "This centrality binning is not defined!" << endl;
		return;
	}
	
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	gROOT->ForceStyle();
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	// input file:
	TFile *inFile_data = new TFile(inFile);

	// Get the histograms from the input file
	TH1D *hResolution_EP1_true, *hResolution_EP1_reco;
	double ResEP1_true[NITER_CENT], ResEP1_exp[NITER_CENT];
	hResolution_EP1_true = (TH1D*) inFile_data->Get("hResolution_EP1_true");
	hResolution_EP1_reco = (TH1D*) inFile_data->Get("hResolution_EP1_reco");
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		ResEP1_true[iter_cent] = hResolution_EP1_true->GetBinContent(iter_cent+1);
		ResEP1_exp[iter_cent] = hResolution_EP1_reco->GetBinContent(iter_cent+1);
		cout << "iter_cent = " << iter_cent << "; ResEP1_true = " << ResEP1_true[iter_cent] << "; ResEP1_exp = " << ResEP1_exp[iter_cent] << endl;
	}
	
	TH1D *hPolarY_Full[NITER_CENT], *hPolarY_Prim[NITER_CENT];
	TH1D *hDeltaPhiRP_Full[NITER_CENT], *hDeltaPhiRP_Prim[NITER_CENT], *hDeltaPhiEP_Full[NITER_CENT], *hDeltaPhiEP_Prim[NITER_CENT];
	TF1 *Fit_hDeltaPhiRP_Full[NITER_CENT], *Fit_hDeltaPhiRP_Prim[NITER_CENT], *Fit_hDeltaPhiEP_Full[NITER_CENT], *Fit_hDeltaPhiEP_Prim[NITER_CENT];
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		hPolarY_Full[iter_cent] = (TH1D*) inFile_data->Get(Form("hPolarY_Full_cent%d",iter_cent));
		hPolarY_Full[iter_cent]->SetYTitle("Counts");
		hPolarY_Full[iter_cent]->SetXTitle("P_{y}");
		hPolarY_Prim[iter_cent] = (TH1D*) inFile_data->Get(Form("hPolarY_Prim_cent%d",iter_cent));	
		hPolarY_Prim[iter_cent]->SetYTitle("Counts");
		hPolarY_Prim[iter_cent]->SetXTitle("P_{y}");	
		
		hDeltaPhiRP_Full[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiRP_Full_cent%d",iter_cent));	
		hDeltaPhiRP_Full[iter_cent]->SetYTitle(title_Yaxis_fit);
		hDeltaPhiRP_Full[iter_cent]->SetXTitle("#Delta #phi^{*}_{p}");
		hDeltaPhiRP_Full[iter_cent]->SetMarkerStyle(21);
		hDeltaPhiRP_Full[iter_cent]->SetLineColor(kBlack);
		hDeltaPhiRP_Full[iter_cent]->SetMarkerSize(1.5);
		
		hDeltaPhiRP_Prim[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiRP_Prim_cent%d",iter_cent));
		hDeltaPhiRP_Prim[iter_cent]->SetYTitle(title_Yaxis_fit);
		hDeltaPhiRP_Prim[iter_cent]->SetXTitle("#Delta #phi^{*}_{p}");
		hDeltaPhiRP_Prim[iter_cent]->SetMarkerStyle(21);
		hDeltaPhiRP_Prim[iter_cent]->SetLineColor(kBlack);
		hDeltaPhiRP_Prim[iter_cent]->SetMarkerSize(1.5);
	
		hDeltaPhiEP_Full[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiEP_Full_cent%d",iter_cent));	
		hDeltaPhiEP_Full[iter_cent]->SetYTitle(title_Yaxis_fit);
		hDeltaPhiEP_Full[iter_cent]->SetXTitle("#Delta #phi^{*}_{p}");
		hDeltaPhiEP_Full[iter_cent]->SetMarkerStyle(21);
		hDeltaPhiEP_Full[iter_cent]->SetLineColor(kBlack);
		hDeltaPhiEP_Full[iter_cent]->SetMarkerSize(1.5);
		
		hDeltaPhiEP_Prim[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiEP_Prim_cent%d",iter_cent));
		hDeltaPhiEP_Prim[iter_cent]->SetYTitle(title_Yaxis_fit);
		hDeltaPhiEP_Prim[iter_cent]->SetXTitle("#Delta #phi^{*}_{p}");
		hDeltaPhiEP_Prim[iter_cent]->SetMarkerStyle(21);
		hDeltaPhiEP_Prim[iter_cent]->SetLineColor(kBlack);
		hDeltaPhiEP_Prim[iter_cent]->SetMarkerSize(1.5);
		
		
	}
	
	// calculate amount of Lambda in each centrality bin
	double n_L_hDeltaPhiRP_Full[NITER_CENT], n_L_hDeltaPhiRP_Prim[NITER_CENT], n_L_hDeltaPhiEP_Full[NITER_CENT], n_L_hDeltaPhiEP_Prim[NITER_CENT];
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		n_L_hDeltaPhiRP_Full[iter_cent] = 0.;
		n_L_hDeltaPhiRP_Prim[iter_cent] = 0.;
		n_L_hDeltaPhiEP_Full[iter_cent] = 0.;
		n_L_hDeltaPhiEP_Prim[iter_cent] = 0.;
		for(int iter_bin = 0; iter_bin < hDeltaPhiRP_Full[iter_cent]->GetNbinsX(); iter_bin++)
		{
			n_L_hDeltaPhiRP_Full[iter_cent] += hDeltaPhiRP_Full[iter_cent]->GetBinContent(iter_bin+1);
			n_L_hDeltaPhiRP_Prim[iter_cent] += hDeltaPhiRP_Prim[iter_cent]->GetBinContent(iter_bin+1);
			n_L_hDeltaPhiEP_Full[iter_cent] += hDeltaPhiEP_Full[iter_cent]->GetBinContent(iter_bin+1);
			n_L_hDeltaPhiEP_Prim[iter_cent] += hDeltaPhiEP_Prim[iter_cent]->GetBinContent(iter_bin+1);
		}
	}
	
	// Get the graphs from the input file
	TGraphErrors *Polar_STAR, *Polar_STAR_nores, *Polar_MC_Full, *Polar_MC_Prim, *Polar_Reco_MCRP_Full, *Polar_Reco_MCRP_Prim, *Polar_Reco_MCEP_Full, *Polar_Reco_MCEP_Prim;
	
	Polar_STAR = (TGraphErrors*) inFile_data->Get("Polar_STAR");
	Polar_STAR_nores = (TGraphErrors*) inFile_data->Get("Polar_STAR_nores");
	Polar_MC_Full = (TGraphErrors*) inFile_data->Get("Polar_MC_Full");
	Polar_MC_Prim = (TGraphErrors*) inFile_data->Get("Polar_MC_Prim");
	Polar_Reco_MCRP_Full = (TGraphErrors*) inFile_data->Get("Polar_Reco_MCRP_Full");
	Polar_Reco_MCRP_Prim = (TGraphErrors*) inFile_data->Get("Polar_Reco_MCRP_Prim");
	Polar_Reco_MCEP_Full = (TGraphErrors*) inFile_data->Get("Polar_Reco_MCEP_Full");
	Polar_Reco_MCEP_Prim = (TGraphErrors*) inFile_data->Get("Polar_Reco_MCEP_Prim");
	
	Polar_MC_Full->GetXaxis()->SetTitle("Centrality [%]");
	Polar_MC_Full->GetYaxis()->SetTitle(title_Yaxis_meanpol);
	Polar_MC_Full->SetLineColor(kRed+1);
	Polar_MC_Full->SetMarkerColor(kRed+1);
	Polar_MC_Full->SetMarkerSize(2);
	Polar_MC_Full->SetMarkerStyle(24);	
	
	Polar_Reco_MCRP_Full->SetLineColor(kBlack);
	Polar_Reco_MCRP_Full->SetMarkerColor(kBlack);
	Polar_Reco_MCRP_Full->SetMarkerSize(2);
	Polar_Reco_MCRP_Full->SetMarkerStyle(20);
	Polar_Reco_MCEP_Full->SetLineColor(kBlack);
	Polar_Reco_MCEP_Full->SetMarkerColor(kBlack);
	Polar_Reco_MCEP_Full->SetMarkerSize(2);
	Polar_Reco_MCEP_Full->SetMarkerStyle(20);
	
	Polar_MC_Prim->SetLineColor(kRed+1);
	Polar_MC_Prim->SetMarkerColor(kRed+1);
	Polar_MC_Prim->SetMarkerSize(2);
	Polar_MC_Prim->SetMarkerStyle(26);	
	
	Polar_Reco_MCRP_Prim->SetLineColor(kBlack);
	Polar_Reco_MCRP_Prim->SetMarkerColor(kBlack);
	Polar_Reco_MCRP_Prim->SetMarkerSize(2);
	Polar_Reco_MCRP_Prim->SetMarkerStyle(22);
	Polar_Reco_MCEP_Prim->SetLineColor(kBlack);
	Polar_Reco_MCEP_Prim->SetMarkerColor(kBlack);
	Polar_Reco_MCEP_Prim->SetMarkerSize(2);
	Polar_Reco_MCEP_Prim->SetMarkerStyle(22);
	
	Polar_STAR->SetLineColor(kBlue+1);
	Polar_STAR->SetMarkerColor(kBlue+1);
	Polar_STAR->SetMarkerSize(2);
	Polar_STAR->SetMarkerStyle(34);
	Polar_STAR_nores->SetLineColor(kBlue+1);
	Polar_STAR_nores->SetMarkerColor(kBlue+1);
	Polar_STAR_nores->SetMarkerSize(2);
	Polar_STAR_nores->SetMarkerStyle(28);	
	
	// Plotting
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		
		if (angle_choice == "RP")
		{
			c1->cd(iter_cent+1);	
			hDeltaPhiRP_Full[iter_cent]->Draw("p9");
			Fit_hDeltaPhiRP_Full[iter_cent] = hDeltaPhiRP_Full[iter_cent]->GetFunction(Form("fitting_fnc_hDeltaPhiRP_Full_cent%d", iter_cent));
			Fit_hDeltaPhiRP_Full[iter_cent]->SetLineWidth(4);
			Fit_hDeltaPhiRP_Full[iter_cent]->SetLineColor(2);
			c1->Update();
			double c1_Ymax = gPad->GetUymax();
			double c1_Ymin = gPad->GetUymin();
			TLatex *title_MCRecpol_Full = new TLatex(3.5,c1_Ymin+0.8*(c1_Ymax-c1_Ymin), Form("P_{H} (MCRec) = %.4f #pm %.4f",Polar_Reco_MCRP_Full->GetPointY(iter_cent),TMath::Abs(Polar_Reco_MCRP_Full->GetErrorY(iter_cent)))); 
			title_MCRecpol_Full->Draw("same");
			TLatex *title_MCpol_Full = new TLatex(3.5,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), Form("P_{H} (MC) = %.4f #pm %.4f",Polar_MC_Full->GetPointY(iter_cent),TMath::Abs(Polar_MC_Full->GetErrorY(iter_cent)))); 
			title_MCpol_Full->Draw("same");
			TLatex *title_chi_Full = new TLatex(3.5,c1_Ymin+0.7*(c1_Ymax-c1_Ymin), Form("#chi^{2}/ndf = %.2f/%d",Fit_hDeltaPhiRP_Full[iter_cent]->GetChisquare(),Fit_hDeltaPhiRP_Full[iter_cent]->GetNDF())); 
			title_chi_Full->Draw("same");
			TLatex *title_NL_Full = new TLatex(3.5,c1_Ymin+0.6*(c1_Ymax-c1_Ymin), Form("N_{#Lambda} =  %1.2e",n_L_hDeltaPhiRP_Full[iter_cent])); 
			title_NL_Full->Draw("same");
			TLatex *title_cent_intervals_Full;
			if(NITER_CENT == 4)
			{
				title_cent_intervals_Full = new TLatex(0.4,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), cent_interval_4bins[iter_cent]);
			}
			else if(NITER_CENT == 7)
			{
				title_cent_intervals_Full = new TLatex(0.4,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), cent_interval_7bins[iter_cent]);			
			}
			else if(NITER_CENT == 10)
			{
				title_cent_intervals_Full = new TLatex(0.4,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), cent_interval_10bins[iter_cent]);
			}
			title_cent_intervals_Full->Draw("same");
			
			c2->cd(iter_cent+1);	
			hDeltaPhiRP_Prim[iter_cent]->Draw("p9");
			Fit_hDeltaPhiRP_Prim[iter_cent] = hDeltaPhiRP_Prim[iter_cent]->GetFunction(Form("fitting_fnc_hDeltaPhiRP_Prim_cent%d", iter_cent));
			Fit_hDeltaPhiRP_Prim[iter_cent]->SetLineWidth(4);
			Fit_hDeltaPhiRP_Prim[iter_cent]->SetLineColor(2);
			c2->Update();
			double c2_Ymax = gPad->GetUymax();
			double c2_Ymin = gPad->GetUymin();
			TLatex *title_MCRecpol_Prim = new TLatex(3.5,c2_Ymin+0.8*(c2_Ymax-c2_Ymin), Form("P_{H} (MCRec) = %.4f #pm %.4f",Polar_Reco_MCRP_Prim->GetPointY(iter_cent),TMath::Abs(Polar_Reco_MCRP_Prim->GetErrorY(iter_cent)))); 
			title_MCRecpol_Prim->Draw("same");
			TLatex *title_MCpol_Prim = new TLatex(3.5,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), Form("P_{H} (MC) = %.4f #pm %.4f",Polar_MC_Prim->GetPointY(iter_cent),TMath::Abs(Polar_MC_Prim->GetErrorY(iter_cent)))); 
			title_MCpol_Prim->Draw("same");
			TLatex *title_chi_Prim = new TLatex(3.5,c2_Ymin+0.7*(c2_Ymax-c2_Ymin), Form("#chi^{2}/ndf = %.2f/%d",Fit_hDeltaPhiRP_Prim[iter_cent]->GetChisquare(),Fit_hDeltaPhiRP_Prim[iter_cent]->GetNDF())); 
			title_chi_Prim->Draw("same");
			TLatex *title_NL_Prim = new TLatex(3.5,c2_Ymin+0.6*(c2_Ymax-c2_Ymin), Form("N_{#Lambda} =  %1.2e",n_L_hDeltaPhiRP_Prim[iter_cent])); 
			title_NL_Prim->Draw("same");
			TLatex *title_cent_intervals_Prim;
			if(NITER_CENT == 4)
			{
				title_cent_intervals_Prim = new TLatex(0.4,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), cent_interval_4bins[iter_cent]);
			}
			else if(NITER_CENT == 7)
			{
				title_cent_intervals_Prim = new TLatex(0.4,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), cent_interval_7bins[iter_cent]);			
			}
			else if(NITER_CENT == 10)
			{
				title_cent_intervals_Prim = new TLatex(0.4,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), cent_interval_10bins[iter_cent]);
			}
			title_cent_intervals_Prim->Draw("same");
		
		}
		else if (angle_choice == "EP")
		{
			c1->cd(iter_cent+1);	
			hDeltaPhiEP_Full[iter_cent]->Draw("p9");
			Fit_hDeltaPhiEP_Full[iter_cent] = hDeltaPhiEP_Full[iter_cent]->GetFunction(Form("fitting_fnc_hDeltaPhiEP_Full_cent%d", iter_cent));
			Fit_hDeltaPhiEP_Full[iter_cent]->SetLineWidth(4);
			Fit_hDeltaPhiEP_Full[iter_cent]->SetLineColor(2);
			c1->Update();
			double c1_Ymax = gPad->GetUymax();
			double c1_Ymin = gPad->GetUymin();
			TLatex *title_MCRecpol_Full = new TLatex(3.5,c1_Ymin+0.8*(c1_Ymax-c1_Ymin), Form("P_{H} (MCRec) = %.4f #pm %.4f",Polar_Reco_MCEP_Full->GetPointY(iter_cent),TMath::Abs(Polar_Reco_MCEP_Full->GetErrorY(iter_cent)))); 
			title_MCRecpol_Full->Draw("same");
			TLatex *title_MCpol_Full = new TLatex(3.5,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), Form("P_{H} (MC) = %.4f #pm %.4f",Polar_MC_Full->GetPointY(iter_cent),TMath::Abs(Polar_MC_Full->GetErrorY(iter_cent)))); 
			title_MCpol_Full->Draw("same");
			TLatex *title_chi = new TLatex(3.5,c1_Ymin+0.7*(c1_Ymax-c1_Ymin), Form("#chi^{2}/ndf = %.2f/%d",Fit_hDeltaPhiEP_Full[iter_cent]->GetChisquare(),Fit_hDeltaPhiEP_Full[iter_cent]->GetNDF())); 
			title_chi->Draw("same");
			TLatex *title_NL = new TLatex(3.5,c1_Ymin+0.6*(c1_Ymax-c1_Ymin), Form("N_{#Lambda} =  %1.2e",n_L_hDeltaPhiEP_Full[iter_cent])); 
			title_NL->Draw("same");
			TLatex *title_cent_intervals;
			if(NITER_CENT == 4)
			{
				title_cent_intervals = new TLatex(0.4,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), cent_interval_4bins[iter_cent]);
			}
			else if(NITER_CENT == 7)
			{
				title_cent_intervals = new TLatex(0.4,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), cent_interval_7bins[iter_cent]);			
			}
			else if(NITER_CENT == 10)
			{
				title_cent_intervals = new TLatex(0.4,c1_Ymin+0.9*(c1_Ymax-c1_Ymin), cent_interval_10bins[iter_cent]);
			}
			title_cent_intervals->Draw("same");
			
			c2->cd(iter_cent+1);	
			hDeltaPhiEP_Prim[iter_cent]->Draw("p9");
			Fit_hDeltaPhiEP_Prim[iter_cent] = hDeltaPhiEP_Prim[iter_cent]->GetFunction(Form("fitting_fnc_hDeltaPhiEP_Prim_cent%d", iter_cent));
			Fit_hDeltaPhiEP_Prim[iter_cent]->SetLineWidth(4);
			Fit_hDeltaPhiEP_Prim[iter_cent]->SetLineColor(2);
			c2->Update();
			double c2_Ymax = gPad->GetUymax();
			double c2_Ymin = gPad->GetUymin();
			TLatex *title_MCRecpol_Prim = new TLatex(3.5,c2_Ymin+0.8*(c2_Ymax-c2_Ymin), Form("P_{H} (MCRec) = %.4f #pm %.4f",Polar_Reco_MCRP_Prim->GetPointY(iter_cent),TMath::Abs(Polar_Reco_MCEP_Prim->GetErrorY(iter_cent)))); 
			title_MCRecpol_Prim->Draw("same");
			TLatex *title_MCpol_Prim = new TLatex(3.5,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), Form("P_{H} (MC) = %.4f #pm %.4f",Polar_MC_Prim->GetPointY(iter_cent),TMath::Abs(Polar_MC_Prim->GetErrorY(iter_cent)))); 
			title_MCpol_Prim->Draw("same");
			TLatex *title_chi_Prim = new TLatex(3.5,c2_Ymin+0.7*(c2_Ymax-c2_Ymin), Form("#chi^{2}/ndf = %.2f/%d",Fit_hDeltaPhiEP_Prim[iter_cent]->GetChisquare(),Fit_hDeltaPhiEP_Prim[iter_cent]->GetNDF())); 
			title_chi_Prim->Draw("same");
			TLatex *title_NL_Prim = new TLatex(3.5,c2_Ymin+0.6*(c2_Ymax-c2_Ymin), Form("N_{#Lambda} =  %1.2e",n_L_hDeltaPhiEP_Prim[iter_cent])); 
			title_NL_Prim->Draw("same");
			TLatex *title_cent_intervals_Prim;
			if(NITER_CENT == 4)
			{
				title_cent_intervals_Prim = new TLatex(0.4,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), cent_interval_4bins[iter_cent]);
			}
			else if(NITER_CENT == 7)
			{
				title_cent_intervals_Prim = new TLatex(0.4,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), cent_interval_7bins[iter_cent]);			
			}
			else if(NITER_CENT == 10)
			{
				title_cent_intervals_Prim = new TLatex(0.4,c2_Ymin+0.9*(c2_Ymax-c2_Ymin), cent_interval_10bins[iter_cent]);
			}
			title_cent_intervals_Prim->Draw("same");
		}
		else
		{
			cout << "This angle choice is not defined! Please choose either RP or EP" << endl;
			return;
		}
		c3->cd(iter_cent+1);	
		hPolarY_Full[iter_cent]->SetLineColor(kBlack);
		hPolarY_Full[iter_cent]->SetLineWidth(3);
		hPolarY_Full[iter_cent]->Draw();
		c3->Update();
		double c3_Ymax = gPad->GetUymax();
		double c3_Ymin = gPad->GetUymin();
		TLatex *title_cent_intervals_Full;
		if(NITER_CENT == 4)
		{
			title_cent_intervals_Full = new TLatex(-0.9,c3_Ymin+0.9*(c3_Ymax-c3_Ymin), cent_interval_4bins[iter_cent]);
		}
		else if(NITER_CENT == 7)
		{
			title_cent_intervals_Full = new TLatex(-0.9,c3_Ymin+0.9*(c3_Ymax-c3_Ymin), cent_interval_7bins[iter_cent]);			
		}
		else if(NITER_CENT == 10)
		{
			title_cent_intervals_Full = new TLatex(-0.9,c3_Ymin+0.9*(c3_Ymax-c3_Ymin), cent_interval_10bins[iter_cent]);
		}
		title_cent_intervals_Full->Draw("same");
		
		c4->cd(iter_cent+1);	
		hPolarY_Prim[iter_cent]->SetLineColor(kBlack);
		hPolarY_Prim[iter_cent]->SetLineWidth(3);
		hPolarY_Prim[iter_cent]->Draw();
		c4->Update();
		double c4_Ymax = gPad->GetUymax();
		double c4_Ymin = gPad->GetUymin();
		TLatex *title_cent_intervals_Prim;
		if(NITER_CENT == 4)
		{
			title_cent_intervals_Prim = new TLatex(-0.9,c4_Ymin+0.9*(c4_Ymax-c4_Ymin), cent_interval_4bins[iter_cent]);
		}
		else if(NITER_CENT == 7)
		{
			title_cent_intervals_Prim = new TLatex(-0.9,c4_Ymin+0.9*(c4_Ymax-c4_Ymin), cent_interval_7bins[iter_cent]);			
		}
		else if(NITER_CENT == 10)
		{
			title_cent_intervals_Prim = new TLatex(-0.9,c4_Ymin+0.9*(c4_Ymax-c4_Ymin), cent_interval_10bins[iter_cent]);
		}
		title_cent_intervals_Prim->Draw("same");
	}
	
	TLine *line = new TLine(0.,0.,80.,0.);
	line->SetLineColor(kBlack);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	TLatex *title_system = new TLatex(32.0,4.5, "Bi+Bi, #sqrt{s_{NN}} = 9.02 GeV"); 
	
	if (angle_choice == "RP")
	{
		c5->cd();	
		Polar_MC_Full->GetYaxis()->SetRangeUser(-1.0,5.);
		Polar_MC_Full->Draw("ap");
		Polar_Reco_MCRP_Full->Draw("psame");
		
		TLegend *legend5_1=new TLegend(0.15,0.65,0.35,0.9);
		legend5_1->SetTextFont(72);
		legend5_1->SetTextSize(0.04);
		legend5_1->SetBorderSize(0);
		legend5_1->AddEntry(Polar_MC_Full,"MC","p");
		legend5_1->AddEntry(Polar_Reco_MCRP_Full,"MC (fit)","p");		
		if (NITER_CENT == 4)
		{
			Polar_STAR->Draw("psame");
			Polar_STAR_nores->Draw("psame");
			legend5_1->AddEntry(Polar_STAR,"STAR (11.5 GeV)","p");	
			legend5_1->AddEntry(Polar_STAR_nores,"STAR (w/o res.)","p");	
		}
		legend5_1->Draw("same");
		line->Draw("same");	
		title_system->Draw("same");
		
		c6->cd();	
		Polar_MC_Full->GetYaxis()->SetRangeUser(-1.0,5.);
		Polar_MC_Full->Draw("ap");
		Polar_Reco_MCRP_Full->Draw("psame");
		Polar_MC_Prim->Draw("psame");
		Polar_Reco_MCRP_Prim->Draw("psame");
		
		TLegend *legend6_1=new TLegend(0.15,0.65,0.35,0.9);
		legend6_1->SetTextFont(72);
		legend6_1->SetTextSize(0.04);
		legend6_1->SetBorderSize(0);
		legend6_1->AddEntry(Polar_MC_Full,"MC","p");
		legend6_1->AddEntry(Polar_Reco_MCRP_Full,"MC (fit)","p");	
		legend6_1->AddEntry(Polar_MC_Prim,"MC (primary)","p");		
		legend6_1->AddEntry(Polar_Reco_MCRP_Prim,"MC (fit, primary)","p");	
		if (NITER_CENT == 4)
		{
			Polar_STAR->Draw("psame");
			Polar_STAR_nores->Draw("psame");
			legend6_1->AddEntry(Polar_STAR,"STAR (11.5 GeV)","p");	
			legend6_1->AddEntry(Polar_STAR_nores,"STAR (w/o res.)","p");	
		}
		legend6_1->Draw("same");
		line->Draw("same");	
		title_system->Draw("same");
	}
	else if(angle_choice == "EP")
	{
		c5->cd();	
		Polar_MC_Full->GetYaxis()->SetRangeUser(-1.0,5.);
		Polar_MC_Full->Draw("ap");
		Polar_Reco_MCEP_Full->Draw("psame");
		
		TLegend *legend5_1=new TLegend(0.15,0.65,0.35,0.9);
		legend5_1->SetTextFont(72);
		legend5_1->SetTextSize(0.04);
		legend5_1->SetBorderSize(0);
		legend5_1->AddEntry(Polar_MC_Full,"MC","p");
		legend5_1->AddEntry(Polar_Reco_MCEP_Full,"MC (fit)","p");		
		if (NITER_CENT == 4)
		{
			Polar_STAR->Draw("psame");
			Polar_STAR_nores->Draw("psame");
			legend5_1->AddEntry(Polar_STAR,"STAR (11.5 GeV)","p");	
			legend5_1->AddEntry(Polar_STAR_nores,"STAR (w/o res.)","p");	
		}
		legend5_1->Draw("same");
		line->Draw("same");	
		title_system->Draw("same");
		
		c6->cd();	
		Polar_MC_Full->GetYaxis()->SetRangeUser(-1.0,5.);
		Polar_MC_Full->Draw("ap");
		Polar_Reco_MCEP_Full->Draw("psame");
		Polar_MC_Prim->Draw("psame");
		Polar_Reco_MCEP_Prim->Draw("psame");
		
		TLegend *legend6_1=new TLegend(0.15,0.65,0.35,0.9);
		legend6_1->SetTextFont(72);
		legend6_1->SetTextSize(0.04);
		legend6_1->SetBorderSize(0);
		legend6_1->AddEntry(Polar_MC_Full,"MC","p");
		legend6_1->AddEntry(Polar_Reco_MCEP_Full,"MC (fit)","p");	
		legend6_1->AddEntry(Polar_MC_Prim,"MC (primary)","p");		
		legend6_1->AddEntry(Polar_Reco_MCEP_Prim,"MC (fit, primary)","p");	
		if (NITER_CENT == 4)
		{
			Polar_STAR->Draw("psame");
			Polar_STAR_nores->Draw("psame");
			legend6_1->AddEntry(Polar_STAR,"STAR (11.5 GeV)","p");	
			legend6_1->AddEntry(Polar_STAR_nores,"STAR (w/o res.)","p");	
		}
		legend6_1->Draw("same");
		line->Draw("same");	
		title_system->Draw("same");
	}
}
