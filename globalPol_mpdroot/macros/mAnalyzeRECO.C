/**
 * @file mAnalyzeRECO.C
 * @author Elizaveta Nazarova
 * @brief Analysis macro for RECO GlobalPol wagon
 * @version 0.1
 * @date 2023-07-15
 * All necessary calculations for the collected histograms using the GlobalPolarizationRECO wagon
 * Make sure that the input values correspond to the ones in the config file used to run the wagon!
 * Calculation of EP resolution in each analyzed centrality bin
 * Distributions of model polarization (P_{y}) in different centrality bins (for primary and full particles), obtained for associated MC tracks for true reconstructed Lambda
 * Fitting of the angular distributions (\Psi_{RP} - \phi or \Psi_{EP} - \phi) 
 * Graphs of average polarization: MC from the model, MC (fit), Reco using EP method (with only statistical uncertainties), values from STAR experiment obtained at energy of 11.5GeV with/without resolution correction
 * Fits of the subtracted background for systematic checks
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
#define pi TMath::Pi()
double ResEventPlane(double chi);
double Chi(double res);
void SidebandsFitInvMass(TH1D *hm0, TString name_c1, TString fullpath, double xmin, double nsig_bckg, double nsig, double &net_lambda_sig_fit, double &net_lambda_full_value_fit, double &net_lambda_bckg_fit);

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
// Function for the fit of background and signal
double fitting_function(double *x, double *par) {
	double vxmix = 1.07;
	double vxmax = 1.17;
	double e = 2.*(x[0] - vxmax)/(vxmax - vxmix) - 1.;
	return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*e + par[5]*0.5*(3.*TMath::Power(e,2) - 1.) + par[6]*0.5*e*(5.*TMath::Power(e,2) - 3.) + par[7]*0.125*(35.*TMath::Power(e,4) - 30.*TMath::Power(e,2) + 3.));
}
//fit function for the angle difference (delta(phi))
double fitting_fnc_delptaPhi(double *x, double *par) 
{
   return par[0]*(1. + 2.*par[1]*TMath::Sin(x[0]) + 2.*par[2]*TMath::Cos(x[0]) + 2.*par[3]*TMath::Sin(2.*x[0]) + 2.*par[4]*TMath::Cos(2.*x[0]));
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
/**
 * @brief Main function for analyzed results of the GlobalPolMC wagon
 * 
 * @param inFile 			input file
 * @param outFile 			output file
 * @param fullpath 			absolute path for the directory to save plots
 * @param NITER_CENT		number of centrality bins 
 * @param NITER_PT		    number of pt bins 
 * @param NITER_ETA		    number of eta bins 
 * @param NITER		        number of deltaphi bins 
 * @param xmin              minimal value of invariant mass (X-axis) for fitting (default: 1.086)
 * @param nsig_bckg         value to obtain the range for sidebands for background fit around the peak (+/- nsig_bckg*sigma) (default: 7.0)
 * @param nsig              value to obtain the range for signal extraction around the peak (+/- nsig*sigma) (default: 3.0)
 * @param particle_pdg 		pdg of analyzed hyperon
 */
void mAnalyzeRECO(TString inFile = "Output_GlopalPolRECO.root", TString outFile = "Output_GlopalPolRECO_Final.root", TString fullpath = "path-for-saved-plots", const int NITER_CENT = 4, const int NITER_PT = 5, const int NITER_ETA = 6, const int NITER = 20, double xmin = 1.086, double nsig_bckg = 7.0, const double nsig = 4.0, double xmin = 1.086, int particle_pdg = 3122)
{
	// histograms and variables
	TH1D *hPolarY_Full[NITER_CENT], *hDeltaPhiRP_Full[NITER_CENT];
	TF1 *fnc_hDeltaPhiRP_Full[NITER_CENT];
	
	TH1D *hm0[NITER_CENT][NITER];
	TH1D *hm0_ptbin[NITER_PT][NITER];
	TH1D *hm0_etabin[NITER_ETA][NITER];
	
	TH1D *Deltaphi_inv_mass[NITER_CENT], *Deltaphi_inv_mass_bckg[NITER_CENT], *Deltaphi_inv_mass_pt[NITER_PT], *Deltaphi_inv_mass_eta[NITER_ETA];
	TF1 *fitting_fnc_deltaphi_reco[NITER_CENT], *fitting_fnc_deltaphi_reco_bckg[NITER_CENT], *fitting_fnc_deltaphi_reco_pt[NITER_PT], *fitting_fnc_deltaphi_reco_eta[NITER_ETA];
	
	TProfile *hPolvsPt, *hPolvsEta;
	
	double mean_pol_MC_Full[NITER_CENT], mean_pol_MC_Full_err[NITER_CENT], mean_pol_recoRP_Full_MC[NITER_CENT], mean_pol_recoRP_Full_MC_err[NITER_CENT];
	double mean_pol_RECO_Full[NITER_CENT], mean_pol_RECO_Full_err[NITER_CENT], mean_pol_RECO_Full_bckg[NITER_CENT], mean_pol_RECO_Full_bckg_err[NITER_CENT];
	double mean_pol_RECO_Full_pt[NITER_CENT], mean_pol_RECO_Full_pt_err[NITER_CENT], mean_pol_RECO_Full_eta[NITER_CENT], mean_pol_RECO_Full_eta_err[NITER_CENT];
	
	double binmin_hDeltaPhiRP_Full[NITER_CENT];
	double net_lambda_cent[NITER_CENT][NITER];
	double net_lambda_full_cent[NITER_CENT][NITER];
	double net_lambda_bckg_cent[NITER_CENT][NITER];
	double net_lambda_pt[NITER_PT][NITER];
	double net_lambda_full_pt[NITER_PT][NITER];
	double net_lambda_bckg_pt[NITER_PT][NITER];
	double net_lambda_eta[NITER_ETA][NITER];
	double net_lambda_full_eta[NITER_ETA][NITER];
	double net_lambda_bckg_eta[NITER_ETA][NITER];
		
	// input file
	TFile *inFile_data = new TFile(inFile);

	// calculate the EP resolution
	TH1D *hNevCentr = (TH1D*) inFile_data->Get("hNevCentr");
	TH1D *hResolution_EP1_true, *hResolution_EP1_reco;
	double NEv_cent[NITER_CENT], ResEP1_true[NITER_CENT], ResEP1_exp[NITER_CENT], SubEvRes1[NITER_CENT];
	
	hResolution_EP1_true = (TH1D*) inFile_data->Get("hResolution_EP1_true");
	hResolution_EP1_reco = (TH1D*) inFile_data->Get("hResolution_EP1_reco");
		
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		NEv_cent[iter_cent] = hNevCentr->GetBinContent(iter_cent+1);
		ResEP1_true[iter_cent] = TMath::Abs(hResolution_EP1_true->GetBinContent(iter_cent+1)/NEv_cent[iter_cent]);
		SubEvRes1[iter_cent] = TMath::Abs(hResolution_EP1_reco->GetBinContent(iter_cent+1)/NEv_cent[iter_cent]);
		SubEvRes1[iter_cent] = TMath::Sqrt(SubEvRes1[iter_cent]);
		ResEP1_exp[iter_cent] = ResEventPlane(TMath::Sqrt(2.)*Chi(SubEvRes1[iter_cent]));
		cout << "iter_cent = " << iter_cent << "; ResEP1_true = " << ResEP1_true[iter_cent] << "; ResEP1_exp = " << ResEP1_exp[iter_cent] << endl;
	}
			
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		hResolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);	
		hResolution_EP1_reco->SetBinContent(iter_cent+1,ResEP1_exp[iter_cent]);	
	}

	// Get the histograms from the input file, calculate model global polarization from P_{y} distributions		
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		hPolarY_Full[iter_cent] = (TH1D*) inFile_data->Get(Form("hPolarY_Full_cent%d",iter_cent));	
		hDeltaPhiRP_Full[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiRP_Full_cent%d",iter_cent));			
		binmin_hDeltaPhiRP_Full[iter_cent] = hDeltaPhiRP_Full[iter_cent]->GetXaxis()->GetXmin();
		mean_pol_MC_Full[iter_cent] = -100.*hPolarY_Full[iter_cent]->GetMean();
		mean_pol_MC_Full_err[iter_cent] = -100.*hPolarY_Full[iter_cent]->GetMeanError();
		for(int iter = 0; iter < NITER; iter++)
		{
			hm0[iter_cent][iter] = (TH1D*) myFile_data->Get(Form("hm0_%d_%d",iter_cent,iter));
		}
		Deltaphi_inv_mass[iter_cent] = new TH1D(Form("Deltaphi_inv_mass_%d",iter_cent),Form("Deltaphi_inv_mass_%d",iter_cent),NITER,0.,2.*pi);
		Deltaphi_inv_mass_bckg[iter_cent] = new TH1D(Form("Deltaphi_inv_mass_bckg_%d",iter_cent),Form("Deltaphi_inv_mass_bckg_%d",iter_cent),NITER,0.,2.*pi);
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter = 0; iter < NITER; iter++)
		{
			double net_lambda_sig_fit = 0.;
			double net_lambda_full_value_fit = 0.;
			double net_lambda_bckg_fit = 0.;
			TString name_canvas = Form("Fit_invmass_cent%d_iter%d.png",iter_cent,iter);
			SidebandsFitInvMass(hm0[iter_cent][iter], name_canvas, fullpath, xmin, nsig_bckg, nsig, net_lambda_sig_fit, net_lambda_full_value_fit, net_lambda_bckg_fit);
			net_lambda_cent[iter_cent][iter] = net_lambda_sig_fit;
			net_lambda_full_cent[iter_cent][iter] = net_lambda_full_value_fit;
			net_lambda_bckg_cent[iter_cent][iter] = net_lambda_bckg_fit;
		}
	}
	
	// set content and errors:
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		int i_cut1 = hDeltaPhiRP_Full[iter_cent]->FindBin(binmin_hDeltaPhiRP_Full[iter_cent]);
		for(int iter_bin = 0; iter_bin < NITER_CENT; iter_bin++)
		{
			hDeltaPhiRP_Full[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(hDeltaPhiRP_Full[iter_cent]->GetBinContent(i_cut1)));
			Deltaphi_inv_mass[iter_cent]->SetBinContent(i_cut,net_lambda_cent[iter_cent][iter]);
			Deltaphi_inv_mass[iter_cent]->SetBinError(i_cut,TMath::Sqrt(net_lambda_full_cent[iter_cent][iter]));
			Deltaphi_bckg_inv_mass[iter_cent]->SetBinContent(i_cut,net_lambda_bckg_cent[iter_cent][iter]);
			Deltaphi_bckg_inv_mass[iter_cent]->SetBinError(i_cut,TMath::Sqrt(net_lambda_bckg_cent[iter_cent][iter]));
			i_cut1++;
		}
	}	
	
	if(NITER_CENT == 4)
	{
		hPolvsPt = (TProfile*) myFile_data->Get("hPolvsPt_2");
		hPolvsEta = (TProfile*) myFile_data->Get("hPolvsEta_2");
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			for(int iter = 0; iter < NITER; iter++)
			{
				hm0_ptbin[iter_pt][iter] = (TH1D*) myFile_data->Get(Form("hm0_ptbin_%d_%d",iter_pt,iter));
				Deltaphi_inv_mass_pt[iter_pt] = new TH1D(Form("Deltaphi_inv_mass_pt_%d",iter_pt),Form("Deltaphi_inv_mass_pt_%d",iter_pt),NITER,0.,2.*pi);
			}
		}
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			for(int iter = 0; iter < NITER; iter++)
			{
				hm0_etabin[iter_eta][iter] = (TH1D*) myFile_data->Get(Form("hm0_etabin_%d_%d",iter_eta,iter));
			}
			Deltaphi_inv_mass_eta[iter_eta] = new TH1D(Form("Deltaphi_inv_mass_eta_%d",iter_eta),Form("Deltaphi_inv_mass_eta_%d",iter_eta),NITER,0.,2.*pi);
		}
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			for(int iter = 0; iter < NITER; iter++)
			{
				double net_lambda_sig_fit = 0.;
				double net_lambda_full_value_fit = 0.;
				double net_lambda_bckg_fit = 0.;
				TString name_canvas = Form("Fit_invmass_pt%d_iter%d.png",iter_pt,iter);
				SidebandsFitInvMass(hm0_ptbin[iter_pt][iter], name_canvas, fullpath, xmin, nsig_bckg, nsig, net_lambda_sig_fit, net_lambda_full_value_fit, net_lambda_bckg_fit);
				net_lambda_pt[iter_pt][iter] = net_lambda_sig_fit;
				net_lambda_full_pt[iter_pt][iter] = net_lambda_full_value_fit;
				net_lambda_bckg_pt[iter_pt][iter] = net_lambda_bckg_fit;
			}
		}
		
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			for(int iter = 0; iter < NITER; iter++)
			{
				double net_lambda_sig_fit = 0.;
				double net_lambda_full_value_fit = 0.;
				double net_lambda_bckg_fit = 0.;
				TString name_canvas = Form("Fit_invmass_eta%d_iter%d.png",iter_eta,iter);
				SidebandsFitInvMass(hm0_etabin[iter_eta][iter], name_canvas, fullpath, xmin, nsig_bckg, nsig, net_lambda_sig_fit, net_lambda_full_value_fit, net_lambda_bckg_fit);
				net_lambda_eta[iter_eta][iter] = net_lambda_sig_fit;
				net_lambda_full_eta[iter_eta][iter] = net_lambda_full_value_fit;
				net_lambda_bckg_eta[iter_eta][iter] = net_lambda_bckg_fit;
			}
		}
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			int i_cut = hDeltaPhiRP_Full[0]->FindBin(hDeltaPhiRP_Full[0]->GetXaxis()->GetXmin());
			for(int iter = 0; iter < NITER; iter++)
			{
				Deltaphi_inv_mass_pt[iter_pt]->SetBinContent(i_cut,net_lambda_pt[iter_pt][iter]);
				Deltaphi_inv_mass_pt[iter_pt]->SetBinError(i_cut,TMath::Sqrt(net_lambda_full_pt[iter_pt][iter]));		
				i_cut++;
			}
		}
		
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			int i_cut = hDeltaPhiRP_Full[0]->FindBin(hDeltaPhiRP_Full[0]->GetXaxis()->GetXmin());
			for(int iter = 0; iter < NITER; iter++)
			{
				Deltaphi_inv_mass_eta[iter_eta]->SetBinContent(i_cut,net_lambda_eta[iter_eta][iter]);
				Deltaphi_inv_mass_eta[iter_eta]->SetBinError(i_cut,TMath::Sqrt(net_lambda_full_eta[iter_eta][iter]));		
				i_cut++;
			}
		}
	}
	
	
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		
		xmin_hDeltaPhiRP_Full[iter_cent] = hDeltaPhiRP_Full[iter_cent]->GetXaxis()->GetXmin();
		xmax_hDeltaPhiRP_Full[iter_cent] = hDeltaPhiRP_Full[iter_cent]->GetXaxis()->GetXmax();
		xmin_deltaphi_reco[iter_cent] = Deltaphi_inv_mass[iter_cent]->GetXaxis()->GetXmin();
		xmax_deltaphi_reco[iter_cent] = Deltaphi_inv_mass[iter_cent]->GetXaxis()->GetXmax();
		
		// create fitting functions for each histogram (delta(phi)_{RP} distributions)		
		fnc_hDeltaPhiRP_Full[iter_cent] = new TF1(Form("fitting_fnc_hDeltaPhiRP_Full_cent%d",iter_cent),fitting_fnc_delptaPhi,xmin_hDeltaPhiRP_Full[iter_cent],xmax_hDeltaPhiRP_Full[iter_cent],5);
		fitting_fnc_deltaphi_reco[iter_cent] = new TF1(Form("fitting_fnc_Deltaphi_inv_mass_cent%d",iter_cent),fitting_fnc_delptaPhi,xmin_deltaphi_reco[iter_cent],xmax_deltaphi_reco[iter_cent],5);
		fitting_fnc_deltaphi_reco_bckg[iter_cent] = new TF1(Form("fitting_fnc_Deltaphi_inv_mass_bckg_cent%d",iter_cent),fitting_fnc_delptaPhi,xmin_deltaphi_reco[iter_cent],xmax_deltaphi_reco[iter_cent],5);
		
		fnc_hDeltaPhiRP_Full[iter_cent]->SetParameters(hDeltaPhiRP_Full[iter_cent]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);	
		fitting_fnc_deltaphi_reco[iter_cent]->SetParameters(Deltaphi_inv_mass[iter_cent]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);
		fitting_fnc_deltaphi_reco_bckg[iter_cent]->SetParameters(Deltaphi_inv_mass_bckg[iter_cent]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);	
		
		// fit the histograms and calculate the global polarization from the fit parameter
		hDeltaPhiRP_Full[iter_cent]->Fit(Form("fitting_fnc_hDeltaPhiRP_Full_cent%d",iter_cent),"q","");
		mean_pol_recoRP_Full_MC[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiRP_Full[iter_cent]->GetParameter(1)));
		mean_pol_recoRP_Full_MC_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiRP_Full[iter_cent]->GetParError(1)));
		
		hDeltaPhiRP_Full[iter_cent]->Fit(Form("fitting_fnc_Deltaphi_inv_mass_cent%d",iter_cent),"q","");
		mean_pol_RECO_Full[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco[iter_cent]->GetParameter(1)));
		mean_pol_RECO_Full_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco[iter_cent]->GetParError(1)));
		
		hDeltaPhiRP_Full_bckg[iter_cent]->Fit(Form("fitting_fnc_Deltaphi_inv_mass_bckg_cent%d",iter_cent),"q","");
		mean_pol_RECO_Full_bckg[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco_bckg[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		mean_pol_RECO_Full_bckg_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco_bckg[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
	}
	
	// Initializing centrality bins, dependent on the number of bins
	int *centrality_min;
	int *centrality_max;
	double *noErr;
	double centrality_bin[NITER_CENT], centrality_bin_shiftleft[NITER_CENT]; 
	if (NITER_CENT == 4)
	{
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		noErr = init_double_array(4, 0, 0., 0., 0., 0.);		
	}
	else if (NITER_CENT == 7)
	{
		centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
		noErr = init_double_array(7, 0, 0., 0., 0., 0., 0., 0., 0.);		
	}
	else if (NITER_CENT == 10)
	{
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		noErr = init_double_array(10, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
	}
	else
	{
		cout << "This centrality binning is not defined!" << endl;
		return;
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		centrality_bin[iter_cent] = centrality_min[iter_cent] + (centrality_max[iter_cent] - centrality_min[iter_cent])/2;
		centrality_bin_shiftleft[iter_cent] = centrality_bin[iter_cent] - 1.5;
	}
	double star_cent_bin[] = {centrality_bin[2]+1.5};
	double star_cent_bin_err[] = {0.};

	// Experimental values from STAR collaboration for Lambda/ALambda at 11.5 GeV
	double star_pol_par[1], star_pol_err[1], star_pol_par_nores[1], star_pol_err_nores[1];
	if (particle_pdg == 3122) // For Lambda:
	{
		star_pol_par[0] = 1.179; // corrected value for 11.5GeV
		star_pol_err[0] = 0.347; // corrected error
		
		star_pol_par_nores[0] = 0.693; // corrected value for 11.5GeV (nores correction)
		star_pol_err_nores[0] = 0.202; // corrected error (nores correction)
	}
	else if (particle_pdg == -3122) // For ALambda:
	{
		star_pol_par[0] = 1.580; // corrected value for 11.5GeV
		star_pol_err[0] = 1.106; // corrected error
		
		star_pol_par_nores[0] = 0.921; // corrected value for 11.5GeV (nores correction)
		star_pol_err_nores[0] = 0.649; // corrected error (nores correction)
	}
	else
	{
		cout << "This particle choice is not defined! Please provide the definition in the code." << endl;
		return;
	}

	// Create the graphs of average global polarization in centrality bins					
	TGraphErrors *Polar_STAR = new TGraphErrors(1, star_cent_bin, star_pol_par, star_cent_bin_err, star_pol_err);	
	Polar_STAR->SetName("Polar_STAR");
	Polar_STAR->SetTitle("Polar_STAR");
	TGraphErrors *Polar_STAR_nores = new TGraphErrors(1, star_cent_bin, star_pol_par_nores, star_cent_bin_err, star_pol_err_nores);
	Polar_STAR_nores->SetName("Polar_STAR_nores");
	Polar_STAR_nores->SetTitle("Polar_STAR_nores");

	TGraphErrors *Polar_MC_Full = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_MC_Full, noErr, mean_pol_MC_Full_err);
	Polar_MC_Full->SetName("Polar_MC_Full");
	Polar_MC_Full->SetTitle("Polar_MC_Full");
	
	TGraphErrors *Polar_Reco_MCRP_Full = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_recoRP_Full_MC, noErr, mean_pol_recoRP_Full_MC_err);
	Polar_Reco_MCRP_Full->SetName("Polar_Reco_MCRP_Full");
	Polar_Reco_MCRP_Full->SetTitle("Polar_Reco_MCRP_Full");
	
	TGraphErrors *Polar_Reco_Cent = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_RECO_Full, noErr, mean_pol_RECO_Full_err);
	Polar_Reco_MCRP_Prim->SetName("Polar_Reco_MCRP_Prim");
	Polar_Reco_MCRP_Prim->SetTitle("Polar_Reco_MCRP_Prim");
	
	TGraphErrors *Polar_Reco_Cent_bckg = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_RECO_Full_bckg, noErr, mean_pol_RECO_Full_bckg_err);
	Polar_Reco_MCEP_Full->SetName("Polar_Reco_MCEP_Full");
	Polar_Reco_MCEP_Full->SetTitle("Polar_Reco_MCEP_Full");

	// Check the amount of Lambda in each case
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		double n_L_hDeltaPhiRP_Full = 0.;
		double n_L_hDeltaPhiReco = 0.;
		for(int iter_bin = 0; iter_bin < hDeltaPhiRP_Full[iter_cent]->GetNbinsX(); iter_bin++)
		{
			n_L_hDeltaPhiRP_Full += hDeltaPhiRP_Full[iter_cent]->GetBinContent(iter_bin+1);
			n_L_hDeltaPhiReco += Deltaphi_inv_mass[iter_cent]->GetBinContent(iter_bin+1);
		}
		cout << "N_L (hDeltaPhiRP_Full) = " << n_L_hDeltaPhiRP_Full << endl;
		cout << "N_L (Deltaphi_inv_mass) = " << n_L_hDeltaPhiReco << endl;
	}
	
	//Calculating polarization in pt and eta bins
	double mean_pol_reco_eta[NITER_ETA], mean_pol_reco_eta_err[NITER_ETA];
	double xmin_deltaphi_reco_eta[NITER_ETA], xmax_deltaphi_reco_eta[NITER_ETA];
	
	if(NITER_CENT == 4)
	{
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			xmin_deltaphi_reco_pt[iter_pt] = Deltaphi_inv_mass_pt[iter_pt]->GetXaxis()->GetXmin();
			xmax_deltaphi_reco_pt[iter_pt] = Deltaphi_inv_mass_pt[iter_pt]->GetXaxis()->GetXmax();
			fitting_fnc_deltaphi_reco_pt[iter_pt] = new TF1(Form("fitting_fnc_Deltaphi_inv_mass_pt%d",iter_pt),fitting_fnc_delptaPhi,xmin_deltaphi_reco_pt[iter_pt],xmax_deltaphi_reco_pt[iter_pt],5);
			fitting_fnc_deltaphi_reco_pt[iter_pt]->SetParameters(Deltaphi_inv_mass_pt[iter_pt]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);	
			Deltaphi_inv_mass_pt[iter_pt]->Fit(Form("fitting_fnc_Deltaphi_inv_mass_pt%d",iter_pt),"q","");
			mean_pol_RECO_Full_pt[iter_pt] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco_pt[iter_pt]->GetParameter(1)));
			mean_pol_RECO_Full_pt_err[iter_pt] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco_pt[iter_pt]->GetParError(1)));
		}
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			xmin_deltaphi_reco_eta[iter_eta] = Deltaphi_inv_mass_eta[iter_eta]->GetXaxis()->GetXmin();
			xmax_deltaphi_reco_eta[iter_eta] = Deltaphi_inv_mass_eta[iter_eta]->GetXaxis()->GetXmax();
			fitting_fnc_deltaphi_reco_eta[iter_eta] = new TF1(Form("fitting_fnc_Deltaphi_inv_mass_eta%d",iter_eta),fitting_fnc_delptaPhi,xmin_deltaphi_reco_eta[iter_eta],xmax_deltaphi_reco_eta[iter_eta],5);
			fitting_fnc_deltaphi_reco_eta[iter_eta]->SetParameters(Deltaphi_inv_mass_eta[iter_eta]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);	
			Deltaphi_inv_mass_eta[iter_eta]->Fit(Form("fitting_fnc_Deltaphi_inv_mass_eta%d",iter_eta),"q","");
			mean_pol_RECO_Full_eta[iter_eta] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco_eta[iter_eta]->GetParameter(1)));
			mean_pol_RECO_Full_eta_err[iter_eta] = 100.*(8./(pi*0.732))*(TMath::Abs(fitting_fnc_deltaphi_reco_eta[iter_eta]->GetParError(1)));
		}
		double n_L_eta[NITER_ETA];
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			n_L_eta[iter_eta] = 0.;
		}
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			for(int iter_bin = 0; iter_bin < NITER; iter_bin++)
			{
				n_L_eta[iter_eta] += Deltaphi_inv_mass_eta[iter_eta]->GetBinContent(iter_bin+1);
			}
			cout << "n_L_eta = " << n_L_eta[iter_eta] << endl;
		}
		double n_L_pt[NITER_PT];
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			n_L_pt[iter_pt] = 0.;
		}
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			for(int iter_bin = 0; iter_bin < NITER; iter_bin++)
			{
				n_L_pt[iter_pt] += Deltaphi_inv_mass_pt[iter_pt]->GetBinContent(iter_bin+1);
			}
			cout << "n_L_pt = " << n_L_pt[iter_pt] << endl;
		}
		double pt_bins[NITER_PT], eta_bins[NITER_ETA];
		double ptbin_min   = 0.0;
		double ptbin_step  = 0.5;       
		double etabin_min  = -1.5; 
		double etabin_step = 0.5;      
		
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			pt_bins[iter_pt] = ptbin_min + ptbin_step*(iter_pt + 0.5);
		}
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			eta_bins[iter_eta] = etabin_min + etabin_step*(iter_eta + 0.5);
		}
		TGraphErrors *Polar_cent_RECO_pt = new TGraphErrors(NITER_PT, pt_bins, mean_pol_RECO_Full_pt, noErr, mean_pol_RECO_Full_pt_err);
		Polar_cent_RECO_pt->SetName("Polar_RECO_pt");
		Polar_cent_RECO_pt->SetTitle("Polar_RECO_pt");
		TGraphErrors *Polar_cent_RECO_eta = new TGraphErrors(NITER_ETA, eta_bins, mean_pol_RECO_Full_eta, noErr, mean_pol_RECO_Full_eta_err);
		Polar_cent_RECO_eta->SetName("Polar_RECO_eta");
		Polar_cent_RECO_eta->SetTitle("Polar_RECO_eta");
		
		hPolvsPt->Scale(-100.);
		hPolvsEta->Scale(-100.);
	}
	
	// Write the results in the output file
	TFile out(outFile,"recreate");
	hResolution_EP1_true->Write();
	hResolution_EP1_reco->Write();
	Polar_STAR->Write();
	Polar_STAR_nores->Write();
	Polar_MC_Full->Write();
	Polar_Reco_MCRP_Full->Write();
	Polar_Reco_Cent->Write();
	Polar_Reco_Cent_bckg->Write();
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		hDeltaPhiRP_Full[iter_cent]->Write();
		Deltaphi_inv_mass[iter_cent]->Write();
		Deltaphi_inv_mass_bckg[iter_cent]->Write();
		hPolarY_Full[iter_cent]->Write();
	}
	if(NITER_CENT == 4)
	{
		Polar_cent_RECO_pt->Write();
		Polar_cent_RECO_eta->Write();
		for(int iter_pt = 0; iter_pt < NITER_PT; iter_pt++)
		{
			Deltaphi_inv_mass_pt[iter_pt]->Write();
		}
		for(int iter_eta = 0; iter_eta < NITER_ETA; iter_eta++)
		{
			Deltaphi_inv_mass_eta[iter_eta]->Write();
		}
	}
	out.Close();
	
}

// plane resolution as function of chi
double ResEventPlane(double chi)
{
	double con = TMath::Sqrt(TMath::Pi()/2)/2 ;   // ~ 0.626657
	double arg = chi * chi / 4.;
	double res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

	return res ;
}

// chi from the event plane resolution, by primitive inversion
double Chi(double res)
{
	double chi   = 2.0;
	double delta = 1.0;
	for(int i = 0; i < 15; i++)// for(int i = 0; i < 50; i++) ---- check!!!
	{
		if(ResEventPlane(chi) < res) { chi = chi + delta ; }
		else                         { chi = chi - delta ; }
		delta = delta / 2.;
	}
	cout << "chi = " << chi << endl;
	return chi ;
}
void SidebandsFitInvMass(TH1D *hm0, TString name_c1, TString fullpath, double xmin, double nsig_bckg, double nsig, double &net_lambda_sig_fit, double &net_lambda_full_value_fit, double &net_lambda_bckg_fit)
{
	TH1D *hm0_signal, *hm0_bckg;
	TF1 *fitting_fnc, *backFcn, *signalFcn;
	
	fullpath += '/';
	fullpath += name_c1;
	
	double sum_sig;
	double sum_full;
	double sum_bckg;
	double entries_new;
	double entries_old;
	double efficiency;
	double ratio_SB;
	double ratio_SSB;
	
	double xmax;
	int bin_left;
	int bin_right;
	
	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
	
	hm0->SetMarkerSize(1);
	hm0->SetMarkerStyle(20);
	hm0->SetLineColor(kBlack);
	hm0->SetMarkerColor(kBlack);
	
	hm0_signal = (TH1D*)hm0->Clone("hm0sig");
	hm0_bckg = (TH1D*)hm0->Clone("hm0bckg");
	xmax = hm0->GetXaxis()->GetXmax();
	
	fitting_fnc = new TF1("int_fitting_fnc",fitting_function,xmin,xmax,8);
	c1->cd();
	
	fitting_fnc->SetNpx(500);
	fitting_fnc->SetLineWidth(4);
	fitting_fnc->SetLineColor(2); //red color for fitting line
	fitting_fnc->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
	int npar = fitting_fnc->GetNumberFreeParameters();
	for (int ip = 0; ip < npar; ++ip) fitting_fnc->SetParameter(ip,0);	
	
	fitting_fnc->SetParameter(0,hm0->GetMaximum());
	fitting_fnc->SetParameter(1,1.116);
	fitting_fnc->SetParameter(2,0.002);
	
	hm0->SetMinimum(1.0E-20);
	hm0->GetXaxis()->SetRangeUser(1.07,1.17);
	hm0->Draw();
	hm0->Fit("int_fitting_fnc","wq0","",xmin,xmax);
	
	double mass = fitting_fnc->GetParameter(npar-7);
	double sigma = TMath::Abs(fitting_fnc->GetParameter(npar-6));
			
	double xmin_cut = mass - nsig_bckg * sigma;
	double xmax_cut = mass + nsig_bckg * sigma;
			
	int imin = hm0->FindBin(xmin_cut);
	int imax = hm0->FindBin(xmax_cut);
	
	//improve the picture:
	backFcn = new TF1("int_backFcn",background,xmin,xmax,5);
	backFcn->SetLineColor(kRed);
	backFcn->SetLineWidth(4);
				
	double errs[200] = {0};
	int nbins = hm0->GetNbinsX();
	for (int ib = 1; ib <= nbins; ++ib) 
	{
		if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg->GetBinContent(ib));
	}
	hm0_bckg->SetError(errs);
	hm0_bckg->SetLineColor(kMagenta);
	hm0_bckg->Fit("int_backFcn","q","same",xmin,xmax);
	
	signalFcn = new TF1("int_signalFcn",gaussian,xmin,xmax,3);
	signalFcn->SetLineColor(kBlue);
	signalFcn->SetNpx(500);
				
	hm0_signal->Sumw2();	
	hm0_signal->Add(backFcn, -1);
	hm0_signal->SetLineColor(kBlack);
	hm0_signal->SetLineWidth(2);
	hm0_signal->SetMarkerStyle(2);
	hm0_signal->SetMarkerSize(1);
	signalFcn->SetParameters(hm0->GetMaximum(),mass,sigma);
	hm0_signal->Fit("int_signalFcn","q0","same",xmin,xmax);
	
	mass = signalFcn->GetParameter(1);
	sigma = signalFcn->GetParameter(2);
				
	double left = mass - nsig * sigma;
	double right = mass + nsig * sigma;
				
	bin_left = hm0_signal->FindBin(left);
	bin_right = hm0_signal->FindBin(right);

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
	TLine *line_backleft = new TLine(xmin_cut,0,xmin_cut,0.6*hm0->GetMaximum());
	line_backleft->SetLineColor(kBlue+1);
	line_backleft->SetLineWidth(3);
	line_backleft->SetLineStyle(9);
	line_backleft->Draw("same");
	TLine *line_backright = new TLine(xmax_cut,0,xmax_cut,0.6*hm0->GetMaximum());
	line_backright->SetLineColor(kBlue+1);
	line_backright->SetLineWidth(3);
	line_backright->SetLineStyle(9);
	line_backright->Draw("same");
				
	entries_new = hm0_signal->Integral(bin_left,bin_right);
	entries_old = hm0->Integral();
	
	sum_sig = 0.;
	sum_full = 0.;
	sum_bckg = 0.;
			
	for (int ib = bin_left; ib <= bin_right; ++ib) 
	{
		sum_sig += hm0_signal->GetBinContent(ib);
		sum_full += hm0->GetBinContent(ib);
	}
	
	for (int ib = imin; ib <= imax; ++ib) 
	{
		sum_bckg += hm0->GetBinContent(ib);
	}
	
	sum_bckg = sum_bckg - sum_sig;
			
	net_lambda_sig_fit = sum_sig;
	net_lambda_full_value_fit = sum_full;
	net_lambda_bckg_fit = sum_bckg;
	
	TLegend *legend=new TLegend(0.6,0.65,0.88,0.9);
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->SetBorderSize(0);
	legend->AddEntry(hm0,"Data","lp");
	legend->AddEntry(backFcn,"Bckg fit","l");	
	legend->AddEntry(line_left,"Cut-off (signal)","l");
	legend->AddEntry(line_backleft,"Cut-off (bckg)","l");
	legend->Draw("same");
	
	ratio_SB = sum_sig/(sum_full - sum_sig);
	ratio_SSB = sum_sig/TMath::Sqrt(sum_full);
	efficiency = 100.*entries_new/entries_old;	
		
	TLatex *title_SSB = new TLatex(1.075,0.92*hm0->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB)); 
	title_SSB->Draw("same");
	TLatex *title_SB = new TLatex(1.075,0.76*hm0->GetMaximum(), Form("S/B = %.2f",ratio_SB)); 
	title_SB->Draw("same");
	TLatex *title_efficiency = new TLatex(1.075,0.66*hm0->GetMaximum(), Form("eff. = %.2f [%%]",efficiency)); 
	title_efficiency->Draw("same");	
	if (c1) 
	{ 
		c1->SaveAs(fullpath);
		c1->Close(); 
		delete c1;
		c1 = 0; 
	}
}
