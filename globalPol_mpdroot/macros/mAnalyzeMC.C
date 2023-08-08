/**
 * @file mAnalyzeMC.C
 * @author Elizaveta Nazarova
 * @brief Analysis macro for MC GlobalPol wagon
 * @version 0.1
 * @date 2023-07-15
 * All necessary calculations for the MCTest using the GlobalPolarizationMC wagon
 * Make sure that the input values correspond to the ones in the config file used to run the wagon!
 * Calculation of EP resolution in each analyzed centrality bin
 * Distributions of model polarization (P_{y}) in different centrality bins (for primary and full particles)
 * Fitting of the angular distributions (\Psi_{RP} - \phi or \Psi_{EP} - \phi) for primary and full particles
 * Graphs of average polarization: MC from the model for primary and full particles, MC (fit) for primary and full particles, values from STAR experiment obtained at energy of 11.5GeV with/without resolution correction
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
 * @param NITER_CENT		number of centrality bins 
 * @param particle_pdg 		pdg of analyzed hyperon
 */
void mAnalyzeMC(TString inFile = "Output_GlopalPolMC.root", TString outFile = "Output_GlopalPolMC_Final.root", const int NITER_CENT = 4, int particle_pdg = 3122)
{
	// histograms and variables
	TH1D *hPolarY_Full[NITER_CENT], *hPolarY_Prim[NITER_CENT];
	TH1D *hDeltaPhiRP_Full[NITER_CENT], *hDeltaPhiRP_Prim[NITER_CENT], *hDeltaPhiEP_Full[NITER_CENT], *hDeltaPhiEP_Prim[NITER_CENT];
	
	TF1 *fnc_hDeltaPhiRP_Full[NITER_CENT], *fnc_hDeltaPhiRP_Prim[NITER_CENT], *fnc_hDeltaPhiEP_Full[NITER_CENT], *fnc_hDeltaPhiEP_Prim[NITER_CENT];
	
	double xmin_hDeltaPhiRP_Full[NITER_CENT], xmax_hDeltaPhiRP_Full[NITER_CENT];
	
	double mean_pol_MC_Full[NITER_CENT], mean_pol_MC_Full_err[NITER_CENT], mean_pol_recoRP_Full_MC[NITER_CENT], mean_pol_recoRP_Full_MC_err[NITER_CENT], mean_pol_recoEP_Full_MC[NITER_CENT], mean_pol_recoEP_Full_MC_err[NITER_CENT];
	double mean_pol_MC_Prim[NITER_CENT], mean_pol_MC_Prim_err[NITER_CENT], mean_pol_recoRP_Prim_MC[NITER_CENT], mean_pol_recoRP_Prim_MC_err[NITER_CENT], mean_pol_recoEP_Prim_MC[NITER_CENT], mean_pol_recoEP_Prim_MC_err[NITER_CENT];
	
	double binmin_hDeltaPhiRP_Full[NITER_CENT];
		
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
		hPolarY_Prim[iter_cent] = (TH1D*) inFile_data->Get(Form("hPolarY_Prim_cent%d",iter_cent));		
		
		hDeltaPhiRP_Full[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiRP_Full_cent%d",iter_cent));	
		hDeltaPhiRP_Prim[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiRP_Prim_cent%d",iter_cent));
	
		hDeltaPhiEP_Full[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiEP_Full_cent%d",iter_cent));	
		hDeltaPhiEP_Prim[iter_cent] = (TH1D*) inFile_data->Get(Form("hDeltaPhiEP_Prim_cent%d",iter_cent));
		
		binmin_hDeltaPhiRP_Full[iter_cent] = hDeltaPhiRP_Full[iter_cent]->GetXaxis()->GetXmin();
		mean_pol_MC_Full[iter_cent] = -100.*hPolarY_Full[iter_cent]->GetMean();
		mean_pol_MC_Full_err[iter_cent] = -100.*hPolarY_Full[iter_cent]->GetMeanError();
		mean_pol_MC_Prim[iter_cent] = -100.*hPolarY_Prim[iter_cent]->GetMean();
		mean_pol_MC_Prim_err[iter_cent] = -100.*hPolarY_Prim[iter_cent]->GetMeanError();
	}
	
	// set errors:
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		int i_cut1 = hDeltaPhiRP_Full[iter_cent]->FindBin(binmin_hDeltaPhiRP_Full[iter_cent]);
		for(int iter_bin = 0; iter_bin < NITER_CENT; iter_bin++)
		{
			hDeltaPhiRP_Full[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(hDeltaPhiRP_Full[iter_cent]->GetBinContent(i_cut1)));
			hDeltaPhiRP_Prim[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(hDeltaPhiRP_Prim[iter_cent]->GetBinContent(i_cut1)));
			hDeltaPhiEP_Full[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(hDeltaPhiEP_Full[iter_cent]->GetBinContent(i_cut1)));
			hDeltaPhiEP_Prim[iter_cent]->SetBinError(i_cut1,TMath::Sqrt(hDeltaPhiEP_Prim[iter_cent]->GetBinContent(i_cut1)));
			i_cut1++;
		}
	}	
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		
		xmin_hDeltaPhiRP_Full[iter_cent] = hDeltaPhiRP_Full[iter_cent]->GetXaxis()->GetXmin();
		xmax_hDeltaPhiRP_Full[iter_cent] = hDeltaPhiRP_Full[iter_cent]->GetXaxis()->GetXmax();
		
		// create fitting functions for each histogram (delta(phi)_{RP} distributions)		
		fnc_hDeltaPhiRP_Full[iter_cent] = new TF1(Form("fitting_fnc_hDeltaPhiRP_Full_cent%d",iter_cent),fitting_fnc_delptaPhi,xmin_hDeltaPhiRP_Full[iter_cent],xmax_hDeltaPhiRP_Full[iter_cent],5);
		fnc_hDeltaPhiRP_Prim[iter_cent] = new TF1(Form("fitting_fnc_hDeltaPhiRP_Prim_cent%d",iter_cent),fitting_fnc_delptaPhi,xmin_hDeltaPhiRP_Full[iter_cent],xmax_hDeltaPhiRP_Full[iter_cent],5);
		fnc_hDeltaPhiEP_Full[iter_cent] = new TF1(Form("fitting_fnc_hDeltaPhiEP_Full_cent%d",iter_cent),fitting_fnc_delptaPhi,xmin_hDeltaPhiRP_Full[iter_cent],xmax_hDeltaPhiRP_Full[iter_cent],5);
		fnc_hDeltaPhiEP_Prim[iter_cent] = new TF1(Form("fitting_fnc_hDeltaPhiEP_Prim_cent%d",iter_cent),fitting_fnc_delptaPhi,xmin_hDeltaPhiRP_Full[iter_cent],xmax_hDeltaPhiRP_Full[iter_cent],5);
		
		fnc_hDeltaPhiRP_Full[iter_cent]->SetParameters(hDeltaPhiRP_Full[iter_cent]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);	
		fnc_hDeltaPhiRP_Prim[iter_cent]->SetParameters(hDeltaPhiRP_Prim[iter_cent]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);
		fnc_hDeltaPhiEP_Full[iter_cent]->SetParameters(hDeltaPhiEP_Full[iter_cent]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);	
		fnc_hDeltaPhiEP_Prim[iter_cent]->SetParameters(hDeltaPhiEP_Prim[iter_cent]->GetMaximum(), 0.01, 0.0003, 0.0002, 0.002);	
		
		// fit the histograms and calculate the global polarization from the fit parameter
		hDeltaPhiRP_Full[iter_cent]->Fit(Form("fitting_fnc_hDeltaPhiRP_Full_cent%d",iter_cent),"q","");
		mean_pol_recoRP_Full_MC[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiRP_Full[iter_cent]->GetParameter(1)));
		mean_pol_recoRP_Full_MC_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiRP_Full[iter_cent]->GetParError(1)));
		
		hDeltaPhiRP_Prim[iter_cent]->Fit(Form("fitting_fnc_hDeltaPhiRP_Prim_cent%d",iter_cent),"q","");
		mean_pol_recoRP_Prim_MC[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiRP_Prim[iter_cent]->GetParameter(1)));
		mean_pol_recoRP_Prim_MC_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiRP_Prim[iter_cent]->GetParError(1)));
		
		hDeltaPhiEP_Full[iter_cent]->Fit(Form("fitting_fnc_hDeltaPhiEP_Full_cent%d",iter_cent),"q","");
		mean_pol_recoEP_Full_MC[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiEP_Full[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		mean_pol_recoEP_Full_MC_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiEP_Full[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
		
		hDeltaPhiEP_Prim[iter_cent]->Fit(Form("fitting_fnc_hDeltaPhiEP_Prim_cent%d",iter_cent),"q","");
		mean_pol_recoEP_Prim_MC[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiEP_Prim[iter_cent]->GetParameter(1)))/ResEP1_exp[iter_cent];
		mean_pol_recoEP_Prim_MC_err[iter_cent] = 100.*(8./(pi*0.732))*(TMath::Abs(fnc_hDeltaPhiEP_Prim[iter_cent]->GetParError(1)))/ResEP1_exp[iter_cent];
	}
	
	// Initializing centrality bins, dependent on the number of bins
	int *centrality_min;
	int *centrality_max;
	double *noErr;
	double centrality_bin[NITER_CENT]; 
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
	
	TGraphErrors *Polar_MC_Prim = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_MC_Prim, noErr, mean_pol_MC_Prim_err);
	Polar_MC_Prim->SetName("Polar_MC_Prim");
	Polar_MC_Prim->SetTitle("Polar_MC_Prim");
	
	TGraphErrors *Polar_Reco_MCRP_Full = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_recoRP_Full_MC, noErr, mean_pol_recoRP_Full_MC_err);
	Polar_Reco_MCRP_Full->SetName("Polar_Reco_MCRP_Full");
	Polar_Reco_MCRP_Full->SetTitle("Polar_Reco_MCRP_Full");
	
	TGraphErrors *Polar_Reco_MCRP_Prim = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_recoRP_Prim_MC, noErr, mean_pol_recoRP_Prim_MC_err);
	Polar_Reco_MCRP_Prim->SetName("Polar_Reco_MCRP_Prim");
	Polar_Reco_MCRP_Prim->SetTitle("Polar_Reco_MCRP_Prim");
	
	TGraphErrors *Polar_Reco_MCEP_Full = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_recoEP_Full_MC, noErr, mean_pol_recoEP_Full_MC_err);
	Polar_Reco_MCEP_Full->SetName("Polar_Reco_MCEP_Full");
	Polar_Reco_MCEP_Full->SetTitle("Polar_Reco_MCEP_Full");
	
	TGraphErrors *Polar_Reco_MCEP_Prim = new TGraphErrors(NITER_CENT, centrality_bin, mean_pol_recoEP_Prim_MC, noErr, mean_pol_recoEP_Prim_MC_err);
	Polar_Reco_MCEP_Prim->SetName("Polar_Reco_MCEP_Prim");
	Polar_Reco_MCEP_Prim->SetTitle("Polar_Reco_MCEP_Prim");
	
	// Check the amount of Lambda in each case
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		double n_L_hDeltaPhiRP_Full = 0.;
		double n_L_hDeltaPhiRP_Prim = 0.;
		double n_L_hDeltaPhiEP_Full = 0.;
		double n_L_hDeltaPhiEP_Prim = 0.;
		for(int iter_bin = 0; iter_bin < hDeltaPhiRP_Full[iter_cent]->GetNbinsX(); iter_bin++)
		{
			n_L_hDeltaPhiRP_Full += hDeltaPhiRP_Full[iter_cent]->GetBinContent(iter_bin+1);
			n_L_hDeltaPhiRP_Prim += hDeltaPhiRP_Prim[iter_cent]->GetBinContent(iter_bin+1);
			n_L_hDeltaPhiEP_Full += hDeltaPhiEP_Full[iter_cent]->GetBinContent(iter_bin+1);
			n_L_hDeltaPhiEP_Prim += hDeltaPhiEP_Prim[iter_cent]->GetBinContent(iter_bin+1);
		}
		cout << "N_L (hDeltaPhiRP_Full) = " << n_L_hDeltaPhiRP_Full << endl;
		cout << "N_L (hDeltaPhiRP_Prim) = " << n_L_hDeltaPhiRP_Prim << endl;
		cout << "N_L (hDeltaPhiEP_Full) = " << n_L_hDeltaPhiEP_Full << endl;
		cout << "N_L (hDeltaPhiEP_Prim) = " << n_L_hDeltaPhiEP_Prim << endl;
	}
	
	// Write the results in the output file
	TFile out(outFile,"recreate");
	hResolution_EP1_true->Write();
	hResolution_EP1_reco->Write();
	Polar_STAR->Write();
	Polar_STAR_nores->Write();
	Polar_MC_Full->Write();
	Polar_MC_Prim->Write();
	Polar_Reco_MCRP_Full->Write();
	Polar_Reco_MCRP_Prim->Write();
	Polar_Reco_MCEP_Full->Write();
	Polar_Reco_MCEP_Prim->Write();
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		hDeltaPhiRP_Full[iter_cent]->Write();
		hDeltaPhiRP_Prim[iter_cent]->Write();
		hDeltaPhiEP_Full[iter_cent]->Write();
		hDeltaPhiEP_Prim[iter_cent]->Write();
		hPolarY_Full[iter_cent]->Write();
		hPolarY_Prim[iter_cent]->Write();
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
