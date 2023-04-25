//creates distribution of invariant mass and performs the fitting procedure (Legendre polinoms + Gaus for Signal), then separate background fit in the sidebands, then cutting off the background, for the full phase space, with different selection criteria omega_2
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

void LegendreFit_selections(TString inFile = "PHSD_omegatest.root", Int_t NITER_CENT = 10, Int_t NITER = 10, Double_t xmin = 1.086, Double_t nsig_bckg = 7.0, Double_t nsig = 3.0)
{
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	
	//starting values for variables for selection: - check!!!
	double chi_pi_start = 0.; // starting value for chi_pi
	double chi_p_start = 0.; // starting value for chi_p
	double chi_V0_start = 0.; // starting value for chi_V0
	double lambda_path_start = 0.; // starting value for lambda_path
	double lambda_angle_start = 0.; // starting value for lambda_angle
	
	//step values for variables for selection: - check!!!
	double chi_pi_step = 1.; // step value for chi_pi
	double chi_p_step = 1.; // step value for chi_p
	double chi_V0_step = 1.; // step value for chi_V0
	double lambda_path_step = 1.; // step value for lambda_path
	double lambda_angle_step = 0.05; // step value for lambda_angle
	
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
	
	TH1D *hm0[NITER_CENT][NITER][NITER][NITER][NITER][NITER], *hm0_signal[NITER_CENT][NITER][NITER][NITER][NITER][NITER], *hm0_bckg[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	TH1D *hm0_before[NITER_CENT];
	TF1 *fitting_fnc[NITER_CENT][NITER][NITER][NITER][NITER][NITER], *backFcn[NITER_CENT][NITER][NITER][NITER][NITER][NITER], *signalFcn[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	char *int_fitting_fnc = new char[NITER];
	char *int_backFcn = new char[NITER];
	char *int_signalFcn = new char[NITER];

	double sum_sig[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	double sum_full[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	double entries_new[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	double entries_old[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	double efficiency[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	double ratio_SB[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	double ratio_SSB[NITER_CENT][NITER][NITER][NITER][NITER][NITER];

	double xmax[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	int bin_left[NITER_CENT][NITER][NITER][NITER][NITER][NITER];
	int bin_right[NITER_CENT][NITER][NITER][NITER][NITER][NITER];

	TCanvas *c1[NITER_CENT];
	char *int_canvas_c1 = new char[NITER_CENT];

	TFile *myFile_data = new TFile(inFile);

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		sprintf(int_canvas_c1,"canvas_c1_%d",iter_cent);
		c1[iter_cent] = new TCanvas(int_canvas_c1,int_canvas_c1,0,0,600,600);
		
		hm0_before[iter_cent] = (TH1D*)myFile_data->Get(Form("hm0_before_%d",iter_cent));
		
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = (TH1D*) myFile_data->Get(Form("hm0_cent%d_pi%d_p%d_V%d_path%d_angle%d",iter_cent,iter_pi,iter_p,iter_V0,iter_path,iter_angle));
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = (TH1D*)hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Clone("hm0sig_cent%d_pi%d_p%d_V%d_path%d_angle%d");
							hm0_bckg[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = (TH1D*)hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Clone("hm0bckg_cent%d_pi%d_p%d_V%d_path%d_angle%d");
							xmax[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetXaxis()->GetXmax();
						}
					}
				}
			}
		}
	}

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		//now it's getting comlicated :)
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							sprintf(int_fitting_fnc,"int_fitting_fnc_cent%d_pi%d_p%d_V%d_path%d_angle%d",iter_cent,iter_pi,iter_p,iter_V0,iter_path,iter_angle);
							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = new TF1(int_fitting_fnc,fitting_function,xmin,xmax[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle],8);
			
							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetNpx(500);
							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetLineWidth(4);
							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetLineColor(2); //red color for fitting line
							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3","pol4");
							Int_t npar = fitting_fncfitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetNumberFreeParameters();
							for (Int_t ip = 0; ip < npar; ++ip) fitting_fncfitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetParameter(ip,0);	

							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetParameter(0,10000);
							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetParameter(1,1.116);
							fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetParameter(2,0.002);
							
							//change to the canvas:
							c1[iter_cent]->cd();	
							c1[iter_cent]->Clear();
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetMinimum(1.0E-20);
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetXaxis()->SetRangeUser(1.07,1.17);
							hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Fit(int_fitting_fnc,"w","",1.086,1.18);
							
							Float_t mass = fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetParameter(npar-7);
							Float_t sigma = TMath::Abs(fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetParameter(npar-6));
						
							Double_t xmin_cut = mass - nsig_bckg * sigma;
							Double_t xmax_cut = mass + nsig_bckg * sigma;
						
							Int_t imin = hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->FindBin(xmin_cut);
							Int_t imax = hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->FindBin(xmax_cut);
							
							//improve the picture:
							sprintf(int_backFcn,"backFcn_cent%d_pi%d_p%d_V%d_path%d_angle%d",iter_cent,iter_pi,iter_p,iter_V0,iter_path,iter_angle);
							backFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = new TF1(int_backFcn,background,xmin,xmax[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle],5);
							backFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetLineColor(kMagenta);
							
							Double_t errs[200] = {0};
							Int_t nbins = hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetNbinsX();
							for (Int_t ib = 1; ib <= nbins; ++ib) 
							{
								if (ib < imin || ib > imax) errs[ib] = TMath::Sqrt(hm0_bckg[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetBinContent(ib));
							}
							hm0_bckg[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetError(errs);
							hm0_bckg[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetLineColor(kMagenta);
							hm0_bckg[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Fit(int_backFcn,"","same",1.086,1.17);
							
							sprintf(int_signalFcn,"signalFcn_cut1_cent%d_pi%d_p%d_V%d_path%d_angle%d",iter_cent,iter_pi,iter_p,iter_V0,iter_path,iter_angle);
							signalFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = new TF1(int_signalFcn,gaussian,xmin,xmax[iter_cent][iter],3);
							signalFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetLineColor(kBlue);
							signalFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetNpx(500);
							
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Sumw2();	
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Add(backFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle], -1);
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetLineColor(kBlack);
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetLineWidth(2);
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetMarkerStyle(2);
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetMarkerSize(1);
							signalFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->SetParameters(hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetMaximum(),1.1157,0.003);
							hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Fit(int_signalFcn,"","same",1.105,1.125);
							
							mass = signalFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetParameter(1);
							sigma = signalFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetParameter(2);
							
							Double_t left = mass - nsig * sigma;
							Double_t right = mass + nsig * sigma;
							left = TMath::Max(left,1.105);
							right = TMath::Min(right,1.13);
							
							bin_left[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->FindBin(left);
							bin_right[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->FindBin(right);
							
							entries_new[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->Integral(bin_left[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle],bin_right[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]);
							entries_old[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetEntries();
							
							int muh1 = hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetMaximum();
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
							legend->AddEntry(hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle],"Data","lp");
							legend->AddEntry(backFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle],"Background","l");
							legend->AddEntry(signalFcn[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle],"Signal","l");
							legend->AddEntry(fitting_fnc[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle],"Global Fit","l");
							legend->AddEntry(line_left,"Cut-off range","l");
							legend->Draw("same");
							
							TLatex *title1 = new TLatex(1.13,0.96*hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetMaximum(), omega_values[iter]); 
							title1->Draw("same");
						}
					}
				}
			}
		}
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							sum_sig[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = 0.;
							sum_full[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = 0.;
						}
					}
				}
			}
		}
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							efficiency[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = 100.*entries_new[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]/entries_old[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle];
			
							for (Int_t ib = bin_left[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]; ib <= bin_right[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]; ++ib) 
							{
								sum_sig[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] += hm0_signal[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetBinContent(ib);
								sum_full[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] += hm0[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]->GetBinContent(ib);
							}
							ratio_SB[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = sum_sig[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]/(sum_full[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] - sum_sig[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]);
							ratio_SSB[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] = sum_sig[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]/TMath::Sqrt(sum_full[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle]);
						}
					}
				}
			}
		}
	}
	
	ofstream output_file_selections;
	output_file_selections.open (Form("Selection_values_Cent_%d.txt",NITER_CENT));

//find the highest values for each centrality interval and the corresponding values of selections parameters:
	double max_value[NITER_CENT];
	int max_value_iter_pi[NITER_CENT];
	int max_value_iter_p[NITER_CENT];
	int max_value_iter_V0[NITER_CENT];
	int max_value_iter_path[NITER_CENT];
	int max_value_iter_angle[NITER_CENT];
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		max_value[iter_cent] = ratio_SSB[iter_cent][0][0][0][0][0];
		max_value_iter_pi[iter_cent] = 0;
		max_value_iter_p[iter_cent] = 0;
		max_value_iter_V0[iter_cent] = 0;
		max_value_iter_path[iter_cent] = 0;
		max_value_iter_angle[iter_cent] = 0;
		
		for(int iter_pi = 0; iter_pi < NITER; iter_pi++)
		{
			for(int iter_p = 0; iter_p < NITER; iter_p++)
			{
				for(int iter_V0 = 0; iter_V0 < NITER; iter_V0++)
				{
					for(int iter_path = 0; iter_path < NITER; iter_path++)
					{
						for(int iter_angle = 0; iter_angle < NITER; iter_angle++)
						{
							if(ratio_SSB[iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle] > max_value[iter_cent])
							{
								max_value[iter_cent] = ratio_SSB[iter_cent][iter_cent][iter_pi][iter_p][iter_V0][iter_path][iter_angle];
								max_value_iter_pi[iter_cent] = iter_pi;
								max_value_iter_p[iter_cent] = iter_p;
								max_value_iter_V0[iter_cent] = iter_V0;
								max_value_iter_path[iter_cent] = iter_path;
								max_value_iter_angle[iter_cent] = iter_angle;
							}
						}
					}
				}
			}
		}
		cout << "centrality interval = " << iter_cent << "; Highest SSB ratio = " << max_value[iter_cent] << endl;
		cout << "iter_pi = " << max_value_iter_pi[iter_cent] << "; chi_pi_value = " << chi_pi_value[max_value_iter_pi[iter_cent]] << endl;
		if(max_value_iter_pi[iter_cent] == 0 || max_value_iter_pi[iter_cent] == NITER-1)
		{
			cout << "iter_pi at the border of the range! : " << endl;
		}
		cout << "iter_p = " << max_value_iter_p[iter_cent] << "; chi_p_value = " << chi_p_value[max_value_iter_p[iter_cent]] << endl;
		if(max_value_iter_p[iter_cent] == 0 || max_value_iter_p[iter_cent] == NITER-1)
		{
			cout << "iter_p at the border of the range! : " << endl;
		}
		cout << "iter_V0 = " << max_value_iter_V0[iter_cent] << "; chi_V0_value = " << chi_V0_value[max_value_iter_V0[iter_cent]] << endl;
		if(max_value_iter_V0[iter_cent] == 0 || max_value_iter_V0[iter_cent] == NITER-1)
		{
			cout << "iter_V0 at the border of the range! : " << endl;
		}
		cout << "iter_path = " << max_value_iter_path[iter_cent] << "; lambda_path_value = " << lambda_path_value[max_value_iter_path[iter_cent]] << endl;
		if(max_value_iter_path[iter_cent] == 0 || max_value_iter_path[iter_cent] == NITER-1)
		{
			cout << "iter_path at the border of the range! : " << endl;
		}
		cout << "iter_angle = " << max_value_iter_angle[iter_cent] << "; lambda_angle_value = " << lambda_angle_value[max_value_iter_angle[iter_cent]] << endl;
		if(max_value_iter_angle[iter_cent] == 0 || max_value_iter_angle[iter_cent] == NITER-1)
		{
			cout << "iter_angle at the border of the range! : " << endl;
		}
		
		output_file_selections << chi_pi_value[max_value_iter_pi[iter_cent]] << chi_p_value[max_value_iter_p[iter_cent]] << chi_V0_value[max_value_iter_V0[iter_cent]] << chi_path_value[max_value_iter_path[iter_cent]] << chi_angle_value[max_value_iter_angle[iter_cent]] << "\n"; //write to file
	}
	
	//now draw the best stuff:
	
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		c1[iter_cent]->cd(iter+1);	
		c1[iter_cent]->Clear();
		hm0[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->Draw();
		hm0_bckg[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->Draw("same");
		hm0_signal[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->Draw("same");
		
		//backFcn[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->Draw("same");
		//signalFcn[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->Draw("same");
		//fitting_fnc[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->Draw("same");
							
		TLatex *title_SSB = new TLatex(1.13,0.86*hm0[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]])); 
		title_SSB->Draw("same");
		TLatex *title_efficiency = new TLatex(1.13,0.76*hm0[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]]->GetMaximum(), Form("eff. = %.2f",efficiency[iter_cent][max_value_iter_pi[iter_cent]][max_value_iter_p[iter_cent]][max_value_iter_V0[iter_cent]][max_value_iter_path[iter_cent]][max_value_iter_angle[iter_cent]])); 
		title_efficiency->Draw("same");
		
	}
	
	output_file_selections.close();
	
}
