/**
 * @file mOmega2_Selection_Mixing_MB.C
 * @author Elizaveta Nazarova
 * @brief Analysis macro for RECO GlobalPol wagon, selection mode, omega_2 selection
 * @version 0.1
 * @date 2023-07-15
 * Using obtained invariant mass distributions from the "selection" mode of GlobalPolarizationRECO wagon, performs background subtraction based on mixed obtained mixed background distribution to extract the signal
 * Make sure that the input values correspond to the ones in the config file used to run the wagon!
 * Finds the optimal value of omega_2 for full MB, corresponding to the maximum significance value
 * Writes the obtained value to the .txt file
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

#include <iostream> // std::cout
#include <fstream>  // std::ifstream
using namespace std;

double massL0 = 1.11568;
float pi = TMath::Pi();

void FitFile(int NITER, TString inFile, double xmin, double nsig, double &max_value_sig, int &max_value_iter);

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
/**
 * @brief Analysis macro for RECO GlobalPol wagon to obtain optimal omega_2 selection values using event mixing for background subtraction
 * 
 * @param inFile               name of merged input file (e.g. Output_GlopalPolRECO_omega2_selection.root)
 * @param configname           name of the config file (e.g. pGlobalPolRECO.txt)
 * @param outFile              name of output file for saving the obtained selection values
 * @param NITER                number of different values of parameters for selection (default: 10)
 * @param xmin                 minimal value of invariant mass (X-axis) for fitting (default: 1.086)
 * @param nsig                 value to obtain the range for signal extraction around the peak (+/- nsig*sigma) (default: 3.0)
 */
void mOmega2_Selection_Mixing_MB(TString inFile = "Output_GlopalPolRECO_omega2_selection.root", TString configname = "pGlobalPolRECO.txt", TString outFile = "Omega2Selection_values_MB_3sigma.txt", int NITER = 30, double xmin = 1.09, double nsig = 3.0)
{
	gROOT->LoadMacro("mpdloadlibs.C");
	gROOT->ProcessLine("mpdloadlibs()");

	cout << "Reading input text file: " << configname << endl;
	double omega_value_start; // starting values 
	double omega_step; // step values

	ifstream ifs(configname);
	if (!ifs) 
	{
		cout << "File " << configname << " can not be opened" << endl;
		return 1;
	}

	string a;
	double b;
	while (ifs.good()) 
	{
		ifs >> a;
		// if comment, skip to the end of line
		if (a.find_first_of('#') == 0 || a.find_first_of("//") == 0) 
		{
			ifs.ignore(999, '\n');
			continue;
		} 
		else if (a == "omega_value_start") 
		{
			ifs >> b;
			omega_value_start = b;
		}
		else if (a == "omega_step") 
		{
			ifs >> b;
			omega_step = b;
		}
		cout << a << "\t" << b << endl;
	}
	ifs.close();

	double omega_value[NITER];
	
	for(int iter = 0; iter < NITER; iter++)
	{
		omega_value[iter] = omega_value_start + omega_step*iter;
	}
	const char *omega_values_30bins[] = {"#omega_{2} > 1.0", "#omega_{2} > 1.1", "#omega_{2} > 1.2", "#omega_{2} > 1.3", "#omega_{2} > 1.4", "#omega_{2} > 1.5", "#omega_{2} > 1.6", "#omega_{2} > 1.7", "#omega_{2} > 1.8", "#omega_{2} > 1.9", "#omega_{2} > 2.0", "#omega_{2} > 2.1", "#omega_{2} > 2.2", "#omega_{2} > 2.3", "#omega_{2} > 2.4", "#omega_{2} > 2.5", "#omega_{2} > 2.6", "#omega_{2} > 2.7", "#omega_{2} > 2.8", "#omega_{2} > 2.9", "#omega_{2} > 3.0", "#omega_{2} > 3.1", "#omega_{2} > 3.2", "#omega_{2} > 3.3", "#omega_{2} > 3.4", "#omega_{2} > 3.5", "#omega_{2} > 3.6", "#omega_{2} > 3.7", "#omega_{2} > 3.8", "#omega_{2} > 3.9"}; 
	for (int i = 0; i < NITER; i++)
	{
		cout << "omega2_cutvalues[i] = " <<  omega_value[i] << "; omega_values[i] = " <<  omega_values[i] << endl;
	} 
	
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);

	ofstream output_file_omega;
	output_file_omega.open(outFile);
	
	double max_value_sig = 0.;
	int max_value_iter = 0;
				
	FitFile(NITER,inFile,xmin,nsig,max_value_sig,max_value_iter);	

	//write out the highest value and the corresponding omega_2 value;
	cout << "; Highest SSB ratio = " << max_value_sig << "; iter = " << max_value_iter << "; Omega_2 = " << omega_value[max_value_iter] << " --- " << omega_values[max_value_iter] << endl;
	if(max_value_iter == NITER-1) 
	{
		cout << "Omega_2 value at the last bin!!! " << "; iter = " << max_value_iter << endl;
	}
	if(max_value_iter == 0) 
	{
		cout << "Omega_2 value at the first bin!!! " << "; iter = " << max_value_iter << endl;
	}
		
	output_file_omega << omega_value[max_value_iter] << "\n"; //write to file
	output_file_omega.close();
	
	//now lets plot the final fits:
	TH1D *hm0_final, *hm0_signal_final, *hm0_bckg_final;
	TF1 *backFcn_final, *signalFcn_final;
	
	double sum_sig_final = 0.;
	double sum_full_final = 0.;
	double entries_new_final;
	double entries_old_final;
	double efficiency_final;
	double ratio_SB_final;
	double ratio_SSB_final;
	
	double xmax_final;
	int bin_left_final;
	int bin_right_final;

	TCanvas *c1_final = new TCanvas("c1_final","c1_final",0,0,600,600);
	TFile *myFile_data = new TFile(inFile);
	
	hm0_final = (TH1D*) myFile_data->Get(Form("hm0_full_%d",max_value_iter));
	hm0_signal_final = (TH1D*)hm0_final->Clone("hm0_fullsig_%d");
	hm0_bckg_final = (TH1D*) myFile_data->Get(Form("hm0_full_mixed_%d",max_value_iter));
	xmax_final = hm0_final->GetXaxis()->GetXmax();
	
	int bin_left_scale = hm0_bckg_final->FindBin(1.14);
	int bin_right_scale = hm0_bckg_final->FindBin(1.16);
	double sum_hm0_scale = 0.;
	double sum_hm0_mixed_scale = 0.;
	for (int ib = bin_left_scale; ib <= bin_right_scale; ++ib) 
	{
	sum_hm0_scale += hm0_final->GetBinContent(ib); 
	sum_hm0_mixed_scale += hm0_bckg_final->GetBinContent(ib);
	}
	hm0_bckg_final->Scale(sum_hm0_scale/sum_hm0_mixed_scale);	
	backFcn_final = new TF1("backFcn_final",background,xmin,xmax_final,5);
	backFcn_final->SetLineColor(kRed);
	backFcn_final->SetLineWidth(4);
	
	//change to the canvas:
	c1_final->cd();	
	hm0_final->SetMinimum(1.0E-20);
	hm0_final->GetXaxis()->SetRangeUser(1.07,1.17);
	hm0_final->SetYTitle("Entries");
	hm0_final->SetXTitle("M_{inv}, GeV/c^{2}");
	hm0_final->Draw();
	hm0_bckg_final->Draw("same");
	hm0_bckg_final->Fit("backFcn_final","q","same",xmin,xmax_final);
	
	signalFcn_final = new TF1("signalFcn_final",gaussian,xmin,xmax_final,3);
	
	signalFcn_final->SetLineColor(kBlue);
	signalFcn_final->SetNpx(500);
	
	hm0_signal_final->Sumw2();	
	hm0_signal_final->Add(backFcn_final, -1);
	hm0_signal_final->SetLineColor(kBlack);
	hm0_signal_final->SetLineWidth(2);
	hm0_signal_final->SetMarkerStyle(2);
	hm0_signal_final->SetMarkerSize(2);
	signalFcn_final->SetParameters(hm0_final->GetMaximum(),massL0,0.003);
	hm0_signal_final->Fit("signalFcn_final","q0","same",xmin,xmax_final);
	
	double mass = signalFcn_final->GetParameter(1);
	double sigma = signalFcn_final->GetParameter(2);
	double left = mass - nsig * sigma;
	double right = mass + nsig * sigma;
	
	cout << "left = " << left << endl;
	cout << "right = " << right << endl;
	
	double int_signal = signalFcn_final->Integral(left,right);
	double int_background = backFcn_final->Integral(left,right);
	
	cout << "signal integral (around peak) = " << int_signal << endl;
	cout << "background integral (around peak) = " << int_background << endl;
	
	bin_left_final = hm0_signal_final->FindBin(left);
	bin_right_final = hm0_signal_final->FindBin(right);
	
	entries_new_final = hm0_signal_final->Integral(bin_left_final,bin_right_final);
	entries_old_final = hm0_final->GetEntries();
	
	cout << "entries_new_final = " << entries_new_final << "; entries_old_final = " << entries_old_final << "; entries_old_final (integral) = " << hm0_final->Integral() << endl;
	
	TLine *line_left = new TLine(left,0,left,0.6*hm0_final->GetMaximum());
	line_left->SetLineColor(kBlack);
	line_left->SetLineWidth(3);
	line_left->SetLineStyle(2);
	line_left->Draw("same");
	TLine *line_right = new TLine(right,0,right,0.6*hm0_final->GetMaximum());
	line_right->SetLineColor(kBlack);
	line_right->SetLineWidth(3);
	line_right->SetLineStyle(2);
	line_right->Draw("same");
	TLegend *legend=new TLegend(0.6,0.65,0.88,0.9);
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->SetBorderSize(0);
	legend->AddEntry(hm0_final,"Data","lp");
	legend->AddEntry(backFcn_final,"Bckg fit","l");
	legend->AddEntry(line_left,"Cut-off (signal)","l");
	legend->Draw("same");	
	
	efficiency_final = 100.*entries_new_final/entries_old_final;
	
	for (int ib = bin_left_final; ib <= bin_right_final; ib++) 
	{
		sum_sig_final += hm0_signal_final->GetBinContent(ib);
		sum_full_final += hm0_final->GetBinContent(ib);
	}
	cout << "Sum of bins; Signal: " << sum_sig_final << "; Full: " << sum_full_final << endl;
	ratio_SB_final = sum_sig_final/(sum_full_final - sum_sig_final);
	ratio_SSB_final = sum_sig_final/TMath::Sqrt(sum_full_final);
	cout << "efficiency = " << efficiency_final << "; S/B = " << ratio_SB_final << "; S/(#sqrt{S+B}) = " << ratio_SSB_final << endl;
	
	TLatex *title_SSB = new TLatex(1.075,0.92*hm0_final->GetMaximum(), Form("#frac{S}{#sqrt{S+B}} = %.2f",ratio_SSB_final)); 
	title_SSB->Draw("same");
	TLatex *title_SB = new TLatex(1.075,0.76*hm0_final->GetMaximum(), Form("S/B = %.2f",ratio_SB_final)); 
	title_SB->Draw("same");
	TLatex *title_efficiency = new TLatex(1.075,0.66*hm0_final->GetMaximum(), Form("eff. = %.2f [%%]",efficiency_final)); 
	title_efficiency->Draw("same");
	
}
void FitFile(int NITER, TString inFile, double xmin, double nsig, double &max_value_sig, int &max_value_iter)
{
	TH1D *hm0_full[NITER], *hm0_signal_full[NITER], *hm0_bckg_full[NITER];
	TF1 *backFcn_full[NITER], *signalFcn_full[NITER];
	char *int_backFcn_full = new char[100];
	char *int_signalFcn_full = new char[100];

	double sum_sig_full[NITER];
	double sum_full_full[NITER];
	double entries_new_full[NITER];
	double entries_old_full[NITER];
	double efficiency_full[NITER];
	double ratio_SB_full[NITER];
	double ratio_SSB_full[NITER];
	
	double xmax_full[NITER];
	int bin_left_full[NITER];
	int bin_right_full[NITER];
	
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

	TFile *myFile_data = new TFile(inFile);

	for(int iter = 0; iter < NITER; iter++)
	{
		hm0_full[iter] = (TH1D*) myFile_data->Get(Form("hm0_full_%d",iter));
		hm0_signal_full[iter] = (TH1D*)hm0_full[iter]->Clone("hm0_fullsig_%d");
		hm0_bckg_full[iter] = (TH1D*) myFile_data->Get(Form("hm0_full_mixed_%d",iter));
		xmax_full[iter] = hm0_full[iter]->GetXaxis()->GetXmax();
		cout << "iter = " << iter << "; xmax = " << xmax_full[iter] << endl;
	}
	
	for(int iter = 0; iter < NITER; iter++)
	{
		int bin_left_scale = hm0_bckg_full[iter]->FindBin(1.14);
		int bin_right_scale = hm0_bckg_full[iter]->FindBin(1.16);
		double sum_hm0_scale = 0.;
		double sum_hm0_mixed_scale = 0.;
		for (int ib = bin_left_scale; ib <= bin_right_scale; ++ib) 
		{
			sum_hm0_scale += hm0_full[iter]->GetBinContent(ib); 
			sum_hm0_mixed_scale += hm0_bckg_full[iter]->GetBinContent(ib);
		}
		hm0_bckg_full[iter]->Scale(sum_hm0_scale/sum_hm0_mixed_scale);	
	
		sprintf(int_backFcn_full,"backFcn_%d",iter);
		backFcn_full[iter] = new TF1(int_backFcn_full,background,xmin,xmax_full[iter],5);
		backFcn_full[iter]->SetLineColor(kRed);
		backFcn_full[iter]->SetLineWidth(4);
		
		//change to the canvas:
		c2->cd(iter+1);	
		hm0_full[iter]->SetMinimum(1.0E-20);
		hm0_full[iter]->GetXaxis()->SetRangeUser(1.07,1.17);
		hm0_full[iter]->Draw();
		hm0_bckg_full[iter]->Draw("same");
		hm0_bckg_full[iter]->Fit(int_backFcn_full,"q","same",xmin,xmax_full[iter]);
		
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
		signalFcn_full[iter]->SetParameters(hm0_full[iter]->GetMaximum(),massL0,0.003);
		hm0_signal_full[iter]->Fit(int_signalFcn_full,"q","same",xmin,xmax_full[iter]);
		
		double mass = signalFcn_full[iter]->GetParameter(1);
		double sigma = signalFcn_full[iter]->GetParameter(2);
		
		double left = mass - nsig * sigma;
		double right = mass + nsig * sigma;
		
		bin_left_full[iter] = hm0_signal_full[iter]->FindBin(left);
		bin_right_full[iter] = hm0_signal_full[iter]->FindBin(right);
		
		entries_new_full[iter] = hm0_signal_full[iter]->Integral(bin_left_full[iter],bin_right_full[iter]);
		entries_old_full[iter] = hm0_full[iter]->GetEntries();
		
		TLine *line_left = new TLine(left,0,left,0.6*hm0_signal_full[iter]->GetMaximum());
		line_left->SetLineColor(kBlack);
		line_left->SetLineWidth(2);
		line_left->SetLineStyle(2);
		line_left->Draw("same");
		TLine *line_right = new TLine(right,0,right,0.6*hm0_signal_full[iter]->GetMaximum());
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
		legend->AddEntry(line_left,"Cut-off range","l");
		legend->Draw("same");
		
	}
	
	for(int iter = 0; iter < NITER; iter++)
	{
		sum_sig_full[iter] = 0.;
		sum_full_full[iter] = 0.;
	}
	for(int iter = 0; iter < NITER; iter++)
	{
		efficiency_full[iter] = 100.*entries_new_full[iter]/entries_old_full[iter];
		
		for (int ib = bin_left_full[iter]; ib <= bin_right_full[iter]; ++ib) {
			sum_sig_full[iter] += hm0_signal_full[iter]->GetBinContent(ib);
			sum_full_full[iter] += hm0_full[iter]->GetBinContent(ib);
		}
		cout << "Sum of bins; Signal: " << sum_sig_full[iter] << "; Full: " << sum_full_full[iter] << endl;
		ratio_SB_full[iter] = sum_sig_full[iter]/(sum_full_full[iter] - sum_sig_full[iter]);
		ratio_SSB_full[iter] = sum_sig_full[iter]/TMath::Sqrt(sum_full_full[iter]);
		cout << "efficiency = " << efficiency_full[iter] << "; S/B = " << ratio_SB_full[iter] << "; S/(#sqrt{S+B}) = " << ratio_SSB_full[iter] << endl;
	}
	
	double max_value = ratio_SSB_full[0];
	int max_value_iter_file = 0;
	
	for(int iter = 0; iter < NITER; iter++)
	{
		if(ratio_SSB_full[iter] > max_value)
		{
			max_value = ratio_SSB_full[iter];
			max_value_iter_file = iter;
		}
	}
	max_value_sig = max_value;
	max_value_iter = max_value_iter_file;
	cout << "File = " << inFile << "; Highest SSB ratio = " << max_value_sig << "; iter = " << max_value_iter << endl;
	if (c2) 
	{ 
		c2->Close(); 
		delete c2;
		c2 = 0; 
	}
	myFile_data->Close();
}
