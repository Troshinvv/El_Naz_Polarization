//Creates different plots for the completed centrality calibration
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

void Plot_Selections_check_histograms(const int NITER_CENT = 10, TString infile = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_DCA_0_5/glauber_qa_7280464_20.root")
{		
	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	c1->Divide(2,2);
	TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);
	c2->Divide(2,2);
	TCanvas *c3 = new TCanvas("c3", "c3",0,0,600,600);
	c3->Divide(2,2);
	TCanvas *c4 = new TCanvas("c4", "c4",0,0,600,600);
	c4->Divide(2,2);
	TCanvas *c5 = new TCanvas("c5", "c5",0,0,600,600);
	c5->Divide(2,2);
	TCanvas *c6 = new TCanvas("c6", "c6",0,0,600,600);
	c6->Divide(2,2);
	TCanvas *c7 = new TCanvas("c7", "c7",0,0,600,600);
	c7->Divide(2,2);
	TCanvas *c8 = new TCanvas("c8", "c8",0,0,600,600);
	c8->Divide(2,2);
	TCanvas *c9 = new TCanvas("c9", "c9",0,0,600,600);
	c9->Divide(2,2);
	TCanvas *c10 = new TCanvas("c10", "c10",0,0,600,600);
	c10->Divide(2,2);
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	const char *cent_interval[4] = {"0-10%","10-20%","20-50%","50-100%"};
	
	TFile *myFile_data = new TFile(infile);
	
	TH1D *Dca_pion[NITER_CENT], *Dca_proton[NITER_CENT], *Chi_pion[NITER_CENT], *Chi_proton[NITER_CENT], *Dca_lambda[NITER_CENT], *Chi_lambda[NITER_CENT], *Dca_v0[NITER_CENT], *Chi_v0[NITER_CENT], *Path_hist[NITER_CENT], *Angle_hist[NITER_CENT];
	TH1D *Dca_pion_true[NITER_CENT], *Dca_proton_true[NITER_CENT], *Chi_pion_true[NITER_CENT], *Chi_proton_true[NITER_CENT], *Dca_lambda_true[NITER_CENT], *Chi_lambda_true[NITER_CENT], *Dca_v0_true[NITER_CENT], *Chi_v0_true[NITER_CENT], *Path_hist_true[NITER_CENT], *Angle_hist_true[NITER_CENT];	
	
	for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Dca_pion[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_pion_%d",iter_cent));	
		Dca_pion[iter_cent]->SetLineWidth(2);
		Dca_proton[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_proton_%d",iter_cent));	
		Dca_proton[iter_cent]->SetLineWidth(2);
		Chi_pion[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_pion_%d",iter_cent));	
		Chi_pion[iter_cent]->SetLineWidth(2);
		Chi_proton[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_proton_%d",iter_cent));	
		Chi_proton[iter_cent]->SetLineWidth(2);
		Dca_lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_lambda_%d",iter_cent));	
		Dca_lambda[iter_cent]->SetLineWidth(2);
		Chi_lambda[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_lambda_%d",iter_cent));	
		Chi_lambda[iter_cent]->SetLineWidth(2);
		Dca_v0[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_v0_%d",iter_cent));	
		Dca_v0[iter_cent]->SetLineWidth(2);
		Chi_v0[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_v0_%d",iter_cent));	
		Chi_v0[iter_cent]->SetLineWidth(2);
		Path_hist[iter_cent] = (TH1D*) myFile_data->Get(Form("Path_hist_%d",iter_cent));	
		Path_hist[iter_cent]->SetLineWidth(2);
		Angle_hist[iter_cent] = (TH1D*) myFile_data->Get(Form("Angle_hist_%d",iter_cent));	
		Angle_hist[iter_cent]->SetLineWidth(2);
		
		Dca_pion_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_pion_true_%d",iter_cent));	
		Dca_pion_true[iter_cent]->SetLineWidth(2);
		Dca_pion_true[iter_cent]->SetLineColor(kRed+1);
		Dca_proton_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_proton_true_%d",iter_cent));	
		Dca_proton_true[iter_cent]->SetLineWidth(2);
		Dca_proton_true[iter_cent]->SetLineColor(kRed+1);
		Chi_pion_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_pion_true_%d",iter_cent));	
		Chi_pion_true[iter_cent]->SetLineWidth(2);
		Chi_pion_true[iter_cent]->SetLineColor(kRed+1);
		Chi_proton_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_proton_true_%d",iter_cent));	
		Chi_proton_true[iter_cent]->SetLineWidth(2);
		Chi_proton_true[iter_cent]->SetLineColor(kRed+1);
		Dca_lambda_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_lambda_true_%d",iter_cent));	
		Dca_lambda_true[iter_cent]->SetLineWidth(2);
		Dca_lambda_true[iter_cent]->SetLineColor(kRed+1);
		Chi_lambda_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_lambda_true_%d",iter_cent));	
		Chi_lambda_true[iter_cent]->SetLineWidth(2);
		Chi_lambda_true[iter_cent]->SetLineColor(kRed+1);
		Dca_v0_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Dca_v0_true_%d",iter_cent));	
		Dca_v0_true[iter_cent]->SetLineWidth(2);
		Dca_v0_true[iter_cent]->SetLineColor(kRed+1);
		Chi_v0_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Chi_v0_true_%d",iter_cent));	
		Chi_v0_true[iter_cent]->SetLineWidth(2);
		Chi_v0_true[iter_cent]->SetLineColor(kRed+1);
		Path_hist_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Path_hist_true_%d",iter_cent));	
		Path_hist_true[iter_cent]->SetLineWidth(2);
		Path_hist_true[iter_cent]->SetLineColor(kRed+1);
		Angle_hist_true[iter_cent] = (TH1D*) myFile_data->Get(Form("Angle_hist_true_%d",iter_cent));	
		Angle_hist_true[iter_cent]->SetLineWidth(2);
		Angle_hist_true[iter_cent]->SetLineColor(kRed+1);
	}
		
	for (int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		//DCA pion histo
		c1->cd(iter_cent+1);	
		c1->SetLogy();
		Dca_pion[iter_cent]->GetXaxis()->SetRangeUser(0.,40.);
		Dca_pion[iter_cent]->Draw();
		Dca_pion_true[iter_cent]->Draw("same");
		//Dca_pion_true[iter_cent]->Draw();
		TPaveText *t1_1 = new TPaveText(.13,.75,.35,.84);
		t1_1->SetTextAlign(22);
		t1_1->SetTextColor(kRed+2);
		t1_1->SetTextFont(72);
		t1_1->SetTextSize(0.04);
		t1_1->Paint("NDC");
		TText *t1_1_1 = t1_1->AddText(cent_interval[iter_cent]);
		t1_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend1=new TLegend(0.75,0.75,0.9,0.85);
			legend1->SetTextFont(42);
			legend1->SetTextSize(0.03);
			legend1->SetBorderSize(0);
			legend1->AddEntry(Dca_pion[iter_cent],"Full","l");
			legend1->AddEntry(Dca_pion_true[iter_cent],"True","l"); 
			legend1->Draw("same");
		}
		//DCA proton histo
		c2->cd(iter_cent+1);	
		c2->SetLogy();
		Dca_proton[iter_cent]->GetXaxis()->SetRangeUser(0.,40.);
		Dca_proton[iter_cent]->Draw();
		Dca_proton_true[iter_cent]->Draw("same");
		//Dca_proton_true[iter_cent]->Draw();
		TPaveText *t2_1 = new TPaveText(.13,.75,.35,.84);
		t2_1->SetTextAlign(22);
		t2_1->SetTextColor(kRed+2);
		t2_1->SetTextFont(72);
		t2_1->SetTextSize(0.04);
		t2_1->Paint("NDC");
		TText *t2_1_1 = t2_1->AddText(cent_interval[iter_cent]);
		t2_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend2=new TLegend(0.75,0.75,0.9,0.85);
			legend2->SetTextFont(42);
			legend2->SetTextSize(0.03);
			legend2->SetBorderSize(0);
			legend2->AddEntry(Dca_proton[iter_cent],"Full","l");
			legend2->AddEntry(Dca_proton_true[iter_cent],"True","l"); 
			legend2->Draw("same");
		}
		//Chi pion histo
		c3->cd(iter_cent+1);	
		c3->SetLogy();
		Chi_pion[iter_cent]->GetXaxis()->SetRangeUser(0.,1000.);
		Chi_pion[iter_cent]->Draw();
		Chi_pion_true[iter_cent]->Draw("same");
		//Chi_pion_true[iter_cent]->Draw();
		TPaveText *t3_1 = new TPaveText(.13,.75,.35,.84);
		t3_1->SetTextAlign(22);
		t3_1->SetTextColor(kRed+2);
		t3_1->SetTextFont(72);
		t3_1->SetTextSize(0.04);
		t3_1->Paint("NDC");
		TText *t3_1_1 = t3_1->AddText(cent_interval[iter_cent]);
		t3_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend3=new TLegend(0.75,0.75,0.9,0.85);
			legend3->SetTextFont(42);
			legend3->SetTextSize(0.03);
			legend3->SetBorderSize(0);
			legend3->AddEntry(Chi_pion[iter_cent],"Full","l");
			legend3->AddEntry(Chi_pion_true[iter_cent],"True","l"); 
			legend3->Draw("same");
		}
		//Chi pion histo
		c4->cd(iter_cent+1);
		c4->SetLogy();	
		Chi_proton[iter_cent]->GetXaxis()->SetRangeUser(0.,1000.);
		Chi_proton[iter_cent]->Draw();
		Chi_proton_true[iter_cent]->Draw("same");
		//Chi_proton_true[iter_cent]->Draw();
		TPaveText *t4_1 = new TPaveText(.13,.75,.35,.84);
		t4_1->SetTextAlign(22);
		t4_1->SetTextColor(kRed+2);
		t4_1->SetTextFont(72);
		t4_1->SetTextSize(0.04);
		t4_1->Paint("NDC");
		TText *t4_1_1 = t4_1->AddText(cent_interval[iter_cent]);
		t4_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend4=new TLegend(0.75,0.75,0.9,0.85);
			legend4->SetTextFont(42);
			legend4->SetTextSize(0.03);
			legend4->SetBorderSize(0);
			legend4->AddEntry(Chi_proton[iter_cent],"Full","l");
			legend4->AddEntry(Chi_proton_true[iter_cent],"True","l"); 
			legend4->Draw("same");
		}
		//DCA lambda histo
		c5->cd(iter_cent+1);	
		c5->SetLogy();
		Dca_lambda[iter_cent]->GetXaxis()->SetRangeUser(0.,40.);
		Dca_lambda[iter_cent]->Draw();
		Dca_lambda_true[iter_cent]->Draw("same");
		//Dca_lambda_true[iter_cent]->Draw();
		TPaveText *t5_1 = new TPaveText(.13,.75,.35,.84);
		t5_1->SetTextAlign(22);
		t5_1->SetTextColor(kRed+2);
		t5_1->SetTextFont(72);
		t5_1->SetTextSize(0.04);
		t5_1->Paint("NDC");
		TText *t5_1_1 = t5_1->AddText(cent_interval[iter_cent]);
		t5_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend5=new TLegend(0.75,0.75,0.9,0.85);
			legend5->SetTextFont(42);
			legend5->SetTextSize(0.03);
			legend5->SetBorderSize(0);
			legend5->AddEntry(Dca_lambda[iter_cent],"Full","l");
			legend5->AddEntry(Dca_lambda_true[iter_cent],"True","l"); 
			legend5->Draw("same");
		}
		//Chi lambda histo
		c6->cd(iter_cent+1);	
		c6->SetLogy();
		Chi_lambda[iter_cent]->GetXaxis()->SetRangeUser(0.,1000.);
		Chi_lambda[iter_cent]->Draw();
		Chi_lambda_true[iter_cent]->Draw("same");
		//Chi_lambda_true[iter_cent]->Draw();
		TPaveText *t6_1 = new TPaveText(.13,.75,.35,.84);
		t6_1->SetTextAlign(22);
		t6_1->SetTextColor(kRed+2);
		t6_1->SetTextFont(72);
		t6_1->SetTextSize(0.04);
		t6_1->Paint("NDC");
		TText *t6_1_1 = t6_1->AddText(cent_interval[iter_cent]);
		t6_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend6=new TLegend(0.75,0.75,0.9,0.85);
			legend6->SetTextFont(42);
			legend6->SetTextSize(0.03);
			legend6->SetBorderSize(0);
			legend6->AddEntry(Chi_lambda[iter_cent],"Full","l");
			legend6->AddEntry(Chi_lambda_true[iter_cent],"True","l"); 
			legend6->Draw("same");
		}
		//DCA V0 histo
		c7->cd(iter_cent+1);
		c7->SetLogy();	
		Dca_v0[iter_cent]->Draw();
		Dca_v0_true[iter_cent]->Draw("same");
		//Dca_v0_true[iter_cent]->Draw();
		TPaveText *t7_1 = new TPaveText(.13,.75,.35,.84);
		t7_1->SetTextAlign(22);
		t7_1->SetTextColor(kRed+2);
		t7_1->SetTextFont(72);
		t7_1->SetTextSize(0.04);
		t7_1->Paint("NDC");
		TText *t7_1_1 = t7_1->AddText(cent_interval[iter_cent]);
		t7_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend7=new TLegend(0.75,0.75,0.9,0.85);
			legend7->SetTextFont(42);
			legend7->SetTextSize(0.03);
			legend7->SetBorderSize(0);
			legend7->AddEntry(Dca_v0[iter_cent],"Full","l");
			legend7->AddEntry(Dca_v0_true[iter_cent],"True","l"); 
			legend7->Draw("same");
		}
		//Chi V0 histo
		c8->cd(iter_cent+1);	
		c8->SetLogy();
		Chi_v0[iter_cent]->Draw();
		Chi_v0_true[iter_cent]->Draw("same");
		//Chi_v0_true[iter_cent]->Draw();
		TPaveText *t8_1 = new TPaveText(.13,.75,.35,.84);
		t8_1->SetTextAlign(22);
		t8_1->SetTextColor(kRed+2);
		t8_1->SetTextFont(72);
		t8_1->SetTextSize(0.04);
		t8_1->Paint("NDC");
		TText *t8_1_1 = t8_1->AddText(cent_interval[iter_cent]);
		t8_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend8=new TLegend(0.75,0.75,0.9,0.85);
			legend8->SetTextFont(42);
			legend8->SetTextSize(0.03);
			legend8->SetBorderSize(0);
			legend8->AddEntry(Chi_v0[iter_cent],"Full","l");
			legend8->AddEntry(Chi_v0_true[iter_cent],"True","l"); 
			legend8->Draw("same");
		}
		//Path histo
		c9->cd(iter_cent+1);
		c9->SetLogy();	
		Path_hist[iter_cent]->GetXaxis()->SetRangeUser(0.,50.);
		Path_hist[iter_cent]->Draw();
		Path_hist_true[iter_cent]->Draw("same");
		//Path_hist_true[iter_cent]->Draw();
		TPaveText *t9_1 = new TPaveText(.13,.75,.35,.84);
		t9_1->SetTextAlign(22);
		t9_1->SetTextColor(kRed+2);
		t9_1->SetTextFont(72);
		t9_1->SetTextSize(0.04);
		t9_1->Paint("NDC");
		TText *t9_1_1 = t9_1->AddText(cent_interval[iter_cent]);
		t9_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend9=new TLegend(0.75,0.75,0.9,0.85);
			legend9->SetTextFont(42);
			legend9->SetTextSize(0.03);
			legend9->SetBorderSize(0);
			legend9->AddEntry(Path_hist[iter_cent],"Full","l");
			legend9->AddEntry(Path_hist_true[iter_cent],"True","l"); 
			legend9->Draw("same");
		}
		//Angle histo
		c10->cd(iter_cent+1);	
		c10->SetLogy();
		Angle_hist[iter_cent]->Draw();
		Angle_hist_true[iter_cent]->Draw("same");
		//Angle_hist_true[iter_cent]->Draw();
		TPaveText *t10_1 = new TPaveText(.13,.75,.35,.84);
		t10_1->SetTextAlign(22);
		t10_1->SetTextColor(kRed+2);
		t10_1->SetTextFont(72);
		t10_1->SetTextSize(0.04);
		t10_1->Paint("NDC");
		TText *t10_1_1 = t10_1->AddText(cent_interval[iter_cent]);
		t10_1->Draw("same");
		if(iter_cent == 0)
		{
			TLegend *legend10=new TLegend(0.75,0.75,0.9,0.85);
			legend10->SetTextFont(42);
			legend10->SetTextSize(0.03);
			legend10->SetBorderSize(0);
			legend10->AddEntry(Angle_hist[iter_cent],"Full","l");
			legend10->AddEntry(Angle_hist_true[iter_cent],"True","l"); 
			legend10->Draw("same");
		}
				
	}	
	
}
