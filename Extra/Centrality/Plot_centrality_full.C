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

void Plot_centrality_full(const int NITER_CENT = 10)
{
	/*TString infile_glaubfit = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_noDCA/glauber_qa_5104340_121.root";
	TString infile_final = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_noDCA/FINAL.root";
	TString infile_histocut = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_noDCA/HistoCutResult.root";
	TString infile_model_b = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_noDCA/B_average_model.root";
	*/
	
	TString infile_glaubfit = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_DCA_0_5/glauber_qa_7280464_20.root";
	TString infile_final = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_DCA_0_5/FINAL.root";
	TString infile_histocut = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_DCA_0_5/HistoCutResult.root";
	TString infile_model_b = "/scratch2/nazarova/CentralityFramework/Framework/results/Dataset2_PHSD_9GeV_b20_DCA_0_5/B_average_model.root";
	
	TH1D *B_vs_Centrality[NITER_CENT], *B_vs_Multiplicity[NITER_CENT];
		
	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
	TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);
	c2->Divide(2,1);
	TCanvas *c3 = new TCanvas("c3", "c3",0,0,600,600);
	TCanvas *c4 = new TCanvas("c4", "c4",0,0,600,300);
	TCanvas *c5 = new TCanvas("c5", "c5",0,0,600,600);
	TCanvas *c6 = new TCanvas("c6", "c6",0,0,600,600);
	gROOT->ForceStyle();
	gStyle->SetOptStat(0000);
	gStyle->SetOptTitle(0);
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
	const char *cent_interval[10] = {"0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
	const double centrality_min[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
	const double centrality_max[] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
	const int coltest[] = {kRed+2, kGreen+1, kBlue+2, kMagenta, kRed, kCyan+2, kBlue, kMagenta+2, kCyan, kGreen+3};
	
	TFile *myFile_data_glaubfit = new TFile(infile_glaubfit);
	TFile *myFile_data_final = new TFile(infile_final);
	TFile *myFile_data_histocut = new TFile(infile_histocut);
	TFile *myFile_data_model_b = new TFile(infile_model_b);
	
	TH1F *Glauber_fit = (TH1F*) myFile_data_glaubfit->Get("glaub_fit_histo");	
	TH1D *Multiplicity = (TH1D*) myFile_data_glaubfit->Get("hRefMultSTAR_DCA");
	
	TH2D *hVertex_Z_vs_Cent = (TH2D*) myFile_data_model_b->Get("hVertex_Z_vs_Cent");	
	
	Glauber_fit->SetYTitle("#frac{dN}{N_{ch}}");
	Glauber_fit->SetXTitle("N_{ch}");
	Glauber_fit->SetMarkerStyle(2);
	Glauber_fit->SetLineColor(kRed+1);
	Glauber_fit->SetMarkerColor(kRed+1);
	Glauber_fit->SetMarkerSize(2);
	Glauber_fit->SetLineWidth(2);
	
	Multiplicity->SetYTitle("#frac{dN}{N_{ch}}");
	Multiplicity->SetXTitle("N_{ch}");
	Multiplicity->SetMarkerStyle(2);
	Multiplicity->SetLineColor(kBlue+1);
	Multiplicity->SetMarkerColor(kBlue+1);
	Multiplicity->SetMarkerSize(2);
	Multiplicity->SetLineWidth(2);
	
	c1->cd();	
	TLegend *legend1=new TLegend(0.35,0.75,0.85,0.85);
	legend1->SetTextFont(42);
	legend1->SetTextSize(0.03);
	legend1->SetBorderSize(0);

	Multiplicity->GetYaxis()->SetRangeUser(1.,10.*Multiplicity->GetMaximum());
	c1->SetLogy();
	Multiplicity->GetXaxis()->SetRangeUser(0.,400.);
	
	Multiplicity->Draw();
	Glauber_fit->Draw("histsame");
	
	legend1->AddEntry(Multiplicity,"PHSD BiBi @9GeV, w/o DCA","l");
	//legend1->AddEntry(Glauber_fit,"glauber fit: f = 0.61, k = 1, M [20,350]","l"); //PHSD_9GeV_b20_noDCA
	legend1->AddEntry(Glauber_fit,"glauber fit: f = 0.0, k = 16, M [20,350]","l"); //Request 25
	legend1->Draw("same");
	
	TH1D *B_vs_Centrality_full = (TH1D*) myFile_data_final->Get("B_VS_CentralityClass 0%-100%");
	B_vs_Centrality_full->SetLineColor(kBlack);
	
	TH1D *B_vs_Multiplicity_full = (TH1D*) Multiplicity->Clone("hRefMultSTAR_2");
	B_vs_Multiplicity_full->SetLineColor(kBlack);
	B_vs_Multiplicity_full->SetYTitle("#frac{dN}{N_{ch}}");
	B_vs_Multiplicity_full->SetXTitle("N_{ch}");
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		B_vs_Centrality[iter_cent] = (TH1D*) myFile_data_final->Get(Form("B_VS_CentralityClass %.1f%%-%.1f%%",centrality_min[iter_cent],centrality_max[iter_cent]));	
		B_vs_Centrality[iter_cent]->SetLineColor(coltest[iter_cent]);
		B_vs_Multiplicity[iter_cent] = (TH1D*) myFile_data_histocut->Get(Form("CentralityClass %.1f%%-%.1f%%",centrality_min[iter_cent],centrality_max[iter_cent]));	
		B_vs_Multiplicity[iter_cent]->SetLineColor(coltest[iter_cent]);
	}
		
	c2->cd(1);
	TLegend *legend2_1 = new TLegend(0.75,0.4,0.85,0.88);
	
	legend2_1->SetTextFont(72);
	legend2_1->SetTextSize(0.04);
	legend2_1->SetBorderSize(0);
	B_vs_Centrality_full->Draw();
	legend2_1->AddEntry(B_vs_Centrality_full,"Full","l");
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		B_vs_Centrality[iter_cent]->Draw("same");
		legend2_1->AddEntry(B_vs_Centrality[iter_cent],cent_interval[iter_cent],"l");
		
	}	
	legend2_1->Draw();
	c2->cd(2);
	TLegend *legend2_2 = new TLegend(0.75,0.4,0.85,0.88);
	
	legend2_2->SetTextFont(72);
	legend2_2->SetTextSize(0.04);
	legend2_2->SetBorderSize(0);
	B_vs_Multiplicity_full->GetYaxis()->SetRangeUser(1.,10.*B_vs_Multiplicity_full->GetMaximum());
	gPad->SetLogy();
	B_vs_Multiplicity_full->GetXaxis()->SetRangeUser(0.,400.);
	B_vs_Multiplicity_full->Draw();
	legend2_2->AddEntry(B_vs_Multiplicity_full,"Full","l");
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		B_vs_Multiplicity[iter_cent]->Draw("same");
		legend2_2->AddEntry(B_vs_Multiplicity[iter_cent],cent_interval[iter_cent],"l");
		
	}	
	legend2_2->Draw();
	
	TH1D *B_average_vs_cent_model = (TH1D*) myFile_data_model_b->Get("B_average_VS_Centrality");
	TH1D *B_average_vs_cent = (TH1D*) myFile_data_final->Get("B_average_VS_Centrality");
	
	B_average_vs_cent_model->SetLineColor(kBlack);
	B_average_vs_cent_model->SetMarkerColor(kBlack);
	B_average_vs_cent_model->SetMarkerStyle(25);
	B_average_vs_cent_model->SetMarkerSize(2);
	B_average_vs_cent->SetLineColor(kRed+1);
	B_average_vs_cent->SetMarkerColor(kRed+1);
	B_average_vs_cent->SetMarkerStyle(23);
	B_average_vs_cent->SetMarkerSize(2);
	B_average_vs_cent->SetYTitle("<b> #pm #sigma_{b}, fm");
	B_average_vs_cent->SetXTitle("Centrality, [%]");
	
	c3->cd();	
	TLegend *legend3_1=new TLegend(0.15,0.72,0.45,0.92);
	legend3_1->SetTextFont(42);
	legend3_1->SetTextSize(0.06);
	legend3_1->SetBorderSize(0);

	B_average_vs_cent->Draw();
	B_average_vs_cent_model->Draw("same");
	
	legend3_1->AddEntry(B_average_vs_cent,"RECO","p");
	legend3_1->AddEntry(B_average_vs_cent_model,"Model","p");
	legend3_1->Draw("same");
	
	TH1D *Divide_hist = (TH1D*) B_average_vs_cent->Clone("Divide_hist");
	Divide_hist->Divide(B_average_vs_cent,B_average_vs_cent_model);
	
	c4->cd();	
	Divide_hist->SetMinimum(0.8);
	Divide_hist->SetMaximum(1.2);
	Divide_hist->SetYTitle("MC-Gl/Model");
	Divide_hist->Draw("HIST p");
	
	TLine *line1 = new TLine(0.,1.,100.,1.);
	line1->SetLineColor(kBlack);
	line1->SetLineWidth(2);
	line1->SetLineStyle(2);
	line1->Draw("same");
	
	TLine *line2 = new TLine(0.,0.95,100.,0.95);
	line2->SetLineColor(kBlack);
	line2->SetLineWidth(2);
	line2->SetLineStyle(2);
	line2->Draw("same");
	
	TLine *line3 = new TLine(0.,1.05,100.,1.05);
	line3->SetLineColor(kBlack);
	line3->SetLineWidth(2);
	line3->SetLineStyle(2);
	line3->Draw("same");
	
	TH1F *NCentr = (TH1F*) myFile_data_model_b->Get("NCentr");
	
	double centrality_bin[NITER_CENT];
	double centr_events[NITER_CENT];
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		centrality_bin[iter_cent] = centrality_min[iter_cent] + (centrality_max[iter_cent] - centrality_min[iter_cent])/2;
		centr_events[iter_cent] = NCentr->GetBinContent(iter_cent+1);
	}
	TGraph *NCentr_graph = new TGraph(NITER_CENT, centrality_bin, centr_events);
	NCentr_graph->SetName("NCentr");
	NCentr_graph->SetTitle("NCentr");
	NCentr_graph->GetXaxis()->SetTitle("Centrality, [%]");
	NCentr_graph->GetYaxis()->SetTitle("Events");
	NCentr_graph->SetLineColor(kRed+1);
	NCentr_graph->SetMarkerColor(kRed+1);
	NCentr_graph->SetMarkerSize(2);
	NCentr_graph->SetMarkerStyle(26);
	NCentr_graph->SetFillColor(40);
	
	c5->cd();		
	NCentr_graph->Draw("AB");
	
	c6->cd();		
	hVertex_Z_vs_Cent->Draw("lego2");
}
