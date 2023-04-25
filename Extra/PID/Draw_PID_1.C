//PID information (dE/dx and M^2)
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


void Draw_PID_1()
{

	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);	
	TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);
	TCanvas *c3 = new TCanvas("c3", "c3",0,0,600,600);
	TCanvas *c4 = new TCanvas("c4", "c4",0,0,600,600);
	TCanvas *c5 = new TCanvas("c5", "c5",0,0,600,600);
	
	//gStyle->SetOptFit();
	//gStyle->SetOptStat(1000000001);
	gStyle->SetOptStat(0000);
	//gStyle->SetOptTitle(0);
	TLatex latex;
	latex.SetNDC();
	TGaxis::SetMaxDigits(3);
	
  
//	TFile *myFile_data = new TFile("m2_TOF_phsd_9GeV.root");
	TFile *myFile_data = new TFile("/home/liza/Nextcloud/Work/polarization/output_req30/m2_TOF_PHSD_9_2GeV.root");

	TH2D *DEDX = (TH2D*)myFile_data->Get("DEDX");
	TH2D *DEDX2 = (TH2D*)myFile_data->Get("DEDX2");
	TH2D *Mass2 = (TH2D*)myFile_data->Get("Mass2");
	TH2D *DEDX_vsM2 = (TH2D*)myFile_data->Get("DEDX_vsM2");
	TH2D *Beta_hist = (TH2D*)myFile_data->Get("Beta_hist");
	
	DEDX->SetYTitle("dE/dx");
	DEDX->SetXTitle("p, GeV");
	DEDX->SetTitle("dE/dx as a function of momentum");
	
	DEDX2->SetYTitle("dE/dx");
	DEDX2->SetXTitle("q*p, GeV");
	DEDX2->SetTitle("dE/dx as a function of momentum");
	
	Mass2->SetYTitle("ToF m^{2}, GeV^{2}");
	Mass2->SetXTitle("p, GeV");
	Mass2->SetTitle("ToF m^{2} as a function of momentum");
	
	DEDX_vsM2->SetYTitle("dE/dx");
	DEDX_vsM2->SetXTitle("ToF m^{2}, GeV^{2}");
	
	Beta_hist->SetYTitle("#beta");
	Beta_hist->SetXTitle("q*p, GeV");
	
c1->cd();	
	gPad->SetLogz();
	
	DEDX->GetYaxis()->SetRangeUser(0.,100.); //put away when you redo the Fillm2 macro
	DEDX->Draw("colz");
	
	latex.SetTextSize(0.025);
	TLatex *pi = new TLatex(0.03,3.5, "#color[2]{#scale[1.2]{#pi}}"); pi->Draw();
	TLatex *kaon = new TLatex(0.2,6., "#color[2]{#scale[1.2]{K}}"); kaon->Draw();
	TLatex *proton = new TLatex(0.4,8., "#color[2]{#scale[1.2]{p}}"); proton->Draw();
	//TLatex *deutron = new TLatex(0.7,15000., "#color[2]{#scale[1.5]{d}}"); deutron->Draw();
c2->cd();	
	gPad->SetLogz();
	Mass2->Draw("colz");
	TLatex *pi_m = new TLatex(0.2,0.0, "#color[2]{#scale[1.5]{#pi}}"); pi_m->Draw();
	TLatex *kaon_m = new TLatex(0.3,0.22, "#color[2]{#scale[1.5]{K}}"); kaon_m->Draw();
	TLatex *proton_m = new TLatex(0.5,0.8, "#color[2]{#scale[1.5]{p}}"); proton_m->Draw();

c3->cd();	
	gPad->SetLogz();
	
	DEDX2->GetYaxis()->SetRangeUser(0.,100.); //put away when you redo the Fillm2 macro
	DEDX2->Draw("colz");
	
	pi->Draw();
	kaon->Draw();
	proton->Draw();
	TLatex *anti_pi = new TLatex(-0.2,3.5, "#color[2]{#scale[1.2]{#pi^{-}}}"); anti_pi->Draw();
	TLatex *anti_kaon = new TLatex(-0.4,6., "#color[2]{#scale[1.2]{K^{-}}}"); anti_kaon->Draw();
	TLatex *anti_proton = new TLatex(-0.6,8., "#color[2]{#scale[1.2]{#barp}}"); anti_proton->Draw();
	
c4->cd();
	gPad->SetLogz();
	
	DEDX_vsM2->GetYaxis()->SetRangeUser(0.,100.); //put away when you redo the Fillm2 macro
	DEDX_vsM2->Draw("colz");
	
c5->cd();
	gPad->SetLogz();
	
	Beta_hist->Draw("colz");
	
}
