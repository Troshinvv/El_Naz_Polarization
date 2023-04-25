#include <iostream>
#include <fstream> // std::ifstream

#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdGlobalPolarizationMC.h"
#include "TFile.h"

ClassImp(MpdGlobalPolarizationMC);

MpdGlobalPolarizationMC::MpdGlobalPolarizationMC(const char *name, const char *outputName) : MpdAnalysisTask2(name, outputName)
{
   readParameters(name);
   param("NITER_CENT", NITER_CENT, 4);
   param("NITER", NITER, 20);
   param("cent_cut_choice", cent_cut_choice, 0);
   param("cent_cut", cent_cut, 70.0);
   param("particle_choice", particle_choice, 3122);
}

void MpdGlobalPolarizationMC::UserInit()
{   
	// Prepare histograms etc.
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

	pdgCodeHyperon = particle_choice;
	if(pdgCodeHyperon == pdgCodeL0)    
	{
		cout << "You have chosen to analyze Lambda hyperons: " << " pdg: " << pdgCodeHyperon << endl;
		pdgCodeDaughter = pdgCodePr;
		massHyperon = massL0;
		massDaughter = massPr;
		cout << "massHyperon: " << massHyperon << endl;
		cout << "massDaughter: " << massDaughter << endl;
	}else if(pdgCodeHyperon == pdgCodeAL0)    
	{
		cout << "You have chosen to analyze anti-Lambda hyperons: " << " pdg: " << pdgCodeHyperon << endl;
		pdgCodeDaughter = pdgCodeAPr;
		massHyperon = massL0;
		massDaughter = massPr;
		cout << "massHyperon: " << massHyperon << endl;
		cout << "massDaughter: " << massDaughter << endl;
	}else
	{
		cout << "This pdg code for particle_choice is not defined! Please provide the definition in the code." << endl;
		exit(0);
	}
	
	if (NITER_CENT == 4)
	{		
		centrality_min = init_int_array(4, 0, 0, 10, 20, 50);
		centrality_max = init_int_array(4, 0, 10, 20, 50, 100);
		_CentrBins = init_double_array(5, 0, 0.,10.,20.,50.,100.);
	}else if (NITER_CENT == 7)
	{
		centrality_min = init_int_array(7, 0, 0, 10, 20, 30, 40, 50, 60);
		centrality_max = init_int_array(7, 0, 10, 20, 30, 40, 50, 60, 70);
		_CentrBins = init_double_array(8, 0, 0., 10., 20., 30., 40., 50., 60., 70.);
	}else if (NITER_CENT == 10)
	{
		centrality_min = init_int_array(10, 0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90);
		centrality_max = init_int_array(10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100);
		_CentrBins = init_double_array(11, 0, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.);
	}else 
	{
		cout << "This values of centrality bins is not defined! Please provide the definition in the code." << endl;
		exit(0);
	}

	// General QA
	mhCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);
	fOutputList->Add(mhCentrality);

	NCentr = new TH1D("NCentr","NCentr",NITER_CENT,_CentrBins);
	fOutputList->Add(NCentr);
	
	Resolution_EP1_true = new TH1D("Resolution_EP1_true","Resolution_EP1_true",NITER_CENT,_CentrBins);
	fOutputList->Add(Resolution_EP1_true);
	
	Resolution_EP1_exp = new TH1D("Resolution_EP1_exp","Resolution_EP1_exp",NITER_CENT,_CentrBins);
	fOutputList->Add(Resolution_EP1_exp);

	Lpolar_y = new TH1D*[NITER_CENT];
	Lpolar_y_prim = new TH1D*[NITER_CENT];
	PstarRP_hist = new TH1D*[NITER_CENT];
	PstarRP_hist_prim = new TH1D*[NITER_CENT];
	PstarEP_hist = new TH1D*[NITER_CENT];
	PstarEP_hist_prim = new TH1D*[NITER_CENT];
	
	ResEP1_true = new double[NITER_CENT];
	SubEvRes1 = new double[NITER_CENT];

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		
		Lpolar_y[iter_cent] = new TH1D(Form("LPolar_y_%d", iter_cent),Form("LPolar_y_%d", iter_cent),100,-1.,1.);
		fOutputList->Add(Lpolar_y[iter_cent]);
		Lpolar_y_prim[iter_cent] = new TH1D(Form("LPolar_y_prim_%d", iter_cent),Form("LPolar_y_prim_%d", iter_cent),100,-1.,1.);
		fOutputList->Add(Lpolar_y_prim[iter_cent]);
		PstarRP_hist[iter_cent] = new TH1D(Form("PstarRP_hist_%d", iter_cent),Form("PstarRP_hist_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarRP_hist[iter_cent]);
		PstarRP_hist_prim[iter_cent] = new TH1D(Form("PstarRP_hist_prim_%d", iter_cent),Form("PstarRP_hist_prim_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarRP_hist_prim[iter_cent]);
		PstarEP_hist[iter_cent] = new TH1D(Form("PstarEP_hist_%d", iter_cent),Form("PstarEP_hist_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarEP_hist[iter_cent]);
		PstarEP_hist_prim[iter_cent] = new TH1D(Form("PstarEP_hist_prim_%d", iter_cent),Form("PstarEP_hist_prim_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarEP_hist_prim[iter_cent]);
		
		ResEP1_true[iter_cent] = 0.;
		SubEvRes1[iter_cent] = 0.;
		
	}
}

void MpdGlobalPolarizationMC::ProcessEvent(MpdAnalysisEvent &event)
{
	if (!selectEvent(event)) 
	{ 
		return;
	}

	phiRP_mc = 0.;
	phiEP_mc = 0.;
	ResEP_mc = 0.;
	ResEPSub_mc = 0.;
	phiRP_mc = event.fMCEventHeader->GetRotZ();
	phiRP_mc = TMath::ATan2(TMath::Sin(phiRP_mc),TMath::Cos(phiRP_mc));
	phiEP_mc = event.fMpdEP.GetPhiEP_FHCal_F_all();
	ResEP_mc = TMath::Cos(phiEP_mc - phiRP_mc);
	ResEPSub_mc = TMath::Cos(event.fMpdEP.GetPhiEP_FHCal_S_all() - event.fMpdEP.GetPhiEP_FHCal_N_all());

	fillHistograms(event);
	
}

void MpdGlobalPolarizationMC::Finish()
{
   cout << "Finish() ..." << endl;

}

//--------------------------------------
bool MpdGlobalPolarizationMC::selectEvent(MpdAnalysisEvent &event)
{
	mMCTracks = event.fMCTrack;
	
	Centrality_tpc = event.getCentrTPC();
	
	if (Centrality_tpc < 0 || Centrality_tpc >= 100) // TPC centrality not defined
	{ 
		return false;
	}

	if (cent_cut_choice == 1)
	{
		cout << "Cutting out events with more than " << cent_cut << "% centrality!" << endl;
		if (Centrality_tpc > cent_cut) 
		{ 
			return false;
		}
	}

	mhCentrality->Fill(Centrality_tpc);

	return true;
}

void MpdGlobalPolarizationMC::fillHistograms(MpdAnalysisEvent &event)
{	
	NCentr->Fill(Centrality_tpc);
	int nMC = mMCTracks->GetEntriesFast();
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		if(Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue; 
		
		ResEP1_true[iter_cent] += ResEP_mc;
		SubEvRes1[iter_cent] += ResEPSub_mc;
		
		for (Int_t j = 0; j < nMC; ++j) 
		{
			TVector3 mom; 
			TVector3 mom_moth; 
			MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(j);
			mcTr->GetMomentum(mom);
			if (mcTr->GetPdgCode() == pdgCodeDaughter) 
			{
				int mcTr_MotherID = mcTr->GetMotherId();
				if (mcTr_MotherID < 0) continue; // if it's a primary particle (e.g. proton), we don't need it													
				MpdMCTrack* mcTr_Mother = (MpdMCTrack*) mMCTracks->UncheckedAt(mcTr_MotherID); // get the track, corresponding to the mother ID
				if(mcTr_Mother->GetPdgCode() == pdgCodeHyperon) // choose only the case, when it's the mother we need (e.g. Lambda)
				{
					mcTr_Mother->GetMomentum(mom_moth);
					int mcTr_Lam_MotherID = mcTr_Mother->GetMotherId(); // ID of the mother of hyperon (e.g. Lambda)
					double polar_x = 0.;
					double polar_y = 0.;
					double polar_z = 0.;
					double pol_length = 0.;
					int sign = 1;
					
					double weight_pol = mcTr_Mother->GetWeight();
					polar_x = mcTr_Mother->GetPolar(0);
					polar_y = mcTr_Mother->GetPolar(1);
					polar_z = mcTr_Mother->GetPolar(2);
					sign = TMath::Sign(1, mcTr_Mother->GetPdgCode());
						
					TVector3 polar_vector(polar_x, polar_y, polar_z); //polarization vector for distributions		
					if (phiRP_mc != 0.) polar_vector.RotateZ(-phiRP_mc);
					polar_x = weight_pol*polar_vector.X();
					polar_y = weight_pol*polar_vector.Y();
					polar_z = weight_pol*polar_vector.Z();
					pol_length = TMath::Sqrt(polar_x*polar_x + polar_y*polar_y + polar_z*polar_z);
					//for primary Lambda
					if (mcTr_Lam_MotherID < 0) 
					{
						Lpolar_y_prim[iter_cent]->Fill(polar_y);
					}
						
					Lpolar_y[iter_cent]->Fill(polar_y);
						
					//daughter (proton) values (MagThetaPhi)
					double p_prot = mom.Mag();
					double eta_prot = mom.Theta();
					double phi_prot = mom.Phi();
					//hyperon (lambda) values (MagThetaPhi)	
					double p_lam = mom_moth.Mag();
					double eta_lam = mom_moth.Theta();
					double phi_lam = mom_moth.Phi();
					TVector3 vPr, vLamb;
					vPr.SetMagThetaPhi(p_prot, eta_prot, phi_prot);
					vLamb.SetMagThetaPhi(p_lam, eta_lam, phi_lam);
						
					//calculate the costheta and phi from dataset:
					double cos_prot = 0.0; 
					double phi_prot_star = 0.0; 
					FindPolarAngle (vPr, vLamb, cos_prot, phi_prot_star);
					double phi_diff_RP = 0.;
					double phi_diff_EP = 0.;
					phi_diff_RP = phiRP_mc - phi_prot_star;
					phi_diff_EP = phiEP_mc - phi_prot_star;
					
					if (phi_diff_RP < 0) phi_diff_RP = phi_diff_RP + 2.*pi;
					if (phi_diff_EP < 0) phi_diff_EP = phi_diff_EP + 2.*pi;

					//fill distribution for all Lambda	
					PstarRP_hist[iter_cent]->Fill(phi_diff_RP);
					PstarEP_hist[iter_cent]->Fill(phi_diff_EP);
					
					//fill distribution for primary Lambda
					if (mcTr_Lam_MotherID < 0) 
					{
						PstarRP_hist_prim[iter_cent]->Fill(phi_diff_RP);
						PstarEP_hist_prim[iter_cent]->Fill(phi_diff_EP);
					}						
				} // mcTr_Mother->GetPdgCode() == pdgCodeHyperon
			} // mcTr->GetPdgCode() == pdgCodeDaughter		
		} // cycle over MC tracks
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		
		Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);
		Resolution_EP1_exp->SetBinContent(iter_cent+1,SubEvRes1[iter_cent]);
	}
}
double* MpdGlobalPolarizationMC::init_double_array (const int n, const double fmt...)
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

int* MpdGlobalPolarizationMC::init_int_array (const int n, const int fmt...)
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

void MpdGlobalPolarizationMC::FindPolarAngle(TVector3 &vPr, TVector3 &vLamb, double &cos_prot, double &phi_prot_star)
{
	cos_prot = 0.;
	phi_prot_star = 0.;
	
	TLorentzVector prLor, lambLor;
	prLor.SetVectM(vPr, massDaughter);
	lambLor.SetVectM(vLamb, massHyperon);
	TVector3 boostV;
	boostV = lambLor.BoostVector();
	boostV *= -1;
  	prLor.Boost(boostV);
	vPr = prLor.Vect();	
	
	cos_prot = vPr.CosTheta(); // cos theta
	phi_prot_star = vPr.Phi(); // phi
}
