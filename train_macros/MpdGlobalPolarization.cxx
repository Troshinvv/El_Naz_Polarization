#include <iostream>
#include <fstream> // std::ifstream

#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdGlobalPolarization.h"
#include "TFile.h"

ClassImp(MpdGlobalPolarization);

MpdGlobalPolarization::MpdGlobalPolarization(const char *name, const char *outputName) : MpdAnalysisTask2(name, outputName)
{
   readParameters(name);
   param("mZvtxCut", mZvtxCut, 130.0);
   param("mNofHitsCut", mNofHitsCut, 10);
   param("mEtaCut", mEtaCut, 0.5);
   param("mPtminCut", mPtminCut, 0.1);
   param("mDcaCut", mDcaCut, 2.0);
   param("NITER_CENT", NITER_CENT, 4);
   param("NITER", NITER, 20);
   param("cent_cut_choice", cent_cut_choice, 0);
   param("cent_cut", cent_cut, 70.0);
   param("particle_choice", particle_choice, "Lambda");
}

void MpdGlobalPolarization::UserInit()
{   
	// Prepare histograms etc.
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting
	
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
		return;
	}

	// General QA
	mhEvents = new TH1D("hEvents", "Number of events", 10, 0., 10.);
	fOutputList->Add(mhEvents);
	mhVertex = new TH1F("hVertex", "Event vertex distribution", 100, -200., 200.);
	fOutputList->Add(mhVertex);
	mhCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);
	fOutputList->Add(mhCentrality);

	NCentr = new TH1D("NCentr","NCentr",NITER_CENT,_CentrBins);
	fOutputList->Add(NCentr);
	
	Resolution_EP1_true = new TH1D("Resolution_EP1_true","Resolution_EP1_true",NITER_CENT,_CentrBins);
	fOutputList->Add(Resolution_EP1_true);
	
	Resolution_EP1_exp = new TH1D("Resolution_EP1_exp","Resolution_EP1_exp",NITER_CENT,_CentrBins);
	fOutputList->Add(Resolution_EP1_exp);

	Lpolar_y = (TProfile**) malloc(sizeof(TProfile*) * NITER_CENT);
	Lpolar_y_prim = (TProfile**) malloc(sizeof(TProfile*) * NITER_CENT);
	PstarRP_hist = (TH1D**) malloc(sizeof(TH1D*) * NITER_CENT);
	PstarRP_hist_prim = (TH1D**) malloc(sizeof(TH1D*) * NITER_CENT);
	PstarEP_hist = (TH1D**) malloc(sizeof(TH1D*) * NITER_CENT);
	PstarEP_hist_prim = (TH1D**) malloc(sizeof(TH1D*) * NITER_CENT);
	
	ResEP1_true = (double*) malloc(sizeof(double) * NITER_CENT);
	SubEvRes1 = (double*) malloc(sizeof(double) * NITER_CENT);

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		//P_{y} component of polarization vector for full hyperons (primary+secondary)  
		Lpolar_y[iter_cent] = new TProfile(Form("LPolar_y_%d", iter_cent),Form("LPolar_y_%d", iter_cent),100,-1.,1.);
		fOutputList->Add(Lpolar_y[iter_cent]);
		//P_{y} component of polarization vector for primary hyperons
		Lpolar_y_prim[iter_cent] = new TProfile(Form("LPolar_y_prim_%d", iter_cent),Form("LPolar_y_prim_%d", iter_cent),100,-1.,1.);
		fOutputList->Add(Lpolar_y_prim[iter_cent]);
		//Angular distributions for full hyperons from MCTracks (w.r.t. reaction plane)
		PstarRP_hist[iter_cent] = new TH1D(Form("PstarRP_hist_%d", iter_cent),Form("PstarRP_hist_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarRP_hist[iter_cent]);
		//Angular distributions for primary hyperons from MCTracks (w.r.t. reaction plane)
		PstarRP_hist_prim[iter_cent] = new TH1D(Form("PstarRP_hist_prim_%d", iter_cent),Form("PstarRP_hist_prim_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarRP_hist_prim[iter_cent]);
		//Angular distributions for full hyperons from MCTracks (w.r.t. event plane)
		PstarEP_hist[iter_cent] = new TH1D(Form("PstarEP_hist_%d", iter_cent),Form("PstarEP_hist_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarEP_hist[iter_cent]);
		//Angular distributions for primary hyperons from MCTracks (w.r.t. event plane)
		PstarEP_hist_prim[iter_cent] = new TH1D(Form("PstarEP_hist_prim_%d", iter_cent),Form("PstarEP_hist_prim_%d", iter_cent),NITER,0.,2.*pi);
		fOutputList->Add(PstarEP_hist_prim[iter_cent]);
		
		ResEP1_true[iter_cent] = 0.;
		SubEvRes1[iter_cent] = 0.;
		
	}
}

void MpdGlobalPolarization::ProcessEvent(MpdAnalysisEvent &event)
{
	if (!selectEvent(event)) 
	{ 
		return;
	}

	phiRP_mc = 0.;
	phiEP_mc = 0.;
	ResEP_mc = 0.;
	ResEPSub_mc = 0.;
	b0 = event.fMCEventHeader->GetB();
	phiRP_mc = event.fMCEventHeader->GetRotZ();
	phiRP_mc = TMath::ATan2(TMath::Sin(phiRP_mc),TMath::Cos(phiRP_mc));
	phiEP_mc = event.fMpdEP.GetPhiEP_FHCal_F_all();
	ResEP_mc = TMath::Cos(phiEP_mc - phiRP_mc);
	ResEPSub_mc = TMath::Cos(event.fMpdEP.GetPhiEP_FHCal_S_all() - event.fMpdEP.GetPhiEP_FHCal_N_all());

	//cout << "phiRP_mc = " << phiRP_mc << "; phiEP_mc = " << phiEP_mc << endl;
	//cout << "ResEP_mc = " << ResEP_mc << "; ResEPSub_mc = " << ResEPSub_mc << endl;

	fillHistograms(event);
	
}

void MpdGlobalPolarization::Finish()
{
   cout << "Finish() ..." << endl;

}

//--------------------------------------
bool MpdGlobalPolarization::selectEvent(MpdAnalysisEvent &event)
{
	mhEvents->Fill(0.5); // Number of full events

	mMCTracks = event.fMCTrack;
	
	int nMC = mMCTracks->GetEntriesFast();
	int nTrMc = 0;
	for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) 
	{
		MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
		if (pr->GetMotherId() == -1) 
		{
			nTrMc++;
		}
	}
	if (nTrMc <= 2*209) 
	{ // Just nucleons of Bi+Bi --> Request30-PHSD
		cout << "nTrMc clause - empty event!" << endl;
		cout << " nTrMc " << nTrMc << "; nMC " << nMC << endl;
		return false;
	}

	mhEvents->Fill(1.5); // Number of events with filled vertex

	if (!event.fVertex) // if even vertex not filled, skip event
	{ 
		return false;
	}
   
	MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
	vertex->Position(mPrimaryVertex);

	if (mPrimaryVertex.Z() == 0) // not reconstructed (==0)
	{ 
		return false;
	}

	if (fabs(mPrimaryVertex.Z()) > mZvtxCut) // beyond the limits
	{ 
		return false;
	}
	
	mhEvents->Fill(2.5); // Number of events after vertex checks (where the vertex was filled, reconstructed and is within the limits (less than the vertex cut))
	
	Centrality_tpc = event.getCentrTPC();
	//cout << "cen = " << cen << endl;
	
	if (Centrality_tpc < 0 || Centrality_tpc >= 100) // TPC centrality not defined
	{ 
		return false;
	}

	mhEvents->Fill(3.5); // Number of events, satisfying all the vertex criteria and having defined centrality

	mhVertex->Fill(mPrimaryVertex.Z());

	mhCentrality->Fill(Centrality_tpc);

	return true;
}

void MpdGlobalPolarization::fillHistograms(MpdAnalysisEvent &event)
{	
	NCentr->Fill(Centrality_tpc);
	int nMC = mMCTracks->GetEntriesFast();

	//cout << "phiRP_mc = " << phiRP_mc << "; phiEP_mc = " << phiEP_mc << endl;
	//cout << "ResEP_mc = " << ResEP_mc << "; ResEPSub_mc = " << ResEPSub_mc << endl;
	
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
			if (particle_choice == "Lambda")
			{	
				if (mcTr->GetPdgCode() == pdgCodePr) 
				{
					int mcTr_MotherID = mcTr->GetMotherId();
					if (mcTr_MotherID < 0) continue; // if it's a primary proton, we don't need it													
					MpdMCTrack* mcTr_Mother = (MpdMCTrack*) mMCTracks->UncheckedAt(mcTr_MotherID); // get the track, corresponding to the ID of protons' mother
					if(mcTr_Mother->GetPdgCode() == pdgCodeL0) // choose only the case, when Lambda is the mother
					{
						mcTr_Mother->GetMomentum(mom_moth);
						int mcTr_Lam_MotherID = mcTr_Mother->GetMotherId(); // ID of the mother of Lambda
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
						
						//cout << "polar_x = " << polar_x << "; polar_y = " << polar_y << "; polar_z = " << polar_z << endl;
							
						TVector3 polar_changed(polar_x, polar_y, polar_z); //polarization vector for distributions		
						if (phiRP_mc != 0.) polar_changed.RotateZ(-phiRP_mc);
						polar_x = weight_pol*polar_changed.X();
						polar_y = weight_pol*polar_changed.Y();
						polar_z = weight_pol*polar_changed.Z();
						pol_length = TMath::Sqrt(polar_x*polar_x + polar_y*polar_y + polar_z*polar_z);
						//for primary Lambda
						if (mcTr_Lam_MotherID < 0) 
						{
							Lpolar_y_prim[iter_cent]->Fill(polar_y,1);
						}
							
						Lpolar_y[iter_cent]->Fill(polar_y,1);
							
						//proton values (MagThetaPhi)
						double p_prot = mom.Mag();
						double eta_prot = mom.Theta();
						double phi_prot = mom.Phi();
						//lambda values (MagThetaPhi)	
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
					} // mcTr_Mother->GetPdgCode() == pdgCodeL0
						
				} // mcTr->GetPdgCode() == pdgCodePr	
			}else if (particle_choice == "ALambda")
			{
				if (mcTr->GetPdgCode() == pdgCodeAPr) 
				{
					int mcTr_MotherID = mcTr->GetMotherId();
					if (mcTr_MotherID < 0) continue; // if it's a primary proton, we don't need it													
					MpdMCTrack* mcTr_Mother = (MpdMCTrack*) mMCTracks->UncheckedAt(mcTr_MotherID); // get the track, corresponding to the ID of protons' mother
					if(mcTr_Mother->GetPdgCode() == pdgCodeAL0) // choose only the case, when Lambda is the mother
					{
						mcTr_Mother->GetMomentum(mom_moth);
						int mcTr_Lam_MotherID = mcTr_Mother->GetMotherId(); // ID of the mother of Lambda
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
							
						TVector3 polar_changed(polar_x, polar_y, polar_z); //polarization vector for distributions		
						if (phiEP_mc != 0.) polar_changed.RotateZ(-phiEP_mc);
						polar_x = weight_pol*polar_changed.X();
						polar_y = weight_pol*polar_changed.Y();
						polar_z = weight_pol*polar_changed.Z();
						pol_length = TMath::Sqrt(polar_x*polar_x + polar_y*polar_y + polar_z*polar_z);
						//for primary Lambda
						if (mcTr_Lam_MotherID < 0) 
						{
							Lpolar_y_prim[iter_cent]->Fill(polar_y,1);
						}
							
						Lpolar_y[iter_cent]->Fill(polar_y,1);
						Lpolar_y[iter_cent]->SetYTitle("Entries");
						Lpolar_y[iter_cent]->SetXTitle("P_{y}");
						Lpolar_y[iter_cent]->SetLineColor(kBlack);	
							
						//proton values (MagThetaPhi)
						double p_prot = mom.Mag();
						double eta_prot = mom.Theta();
						double phi_prot = mom.Phi();
						//lambda values (MagThetaPhi)	
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
					} // mcTr_Mother->GetPdgCode() == pdgCodeAL0
						
				} // mcTr->GetPdgCode() == pdgCodeAPr	
			}else
			{
				cout << "This particle choice is not defined! Please provide the definition in the code." << endl;
			}
		}
	}
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{		
		
		Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);
		Resolution_EP1_exp->SetBinContent(iter_cent+1,SubEvRes1[iter_cent]);
	}
}
double* MpdGlobalPolarization::init_double_array (const int n, const double fmt...)
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

int* MpdGlobalPolarization::init_int_array (const int n, const int fmt...)
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

//Get costheta and phi of the daughter particle
void MpdGlobalPolarization::FindPolarAngle(TVector3 &vPr, TVector3 &vLamb, double &cos_prot, double &phi_prot_star)
{
	//Compute azimuthal angle of proton in the lambda frame
	cos_prot = 0.;
	phi_prot_star = 0.;
	
	TLorentzVector prLor, lambLor;
      
	prLor.SetVectM(vPr, 0.938272);
	lambLor.SetVectM(vLamb, 1.11568);

	//standard Lorentz boost
	TVector3 boostV;
	boostV = lambLor.BoostVector();
	boostV *= -1;
  	
  	prLor.Boost(boostV);
	vPr = prLor.Vect();	
	
	//calculating the azimuthal angle of proton in the lambda frame (phi*):
	
	cos_prot = vPr.CosTheta(); // cos theta
	phi_prot_star = vPr.Phi();
}

void runit (TString output){
   gROOT->LoadMacro("mpdloadlibs.C");
   gROOT->ProcessLine("mpdloadlibs()");

   MpdAnalysisManager man("ManagerAnal",-1) ;
   man.InputFileList("list.txt") ;
   man.ReadBranches("*") ; 
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;

   MpdEventPlaneAll pEP("pEP","pEP") ;
   man.AddTask(&pEP) ;
   
   MpdGlobalPolarization pGlobalPol("pGlobalPol",output) ;
   man.AddTask(&pGlobalPol) ;

   man.Process() ;
}
