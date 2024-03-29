#include <iostream>
#include <fstream> // std::ifstream

#include "MpdVertex.h"
#include "MpdEvent.h"
#include "TFile.h"
#include "MpdGlobalPolarization_RECO.h"
#include "Struct_L0_v1.h"
ClassImp(MpdGlobalPolarization_RECO);

MpdGlobalPolarization_RECO::MpdGlobalPolarization_RECO(const char *name, const char *outputName, const char *analysis_choice, const char *selection_choice) : MpdAnalysisTask2(name, outputName)
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
   param("MCFile", MCFile, "");
   param("sigM", sigM, 4.0);
   param("sigE", sigE, 4.0);
   param("energy", energy, 9.0);
   param("coef", coef, 1.0);
   param("generator", generator, "PHSD");
   param("tracking", tracking, "CFHM");
   param("NITER_Selections", NITER_Selections, 30);
   param("omega_start", omega_start, 1.4);
   param("omega_step", omega_step, 0.1);

   param("chi_pi_start", chi_pi_start, 2.0);
   param("chi_p_start", chi_p_start, 2.0);
   param("chi_V0_start", chi_V0_start, 2.0);
   param("lambda_path_start", lambda_path_start, 1.6);
   param("lambda_angle_start", lambda_angle_start, 0.06);
   param("chi_pi_step", chi_pi_step, 0.4);
   param("chi_p_step", chi_p_step, 0.4);
   param("chi_V0_step", chi_V0_step, 0.4);
   param("lambda_path_step", lambda_path_step, 0.4);
   param("lambda_angle_step", lambda_angle_step, 0.02);

   param("selections_values", selections_values, "");

   this->analysis_choice = analysis_choice;
   this->selection_choice = selection_choice;
}

void MpdGlobalPolarization_RECO::UserInit()
{   
	cout << analysis_choice << endl;
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
	
	cout<<"---------- BEFORE PID ----------"<<endl;
	newPid = new MpdPid(sigM, sigE, energy, coef, generator, tracking, "pikapr");

	// Creating the necessary histograms:
	// Only the ones required for "selection" or "analysis" will be saved
	mhEvents = new TH1D("hEvents", "Number of events", 10, 0., 10.);
	mhVertex = new TH1F("hVertex", "Event vertex distribution", 100, -200., 200.);
	mhCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);

	NCentr = new TH1D("NCentr","NCentr",NITER_CENT,_CentrBins);
	Resolution_EP1_true = new TH1D("Resolution_EP1_true","Resolution_EP1_true",NITER_CENT,_CentrBins);
	Resolution_EP1_exp = new TH1D("Resolution_EP1_exp","Resolution_EP1_exp",NITER_CENT,_CentrBins);

	hMassL = new TH1D("hMassL", "Lambda mass", 50, 1.070, 1.170);
	hMassLsig = new TH1D("hMassLsig", "Lambda mass (signal)", 50, 1.070, 1.170);
	hMassLbkg = new TH1D("hMassLbkg", "Lambda mass (bckg.)", 50, 1.070, 1.170);
	Lpolar_full = new TH1D("Lpolar_full", "Lpolar_full", 100, -1., 1.);
	hPIDflag = new TH1D("hPIDflag", "PID flags", 12, 0, 12);
	hPdg = new TH1D("hPdg", "PdgCode if is not Pion", 1000, -2500, 2500);
	hLambFlag = new TH1D("hLambFlag","Flags for lambda", 12, 0, 12);
	hRecognitionFlag = new TH1D("hRecognitionFlag","Flags for Recognition", 10, 0, 10);
	hLambPTH = new TH1D("hLambPTH","Flags for lambdaPTH", 12, 0, 12);
	hAngle = new TH2D("hAngle","Acollinearity angle vs Pt", 50, 0, 1, 60, 0, 90);
	hProbTrueP = new TH2D("hProbTrueP","Probability for true Protons", 50, 0, 1.1, 50, 0, 1.1);
	hProbP = new TH2D("hProbfalseP","Probability for Pions and identification Protons", 50, 0, 1.1, 50, 0, 1.1);
	hProbTruePi = new TH2D("hProbTruePi","Probability for true Pions", 50, 0, 1.1, 50, 0, 1.1);
	hProbPi = new TH2D("hProbfalsePi","Probability for Protons and identification Pions", 50, 0, 1.1, 50, 0, 1.1);
	hProbTrueK = new TH2D("hProbTrueK","Probability for true Kaons", 50, 0, 1.1, 50, 0, 1.1);


	results = new TTree("event","Event");
	results->Branch("b0",&b0,"b0/F"); //impact parameter
	results->Branch("ntr",&ntr,"ntr/I"); // number of tracks selected for analysis
	TBranch *br = results->Branch("l0","std::vector<L0>", &vLambdas); //lambda candidates
	results->Branch("ptetayl0","std::vector<tuple<float,float,float> >", &vvvLpt); //lambda phase space (MC)
	results->Branch("nLamb",&nLamb,"nLamb/I"); //number of Lambda (in collection)
	results->Branch("nLamb_MC",&nLamb_MC,"nLamb_MC/I"); //number of Lambda (MC)


	if(analysis_choice == "analysis")
	{
		fOutputList->Add(mhEvents);
		fOutputList->Add(mhVertex);
		fOutputList->Add(mhCentrality);
		fOutputList->Add(NCentr);
		fOutputList->Add(Resolution_EP1_true);
		fOutputList->Add(Resolution_EP1_exp);
		fOutputList->Add(hMassL);
		fOutputList->Add(hMassLsig);
		fOutputList->Add(hMassLbkg);
		fOutputList->Add(Lpolar_full);
		fOutputList->Add(hPIDflag);
		fOutputList->Add(hPdg);
		fOutputList->Add(hLambFlag);
		fOutputList->Add(hRecognitionFlag);
		fOutputList->Add(hLambPTH);
		fOutputList->Add(hAngle);
		fOutputList->Add(hProbTrueP);
		fOutputList->Add(hProbP);
		fOutputList->Add(hProbTruePi);
		fOutputList->Add(hProbPi);
		fOutputList->Add(hProbTrueK);

		hm0_Full = new TH1D("hm0_Full", "hm0_Full", 100, 1.07, 1.17);
		fOutputList->Add(hm0_Full);
		hm0_before_full = new TH1D("hm0_before_full", "hm0_before_full", 100, 1.07, 1.17);
		fOutputList->Add(hm0_before_full);
		hm0_before = new TH1D*[NITER_CENT];
		hm0_after = new TH1D*[NITER_CENT];
		Lpolar = new TH1D*[NITER_CENT];
		Lpolar_prim = new TH1D*[NITER_CENT];
		PstarEP_hist = new TH1D*[NITER_CENT];
		PstarEP_hist_prim = new TH1D*[NITER_CENT];
		PstarRP_hist = new TH1D*[NITER_CENT];
		PstarRP_hist_prim = new TH1D*[NITER_CENT];
		PstarRP_hist_MC = new TH1D*[NITER_CENT];
		PstarRP_hist_MC_prim = new TH1D*[NITER_CENT];
		Dca_pion = new TH1D*[NITER_CENT];
		Dca_proton = new TH1D*[NITER_CENT];
		Chi_pion = new TH1D*[NITER_CENT];
		Chi_proton = new TH1D*[NITER_CENT];
		Dca_lambda = new TH1D*[NITER_CENT];
		Chi_lambda = new TH1D*[NITER_CENT];
		Dca_v0 = new TH1D*[NITER_CENT];
		Chi_v0 = new TH1D*[NITER_CENT];
		Path_hist = new TH1D*[NITER_CENT];
		Angle_hist = new TH1D*[NITER_CENT];
		hm0 = new TH1D**[NITER_CENT];
		angle_min = (double*) malloc(sizeof(double) * NITER);
		angle_max = (double*) malloc(sizeof(double) * NITER);

		double step_angle = (xmax_anglemax - xmin_anglemin)/NITER;
		cout << "xmin_anglemin = " << xmin_anglemin << "; xmax_anglemax = " << xmax_anglemax << endl;
		for(int iter = 0; iter < NITER; iter++)
		{
			angle_min[iter] = xmin_anglemin + step_angle*iter;
			angle_max[iter] = xmin_anglemin + step_angle*(iter+1);
			cout << "iter = " << iter << "; angle_min = " << angle_min[iter] << "; angle_max = " << angle_max[iter] << endl;
		}

		if(selection_choice == "omega2")
		{
			omega_value = (double*) malloc(sizeof(double) * NITER_CENT);
			//reading the omega_2 values from the file:
			cout << "Topology selection using omega_2 parameter" << endl;
			ifstream selections_file;
			selections_file.open(selections_values);
			if(selections_file.fail())
			{
				cout << "File with selection values does not exist! Please run the 'selection' choice first! Exiting... " << endl;
				exit(0);
			}
			selections_file >> omega_value_full;
			Int_t nlines = 0;
			while (true) 
			{
				if(nlines > NITER_CENT)
				{
					cout << "Too many lines in the selections file! Check it! " << endl;
					exit(0);
				}
				selections_file >> omega_value[nlines];
				if (!selections_file.good()) break;
				nlines++;
			}
			printf("found %d points\n",nlines);
			if(nlines != NITER_CENT)
			{
				cout << "Amount of omega_2 selection cuts is not equal to the centrality bins: " << "; nlines = " << nlines << "; NITER_CENT = " << NITER_CENT << endl;
				exit(0);
			}
			selections_file.close();
			cout << "omega_value_full =  " << omega_value_full << endl;
			for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
			{
				cout << "omega_value[iter_cent] =  " << omega_value[iter_cent] << endl;
			}
		}else if(selection_choice == "chi")
		{
			chi_pi_value = (double*) malloc(sizeof(double) * NITER_CENT);
			chi_p_value = (double*) malloc(sizeof(double) * NITER_CENT);
			chi_V0_value = (double*) malloc(sizeof(double) * NITER_CENT);
			lambda_path_value = (double*) malloc(sizeof(double) * NITER_CENT);
			lambda_angle_value = (double*) malloc(sizeof(double) * NITER_CENT);
			//reading the chi selection values from the file:
			cout << "Topology selection using chi selection parameters" << endl;
			ifstream selections_file;
			selections_file.open(selections_values);
			if(selections_file.fail())
			{
				cout << "File with selection values does not exist! Please run the 'selection' choice first! Exiting... " << endl;
				exit(0);
			}
			//selections_file >> omega_value_full; // should i have a set for "full"?
			Int_t nlines = 0;
			while (true) 
			{
				if(nlines > NITER_CENT)
				{
					cout << "Too many lines in the selections file! Check it! " << endl;
					exit(0);
				}
				selections_file >> chi_pi_value[nlines] >> chi_p_value[nlines] >> chi_V0_value[nlines] >> lambda_path_value[nlines] >> lambda_angle_value[nlines];
				if (!selections_file.good()) break;
				nlines++;
			}
			printf("found %d points\n",nlines);
			if(nlines != NITER_CENT)
			{
				cout << "Amount of omega_2 selection cuts is not equal to the centrality bins: " << "; nlines = " << nlines << "; NITER_CENT = " << NITER_CENT << endl;
				exit(0);
			}
			selections_file.close();
			//cout << "omega_value_full =  " << omega_value_full << endl;
			for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
			{
				cout << "iter_cent = " << iter_cent << "; chi_pi_value[iter_cent] =  " << chi_pi_value[iter_cent] << "; chi_p_value[iter_cent] =  " << chi_p_value[iter_cent] << "; chi_V0_value[iter_cent] =  " << chi_V0_value[iter_cent] << "; lambda_path_value[iter_cent] =  " << lambda_path_value[iter_cent] << "; lambda_angle_value[iter_cent] =  " << lambda_angle_value[iter_cent] << endl;
			}
		}else 
		{
			cout << "No such analysis_choice defined yet! Please choose either 'omega2' or 'chi'! Exiting... " << endl;
			exit(0);
		}

		for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
		{
			hm0[iter_cent] = new TH1D*[NITER];
			hm0_before[iter_cent] = new TH1D(Form("hm0_before_%d", iter_cent),Form("hm0_before_%d", iter_cent), 100, 1.07, 1.17);
			fOutputList->Add(hm0_before[iter_cent]);
			hm0_after[iter_cent] = new TH1D(Form("hm0_after_%d", iter_cent),Form("hm0_after_%d", iter_cent), 100, 1.07, 1.17);
			fOutputList->Add(hm0_after[iter_cent]);
			Lpolar[iter_cent] = new TH1D(Form("Lpolar_%d", iter_cent),Form("Lpolar_%d", iter_cent), 100, -1., 1.);
			fOutputList->Add(Lpolar[iter_cent]);
			Lpolar_prim[iter_cent] = new TH1D(Form("Lpolar_prim_%d", iter_cent),Form("Lpolar_prim_%d", iter_cent), 100, -1., 1.);
			fOutputList->Add(Lpolar_prim[iter_cent]);
			PstarEP_hist[iter_cent] = new TH1D(Form("PstarEP_hist_%d", iter_cent),Form("PstarEP_hist_%d", iter_cent), NITER, 0., 2.*pi);
			fOutputList->Add(PstarEP_hist[iter_cent]);
			PstarEP_hist_prim[iter_cent] = new TH1D(Form("PstarEP_hist_prim_%d", iter_cent),Form("PstarEP_hist_prim_%d", iter_cent), NITER, 0., 2.*pi);
			fOutputList->Add(PstarEP_hist_prim[iter_cent]);
			PstarRP_hist[iter_cent] = new TH1D(Form("PstarRP_hist_%d", iter_cent),Form("PstarRP_hist_%d", iter_cent), NITER, 0., 2.*pi);
			fOutputList->Add(PstarRP_hist[iter_cent]);
			PstarRP_hist_prim[iter_cent] = new TH1D(Form("PstarRP_hist_prim_%d", iter_cent),Form("PstarRP_hist_prim_%d", iter_cent), NITER, 0., 2.*pi);
			fOutputList->Add(PstarRP_hist_prim[iter_cent]);
			PstarRP_hist_MC[iter_cent] = new TH1D(Form("PstarRP_hist_MC_%d", iter_cent),Form("PstarRP_hist_MC_%d", iter_cent), NITER, 0., 2.*pi);
			fOutputList->Add(PstarRP_hist_MC[iter_cent]);
			PstarRP_hist_MC_prim[iter_cent] = new TH1D(Form("PstarRP_hist_MC_prim_%d", iter_cent),Form("PstarRP_hist_MC_prim_%d", iter_cent), NITER, 0., 2.*pi);
			fOutputList->Add(PstarRP_hist_MC_prim[iter_cent]);
			Dca_pion[iter_cent] = new TH1D(Form("Dca_pion_%d", iter_cent),Form("Dca_pion_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Dca_pion[iter_cent]);
			Dca_proton[iter_cent] = new TH1D(Form("Dca_proton_%d", iter_cent),Form("Dca_proton_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Dca_proton[iter_cent]);
			Chi_pion[iter_cent] = new TH1D(Form("Chi_pion_%d", iter_cent),Form("Chi_pion_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Chi_pion[iter_cent]);
			Chi_proton[iter_cent] = new TH1D(Form("Chi_proton_%d", iter_cent),Form("Chi_proton_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Chi_proton[iter_cent]);
			Dca_lambda[iter_cent] = new TH1D(Form("Dca_lambda_%d", iter_cent),Form("Dca_lambda_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Dca_lambda[iter_cent]);
			Chi_lambda[iter_cent] = new TH1D(Form("Chi_lambda_%d", iter_cent),Form("Chi_lambda_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Chi_lambda[iter_cent]);
			Dca_v0[iter_cent] = new TH1D(Form("Dca_v0_%d", iter_cent),Form("Dca_v0_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Dca_v0[iter_cent]);
			Chi_v0[iter_cent] = new TH1D(Form("Chi_v0_%d", iter_cent),Form("Chi_v0_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Chi_v0[iter_cent]);
			Path_hist[iter_cent] = new TH1D(Form("Path_hist_%d", iter_cent),Form("Path_hist_%d", iter_cent), NITER, 0., 100.);
			fOutputList->Add(Path_hist[iter_cent]);
			Angle_hist[iter_cent] = new TH1D(Form("Angle_hist_%d", iter_cent),Form("Angle_hist_%d", iter_cent), NITER, 0., 1.6);
			fOutputList->Add(Angle_hist[iter_cent]);
			
			for(int iter = 0; iter < NITER; iter++)
			{
				hm0[iter_cent][iter] = new TH1D(Form("hm0_%d_%d", iter_cent, iter),Form("hm0_%d_%d", iter_cent, iter), 100, 1.07, 1.17);
				fOutputList->Add(hm0[iter_cent][iter]);
			}
			
		}

	}else if(analysis_choice == "selection")
	{
		if(selection_choice == "omega2")
		{
			cout << "Topology selection using omega_2 parameter" << endl;
			hm0_full = new TH1D*[NITER_Selections];
			hm0_before_full = new TH1D("hm0_before_full", "hm0_before_full", 100, 1.07, 1.17);
			fOutputList->Add(hm0_before_full);
			hm0_before = new TH1D*[NITER_CENT];
			hm0 = new TH1D**[NITER_CENT];
			omega_value = (double*) malloc(sizeof(double) * NITER_Selections);
			for(int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++)
			{
				omega_value[iter_sel] = omega_start + omega_step*iter_sel;
				cout << "iter_sel = " << iter_sel << "; omega_value = " << omega_value[iter_sel] << endl;
			}
			for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
			{
				hm0[iter_cent] = new TH1D*[NITER_Selections];
				hm0_before[iter_cent] = new TH1D(Form("hm0_before_%d", iter_cent),Form("hm0_before_%d", iter_cent), 100, 1.07, 1.17);
				fOutputList->Add(hm0_before[iter_cent]);
				
				for(int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++)
				{
					hm0[iter_cent][iter_sel] = new TH1D(Form("hm0_%d_%d", iter_cent, iter_sel),Form("hm0_%d_%d", iter_cent, iter_sel), 100, 1.07, 1.17);
					fOutputList->Add(hm0[iter_cent][iter_sel]);
				}
				
			}
			for(int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++)
			{
				hm0_full[iter_sel] = new TH1D(Form("hm0_full_%d", iter_sel),Form("hm0_full_%d", iter_sel), 100, 1.07, 1.17);
				fOutputList->Add(hm0_full[iter_sel]);
			}
		}else if(selection_choice == "chi")
		{
			cout << "Topology selection using chi parameters" << endl;
			hm0_before_full = new TH1D("hm0_before_full", "hm0_before_full", 100, 1.07, 1.17);
			fOutputList->Add(hm0_before_full);
			hm0_before = new TH1D*[NITER_CENT];
			chi_pi_value = (double*) malloc(sizeof(double) * NITER_Selections);
			chi_p_value = (double*) malloc(sizeof(double) * NITER_Selections);
			chi_V0_value = (double*) malloc(sizeof(double) * NITER_Selections);
			lambda_path_value = (double*) malloc(sizeof(double) * NITER_Selections);
			lambda_angle_value = (double*) malloc(sizeof(double) * NITER_Selections);
			for(int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++)
			{
				chi_pi_value[iter_sel] = chi_pi_start + chi_pi_step*iter_sel;
				chi_p_value[iter_sel] = chi_p_start + chi_p_step*iter_sel;
				chi_V0_value[iter_sel] = chi_V0_start + chi_V0_step*iter_sel;
				lambda_path_value[iter_sel] = lambda_path_start + lambda_path_step*iter_sel;
				lambda_angle_value[iter_sel] = lambda_angle_start + lambda_angle_step*iter_sel;
				cout << "iter_sel = " << iter_sel << "; chi_pi_value = " << chi_pi_value[iter_sel] << "; chi_p_value = " << chi_p_value[iter_sel] << "; chi_V0_value = " << chi_V0_value[iter_sel] << "; lambda_path_value = " << lambda_path_value[iter_sel] << "; lambda_angle_value = " << lambda_angle_value[iter_sel] << endl;
			}
			for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
			{
				hm0_before[iter_cent] = new TH1D(Form("hm0_before_%d", iter_cent),Form("hm0_before_%d", iter_cent), 100, 1.07, 1.17);
				fOutputList->Add(hm0_before[iter_cent]);				
			}
		}
	}else 
	{
		cout << "No such analysis_choice defined yet! Please choose either 'selection' or 'analysis'! Exiting... " << endl;
		exit(0);
	}
	
	ResEP1_true = new double[NITER_CENT];
	SubEvRes1 = new double[NITER_CENT];

	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		ResEP1_true[iter_cent] = 0.;
		SubEvRes1[iter_cent] = 0.;
	}
	
	
}

void MpdGlobalPolarization_RECO::ProcessEvent(MpdAnalysisEvent &event)
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

	ParticleMCProperties(event);

	// Find TPC track IDs 
	int idMax = 0;
	for (int j = 0; j < mKalmanTracks->GetEntriesFast(); j++) 
	{
		MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(j);
		idMax = TMath::Max(idMax,tr->GetTrackID());
	}
	lays = new int [idMax+1];
	kProb = new double [idMax+1];
	piProb = new double [idMax+1];
	pProb = new double [idMax+1];
	eProb = new double [idMax+1];
	id2dst = new int [idMax+1];

	CalculateLambdaAcceptance(event, idMax);

	// Collect "good" pions, kaons and protons		
	vector<int> vecPi, vecK, vecP;
	CollectParticles(vecPi, vecK, vecP);
	RecoEff(vecP, vecPi, 1);
	ApplyPid(newPid, vecP, vecPi);  	
	vector<MpdParticle*> vecL;
	vecL.clear();
	vLambdas.clear();
	BuildLambda(vecP, vecPi, vecL, phiEP_mc);
	nLamb = vecL.size();
	//cout << "Centrality_tpc = " << Centrality_tpc << "; nLamb_MC = " << nLamb_MC << "; vecL.size() = " << vecL.size() << "; vLambdas.size() = " << vLambdas.size() << endl;
	results->Fill(); 
	
	fillHistograms(event);

	for (int ipart = 0; ipart < nLamb; ipart++) delete vecL[ipart];

	delete [] lays;
	delete [] kProb;
	delete [] piProb;
	delete [] pProb;
	delete [] eProb;
	delete [] id2dst;

}

void MpdGlobalPolarization_RECO::Finish()
{
   cout << "Finish() ..." << endl;

}

bool MpdGlobalPolarization_RECO::selectEvent(MpdAnalysisEvent &event)
{
	static bool isInitialized = false;

	if (!isInitialized) 
	{
		if(MCFile == "")
		{
			cout << "Parameter MCFile must be set in the config file!" << endl;
			exit(0);
		}
		simMC = new TChain("mpdsim");
		simMC->AddFile(TString(MCFile)); //using the one from Request 25 (should be same geometry)
		simMC->SetName("mpdsim1");
		TFile fileMC(simMC->GetListOfFiles()->First()->GetTitle());
		fileMC.Get("FairGeoParSet");
		TClonesArray *tpcPoints = (TClonesArray*) fileMC.FindObjectAny("TpcPoint");
		simMC->SetBranchAddress("TpcPoint",&tpcPoints);
		TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
		secGeo = new TpcSectorGeoAZ();
		recoTpc = new MpdTpcKalmanFilter(*secGeo,"TPC Kalman filter");
		recoTpc->FillGeoScheme();
		isInitialized = true;
	}
	
	mhEvents->Fill(0.5); // Number of full events

	mMCTracks = event.fMCTrack;
	mKalmanTracks = event.fTPCKalmanTrack;
	mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();

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
	if (nTrMc <= 2*209) // Just nucleons of Bi+Bi --> Request30-PHSD
	{ 
		cout << "nTrMc clause - empty event!" << endl;
		cout << " nTrMc " << nTrMc << "; nMC " << nMC << endl;
		return false;
	}

	mhEvents->Fill(1.5); // Number of events with filled vertex

	if (!event.fVertex) // if even vertex not filled, skip event
	{ 
		return false;
	}
   
	vertex = (MpdVertex *)event.fVertex->First();
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
	
	if (Centrality_tpc < 0 || Centrality_tpc >= 100) // TPC centrality not defined
	{ 
		return false;
	}

	mhEvents->Fill(3.5); // Number of events, satisfying all the vertex criteria and having defined centrality

	mhVertex->Fill(mPrimaryVertex.Z());

	mhCentrality->Fill(Centrality_tpc);

	return true;
}

void MpdGlobalPolarization_RECO::ParticleMCProperties(MpdAnalysisEvent &event)
{
	vLambPtEtaY.clear();
	TVector3 genVert;
	event.fMCEventHeader->GetVertex(genVert); 
	nLamb_MC = 0;
	for (int j = 0; j < mMCTracks->GetEntriesFast(); j++) 
	{
		TVector3 mom; 
		MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(j);
		mcTr->GetMomentum(mom);	
		if (mcTr->GetPdgCode() == pdgCodeL0) 
		{
			// Check production vertex
			TVector3 pos;
			Double_t r = 0.0;
			if (mcTr->GetMotherId() >= 0) mcTr->GetStartVertex(pos);
			pos -= genVert;
			r = pos.Mag();
			if (r < 50.0) 
			{
				// Production vertex constraint 50 cm
				hLambFlag->Fill(0);
				hLambPTH->Fill(0);
				double pt = mom.Pt();
				if (pt < 0.5) hLambPTH->Fill(2);
				if (pt > 0.5 && pt < 1.0) hLambPTH->Fill(4);
				if (pt > 1.0 && pt < 1.5) hLambPTH->Fill(6);
				if (pt > 1.5 && pt < 2.0) hLambPTH->Fill(8);
				if (pt > 2.0) hLambPTH->Fill(10);					
				vLambPtEtaY.push_back(make_tuple(pt,mom.Eta(),mcTr->GetRapidity()));
			}
			nLamb_MC++;
		}
		
		if (mcTr->GetMotherId() < 0) continue;
		TVector3 pos;
		mcTr->GetStartVertex(pos);
		if (mom.Pt() < 0.001) continue;
		if (TMath::Abs(mom.Eta()) < 1.3) 
		{
			pdg = mcTr->GetPdgCode();
			moth = ((MpdMCTrack*) mMCTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
			if (moth == pdgCodeL0) 
			{
				hAngle->Fill(mom.Pt(),mom.Angle(pos)*TMath::RadToDeg());
			}
		}
	}
}

void MpdGlobalPolarization_RECO::CalculateLambdaAcceptance(MpdAnalysisEvent &event, int idMax)
{
	int *ids = new int [idMax+1];
	int *moths = new int [idMax+1];
	int *pdgs = new int [idMax+1];
	double *pt = new double [idMax+1];
	double *th = new double [idMax+1];
	double *rad = new double [idMax+1];
	double *dZ = new double [idMax+1];
	FairMCPoint **point = new FairMCPoint* [idMax+1];
	MpdTpcKalmanTrack **track = new MpdTpcKalmanTrack* [idMax+1];
	
	for (int j = 0; j <= idMax; j++) 
	{ 
		ids[j] = lays[j] = 0; 
		point[j] = 0x0;
		track[j] = 0x0;
		dZ[j] = 999999;
	} 
	
	// Get max. reached layer No.
	for (int j = 0; j < mKalmanTracks->GetEntriesFast(); j++) 
	{
		MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(j);
		int id = tr->GetTrackID();
		ids[id]++;		
		MpdKalmanHit *hit = (MpdKalmanHit*) tr->GetTrHits()->First();
		lays[id] = TMath::Max (hit->GetLayer(), lays[id]);
   
	}
	
	// Exclude "clones" (multiple loops)
	for (int j = 0; j < mKalmanTracks->GetEntriesFast(); j++) 
	{
		TVector3 mom; 
		MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(j);
		int id = tr->GetTrackID();
		if (track[id] == 0x0) track[id] = tr;		
		// MC track
		MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(id);
		mcTr->GetMomentum(mom);
		pt[id] = mom.Pt();
		th[id] = mom.Theta();
	}

	// Loop over DST tracks
	for (int j = 0; j < mMpdGlobalTracks->GetEntriesFast(); j++) 
	{
		MpdTrack *mpdTr = (MpdTrack*) mMpdGlobalTracks->UncheckedAt(j);
		
		int id = mpdTr->GetID();
		if (id > idMax || track[id] == 0x0) continue;
		if (ids[id] == 1) 
		{
			kProb[id] = mpdTr->GetPidProbKaon();
			piProb[id] = mpdTr->GetPidProbPion();
			pProb[id] = mpdTr->GetPidProbProton();
			eProb[id] = mpdTr->GetPidProbElectron();
			id2dst[id] = j;
		} else 
		{
			if (TMath::Abs(mpdTr->GetFirstPointZ()-track[id]->GetParam(1)) < dZ[id]) 
			{
				dZ[id] = TMath::Abs(mpdTr->GetFirstPointZ()-track[id]->GetParam(1));
				kProb[id] = mpdTr->GetPidProbKaon();
				piProb[id] = mpdTr->GetPidProbPion();
				pProb[id] = mpdTr->GetPidProbProton();
				eProb[id] = mpdTr->GetPidProbElectron();
				id2dst[id] = j;
			}
		}
		
	}

	// Lambda acceptance
	multimap<int,int> mapLamb, mapXi;
	for (int j = 0; j <= idMax; j++) 
	{
		TVector3 mom; 
		MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(j);
		mcTr->GetMomentum(mom);
		int mothID = mcTr->GetMotherId();
		if (mothID == -1 && lays[j] != 0) 
		{
			lays[j] = -lays[j]; // flag primary tracks
		}
		TVector3 pos;
		mcTr->GetStartVertex(pos);
		rad[j] = pos.Pt();
		moths[j] = 0;
		pdgs[j] = mcTr->GetPdgCode();
		if (mothID >= 0) {
			// Check lambda production vertex ( < 50 cm)
			MpdMCTrack* moth = (MpdMCTrack*) mMCTracks->UncheckedAt(mothID);
			moth->GetStartVertex(pos);
			if (pos.Mag() < 50.0) 
			{
				moths[j] = moth->GetPdgCode();
				if (moths[j] == pdgCodeL0 && (pdgs[j] == pdgCodePr || pdgs[j] == pdgCodeNeg)) mapLamb.insert(pair<int,int>(mothID,j));
			}
			if (moths[j] == pdgCodeXi && (pdgs[j] == pdgCodeL0 || pdgs[j] == pdgCodeNeg)) mapXi.insert(pair<int,int>(mothID,j));
		}
	}
	multimap<int,int>::iterator mit, mit1;
	pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;
	
	mit = mapLamb.begin();
	while (mit != mapLamb.end()) 
	{
		int mothID = mit->first;
		if (mapLamb.count(mothID) != 2) {mit = mapLamb.upper_bound(mothID); continue; } // only one decay particle
		ret = mapLamb.equal_range(mothID);
		int nppi[2] = {0}, nok = 0;
		int nok1 = 0, nok2 = 0, nok3 = 0;

		for (mit1 = ret.first; mit1 != ret.second; ++mit1) 
		{
			TVector3 mom; 
			MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(mit1->second);
			if (mcTr->GetPdgCode() == pdgCodePr) nppi[0] = 1; 
			else if (mcTr->GetPdgCode() == pdgCodeNeg) nppi[1] = 1;
			mcTr->GetMomentum(mom);
			if (mom.Pt() < 0.001) continue;
			if (TMath::Abs(mom.Eta()) < 1.3) ++nok;
			if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.05) ++nok1;
			if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.1) ++nok2;
			if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.2) ++nok3;
		}
		if (nppi[0] != 1 || nppi[1] != 1) { 
			// not p - p- decay   
			cout << " Wrong decay mode !!! " << endl; 
			mit = mapLamb.upper_bound(mothID);
			continue; 
		}
		if (nppi[0] == 1 && nppi[1] == 1) hLambFlag->Fill(1);
		if (nok == 2) hLambFlag->Fill(2); 
		if (nok1 == 2) hLambFlag->Fill(4); 
		if (nok2 == 2) hLambFlag->Fill(6); 
		if (nok3 == 2) hLambFlag->Fill(8); 

		// Check Xi-
		MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(mothID);
		int gmID = mcTr->GetMotherId();

		if (mapXi.find(gmID) != mapXi.end()) 
		{
			ret = mapXi.equal_range(gmID);
			for (mit1 = ret.first; mit1 != ret.second; ++mit1) 
			{
				TVector3 mom; 
				MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(mit1->second);
				if (mcTr->GetPdgCode() != pdgCodeNeg) continue;
				mcTr->GetMomentum(mom);
				if (mom.Pt() < 0.001) continue;
			}
		}
		mit = mapLamb.upper_bound(mothID);
	} // while (mit != mapLamb.end())

	delete [] ids;
	delete [] moths;
	delete [] pdgs;
	delete [] pt;
	delete [] th;
	delete [] point;
	delete [] rad;
	delete [] dZ;
	delete [] track;
}

void MpdGlobalPolarization_RECO::fillHistograms(MpdAnalysisEvent &event)
{
	NCentr->Fill(Centrality_tpc);
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		if(Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue; 
		
		ResEP1_true[iter_cent] += ResEP_mc;
		SubEvRes1[iter_cent] += ResEPSub_mc;		
	}
	
	for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
	{
		Resolution_EP1_true->SetBinContent(iter_cent+1,ResEP1_true[iter_cent]);
		Resolution_EP1_exp->SetBinContent(iter_cent+1,SubEvRes1[iter_cent]);
	}
	if(analysis_choice == "selection")
	{
		for(int i = 0; i < vLambdas.size(); i++) //cycle for reco Lambda
		{
			L0* lamb = (L0*) &vLambdas.at(i);
			hm0_before_full->Fill(lamb->massh);
			if(selection_choice == "omega2")
			{
				for(int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++)
				{	
					if(lamb->omega2 > omega_value[iter_sel]) 
						hm0_full[iter_sel]->Fill(lamb->massh);
				}	
				for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
				{
					if(Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue; 
					hm0_before[iter_cent]->Fill(lamb->massh);
					for(int iter_sel = 0; iter_sel < NITER_Selections; iter_sel++)
					{	
						if(lamb->omega2 > omega_value[iter_sel]) 
							hm0[iter_cent][iter_sel]->Fill(lamb->massh);
					}		
				}
			}
		}
		
	}else if(analysis_choice == "analysis")
	{
		for(int i = 0; i < vLambdas.size(); i++) //cycle for reco Lambda
		{
			L0* lamb = (L0*) &vLambdas.at(i);
			hm0_before_full->Fill(lamb->massh);
			if(selection_choice == "omega2")
			{
				if(lamb->omega2 > omega_value_full) 
					hm0_Full->Fill(lamb->massh); 
			}
			for(int iter_cent = 0; iter_cent < NITER_CENT; iter_cent++)
			{
				if(Centrality_tpc >= centrality_max[iter_cent] || Centrality_tpc < centrality_min[iter_cent]) continue; 

				if (lamb->origs[0] > 0)
				{							
					Lpolar[iter_cent]->Fill(lamb->polarhy);
					
					double phi_diff_hist = phiEP_mc - lamb->phi_star;
					double phi_diff_histRP = phiRP_mc - lamb->phi_star;
					double phi_diff_MC = phiRP_mc - lamb->phi_star_MC;
					if (phi_diff_hist < 0) phi_diff_hist = phi_diff_hist + 2.*pi;
					if (phi_diff_histRP < 0) phi_diff_histRP = phi_diff_histRP + 2.*pi;	
					if (phi_diff_MC < 0) phi_diff_MC = phi_diff_MC + 2.*pi;	
						
					PstarEP_hist[iter_cent]->Fill(phi_diff_hist); 
					PstarRP_hist[iter_cent]->Fill(phi_diff_histRP);
					PstarRP_hist_MC[iter_cent]->Fill(phi_diff_MC);				
				}
				if (lamb->origs[0] == 1) //true lambda (primary??? how now???)
				{	
					Lpolar_prim[iter_cent]->Fill(lamb->polarhy);
					
					float phi_diff_hist = phiEP_mc - lamb->phi_star;
					float phi_diff_histRP = phiRP_mc - lamb->phi_star;
					float phi_diff_MC = phiRP_mc - lamb->phi_star_MC;
					if (phi_diff_hist < 0) phi_diff_hist = phi_diff_hist + 2.*pi;
					if (phi_diff_histRP < 0) phi_diff_histRP = phi_diff_histRP + 2.*pi;	
					if (phi_diff_MC < 0) phi_diff_MC = phi_diff_MC + 2.*pi;	
						
					PstarEP_hist_prim[iter_cent]->Fill(phi_diff_hist); 
					PstarRP_hist_prim[iter_cent]->Fill(phi_diff_histRP);
					PstarRP_hist_MC_prim[iter_cent]->Fill(phi_diff_MC);
				}
				hm0_before[iter_cent]->Fill(lamb->massh);

				if(selection_choice == "omega2")
				{
					if(lamb->omega2 <= omega_value[iter_cent]) continue; 
				}else if(selection_choice == "chi")
				{
					if(lamb->chi2s[0] <= chi_pi_value[iter_cent]) continue; 
					if(lamb->chi2s[1] <= chi_p_value[iter_cent]) continue; 
					if(lamb->chi2h >= chi_pi_value[iter_cent]) continue; 
					if(lamb->path <= lambda_path_value[iter_cent]) continue; 
					if(lamb->angle >= lambda_angle_value[iter_cent]) continue;	
				}

				hm0_after[iter_cent]->Fill(lamb->massh);
				Dca_pion[iter_cent]->Fill(lamb->dcas[0]);
				Dca_proton[iter_cent]->Fill(lamb->dcas[1]);
				Chi_pion[iter_cent]->Fill(lamb->chi2s[0]);
				Chi_proton[iter_cent]->Fill(lamb->chi2s[1]);
				Dca_lambda[iter_cent]->Fill(lamb->dca);
				Chi_lambda[iter_cent]->Fill(lamb->c2pv);
				Dca_v0[iter_cent]->Fill(lamb->disth);
				Chi_v0[iter_cent]->Fill(lamb->chi2h);
				Path_hist[iter_cent]->Fill(lamb->path);
				Angle_hist[iter_cent]->Fill(lamb->angle);	
				double phi_diff = phiEP_mc - lamb->phi_star;
				if (phi_diff < 0) phi_diff = phi_diff + 2.*pi;
										
				for(int iter = 0; iter < NITER; iter++)
				{
					if(phi_diff < angle_max[iter] && phi_diff >= angle_min[iter])
					{
						hm0[iter_cent][iter]->Fill(lamb->massh);
					}
				}
			}			
		}
	}
}

// Collect "good" pions, kaons and protons		
void MpdGlobalPolarization_RECO::CollectParticles(vector<int> &vecPi, vector<int> &vecK, vector<int> &vecP)
{
	for (Int_t j = 0; j < mKalmanTracks->GetEntriesFast(); ++j) 
	{
		MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(j);
		if (tr->GetChi2() < -8) continue;
		Int_t id = tr->GetTrackID();
		MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(id);
		// !!!
		if (mcTr->GetMotherId() == 0 &&
		((MpdMCTrack*)mMCTracks->UncheckedAt(0))->GetPdgCode() == 1010010030) continue; // !!! decay product of artificial H3L (do i need it?)
		// !!! 
		if (mcTr->GetPdgCode() == pdgCodePr && tr->Charge() == 
		TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodePr)->Charge()/3)) vecP.push_back(j);
		else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == 
		TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeNeg)->Charge()/3)) vecPi.push_back(j); 	
		if (mcTr->GetPdgCode() == pdgCodePos && tr->Charge() == 1) hRecognitionFlag->Fill(1);
		else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) hRecognitionFlag->Fill(5);	
		
		if (tr->Charge() == 1 && pProb[id] > piProb[id] && pProb[id] > 0.25) 
		{
			// Fill if Proton
			if (mcTr->GetPdgCode() == pdgCodePr)
			{
				hRecognitionFlag->Fill(2); //true proton
				hProbTrueP->Fill(pProb[id],piProb[id]);
			}
			// Fill if not Proton
			if (mcTr->GetPdgCode() != pdgCodePr) hRecognitionFlag->Fill(3); //false proton
			hProbP->Fill(pProb[id],piProb[id]);
		}	
		else if (tr->Charge() == -1 && piProb[id] > pProb[id] && piProb[id] > kProb[id] && piProb[id] > eProb[id] && piProb[id] > 0.25) 
		{			
			hProbTrueK->Fill(kProb[id],piProb[id]);
			// Fill if Pion
			if (mcTr->GetPdgCode() == pdgCodeNeg)
			{
				hRecognitionFlag->Fill(6); // true pion
				hProbTruePi->Fill(pProb[id],piProb[id]);
			}
			// Fill if not Pion
			if (mcTr->GetPdgCode() != pdgCodeNeg) 
			{
				hRecognitionFlag->Fill(7); // false pion
				hPdg->Fill(mcTr->GetPdgCode());
				hProbPi->Fill(pProb[id],piProb[id]);
			}
		}	

	}
}
MpdHelix MpdGlobalPolarization_RECO::MakeHelix(const MpdTpcKalmanTrack *tr) 
{
	double r = tr->GetPosNew();
	double phi = tr->GetParam(0) / r;
	double x = r * TMath::Cos(phi);
	double y = r * TMath::Sin(phi);
	double dip = tr->GetParam(3);
	double cur = 0.3 * 0.01 * 5 / 10; // 5 kG
	cur *= TMath::Abs (tr->GetParam(4));
	TVector3 o(x, y, tr->GetParam(1));
	int h = (int) TMath::Sign(1.1,tr->GetParam(4));
	MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
	return helix;
}

MpdHelix MpdGlobalPolarization_RECO::MakeHelix(const MpdParticle *part) 
{
	double dip = TMath::PiOver2() - part->Theta();
	double cur = TMath::Abs (part->GetMeas(4));
	if (part->GetCharge() == 0) cur = numeric_limits<double>::epsilon();
	int h = (int) TMath::Sign(1.1,part->GetMeas(4));
	double phase = part->GetMeas(2) - TMath::PiOver2() * h;
	double x = part->GetXY(0);
	double y = part->GetXY(1);
	TVector3 o(x, y, part->GetMeas(1));
	MpdHelix helix(cur, dip, phase, o, h);
	return helix;
}  

void MpdGlobalPolarization_RECO::ApplyPid(MpdPid *pid, vector<int> &vecP, vector<int> &vecPi)
{
	map<int,set<int> > mapL, mapXi;

	int nP = vecP.size();
	int nPi = vecPi.size();

	for (int ip = 0; ip < nP; ip++) 
	{
		MpdTpcKalmanTrack *trP = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecP[ip]);
		if (trP->GetUniqueID() > 0) 
		{
			// Lambda decay product
      			int mid = trP->GetUniqueID();
      			if (mapL.find(mid) == mapL.end()) { set<int> aaa; mapL[mid] = aaa; }
      			mapL[mid].insert(vecP[ip]);
      			trP->SetUniqueID(0); // reset 
   		}
		if (trP->GetVertex().GetUniqueID() > 0) 
		{
		// Xi- decay product
		int mid = trP->GetVertex().GetUniqueID();
		if (mapXi.find(mid) == mapXi.end()) { set<int> aaa; mapXi[mid] = aaa; }
		mapXi[mid].insert(vecP[ip]);
		}
	}

	for (int ip = 0; ip < nPi; ip++) 
	{
		MpdTpcKalmanTrack *trP = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecPi[ip]);
		if (trP->GetUniqueID() > 0) 
		{
			// Lambda decay product
			int mid = trP->GetUniqueID();
			if (mapL.find(mid) == mapL.end()) { set<int> aaa; mapL[mid] = aaa; }
			mapL[mid].insert(vecPi[ip]);
			trP->SetUniqueID(0); // reset 
		}
		if (trP->GetVertex().GetUniqueID() > 0) 
		{
			// Xi- decay product
			int mid = trP->GetVertex().GetUniqueID();
			if (mapXi.find(mid) == mapXi.end()) { set<int> aaa; mapXi[mid] = aaa; }
			mapXi[mid].insert(vecPi[ip]);
		}
	}

	vecP.clear();
	vecPi.clear();

	int nTracks = mKalmanTracks->GetEntriesFast();
	for (int j = 0; j < nTracks; j++) 
	{
		MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(j);
		if (tr->GetChi2() < -8) continue;
		int id = tr->GetTrackID();
		MpdMCTrack* mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(id);
		int mothId = mcTr->GetMotherId();
		int uid = tr->GetVertex().GetUniqueID();
   
		MpdTrack* mpdTrack = (MpdTrack*) mMpdGlobalTracks->UncheckedAt(j);
		if (mpdTrack->GetID() != id) { cout << id << " " << mpdTrack->GetID() << endl; Fatal("ApplyPid"," Different ID"); }
		int ret = 0, charge = tr->Charge(), tofFlag = mpdTrack->GetTofFlag();
		double dedx = tr->GetDedx(), m2 = mpdTrack->GetTofMass2();
		
		if (tofFlag == 2 || tofFlag == 6)          // dE/dx+TOF
		ret = pid->FillProbs(tr->Momentum(), dedx, m2, charge);
		if (ret == 0) ret = pid->FillProbs(tr->Momentum(), dedx, charge);
		if (ret == 0) // No PID
		{
			if (mcTr->GetPdgCode() == pdgCodeNeg) hPIDflag->Fill(2.1); // lost pion
			if (mcTr->GetPdgCode() == pdgCodePr) hPIDflag->Fill(6.1); // lost proton
			continue;
		}
		double piThr = -0.75;
		double probThr = -0.60;
		if (pdgCodeL0 * tr->Charge() < 0) 
		{
			double prob = pid->GetProbPi();
			if (prob > piThr && prob > pid->GetProbKa() && prob > pid->GetProbEl() && prob > pid->GetProbPr() && prob > pid->GetProbMu()) 
			{
				// "pion"
				if (mcTr->GetPdgCode() == pdgCodeNeg) hPIDflag->Fill(0.1); // correct pion
				else if (mcTr->GetPdgCode() != pdgCodeNeg) hPIDflag->Fill(1.1); // false pion
				//
				if (mapL.find(mothId+1) != mapL.end() && mapL[mothId+1].find(j) != mapL[mothId+1].end())
				mapL[mothId+1].erase(j);
				if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end())
				mapXi[uid].erase(j);
				//
				double chi2 = TMath::Min (tr->GetChi2Vertex(),999.);
				if (chi2 < gC2pi) continue;
				vecPi.push_back(j);
			} else if (mcTr->GetPdgCode() == pdgCodeNeg) hPIDflag->Fill(2.1); // lost pion
		} else 
		{
			double prob = pid->GetProbPr();
			if (prob > probThr && prob > pid->GetProbKa() && prob > pid->GetProbPi() && prob > pid->GetProbDe()) 
			{
				// "proton"
				if (mcTr->GetPdgCode() == pdgCodePr) hPIDflag->Fill(4.1); // correct proton
				else if (mcTr->GetPdgCode() != pdgCodePr) hPIDflag->Fill(5.1); // false proton
				//
				if (mapL.find(mothId+1) != mapL.end() && mapL[mothId+1].find(j) != mapL[mothId+1].end())
				mapL[mothId+1].erase(j);
				if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end())
				mapXi[uid].erase(j);
				//
				MpdTpcKalmanTrack trCor = *tr;
				trCor.SetDirection(MpdKalmanTrack::kInward);
				recoTpc->Refit(&trCor, 0.93827, 1); // refit
				//~ ShowTrackStats(tr);
				MpdParticle prot(trCor, 0);
				prot.SetPdg(pdgCodePr);
				prot.SetMass();
				double chi2 = TMath::Min (prot.Chi2Vertex(vertex),999.);
				if (chi2 < gC2p) continue;
				vecP.push_back(j);
			} else if (mcTr->GetPdgCode() == pdgCodePr) hPIDflag->Fill(6.1); // lost proton
		}
	}    
	//
	int nLok = 0;
	int nXiok = 0;
	for (map<int,set<int> >::iterator mit = mapL.begin(); mit != mapL.end(); mit++) 
	{
		if (mit->second.size() == 0) nLok++;
	}
	for (map<int,set<int> >::iterator mit = mapXi.begin(); mit != mapXi.end(); mit++) 
	{
		if (mit->second.size() == 0) nXiok++;
	}
}

void MpdGlobalPolarization_RECO::RecoEff(vector<int> &vecP, vector<int> &vecPi, bool use_pid)
{
	//TODO: need to carefully look through and rewrite anything which I don't like/need

	int nPi = vecPi.size();
	int nP = vecP.size();

	for (int ip = nP - 1; ip >= 0; ip--) // Proton
	{
		MpdTpcKalmanTrack *trP = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecP[ip]);
		MpdMCTrack *mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(trP->GetTrackID());
		int mothId = mcTr->GetMotherId();
		if (mothId < 0) continue;
		MpdMCTrack *moth = (MpdMCTrack*) mMCTracks->UncheckedAt(mothId);  
		if (moth->GetPdgCode() == pdgCodeL0) 
		{
			Int_t mp = mothId;
			// Proton from Lambda
			for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) // Pion
			{	
				MpdTpcKalmanTrack *trPi = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecPi[jpi]);
				MpdMCTrack *mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(trPi->GetTrackID());
				int mothId = mcTr->GetMotherId();
				if (mothId < 0) continue;
				MpdMCTrack *moth = (MpdMCTrack*) mMCTracks->UncheckedAt(mothId);  
				if (moth->GetPdgCode() == pdgCodeL0 && mp == mothId) 
				{
					//AZ - flag decay tracks to check PID influence later
					trP->SetUniqueID(mothId+1);
					trPi->SetUniqueID(mothId+1);
					//
					int gmId = moth->GetMotherId();
					if (gmId >= 0) 
					{
						MpdMCTrack *gmoth = (MpdMCTrack*) mMCTracks->UncheckedAt(gmId);
						if (gmoth->GetPdgCode() == pdgCodeXi) 
						{
							for (int kpi = nPi - 1; kpi >= 0; kpi--) 
							{
								// Pion
								MpdTpcKalmanTrack *trK = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecPi[kpi]);
								MpdMCTrack *mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(trK->GetTrackID());
								int mothId = mcTr->GetMotherId();
								if (mothId < 0) continue;
								MpdMCTrack *moth = (MpdMCTrack*) mMCTracks->UncheckedAt(mothId);  
								if (moth->GetPdgCode() == pdgCodeXi && gmId == mothId) 
								{
									//AZ - flag decay tracks to check PID influence later
									trK->GetVertex().SetUniqueID(mothId+1); 
									trP->GetVertex().SetUniqueID(mothId+1);
									trPi->GetVertex().SetUniqueID(mothId+1);
									break;
								}
							}
						} // if (gmoth->GetPdgCode() == pdgCodeXi)
					}
					break;
				}
			}
		}
	}

	if (use_pid) return; // skip the rest if PID is used

	for (int ip = nP - 1; ip >= 0; ip--) // Proton
	{
		MpdTpcKalmanTrack *trP = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecP[ip]);
		MpdTpcKalmanTrack trCor = *trP;
		trCor.SetDirection(MpdKalmanTrack::kInward);
		recoTpc->Refit(&trCor, 0.93827, 1); // refit
		MpdParticle prot(trCor, vecP[ip]);
		prot.SetPdg(pdgCodePr);
		prot.SetMass();

		double chi2 = TMath::Min (prot.Chi2Vertex(vertex),999.);
		if (chi2 < gC2p) vecP.erase(vecP.begin()+ip);
	}

	if (nP) 
	{
		for (Int_t jpi = nPi - 1; jpi >= 0; jpi--) 
		{
			// Pion
			MpdTpcKalmanTrack *trPi = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecPi[jpi]);
			MpdTpcKalmanTrack trCor = *trPi;
 
			double chi2 = TMath::Min (trPi->GetChi2Vertex(),999.);
			if (chi2 < gC2pi) vecPi.erase(vecPi.begin()+jpi);
		}
	}
}

void MpdGlobalPolarization_RECO::BuildLambda(vector<int> &vecP, vector<int> &vecPi, vector<MpdParticle*> &vecL, double &phiEP) 
{
	int nPi = vecPi.size();
	int nP = vecP.size();
	vector<MpdParticle*> vPart;
	vecL1.clear(); 

	for (int ip = 0; ip < nP; ip++) // Proton
	{
		MpdTpcKalmanTrack *trP = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecP[ip]); 
		MpdMCTrack *mcTr = (MpdMCTrack*) mMCTracks->UncheckedAt(trP->GetTrackID()); 
		TVector3 mcmom1;
		mcTr->GetMomentum(mcmom1);
		mcps[1] = mcmom1.Mag();
		mcphis[1] = mcmom1.Phi();
		mcthetas[1] = mcmom1.Theta();
		Int_t mothId = mcTr->GetMotherId();
		MpdTpcKalmanTrack trCor = *trP;
		trCor.SetDirection(MpdKalmanTrack::kInward);
		recoTpc->Refit(&trCor, 0.93827, 1); // refit
		MpdParticle prot(trCor, vecP[ip]);
		prot.SetPdg(pdgCodePr);
		prot.SetMass();
		qs[1] = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodePr)->Charge()/3); //again the same error
		etas[1] = trP->Momentum3().Eta();
		
		chi2s[1] = TMath::Min (prot.Chi2Vertex(vertex),9999.);
		layMx[1] = TMath::Abs (lays[trP->GetTrackID()]);
		MpdHelix helix = MakeHelix(trP);
		//Get 3-D DCA to primary vertex
		TVector3 pca;
		double s = helix.pathLength(mPrimaryVertex);
		pca = helix.at(s);
		pca -= mPrimaryVertex;
		dcas[1] = pca.Mag();
		origs[1] = 0;
		if (mothId >= 0 && ((MpdMCTrack*) mMCTracks->UncheckedAt(mothId))->GetPdgCode() == pdgCodeL0)
			origs[1] = -1; // from lambda

		for (Int_t jpi = 0; jpi < nPi; ++jpi) // Pion
		{		
			MpdTpcKalmanTrack *trPi = (MpdTpcKalmanTrack*) mKalmanTracks->UncheckedAt(vecPi[jpi]);
			MpdMCTrack *mcTr1 = (MpdMCTrack*) mMCTracks->UncheckedAt(trPi->GetTrackID());
			TVector3 mcmom2;
			mcTr1->GetMomentum(mcmom2);
			mcps[0] = mcmom2.Mag();
			mcphis[0] = mcmom2.Phi();
			mcthetas[0] = mcmom2.Theta();
			int mothId1 = mcTr1->GetMotherId();
			origs[0] = 0;
			if (mothId1 >= 0 && ((MpdMCTrack*) mMCTracks->UncheckedAt(mothId1))->GetPdgCode() == pdgCodeL0)
				origs[0] = -1; // from lambda
			MpdParticle *pion = new MpdParticle(*trPi, vecPi[jpi]);
			pion->SetPdg(pdgCodeNeg);
			pion->SetMass();

			vPart.clear();
			vPart.push_back(new MpdParticle(prot));
			vPart.push_back(pion);


			MpdParticle lambPart;
			Double_t chi2 = lambPart.BuildMother(vPart);
			TVector3 v0(lambPart.Getx()(0,0), lambPart.Getx()(1,0), lambPart.Getx()(2,0));
			v0 -= mPrimaryVertex;
			Double_t decay = v0.Mag();
			path = TMath::Sign (decay, v0*lambPart.Momentum3());

			if (chi2 >= 0 && chi2 < gC2L && path > gPathL) 
			{
				if (origs[1] > 0) origs[1] = -1;
				MpdMCTrack *moth = NULL;
				hMassL->Fill(lambPart.GetMass());
				if (mothId != mothId1 || mothId < 0) 
				{
					hMassLbkg->Fill(lambPart.GetMass());
				}else 
				{
					if (origs[0] == -1) 
					{
						hMassLsig->Fill(lambPart.GetMass());
						origs[0] = origs[1] = 1;
						moth = (MpdMCTrack*) mMCTracks->UncheckedAt(mothId);
					}
					else hMassLbkg->Fill(lambPart.GetMass());
				}

				// Filling results branch:
	
				qs[0] = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeNeg)->Charge()/3); // pion
				etas[0] = trPi->Momentum3().Eta();
	
	
				chi2s[0] = TMath::Min (pion->Chi2Vertex(vertex),9999.);
				layMx[0] = TMath::Abs (lays[trPi->GetTrackID()]);
				MpdHelix helix1 = MakeHelix(trPi);
				// Get 3-D DCA to primary vertex
				s = helix1.pathLength(mPrimaryVertex);
				pca = helix1.at(s);
				pca -= mPrimaryVertex;
				dcas[0] = pca.Mag();

				massh = lambPart.GetMass();
				chi2h = chi2;
				angle = v0.Angle(lambPart.Momentum3());
				pth = lambPart.Pt(); // reconstructed
				ph = lambPart.Momentum(); // reconstructed
				phi_Lam = lambPart.Phi(); // reconstructed 
				if (pth > 0.001) 
					etah = lambPart.Momentum3().Eta(); 
				else 
					etah = TMath::Sign(100.,lambPart.Momentum3().Z()); 
				pair<Double_t,Double_t> paths = helix.pathLengths(helix1);
				TVector3 p1 = helix.at(paths.first);
				TVector3 p2 = helix1.at(paths.second);
				p1 -= p2;
				disth = p1.Mag(); // closest distance between daughters
				// Get 3-D DCA of lambda to primary vertex
				MpdHelix helix2 = MakeHelix(&lambPart);
				s = helix2.pathLength(mPrimaryVertex);
				pca = helix2.at(s);
				pca -= mPrimaryVertex;
				dca = pca.Mag();
				c2pv = TMath::Min (lambPart.Chi2Vertex(vertex),9999.);
				omega1 = dcas[0] * dcas[1] / (dca * dca + disth * disth);
				omega2 = TMath::Sqrt (chi2s[0] * chi2s[1]) / (c2pv + chi2h);

				dstNo[0] = vecPi[jpi]; // pion index
				dstNo[1] = vecP[ip]; // proton index

				
				if (lambPart.GetMass() >= 1.10518 && lambPart.GetMass() <= 1.12668) // lambda mass +- 5*2.15 MeV
				{ 
					vecL.push_back(new MpdParticle(lambPart));
					vector<Double_t> lambPars(6);
					lambPars[0] = disth;
					lambPars[1] = angle;
					for (Int_t jl = 0; jl < 2; ++jl) {
						lambPars[jl+2] = chi2s[jl];
						lambPars[jl+4] = dcas[jl];
					}
					vecL1.push_back(lambPars);
				}
	
				if (origs[0] == 1) // True lambda
				{
					lambPart.SetMass(1.11568); // set true mass
					yh = lambPart.Rapidity();
					// Check mother of lambda
					Int_t gMothId = moth->GetMotherId();
					if (gMothId >= 0) origs[0] = origs[1] = 2; // secondary lambda
				}


				polarhx = 0.0;
				polarhy = 0.0;
				polarhz = 0.0;
				if (origs[0] > 0) 
				{
					double weight_pol = ((MpdMCTrack*)mMCTracks->UncheckedAt(mothId))->GetWeight();
					polarhx = ((MpdMCTrack*)mMCTracks->UncheckedAt(mothId))->GetPolar(0);
					polarhy = ((MpdMCTrack*)mMCTracks->UncheckedAt(mothId))->GetPolar(1);
					polarhz = ((MpdMCTrack*)mMCTracks->UncheckedAt(mothId))->GetPolar(2);
					TVector3 polar_changed(polarhx, polarhy, polarhz);
					if (phiEP != 0.) 
						polar_changed.RotateZ(-phiEP);
					polarhx = weight_pol*polar_changed.X();
					polarhy = weight_pol*polar_changed.Y();
					polarhz = weight_pol*polar_changed.Z();
					Lpolar_full->Fill(polarhy);
				}

				FindPolarAngle (lambPart, vPart);
				L0 l0(massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas, mcthetas, thetas, mcphis, phis, mcps, ps, pts, chi2s, dcas, dca, c2pv, omega1, omega2, cosA, cosAmc, polarhx, polarhy, polarhz, phi_star, phi_star_MC, phi_Lam, origs, qs, layMx, evNo);
				vLambdas.push_back(l0);
			} // if (chi2 >= 0 && chi2 < gC2L...

			Int_t nPart = vPart.size();
			for (Int_t ipart = 0; ipart < nPart; ++ipart) delete vPart[ipart];
		} // for (Int_t jpi = 0; jpi < nPi;
	} // for (Int_t ip = 0; ip < nP;

}

void MpdGlobalPolarization_RECO::FindPolarAngle(MpdParticle &lamb, vector<MpdParticle*> &vPart)
{
	// Compute decay proton angle w.r.t. lambda decay plane

	TVector3 vPr, vPi, vLamb;
	TLorentzVector prLor, piLor, lambLor; 

	// True (exact) parameters
	vPr.SetMagThetaPhi(mcps[1], mcthetas[1], mcphis[1]);
	vPi.SetMagThetaPhi(mcps[0], mcthetas[0], mcphis[0]);
      
	prLor.SetVectM(vPr, 0.938272);
	piLor.SetVectM(vPi, 0.139570);

	lambLor = prLor + piLor;

	TVector3 boostV_MC;
	boostV_MC = lambLor.BoostVector();
	boostV_MC *= -1;
	
	prLor.Boost(boostV_MC);
	vPr = prLor.Vect();

	cosAmc = vPr.CosTheta();
	phi_star_MC = vPr.Phi();

	// Reco parameters
	vLamb.SetMagThetaPhi(lamb.Momentum(), lamb.Theta(), lamb.Phi());
	lambLor.SetVectM(vLamb, 1.11568);

	MpdParticle *prot = vPart[0];
	vPr.SetMagThetaPhi(prot->Momentum(), prot->Theta(), prot->Phi());
	
	prLor.SetVectM(vPr, 0.938272);
	
	TVector3 boostV;
	boostV = lambLor.BoostVector();
	boostV *= -1;
	  	
  	prLor.Boost(boostV);
	vPr = prLor.Vect();
	
	cosA = vPr.CosTheta();
	//calculating the azimuthal angle of proton in the lambda frame (phi*):
	phi_star = vPr.Phi();
	
}
double* MpdGlobalPolarization_RECO::init_double_array (const int n, const double fmt...)
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

int* MpdGlobalPolarization_RECO::init_int_array (const int n, const int fmt...)
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
