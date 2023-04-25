//Main analyzer for the study of global polarization of Lambda using PHSD data
//Collection of Lambda hyperons (using PID), multiplicity in TPC for centrality classes, EP angles and resolutions (1st order) using FHCal, global polarization values extracted from MCTracks 
//compile with: make Cluster_analyzer_Lambda, start with the corresponding sge script phsd_analysis_Lambda_v1.sge
//e.g. ./Cluster_analyzer_Lambda -gen_choice PHSD -inname $INFILE -outname $OUTFILE
//only for Lambda hyperons!

///check whether there is anything missing/obsolete!
///generalize the code as much as possible!
///recheck what you need in Stuct_L0_v1.h

// MPD includes: 
#include "TpcPoint.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdHelix.h"
#include "MpdVertex.h"
#include "MpdParticle.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdMotherFitterPart.h"
#include "MpdTrackFinderIts5spd.h"
#include "MpdMCEventHeader.h"
#include "MpdZdcDigi.h"
#include "MpdPid.h"
#include "MpdMCTrack.h"

// CBM includes:
#include "FairMCPoint.h"
#include "FairRunAna.h"

// ROOT includes:
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <Riostream.h>

#include <set>
#include <map>
#include <tuple>
#include <vector>

#include <iostream>

using namespace std;

#include "TStopwatch.h"
// header, which includes all three MC projections of polarization:
#include <Struct_L0_v1.h> 

// Global parameters:
#define _N_MODULES_TOTAL 90
#define _N_ARM 2 
//pdg numbers 
Int_t pdgCodePr= 2212; // proton
Int_t pdgCodeNeg= -211; // pi-
Int_t pdgCodePos= 211; // pi+
Int_t pdgCodeL0 = 3122; // lambda (1.11568)
Int_t pdgCodeXi = 3312; // Xi- (1.3213 GeV/c)
Int_t pdgCodeKm = -321; // K-

Int_t *lays = 0x0;
MpdHelix trC(0,0,0,TVector3(0,0,0),0);
MpdVertex *mpdVert;
TVector3 vtxN, momN, primVert;
TClonesArray *itsTracks, *mcTracks, *mpdTracks;

TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params 

//defining the parameters for the analysis (should check this at some point):
const Double_t gC2p      = 3.; //4.; //9.;           //chi2 of p to PV
const Double_t gC2pi     = 5.; //5.; //11.;          //chi2 of pion to PV
const Double_t gC2Lpv    = 0.;           //chi2 of L to PV
const Double_t gC2Xipv   = 0.;           //chi2 of Xi to PV

const Double_t gDCAp     = 0.;           //cm - DCA of p to PV
const Double_t gDCApi    = 0.;           //cm - DCA of pion to PV
const Double_t gDCAL     = 0.;           //cm - DCA of L to PV
const Double_t gDCAXi    = 0.;           //cm - DCA of Xi to PV

const Double_t gDistL    = 9999.;        //cm - DCA between pion & p in V0
const Double_t gDistXi   = 9999.;        //cm - DCA between pion & L in Xi
//const Double_t gPathL    = 2.0;          //cm - path to Lambda decay
const Double_t gPathL    = 0.0;          //cm - path to Lambda decay
const Double_t gPathXi   = 0.;           //cm - path to Xi decay
const Double_t gC2L      = 25.; //9999.;  //chi2 between pion & p in V0
const Double_t gC2Xi     = 25.; //16.;//9999.;  //chi2 between pion & L in Xi

const Double_t gDcaL0 = 0.; //0.15; 
const Double_t gDcaK = 0.; //0.3;
const Double_t gChi2K = 0.; //100; //50;
const Double_t gDecayOm = 0.; //1.0;
const Double_t gDistLK = 9999.; //1.15; //0.2;
const Double_t gDcaOm = 9999.; //0.1; //0.15;

//#define ITS
#ifdef ITS
typedef MpdItsKalmanTrack AzTrack;
#else
typedef MpdTpcKalmanTrack AzTrack;
#endif
MpdTpcKalmanFilter* recoTpc = NULL;
MpdTrackFinderIts5spd* recoIts = NULL;

//define the functions:
TChain* Chain(Int_t nFiles, TString firstFile);
TChain* ChainFile(Int_t nFiles, TString fileNameList, Int_t skipLines);
void RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi, Int_t pid = 0);
void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL, Float_t &phiEP_mc);
void ApplyPid(MpdPid *pid, vector<Int_t> &vecP, vector<Int_t> &vecPi);
MpdHelix MakeHelix(const MpdKalmanTrack *tr);
MpdHelix MakeHelix(const MpdParticle *part);
Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu);
void FindPolarAngle(MpdParticle &lamb, vector<MpdParticle*> &vPart);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

Double_t GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID);
Double_t GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy);
Double_t GetPsiFullZdc(Float_t* zdc_energy, Int_t n);
void GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy);
Double_t* GetAngles();

//define the variables for the tree:
Float_t massh, pth, ph, etah, yh, chi2h, disth, path, c2pv, d2pv, masshL, chi2hL, disthL, pathL, angL;
Float_t etas[2], mcps[2], ps[2], pts[2], chi2s[2], dcas[2], probs[2], chi2sL[2], dcasL[2],  angle, masshL1, chi2hL1;
Float_t b0, mcthetas[2], thetas[2], mcphis[2], phis[2];
Float_t dca, omega1, omega2, omegaL[3], cosA, cosAmc, polarhx, polarhy, polarhz, phi_star, phi_star_MC, phi_Lam;
Double_t *kProb, *piProb, *pProb, *eProb;
Int_t evNo, origs[2], qs[2], dstNo[2], layMx[2], ntr13, n_tracks_mpd, n_tracks_mpd_DCA;
	
Int_t *id2dst;
vector<vector<Double_t> > vecL1;
TRandom *RNG = new TRandom();
Float_t phiEP_mc, phiEP_mc_or, phiEP_1_ZDC, phiEP_2_ZDC;
Float_t ResEP_1_ZDC, ResEP_2_ZDC, ResEPSub_1_ZDC, ResEPSub_2_ZDC;
Float_t ZDC_energy_mpd[_N_MODULES_TOTAL];
Int_t nLamb, nLamb_MC;

TClonesArray *ZDCHits;
MpdZdcDigi* ZDCHit;

std::vector<L0> vLambdas;
std::vector<L0> *vvvL = &vLambdas;
L0 aaaa;
std::vector<tuple<float,float,float> > vLambPtEtaY;
std::vector<tuple<float,float,float> > *vvvLpt = &vLambPtEtaY;

//function to check input file:
bool checkfile (char* filename) 
{
	TFile *file_for_check = TFile::Open(filename);
	// if exist
    if (!file_for_check)
    {
        cout << "no specified file: " << filename << endl;
        return false;
    } 
    // if zombie
	if (file_for_check->IsZombie())
	{
        cout << "file: " << filename << " is a Zombie!" << endl;
        return false;
    }
	file_for_check->Close();
	return true;
}
//starting the main Analyzer:
int main(int argc, char** argv)
{
	TStopwatch timer;
    timer.Start();
    
	Int_t n1 = 0;
	Int_t n2 = 0;
	Int_t skipFiles = 0;
	Int_t iset = 1;
	
	Int_t _N_Hits = 16; // default value
	Double_t cut_pt = 0.15; // default value (GeV/c)
	Double_t cut_eta = 0.5; // default value
	Double_t dca_cut = 0.5; // default: 0.5 cm
	
	TChain chain("mpdsim");
	TString outname = ""; 
	TString generator = ""; 
	// generating dictionaries for vector classes (for special structures it needs to be done separately)
	
	// put the definition of the structure classes L0 and Xi, as well as the linking in a separate header file (e.g. Struct_L0.h), then:
	// rootcling -v4 -f mydict.cxx -rmf libmydict.rootmap -rml libmydict.so Struct_L0_v1.h
	// g++ -shared -o libmydict.so mydict.cxx `root-config --cflags --libs` -fPIC
	
//	possible choices for the analysis:
//	cout << "Usage: ./Cluster_analyzer -events 100 -gen_choice PHSD(URQMD) -inname path-to-infile -outname path-to-outfile << endl;
	
	for (int i = 0 ; i < argc ; i++) 
	{
		string tmp(argv[i]);
		if (tmp.compare("-events") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				n2 = number;
			}
		}
		if (tmp.compare("-cut_pt") == 0 && i+1 < argc) 
		{
			float number = atof(argv[i+1]);
			if (number > 0) 
			{
				cut_pt = number;
			}
		}
		if (tmp.compare("-cut_eta") == 0 && i+1 < argc) 
		{
			float number = atof(argv[i+1]);
			if (number > 0) 
			{
				cut_eta = number;
			}
		}
		if (tmp.compare("-_N_Hits") == 0 && i+1 < argc) 
		{
			int number = atoi(argv[i+1]);
			if (number > 0) 
			{
				_N_Hits = number;
			}
		}		
		if (tmp.compare("-gen_choice") == 0 && i+1 < argc) 
		{
			generator = argv[i+1];
			cout << "generator: " << argv[i+1] << endl;
		}
		if (tmp.compare("-inname") == 0 && i+1 < argc) 
		{
			int count = 0;
			while(i+count+1 < argc && argv[i+count+1][0] != '-')
			{
				//add the file to the chain after checking:
				if (checkfile(argv[i+count+1]))
				{
					cout << "adding to the chain: " << argv[i+count+1] << endl;
					chain.Add(argv[i+count+1]);
				}
				count++;
			}
			cout << "number of files = " << count << endl;
		}
		if (tmp.compare("-outname") == 0 && i+1 < argc) 
		{
			outname = argv[i+1];
		}
	}
	
	if(outname.Length() == 0)
	{
		cerr << "!! Please provide the path for the output file !!" << endl;
		return 1;
	}
	
	cout << "events chosen = " << n2 << "; N_Hits = " << _N_Hits << "; cut_pt = " << cut_pt << "; cut_eta = " << cut_eta << endl;
	
	TChain *simMC = new TChain("mpdsim");
	//always check which one do you have there!
	//simMC->AddFile("/eos/nica/mpd/users/nazarova/extra/evetest-1.root"); //have the updated one there
	simMC->AddFile("/eos/nica/mpd/sim/data/mc/urqmd-BiBi-09.2GeV-mb-eos0-500-0.dst.root"); //using the one from Request 25 (should be same geometry)
	
	simMC->SetName("mpdsim1");
	
	cout << "entries = " << chain.GetEntries() << endl;
	
	//return(0); //for debugging
	
	TChain *simITS = (TChain*) gROOT->FindObject("mpdsim");
	TFile fileITS(simITS->GetListOfFiles()->First()->GetTitle());
	TFile fileMC(simMC->GetListOfFiles()->First()->GetTitle());
	fileMC.Get("FairGeoParSet");
	TClonesArray *vtxs = (TClonesArray*) fileITS.FindObjectAny("Vertex");
	simITS->SetBranchAddress("Vertex",&vtxs);
	TBranch *vtxB = simITS->GetBranch("Vertex");
	
	itsTracks = (TClonesArray*) fileITS.FindObjectAny("TpcKalmanTrack");
	simITS->SetBranchAddress("TpcKalmanTrack",&itsTracks);
	TBranch *itsRecoB = simITS->GetBranch("TpcKalmanTrack");
	
	MpdEvent *event = 0x0;
	simITS->SetBranchAddress("MPDEvent.", &event);
	MpdMCEventHeader *mcHeader = 0x0;
	simITS->SetBranchAddress("MCEventHeader.", &mcHeader); 
	
	TClonesArray *tpcPoints = (TClonesArray*) fileMC.FindObjectAny("TpcPoint");
	simMC->SetBranchAddress("TpcPoint",&tpcPoints);
	TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
	
	mcTracks = (TClonesArray*) fileITS.FindObjectAny("MCTrack");
	simITS->SetBranchAddress("MCTrack",&mcTracks);
	TBranch *mcBranch = simITS->GetBranch("MCTrack");
	
	ZDCHits = 0;
	simITS->SetBranchAddress("ZdcDigi",&ZDCHits);
	
	
	//output file:
	TFile out(outname,"recreate");
	
	FairRunAna ana;
	MpdKalmanFilter::Instance("KF")->Init();
	recoTpc = new MpdTpcKalmanFilter("TPC Kalman filter");
	recoTpc->SetSectorGeo(MpdTpcSectorGeo::Instance());
	recoTpc->FillGeoScheme();
	
	Double_t sigM = 4.0, sigE = 4.0, energy = 9.0, coef = 1.0; // n-sigma bands for PID selection
	//TString generator = "PHSD", tracking = "CF"; // - old
	//TString tracking = "CF"; // - old
	TString tracking = "CFHM";
	cout<<"---------- BEFORE PID ----------"<<endl;
	MpdPid *newPid = new MpdPid(sigM, sigE, energy, coef, generator, tracking, "pikapr");
	
	//whatever histograms you want to have in the output file
	TH1D *hLambFlag = new TH1D("hLambFlag","Flags for lambda",12,0,12);
	TH1D *hRecognitionFlag = new TH1D("hRecognitionFlag","Flags for Recognition",10,0,10);
	TH1D *hLambPTH = new TH1D("hLambPTH","Flags for lambdaPTH",12,0,12);
	
	TH1D *hMassL = new TH1D("hMassL","Lambda mass",50,1.070,1.170);
	TH1D *hMassLsig = new TH1D("hMassLsig","Lambda mass (signal)",50,1.070,1.170);
	TH1D *hMassLbkg = new TH1D("hMassLbkg","Lambda mass (bkg.)",50,1.070,1.170);
	
	TH1D *hPdg = new TH1D("hPdg","PdgCode if is not Pion",1000,-2500,2500);
	new TH1D("hPIDflag","PID flags",12,0,12);
	TH2D *hAngle = new TH2D("hAngle","Acollinearity angle vs Pt",50,0,1,60,0,90);
	TH1D *hPolarL = new TH1D("hPolarL","Lambda polarization",100,-1.0,1.0);

	TH2D *hProbTrueP = new TH2D("hProbTrueP","Probability for true Protons",50,0,1,50,0,1.1);
	TH2D *hProbP = new TH2D("hProbfalseP","Probability for Pions and identification Protons",50,0,1.1,50,0,1.1);
	TH2D *hProbTruePi = new TH2D("hProbTruePi","Probability for true Pions",50,0,1.1,50,0,1.1);
	TH2D *hProbPi = new TH2D("hProbfalsePi","Probability for Protons and identification Pions",50,0,1.1,50,0,1.1);
	TH2D *hProbTrueK = new TH2D("hProbTrueK","Probability for true Kaons",50,0,1.1,50,0,1.1);
	TH2D *hVertex = new TH2D("hVertex","Distribution of Event Vertex",100,-100,100,100,-2,2);
	
	Double_t pmom, eta1, dpp, rorig, ptt;
	Int_t prim, idtr, np, moth, pdg;
	
	TTree *tree = new TTree("event","Event");
	tree->Branch("b0",&b0,"b0/F"); //impact parameter
	tree->Branch("ntr13",&ntr13,"ntr13/I"); // number of tracks selected for analysis
	tree->Branch("n_tracks_mpd",&n_tracks_mpd,"n_tracks_mpd/I"); // number of tracks (multiplicity) in TPC for centrality determination
	tree->Branch("n_tracks_mpd_DCA",&n_tracks_mpd_DCA,"n_tracks_mpd_DCA/I"); // number of tracks (multiplicity) in TPC for centrality determination, using DCA cut
	tree->Branch("phiEP_mc_or",&phiEP_mc_or,"phiEP_mc_or/F"); // Reaction plane angle from MC
	tree->Branch("phiEP_mc",&phiEP_mc,"phiEP_mc/F"); // Reaction plane angle from MC
	tree->Branch("phiEP_1_ZDC",&phiEP_1_ZDC,"phiEP_1_ZDC/F"); // 1st-order Event plane angle (ZDC)
	tree->Branch("phiEP_2_ZDC",&phiEP_2_ZDC,"phiEP_2_ZDC/F"); // 2st-order Event plane angle (ZDC)
	tree->Branch("ResEP_1_ZDC",&ResEP_1_ZDC,"ResEP_1_ZDC/F"); // 1st-order Event plane resolution (ZDC)
	tree->Branch("ResEP_2_ZDC",&ResEP_2_ZDC,"ResEP_2_ZDC/F"); // 2st-order Event plane resolution (ZDC)
	tree->Branch("ResEPSub_1_ZDC",&ResEPSub_1_ZDC,"ResEPSub_1_ZDC/F"); // 1st-order sub-Event plane resolution (ZDC)
	tree->Branch("ResEPSub_2_ZDC",&ResEPSub_2_ZDC,"ResEPSub_2_ZDC/F"); // 2st-order sub-Event plane resolution (ZDC)
	TBranch *br = tree->Branch("l0","std::vector<L0>", &vvvL); 
	tree->Branch("ptetayl0","std::vector<tuple<float,float,float> >", &vvvLpt);
	tree->Branch("nLamb",&nLamb,"nLamb/I");
	tree->Branch("nLamb_MC",&nLamb_MC,"nLamb_MC/I");
	
	Int_t events = simITS->GetEntries();
	cout << " Number of events = " << events << " Number n2 = " << n2 << endl;
	if (n2 != 0) events = TMath::Min (events, n2);
	cout << " Number of events = " << events << endl;
	
	for (Int_t i = 0; i < events; ++i) 
	{
		if (i < n1) continue;
		simITS->GetEntry(i);
		evNo = i + 1;
		cout << " Event " << evNo << " out of " << events << endl;
		TVector3 genVert;
		mcHeader->GetVertex(genVert); 
		TVector3 mom; 
		set<Int_t> idxs;
		
		Int_t nMC = mcTracks->GetEntriesFast();
		Int_t skip = 0;
		vLambPtEtaY.clear();
		
		nLamb_MC = 0;
		for (Int_t j = 0; j < nMC; ++j) 
		{
			MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(j);
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
					Double_t pt = mom.Pt();
					if (pt < 0.5) hLambPTH->Fill(2);
					if (pt > 0.5 && pt < 1.0) hLambPTH->Fill(4);
					if (pt > 1.0 && pt < 1.5) hLambPTH->Fill(6);
					if (pt > 1.5 && pt < 2.0) hLambPTH->Fill(8);
					if (pt > 2.0) hLambPTH->Fill(10);					
					vLambPtEtaY.push_back(make_tuple(pt,mom.Eta(),mcTr->GetRapidity()));
				}
				nLamb_MC++;
			}
			
			//AZ if (mcTr->GetMotherId() == -1) continue;
			if (mcTr->GetMotherId() < 0) continue;
			TVector3 pos;
			mcTr->GetStartVertex(pos);
			if (mom.Pt() < 0.001) continue;
			if (TMath::Abs(mom.Eta()) < 1.3) 
			{
				pdg = mcTr->GetPdgCode();
				moth = ((MpdMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
				if (moth == 3122) 
				{
					hAngle->Fill(mom.Pt(),mom.Angle(pos)*TMath::RadToDeg());
				}
      		}
    	}
    	if (skip) continue;
    	
		Int_t nMpdTr = 0;
		Int_t nITS = itsTracks->GetEntriesFast();
		Int_t nVert = vtxs->GetEntriesFast();
		if (event) mpdTracks = event->GetGlobalTracks();
		if (mpdTracks) nMpdTr = mpdTracks->GetEntriesFast();
		cout << " *** Event No: " << i << ", reco tracks in TPC (ITS), global: " << " " << nITS << " " << nMpdTr << ", vertices: " << nVert << endl;
		
		MpdVertex *vtx = (MpdVertex*) vtxs->First();
		mpdVert = vtx;
		vtx->Position(primVert);
		TArrayI *indxs = vtx->GetIndices();
		Int_t nPrim = indxs->GetSize();
		set<int> indxVert;
		for (Int_t k = 0; k < nPrim; ++k) indxVert.insert((*indxs)[k]);
		cout << " Number of primary (used for vertex reco) tracks: " << indxVert.size() << endl;
		
		Float_t Z_vertex = vtx->GetZ();
		Float_t Y_vertex = vtx->GetY();
		hVertex->Fill(Z_vertex,Y_vertex);
		hVertex->SetYTitle("Y_{Vertex}, [mm]");
		hVertex->SetXTitle("Z_{Vertex}, [mm]");
		
		// Find TPC track IDs 
		Int_t nPoints = 0, idMax = 0;
		for (Int_t j = 0; j < nITS; ++j) 
		{
			MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
			idMax = TMath::Max(idMax,tr->GetTrackID());
		}
		cout << " Max ID: " << idMax << endl;
		Int_t *ids = new Int_t [idMax+1];
		lays = new Int_t [idMax+1];
		Int_t *moths = new Int_t [idMax+1];
		Int_t *pdgs = new Int_t [idMax+1];
		Double_t *pt = new Double_t [idMax+1];
		Double_t *th = new Double_t [idMax+1];
		Double_t *rad = new Double_t [idMax+1];
		kProb = new Double_t [idMax+1];
		piProb = new Double_t [idMax+1];
		pProb = new Double_t [idMax+1];
		eProb = new Double_t [idMax+1];
		id2dst = new Int_t [idMax+1];
		Double_t *dZ = new Double_t [idMax+1];
		FairMCPoint **point = new FairMCPoint* [idMax+1];
		AzTrack **track = new AzTrack* [idMax+1];
		
		for (Int_t j = 0; j <= idMax; ++j) 
		{ 
			ids[j] = lays[j] = 0; 
			point[j] = 0x0;
			track[j] = 0x0;
			dZ[j] = 999999;
		} 
		
		// Get max. reached layer No.
		for (Int_t j = 0; j < nITS; ++j) 
		{
			AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
			Int_t id = tr->GetTrackID();
			ids[id]++;
			MpdKalmanHit *hit = (MpdKalmanHit*) tr->GetTrHits()->First();
			lays[id] = TMath::Max (hit->GetLayer(), lays[id]);
		}
		
		// Exclude "clones" (multiple loops)
		for (Int_t j = 0; j < nITS; ++j) 
		{
			AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
			Int_t id = tr->GetTrackID();
			if (track[id] == 0x0) track[id] = tr;		
			// MC track
			MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(id);
			mcTr->GetMomentum(mom);
			pt[id] = mom.Pt();
			th[id] = mom.Theta();
		}
		
		// Loop over DST tracks
		//calculating centrality in TPC (new calibration)
		Int_t k_check = 0;
		Int_t k_check_DCA = 0;
		for (Int_t j = 0; j < nMpdTr; ++j) 
		{
			MpdTrack *mpdTr = (MpdTrack*) mpdTracks->UncheckedAt(j);
			
			if (TMath::Abs(mpdTr->GetEta())<cut_eta && TMath::Abs(mpdTr->GetPt())>cut_pt && mpdTr->GetNofHits()>_N_Hits) 
			{
				k_check++;
				if (TMath::Sqrt(TMath::Power(mpdTr->GetDCAX(),2) + TMath::Power(mpdTr->GetDCAY(),2) + TMath::Power(mpdTr->GetDCAZ(),2)) < dca_cut) {k_check_DCA++;}
			}
			
			Int_t id = mpdTr->GetID();
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
		n_tracks_mpd = k_check;
		n_tracks_mpd_DCA = k_check_DCA;
		
		// calculating ZDC energy loss:
		for (long int j = 0; j < 90; ++j)
        {
			ZDC_energy_mpd[j] = 0;
		}
		Int_t number_of_zdchits = ZDCHits->GetEntries();
		for (Int_t zdchit_number = 0; zdchit_number < number_of_zdchits; ++zdchit_number)
		{
			ZDCHit = (MpdZdcDigi*) ZDCHits->At(zdchit_number);
		    Int_t detector_ID = ZDCHit->GetDetectorID();//1,2
	        Int_t module_ID = ZDCHit->GetModuleID();//1-45
			Double_t energy_deposit_per_hit = ZDCHit->GetELoss();
			
			ZDC_energy_mpd[(detector_ID - 1)*45 + module_ID - 1] += energy_deposit_per_hit;
		}
		
		// calculating EP angles:
		b0 = mcHeader->GetB();
		phiEP_mc = mcHeader->GetRotZ();
		phiEP_mc_or = phiEP_mc;
		phiEP_mc = TMath::ATan2(TMath::Sin(phiEP_mc),TMath::Cos(phiEP_mc));
		
		//if both ZDCs have some signal then calculate EP and its resolutions using ZDC
		if ((GetTotalEnergy(ZDC_energy_mpd, 0) != 0)&&(GetTotalEnergy(ZDC_energy_mpd, 1) != 0))
		{
			phiEP_1_ZDC = GetPsiFullZdc(ZDC_energy_mpd, 1);
			phiEP_2_ZDC = GetPsiFullZdc(ZDC_energy_mpd, 2);
			ResEP_1_ZDC = TMath::Cos(phiEP_1_ZDC - phiEP_mc);
			ResEP_2_ZDC = TMath::Cos(2*(phiEP_2_ZDC - phiEP_mc));
			
			Double_t qx,qy;
			Double_t psi_N_R_1 = GetPsiHalfZdc(ZDC_energy_mpd, 0, 1, qx, qy);
			Double_t psi_N_L_1 = GetPsiHalfZdc(ZDC_energy_mpd, 1, 1, qx, qy);
			ResEPSub_1_ZDC = TMath::Cos(psi_N_R_1 - psi_N_L_1);
			Double_t psi_N_R_2 = GetPsiHalfZdc(ZDC_energy_mpd, 0, 2, qx, qy);
			Double_t psi_N_L_2 = GetPsiHalfZdc(ZDC_energy_mpd, 1, 2, qx, qy);
			ResEPSub_2_ZDC = TMath::Cos(2*(psi_N_R_2 - psi_N_L_2));
		}	
		
		// Lambda acceptance
		multimap<Int_t,Int_t> mapLamb, mapXi;
		for (Int_t j = 0; j <= idMax; ++j) 
		{
			MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(j);
			mcTr->GetMomentum(mom);
			Int_t mothID = mcTr->GetMotherId();
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
				MpdMCTrack* moth = (MpdMCTrack*) mcTracks->UncheckedAt(mothID);
				moth->GetStartVertex(pos);
				if (pos.Mag() < 50.0) 
				{
					moths[j] = moth->GetPdgCode();
					if (moths[j] == pdgCodeL0 && (pdgs[j] == pdgCodePr || pdgs[j] == pdgCodeNeg)) mapLamb.insert(pair<Int_t,Int_t>(mothID,j));
				}
				if (moths[j] == pdgCodeXi && (pdgs[j] == pdgCodeL0 || pdgs[j] == pdgCodeNeg)) mapXi.insert(pair<Int_t,Int_t>(mothID,j));
			}
		}
		multimap<int,int>::iterator mit, mit1;
		pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;
		
		mit = mapLamb.begin();
		while (mit != mapLamb.end()) 
		{
			Int_t mothID = mit->first;
			if (mapLamb.count(mothID) != 2) {mit = mapLamb.upper_bound(mothID); continue; } // only one decay particle
			ret = mapLamb.equal_range(mothID);
			Int_t nppi[2] = {0}, nok = 0;
			Int_t nok1 = 0, nok2 = 0, nok3 = 0;

			for (mit1 = ret.first; mit1 != ret.second; ++mit1) 
			{
				MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(mit1->second);
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
			MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(mothID);
			Int_t gmID = mcTr->GetMotherId();

			if (mapXi.find(gmID) != mapXi.end()) 
			{
				ret = mapXi.equal_range(gmID);
				for (mit1 = ret.first; mit1 != ret.second; ++mit1) 
				{
					MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(mit1->second);
					if (mcTr->GetPdgCode() != pdgCodeNeg) continue;
					mcTr->GetMomentum(mom);
					if (mom.Pt() < 0.001) continue;
				}
			}
			mit = mapLamb.upper_bound(mothID);
		} // while (mit != mapLamb.end())
		
		// Track selection 
		ntr13 = 0;
		for (Int_t j = 0; j < nITS; ++j) 
		{
			AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
			if (tr->GetChi2() < -8) continue;
			Int_t id = tr->GetTrackID();
			Double_t thRec = tr->Theta();
			Double_t etaRec = tr->Momentum3().Eta();
			Double_t ptRec = tr->Momentum3().Pt();
			if (TMath::Abs(lays[id]) < -41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
			Int_t iQ = tr->Charge(); 
      
			if (tr->GetNofHits() < 10) tr->SetChi2(-9.);
			if (tr->GetChi2() < -8) continue;
			// Create MpdHelix
			MpdHelix helix = MakeHelix(tr); 
     
			// Get 3-D DCA to primary vertex
			TVector3 pca;
			Double_t s = helix.pathLength(primVert);
			pca = helix.at(s);
			pca -= primVert;
			if (iQ < 0) 
			{
				if (pdgs[id] != pdgCodeKm && pca.Mag() < gDCApi) tr->SetChi2(-9.);
			}
			else if (iQ > 0 && pca.Mag() < gDCAp) tr->SetChi2(-9.);			
			
			++ntr13;
		}   
		
// Collect "good" pions, kaons and protons		
		vector<Int_t> vecPi, vecK, vecP;
		for (Int_t j = 0; j < nITS; ++j) 
		{
			MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
			if (tr->GetChi2() < -8) continue;
			Int_t id = tr->GetTrackID();
			MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(id);
			// !!!
			if (mcTr->GetMotherId() == 0 &&
			((MpdMCTrack*)mcTracks->UncheckedAt(0))->GetPdgCode() == 1010010030) continue; // !!! decay product of artificial H3L (do i need it?)
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
		
		//~ cout << "nPi = " << vecPi.size() << "; nP = " << vecP.size() << endl;
		//~ cout << "------- Doing recoeff -------"<< endl;
		RecoEff(vecP, vecPi, 1);
		//~ cout << "nPi = " << vecPi.size() << "; nP = " << vecP.size() << endl;
		
		//~ cout << "Doing applypid"<< endl;
		ApplyPid(newPid, vecP, vecPi);  	
		//~ cout << "nPi (pid) = " << vecPi.size() << "; nP (pid) = " << vecP.size() << endl;
		
		vector<MpdParticle*> vecL;
		vecL.clear();
		vLambdas.clear();
		//~ cout << "Building lambda"<< endl;
		BuildLambda(vecP, vecPi, vecL, phiEP_mc);
		nLamb = vecL.size();
		
		//cout << "nPi = " << vecPi.size() << "; nP = " << vecP.size() << "; nL = " << vecL.size() << endl;
		//~ cout << "Filling tree"<< endl;
		tree->Fill();
		
		for (Int_t ipart = 0; ipart < nLamb; ++ipart) delete vecL[ipart];
		
		delete [] lays;
		delete [] ids;
		delete [] moths;
		delete [] pdgs;
		delete [] pt;
		delete [] th;
		delete [] point;
		delete [] rad;
		delete [] kProb;
		delete [] piProb;
		delete [] pProb;
		delete [] eProb;
		delete [] dZ;
		delete [] id2dst;
		delete [] track;
	} // for (Int_t i = 0; i < events;
	
	out.Write();
	out.Close();
	
	timer.Stop();
    Double_t rtime = timer.RealTime(), ctime = timer.CpuTime();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
	cout << "End of the program" << endl;
	
	return(0);
}

//__________________________________________________________________________
TChain* Chain(Int_t nFiles, TString firstFile)
{
	// File name parsing

	Int_t ndash = 0, ndashOK = 2;
  
	// Get first file number
	Int_t leng = firstFile.Length(), i1 = 0, i2 = 0;
	cout << leng << endl;
	TString numb, prefix, suffix, symb, symb0;
	cout << numb.Length() << endl;
	for (Int_t i = leng-1; i > -1; --i) {
		symb = firstFile(i,1);
		if (symb == "_" || symb == "-") {
			++ndash;
			if (ndash < ndashOK) continue;
			prefix = firstFile(0,i+1);
			i1 = i + 1;
			break;
		} else if (symb == ".") {
			suffix = firstFile(i,leng-i);
			i2 = i - 1;
		}
	}
	numb = firstFile(i1,i2-i1+1);

	Int_t numb0 = numb.Atoi();
	cout << numb << endl;
	cout << numb0 << endl;
	cout << prefix << endl;
	cout << suffix << endl;

	TChain *chain = new TChain("mpdsim");
	TString fileName;
	nFiles += numb0;
	for (Int_t i = numb0; i < nFiles; ++i) {
		fileName = prefix;
		fileName += TString::Format("%04d",i);
		fileName += suffix;
		if (!fileName.Contains("root:")) {
			// local file
			if (!gSystem->FindFile("./",fileName)) break;
		} else {
			// xrootd
			TFile *f0 = TFile::Open(fileName);
			if (f0 == NULL) break;
			f0->Close();
		}
		chain->AddFile(fileName);
	}
	chain->ls();
	return chain;
}
//__________________________________________________________________________
TChain* ChainFile(Int_t nFiles, TString fileNameList, Int_t skipLines)
{
	// Read "nFiles" lines from the file "fileNameList", containing input file names,
	// skipping "skipLines" lines

	ifstream fin(fileNameList.Data());
	string chline;
	TString fileName;

	for (Int_t line = 0 ; line < skipLines; ++line) getline(fin,chline);

	TChain *chain = new TChain("mpdsim");
	for (Int_t line = 0; line < nFiles; ++line) {
		getline(fin,chline);
		Int_t i = chline.rfind(" ");
		fileName = chline.substr(i+1,string::npos);
		chain->AddFile(fileName);
		cout << "Entries = " << chain->GetEntries() << endl; //to check how many events in each file
	}
	chain->ls();
	return chain;
}
//--------------------
void RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi, Int_t pid)
{
	// Check reco efficiency

	Int_t nPi = vecPi.size(), nP = vecP.size();

	for (Int_t ip = nP - 1; ip >= 0; --ip) 
	{
		// AntiProton
		AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
		MpdMCTrack *mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(trP->GetTrackID());
		Int_t mothId = mcTr->GetMotherId();
		if (mothId < 0) continue;
		MpdMCTrack *moth = (MpdMCTrack*) mcTracks->UncheckedAt(mothId);  
		if (moth->GetPdgCode() == pdgCodeL0) 
		{
			Int_t mp = mothId;
			// Proton from Lambda
			for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) 
			{
				// Pion
				AzTrack *trPi = (AzTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
				MpdMCTrack *mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
				Int_t mothId = mcTr->GetMotherId();
				if (mothId < 0) continue;
				MpdMCTrack *moth = (MpdMCTrack*) mcTracks->UncheckedAt(mothId);  
				if (moth->GetPdgCode() == pdgCodeL0 && mp == mothId) 
				{
					//AZ - flag decay tracks to check PID influence later
					trP->SetUniqueID(mothId+1);
					trPi->SetUniqueID(mothId+1);
					//
					Int_t gmId = moth->GetMotherId();
					if (gmId >= 0) {
						MpdMCTrack *gmoth = (MpdMCTrack*) mcTracks->UncheckedAt(gmId);
						if (gmoth->GetPdgCode() == pdgCodeXi) 
						{
							for (Int_t kpi = nPi - 1; kpi >= 0; --kpi) 
							{
								// Pion
								AzTrack *trK = (AzTrack*) itsTracks->UncheckedAt(vecPi[kpi]);
								MpdMCTrack *mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(trK->GetTrackID());
								Int_t mothId = mcTr->GetMotherId();
								if (mothId < 0) continue;
								MpdMCTrack *moth = (MpdMCTrack*) mcTracks->UncheckedAt(mothId);  
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

	if (pid) return; // skip the rest if PID is used

	for (Int_t ip = nP - 1; ip >= 0; --ip) 
	{
		// Proton
		AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
		AzTrack trCor = *trP;
		trCor.SetDirection(MpdKalmanTrack::kInward);
		if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
		else recoTpc->Refit(&trCor, 0.93827, 1); // refit
		MpdParticle prot(trCor, vecP[ip]);
		prot.SetPdg(pdgCodePr);
		prot.SetMass();

		Double_t chi2 = TMath::Min (prot.Chi2Vertex(mpdVert),999.);
		if (chi2 < gC2p) vecP.erase(vecP.begin()+ip);
	}

	if (nP) 
	{
		for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) 
		{
			// Pion
			AzTrack *trPi = (AzTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
			AzTrack trCor = *trPi;
 
			Double_t chi2 = TMath::Min (trPi->GetChi2Vertex(),999.);
			if (chi2 < gC2pi) vecPi.erase(vecPi.begin()+jpi);
		}
	}
}
//--------------------___________________________________________
MpdHelix MakeHelix(const MpdKalmanTrack *tr) 
{
  Double_t r = tr->GetPosNew();
  Double_t phi = tr->GetParam(0) / r;
  Double_t x = r * TMath::Cos(phi);
  Double_t y = r * TMath::Sin(phi);
  Double_t dip = tr->GetParam(3);
  Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
  cur *= TMath::Abs (tr->GetParam(4));
  TVector3 o(x, y, tr->GetParam(1));
  Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
  MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
  return helix;
}
//--------------------___________________________________________
MpdHelix MakeHelix(const MpdParticle *part) 
{
	Double_t dip = TMath::PiOver2() - part->Theta();
	Double_t cur = TMath::Abs (part->GetMeas(4));
	if (part->GetCharge() == 0) cur = numeric_limits<double>::epsilon();
	Int_t h = (Int_t) TMath::Sign(1.1,part->GetMeas(4));
	Double_t phase = part->GetMeas(2) - TMath::PiOver2() * h;
	Double_t x = part->GetXY(0);
	Double_t y = part->GetXY(1);
	TVector3 o(x, y, part->GetMeas(1));
	MpdHelix helix(cur, dip, phase, o, h);
	return helix;
}  
//--------------------________________________________________
void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL, Float_t &phiEP) 
{
	// Make antilambdas

	Int_t nPi = vecPi.size(), nP = vecP.size();
	vector<MpdParticle*> vPart;
	vecL1.clear();
	//vecL2.clear();

	for (Int_t ip = 0; ip < nP; ++ip) {
		// Proton
		AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
		MpdMCTrack *mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(trP->GetTrackID());
		TVector3 mcmom1;
		mcTr->GetMomentum(mcmom1);
		mcps[1] = mcmom1.Mag();
		mcphis[1] = mcmom1.Phi();
		mcthetas[1] = mcmom1.Theta();
		Int_t mothId = mcTr->GetMotherId();
		AzTrack trCor = *trP;
		trCor.SetDirection(MpdKalmanTrack::kInward);
		if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
		else recoTpc->Refit(&trCor, 0.93827, 1); // refit
		MpdParticle prot(trCor, vecP[ip]);
		prot.SetPdg(pdgCodePr);
		prot.SetMass();
		qs[1] = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodePr)->Charge()/3); //again the same error
		etas[1] = trP->Momentum3().Eta();
		
    
		chi2s[1] = TMath::Min (prot.Chi2Vertex(mpdVert),9999.);
		layMx[1] = TMath::Abs (lays[trP->GetTrackID()]);
		MpdHelix helix = MakeHelix(trP);
		//Get 3-D DCA to primary vertex
		TVector3 pca;
		Double_t s = helix.pathLength(primVert);
		pca = helix.at(s);
		pca -= primVert;
		dcas[1] = pca.Mag();
		origs[1] = 0;
		if (mothId >= 0 && ((MpdMCTrack*) mcTracks->UncheckedAt(mothId))->GetPdgCode() == pdgCodeL0)
		origs[1] = -1; // from lambda

		for (Int_t jpi = 0; jpi < nPi; ++jpi) {
			// Pion
			MpdKalmanTrack *trPi = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
			MpdMCTrack *mcTr1 = (MpdMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
			TVector3 mcmom2;
			mcTr1->GetMomentum(mcmom2);
			mcps[0] = mcmom2.Mag();
			mcphis[0] = mcmom2.Phi();
			mcthetas[0] = mcmom2.Theta();
			Int_t mothId1 = mcTr1->GetMotherId();
			origs[0] = 0;
			if (mothId1 >= 0 && ((MpdMCTrack*) mcTracks->UncheckedAt(mothId1))->GetPdgCode() == pdgCodeL0)
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
			v0 -= primVert;
			Double_t decay = v0.Mag();
			path = TMath::Sign (decay, v0*lambPart.Momentum3());

			if (chi2 >= 0 && chi2 < gC2L && path > gPathL) {
				if (origs[1] > 0) origs[1] = -1;
				MpdMCTrack *moth = NULL;
				((TH1D*)gROOT->FindObjectAny("hMassL"))->Fill(lambPart.GetMass());
				if (mothId != mothId1 || mothId < 0) {
					((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
				} else {
					if (origs[0] == -1) {
						((TH1D*)gROOT->FindObjectAny("hMassLsig"))->Fill(lambPart.GetMass());
						origs[0] = origs[1] = 1;
						moth = (MpdMCTrack*) mcTracks->UncheckedAt(mothId);
					}
					else ((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
				}

				// Fill tree
	
				qs[0] = TMath::Nint(TDatabasePDG::Instance()->GetParticle(pdgCodeNeg)->Charge()/3); // pion
				etas[0] = trPi->Momentum3().Eta();
	
	
				chi2s[0] = TMath::Min (pion->Chi2Vertex(mpdVert),9999.);
				layMx[0] = TMath::Abs (lays[trPi->GetTrackID()]);
				MpdHelix helix1 = MakeHelix(trPi);
				// Get 3-D DCA to primary vertex
				s = helix1.pathLength(primVert);
				pca = helix1.at(s);
				pca -= primVert;
				dcas[0] = pca.Mag();

				massh = lambPart.GetMass();
				chi2h = chi2;
				angle = v0.Angle(lambPart.Momentum3());
				pth = lambPart.Pt(); // reconstructed
				ph = lambPart.Momentum(); // reconstructed
				phi_Lam = lambPart.Phi(); // reconstructed (new addition!)
				if (pth > 0.001) etah = lambPart.Momentum3().Eta(); 
				else etah = TMath::Sign(100.,lambPart.Momentum3().Z()); 
				pair<Double_t,Double_t> paths = helix.pathLengths(helix1);
				TVector3 p1 = helix.at(paths.first);
				TVector3 p2 = helix1.at(paths.second);
				p1 -= p2;
				disth = p1.Mag(); // closest distance between daughters
				// Get 3-D DCA of lambda to primary vertex
				MpdHelix helix2 = MakeHelix(&lambPart);
				s = helix2.pathLength(primVert);
				pca = helix2.at(s);
				pca -= primVert;
				dca = pca.Mag();
				c2pv = TMath::Min (lambPart.Chi2Vertex(mpdVert),9999.);
				omega1 = dcas[0] * dcas[1] / (dca * dca + disth * disth);
				omega2 = TMath::Sqrt (chi2s[0] * chi2s[1]) / (c2pv + chi2h);

				dstNo[0] = vecPi[jpi]; // pion index
				dstNo[1] = vecP[ip]; // proton index

				
				if (lambPart.GetMass() >= 1.10518 && lambPart.GetMass() <= 1.12668) { // lambda mass +- 5*2.15 MeV
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
	
				if (origs[0] == 1) {
					// True lambda
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
					Float_t weight_pol = ((MpdMCTrack*)mcTracks->UncheckedAt(mothId))->GetWeight();
					polarhx = ((MpdMCTrack*)mcTracks->UncheckedAt(mothId))->GetPolar(0);
					polarhy = ((MpdMCTrack*)mcTracks->UncheckedAt(mothId))->GetPolar(1);
					polarhz = ((MpdMCTrack*)mcTracks->UncheckedAt(mothId))->GetPolar(2);
					TVector3 polar_changed(polarhx, polarhy, polarhz);
					if (phiEP != 0.) polar_changed.RotateZ(-phiEP);
					polarhx = weight_pol*polar_changed.X();
					polarhy = weight_pol*polar_changed.Y();
					polarhz = weight_pol*polar_changed.Z();
					
					((TH1D*)gROOT->FindObjectAny("hPolarL"))->GetXaxis()->SetTitle("Global Polarization (Full)");
					((TH1D*)gROOT->FindObjectAny("hPolarL"))->GetYaxis()->SetTitle("Entries ");
					((TH1D*)gROOT->FindObjectAny("hPolarL"))->Fill(polarhy);
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
//--------------------_____________________________________
void FindPolarAngle(MpdParticle &lamb, vector<MpdParticle*> &vPart)
{
	// Compute decay proton angle w.r.t. lambda decay plane

	TVector3 vPr, vPi, vLamb;
	TLorentzVector prLor, piLor, lambLor; 

	// True (exact) parameters
	vPr.SetMagThetaPhi(mcps[1], mcthetas[1], mcphis[1]);
	vPi.SetMagThetaPhi(mcps[0], mcthetas[0], mcphis[0]);
	//cout << "mcps[1] = " << mcps[1] << "; mcthetas[1] = " << mcthetas[1] << "; mcphis[1] = " << mcphis[1] << endl;
	//cout << "mcps[0] = " << mcps[0] << "; mcthetas[0] = " << mcthetas[0] << "; mcphis[0] = " << mcphis[0] << endl;
      
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

	//cout << "cosAmc = " << cosAmc << "; phi_star_MC = " << phi_star_MC << endl;

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
	
	//cout << "cosA = " << cosA << "; phi_star = " << phi_star << endl;
}
//--------------------_____________________________________
Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu)
{
	// Compute distance between helix and straight line

	// Create MpdHelix
	trC = MakeHelix(helix);
	vtxN.SetXYZ(neu->Getx()(0,0), neu->Getx()(1,0), neu->Getx()(2,0));
	momN = neu->Momentum3();
	momN *= (1. / momN.Mag());

	gMinuit->SetFCN(fcn);

	Double_t arglist[10];
	Int_t ierflg = 0;
	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
	arglist[0] = -1; //1; //-1;
	gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);

	Double_t vstart[2] = {-0.1,-0.1};
	static Double_t step[2] = {0.1, 0.1};	
	gMinuit->mnparm(0, "lengN", vstart[0], step[0], 0,0,ierflg);
	gMinuit->mnparm(1, "lengC", vstart[1], step[1], 0,0,ierflg);
    
	// Now ready for minimization step
	arglist[0] = 500;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
	// Get results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	return amin;
}
//--------------------___________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	// Compute distance between straight line and helix

	TVector3 mom = momN;
	mom *= par[0];
	TVector3 posN = vtxN;
	posN += mom;
  
	TVector3 posC = trC.at(par[1]);
	posC -= posN;
	f = posC.Mag();
}
//--------------------
void ApplyPid(MpdPid *pid, vector<Int_t> &vecP, vector<Int_t> &vecPi)
{
	// Apply PID
	//~ cout << "Apply PID started"<< endl;
	//AZ - get information on hyperon decay products
	map<Int_t,set<Int_t> > mapL, mapXi;

	Int_t nP = vecP.size(), nPi = vecPi.size();
	
	//cout << "PID started ................" << endl; 
	//cout << "nP = "<< nP << "; nPi = " << nPi << endl;
	
	for (Int_t ip = 0; ip < nP; ++ip) 
	{
		//cout << "First cycle started"<< endl;
		AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
		if (trP->GetUniqueID() > 0) 
		{
			// Lambda decay product
      			Int_t mid = trP->GetUniqueID();
      			if (mapL.find(mid) == mapL.end()) { set<Int_t> aaa; mapL[mid] = aaa; }
      			mapL[mid].insert(vecP[ip]);
      			trP->SetUniqueID(0); // reset 
   		}
    		if (trP->GetVertex().GetUniqueID() > 0) 
    		{
      		// Xi- decay product
      		Int_t mid = trP->GetVertex().GetUniqueID();
      		if (mapXi.find(mid) == mapXi.end()) { set<Int_t> aaa; mapXi[mid] = aaa; }
      		mapXi[mid].insert(vecP[ip]);
    		}
	}

	for (Int_t ip = 0; ip < nPi; ++ip) 
	{
		//cout << "Second cycle started"<< endl;
		AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecPi[ip]);
		if (trP->GetUniqueID() > 0) 
		{
			// Lambda decay product
			Int_t mid = trP->GetUniqueID();
			if (mapL.find(mid) == mapL.end()) { set<Int_t> aaa; mapL[mid] = aaa; }
			mapL[mid].insert(vecPi[ip]);
			trP->SetUniqueID(0); // reset 
		}
		if (trP->GetVertex().GetUniqueID() > 0) 
		{
			// Xi- decay product
			Int_t mid = trP->GetVertex().GetUniqueID();
			if (mapXi.find(mid) == mapXi.end()) { set<Int_t> aaa; mapXi[mid] = aaa; }
			mapXi[mid].insert(vecPi[ip]);
		}
	}


	vecP.clear();
	vecPi.clear();

	Int_t nITS = itsTracks->GetEntriesFast();
	//~ cout << "nITS = "<< nITS << endl;
	for (Int_t j = 0; j < nITS; ++j) 
	{
		//~ cout << "Third cycle started"<< endl;
		AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
		if (tr->GetChi2() < -8) continue;
		Int_t id = tr->GetTrackID();
		MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(id);
		Int_t mothId = mcTr->GetMotherId();
		Int_t uid = tr->GetVertex().GetUniqueID();
		//cout << "blub1"<< endl;
   
		MpdTrack* mpdTrack = (MpdTrack*) mpdTracks->UncheckedAt(j);
		if (mpdTrack->GetID() != id) { cout << id << " " << mpdTrack->GetID() << endl; Fatal("ApplyPid"," Different ID"); }
		//cout << "blub2"<< endl;
		Int_t ret = 0, charge = tr->Charge(), tofFlag = mpdTrack->GetTofFlag();
		//Double_t dedx = tr->GetPartID(), m2 = mpdTrack->GetTofMass2();
		Double_t dedx = tr->GetDedx(), m2 = mpdTrack->GetTofMass2();
		
		if (tofFlag == 2 || tofFlag == 6)          // dE/dx+TOF
		ret = pid->FillProbs(tr->Momentum(), dedx, m2, charge);
		if (ret == 0) ret = pid->FillProbs(tr->Momentum(), dedx, charge);
		//cout << "blub3"<< endl;
		TH1D *hFlag = (TH1D*) gROOT->FindObjectAny("hPIDflag");
		if (ret == 0) {
			// No PID
			if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(2.1); // lost pion
			if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(6.1); // lost proton
			continue;
		}
		//~ cout << "blub4"<< endl;
		Double_t piThr = -0.75, probThr = -0.60;
		//~ cout << "step = " << j << " out of " << nITS << endl;
		if (pdgCodeL0 * tr->Charge() < 0) 
		{
			Double_t prob = pid->GetProbPi();
			//cout << "probPi = " << prob << endl;
			if (prob > piThr && prob > pid->GetProbKa() && prob > pid->GetProbEl() && prob > pid->GetProbPr() && prob > pid->GetProbMu()) 
			{
				// "pion"
				//~ cout << "Found pion! " << endl;
				if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(0.1); // correct pion
				else if (mcTr->GetPdgCode() != pdgCodeNeg) hFlag->Fill(1.1); // false pion
				//
				if (mapL.find(mothId+1) != mapL.end() && mapL[mothId+1].find(j) != mapL[mothId+1].end())
				mapL[mothId+1].erase(j);
				if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end())
				mapXi[uid].erase(j);
				//
				Double_t chi2 = TMath::Min (tr->GetChi2Vertex(),999.);
				//~ cout << "chi2 = " << chi2 << "; gC2pi = " << gC2pi << endl;
				if (chi2 < gC2pi) continue;
				vecPi.push_back(j);
			} else if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(2.1); // lost pion
		} else 
		{
			Double_t prob = pid->GetProbPr();
			//cout << "probPr = " << prob << endl;
			if (prob > probThr && prob > pid->GetProbKa() && prob > pid->GetProbPi() && prob > pid->GetProbDe()) 
			{
				// "proton"
				//~ cout << "Found proton! " << endl;
				if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(4.1); // correct proton
				else if (mcTr->GetPdgCode() != pdgCodePr) hFlag->Fill(5.1); // false proton
//
				if (mapL.find(mothId+1) != mapL.end() && mapL[mothId+1].find(j) != mapL[mothId+1].end())
				mapL[mothId+1].erase(j);
				if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end())
				mapXi[uid].erase(j);
				//
				AzTrack trCor = *tr;
				trCor.SetDirection(MpdKalmanTrack::kInward);
				//~ cout << "Starting refit " << endl;
				if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
				else recoTpc->Refit(&trCor, 0.93827, 1); // refit
				//~ cout << "Finished refit " << endl;
				//~ ShowTrackStats(tr);
				MpdParticle prot(trCor, 0);
				//~ cout << "Setting pdg" << endl;
				prot.SetPdg(pdgCodePr);
				//~ cout << "Setting mass" << endl;
				prot.SetMass();
				//~ cout << "Calculating chi2" << endl;
				Double_t chi2 = TMath::Min (prot.Chi2Vertex(mpdVert),999.);
				//~ cout << "chi2 = " << chi2 << "; gC2p = " << gC2p << endl;
				if (chi2 < gC2p) continue;
				vecP.push_back(j);
			} else if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(6.1); // lost proton
		}
		//cout << "blub5"<< endl;
		// GetProbKa(), GetProbEl(), GetProbPr(), GetProbDe(), GetProbHe3
	}    
	//~ cout << "Almost done"<< endl;
	//
	Int_t nLok = 0, nXiok = 0;
	for (map<Int_t,set<Int_t> >::iterator mit = mapL.begin(); mit != mapL.end(); ++mit) 
	{
		if (mit->second.size() == 0) ++nLok;
	}
	for (map<Int_t,set<Int_t> >::iterator mit = mapXi.begin(); mit != mapXi.end(); ++mit) 
	{
		if (mit->second.size() == 0) ++nXiok;
	}
	//~ cout << "Apply PID finished"<< endl;
}

//
Double_t GetPsiHalfZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t n, Double_t &qx, Double_t &qy)
{
	Double_t Qcos=0., Qsin=0.;
	GetQsZdc(zdc_energy,zdc_ID,n,Qcos,Qsin);

	qx = Qcos;
	qy = Qsin;
	Double_t PsiEP = (1/(Double_t)n) * TMath::ATan2(Qsin,Qcos);
	return PsiEP;
}

//
Double_t GetPsiFullZdc(Float_t* zdc_energy, Int_t n)
{
	Double_t QcosR=0., QsinR=0.;
	GetQsZdc(zdc_energy,0,n,QcosR, QsinR);

	Double_t QcosL=0., QsinL=0.;
	GetQsZdc(zdc_energy,1,n,QcosL,QsinL);

	Double_t psiEP = TMath::ATan2(QsinR + QsinL,QcosR + QcosL)/(Double_t) n; // (-pi/n,pi/n]
    return psiEP;
}
//
void GetQsZdc(Float_t* zdc_energy, Int_t zdc_ID, Int_t harm, Double_t &Qx, Double_t &Qy)
{
	Double_t *phi_angle_of_modules = GetAngles();
	Double_t Qcos = 0 , Qsin = 0;
	Double_t w_sign = 0;
	Double_t total_energy = GetTotalEnergy(zdc_energy,zdc_ID);

	if (harm == 2) w_sign = 1;
	else if (harm == 1)
	{
		if (zdc_ID == 0) w_sign = 1;
		else if (zdc_ID == 1) w_sign = -1;
	}


	for (int module = _N_MODULES_TOTAL/2 * zdc_ID; module < _N_MODULES_TOTAL/2 * (zdc_ID + 1); ++module)
	{
		if ((module==22) || (module==67)) continue;
//		if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
//				(i==66)||(i==68)||(i==74)) continue;
		Qcos += w_sign*zdc_energy[module] / total_energy * TMath::Cos(harm*phi_angle_of_modules[module]);
        Qsin += w_sign*zdc_energy[module] / total_energy * TMath::Sin(harm*phi_angle_of_modules[module]);
	}

	Qx = Qcos; Qy = Qsin;
}
//
Double_t GetTotalEnergy(Float_t* zdc_energy, Int_t zdc_ID)
{
	Double_t total_energy = 0.;
	for (int i = _N_MODULES_TOTAL/2 * zdc_ID; i < _N_MODULES_TOTAL/2 * (zdc_ID + 1); ++i)
	{
		if ((i==22) || (i==67)) continue;
//		if ((i ==15)||(i ==21)||(i==23)||(i ==29)||(i==60)||
//				(i==66)||(i==68)||(i==74)) continue;
		total_energy += zdc_energy[i];
	}
	return total_energy;
}
Double_t* GetAngles()
{
	Double_t *phi_angle_of_module = new Double_t[_N_MODULES_TOTAL];

	for (int i = 0; i < _N_ARM; ++i)
	{
		Int_t x_axis_switch;
		if (i == 0) x_axis_switch = 1;
		else if (i == 1) x_axis_switch = -1;

		for (Int_t j = 0; j < _N_MODULES_TOTAL/2; ++j)
		{
			Double_t y = 0, x = 0;

			if ((j>=0) && (j<=4))
			{
				y = 45., x = (j-2)*15.;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = TMath::ATan2(y,x_axis_switch*x);
			}
			else if ((j>=5) && (j<=39))
			{
				y = (3-(j+2)/7)*15, x = (3-(j+2)%7)*15;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = TMath::ATan2(y,x_axis_switch*x);
			}
			else if ((j>=40) && (j<=44))
			{
				y = -45. , x = (j-42)*15.;
				phi_angle_of_module[j + i*_N_MODULES_TOTAL/2] = TMath::ATan2(y,x_axis_switch*x);
			}
		}
	}

	return phi_angle_of_module;
}
