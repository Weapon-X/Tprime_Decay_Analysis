
#ifndef lvb_qqqq_new
#define lvb_qqqq_new

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//Header file for the classes stored in the TTree if any.
#include <TVector.h>
#include "TString.h"
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>

#include <TH2D.h>
#include <TH1I.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
//#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>
#include <map>
//#include "TRFIOFile.h"
#include "TMath.h"
#include <vector>
#include <TList.h>
#include <TLatex.h>
#include <Riostream.h>
#include <set>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TKDE.h"

#include <iostream>
#include <array>

using namespace std;
using namespace TMath; // To include TMath class functions
using namespace ROOT;


class lvb_qqqq_dCSV {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
		Int_t           v_event;
		Int_t           n_Vtx;
		Int_t           n_GoodVtx;
		Int_t           n_TrksPV;
		Bool_t          is_PVGood;
		Float_t         v_vtx;
		Float_t         v_vty;
		Float_t         v_vtz;
		Float_t         rho_Central;
		Int_t           n_MC;
		vector<int>     *mc_PID;
		vector<float>   *mc_Vtx;
		vector<float>   *mc_Vty;
		vector<float>   *mc_Vtz;
		vector<float>   *mc_Pt;
		vector<float>   *mc_Mass;
		vector<float>   *mc_Eta;
		vector<float>   *mc_Phi;
		vector<float>   *mc_E;
		vector<float>   *mc_Et;
		vector<int>     *mc_GMomPID;
		vector<int>     *mc_MomPID;
		vector<float>   *mc_MomPt;
		vector<float>   *mc_MomMass;
		vector<float>   *mc_MomEta;
		vector<float>   *mc_MomPhi;
		vector<unsigned short> *mc_StatusFlag;
		vector<int>     *mc_Parentage;
		vector<int>     *mc_Status;
		vector<float>   *mc_CalIsoDR03;
		vector<float>   *mc_TrkIsoDR03;
		vector<float>   *mc_CalIsoDR04;
		vector<float>   *mc_TrkIsoDR04;
		Float_t         gen_MET;
		Float_t         gen_METPhi;
		Float_t         pf_MET;
		Float_t         pf_METPhi;
		Float_t         pf_METsumEt;
		Float_t         pf_METmEtSig;
		Float_t         pf_METSig;
		Int_t           n_Ele;
		vector<int>     *ele_Charge;
		vector<float>   *ele_En;
		vector<float>   *ele_D0;
		vector<float>   *ele_Dz;
		vector<float>   *ele_Pt;
		vector<float>   *ele_Eta;
		vector<float>   *ele_Phi;
		vector<float>   *ele_R9;
		vector<float>   *ele_SCEta;
		vector<float>   *ele_SCPhi;
		vector<float>   *ele_HoverE;
		vector<float>   *ele_EoverP;
		vector<float>   *ele_EoverPout;
		vector<float>   *ele_EoverPInv;
		vector<float>   *ele_dEtaAtVtx;
		vector<float>   *ele_dPhiAtVtx;
		vector<float>   *ele_SigmaIEtaIEtaFull5x5;
		vector<float>   *ele_SigmaIPhiIPhiFull5x5;
		vector<int>     *ele_ConvVeto;
		vector<int>     *ele_MissHits;
		vector<float>   *ele_PFChIso;
		vector<float>   *ele_PFPhoIso;
		vector<float>   *ele_PFNeuIso;
		vector<float>   *ele_PFMiniIso;
		vector<float>   *ele_dEtaseedAtVtx;
		Int_t           N_Mu;
		vector<float>   *mu_Pt;
		vector<float>   *mu_En;
		vector<float>   *mu_Eta;
		vector<float>   *mu_Phi;
		vector<int>     *mu_Charge;
		vector<unsigned short> *mu_IDbit;
		vector<float>   *mu_D0;
		vector<float>   *mu_Dz;
		vector<float>   *mu_Chi2NDF;
		vector<float>   *mu_InnerD0;
		vector<float>   *mu_InnerDz;
		vector<float>   *mu_InnervalidFraction;
		vector<float>   *mu_segmentCompatibility;
		vector<float>   *mu_chi2LocalPosition;
		vector<float>   *mu_trkKink;
		vector<float>   *mu_PFChIso;
		vector<float>   *mu_PFPhoIso;
		vector<float>   *mu_PFNeuIso;
		vector<float>   *mu_PFMiniIso;
		vector<int>     *mu_TrkLayers;
		vector<int>     *mu_PixelLayers;
		vector<int>     *mu_PixelHits;
		vector<int>     *mu_MuonHits;
		vector<int>     *mu_Stations;
		vector<int>     *mu_Matches;
		vector<int>     *mu_TrkQuality;
		vector<float>   *mu_IsoTrk;

		Int_t           n_Jet;
		vector<float>   *jet_Pt;
		vector<float>   *jet_En;
		vector<float>   *jet_Eta;
		vector<float>   *jet_Phi;
		vector<float>   *jet_Area;
		vector<float>   *jet_Mt;
		vector<float>   *jet_CSV2BJetTags;
		vector<float>   *jet_JetProbabilityBJetTags;
		vector<float>   *jet_pfCombinedMVAV2BJetTags;
		vector<float>   *jet_DeepCSVTags_b;
		vector<float>   *jet_DeepCSVTags_bb;
		vector<float>   *jet_DeepCSVTags_c;
		vector<float>   *jet_DeepCSVTags_cc;
		vector<float>   *jet_DeepCSVTags_udsg;
		vector<int>     *jet_PartonID;
		vector<int>     *jet_HadFlvr;
		vector<float>   *jet_GenJetEn;
		vector<float>   *jet_GenJetPt;
		vector<float>   *jet_GenJetEta;
		vector<float>   *jet_GenJetPhi;
		vector<int>     *jet_GenPartonID;
		vector<float>   *jet_GenEn;
		vector<float>   *jet_GenPt;
		vector<float>   *jet_GenEta;
		vector<float>   *jet_GenPhi;
		vector<int>     *jet_GenPartonMomID;
		vector<bool>    *jet_PFLooseId;
		vector<int>     *jet_ID;
		vector<float>   *jet_PUID;
		vector<int>     *jet_PUFullID;
		vector<float>   *jet_CHF;
		vector<float>   *jet_NHF;
		vector<float>   *jet_CEF;
		vector<float>   *jet_NEF;
		vector<int>     *jet_NCH;
		vector<int>     *jet_NNP;
		vector<float>   *jet_MUF;
		Int_t           N_AK8Jet;
		vector<float>   *AK8_JetPt;
		vector<float>   *AK8_JetEn;
		vector<float>   *AK8_JetEta;
		vector<float>   *AK8_JetPhi;
		vector<float>   *AK8_JetMass;
		vector<float>   *AK8_Jet_tau1;
		vector<float>   *AK8_Jet_tau2;
		vector<float>   *AK8_Jet_tau3;
		vector<float>   *AK8_Jet_tau4;
		vector<float>   *AK8_JetCHF;
		vector<float>   *AK8_JetNHF;
		vector<float>   *AK8_JetCEF;
		vector<float>   *AK8_JetNEF;
		vector<int>     *AK8_JetNCH;
		vector<int>     *AK8_JetNNP;
		vector<float>   *AK8_JetMUF;
		vector<int>     *AK8_Jetnconstituents;
		vector<bool>    *AK8_JetPFLooseId;
		vector<bool>    *AK8_JetPFTightLepVetoId;
		vector<float>   *AK8_JetSoftDropMass;
		vector<float>   *AK8_JetPrunedMass;
		vector<float>   *AK8_JetpfBoostedDSVBTag;
		vector<float>   *AK8_JetCSV;
		vector<float>   *AK8_JetDSVnewV4;
		vector<float>   *AK8_puppiPt;
		vector<float>   *AK8_puppiMass;
		vector<float>   *AK8_puppiEta;
		vector<float>   *AK8_puppiPhi;
		vector<float>   *AK8_puppiTau1;
		vector<float>   *AK8_puppiTau2;
		vector<float>   *AK8_puppiTau3;
		vector<float>   *AK8_puppiSDMass;
		vector<int>     *AK8_JetPartonID;
		vector<int>     *AK8_JetHadFlvr;
		vector<int>     *AK8_JetGenJetIndex;
		vector<float>   *AK8_JetGenJetEn;
		vector<float>   *AK8_JetGenJetPt;
		vector<float>   *AK8_JetGenJetEta;
		vector<float>   *AK8_JetGenJetPhi;
		vector<int>     *AK8_JetGenPartonID;
		vector<float>   *AK8_JetGenEn;
		vector<float>   *AK8_JetGenPt;
		vector<float>   *AK8_JetGenEta;
		vector<float>   *AK8_JetGenPhi;
		vector<int>     *AK8_JetGenPartonMomID;
		vector<int>     *n_AK8SDSJ;
		vector<vector<float> > *AK8_SDSJPt;
		vector<vector<float> > *AK8_SDSJEta;
		vector<vector<float> > *AK8_SDSJPhi;
		vector<vector<float> > *AK8_SDSJMass;
		vector<vector<float> > *AK8_SDSJE;
		vector<vector<int> > *AK8_SDSJCharge;
		vector<vector<int> > *AK8_SDSJFlavour;
		vector<vector<float> > *AK8_SDSJCSV;
		vector<int>     *n_AK8puppiSDSJ;
		vector<vector<float> > *AK8_puppiSDSJPt;
		vector<vector<float> > *AK8_puppiSDSJEta;
		vector<vector<float> > *AK8_puppiSDSJPhi;
		vector<vector<float> > *AK8_puppiSDSJMass;
		vector<vector<float> > *AK8_puppiSDSJE;
		vector<vector<int> > *AK8_puppiSDSJCharge;
		vector<vector<int> > *AK8_puppiSDSJFlavour;
		vector<vector<float> > *AK8_puppiSDSJCSV;

		// List of branches
		TBranch        *b_v_event;   //!
		TBranch        *b_n_Vtx;   //!
		TBranch        *b_n_GoodVtx;   //!
		TBranch        *b_n_TrksPV;   //!
		TBranch        *b_is_PVGood;   //!
		TBranch        *b_v_vtx;   //!
		TBranch        *b_v_vty;   //!
		TBranch        *b_v_vtz;   //!
		TBranch        *b_rho_Central;   //!
		TBranch        *b_n_MC;   //!
		TBranch        *b_mc_PID;   //!
		TBranch        *b_mc_Vtx;   //!
		TBranch        *b_mc_Vty;   //!
		TBranch        *b_mc_Vtz;   //!
		TBranch        *b_mc_Pt;   //!
		TBranch        *b_mc_Mass;   //!
		TBranch        *b_mc_Eta;   //!
		TBranch        *b_mc_Phi;   //!
		TBranch        *b_mc_E;   //!
		TBranch        *b_mc_Et;   //!
		TBranch        *b_mc_GMomPID;   //!
		TBranch        *b_mc_MomPID;   //!
		TBranch        *b_mc_MomPt;   //!
		TBranch        *b_mc_MomMass;   //!
		TBranch        *b_mc_MomEta;   //!
		TBranch        *b_mc_MomPhi;   //!
		TBranch        *b_mc_StatusFlag;   //!
		TBranch        *b_mc_Parentage;   //!
		TBranch        *b_mc_Status;   //!
		TBranch        *b_mc_CalIsoDR03;   //!
		TBranch        *b_mc_TrkIsoDR03;   //!
		TBranch        *b_mc_CalIsoDR04;   //!
		TBranch        *b_mc_TrkIsoDR04;   //!
		TBranch        *b_gen_MET;   //!
		TBranch        *b_gen_METPhi;   //!
		TBranch        *b_pf_MET;   //!
		TBranch        *b_pf_METPhi;   //!
		TBranch        *b_pf_METsumEt;   //!
		TBranch        *b_pf_METmEtSig;   //!
		TBranch        *b_pf_METSig;   //!
		TBranch        *b_n_Ele;   //!
		TBranch        *b_ele_Charge;   //!
		TBranch        *b_ele_En;   //!
		TBranch        *b_ele_D0;   //!
		TBranch        *b_ele_Dz;   //!
		TBranch        *b_ele_Pt;   //!
		TBranch        *b_ele_Eta;   //!
		TBranch        *b_ele_Phi;   //!
		TBranch        *b_ele_R9;   //!
		TBranch        *b_ele_SCEta;   //!
		TBranch        *b_ele_SCPhi;   //!
		TBranch        *b_ele_HoverE;   //!
		TBranch        *b_ele_EoverP;   //!
		TBranch        *b_ele_EoverPout;   //!
		TBranch        *b_ele_EoverPInv;   //!
		TBranch        *b_ele_dEtaAtVtx;   //!
		TBranch        *b_ele_dPhiAtVtx;   //!
		TBranch        *b_ele_SigmaIEtaIEtaFull5x5;   //!
		TBranch        *b_ele_SigmaIPhiIPhiFull5x5;   //!
		TBranch        *b_ele_ConvVeto;   //!
		TBranch        *b_ele_MissHits;   //!
		TBranch        *b_ele_PFChIso;   //!
		TBranch        *b_ele_PFPhoIso;   //!
		TBranch        *b_ele_PFNeuIso;   //!
		TBranch        *b_ele_PFMiniIso;   //!
		TBranch        *b_ele_dEtaseedAtVtx;   //!
		TBranch        *b_N_Mu;   //!
		TBranch        *b_mu_Pt;   //!
		TBranch        *b_mu_En;   //!
		TBranch        *b_mu_Eta;   //!
		TBranch        *b_mu_Phi;   //!
		TBranch        *b_mu_Charge;   //!
		TBranch        *b_mu_IDbit;   //!
		TBranch        *b_mu_D0;   //!
		TBranch        *b_mu_Dz;   //!
		TBranch        *b_mu_Chi2NDF;   //!
		TBranch        *b_mu_InnerD0;   //!
		TBranch        *b_mu_InnerDz;   //!
		TBranch        *b_mu_InnervalidFraction;   //!
		TBranch        *b_mu_segmentCompatibility;   //!
		TBranch        *b_mu_chi2LocalPosition;   //!
		TBranch        *b_mu_trkKink;   //!
		TBranch        *b_mu_PFChIso;   //!
		TBranch        *b_mu_PFPhoIso;   //!
		TBranch        *b_mu_PFNeuIso;   //!
		TBranch        *b_mu_PFMiniIso;   //!
		TBranch        *b_mu_TrkLayers;
		TBranch        *b_mu_PixelLayers;
		TBranch        *b_mu_PixelHits;
		TBranch        *b_mu_MuonHits;
		TBranch        *b_mu_Stations;
		TBranch        *b_mu_Matches;
		TBranch        *b_mu_TrkQuality;
		TBranch        *b_mu_IsoTrk;   
		TBranch        *b_n_Jet;   //!
		TBranch        *b_jet_Pt;   //!
		TBranch        *b_jet_En;   //!
		TBranch        *b_jet_Eta;   //!
		TBranch        *b_jet_Phi;   //!
		TBranch        *b_jet_Area;   //!
		TBranch        *b_jet_Mt;   //!
		TBranch        *b_jet_CSV2BJetTags;   //!
		TBranch        *b_jet_JetProbabilityBJetTags;   //!
		TBranch        *b_jet_pfCombinedMVAV2BJetTags;   //!
		TBranch        *b_jet_DeepCSVTags_b;   //!
		TBranch        *b_jet_DeepCSVTags_bb;   //!
		TBranch        *b_jet_DeepCSVTags_c;   //!
		TBranch        *b_jet_DeepCSVTags_cc;   //!
		TBranch        *b_jet_DeepCSVTags_udsg;   //!
		TBranch        *b_jet_PartonID;   //!
		TBranch        *b_jet_HadFlvr;   //!
		TBranch        *b_jet_GenJetEn;   //!
		TBranch        *b_jet_GenJetPt;   //!
		TBranch        *b_jet_GenJetEta;   //!
		TBranch        *b_jet_GenJetPhi;   //!
		TBranch        *b_jet_GenPartonID;   //!
		TBranch        *b_jet_GenEn;   //!
		TBranch        *b_jet_GenPt;   //!
		TBranch        *b_jet_GenEta;   //!
		TBranch        *b_jet_GenPhi;   //!
		TBranch        *b_jet_GenPartonMomID;   //!
		TBranch        *b_jet_PFLooseId;   //!
		TBranch        *b_jet_ID;   //!
		TBranch        *b_jet_PUID;   //!
		TBranch        *b_jet_PUFullID;   //!
		TBranch        *b_jet_CHF;   //!
		TBranch        *b_jet_NHF;   //!
		TBranch        *b_jet_CEF;   //!
		TBranch        *b_jet_NEF;   //!
		TBranch        *b_jet_NCH;   //!
		TBranch        *b_jet_NNP;   //!
		TBranch        *b_jet_MUF;   //!
		TBranch        *b_N_AK8Jet;   //!
		TBranch        *b_AK8_JetPt;   //!
		TBranch        *b_AK8_JetEn;   //!
		TBranch        *b_AK8_JetEta;   //!
		TBranch        *b_AK8_JetPhi;   //!
		TBranch        *b_AK8_JetMass;   //!
		TBranch        *b_AK8_Jet_tau1;   //!
		TBranch        *b_AK8_Jet_tau2;   //!
		TBranch        *b_AK8_Jet_tau3;   //!
		TBranch        *b_AK8_Jet_tau4;   //!
		TBranch        *b_AK8_JetCHF;   //!
		TBranch        *b_AK8_JetNHF;   //!
		TBranch        *b_AK8_JetCEF;   //!
		TBranch        *b_AK8_JetNEF;   //!
		TBranch        *b_AK8_JetNCH;   //!
		TBranch        *b_AK8_JetNNP;   //!
		TBranch        *b_AK8_JetMUF;   //!
		TBranch        *b_AK8_Jetnconstituents;   //!
		TBranch        *b_AK8_JetPFLooseId;   //!
		TBranch        *b_AK8_JetPFTightLepVetoId;   //!
		TBranch        *b_AK8_JetSoftDropMass;   //!
		TBranch        *b_AK8_JetPrunedMass;   //!
		TBranch        *b_AK8_JetpfBoostedDSVBTag;   //!
		TBranch        *b_AK8_JetCSV;   //!
		TBranch        *b_AK8_JetDSVnewV4;   //!
		TBranch        *b_AK8_puppiPt;   //!
		TBranch        *b_AK8_puppiMass;   //!
		TBranch        *b_AK8_puppiEta;   //!
		TBranch        *b_AK8_puppiPhi;   //!
		TBranch        *b_AK8_puppiTau1;   //!
		TBranch        *b_AK8_puppiTau2;   //!
		TBranch        *b_AK8_puppiTau3;   //!
		TBranch        *b_AK8_puppiSDMass;   //!
		TBranch        *b_AK8_JetPartonID;   //!
		TBranch        *b_AK8_JetHadFlvr;   //!
		TBranch        *b_AK8_JetGenJetIndex;   //!
		TBranch        *b_AK8_JetGenJetEn;   //!
		TBranch        *b_AK8_JetGenJetPt;   //!
		TBranch        *b_AK8_JetGenJetEta;   //!
		TBranch        *b_AK8_JetGenJetPhi;   //!
		TBranch        *b_AK8_JetGenPartonID;   //!
		TBranch        *b_AK8_JetGenEn;   //!
		TBranch        *b_AK8_JetGenPt;   //!
		TBranch        *b_AK8_JetGenEta;   //!
		TBranch        *b_AK8_JetGenPhi;   //!
		TBranch        *b_AK8_JetGenPartonMomID;   //!
		TBranch        *b_n_AK8SDSJ;   //!
		TBranch        *b_AK8_SDSJPt;   //!
		TBranch        *b_AK8_SDSJEta;   //!
		TBranch        *b_AK8_SDSJPhi;   //!
		TBranch        *b_AK8_SDSJMass;   //!
		TBranch        *b_AK8_SDSJE;   //!
		TBranch        *b_AK8_SDSJCharge;   //!
		TBranch        *b_AK8_SDSJFlavour;   //!
		TBranch        *b_AK8_SDSJCSV;   //!
		TBranch        *b_n_AK8puppiSDSJ;   //!
		TBranch        *b_AK8_puppiSDSJPt;   //!
		TBranch        *b_AK8_puppiSDSJEta;   //!
		TBranch        *b_AK8_puppiSDSJPhi;   //!
		TBranch        *b_AK8_puppiSDSJMass;   //!
		TBranch        *b_AK8_puppiSDSJE;   //!
		TBranch        *b_AK8_puppiSDSJCharge;   //!
		TBranch        *b_AK8_puppiSDSJFlavour;   //!
		TBranch        *b_AK8_puppiSDSJCSV;   //!

		// lvb_qqqq_dCSV(TTree *tree=0);
		lvb_qqqq_dCSV(TString inputFile);
		virtual ~lvb_qqqq_dCSV();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		//  virtual void     Init(TTree *tree);
		virtual void     Init(TChain *tree);
		virtual void     Loop(TString OutputFileName);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		//===================User Defined Function List=========================
		virtual  void    Clear_Vector() ;
		virtual  float  delta_phi(float phi1,float phi2);

		// Histogram Defining & Filling Functions
		//1. For Gen Level Objects
		virtual void      DefineMC_NPtEta_Histo();
		virtual void      dRHisto_MCObject() ;
		virtual void      Fill_MC_PtEta_Histo(int id, int en);
		// virtual  float    Fill_mcMu_recoMu_plots(int Y, int Z, int lvl);  
		virtual void     dR_Plots_Genlvl() ;

		// 2. For Reco Objects  at Preselection
		virtual void      Define_NPtEta_Histo();  // for reco objects definitions at preselection level
		virtual void      dRHisto_RecoObject();  // for reco objects dR Plots  at preselection level
		virtual void      dRHisto_MCRecoObject();
		virtual void      Fill_NPtEta_Histo(int id, int en, int idx);
		virtual void      RecoPlots_dRHisto();
		virtual void       Fill_RecoObject();
		virtual float    Fill_MET_var(int lep, int lvl);  // X dlt
		//3. For Tag jets before Category Selection
		virtual void      Define_Tag_Jet_Histo() ;       // for  tag jets before Category Selection
		virtual void      dR_tagjetHisto();
		virtual void      Define_Reco_tagjetHisto() ; 
		virtual void       Fill_Puppi_jet(int jet, int idx);

		//---- Object Selection Functions------
		virtual  bool   Cut_Muon(int c_muon);
		virtual  void   Cut_bjet(int c_jet, int wp);
		virtual  void   Cut_AK8jet(int c_jet);

		// -----dR & dPt Calculating Functions--------
		//1. For MC Objects
		virtual  float   dR_mcbjet(int Y, int Z, int X, int idx, int lvl);
		virtual  float   dR_mc_mub(int Y, int Z, int A);
		virtual  void   dR_mc_bqi(int Y, int Z, int A);
		virtual  void   dR_mc_qiqj(int Y, int Z, int A);
		virtual  void   dR_mcAK8(int Y, int Z, int X, int idx);
		virtual  float   dPt_Gen_jetlep(int jet, int lep, int idx);

		//2. For Reco Objects
		virtual  float   dR_mu(int Y, int Z_jet, int idx, int lvl);
		virtual  float   dR_mu_AK8(int Y, int Z_jet, int idx, int lvl );
		virtual  float   dRPlots_bjet(int Y, int Z, int lvl) ;
		virtual  float   dR_AK8jet(int Y, int Z, int id8, int lvl);
		virtual  float   dR_AK8bjet(int Y, int Z, int id8, int idb, int lvl);
		virtual  void   dR_qjet_objects(int qjet, int obj, int idq, int idob);  
		virtual  float   dPt_lep(int jet, int lep, int idx, int lvl);

		//---- Jet Tagging Selections & Plots Functions------
		//1. Selections
		virtual  void    Higgs_selection() ;
		virtual  void    Top_selection() ;
		virtual  void    Wjet_selection() ;
		virtual  void    Fatjet_selection() ;
		// 2. Plot Function
		virtual  void     Wjet_Plots( int bjet );
		virtual  void     Topjet_Plots( ) ;
		virtual  void      Higgsjet_Plots( int bjet) ;
		virtual  void     Fatjet_Plots(int topbjet) ; 
		virtual void      TagJets_dRPlots() ;

		// ---------Signal Categorization , Histogram & Plot Function----------
		//1. Categorization
		virtual  void     Wtag0_Category(); 
		virtual  void     Wtag1_Category(); 
		// 2. Histogram definition
		virtual void      Category_Object_Histo() ;
		virtual void      Category_Object_MtHisto() ;  
		virtual void      Category_Object_dRHisto() ;  
		//3. Plot Functions  
		virtual  void     top_fatjet_Plots();
		virtual  void     top_Wjet_Plots() ; 
		virtual  void     W_fatjet_Plots() ;
		virtual  void     WW_lvbjet_Plots() ;
		virtual void      Higgs_lbjet_Plots() ;
		virtual void      CatI_Objects_Plots() ;  
		virtual void      CatII_Objects_Plots() ;    
		virtual void      CatIII_Objects_Plots() ;
		virtual void      CatIV_Objects_Plots() ;  
		virtual void      CatV_Objects_Plots() ;    
		virtual void      CatVI_Objects_Plots() ;
		virtual void      Category_Wjet_Plot( int cat, int W, int idx) ;
		virtual void      Category_Top_Plot( int cat,int top, int idx) ;
		virtual void      Category_Fat_Plot( int cat,int fat, int idx) ;
		virtual void      Category_Higgs_Plot( int cat, int higg, int idx) ;
		// 3. Transverse Mass Calculation   
		virtual  void    toptag_MTCalculation( ) ;
		virtual  void    Higgstag_MTCalculation( int topbjet)  ;
		virtual  void     genlvl_MTCalculation(int type ) ;

		//////////////
		virtual void      Define_Mt_Histo();  // for transverse mass histo only...delete it as soon for new histo

		// ===================create an array of Histograms======================
		// Category wise tag object Histogram

		std::array< std::array< TH1F*, 5> , 6> h_Histo_Pt;
		std::array< std::array< TH1F*, 5> , 6> h_Histo_Mt ;  
		std::array< std::array< TH1F*, 5> , 6> h_Histo_Eta ;                          
		std::array< std::array< TH1F*, 2> , 6> h_Histo_Mass ;  
		std::array< std::array< TH1F*, 2> , 6> h_Histo_SD ;  
		std::array< std::array< TH1F*, 2> , 6> h_Histo_tau ; 

		std::array< std::array< TH1F*, 6> , 6> h_Histo_dR ;   
		std::array< std::array< TH1F*, 6> , 6> h_Histo_dPt ;   

		// tag object histogram before Category selection
		std::array< TH1F*, 4>  h_tag_N ;  
		std::array< TH1F*, 12> h_tag_Pt ;  
		std::array< TH1F*, 12> h_tag_Eta ;                          
		std::array< TH1F*, 12> h_tag_Mass ;  
		std::array< TH1F*, 12> h_tag_SD ;  
		std::array< TH1F*, 12> h_tag_tau ;  
		std::array< std::array< TH1F*, 12> , 12> h_dR_tagjet ;

		// dR Histogram for tagjet - reco object
		std::array< TH1F*, 12>  h_dR_Recomu_tagjet ;
		std::array< TH1F*, 12>  h_dPt_lep_tagjet ;
		std::array< TH1F*, 12>  h_dR_Recob1_tagjet ;
		std::array< TH1F*, 12>  h_dR_Recob2_tagjet ;
		std::array< TH1F*, 12>  h_dR_Recob3_tagjet ;
		std::array< TH1F*, 12>  h_dR_Recob4_tagjet ;

		// mc level object histogram
		std::array< TH1F*, 7>  h_mcobject_pt;
		std::array< TH1F*, 7>  h_mcobject_eta;
		std::array< std::array< TH1F*, 6> , 6>  h_dR_MC ;
		std::array< std::array< TH1F*, 6> , 6>  h_dPt_MC ;
		std::array< std::array< TH1F*, 8> , 5>  h_dR_MCReco ;

		// reco object variables histogram
		std::array< TH1F*, 9>  h_object_pt ;
		std::array< TH1F*, 9>  h_object_eta ;
		std::array< TH1F*, 9>  h_object_no ;
		std::array< TH1F*, 3>  h_AK8_PUPPImass ;
		std::array< TH1F*, 3>  h_AK8_PUPPItau21 ; 
		std::array< TH1F*, 3>  h_AK8_PUPPItau32 ; 
		std::array< TH1F*, 3>  h_AK8_CHStau42 ; 
		std::array< std::array< TH1F*, 8> , 8>  h_dR_Reco ;
		std::array< std::array< TH1F*, 8> , 8>  h_dPt_Reco ;
		std::array< TH1F*, 3>  h_MET_var ;
		std::array< TH1F*, 9>  h_object_MT ;


		//======Global Variables =================

		int T_top = -1;
		int T_higgs = -1;
		int  Higgs_W = -1 ;
		int  eventI = -1 ;
		int  eventII = -1 ;
		int  eventIII = -1 ;
		int  eventIV = -1 ;
		int  eventV = -1 ;
		int  eventVI = -1 ;
		int  b_top = -1 ;

		vector<int>   puppi_jet ;
		vector <int> n_ele;
		vector <int> b_jet_loose;
		vector <int> b_jet;
		vector <int> n_Mu;
		vector <int> top_Mu;
		vector <int> Higgs_Mu;                
		vector <int> b_jet_tight;
		vector <int> n_jet;
		vector <int> n_forwjet;
		vector <int> n_AK8Jet;
		vector <int> b_jet_medium;
		vector <int> fat_jet;
		vector <int> Higgs_Jet;
		vector <int> Higgsjets;                
		vector <int> Higgs_qJet;
		vector <int> topW_q;
		vector <int> W_boson;
		vector <int> bjet_match ; 
		vector <int> topjet ;

		vector <int> CatI_Objects;              // for top, W, muon & MET
		vector <int> CatII_Objects;             // for top, fat, muon & MET
		vector <int> CatIII_Objects;            // for W, b, fat, muon & MET
		vector <int> CatIV_Objects;            // for  W, b, W, muon & MET
		vector <int> CatV_Objects;             // for  muon, MET, b, W & fat
		vector <int> CatVI_Objects;            // for Higgs, b, l & MET 
		//   map < int, string> CategoryII;            // for keeping track of selected elements.


};


#endif

#ifdef lvb_qqqq_new_cxx

/*lvb_qqqq_dCSV::lvb_qqqq_dCSV(TTree *tree) : fChain(0) 
  {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
cout << "Ayee babuai!! ";
if (tree == 0) {
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("WW_LH_M1500_Higgs_CHS.root");
if (!f || !f->IsOpen()) {
f = new TFile("WW_LH_M1500_Higgs_CHS.root");
}
f->GetObject("v_tree",tree);

}
Init(tree);
}
 */

lvb_qqqq_dCSV::lvb_qqqq_dCSV(TString inputFile)
{
	TChain *tree = new TChain("v_tree");
	//        add input files here


	ifstream datafile;
	datafile.open(inputFile.Data(), ifstream::in );

	TString datafilename;

	//               for signal
	string location ="/home/logan/VLQAnalysis/New_Method/";
	while (true) {
		datafile >> datafilename;
		string fname = datafilename.Data();
		string nameFile = fname;
		string string_search (".root");
		size_t found =nameFile.find(string_search);

		if(found != string::npos){
			tree->Add((location+fname).c_str());
			cout<<" FileName =" <<fname<<" , Events= "<<tree->GetEntries()<<"  file= "<<((location+fname).c_str())<<endl;
		}
		if( datafile.eof() ) break;
	}
	Init(tree);
}


lvb_qqqq_dCSV::~lvb_qqqq_dCSV()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t lvb_qqqq_dCSV::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t lvb_qqqq_dCSV::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

//====================User Defined Functions' Definitions====================

void  lvb_qqqq_dCSV::Clear_Vector()
{
	//   for clearing objects vector
	puppi_jet.clear();
	n_ele.clear();
	n_Mu.clear();
	top_Mu.clear();
	topjet.clear();        
	Higgs_Mu.clear();        
	n_jet.clear();
	n_forwjet.clear();
	b_jet_loose.clear();
	b_jet.clear();
	b_jet_medium.clear();
	b_jet_tight.clear();
	fat_jet.clear();
	Higgs_Jet.clear();
	Higgsjets.clear();        
	Higgs_qJet.clear();
	n_AK8Jet.clear();
	topW_q.clear();
	W_boson.clear();
	bjet_match.clear();
	CatI_Objects.clear() ;    
	CatII_Objects.clear() ;
	CatIII_Objects.clear() ;
	CatIV_Objects.clear() ;
	CatV_Objects.clear() ;
	CatVI_Objects.clear() ;        
	//     CategoryII.clear();             
}


//========== Histogram Definations======================================

void lvb_qqqq_dCSV::Define_Mt_Histo()
{
	//   New defination of histogram
	char  dr_MTname[100] , dr_MTtitle[100] ;

	int s = -1;
	int rk = -1 ;
	string MT_comp[9] = {"Mt_muMeT", "Mt_WmuMeT", "Mt_Tprime", "Mt_topjet", "Mt_Wjet", "Mt_Wmc", "MT_HiggsMC", "MT_topmc", "MT_Tmc"} ;
	//cout<<"\nFunc 1" ;
	// for Transverse mass plots, we have 6 combinations for objects
	for( int Sk = 0 ; Sk < 9 ; Sk ++ )
	{
		sprintf( dr_MTname, "%s",  MT_comp[Sk].c_str() ) ;
		sprintf( dr_MTtitle, "%s Distribution",  MT_comp[Sk].c_str() ) ;		
		h_object_MT.at(Sk) = new TH1F(dr_MTname,dr_MTtitle, 300, 0.0 , 3000.0 );
		h_object_MT.at(Sk)->GetYaxis()->SetTitle("Events");
		//cout<< " , h_object_Mt["<<Sk<<"].name() : "<< h_object_MT.at(Sk) ->GetName()<<endl; 
	}      
}

// ---------------------------------------------------mc level objects ----------------------------------------------

void lvb_qqqq_dCSV::DefineMC_NPtEta_Histo()
{
	//   New defination of histogram
	//        TString   mc_pT , mc_pT_title, gh_name, gh_title;
	char mc_pT[100], mc_pT_title[100], gh_name[100], gh_title[100] ;

	std::string mc_objectI[7] = {"topb", "topWl","higgWq1", "higgWq2", "higgWl", "assob", "forwq" };
	//std::string mc_objectI[3] = {"topb", "topWl","higgWq1"};
	// for object pT histograms
	//cout << "\nFunc 2" ;
	for(Int_t k = 0 ; k < 7; ++k) {
		sprintf (mc_pT, "Pt(%s)", mc_objectI[k].c_str() )  ;
		sprintf (mc_pT_title, "Pt(%s) Distribution", mc_objectI[k].c_str() )  ;               
		h_mcobject_pt.at(k) = new TH1F(mc_pT,mc_pT_title, 200, 0, 1000.0);
		h_mcobject_pt.at(k) ->GetYaxis()->SetTitle("Events");        
		//	cout<< " , h_mcobject_pt["<<k<<"].name() : "<< h_mcobject_pt.at(k) ->GetName()<<endl; 	
	}


	for(Int_t h = 0 ; h < 7; ++h) {
		sprintf (gh_name, "Eta(%s)", mc_objectI[h].c_str() )  ;
		sprintf (gh_title, "Eta(%s) Distribution", mc_objectI[h].c_str() )  ;                            
		h_mcobject_eta.at(h) = new  TH1F(gh_name,gh_title, 200, -5.0, 5.0);
		h_mcobject_eta.at(h) ->GetYaxis()->SetTitle("Events");
		//cout<<" , h_mcobject_eta["<<h<<"].name() : "<< h_mcobject_eta.at(h) ->GetName()<<endl; 
	}

}



void lvb_qqqq_dCSV::dRHisto_MCObject()
{
	// New defination of histogram
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string mc_object[7] = {"topb", "topWl","higgWq1", "higgWq2", "higgWl", "assob", "forwq" };
	// Here, i = 0-1 top decayed, i = 2-4 Higgs decayed, i = 5 assob, i = 6 forwq
	Int_t s = -1;
	//cout << "\nFunc 3" ;	
	for (Int_t i = 0; i < 6; i++) {     
		for (Int_t k = i+1; k < 6; k++) {
			// for mc  objects dR histograms
			sprintf(dR_name, "DeltaR(%s_%s)", mc_object[i].c_str(), mc_object[k].c_str()) ;
			sprintf(dR_title, "DeltaR(%s_%s) Distribution", mc_object[i].c_str(), mc_object[k].c_str()) ;
			h_dR_MC.at(i).at(k) = new TH1F(dR_name,dR_title, 500, 0, 5.0);
			h_dR_MC.at(i).at(k) ->GetYaxis()->SetTitle("Events");			
			////cout<<" , h_dR_MC["<<i<<"]["<<k<<"].name() : "<< h_dR_MC.at(i).at(k) ->GetName()<<endl; 
			// for mc Object dPt histograms     
			sprintf(dPt_name, "DeltaPt(%s_%s)", mc_object[i].c_str(), mc_object[k].c_str()) ;
			sprintf(dPt_title, "DeltaPt(%s_%s) Distribution", mc_object[i].c_str(), mc_object[k].c_str()) ;
			h_dPt_MC.at(i).at(k) = new TH1F(dPt_name,dPt_title, 250, 0, 500.0);
			h_dPt_MC.at(i).at(k) ->GetYaxis()->SetTitle("Events");			
			//cout<<" h_dPt_MC["<<i<<"]["<<k<<"].name() : "<< h_dPt_MC.at(i).at(k) ->GetName()<<endl;			
		}
	}
}     


//------------------------------------------reco level objects - Preselection-----------------------------------------------

void lvb_qqqq_dCSV::Define_NPtEta_Histo()
{
	//   New defination of histogram
	char   pT_name[100] , pT_title[100] ; 
	string object_no[5] = {"mu", "bjet", "AK8jet","qjet","forwjet"};
	string reco_object[9] = {"MET", "mu1", "mu2","bjet1", "bjet2", "bjet3","AK8jet1", "AK8jet2", "AK8jet3"};
	string variable[6] = {"Pt", "Eta", "Puppimass", "Tau21", "Tau32", "Tau42"} ;
	string ar;
	//cout << "\nFunc 4" ;
	int s = -1;

	// for object pT histograms
	for(Int_t k = 0 ; k < 9; k++) {
		for(Int_t j = 0; j < 6; j++) {

			if( k < 6 && j > 1) continue ;
			sprintf( pT_name , "%s(%s)" , variable[j].c_str(), reco_object[k].c_str() ) ;
			sprintf( pT_title , "%s(%s) Distribution" , variable[j].c_str(), reco_object[k].c_str() ) ;			

			if ( j == 0) {
				h_object_pt.at(k) = new TH1F(pT_name,pT_title, 200, 0, 1000.0);
				h_object_pt.at(k) ->GetYaxis()->SetTitle("Events");
				//cout<< " h_object_pt["<<k<<"].name() : "<< h_object_pt.at(k) ->GetName()<<endl; 
			}
			if ( j == 1) {
				h_object_eta.at(k) = new TH1F(pT_name,pT_title, 200, -5.0, 5.0);
				h_object_eta.at(k) ->GetYaxis()->SetTitle("Events");
				//cout<< " h_object_eta["<<k<<"].name() : "<< h_object_eta.at(k) ->GetName()<<endl; 	
			}
			if( j == 2) {
				s = k - 3*j ;
				h_AK8_PUPPImass.at(s) =  new TH1F(pT_name,pT_title, 100, 0.0, 1000.0) ;
				h_AK8_PUPPImass.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPImass["<<s<<"].name() : "<< h_AK8_PUPPImass.at(s)->GetName()<<endl; 	
			}
			if( j == 3) {
				s = k - 2*j ;
				h_AK8_PUPPItau21.at(s) =  new TH1F(pT_name,pT_title, 100, 0.0, 1.0) ;
				h_AK8_PUPPItau21.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau21["<<s<<"].name() : "<< h_AK8_PUPPItau21.at(s)->GetName()<<endl; 					
			}
			if( j == 4) {
				s = k - j - 2 ;
				h_AK8_PUPPItau32.at(s) =  new TH1F(pT_name,pT_title, 100, 0.0, 1.0) ;
				h_AK8_PUPPItau32.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_PUPPItau32["<<s<<"].name() : "<< h_AK8_PUPPItau32.at(s)->GetName()<<endl; 									
			}
			if( j == 5) {
				s  = k - j -1;
				h_AK8_CHStau42.at(s) =  new TH1F(pT_name,pT_title, 100, 0.0, 1.0) ;
				h_AK8_CHStau42.at(s) -> GetYaxis()->SetTitle("Events");
				//cout<<"h_AK8_CHStau42["<<s<<"].name() : "<< h_AK8_CHStau42.at(s)->GetName()<<endl; 									
			}

		}
	}

	// for No. object  histograms
	for(Int_t j = 0; j < 5; j++) {	
		sprintf( pT_name , "#N(%s)" , object_no[j].c_str() ) ;
		sprintf( pT_title , "#N(%s) Distribution" , object_no[j].c_str() ) ;	
		h_object_no.at(j) = new TH1F(pT_name,pT_title, 10, 0, 10);  
		h_object_no.at(j) ->GetYaxis()->SetTitle("Events");
		//cout<<"h_object_no["<<j<<"].name() : "<< h_object_no.at(j)->GetName()<<endl; 			
	}

}


void lvb_qqqq_dCSV::dRHisto_RecoObject() 
{
	// New defination of histogram for reco level at Preselection 
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string reco_object[8] = {"mu0", "mu1", "bjet1", "bjet2", "bjet3", "AK8jet1", "AK8jet2", "AK8jet3" };
	// Here i -> 0-1 muon, i -> 2-4 bjet & i -> 5-7 AK8jet
	//cout<<"\nFunc 5" ;	
	for (Int_t i = 0; i < 8; i++) {          
		for (Int_t k = i+1; k < 8; k++) {      

			// for reco preselected objects dR histograms
			sprintf( dR_name , "DeltaR(%s_%s)" , reco_object[i].c_str(), reco_object[k].c_str() ) ;
			sprintf( dR_title ,"DeltaR(%s_%s) Diatribution",reco_object[i].c_str(),reco_object[k].c_str()) ;
			h_dR_Reco.at(i).at(k) = new TH1F(dR_name,dR_title, 500, 0, 5.0);
			h_dR_Reco.at(i).at(k) ->GetYaxis()->SetTitle("Events");
			//cout<<"h_dR_Reco[" <<i<<"]["<<k<<"].name() : "<< h_dR_Reco.at(i).at(k)->GetName()<<endl; 									


			// for lep-bjet dPt histograms  
			sprintf( dPt_name , "DeltaPt(%s_%s)" , reco_object[i].c_str(), reco_object[k].c_str() ) ;
			sprintf( dPt_title ,"DeltaPt(%s_%s) Diatribution",reco_object[i].c_str(),reco_object[k].c_str());
			h_dPt_Reco.at(i).at(k) = new TH1F(dPt_name,dPt_title, 250, 0, 500.0);
			h_dPt_Reco.at(i).at(k) ->GetYaxis()->SetTitle("Events");
			//cout<<"h_dPt_Reco[" <<i<<"]["<<k<<"].name() : "<< h_dPt_Reco.at(i).at(k)->GetName()<<endl; 			
		}
	}
}     

//------------------------------------------mc-Reco dR Histogram--------------------------------------

void lvb_qqqq_dCSV::dRHisto_MCRecoObject()
{
	// New defination of histogram
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string mc_object[5] = {"topb", "topWl","higgWq1", "higgWq2", "higgWl" };
	string reco_object[8] = {"mu0", "mu1", "bjet1", "bjet2", "bjet3", "AK8jet1", "AK8jet2", "AK8jet3" };
	//cout<<"\nFunc 6" ;
	for (Int_t i = 0; i < 5; i++) {          
		for (Int_t k = 0; k < 8; k++) {
			// for mc  objects dR histograms
			sprintf( dR_name, "DeltaR(%s_%s)" , mc_object[i].c_str(), reco_object[k].c_str()) ;
			sprintf( dR_title,"DeltaR(%s_%s) Distribution" , mc_object[i].c_str(), reco_object[k].c_str()) ;
			h_dR_MCReco.at(i).at(k) = new TH1F(dR_name,dR_title, 500, 0, 5.0);
			h_dR_MCReco.at(i).at(k) ->GetYaxis()->SetTitle("Events");
			//cout<<"h_dR_MCReco[" <<i<<"]["<<k<<"].name() : "<< h_dR_MCReco.at(i).at(k)->GetName()<<endl; 			
		}
	}     
}
//-----------------------------Tagged jets Histogram ---------------------------------------------------

void  lvb_qqqq_dCSV::Define_Tag_Jet_Histo()
{
	char tag_name[100] , tag_title[100] ;
	int x1 = -1 ;
	int rk = -1 ;
	string jet_tag[4] = { "Wjet", "topjet", "Higgsjet", "fatjet" } ; 
	string Wjet_var[5] = { "Pt", "Eta", "TransMass", "PuppiSDMass",  "PUPPItau21"} ;
	string top_var[5] = { "Pt", "Eta", "TransMass", "PuppiSDMass",  "PUPPItau32"} ;
	string Higgs_var[5] = { "Pt", "Eta", "TransMass", "CHSPrunedMass",  "CHStau42"} ;
	//cout<<"\nFunc 7" ;
	for( int h =0 ; h < 4 ; h ++) {       
		sprintf( tag_name, "N(%s)" , jet_tag[h].c_str() ) ;
		sprintf( tag_title, "N(%s) Distribution" , jet_tag[h].c_str() ) ;
		h_tag_N.at(h)  =  new TH1F(tag_name,tag_title, 10, 0, 10);
		h_tag_N.at(h)  -> GetYaxis()->SetTitle("Events");
		//cout<<"h_tag_N[" <<h<<"].name() : "<< h_tag_N.at(h)->GetName()<<endl; 		
	}

	// Here i = 0-3 Wjets, i -> 4-5 topjet, i = 6-7 Higgsjet, i = 8-11 fatjet,Use this index cycle for filling plots  
	for(     int Sj = 0 ; Sj < 12; Sj ++) {
		for( int Tk = 0 ; Tk < 5 ; Tk ++ ) {

			if ( Sj < 4  ) {
				x1 = Sj + 1 ;                // for Wjets index   0-3     
				sprintf( tag_name, "%s(%s%d)",  Wjet_var[Tk].c_str(), jet_tag[0].c_str(), x1 ) ;  
				sprintf( tag_title,"%s(%s%d) Distribution", Wjet_var[Tk].c_str(), jet_tag[0].c_str(), x1);  				    
			}
			if ( Sj == 4 || Sj == 5 ) {
				x1 = Sj - 3 ;           // for topjet index  4-5
				sprintf( tag_name, "%s(%s%d)",  top_var[Tk].c_str(), jet_tag[1].c_str(), x1 ) ;  
				sprintf( tag_title,"%s(%s%d) Distribution", top_var[Tk].c_str(), jet_tag[1].c_str(), x1); 
			}
			if ( Sj == 6 || Sj == 7 ) {
				x1 = Sj - 5 ;           // for Higgsjet index   6-7 
				sprintf( tag_name, "%s(%s%d)",  Higgs_var[Tk].c_str(), jet_tag[2].c_str(), x1 ) ;     	
				sprintf( tag_title,"%s(%s%d) Distribution", Higgs_var[Tk].c_str(), jet_tag[2].c_str(), x1); 	
			}          
			if ( Sj > 7 )      {
				x1 = Sj - 7 ;            // for fatjet index 8-11
				sprintf( tag_name, "%s(%s%d)",  Wjet_var[Tk].c_str(), jet_tag[3].c_str(), x1 ) ;  		
				sprintf( tag_title,"%s(%s%d) Distribution", Wjet_var[Tk].c_str(), jet_tag[3].c_str(), x1); 	
			}         


			if(Tk == 0 ) {
				h_tag_Pt.at(Sj)             = new TH1F(tag_name,tag_title, 400, 0, 1000.0);
				h_tag_Pt.at(Sj)             ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_Pt[" <<Sj<<"].name() : "<< h_tag_Pt.at(Sj)->GetName()<<endl; 						
			}      
			if(Tk == 1 )  {
				h_tag_Eta.at(Sj)           = new TH1F(tag_name,tag_title, 200, -5.0, 5.0);                   
				h_tag_Eta.at(Sj)            ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_Eta[" <<Sj<<"].name() : "<< h_tag_Eta.at(Sj)->GetName()<<endl; 			
			}      
			if(Tk == 2 ) {
				h_tag_Mass.at(Sj)         = new TH1F(tag_name,tag_title, 300, 0, 1000.0);         
				h_tag_Mass.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_Mass[" <<Sj<<"].name() : "<< h_tag_Mass.at(Sj)->GetName()<<endl; 						
			}                                
			if(Tk == 3 )  {
				h_tag_SD.at(Sj)      = new TH1F(tag_name,tag_title, 300, 0, 1000.0);       
				h_tag_SD.at(Sj)       ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_SD[" <<Sj<<"].name() : "<< h_tag_SD.at(Sj)->GetName()<<endl;				
			}            
			if(Tk == 4 )  {
				h_tag_tau.at(Sj)         = new TH1F(tag_name,tag_title, 50, 0, 1.0);               
				h_tag_tau.at(Sj)         ->GetYaxis()->SetTitle("Events");
				//cout<<"h_tag_tau[" <<Sj<<"].name() : "<< h_tag_tau.at(Sj)->GetName()<<endl;				
			}                        

		}      
	}
}


void lvb_qqqq_dCSV::dR_tagjetHisto()
{
	// New defination of histogram
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string tag_jet[12] = {"Wjet1", "Wjet2", "Wjet3", "Wjet4", "topjet1", "topjet2", "Higgsjet1", "Higgsjet2", "fatjet1", "fatjet2", "fatjet3", "fatjet4"} ;                   
	//for tag jet   dR histograms
	// Here i = 0-3 Wjets, i -> 4-5 topjet, i = 6-7 Higgsjet, i = 8-11 fatjet, Use this index cycle for illing plots
	//cout<<"\nFunc 8" ;
	for (int i = 0; i < 12; i++) {
		for(int j = i+1; j < 12; j++) {

			sprintf(dR_name,"DeltaR(%s_%s)" , tag_jet[i].c_str(), tag_jet[j].c_str() ) ;
			sprintf(dR_title,"DeltaR(%s_%s) Distribution", tag_jet[i].c_str(), tag_jet[j].c_str() ) ;	
			h_dR_tagjet.at(i).at(j)   =    new TH1F(dR_name,dR_title, 500, 0, 5.0);
			h_dR_tagjet.at(i).at(j)   ->   GetYaxis()->SetTitle("Events");         
			//cout<<"h_dR_tagjet["<<i<<"]["<<j<<"].name() : "<< h_dR_tagjet.at(i).at(j)->GetName()<<endl;			
		}
	}
}


void lvb_qqqq_dCSV::Define_Reco_tagjetHisto()
{
	// New defination of histogram
	char dR_name[100], dR_title[100], dPt_name[100], dPt_title[100] ;
	string reco_object[5]  =  {"muon","bjet1", "bjet2", "bjet3", "bjet4"} ;
	string tag_jet[12] = {"Wjet1", "Wjet2", "Wjet3", "Wjet4", "topjet1", "topjet2", "Higgsjet1", "Higgsjet2", "fatjet1", "fatjet2", "fatjet3", "fatjet4"} ;  

	//for tag jet - muon  dR histograms
	// Here i = 0-3 Wjets, i -> 4-5 topjet, i = 6-7 Higgsjet, i = 8-11 fatjet, Use this index cycle for illing plots	
	//cout<<"\nFunc 9" ;	
	int s = -1 ;
	for (Int_t i = 0; i < 5; i++) {
		for(Int_t j = 0; j < 12; j++) {
			//cout << "\n value i = " << i << " && value j << " << j ;        

			if( i == 0){			
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str()) ;
				h_dR_Recomu_tagjet.at(j)= new TH1F(dR_name,dR_title, 500, 0, 5.0);
				h_dR_Recomu_tagjet.at(j)->GetYaxis()->SetTitle("Events");
				//cout<<"h_dR_Recomu_tagjet]"<<j<<"].name() : "<< h_dR_Recomu_tagjet.at(j)->GetName()<<endl;					

				sprintf( dPt_name,"DeltaPt(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dPt_title,"DeltaPt(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str());
				h_dPt_lep_tagjet.at(j)= new TH1F(dPt_name,dPt_title, 250, 0, 500.0);
				h_dPt_lep_tagjet.at(j)->GetYaxis()->SetTitle("Events");				
				//cout<<"h_dPt_lep_tagjet]"<<j<<"].name() : "<< h_dPt_lep_tagjet.at(j)->GetName()<<endl;						
			}
			if( i == 1){
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str());	
				h_dR_Recob1_tagjet.at(j)= new TH1F(dR_name,dR_title, 500, 0, 5.0);
				h_dR_Recob1_tagjet.at(j)->GetYaxis()->SetTitle("Events");      
				//cout<<"h_dR_Recob1_tagjet]"<<j<<"].name() : "<< h_dR_Recob1_tagjet.at(j)->GetName()<<endl;										
			}
			if( i == 2){
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str() );
				h_dR_Recob2_tagjet.at(j)= new TH1F(dR_name,dR_title, 500, 0, 5.0);
				h_dR_Recob2_tagjet.at(j)->GetYaxis()->SetTitle("Events");
				//cout<<"h_dR_Recob2_tagjet]"<<j<<"].name() : "<< h_dR_Recob2_tagjet.at(j)->GetName()<<endl;				
			}
			if( i == 3){
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str() );
				h_dR_Recob3_tagjet.at(j)= new TH1F(dR_name,dR_title, 500, 0, 5.0);
				h_dR_Recob3_tagjet.at(j)->GetYaxis()->SetTitle("Events");
				//cout<<"h_dR_Recob3_tagjet]"<<j<<"].name() : "<< h_dR_Recob3_tagjet.at(j)->GetName()<<endl;				
			}
			if( i == 4){
				sprintf( dR_name,"DeltaR(%s_%s)",  reco_object[i].c_str(), tag_jet[j].c_str() ) ;
				sprintf( dR_title,"DeltaR(%s_%s) Distribution", reco_object[i].c_str(), tag_jet[j].c_str() );
				h_dR_Recob4_tagjet.at(j)= new TH1F(dR_name,dR_title, 500, 0, 5.0);
				h_dR_Recob4_tagjet.at(j)->GetYaxis()->SetTitle("Events");
				//cout<<"h_dR_Recob4_tagjet]"<<j<<"].name() : "<< h_dR_Recob4_tagjet.at(j)->GetName()<<endl;
			}
		}//END OF J
	}// END OF I
	//	cout << "Size in loop = " <<  sizeof(h_dR_Recomu_tagjet) <<endl;
} // END OF FUNC


//================Signal Category Objects Histogram==================


void lvb_qqqq_dCSV::Category_Object_Histo()
{
	char Histo_name[100], Histo_title[100] ;
	string  obj[5], var  ;
	string  reco_var[2] = { "Pt", "Eta" } ;
	string Wjet_var[5] = { "Pt", "Eta", "TransMass", "PuppiSDMass",  "PUPPItau21"} ;
	string top_var[5] = { "Pt", "Eta", "TransMass", "PuppiSDMass",  "PUPPItau32"} ;
	string Higgs_var[5] = { "Pt", "Eta", "TransMass", "CHSPrunedMass",  "CHStau42"} ;
	int x ;
	int f = 0 ;
	//cout<<"\nFunc 10" ;	
	//For category  plots
	for ( int u = 0 ; u < 6 ; u ++ ) 
	{ 
		x = u +1 ;
		// Category Decision
		if ( u == 0 ) {
			obj[0] = "Topjet" ; obj[1] = "HWjet"; obj[2] = "muon" ; obj[3] = "MET";
			f = 4;
		}
		if ( u == 1 )  {
			obj[0] = "Topjet" ; obj[1] = "HFatjet"; obj[2] = "muon" ; obj[3] = "MET";
			f =4 ;
		}
		if ( u == 2 )  {
			obj[0] = "TWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 3 )  {
			obj[0] = "TWjet" ; obj[1] = "HWjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 4 )  {
			obj[0] = "HWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 5 ) {
			obj[0] = "no" ; obj[1] = "Higgsjet" ; obj[2] = "bjet"; obj[3] = "muon" ; obj[4] = "MET";
			f = 5;
		}

		for ( int j = 0; j < f; j ++ )
		{
			if ( obj[j] == "no" ) continue;
			for ( int k = 0; k < 5 ; k++) 
			{
				if ( j > 1 && k > 1) continue ;

				if( obj[j] == "Topjet" )                                var =  top_var[k]  ;
				if(obj[j] == "Higgsjet")                                var =  Higgs_var[k]  ;
				if( obj[j] == "muon" )                                 var =  reco_var[k]  ;
				if( obj[j] == "MET" )                                 var =  reco_var[k]  ;
				if( obj[j] == "bjet" )                                   var =  reco_var[k] ;
				if( obj[j] == "TWjet" || obj[j] == "HWjet" )		var  =  Wjet_var[k]  ;				
				if( obj[j] == "TFatjet" || obj[j] == "HFatjet")     var  =  Wjet_var[k]  ;				

				sprintf( Histo_name, "Category%d_%s(%s)", x, obj[j].c_str(), var.c_str() ) ; 
				sprintf( Histo_title, "Category%d_%s(%s) Distribution", x, obj[j].c_str(), var.c_str() ) ; 				
				if(k == 0 ) {
					h_Histo_Pt.at(u).at(j)       = new TH1F(Histo_name,Histo_title, 400, 0, 1000.0);
					h_Histo_Pt.at(u).at(j)             ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_Pt["<<u<<"]["<<j<<"].name() : "<< h_Histo_Pt.at(u).at(j)->GetName()<<endl;     				
				}      
				if(k == 1 )  {
					h_Histo_Eta.at(u).at(j)      = new TH1F(Histo_name,Histo_title, 200, -5.0, 5.0);    
					h_Histo_Eta.at(u).at(j)      ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_Eta["<<u<<"]["<<j<<"].name() : "<< h_Histo_Eta.at(u).at(j)->GetName()<<endl;					
				}      
				if(k == 2 ) {
					h_Histo_Mass.at(u).at(j)     = new TH1F(Histo_name,Histo_title, 300, 0, 1000.0);   
					h_Histo_Mass.at(u).at(j)      ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_Mass["<<u<<"]["<<j<<"].name() : "<< h_Histo_Mass.at(u).at(j)->GetName()<<endl;						
				}                                
				if(k == 3 )  {
					h_Histo_SD.at(u).at(j)        = new TH1F(Histo_name,Histo_title, 300, 0, 1000.0);
					h_Histo_SD.at(u).at(j)        ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_SD["<<u<<"]["<<j<<"].name() : "<< h_Histo_SD.at(u).at(j)->GetName()<<endl;						
				}            
				if(k == 4 )  {
					h_Histo_tau.at(u).at(j)         = new TH1F(Histo_name,Histo_title, 50, 0, 1.0);  
					h_Histo_tau.at(u).at(j)         ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_tau["<<u<<"]["<<j<<"].name() : "<< h_Histo_tau.at(u).at(j)->GetName()<<endl;
				}  // end of k ==4                       

			}  //end of kth loop
		}    //end of jth loop

	}    // end of uth loop
}  //  end of function


void lvb_qqqq_dCSV::Category_Object_MtHisto()
{
	char Histo_name[100], Histo_title[100];
	string TransM ;
	string  obj[5]  ;
	int x ;
	int f = 0 ;
	//cout<<"\nFunc 11" ;		
	//For category  plots
	for ( int u = 0 ; u < 6 ; u ++ ) 
	{ 
		x = u +1 ;
		// Category Decision
		if ( u == 0 ) {
			obj[0] = "Top" ; obj[1] = "HW"; obj[2] = "mu" ; obj[3] = "MET";
			f = 4;
		}
		if ( u == 1 )  {
			obj[0] = "Top" ; obj[1] = "HFat"; obj[2] = "mu" ; obj[3] = "MET";
			f =4 ;
		}
		if ( u == 2 )  {
			obj[0] = "TW" ; obj[1] = "b"; obj[2] = "HFat" ; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 3 )  {
			obj[0] = "TW" ; obj[1] = "b"; obj[2] = "HW" ; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 4 )  {
			obj[0] = "HW" ; obj[1] = "HFat"; obj[2] = "b" ; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}
		if ( u == 5 ) {
			obj[0] = "no" ; obj[1] = "Higgs" ; obj[2] = "b"; obj[3] = "mu" ; obj[4] = "MET";
			f = 5;
		}

		TransM = TransM + "MET" ;
		for ( int j = f-2; j > -1; j -- )
		{
			if ( obj[j] == "no" )                                continue;    
			if(  obj[j] == "b" && ( u == 2 || u == 3)  )  continue; 
			if( obj[j] == "HFat" && u == 4 )               continue ;     
			if( obj[j] == "Top" )                               TransM = TransM + obj[j] ;
			if(obj[j] == "Higgs")                               TransM = TransM + obj[j] ;
			if( obj[j] == "mu" )                                TransM = TransM + obj[j] ;   
			if(  obj[j] == "b" )                                 TransM = TransM + obj[j] ;
			if(  obj[j] == "HFat" )                             TransM = TransM + obj[j] ;     
			if( obj[j] == "TW" || obj[j] == "HW" ) {
				if ( u == 2 || u == 3 || u == 4 )   TransM = TransM + obj[j+1] + obj[j];
				else {
					TransM = TransM + obj[j] ;
				}
			}

			sprintf( Histo_name, "Cat%d_Mt(%s)", x , TransM.c_str() ) ;
			sprintf( Histo_title, "Cat%d_Mt(%s) Distribution", x , TransM.c_str() ) ;			
			h_Histo_Mt.at(u).at(j)         = new TH1F(Histo_name,Histo_title, 500, 0, 1000.0);
			h_Histo_Mt.at(u).at(j)         ->GetYaxis()->SetTitle("Events");   
			//cout<<"h_Histo_Mt["<<u<<"]["<<j<<"].name() : "<< h_Histo_Mt.at(u).at(j)->GetName()<<endl;			
		} 
		TransM = "" ;
	}
}


void lvb_qqqq_dCSV::Category_Object_dRHisto()
{
	char Histo_name[100], Histo_title[100] , dp_name[100], dp_title[100];
	string  obj[5]  ;
	int x ;
	int f = 0 , g  ;
	//cout<<"\nFunc 12" ;		
	//For category  plots
	for ( int u = 0 ; u < 6 ; u ++ ) 
	{ 
		g = -1 ; 
		x = u +1 ;
		// Category Decision
		if ( u == 0 ) {
			obj[0] = "Topjet" ; obj[1] = "HWjet"; obj[2] = "muon" ; obj[3] = "MET";
			f = 3;
		}
		if ( u == 1 )  {
			obj[0] = "Topjet" ; obj[1] = "HFatjet"; obj[2] = "muon" ; obj[3] = "MET";
			f =3 ;
		}
		if ( u == 2 )  {
			obj[0] = "TWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}
		if ( u == 3 )  {
			obj[0] = "TWjet" ; obj[1] = "HWjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}
		if ( u == 4 )  {
			obj[0] = "HWjet" ; obj[1] = "HFatjet"; obj[2] = "bjet" ; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}
		if ( u == 5 ) {
			obj[0] = "no" ; obj[1] = "Higgsjet" ; obj[2] = "bjet"; obj[3] = "muon" ; obj[4] = "MET";
			f = 4;
		}

		for ( int j = 0; j < f; j ++ )
		{
			if ( obj[j] == "no" ) continue;
			for ( int k = j+1; k < f ; k++) 
			{
				g ++ ;
				if( obj[j] == "Topjet" ){
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if(obj[j] == "Higgsjet") {
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if( obj[j] == "bjet" ) {
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if( obj[j] == "TWjet" || obj[j] == "HWjet" ){
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				}
				if( obj[j] == "TFatjet" || obj[j] == "HFatjet"){
					sprintf( Histo_name, "Category%d_dR(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;
					sprintf( dp_name, "Category%d_dPt(%s_%s)",x ,obj[j].c_str(), obj[k].c_str() ) ;		
				} 

				sprintf(Histo_title,"Category%d_dR(%s_%s) Distribution",x ,obj[j].c_str(), obj[k].c_str() );
				sprintf(dp_title, "Category%d_dPt(%s_%s) Distribution",x ,obj[j].c_str(), obj[k].c_str() ) ;

				h_Histo_dR.at(u).at(g)             = new TH1F(Histo_name,Histo_title, 500, 0.0, 5.0);
				h_Histo_dR.at(u).at(g)             ->GetYaxis()->SetTitle("Events");
				//cout<<"h_Histo_dR["<<u<<"]["<<g<<"].name() : "<< h_Histo_dR.at(u).at(g)->GetName()<<endl;				

				if( obj[k] == "muon"){
					h_Histo_dPt.at(u).at(g)        = new TH1F(dp_name,dp_title, 500, 0.0, 500.0);
					h_Histo_dPt.at(u).at(g)        ->GetYaxis()->SetTitle("Events");
					//cout<<"h_Histo_dPt["<<u<<"]["<<g<<"].name() : "<< h_Histo_dPt.at(u).at(g)->GetName()<<endl;					
				}

			} // end of kth loop 
		} // end of jth loop
	}  //  end of uth loop
}    


//===============================================================


void  lvb_qqqq_dCSV::dR_Plots_Genlvl()
{
	//  for pt & eta topb= 0, topWl = 1, higgWq1 = 2, higgWq2= 3, higgWl=4, assob = 5, forwq = 6
	//  for dR , Here, i = 0-1 top decayed, i = 2-4 Higgs decayed, i = 5 assob, i = 6 forwq

	// int it = topW_q[0] ;
	//  int it = T_top ;
	//  int jet1 = topW_q[0] ;
	//   int j =  Higgs_Mu[0] ;    

	int it = 0 ;  
	int jet1 = 1  ;      
	int j =  2 ;

	int A = -1 ; // for mom ID
	/*    
	      h_mcobject_pt.at(1) ->Fill((*mc_Pt)[j]);       
	      h_mcobject_eta.at(1) ->Fill((*mc_Eta)[j]);

	      float dR = dR_mc_mub(it , j, 1);   
	      h_dR_MC.at(1).at(2)  -> Fill(dR) ;  // dR & dPt (Wtop, muon)
	      dR = dPt_Gen_jetlep(it , j, A);
	      h_dPt_MC.at(1)(2) -> Fill(dR) ;
	//  h_mcobject_pt.at(2) ->Fill((*mc_Pt)[it]);       
	h_mcobject_eta.at(2) ->Fill((*mc_Eta)[it]); 

	dR = dR_mc_mub(b_top , j, 1);
	h_dR_MC.at(0)(2)  -> Fill(dR) ;   // dR & dPt (muon, btop)
	dR = dPt_Gen_jetlep(b_top , j, 1);
	h_dPt_MC.at(0)(2) -> Fill(dR) ;   
	h_mcobject_pt.at(0) ->Fill((*mc_Pt)[b_top]);       
	h_mcobject_eta.at(0) ->Fill((*mc_Eta)[b_top]);        
	dR = dR_mc_mub(it , b_top, 1);		 // dR (Whiggs, btop)
	h_dR_MC.at(0)(3)  -> Fill(dR) ; 
	 */
}



void   lvb_qqqq_dCSV::Fill_NPtEta_Histo(int id, int en, int idx)
{
	float ptau1 = 0.0 ;

	//   id for defining particle type, ele =1, mu= 3, bjet = 05, AK8jet= 24(W) , qjet & forward jets = 13
	if (id == 0) {
		h_object_no.at(0) ->Fill(n_Mu.size());
		h_object_no.at(1) ->Fill(b_jet.size());
		h_object_no.at(2) ->Fill(n_AK8Jet.size()) ;
		h_object_no.at(3) ->Fill(n_jet.size()) ;
		h_object_no.at(4) ->Fill(n_forwjet.size()) ;
		h_object_pt.at(0) ->Fill(pf_MET); 
	}

	if ( id == 3){
		h_object_pt.at(idx+1)    ->Fill((*mu_Pt)[en]); 
		h_object_eta.at(idx+1)   ->Fill((*mu_Eta)[en]);
	}

	if( id == 2){

		//h_object_eta.at(id) ->Fill(pf_METEta);
	}
	if( id == 5) {
		h_object_pt.at(idx+3) ->Fill((*jet_Pt)[en]); 
		h_object_eta.at(idx+3) ->Fill((*jet_Eta)[en]);
	}
	if (id == 24 ) {
		h_object_pt.at(idx+6) ->Fill((*AK8_JetPt)[en]); 
		h_object_eta.at(idx+6) ->Fill((*AK8_JetEta)[en]);
		h_AK8_PUPPImass.at(idx) -> Fill((*AK8_puppiSDMass)[en]) ;

		ptau1 = ((*AK8_puppiTau2)[en]/(*AK8_puppiTau1)[en]) ;
		h_AK8_PUPPItau21.at(idx) -> Fill(ptau1);
		ptau1 = ((*AK8_puppiTau3)[en]/(*AK8_puppiTau2)[en]) ;
		h_AK8_PUPPItau32.at(idx) -> Fill(ptau1);
		ptau1 = ((*AK8_Jet_tau4)[en]/(*AK8_Jet_tau2)[en]) ;
		h_AK8_CHStau42.at(idx) -> Fill(ptau1);


	}  
	if(id == 14){
		h_object_pt.at(id+idx) ->Fill((*jet_Pt)[en]);
		h_object_eta.at(id+idx) ->Fill((*jet_Eta)[en]);
	}
}


void   lvb_qqqq_dCSV::Fill_MC_PtEta_Histo(int id, int en)
{
	//   id for defining particle type, 6 = asso bquark, 7= top bquark, 8-9 Wdeacay(Higgs) qi, 10-11 Wdeacay(top) qi, 12 = muon, 13= forw quark
	h_object_pt.at(id) ->Fill((*mc_Pt)[en]); 
	//cout<<", ID = "<< id;
	h_object_eta.at(id) ->Fill((*mc_Eta)[en]);
}


float lvb_qqqq_dCSV::Fill_MET_var(int lep, int lvl)
{
	float pass = 1.0 ;
	//bool pass = true ;
	float dPhi = 0.0 ;
	for(int h =0; h< n_MC; h++){
		if( abs((*mc_PID)[h]) == 14  && abs((*mc_MomPID)[h]) == 24 && abs((*mc_GMomPID)[h]) == 6) {
			dPhi = delta_phi((*mc_Phi)[h] , pf_METPhi);
		}
	}
	if (lvl == 1){
		h_MET_var.at(0) ->Fill(pf_MET);
		float sum_HT = pf_MET + (*mu_Pt)[lep] ;
		h_MET_var.at(1) ->Fill(sum_HT);
		h_MET_var.at(2) ->Fill(dPhi);
	}
	if ( dPhi < 0.6 ) pass = dPhi;
	return pass ;
}

void   lvb_qqqq_dCSV::Fill_Puppi_jet(int jet, int idx)
{// here idx = 0 for toptag, idx= 1,2 for Wtag, idx =3 for fatjet, idx =4,5 for Higgtag

	float ptau3 =   ((*AK8_puppiTau3)[jet]/(*AK8_puppiTau2)[jet]);
	float ptau2 =   ((*AK8_puppiTau2)[jet]/(*AK8_puppiTau1)[jet]) ;
	float ptau4 =   ((*AK8_Jet_tau4)[jet]/(*AK8_Jet_tau2)[jet]) ;
	float mass_puppi = 0.0;
	TLorentzVector v1 ;

	v1.SetPtEtaPhiE((*AK8_JetPt)[jet],(*AK8_JetEta)[jet],(*AK8_JetPhi)[jet],(*AK8_JetEn)[jet]);
	mass_puppi = v1.M();

}


bool lvb_qqqq_dCSV::Cut_Muon(int c_muon){

	bool pass_muon = true;
	if ( (*mu_Pt)[c_muon] <= 40.0 ) pass_muon = false;
	if ( fabs((*mu_Eta)[c_muon]) >= 2.1 ) pass_muon = false;
	bool tight_muon = (*mu_IDbit)[c_muon]>>2 & 1;    //for tight muon ID
	if ( tight_muon == 0 ) pass_muon = false;

	//cout << "Tight mu = " << tight_muon << endl ;
	// for medium ID        
	/*   //      if ((*mu_Chi2NDF)[c_muon] >= 3) pass_muon = false;
	     if ((*mu_chi2LocalPosition)[c_muon] >= 12 ) pass_muon = false;
	     if ((*mu_trkKink)[c_muon] >= 20 ) pass_muon = false;
	     if ((*mu_InnervalidFraction)[c_muon] <= 0.8) pass_muon = false;
	     if ((*mu_segmentCompatibility)[c_muon] <= 0.303 ) pass_muon = false;

	// for tight id
	//     if ((*mu_Chi2NDF)[c_muon] >= 10) pass_muon = false   ;
	if ((*mu_TrkLayers)[c_muon] <= 5  ) pass_muon = false ;
	if ((*mu_PixelHits)[c_muon] <= 0) pass_muon = false ;
	if ((*mu_MuonHits)[c_muon] <= 0 ) pass_muon = false ;
	if ((*mu_Stations)[c_muon] <= 1) pass_muon = false ;
	if((*mu_D0)[c_muon] >= 0.2) pass_muon = false  ;
	if((*mu_Dz)[c_muon] >= 0.5) pass_muon = false   ;  
	 */

	return pass_muon;

}


void lvb_qqqq_dCSV::Cut_bjet(int c_jet, int wp)
{     
	if( (*jet_Pt)[c_jet] >= 30.0){
		if (fabs((*jet_Eta)[c_jet]) < 2.4 ){ 
			if (wp == 1 ){
				if ( ((*jet_DeepCSVTags_b)[c_jet] + (*jet_DeepCSVTags_bb)[c_jet]) >  0.2219) {
					b_jet_loose.push_back(c_jet) ;	
					b_jet.push_back(c_jet) ;	  
				}
				else {
					n_jet.push_back(c_jet) ;
				}}

			if (wp == 2 ){
				if ( ((*jet_DeepCSVTags_b)[c_jet] + (*jet_DeepCSVTags_bb)[c_jet]) >  0.6324) {
					b_jet_medium.push_back(c_jet) ;	
					b_jet.push_back(c_jet) ;	  
				}
				else {
					n_jet.push_back(c_jet) ;
				}}

			if (wp == 3 ){
				if ( ((*jet_DeepCSVTags_b)[c_jet] + (*jet_DeepCSVTags_bb)[c_jet]) >  0.8958) {
					b_jet_tight.push_back(c_jet) ;	
					b_jet.push_back(c_jet) ;	  
				}
				else {
					n_jet.push_back(c_jet) ;
				}}

		}
		else{
			n_forwjet.push_back(c_jet) ;

		}}}


void   lvb_qqqq_dCSV::Cut_AK8jet(int c_jet)
{
	if((*AK8_JetPt)[c_jet] >= 150.0 && fabs((*AK8_JetEta)[c_jet]) <= 2.4) n_AK8Jet.push_back(c_jet) ;
}


//===================dR distribution functions===================

float  lvb_qqqq_dCSV::dR_AK8jet(int Y, int Z, int id8, int lvl)
{ // here lvl value -1 for Puppi & CHS jets dR calculation & lvl = 1 for  dR with tagged jets
	bool pass_dr = false;
	float dPhi ; 
	float dEta;
	float dR   ;        
	if(lvl == 1){
		dPhi = delta_phi((*AK8_JetPhi)[Y], (*AK8_JetPhi)[Z]) ;
		dEta = fabs((*AK8_JetEta)[Y] - (*AK8_JetEta)[Z]);
		dR = Sqrt(dPhi*dPhi + dEta*dEta) ;
	}       
	else {
		dPhi = delta_phi((*AK8_JetPhi)[Y], (*AK8_puppiPhi)[Z]) ;
		dEta = fabs((*AK8_JetEta)[Y] - (*AK8_puppiEta)[Z]);
		dR = Sqrt(dPhi*dPhi + dEta*dEta) ;
	}
	return dR;
}


float  lvb_qqqq_dCSV::dRPlots_bjet(int Y, int Z, int lvl)
{
	// for bjets only  
	float dPhi = delta_phi((*jet_Phi)[Y] , (*jet_Phi)[Z]) ;
	float dEta = (*jet_Eta)[Y] - (*jet_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;     
	return dR;
}


float  lvb_qqqq_dCSV::dR_AK8bjet(int Y, int Z, int id8, int idb, int lvl)
{
	// for topjet id8 = 0, fatjet id8 = 1, Wtag id8 = 2,3
	bool pass_dr = false;
	int d =  id8 +4*idb;

	float dPhi = delta_phi((*AK8_JetPhi)[Y] , (*jet_Phi)[Z]) ;
	float dEta = (*AK8_JetEta)[Y] - (*jet_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;     
	return dR;
}

float  lvb_qqqq_dCSV::dR_mc_mub(int Y, int Z, int A)
{
	// here Z for muon/bquark & Y for top/Wtop & Higgs/Whiggs
	float phi, eta ;
	if ( A == -1 ) {
		phi = (*mc_MomPhi)[Y] ;
		eta = (*mc_MomEta)[Y] ;
	}  
	else {
		phi = (*mc_Phi)[Y] ;
		eta = (*mc_Eta)[Y] ;  
	}
	float dPhi = delta_phi(phi , (*mc_Phi)[Z]) ;
	float dEta = eta - (*mc_Eta)[Z];
	float dR   = Sqrt(dPhi*dPhi + dEta*dEta) ;               
	return dR ; 
}


void  lvb_qqqq_dCSV::dR_mc_qiqj(int Y, int Z,int A)
{
	bool pass_dr = true;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*mc_Phi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*mc_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;       
	pass_dr = false;
}

void  lvb_qqqq_dCSV::dR_mc_bqi(int Y, int Z,int A)
{
	bool pass_dr = true;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*mc_Phi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*mc_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;       
	pass_dr = false;
}


float  lvb_qqqq_dCSV::dR_mcbjet(int Y, int Z, int X, int idx, int lvl)
{
	float pass_dr = 1.0;
	int d = 2*X + Z;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*jet_Phi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*jet_Eta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;
	return dR; 
}


void  lvb_qqqq_dCSV::dR_mcAK8(int Y, int Z, int X, int idx)
{

	bool pass_dr = true;
	int d = 2*X + idx;
	float dPhi = delta_phi((*mc_Phi)[Y] , (*AK8_JetPhi)[Z]) ;
	float dEta = (*mc_Eta)[Y] - (*AK8_JetEta)[Z];
	float dR = Sqrt(dPhi*dPhi + dEta*dEta) ;       
}

float  lvb_qqqq_dCSV::dR_mu(int Y, int Z_jet, int idx, int lvl)
{
	float pass_dr = 2.0;
	double dPhi = delta_phi((*mu_Phi)[Y] , (*jet_Phi)[Z_jet]) ;
	double dEta = (*mu_Eta)[Y] - (*jet_Eta)[Z_jet];
	double dR   = sqrt(dPhi*dPhi + dEta*dEta) ;
	return dR;
}

float lvb_qqqq_dCSV::dR_mu_AK8(int Y, int Z_jet, int idx, int lvl)
{
	bool pass_dr = true;
	double dPhi = delta_phi((*mu_Phi)[Y] , (*AK8_JetPhi)[Z_jet]) ;
	double dEta = (*mu_Eta)[Y] - (*AK8_JetEta)[Z_jet];
	double dR   = sqrt(dPhi*dPhi + dEta*dEta) ;
	return dR;
}


float lvb_qqqq_dCSV::dPt_lep(int jet, int lep, int idx, int lvl)
{
	float  pass_dPt =  0.0;
	float  vec_mag = 0.0 ;
	TVector3 vec_Xprod, vec_lep, vec_jet;
	float dPt = 0.0 ;

	vec_lep.SetPtEtaPhi( (*mu_Pt)[lep], (*mu_Eta)[lep] , (*mu_Phi)[lep] ) ;
	if(idx == 0) vec_jet.SetPtEtaPhi( (*jet_Pt)[jet], (*jet_Eta)[jet], (*jet_Phi)[jet] ) ;
	else {
		vec_jet.SetPtEtaPhi( (*AK8_JetPt)[jet], (*AK8_JetEta)[jet], (*AK8_JetPhi)[jet] ) ; 
	}    

	vec_Xprod = vec_lep. Cross(vec_jet);
	vec_mag = vec_jet.Mag();
	dPt =  (vec_Xprod.Mag() ) / vec_mag ;
	return dPt ;
}


float lvb_qqqq_dCSV::dPt_Gen_jetlep(int jet, int lep, int idx)
{
	// idx = -1 for Mom & idx = 1 for daughter
	float  pass_dPt =  0.0;
	float  vec_mag = 0.0 ;
	TVector3 vec_Xprod, vec_lep, vec_jet;
	float dPt = 0.0 ;

	vec_lep.SetPtEtaPhi( (*mc_Pt)[lep], (*mc_Eta)[lep] , (*mc_Phi)[lep] ) ;
	if(idx == -1) vec_jet.SetPtEtaPhi( (*mc_MomPt)[jet], (*mc_MomEta)[jet] , (*mc_MomPhi)[jet] ) ;
	if(idx == 1 ) vec_jet.SetPtEtaPhi( (*mc_Pt)[jet], (*mc_Eta)[jet] , (*mc_Phi)[jet] ) ;

	vec_Xprod = vec_lep. Cross(vec_jet);
	vec_mag = vec_jet.Mag();
	dPt =  (vec_Xprod.Mag() ) / vec_mag ;                      
	return dPt ;
}


void  lvb_qqqq_dCSV::dR_qjet_objects(int qjet, int obj, int idq, int idob)   
{
	// for reco, mu = 0, bjets = 1&2, ak8jet = 3&4, qjet = 13. for gen, assob =5, topb = 6, q from W(higgs) = 7 & 8, q from W(top) = 9 & 10, mu = 11, forw q = 12,
	bool pass_dr = false;
	float dR = 0.0;
	float dPhi = 0.0 , dEta  = 0.0 ;        
	if( idob == 0){         
		dPhi = delta_phi((*mu_Phi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*mu_Eta)[obj] - (*jet_Eta)[qjet];
	}   

	if( idob >= 5 && idob <=12 ){
		dPhi = delta_phi((*mc_Phi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*mc_Eta)[obj] - (*jet_Eta)[qjet];       
	}

	if( idob ==1 || idob == 2){
		dPhi = delta_phi((*jet_Phi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*jet_Eta)[obj] - (*jet_Eta)[qjet];
	}

	if( idob ==3|| idob == 4){
		dPhi = delta_phi((*AK8_JetPhi)[obj] , (*jet_Phi)[qjet]) ;
		dEta = (*AK8_JetEta)[obj] - (*jet_Eta)[qjet];
	}

	if( idob <= 12 ) {
		dR = Sqrt(dPhi*dPhi + dEta*dEta) ;

	}

	//        return pass_dr;
}


float lvb_qqqq_dCSV::delta_phi(float phi1,float phi2)
{
	const float PI=2.0*acos(0.);
	const float TWOPI=2.0*PI;
	float PHI=fabs( phi1 - phi2 ) ;
	return (PHI<=PI)? PHI : TWOPI-PHI;

}
// ---------------------------------------Reco Object Plots----------------------------------------


void lvb_qqqq_dCSV::RecoPlots_dRHisto()
{
	float dR = 0.0;
	int b1 = -1, b2 =-1 ;
	int Msize  = ( n_Mu.size() >= 2) ? 2 : n_Mu.size() ;
	int Asize  = ( n_AK8Jet.size() >= 3) ? 3 : n_AK8Jet.size() ;
	int Bsize  = ( b_jet.size() >= 3) ? 3 : b_jet.size() ;

	// for muon w.r.t other jets dR plots
	for ( int h = 0 ; h < Msize ; h ++ ){
		b1 = n_Mu[h];   
		for( int g = 0 ; g < Bsize ; g ++ ) {
			b2 = b_jet[g];
			dR = dR_mu(b1, b2, b1, b1);
			h_dR_Reco.at(h).at(g+2) ->Fill(dR) ;
			dR = dPt_lep(b2, b1, 0, b1);
			h_dPt_Reco.at(h).at(g+2) ->Fill(dR) ;
		}
		for( int g = 0 ; g < Asize ; g ++ ) {
			b2 = n_AK8Jet[g];
			dR = dR_mu_AK8(b1, b2, b1, b2);
			h_dR_Reco.at(h).at(g+5) ->Fill(dR) ;
			dR = dPt_lep(b2, b1, 1, b1);
			h_dPt_Reco.at(h).at(g+5) ->Fill(dR) ;
		}  
	}

	// for bjet w.r.t other jets dR plots
	for( int h = 0 ; h < Bsize ; h ++ ) {
		b1 = b_jet[h];
		for( int g = h+1 ; g < Bsize ; g ++ ) {
			b2 = b_jet[g];  
			dR = dRPlots_bjet(b1, b2, -1);
			h_dR_Reco.at(h+2).at(g+2) ->Fill(dR) ;  
		}
		for( int g = 0 ; g < Asize ; g ++ ) {
			b2 = n_AK8Jet[g];
			dR = dR_AK8bjet(b2, b1, b1, 0, 0);
			h_dR_Reco.at(h+2).at(g+5) ->Fill(dR) ; 
		}  
	}

	// for AK8jet dR Plots
	for( int h = 0 ; h < Asize ; h ++ ) {
		b1 = n_AK8Jet[h];
		for( int g = h+1 ; g < Asize ; g ++ ) {
			b2 = n_AK8Jet[g];
			dR = dR_AK8jet(b1, b2, 0, 1);
			h_dR_Reco.at(h+5).at(g+5) ->Fill(dR) ; 
		}  
	}
}



void lvb_qqqq_dCSV::Fill_RecoObject()
{
	int Msize  = ( n_Mu.size() >= 2) ? 2 : n_Mu.size() ;             // id = 3 for muon
	int Asize  = ( n_AK8Jet.size() >= 3) ? 3 : n_AK8Jet.size() ;    // id = 24 for AK8jet
	int Bsize  = ( b_jet.size() >= 3) ? 3 : b_jet.size() ;                // id = 5 for bjet
	int b2 = -1 ;

	Fill_NPtEta_Histo(0, 0, 0) ;                                     // for MET & object population

	for ( int i = 0; i < Msize; i ++)
	{    b2 = n_Mu[i] ;
		Fill_NPtEta_Histo(3, b2, i) ;
	}

	for ( int i = 0; i < Bsize; i ++)
	{    b2 = b_jet[i] ;
		Fill_NPtEta_Histo(5, b2, i) ;
	}

	for ( int i = 0; i < Asize; i ++)
	{    b2 = n_AK8Jet[i] ;
		Fill_NPtEta_Histo(24, b2, i) ;
	}

	RecoPlots_dRHisto() ;  
}

//===============Jet Category Selection=================================
void  lvb_qqqq_dCSV::Higgs_selection( )  // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1 ; 
	int sel = -1 ;
	float Pruned_mass;
	float tau42  ;

	for (int g = 0 ; g < n_AK8Jet.size() ; g ++ ){
		if ( n_AK8Jet[g] == -1 )   continue;
		jt = n_AK8Jet[g] ;
		Pruned_mass  = (*AK8_JetPrunedMass)[jt] ;
		tau42 = ((*AK8_Jet_tau4)[jt]/(*AK8_Jet_tau2)[jt]) ;

		// for higgs jet selection     
		if ((*AK8_JetPt)[jt] < 300.0) continue ;
		if (tau42 > 0.55 ) continue ;
		if (tau42 < 0.02 ) continue ;
		if(Pruned_mass < 110.0 || Pruned_mass > 135.0 ) continue ;
		Higgsjets.push_back(jt);
		n_AK8Jet[g] = -1 ;    
	}
}


void  lvb_qqqq_dCSV::Top_selection( )  // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1 ; 
	int sel = -1 ;
	float tau32  ;
	float SDmass ;

	for (int g = 0 ; g < n_AK8Jet.size() ; g ++ ){
		if ( n_AK8Jet[g] == -1 )   continue;
		sel = -1 ;
		jt = n_AK8Jet[g] ;
		tau32 = ((*AK8_puppiTau3)[jt]/(*AK8_puppiTau2)[jt]) ;
		SDmass = (*AK8_puppiSDMass)[jt];

		// for top jet selection     
		if ((*AK8_JetPt)[jt] < 400.0) continue ;
		if (tau32 > 0.65 ) continue ;
		if (tau32 < 0.02 ) continue ;
		if ( SDmass <140.0 || SDmass > 210.0 ) continue ;

		for(int jg = 0 ; jg < (*n_AK8puppiSDSJ)[jt] ; jg ++ )
		{
			if ((*AK8_puppiSDSJCSV)[jt][jg] > 0.5426) sel = 1;
		}
		if (sel == -1) continue ;
		topjet.push_back(jt) ;
		n_AK8Jet[g] = -1 ;
	}
}


void  lvb_qqqq_dCSV::Wjet_selection()   // tag = 0 for Higgs & tag = 1 for top, W & fat
{
	int jt = -1  ; 
	int sel = -1 ;
	float tau21 ;
	float SDmass ;

	for (int g = 0 ; g < n_AK8Jet.size() ; g ++ ){
		jt = n_AK8Jet[g] ;
		tau21    =  ((*AK8_puppiTau2)[jt]/(*AK8_puppiTau1)[jt]) ;
		SDmass  =  (*AK8_puppiSDMass)[jt];

		// for Wjet selection     
		if ((*AK8_JetPt)[jt] < 220.0) continue ;
		if (tau21 > 0.55 ) continue ;
		if (tau21 < 0.02) continue ;
		if ( SDmass < 65.0 ) continue ;
		if (SDmass > 105.0 ) continue ;
		W_boson.push_back(jt) ;
		n_AK8Jet[g] = -1 ;
	}

}


void  lvb_qqqq_dCSV::Fatjet_selection()
{
	int g3 = -1 ;
	float dR2 = 0.0 ;  
	for (int g = 0 ; g < n_AK8Jet.size() ; g ++ )
	{
		g3 = n_AK8Jet[g]  ;
		if( g3 == -1 ) continue ;
		dR2 =((*AK8_puppiTau2)[g3]/(*AK8_puppiTau1)[g3]) ;
		if (dR2 > 0.80) continue ;
		if (dR2 < 0.02) continue ;     
		if ( (*AK8_puppiSDMass)[g3] > 65.0 && (*AK8_puppiSDMass)[g3] < 210.0 ) continue ;
		//if ( (*AK8_JetPrunedMass)[g3] > 65.0 && (*AK8_puppiSDMass)[g3] < 210.0 ) continue ;
		fat_jet.push_back(g3) ;
	}

}


// ============================Category Study =================
void     lvb_qqqq_dCSV::Wtag0_Category()
{
	if (topjet.size() != 0) {
		if ( fat_jet.size() != 0 ) {
			top_fatjet_Plots() ;	
			eventII ++ ;		                 
		}
	}

	else{
		if (Higgsjets.size() != 0 ) {
			if (b_jet.size() != 0 ) {
				Higgs_lbjet_Plots() ;	
				eventVI ++ ;		                 
			}
		}
	}
}

void     lvb_qqqq_dCSV::Wtag1_Category()
{
	if (topjet.size() == 1 ) {
		top_Wjet_Plots() ; 
		eventI ++ ;		                 		
	}
	else{
		if (b_jet.size() != 0 ) {
			if ( fat_jet.size() != 0 ) W_fatjet_Plots() ;   
		}

	}
}
// ======================Signal Profile Plots ==================
void  lvb_qqqq_dCSV::Wjet_Plots(int bjet ) 
{
	int j = n_Mu[0] ;
	int fill = -1 ;
	int Pt_idx = 0; // for dPt w.r.t muon plots
	int t2 = -1 ;
	int t3 = -1 ;
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;
	float mass_puppi = 0.0;
	TLorentzVector  Ws ;
	int Wsize        =   W_boson.size()  ;
	if ( Wsize != 0 )  h_tag_N.at(0)    ->Fill(Wsize) ;
	Wsize  = ( W_boson.size() >= 4) ? 4 : W_boson.size() ;

	for (int k = 0 ; k < Wsize ; k ++ ) {
		t2   =   W_boson[k] ;

		h_tag_Pt.at(k)   -> Fill((*AK8_JetPt)[t2]) ;
		h_tag_Eta.at(k) -> Fill((*AK8_JetEta)[t2]) ; 
		h_tag_SD.at(k)  ->Fill((*AK8_puppiSDMass)[t2]);

		tau =  ((*AK8_puppiTau2)[t2]/(*AK8_puppiTau1)[t2]) ;
		h_tag_tau.at(k)  -> Fill(tau) ;
		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.at(k) ->Fill(mass_puppi);       

		tau =  dPt_lep( t2, j, Pt_idx, fill);
		tau1 = dR_mu_AK8(j, t2, k, fill) ;
		h_dR_Recomu_tagjet.at(k) ->Fill(tau1);
		h_dPt_lep_tagjet.at(k) ->Fill(tau);


		// for bjet dR plots
		if ( bjet != -1 ) {
			tau = dR_AK8bjet( t2, bjet, k, 0, fill);
			h_dR_Recob1_tagjet.at(k) -> Fill(tau) ;
		}
		else{
			for(int l = 0 ; l < b_jet.size() ; l ++){
				j = b_jet[l] ;
				tau = dR_AK8bjet( t2, j, k, l, fill);
				if( l == 0) h_dR_Recob1_tagjet.at(k) -> Fill(tau) ;
				if( l == 1) h_dR_Recob2_tagjet.at(k) -> Fill(tau) ;
				if( l == 2) h_dR_Recob3_tagjet.at(k) -> Fill(tau) ;
				if( l == 3) h_dR_Recob4_tagjet.at(k) -> Fill(tau) ;
			}
		}
	}

	// for dR plots in Wjets
	for ( int h = 0 ; h < W_boson.size() ; h ++ ){
		t2 = W_boson[h];
		for ( int g = h+1; g < W_boson.size(); g++) {
			t3 = W_boson[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
			h_dR_tagjet.at(h).at(g) ->Fill(tau) ;
		}  
	}
}

void  lvb_qqqq_dCSV::Topjet_Plots( ) 
{
	int j = n_Mu[0] ;
	int fill = -1 ;
	int Pt_idx = 2; //for dPt w.r.t muon plots
	int t2 = -1 ;
	int t3 = -1 ;     
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;     
	float mass_puppi = 0.0;
	TLorentzVector  Ws ;
	int Tsize        =   topjet.size()  ;
	if ( Tsize != 0 )  h_tag_N.at(1)     ->Fill(Tsize) ;
	Tsize  = ( topjet.size() >= 2) ? 2 : topjet.size() ;

	for (int k = 0 ; k < Tsize ; k ++ ) {

		t2   =   topjet[k] ;
		h_tag_Pt.at(k+4)  -> Fill((*AK8_JetPt)[t2]) ;
		h_tag_Eta.at(k+4) -> Fill((*AK8_JetEta)[t2]) ; 
		h_tag_SD.at(k+4)  ->Fill((*AK8_puppiSDMass)[t2]);

		tau =  ((*AK8_puppiTau2)[t2]/(*AK8_puppiTau1)[t2]) ;
		h_tag_tau.at(k+4)  -> Fill(tau) ;
		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.at(k+4) ->Fill(mass_puppi);

		if(Pt_idx == 2) {
			tau =  dPt_lep( t2, j, Pt_idx, fill);
			tau1 = dR_mu_AK8(j, t2, k, fill) ;
			h_dPt_lep_tagjet.at(k+4) ->Fill(tau);
			h_dR_Recomu_tagjet.at(k+4) ->Fill(tau1);
		}     

		// for bjet dR plots
		for(int l = 0 ; l < b_jet.size() ; l ++){
			j = b_jet[l] ;
			tau = dR_AK8bjet( t2, j, k, l, fill);
			if( l == 0) h_dR_Recob1_tagjet.at(k+4) -> Fill(tau) ;
			if( l == 1) h_dR_Recob2_tagjet.at(k+4) -> Fill(tau) ;
			if( l == 2) h_dR_Recob3_tagjet.at(k+4) -> Fill(tau) ;
			if( l == 3) h_dR_Recob4_tagjet.at(k+4) -> Fill(tau) ;
		}    
	}  

	// for dR plots in Wjets
	for ( int h = 0 ; h < Tsize ; h ++ ){
		t2 = topjet[h];
		for ( int g = h+1; g < Tsize; g++) {
			t3 = topjet[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
			h_dR_tagjet.at(h+4).at(g+4) ->Fill(tau) ;
		}  
	}
}


void  lvb_qqqq_dCSV::Higgsjet_Plots( int bjet) 
{
	int j = n_Mu[0] ;
	int fill = -1 ;
	int Pt_idx = 3; //for dPt w.r.t muon plots
	int t2 = -1 ;
	int t3 = -1 ;     
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;
	float mass_puppi = 0.0;
	TLorentzVector  Ws ;
	int Hsize        =   Higgsjets.size()  ;
	if ( Hsize != 0 )  h_tag_N.at(2)     ->Fill(Hsize) ;
	Hsize  = ( Higgsjets.size() >= 2) ? 2 : Higgsjets.size() ;

	for (int k = 0 ; k < Hsize ; k ++ ) {

		t2   =   Higgsjets[k] ;
		h_tag_Pt.at(k+6)  -> Fill((*AK8_JetPt)[t2]) ;
		h_tag_Eta.at(k+6) -> Fill((*AK8_JetEta)[t2]) ; 
		h_tag_SD.at(k+6)  ->Fill((*AK8_JetPrunedMass)[t2]);

		tau =  ((*AK8_Jet_tau4)[t2]/(*AK8_Jet_tau2)[t2]) ;
		h_tag_tau.at(k+6)  -> Fill(tau) ;

		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.at(k+6) ->Fill(mass_puppi);

		if(Pt_idx == 3) {
			tau =  dPt_lep( t2, j, Pt_idx, fill);
			tau1 = dR_mu_AK8(j, t2, k, fill) ;
			h_dPt_lep_tagjet.at(k+6) ->Fill(tau);
			h_dR_Recomu_tagjet.at(k+6) ->Fill(tau1);

		}   

		if ( bjet != -1 ) {
			// for bjet dR plots
			for(int l = 0 ; l < b_jet.size() ; l ++){
				j = b_jet[l] ;
				tau = dR_AK8bjet( t2, j, k, l, fill);
				if( l == 0) h_dR_Recob1_tagjet.at(k+6) -> Fill(tau) ;
				if( l == 1) h_dR_Recob2_tagjet.at(k+6) -> Fill(tau) ;
				if( l == 2) h_dR_Recob3_tagjet.at(k+6) -> Fill(tau) ;
				if( l == 3) h_dR_Recob4_tagjet.at(k+6) -> Fill(tau) ;
			}
		} 
	}

	// for dR plots in Wjets
	for ( int h = 0 ; h < Hsize ; h ++ ){
		t2 = Higgsjets[h];
		for ( int g = h+1; g < Hsize; g++) {
			t3 = Higgsjets[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
			h_dR_tagjet.at(h+6).at(g+6) ->Fill(tau) ;
		}  
	}
}


void  lvb_qqqq_dCSV::Fatjet_Plots(int topbjet ) 
{
	int j = n_Mu[0] ;
	int fill = -1 ;
	int Pt_idx = 4; //for dPt w.r.t muon plots     
	int t2 = -1 ;
	int t3 = -1 ;     
	float  tau = 0.0 ;
	float  tau1 = 0.0 ;
	float mass_puppi = 0.0;
	TLorentzVector  Ws ;
	int Fsize        =   fat_jet.size()  ;
	if ( Fsize != 0 )  h_tag_N.at(3)     ->Fill(Fsize) ;
	Fsize  = ( fat_jet.size() >= 4) ? 4 : fat_jet.size() ;
	for (int k = 0 ; k < Fsize ; k ++ ) {

		t2   =   fat_jet[k] ;
		h_tag_Pt.at(k+8)  -> Fill((*AK8_JetPt)[t2]) ;
		h_tag_Eta.at(k+8) -> Fill((*AK8_JetEta)[t2]) ; 
		h_tag_SD.at(k+8)  ->Fill((*AK8_puppiSDMass)[t2]);

		tau =  ((*AK8_puppiTau2)[t2]/(*AK8_puppiTau1)[t2]) ;
		h_tag_tau.at(k+8)  -> Fill(tau) ;
		Ws.SetPtEtaPhiE((*AK8_JetPt)[t2],(*AK8_JetEta)[t2],(*AK8_JetPhi)[t2],(*AK8_JetEn)[t2]);
		mass_puppi = Ws.M();
		h_tag_Mass.at(k+8) ->Fill(mass_puppi);

		if(Pt_idx < 6) {
			tau =  dPt_lep( t2, j, Pt_idx, fill);
			tau1 = dR_mu_AK8(j, t2, k, fill) ;
			h_dPt_lep_tagjet.at(k+8) ->Fill(tau);
			h_dR_Recomu_tagjet.at(k+8) ->Fill(tau1);
		}     

		// for bjet dR plots
		if ( topbjet != -1 ) {
			tau = dR_AK8bjet( t2, topbjet, k, 0, fill);
			h_dR_Recob1_tagjet.at(k+8) -> Fill(tau) ;
		}
		else {
			for(int l = 0 ; l < b_jet.size() ; l ++){
				j = b_jet[l] ;
				tau = dR_AK8bjet( t2, j, k, l, fill);
				if( l == 0) h_dR_Recob1_tagjet.at(k+8) -> Fill(tau) ;
				if( l == 1) h_dR_Recob2_tagjet.at(k+8) -> Fill(tau) ;
				if( l == 2) h_dR_Recob3_tagjet.at(k+8) -> Fill(tau) ;
				if( l == 3) h_dR_Recob4_tagjet.at(k+8) -> Fill(tau) ;
			}
		} 
	}

	// for dR plots in fatjet
	for ( int h = 0 ; h < Fsize ; h ++ ){
		t2 = fat_jet[h];
		for ( int g = h+1; g < Fsize; g++) {
			t3 = fat_jet[g];
			tau = dR_AK8jet(t3, t2, t2, 1);
			h_dR_tagjet.at(h+8).at(g+8) ->Fill(tau) ;
		}  
	}
}


void  lvb_qqqq_dCSV::Category_Wjet_Plot( int cat, int W, int idx) 
{    
	// II is for whether, W is subleading Wboson belonging to category IV
	int fill = 0 ;
	float  tau = 0.0 ;
	TLorentzVector  Ws ;
	int Wsize        =   W_boson.size()  ;     

	h_Histo_Pt.at(cat).at(idx)  -> Fill((*AK8_JetPt)[W]) ;
	h_Histo_Eta.at(cat).at(idx) -> Fill((*AK8_JetEta)[W]) ; 
	h_Histo_SD.at(cat).at(idx)  ->Fill((*AK8_puppiSDMass)[W]);

	tau =  ((*AK8_puppiTau2)[W]/(*AK8_puppiTau1)[W]) ;
	h_Histo_tau.at(cat).at(idx)  -> Fill(tau) ;
	Ws.SetPtEtaPhiE((*AK8_JetPt)[W],(*AK8_JetEta)[W],(*AK8_JetPhi)[W],(*AK8_JetEn)[W]);
	float   Et_W        =  Ws.Et() ;
	float   pT_W  =  (*AK8_JetPt)[W] ;
	double mass_puppi  = Sqrt (fabs ((Et_W * Et_W) - (pT_W * pT_W) ) );    
	h_Histo_Mass.at(cat).at(idx) ->Fill(mass_puppi);
}


void  lvb_qqqq_dCSV::Category_Top_Plot( int cat,int top, int idx) 
{    
	int fill = -1 ;
	float  tau = 0.0 ;     
	TLorentzVector  Ws ;
	int Tsize        =   topjet.size()  ;

	h_Histo_Pt.at(cat).at(idx)  -> Fill((*AK8_JetPt)[top]) ;
	h_Histo_Eta.at(cat).at(idx) -> Fill((*AK8_JetEta)[top]) ; 
	h_Histo_SD.at(cat).at(idx)  ->Fill((*AK8_puppiSDMass)[top]);

	tau =  ((*AK8_puppiTau3)[top]/(*AK8_puppiTau2)[top]) ;
	h_Histo_tau.at(cat).at(idx)  -> Fill(tau) ;
	Ws.SetPtEtaPhiE((*AK8_JetPt)[top],(*AK8_JetEta)[top],(*AK8_JetPhi)[top],(*AK8_JetEn)[top]);

	float   Et_top        =  Ws.Et() ;
	float   pT_top  =  (*AK8_JetPt)[top] ;
	double mass_puppi  = Sqrt (fabs ((Et_top * Et_top) - (pT_top * pT_top) ) );    
	h_Histo_Mass.at(cat).at(idx) ->Fill(mass_puppi);

}


void  lvb_qqqq_dCSV::Category_Higgs_Plot( int cat, int higg, int idx) 
{      
	float  tau = 0.0 ;
	TLorentzVector  Ws ;
	int Hsize        =   Higgsjets.size()  ;

	h_Histo_Pt.at(cat).at(idx)  -> Fill((*AK8_JetPt)[higg]) ;
	h_Histo_Eta.at(cat).at(idx) -> Fill((*AK8_JetEta)[higg]) ; 
	h_Histo_SD.at(cat).at(idx)  ->Fill((*AK8_JetPrunedMass)[higg]);

	tau =  ((*AK8_Jet_tau4)[higg]/(*AK8_Jet_tau2)[higg]) ;
	h_Histo_tau.at(cat).at(idx)  -> Fill(tau) ;

	Ws.SetPtEtaPhiE((*AK8_JetPt)[higg],(*AK8_JetEta)[higg],(*AK8_JetPhi)[higg],(*AK8_JetEn)[higg]);
	float   Et_top        =  Ws.Et() ;
	float   pT_top  =  (*AK8_JetPt)[higg] ;
	double mass_puppi  = Sqrt (fabs ((Et_top * Et_top) - (pT_top * pT_top) ) );    
	h_Histo_Mass.at(cat).at(idx) ->Fill(mass_puppi);      
}        


void  lvb_qqqq_dCSV::Category_Fat_Plot( int cat,int fat, int idx) 
{    
	int fill = -1 ;
	float  tau = 0.0 ;
	TLorentzVector  Ws ;
	int Tsize        =   fat_jet.size()  ;

	h_Histo_Pt.at(cat).at(idx)  -> Fill((*AK8_JetPt)[fat]) ;
	h_Histo_Eta.at(cat).at(idx) -> Fill((*AK8_JetEta)[fat]) ; 
	h_Histo_SD.at(cat).at(idx)  ->Fill((*AK8_puppiSDMass)[fat]);

	tau =  ((*AK8_puppiTau2)[fat]/(*AK8_puppiTau1)[fat]) ;
	h_Histo_tau.at(cat).at(idx)  -> Fill(tau) ;
	Ws.SetPtEtaPhiE((*AK8_JetPt)[fat],(*AK8_JetEta)[fat],(*AK8_JetPhi)[fat],(*AK8_JetEn)[fat]);
	float   Et_top        =  Ws.Et() ;
	float   pT_top  =  (*AK8_JetPt)[fat] ;
	double mass_puppi  = Sqrt (fabs ((Et_top * Et_top) - (pT_top * pT_top) ) );    
	h_Histo_Mass.at(cat).at(idx) ->Fill(mass_puppi);

}


void lvb_qqqq_dCSV::TagJets_dRPlots()
{
	float dR = 0.0;
	int b1 = -1, b2 =-1 ;
	int Fsize  = ( fat_jet.size() >= 4) ? 4 : fat_jet.size() ;
	int Wsize  = ( W_boson.size() >= 4) ? 4 : W_boson.size() ;
	int Tsize  = ( topjet.size() >= 2) ? 2 : topjet.size() ;
	int Hsize  = ( Higgsjets.size() >= 2 ) ? 2 : Higgsjets.size() ;
	// for Wjet w.r.t tag jets dR plots
	for ( int h = 0 ; h < Wsize ; h ++ ){
		b1 = W_boson[h];  
		for( int g = 0 ; g < Tsize ; g ++ ) {
			b2 = topjet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h).at(g+4) ->Fill(dR) ;
		}
		for( int g = 0 ; g < Hsize ; g ++ ) {
			b2 = Higgsjets[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h).at(g+6) ->Fill(dR) ;
		}
		for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = fat_jet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h).at(g+8) ->Fill(dR) ;
		}  
	} 
	// for topjet w.r.t tag jets dR plots
	for( int h = 0 ; h < Tsize ; h ++ ) {
		b1 = topjet[h]; 
		for( int g = 0 ; g < Hsize ; g ++ ) {
			b2 = Higgsjets[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h+4).at(g+6) ->Fill(dR) ;
		}
		for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = fat_jet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h+4).at(g+8) ->Fill(dR) ;
		}  
	}
	// for Higgsjet w.r.t tag jets dR plots  
	for( int h = 0 ; h < Hsize ; h ++ ) {
		b1 = Higgsjets[h];  
		for( int g = 0 ; g < Fsize ; g ++ ) {
			b2 = fat_jet[g];
			dR = dR_AK8jet(b2, b1, b1, 1);
			h_dR_tagjet.at(h+6).at(g+8) ->Fill(dR) ;
		}    
	}

}

// ============Mt Calculation ========================

void  lvb_qqqq_dCSV::toptag_MTCalculation( ) 
{
	int top = topjet[0] ;
	int W;
	if ( W_boson.size() == 0 )    W = fat_jet[0] ;
	if ( W_boson.size() == 1 )    W  = W_boson[0] ;

	int mu = n_Mu[0] ;
	TLorentzVector v_top, v_mu, v_W ;

	v_top.SetPtEtaPhiE((*AK8_JetPt)[top], (*AK8_JetEta)[top], (*AK8_JetPhi)[top], (*AK8_JetEn)[top] ) ;
	v_W.SetPtEtaPhiE((*AK8_JetPt)[W], (*AK8_JetEta)[W], (*AK8_JetPhi)[W], (*AK8_JetEn)[W] ) ;
	v_mu.SetPtEtaPhiE( (*mu_Pt)[mu], (*mu_Eta)[mu], (*mu_Phi)[mu], (*mu_En)[mu]);

	float Et_top        =  v_top.Et() ;
	float Et_W         = v_W.Et() ;
	float Et_mu        = v_mu.Et() ;
	double Et_sum    = 0.0 ;

	float pT_top  =  (*AK8_JetPt)[top] ;
	float pT_W  =  (*AK8_JetPt)[W] ;
	float pT_mu   =  (*mu_Pt)[mu] ;
	double pT_sum  = 0.0 ;

	double  W_Trans                = 0.0 ;
	double  Higgs_Trans            = 0.0 ;
	double  TprimEt_Trans          = 0.0 ;
	double  trans_mass            = 0.0 ;

	// for cos(deltaPhi) calculation  
	float  dPhi           =  delta_phi((*mu_Phi)[mu] , pf_METPhi) ;   // for muon & MET dphi
	float  MET_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[W] , pf_METPhi)  ;   //  for W & MET dphi
	float  MET2_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[top] , pf_METPhi)  ;   //  for top & MET dphi
	float  MET3_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[W] , (*mu_Phi)[mu])  ;   //  for W & muon dphi
	float  mu_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[top] , (*mu_Phi)[mu])  ;   //  for top & muon dphi
	float  mu2_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[W] , (*AK8_JetPhi)[top])  ;   //  for W & top dphi
	float  W_phi       =   pT_W * Cos(dPhi) ;

	float  Dot_prod = 0.0   ;
	//  for W transverse mass    
	Et_sum = Et_sum + Et_mu  + pf_MET  ;
	pT_sum  = pT_sum + ( pf_MET * pf_MET ) + (pT_mu * pT_mu) ;  
	Dot_prod  =  Dot_prod +  2.0 * pT_mu * MET_phi ;
	W_Trans = fabs (( Et_sum * Et_sum) - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	h_object_MT.at(0) -> Fill( trans_mass) ;

	// for Higgs Transverse mass
	Et_sum               =     Et_sum +   Et_W ;
	pT_sum              =     pT_sum + ( pT_W* pT_W ) ;
	Dot_prod            =     Dot_prod + 2.0 * pT_W *( MET2_phi + mu_phi ) ;
	Higgs_Trans        =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass          =    Sqrt ( Higgs_Trans) ;
	h_object_MT.at(1) -> Fill( trans_mass) ;

	// for Tprime Transverse mass
	Et_sum          = Et_sum + Et_top ;
	pT_sum        =    pT_sum + (pT_top * pT_top) ;    
	Dot_prod       =   Dot_prod + 2.0 * pT_top * ( MET3_phi + mu2_phi + W_phi ) ;
	TprimEt_Trans =   fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( TprimEt_Trans ) ;
	h_object_MT.at(2) -> Fill(trans_mass) ;

	//  for top tag transverse mass              
	double  MT_vec  = fabs ((Et_top * Et_top) - (pT_top * pT_top) ) ;
	trans_mass  = Sqrt ( MT_vec) ;
	h_object_MT.at(3) -> Fill(trans_mass) ;
	MT_vec = 0.0 ;

	// for Wboson transverse mass 
	MT_vec  = fabs ((Et_W * Et_W) - (pT_W * pT_W) ) ;   
	trans_mass  = Sqrt ( MT_vec) ;
	h_object_MT.at(4) -> Fill(trans_mass) ;

}


// =========for Higgs tag Mt calculation =======================


void  lvb_qqqq_dCSV::Higgstag_MTCalculation( int topbjet) 
{
	int higgs = Higgsjets[0] ;  
	int mu = n_Mu[0] ;
	TLorentzVector v_higgs, v_mu, v_topbjet ;

	v_higgs.SetPtEtaPhiE((*AK8_JetPt)[higgs], (*AK8_JetEta)[higgs], (*AK8_JetPhi)[higgs], (*AK8_JetEn)[higgs] ) ;
	v_topbjet.SetPtEtaPhiE((*jet_Pt)[topbjet], (*jet_Eta)[topbjet], (*jet_Phi)[topbjet], (*jet_En)[topbjet] ) ;
	v_mu.SetPtEtaPhiE( (*mu_Pt)[mu], (*mu_Eta)[mu], (*mu_Phi)[mu], (*mu_En)[mu]);

	float Et_higgs        =  v_higgs.Et() ;
	float Et_bjet          = v_topbjet.Et() ;
	float Et_mu          = v_mu.Et() ;
	double Et_sum       = 0.0 ;

	float pT_higgs       =  (*AK8_JetPt)[higgs] ;
	float pT_bjet         =  (*jet_Pt)[topbjet] ;
	float pT_mu         =  (*mu_Pt)[mu] ;
	double pT_sum      =  0.0 ;

	double  W_Trans                = 0.0 ;
	double  top_Trans                = 0.0 ;  
	double  Higgs_Trans            = 0.0 ;
	double  TprimEt_Trans        = 0.0 ;
	double  trans_mass             = 0.0 ;

	// for cos(deltaPhi) calculation  
	float  dPhi           =  delta_phi((*mu_Phi)[mu] , pf_METPhi) ;   // for muon & MET dphi
	float  MET_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet] , pf_METPhi)  ;   //  for topbjet & MET dphi
	float  MET2_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs] , pf_METPhi)  ;   //  for higgs & MET dphi
	float  MET3_phi    =  pf_MET * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*mu_Phi)[mu])  ;   //  for topbjet & muon dphi
	float  mu_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*AK8_JetPhi)[higgs], (*mu_Phi)[mu]) ; //  for higgs & muon dphi
	float  mu2_phi       =   pT_mu * Cos(dPhi) ;

	dPhi                   =  delta_phi((*jet_Phi)[topbjet], (*AK8_JetPhi)[higgs])  ; //  for W & higgs dphi
	float  bjet_phi       =   pT_bjet * Cos(dPhi) ;

	float  Dot_prod = 0.0   ;
	//  for W transverse mass    
	Et_sum = Et_sum + Et_mu  + pf_MET  ;
	pT_sum  = pT_sum + ( pf_MET * pf_MET ) + (pT_mu * pT_mu) ;  
	Dot_prod  =  Dot_prod +  2.0 * pT_mu * MET_phi ;
	W_Trans = fabs (( Et_sum * Et_sum) - ( pT_sum + Dot_prod) ) ;
	trans_mass  = Sqrt ( W_Trans) ;
	h_object_MT.at(0) -> Fill( trans_mass) ;

	//  for Top transverse mass 
	Et_sum               =     Et_sum +   Et_bjet ;
	pT_sum              =     pT_sum + ( pT_bjet* pT_bjet ) ;
	Dot_prod            =     Dot_prod + 2.0 * pT_bjet *( MET2_phi + mu_phi ) ;
	top_Trans           =     fabs (( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) )  ;
	trans_mass          =    Sqrt ( top_Trans) ;
	h_object_MT.at(1) -> Fill( trans_mass) ;

	// for Tprime Transverse mass
	Et_sum              = Et_sum + Et_higgs ;
	pT_sum             =    pT_sum + (pT_higgs * pT_higgs) ;    
	Dot_prod            =   Dot_prod + 2.0 * pT_higgs * ( MET3_phi + mu2_phi + bjet_phi ) ;
	TprimEt_Trans     =   fabs(( Et_sum * Et_sum)  - ( pT_sum + Dot_prod) ) ;
	trans_mass         = Sqrt ( TprimEt_Trans ) ;
	h_object_MT.at(2) -> Fill(trans_mass) ;

	//  for higgs tag transverse mass              
	double  MT_vec  = fabs ((Et_higgs * Et_higgs) - (pT_higgs * pT_higgs) ) ;
	trans_mass  = Sqrt ( MT_vec) ;
	h_object_MT.at(3) -> Fill(trans_mass) ;
	MT_vec = 0.0 ;

	/* for Wboson transverse mass 
	   MT_vec  = fabs ((Et_W * Et_W) - (pT_W * pT_W) ) ;   
	   trans_mass  = Sqrt ( MT_vec) ;
	   h_object_MT.at(4) -> Fill(trans_mass) ;
	 */                
}


void  lvb_qqqq_dCSV::genlvl_MTCalculation(int type ) 
{
	TLorentzVector mc_Tprime, mc_mu ;
	int ID = -1 ;
	float Et = 0.0 ;
	float Pt2 = 0.0 ;
	float MT = 0.0 ;

	float  dPhi      =  delta_phi((*mc_Phi)[T_higgs] , (*mc_Phi)[T_top]) ;   // for muon & MET dphi

	mc_Tprime.SetPtEtaPhiM((*mc_MomPt)[T_higgs],(*mc_MomEta)[T_higgs],(*mc_MomPhi)[T_higgs], (*mc_MomMass)[T_higgs] ) ;

	Et  =  (*mc_Et)[T_higgs] + (*mc_Et)[T_top] ;
	Pt2  = Power( (*mc_Pt)[T_higgs], 2)  + Power ((*mc_Pt)[T_top], 2 ) + 2.0 * (*mc_Pt)[T_higgs] * (*mc_Pt)[T_top] * Cos(dPhi) ;
	MT =  Sqrt (fabs( Power( Et , 2 ) -  Pt2  ) ) ;
	h_object_MT.at(8) -> Fill (MT );

	MT = Sqrt (fabs( Power( (*mc_Et)[T_higgs] , 2 ) - Power ( (*mc_Pt)[T_higgs] , 2 ) ) );
	h_object_MT.at(6) -> Fill ( MT );

	MT = Sqrt( fabs( Power( (*mc_Et)[T_top] , 2 ) - Power ( (*mc_Pt)[T_top] , 2 ) ) );
	h_object_MT.at(7) -> Fill ( MT );	

	if ( type  == -1 ) {
		if ( Higgs_Mu.size() != 0 ) ID = Higgs_Mu[0] ;
	}
	if ( type  == 1 ) {
		if ( top_Mu.size() != 0 )  ID = top_Mu[0] ;
	}
	if ( ID != -1 ) {  
		mc_mu.SetPtEtaPhiM((*mc_MomPt)[ID],(*mc_MomEta)[ID],(*mc_MomPhi)[ID], (*mc_MomMass)[ID] ) ;
		Et  =  mc_mu.Et();
		Pt2  = (*mc_MomPt)[ID] ;
		MT =  Sqrt (fabs( Power( Et , 2 ) - Power( Pt2 , 2 ) ) ) ;
		h_object_MT.at(5) -> Fill (MT );  
	}   
}

//===========================Category Plots ===============

void     lvb_qqqq_dCSV::top_fatjet_Plots()
{
	TString Cat_pass  = "no" ;
	int jet1 = topjet[0] ;
	int jet2 = fat_jet[0] ;
	int j2 = n_Mu[0] ;
	int lvl = -1 ;
	int mu = -1 ;
	int top = -1 ;
	int W = -1;                 
	float dR  =   dR_mu_AK8(j2, jet2, j2, lvl) ;
	float dPt =    dPt_lep( jet2, j2, 1, lvl);

	if ( dR_mu_AK8(j2, jet1, j2, lvl) >  2.0)  mu = 1 ;
	if ( dR_AK8jet(jet2, jet1, jet1, 1) > 2.0) top = 1 ;
	// if ( ( (dR > 0.2 && dR < 2.0 ) || dPt < 20.0 ) ) W = 1 ;     
	if ( dR < 2.5 ) W = 1 ;
	if ( mu != -1 && top != -1 && W != -1 ) {
		Cat_pass  = "yes" ;
		CatII_Objects.push_back(jet1);               // topjet at 0th place
		CatII_Objects.push_back(jet2);			  //  fatjet at 1st place	
		CatII_Objects.push_back(j2);				  //  muon at 2nd place & MET is global
		CatII_Objects_Plots();
	}				

	//	if ( Cat_pass  == "yes" ) { }

	// genlvl_MTCalculation(1) ;    			// for muon from Higgs decayed W type
}			


void  lvb_qqqq_dCSV::top_Wjet_Plots() {
	TString Cat_pass  = "no" ;
	int jet1 = topjet[0] ;
	int jet2 = W_boson[0] ;
	int j = n_Mu[0] ;
	int lvl = -1 ;
	int mu = -1 ;
	int top = -1 ;
	int W = -1;                 
	float dR  =   dR_mu_AK8(j, jet2, j, lvl) ;
	float dPt =    dPt_lep( jet2, j, 1, lvl);

	if ( dR_mu_AK8(j, jet1, j, lvl) >  2.0)  mu = 1 ;
	if ( dR_AK8jet(jet2, jet1, jet1, 1) > 2.0) top = 1 ;
	//if ( ( (dR > 0.2 && dR < 2.0 ) || dPt < 20.0 ) ) W = 1 ;     
	if ( dR < 2.5 ) W =1 ;
	if ( mu != -1 && top != -1 & W != -1  )  {
		Cat_pass  = "yes" ;
		CatI_Objects.push_back(jet1);               // topjet at 0th place
		CatI_Objects.push_back(jet2);			  //  Wjet at 1st place
		CatI_Objects.push_back(j);				  //  muon at 2nd place & MET is global
		CatI_Objects_Plots();   
	}
	if ( Cat_pass  == "yes" ) {
		/*   
		     toptag_MTCalculation( ) ;			
		     genlvl_MTCalculation( 1) ;      // for muon from Higgs decayed W type
		 */
	}
}


void  lvb_qqqq_dCSV::Higgs_lbjet_Plots() {
	TString Cat_pass  = "no" ;
	int jet1 = Higgsjets[0] ;
	int jet2 = -1 ; 
	int b2  = -1 ;
	int j = n_Mu[0] ;
	float dR = 0.0, dPt = 0.0 ;
	int lvl = -1 ;
	bool pass = true ;
	int mu = -1 , top = -1, higg = -1 , bjet = -1 ;
	if ( dR_mu_AK8(j, jet1, j, lvl) < 2.0 )  pass  = false ;
	int size4  = ( b_jet.size() >= 4) ? 4 : b_jet.size() ;
	for( int h =0 ; h < size4 ; h ++ ){
		if ( jet2 != -1 ) continue ;
		b2  = b_jet[h] ;                
		dR = dR_mu(j, b2, 0, lvl);
		dPt = dPt_lep(j, b2, 0, lvl);   
		if ( dR < 2.5 ) top = 1  ;    
		if(((dR > 0.4 && dR < 2.0)  || dPt > 20.0))  mu = 1 ;                                            
		if ( (*jet_Pt)[b2] > 40.0 ) bjet = 1 ;
		if ( dR_AK8bjet( jet1, b2, 6, 0, lvl) > 2.0 )  higg = 1 ;           
		if ( top == 1&& mu == 1 && higg == 1 && bjet == 1)  jet2  = b2 ;
		mu = -1;
		top = -1;
		higg = -1;
	}    
	if ( pass == true && jet2 != -1 ) {
		Cat_pass  = "yes"  ; 
		CatVI_Objects.push_back(-1);    //vague value as all histogram defined from index =1                 
		CatVI_Objects.push_back(jet1);               // Higgsjet at 0th place
		CatVI_Objects.push_back(jet2);			  //  bjet at 1st place	
		CatVI_Objects.push_back(j);				  //  muon at 2nd place & MET is global
		CatVI_Objects_Plots();    
	}
	if ( Cat_pass  == "yes" ) {   
		/*     
		       Higgstag_MTCalculation( jet2) ;
		       genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
		 */
	}

} 



void  lvb_qqqq_dCSV::W_fatjet_Plots() {
	TString Cat_pass  = "no" ;
	int jet1 = W_boson[0] ;
	int jet2 = -1 ;  // jet2 for Wfatjet
	int jet3 = -1 ;  // jet3 for lvfatjet
	int j    = n_Mu[0] ;
	float dR, dPt ;
	int lvl = -1 ;
	int mu = -1 , top = -1 , higg = -1 , topbjet = -1 , hW = -1, hfat = -1, tb = -1, tmu = -1 ;                   
	int b2  = -1;     

	for ( int g = 0; g < fat_jet.size() ; g ++ ) {
		b2  = fat_jet[g] ;
		dPt = dR_AK8jet(b2, jet1, jet1, 1) ;
		dR =  dR_mu_AK8(j, b2, j, lvl) ;
		if (  dR < 2.0 ) jet2 =  b2 ;
		if (  dR > 2.0 ) jet3 =  b2 ;
	}		


	int size4  = ( b_jet.size() >= 4) ? 4 : b_jet.size() ;
	for( int h =0 ; h < size4 ; h ++ ){
		b2  = b_jet[h] ;         
		//         if ( dR_AK8bjet( jet2, b2, h, 1, lvl) > 2.0 ) higg = 1 ; 
		//         if ( dR_AK8bjet( jet3, b2, h, 1, lvl) > 2.0 ) hfat = 1 ;          // for Cat V
		if ( dR_AK8bjet( jet1, b2, h, 0, lvl) < 2.0 ) top = 1 ; 
		if ( dR_AK8bjet( jet1, b2, h, 0, lvl) > 2.0 ) hW = 1 ;                   
		if ( dR_mu(j, b2, h, lvl) > 2.5 ) mu = 1 ;  
		if ( dR_mu(j, b2, h, lvl) < 2.5 ) tmu = 1 ;           
		if ( tmu == 1 && hW == 1 ) tb = b2 ;               
		if ( mu == 1 && top == 1  ) topbjet = b2 ;    
		mu = -1;  tmu = -1;          
		higg = -1 ;  top = - 1 ;    
		hW = -1;  hfat = -1 ;
	}

	dR  =   dR_mu_AK8(j, jet2, j, lvl) ;
	if ( dR_mu_AK8(j, jet1, j, lvl) > 2.0) mu = 1;
	//    if ( dR_AK8jet(jet2, jet1, jet1, 1) > 2.5 )higg = 1;
	if ( topbjet != -1 && mu == 1 && jet2 != -1 )  {
		Cat_pass  = "yes"  ;
		CatIII_Objects.push_back(jet1) ;
		CatIII_Objects.push_back(jet2) ;
		CatIII_Objects.push_back(topbjet) ;
		CatIII_Objects.push_back(j) ;        
	}               
	//   if ( dR_AK8jet(jet3, jet1, jet1, 1) < 2.5 ) hW = 1;
	if ( tb != -1 && mu == 1 && jet3 != -1 )  {
		CatV_Objects.push_back(jet1) ;
		CatV_Objects.push_back(jet3) ;
		CatV_Objects.push_back(tb) ;
		CatV_Objects.push_back(j) ;        
	}

	if (  CatV_Objects.size() != 0 )  	CatV_Objects_Plots() ;       
	if( CatIII_Objects.size() != 0 && CatV_Objects.size() == 0)     CatIII_Objects_Plots() ;

	if ( Cat_pass  == "no" ) {          		
		//		       genlvl_MTCalculation(1) ;    // for muon from Higgs decayed W type	
	}    
}


void  lvb_qqqq_dCSV::WW_lvbjet_Plots() 
{
	TString Cat_pass  = "no" ;
	int jet1 = W_boson[0] ;   // jet1 = top decayed jet
	int jet2 = -1 ;   // jet2 = higgs decayed jet
	int j = n_Mu[0] ;
	float dR, dPt ;
	int lvl = -1 ;
	int b2 = -1 ;   
	int topbjet  = -1 ;    
	int mu = -1, top = -1 , higg = -1, bj = -1 ;
	//  for differentiating top W from higgs W  
	dR  =   dR_mu_AK8(j, jet1, j, lvl) ;
	b2   =   W_boson[1];  
	dPt  =  dR_mu_AK8(j, b2, j, lvl)  ;       
	if ( dPt > dR) {
		jet2 = W_boson[0];
		jet1 = W_boson[1] ;
	}             
	else {
		jet2 = W_boson[1] ; 
	}                     
	W_boson.clear();
	W_boson.push_back(jet1) ;
	W_boson.push_back(jet2) ; 
	jet1 = W_boson[0];
	jet2 = W_boson[1] ;                                     

	//  for selecting top bjet                   
	int size4  = ( b_jet.size() >= 4) ? 4 : b_jet.size() ;
	for( int h =0 ; h < size4 ; h ++ ){
		if ( topbjet != -1) continue ;
		b2  = b_jet[h] ;               
		if ( dR_mu(j, b2, h, lvl) > 2.0 ) mu = 1 ;
		if ( dR_AK8bjet( jet2, b2, h, 1, lvl) > 2.0 ) higg = 1 ; 
		if ( dR_AK8bjet( jet1, b2, h, 0, lvl) < 2.0 ) top = 1 ;  
		if ( (*jet_Pt)[b2] > 40.0 ) bj = 1;           
		if ( mu == 1 && higg == 1 && top == 1 && bj == 1) topbjet = b2 ;
		mu = -1; 
		higg = -1;
		top = -1 ; 
		bj = -1 ;
	}

	if ( dR_mu_AK8(j, jet1, j, lvl) > 2.0 ) top = 1 ;
	dR  =   dR_mu_AK8(j, jet2, j, lvl) ;    
	if ( dR < 2.0 ) higg = 1 ;
	if ( topbjet != -1 && higg != -1 && top != -1 ) {
		Cat_pass  = "yes" ;
		CatIV_Objects.push_back(jet1);               // Wjet at 0th place
		CatIV_Objects.push_back(jet2);			  //  Wjet at 1st place
		CatIV_Objects.push_back(topbjet);				//	bjet at 2nd place	
		CatIV_Objects.push_back(j);				  //  muon at 3rd place & MET is global
		CatIV_Objects_Plots();              
	}
	if ( Cat_pass  == "yes" ) {
		/*        
			  genlvl_MTCalculation(1) ;    // for muon from Higgs decayed W type 
		 */
	}
}
//=======================Plotting Functions=====================

void     lvb_qqqq_dCSV::CatI_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(0).at(3)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     

	for(int g = 0 ; g < CatI_Objects.size(); g ++ )
	{   
		prtcl = CatI_Objects[g] ;
		if (g == 0)       	  Category_Top_Plot( 0,prtcl, g) ; 
		if (g == 1)          Category_Wjet_Plot( 0,prtcl, g) ; 
		if (g == 2)  {
			h_Histo_Pt.at(0).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(0).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}        
		for( int h = g+1; h < CatI_Objects.size() ; h++)
		{ 
			j ++ ;
			b2 = CatI_Objects[h] ;          
			if ( j == 0)   {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(0).at(j) ->Fill(dR) ;
			}
			if ( j >= 1 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(0).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(0).at(j) ->Fill(dR) ;         
			}             				  
		}
	}
	//         genlvl_MTCalculation(1) ;
	//     	toptag_MTCalculation() ;
}     	


void     lvb_qqqq_dCSV::CatII_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(1).at(3)  -> Fill(pf_MET) ;    // for MET being  global & [cat)(histono.]				     

	for(int g = 0 ; g < CatII_Objects.size(); g ++ )
	{   
		prtcl = CatII_Objects[g] ;
		if (g == 0)       	  Category_Top_Plot( 1,prtcl, g) ; 
		if (g == 1)          Category_Fat_Plot( 1,prtcl, g) ; 
		if (g == 2)  {
			h_Histo_Pt.at(1).at(g) -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(1).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatII_Objects.size() ; h++)
		{ 
			j ++ ;
			b2 = CatII_Objects[h] ;          
			if ( j == 0)   {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(1).at(j) ->Fill(dR) ;
			}
			if ( j >= 1 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(1).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(1).at(j) ->Fill(dR) ;         
			}             

		}
	}
	//     	toptag_MTCalculation() ;
}     	


void     lvb_qqqq_dCSV::CatIII_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(2).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     

	for(int g = 0 ; g < CatIII_Objects.size(); g ++ )
	{   
		prtcl = CatIII_Objects[g] ;
		if (g == 0)       	  Category_Wjet_Plot( 2,prtcl, g) ; 
		if (g == 1)       	  Category_Fat_Plot( 2,prtcl, g) ;         
		if (g == 2)    {
			h_Histo_Pt.at(2).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(2).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(2).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(2).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatIII_Objects.size() ; h++)
		{ 
			j  ++ ;
			b2 = CatIII_Objects[h] ;   
			if ( j == 0 ) {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;                             
			}       
			if ( j == 3 || j == 1 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;
			}
			if ( j == 4 || j == 2 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(2).at(j) ->Fill(dR) ;         
			}
			if ( j == 5 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(2).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(2).at(j) ->Fill(dR) ;   
			}             
		}
	}

	//          Higgstag_MTCalculation( jet2) ;
	//          genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}  


void     lvb_qqqq_dCSV::CatIV_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(3).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     

	for(int g = 0 ; g < CatIV_Objects.size(); g ++ )
	{   
		prtcl = CatIV_Objects[g] ;
		if (g == 0)       	  Category_Wjet_Plot( 3,prtcl, g) ; 
		if (g == 1)       	  Category_Wjet_Plot( 3,prtcl, g) ;         
		if (g == 2)    {
			h_Histo_Pt.at(3).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(3).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(3).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(3).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatIV_Objects.size() ; h++)
		{ 
			j  ++ ;
			b2 = CatIV_Objects[h] ;   
			if ( j == 0 ) {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;                             
			}       
			if ( j == 3 || j == 1 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;
			}
			if ( j == 4 || j == 2 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(3).at(j) ->Fill(dR) ;         
			}
			if ( j == 5 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(3).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(3).at(j) ->Fill(dR) ;   
			}             
		}
	}

	//          Higgstag_MTCalculation( jet2) ;
	//          genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}     	


void     lvb_qqqq_dCSV::CatV_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(4).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     

	for(int g = 0 ; g < CatV_Objects.size(); g ++ )
	{   
		prtcl = CatV_Objects[g] ;
		if (g == 0)       	  Category_Wjet_Plot( 4,prtcl, g) ; 
		if (g == 1)       	  Category_Fat_Plot( 4,prtcl, g) ;         
		if (g == 2)    {
			h_Histo_Pt.at(4).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(4).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(4).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(4).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatV_Objects.size() ; h++)
		{ 
			j  ++ ;
			b2 = CatV_Objects[h] ;   
			if ( j == 0 ) {
				dR = dR_AK8jet(b2, prtcl, prtcl, 1);
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;                             
			}       
			if ( j == 3 || j == 1 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;
			}
			if ( j == 4 || j == 2 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(4).at(j) ->Fill(dR) ;         
			}
			if ( j == 5 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(4).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(4).at(j) ->Fill(dR) ;   
			}             
		}
	}

	//          Higgstag_MTCalculation( jet2) ;
	//          genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}  


void     lvb_qqqq_dCSV::CatVI_Objects_Plots()
{
	int prtcl = -1, b2 = -1, j = -1 ;     
	float dR = 0.0 ;                    
	h_Histo_Pt.at(5).at(4)  -> Fill(pf_MET) ;    // for MET being  global & [cat).at(histono.]				     

	for(int g = 1 ; g < CatVI_Objects.size(); g ++ )
	{   
		prtcl = CatVI_Objects[g] ;
		if (g == 1)       	  Category_Higgs_Plot( 5,prtcl, g) ; 
		if (g == 2)    {
			h_Histo_Pt.at(5).at(g)  -> Fill((*jet_Pt)[prtcl]) ;
			h_Histo_Eta.at(5).at(g) -> Fill((*jet_Eta)[prtcl]) ; 
		}
		if (g == 3)  {
			h_Histo_Pt.at(5).at(g)  -> Fill((*mu_Pt)[prtcl]) ;
			h_Histo_Eta.at(5).at(g) -> Fill((*mu_Eta)[prtcl]) ; 
		}

		for( int h = g+1; h < CatVI_Objects.size() ; h++)
		{ 
			j ++ ;
			b2 = CatVI_Objects[h] ;          
			if ( j == 0 )   {               
				dR = dR_AK8bjet( prtcl, b2 , g, h, j) ;
				h_Histo_dR.at(5).at(j) ->Fill(dR) ;
			}
			if ( j == 1 )   {
				dR = dR_mu_AK8(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(5).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, h) ;
				h_Histo_dPt.at(5).at(j) ->Fill(dR) ;         
			}
			if ( j == 2 ) {
				dR = dR_mu(b2, prtcl, prtcl , h) ;
				h_Histo_dR.at(5).at(j) ->Fill(dR) ;
				dR =  dPt_lep(prtcl, b2, j, 0) ;
				h_Histo_dPt.at(5).at(j) ->Fill(dR) ;   
			}             
		}
	}

	//          Higgstag_MTCalculation( jet2) ;
	//          genlvl_MTCalculation(-1) ;  // for muon from top decayed W type
}     	

//=========================Vague Functions need to be dump=================



//=======================Default Funtions=========================
//void lvb_qqqq_dCSV::Init(TTree *tree)
void lvb_qqqq_dCSV::Init(TChain *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
	mc_PID = 0;
	mc_Vtx = 0;
	mc_Vty = 0;
	mc_Vtz = 0;
	mc_Pt = 0;
	mc_Mass = 0;
	mc_Eta = 0;
	mc_Phi = 0;
	mc_E = 0;
	mc_Et = 0;
	mc_GMomPID = 0;
	mc_MomPID = 0;
	mc_MomPt = 0;
	mc_MomMass = 0;
	mc_MomEta = 0;
	mc_MomPhi = 0;
	mc_StatusFlag = 0;
	mc_Parentage = 0;
	mc_Status = 0;
	mc_CalIsoDR03 = 0;
	mc_TrkIsoDR03 = 0;
	mc_CalIsoDR04 = 0;
	mc_TrkIsoDR04 = 0;
	ele_Charge = 0;
	ele_En = 0;
	ele_D0 = 0;
	ele_Dz = 0;
	ele_Pt = 0;
	ele_Eta = 0;
	ele_Phi = 0;
	ele_R9 = 0;
	ele_SCEta = 0;
	ele_SCPhi = 0;
	ele_HoverE = 0;
	ele_EoverP = 0;
	ele_EoverPout = 0;
	ele_EoverPInv = 0;
	ele_dEtaAtVtx = 0;
	ele_dPhiAtVtx = 0;
	ele_SigmaIEtaIEtaFull5x5 = 0;
	ele_SigmaIPhiIPhiFull5x5 = 0;
	ele_ConvVeto = 0;
	ele_MissHits = 0;
	ele_PFChIso = 0;
	ele_PFPhoIso = 0;
	ele_PFNeuIso = 0;
	ele_PFMiniIso = 0;
	ele_dEtaseedAtVtx = 0;
	mu_Pt = 0;
	mu_En = 0;
	mu_Eta = 0;
	mu_Phi = 0;
	mu_Charge = 0;
	mu_IDbit = 0;
	mu_D0 = 0;
	mu_Dz = 0;
	mu_Chi2NDF = 0;
	mu_InnerD0 = 0;
	mu_InnerDz = 0;
	mu_InnervalidFraction = 0;
	mu_segmentCompatibility = 0;
	mu_chi2LocalPosition = 0;
	mu_trkKink = 0;
	mu_PFChIso = 0;
	mu_PFPhoIso = 0;
	mu_PFNeuIso = 0;
	mu_PFMiniIso = 0;
	mu_TrkLayers = 0;
	mu_PixelLayers = 0;
	mu_PixelHits = 0;
	mu_MuonHits = 0;
	mu_Stations = 0;
	mu_Matches = 0;
	mu_TrkQuality = 0;
	mu_IsoTrk = 0;


	jet_Pt = 0;
	jet_En = 0;
	jet_Eta = 0;
	jet_Phi = 0;
	jet_Area = 0;
	jet_Mt = 0;
	jet_CSV2BJetTags = 0;
	jet_JetProbabilityBJetTags = 0;
	jet_pfCombinedMVAV2BJetTags = 0;
	jet_DeepCSVTags_b = 0;
	jet_DeepCSVTags_bb = 0;
	jet_DeepCSVTags_c = 0;
	jet_DeepCSVTags_cc = 0;
	jet_DeepCSVTags_udsg = 0;
	jet_PartonID = 0;
	jet_HadFlvr = 0;
	jet_GenJetEn = 0;
	jet_GenJetPt = 0;
	jet_GenJetEta = 0;
	jet_GenJetPhi = 0;
	jet_GenPartonID = 0;
	jet_GenEn = 0;
	jet_GenPt = 0;
	jet_GenEta = 0;
	jet_GenPhi = 0;
	jet_GenPartonMomID = 0;
	jet_PFLooseId = 0;
	jet_ID = 0;
	jet_PUID = 0;
	jet_PUFullID = 0;
	jet_CHF = 0;
	jet_NHF = 0;
	jet_CEF = 0;
	jet_NEF = 0;
	jet_NCH = 0;
	jet_NNP = 0;
	jet_MUF = 0;
	AK8_JetPt = 0;
	AK8_JetEn = 0;
	AK8_JetEta = 0;
	AK8_JetPhi = 0;
	AK8_JetMass = 0;
	AK8_Jet_tau1 = 0;
	AK8_Jet_tau2 = 0;
	AK8_Jet_tau3 = 0;
	AK8_Jet_tau4 = 0;
	AK8_JetCHF = 0;
	AK8_JetNHF = 0;
	AK8_JetCEF = 0;
	AK8_JetNEF = 0;
	AK8_JetNCH = 0;
	AK8_JetNNP = 0;
	AK8_JetMUF = 0;
	AK8_Jetnconstituents = 0;
	AK8_JetPFLooseId = 0;
	AK8_JetPFTightLepVetoId = 0;
	AK8_JetSoftDropMass = 0;
	AK8_JetPrunedMass = 0;
	AK8_JetpfBoostedDSVBTag = 0;
	AK8_JetCSV = 0;
	AK8_JetDSVnewV4 = 0;
	AK8_puppiPt = 0;
	AK8_puppiMass = 0;
	AK8_puppiEta = 0;
	AK8_puppiPhi = 0;
	AK8_puppiTau1 = 0;
	AK8_puppiTau2 = 0;
	AK8_puppiTau3 = 0;
	AK8_puppiSDMass = 0;
	AK8_JetPartonID = 0;
	AK8_JetHadFlvr = 0;
	AK8_JetGenJetIndex = 0;
	AK8_JetGenJetEn = 0;
	AK8_JetGenJetPt = 0;
	AK8_JetGenJetEta = 0;
	AK8_JetGenJetPhi = 0;
	AK8_JetGenPartonID = 0;
	AK8_JetGenEn = 0;
	AK8_JetGenPt = 0;
	AK8_JetGenEta = 0;
	AK8_JetGenPhi = 0;
	AK8_JetGenPartonMomID = 0;
	n_AK8SDSJ = 0;
	AK8_SDSJPt = 0;
	AK8_SDSJEta = 0;
	AK8_SDSJPhi = 0;
	AK8_SDSJMass = 0;
	AK8_SDSJE = 0;
	AK8_SDSJCharge = 0;
	AK8_SDSJFlavour = 0;
	AK8_SDSJCSV = 0;
	n_AK8puppiSDSJ = 0;
	AK8_puppiSDSJPt = 0;
	AK8_puppiSDSJEta = 0;
	AK8_puppiSDSJPhi = 0;
	AK8_puppiSDSJMass = 0;
	AK8_puppiSDSJE = 0;
	AK8_puppiSDSJCharge = 0;
	AK8_puppiSDSJFlavour = 0;
	AK8_puppiSDSJCSV = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("v_event", &v_event, &b_v_event);
	fChain->SetBranchAddress("n_Vtx", &n_Vtx, &b_n_Vtx);
	fChain->SetBranchAddress("n_GoodVtx", &n_GoodVtx, &b_n_GoodVtx);
	fChain->SetBranchAddress("n_TrksPV", &n_TrksPV, &b_n_TrksPV);
	fChain->SetBranchAddress("is_PVGood", &is_PVGood, &b_is_PVGood);
	fChain->SetBranchAddress("v_vtx", &v_vtx, &b_v_vtx);
	fChain->SetBranchAddress("v_vty", &v_vty, &b_v_vty);
	fChain->SetBranchAddress("v_vtz", &v_vtz, &b_v_vtz);
	fChain->SetBranchAddress("rho_Central", &rho_Central, &b_rho_Central);
	fChain->SetBranchAddress("n_MC", &n_MC, &b_n_MC);
	fChain->SetBranchAddress("mc_PID", &mc_PID, &b_mc_PID);
	fChain->SetBranchAddress("mc_Vtx", &mc_Vtx, &b_mc_Vtx);
	fChain->SetBranchAddress("mc_Vty", &mc_Vty, &b_mc_Vty);
	fChain->SetBranchAddress("mc_Vtz", &mc_Vtz, &b_mc_Vtz);
	fChain->SetBranchAddress("mc_Pt", &mc_Pt, &b_mc_Pt);
	fChain->SetBranchAddress("mc_Mass", &mc_Mass, &b_mc_Mass);
	fChain->SetBranchAddress("mc_Eta", &mc_Eta, &b_mc_Eta);
	fChain->SetBranchAddress("mc_Phi", &mc_Phi, &b_mc_Phi);
	fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
	fChain->SetBranchAddress("mc_Et", &mc_Et, &b_mc_Et);
	fChain->SetBranchAddress("mc_GMomPID", &mc_GMomPID, &b_mc_GMomPID);
	fChain->SetBranchAddress("mc_MomPID", &mc_MomPID, &b_mc_MomPID);
	fChain->SetBranchAddress("mc_MomPt", &mc_MomPt, &b_mc_MomPt);
	fChain->SetBranchAddress("mc_MomMass", &mc_MomMass, &b_mc_MomMass);
	fChain->SetBranchAddress("mc_MomEta", &mc_MomEta, &b_mc_MomEta);
	fChain->SetBranchAddress("mc_MomPhi", &mc_MomPhi, &b_mc_MomPhi);
	fChain->SetBranchAddress("mc_StatusFlag", &mc_StatusFlag, &b_mc_StatusFlag);
	fChain->SetBranchAddress("mc_Parentage", &mc_Parentage, &b_mc_Parentage);
	fChain->SetBranchAddress("mc_Status", &mc_Status, &b_mc_Status);
	fChain->SetBranchAddress("mc_CalIsoDR03", &mc_CalIsoDR03, &b_mc_CalIsoDR03);
	fChain->SetBranchAddress("mc_TrkIsoDR03", &mc_TrkIsoDR03, &b_mc_TrkIsoDR03);
	fChain->SetBranchAddress("mc_CalIsoDR04", &mc_CalIsoDR04, &b_mc_CalIsoDR04);
	fChain->SetBranchAddress("mc_TrkIsoDR04", &mc_TrkIsoDR04, &b_mc_TrkIsoDR04);
	fChain->SetBranchAddress("gen_MET", &gen_MET, &b_gen_MET);
	fChain->SetBranchAddress("gen_METPhi", &gen_METPhi, &b_gen_METPhi);
	fChain->SetBranchAddress("pf_MET", &pf_MET, &b_pf_MET);
	fChain->SetBranchAddress("pf_METPhi", &pf_METPhi, &b_pf_METPhi);
	fChain->SetBranchAddress("pf_METsumEt", &pf_METsumEt, &b_pf_METsumEt);
	fChain->SetBranchAddress("pf_METmEtSig", &pf_METmEtSig, &b_pf_METmEtSig);
	fChain->SetBranchAddress("pf_METSig", &pf_METSig, &b_pf_METSig);
	fChain->SetBranchAddress("n_Ele", &n_Ele, &b_n_Ele);
	fChain->SetBranchAddress("ele_Charge", &ele_Charge, &b_ele_Charge);
	fChain->SetBranchAddress("ele_En", &ele_En, &b_ele_En);
	fChain->SetBranchAddress("ele_D0", &ele_D0, &b_ele_D0);
	fChain->SetBranchAddress("ele_Dz", &ele_Dz, &b_ele_Dz);
	fChain->SetBranchAddress("ele_Pt", &ele_Pt, &b_ele_Pt);
	fChain->SetBranchAddress("ele_Eta", &ele_Eta, &b_ele_Eta);
	fChain->SetBranchAddress("ele_Phi", &ele_Phi, &b_ele_Phi);
	fChain->SetBranchAddress("ele_R9", &ele_R9, &b_ele_R9);
	fChain->SetBranchAddress("ele_SCEta", &ele_SCEta, &b_ele_SCEta);
	fChain->SetBranchAddress("ele_SCPhi", &ele_SCPhi, &b_ele_SCPhi);
	fChain->SetBranchAddress("ele_HoverE", &ele_HoverE, &b_ele_HoverE);
	fChain->SetBranchAddress("ele_EoverP", &ele_EoverP, &b_ele_EoverP);
	fChain->SetBranchAddress("ele_EoverPout", &ele_EoverPout, &b_ele_EoverPout);
	fChain->SetBranchAddress("ele_EoverPInv", &ele_EoverPInv, &b_ele_EoverPInv);
	fChain->SetBranchAddress("ele_dEtaAtVtx", &ele_dEtaAtVtx, &b_ele_dEtaAtVtx);
	fChain->SetBranchAddress("ele_dPhiAtVtx", &ele_dPhiAtVtx, &b_ele_dPhiAtVtx);
	fChain->SetBranchAddress("ele_SigmaIEtaIEtaFull5x5", &ele_SigmaIEtaIEtaFull5x5, &b_ele_SigmaIEtaIEtaFull5x5);
	fChain->SetBranchAddress("ele_SigmaIPhiIPhiFull5x5", &ele_SigmaIPhiIPhiFull5x5, &b_ele_SigmaIPhiIPhiFull5x5);
	fChain->SetBranchAddress("ele_ConvVeto", &ele_ConvVeto, &b_ele_ConvVeto);
	fChain->SetBranchAddress("ele_MissHits", &ele_MissHits, &b_ele_MissHits);
	fChain->SetBranchAddress("ele_PFChIso", &ele_PFChIso, &b_ele_PFChIso);
	fChain->SetBranchAddress("ele_PFPhoIso", &ele_PFPhoIso, &b_ele_PFPhoIso);
	fChain->SetBranchAddress("ele_PFNeuIso", &ele_PFNeuIso, &b_ele_PFNeuIso);
	fChain->SetBranchAddress("ele_PFMiniIso", &ele_PFMiniIso, &b_ele_PFMiniIso);
	fChain->SetBranchAddress("ele_dEtaseedAtVtx", &ele_dEtaseedAtVtx, &b_ele_dEtaseedAtVtx);
	fChain->SetBranchAddress("N_Mu", &N_Mu, &b_N_Mu);
	fChain->SetBranchAddress("mu_Pt", &mu_Pt, &b_mu_Pt);
	fChain->SetBranchAddress("mu_En", &mu_En, &b_mu_En);
	fChain->SetBranchAddress("mu_Eta", &mu_Eta, &b_mu_Eta);
	fChain->SetBranchAddress("mu_Phi", &mu_Phi, &b_mu_Phi);
	fChain->SetBranchAddress("mu_Charge", &mu_Charge, &b_mu_Charge);
	fChain->SetBranchAddress("mu_IDbit", &mu_IDbit, &b_mu_IDbit);
	fChain->SetBranchAddress("mu_D0", &mu_D0, &b_mu_D0);
	fChain->SetBranchAddress("mu_Dz", &mu_Dz, &b_mu_Dz);
	fChain->SetBranchAddress("mu_Chi2NDF", &mu_Chi2NDF, &b_mu_Chi2NDF);
	fChain->SetBranchAddress("mu_InnerD0", &mu_InnerD0, &b_mu_InnerD0);
	fChain->SetBranchAddress("mu_InnerDz", &mu_InnerDz, &b_mu_InnerDz);
	fChain->SetBranchAddress("mu_InnervalidFraction", &mu_InnervalidFraction, &b_mu_InnervalidFraction);
	fChain->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility, &b_mu_segmentCompatibility);
	fChain->SetBranchAddress("mu_chi2LocalPosition", &mu_chi2LocalPosition, &b_mu_chi2LocalPosition);
	fChain->SetBranchAddress("mu_trkKink", &mu_trkKink, &b_mu_trkKink);
	fChain->SetBranchAddress("mu_PFChIso", &mu_PFChIso, &b_mu_PFChIso);
	fChain->SetBranchAddress("mu_PFPhoIso", &mu_PFPhoIso, &b_mu_PFPhoIso);
	fChain->SetBranchAddress("mu_PFNeuIso", &mu_PFNeuIso, &b_mu_PFNeuIso);
	fChain->SetBranchAddress("mu_PFMiniIso", &mu_PFMiniIso, &b_mu_PFMiniIso);
	fChain->SetBranchAddress("mu_TrkLayers", &mu_TrkLayers , &b_mu_TrkLayers);
	fChain->SetBranchAddress("mu_PixelLayers", &mu_PixelLayers , &b_mu_PixelLayers);
	fChain->SetBranchAddress("mu_PixelHits", &mu_PixelHits , &b_mu_PixelHits);
	fChain->SetBranchAddress("mu_MuonHits", &mu_MuonHits , &b_mu_MuonHits);
	fChain->SetBranchAddress("mu_Stations", &mu_Stations , &b_mu_Stations);
	fChain->SetBranchAddress("mu_Matches", &mu_Matches , &b_mu_Matches);
	fChain->SetBranchAddress("mu_TrkQuality", &mu_TrkQuality , &b_mu_TrkQuality);
	fChain->SetBranchAddress("mu_IsoTrk", &mu_IsoTrk , &b_mu_IsoTrk);

	fChain->SetBranchAddress("n_Jet", &n_Jet, &b_n_Jet);
	fChain->SetBranchAddress("jet_Pt", &jet_Pt, &b_jet_Pt);
	fChain->SetBranchAddress("jet_En", &jet_En, &b_jet_En);
	fChain->SetBranchAddress("jet_Eta", &jet_Eta, &b_jet_Eta);
	fChain->SetBranchAddress("jet_Phi", &jet_Phi, &b_jet_Phi);
	fChain->SetBranchAddress("jet_Area", &jet_Area, &b_jet_Area);
	fChain->SetBranchAddress("jet_Mt", &jet_Mt, &b_jet_Mt);
	fChain->SetBranchAddress("jet_CSV2BJetTags", &jet_CSV2BJetTags, &b_jet_CSV2BJetTags);
	fChain->SetBranchAddress("jet_JetProbabilityBJetTags", &jet_JetProbabilityBJetTags, &b_jet_JetProbabilityBJetTags);
	fChain->SetBranchAddress("jet_pfCombinedMVAV2BJetTags", &jet_pfCombinedMVAV2BJetTags, &b_jet_pfCombinedMVAV2BJetTags);
	fChain->SetBranchAddress("jet_DeepCSVTags_b", &jet_DeepCSVTags_b, &b_jet_DeepCSVTags_b);
	fChain->SetBranchAddress("jet_DeepCSVTags_bb", &jet_DeepCSVTags_bb, &b_jet_DeepCSVTags_bb);
	fChain->SetBranchAddress("jet_DeepCSVTags_c", &jet_DeepCSVTags_c, &b_jet_DeepCSVTags_c);
	fChain->SetBranchAddress("jet_DeepCSVTags_cc", &jet_DeepCSVTags_cc, &b_jet_DeepCSVTags_cc);
	fChain->SetBranchAddress("jet_DeepCSVTags_udsg", &jet_DeepCSVTags_udsg, &b_jet_DeepCSVTags_udsg);
	fChain->SetBranchAddress("jet_PartonID", &jet_PartonID, &b_jet_PartonID);
	fChain->SetBranchAddress("jet_HadFlvr", &jet_HadFlvr, &b_jet_HadFlvr);
	fChain->SetBranchAddress("jet_GenJetEn", &jet_GenJetEn, &b_jet_GenJetEn);
	fChain->SetBranchAddress("jet_GenJetPt", &jet_GenJetPt, &b_jet_GenJetPt);
	fChain->SetBranchAddress("jet_GenJetEta", &jet_GenJetEta, &b_jet_GenJetEta);
	fChain->SetBranchAddress("jet_GenJetPhi", &jet_GenJetPhi, &b_jet_GenJetPhi);
	fChain->SetBranchAddress("jet_GenPartonID", &jet_GenPartonID, &b_jet_GenPartonID);
	fChain->SetBranchAddress("jet_GenEn", &jet_GenEn, &b_jet_GenEn);
	fChain->SetBranchAddress("jet_GenPt", &jet_GenPt, &b_jet_GenPt);
	fChain->SetBranchAddress("jet_GenEta", &jet_GenEta, &b_jet_GenEta);
	fChain->SetBranchAddress("jet_GenPhi", &jet_GenPhi, &b_jet_GenPhi);
	fChain->SetBranchAddress("jet_GenPartonMomID", &jet_GenPartonMomID, &b_jet_GenPartonMomID);
	fChain->SetBranchAddress("jet_PFLooseId", &jet_PFLooseId, &b_jet_PFLooseId);
	fChain->SetBranchAddress("jet_ID", &jet_ID, &b_jet_ID);
	fChain->SetBranchAddress("jet_PUID", &jet_PUID, &b_jet_PUID);
	fChain->SetBranchAddress("jet_PUFullID", &jet_PUFullID, &b_jet_PUFullID);
	fChain->SetBranchAddress("jet_CHF", &jet_CHF, &b_jet_CHF);
	fChain->SetBranchAddress("jet_NHF", &jet_NHF, &b_jet_NHF);
	fChain->SetBranchAddress("jet_CEF", &jet_CEF, &b_jet_CEF);
	fChain->SetBranchAddress("jet_NEF", &jet_NEF, &b_jet_NEF);
	fChain->SetBranchAddress("jet_NCH", &jet_NCH, &b_jet_NCH);
	fChain->SetBranchAddress("jet_NNP", &jet_NNP, &b_jet_NNP);
	fChain->SetBranchAddress("jet_MUF", &jet_MUF, &b_jet_MUF);
	fChain->SetBranchAddress("N_AK8Jet", &N_AK8Jet, &b_N_AK8Jet);
	fChain->SetBranchAddress("AK8_JetPt", &AK8_JetPt, &b_AK8_JetPt);
	fChain->SetBranchAddress("AK8_JetEn", &AK8_JetEn, &b_AK8_JetEn);
	fChain->SetBranchAddress("AK8_JetEta", &AK8_JetEta, &b_AK8_JetEta);
	fChain->SetBranchAddress("AK8_JetPhi", &AK8_JetPhi, &b_AK8_JetPhi);
	fChain->SetBranchAddress("AK8_JetMass", &AK8_JetMass, &b_AK8_JetMass);
	fChain->SetBranchAddress("AK8_Jet_tau1", &AK8_Jet_tau1, &b_AK8_Jet_tau1);
	fChain->SetBranchAddress("AK8_Jet_tau2", &AK8_Jet_tau2, &b_AK8_Jet_tau2);
	fChain->SetBranchAddress("AK8_Jet_tau3", &AK8_Jet_tau3, &b_AK8_Jet_tau3);
	fChain->SetBranchAddress("AK8_Jet_tau4", &AK8_Jet_tau4, &b_AK8_Jet_tau4);
	fChain->SetBranchAddress("AK8_JetCHF", &AK8_JetCHF, &b_AK8_JetCHF);
	fChain->SetBranchAddress("AK8_JetNHF", &AK8_JetNHF, &b_AK8_JetNHF);
	fChain->SetBranchAddress("AK8_JetCEF", &AK8_JetCEF, &b_AK8_JetCEF);
	fChain->SetBranchAddress("AK8_JetNEF", &AK8_JetNEF, &b_AK8_JetNEF);
	fChain->SetBranchAddress("AK8_JetNCH", &AK8_JetNCH, &b_AK8_JetNCH);
	fChain->SetBranchAddress("AK8_JetNNP", &AK8_JetNNP, &b_AK8_JetNNP);
	fChain->SetBranchAddress("AK8_JetMUF", &AK8_JetMUF, &b_AK8_JetMUF);
	fChain->SetBranchAddress("AK8_Jetnconstituents", &AK8_Jetnconstituents, &b_AK8_Jetnconstituents);
	fChain->SetBranchAddress("AK8_JetPFLooseId", &AK8_JetPFLooseId, &b_AK8_JetPFLooseId);
	fChain->SetBranchAddress("AK8_JetPFTightLepVetoId", &AK8_JetPFTightLepVetoId, &b_AK8_JetPFTightLepVetoId);
	fChain->SetBranchAddress("AK8_JetSoftDropMass", &AK8_JetSoftDropMass, &b_AK8_JetSoftDropMass);
	fChain->SetBranchAddress("AK8_JetPrunedMass", &AK8_JetPrunedMass, &b_AK8_JetPrunedMass);
	fChain->SetBranchAddress("AK8_JetpfBoostedDSVBTag", &AK8_JetpfBoostedDSVBTag, &b_AK8_JetpfBoostedDSVBTag);
	fChain->SetBranchAddress("AK8_JetCSV", &AK8_JetCSV, &b_AK8_JetCSV);
	fChain->SetBranchAddress("AK8_JetDSVnewV4", &AK8_JetDSVnewV4, &b_AK8_JetDSVnewV4);
	fChain->SetBranchAddress("AK8_puppiPt", &AK8_puppiPt, &b_AK8_puppiPt);
	fChain->SetBranchAddress("AK8_puppiMass", &AK8_puppiMass, &b_AK8_puppiMass);
	fChain->SetBranchAddress("AK8_puppiEta", &AK8_puppiEta, &b_AK8_puppiEta);
	fChain->SetBranchAddress("AK8_puppiPhi", &AK8_puppiPhi, &b_AK8_puppiPhi);
	fChain->SetBranchAddress("AK8_puppiTau1", &AK8_puppiTau1, &b_AK8_puppiTau1);
	fChain->SetBranchAddress("AK8_puppiTau2", &AK8_puppiTau2, &b_AK8_puppiTau2);
	fChain->SetBranchAddress("AK8_puppiTau3", &AK8_puppiTau3, &b_AK8_puppiTau3);
	fChain->SetBranchAddress("AK8_puppiSDMass", &AK8_puppiSDMass, &b_AK8_puppiSDMass);
	fChain->SetBranchAddress("AK8_JetPartonID", &AK8_JetPartonID, &b_AK8_JetPartonID);
	fChain->SetBranchAddress("AK8_JetHadFlvr", &AK8_JetHadFlvr, &b_AK8_JetHadFlvr);
	fChain->SetBranchAddress("AK8_JetGenJetIndex", &AK8_JetGenJetIndex, &b_AK8_JetGenJetIndex);
	fChain->SetBranchAddress("AK8_JetGenJetEn", &AK8_JetGenJetEn, &b_AK8_JetGenJetEn);
	fChain->SetBranchAddress("AK8_JetGenJetPt", &AK8_JetGenJetPt, &b_AK8_JetGenJetPt);
	fChain->SetBranchAddress("AK8_JetGenJetEta", &AK8_JetGenJetEta, &b_AK8_JetGenJetEta);
	fChain->SetBranchAddress("AK8_JetGenJetPhi", &AK8_JetGenJetPhi, &b_AK8_JetGenJetPhi);
	fChain->SetBranchAddress("AK8_JetGenPartonID", &AK8_JetGenPartonID, &b_AK8_JetGenPartonID);
	fChain->SetBranchAddress("AK8_JetGenEn", &AK8_JetGenEn, &b_AK8_JetGenEn);
	fChain->SetBranchAddress("AK8_JetGenPt", &AK8_JetGenPt, &b_AK8_JetGenPt);
	fChain->SetBranchAddress("AK8_JetGenEta", &AK8_JetGenEta, &b_AK8_JetGenEta);
	fChain->SetBranchAddress("AK8_JetGenPhi", &AK8_JetGenPhi, &b_AK8_JetGenPhi);
	fChain->SetBranchAddress("AK8_JetGenPartonMomID", &AK8_JetGenPartonMomID, &b_AK8_JetGenPartonMomID);
	fChain->SetBranchAddress("n_AK8SDSJ", &n_AK8SDSJ, &b_n_AK8SDSJ);
	fChain->SetBranchAddress("AK8_SDSJPt", &AK8_SDSJPt, &b_AK8_SDSJPt);
	fChain->SetBranchAddress("AK8_SDSJEta", &AK8_SDSJEta, &b_AK8_SDSJEta);
	fChain->SetBranchAddress("AK8_SDSJPhi", &AK8_SDSJPhi, &b_AK8_SDSJPhi);
	fChain->SetBranchAddress("AK8_SDSJMass", &AK8_SDSJMass, &b_AK8_SDSJMass);
	fChain->SetBranchAddress("AK8_SDSJE", &AK8_SDSJE, &b_AK8_SDSJE);
	fChain->SetBranchAddress("AK8_SDSJCharge", &AK8_SDSJCharge, &b_AK8_SDSJCharge);
	fChain->SetBranchAddress("AK8_SDSJFlavour", &AK8_SDSJFlavour, &b_AK8_SDSJFlavour);
	fChain->SetBranchAddress("AK8_SDSJCSV", &AK8_SDSJCSV, &b_AK8_SDSJCSV);
	fChain->SetBranchAddress("n_AK8puppiSDSJ", &n_AK8puppiSDSJ, &b_n_AK8puppiSDSJ);
	fChain->SetBranchAddress("AK8_puppiSDSJPt", &AK8_puppiSDSJPt, &b_AK8_puppiSDSJPt);
	fChain->SetBranchAddress("AK8_puppiSDSJEta", &AK8_puppiSDSJEta, &b_AK8_puppiSDSJEta);
	fChain->SetBranchAddress("AK8_puppiSDSJPhi", &AK8_puppiSDSJPhi, &b_AK8_puppiSDSJPhi);
	fChain->SetBranchAddress("AK8_puppiSDSJMass", &AK8_puppiSDSJMass, &b_AK8_puppiSDSJMass);
	fChain->SetBranchAddress("AK8_puppiSDSJE", &AK8_puppiSDSJE, &b_AK8_puppiSDSJE);
	fChain->SetBranchAddress("AK8_puppiSDSJCharge", &AK8_puppiSDSJCharge, &b_AK8_puppiSDSJCharge);
	fChain->SetBranchAddress("AK8_puppiSDSJFlavour", &AK8_puppiSDSJFlavour, &b_AK8_puppiSDSJFlavour);
	fChain->SetBranchAddress("AK8_puppiSDSJCSV", &AK8_puppiSDSJCSV, &b_AK8_puppiSDSJCSV);
	Notify();
}

Bool_t lvb_qqqq_dCSV::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void lvb_qqqq_dCSV::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t lvb_qqqq_dCSV::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef lvb_qqqq_dCSV_cxx
