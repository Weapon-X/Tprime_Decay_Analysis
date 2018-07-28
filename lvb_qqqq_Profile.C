
#define    lvb_qqqq_new_cxx
#include  "lvb_qqqq_new.h"
#include  <TH2.h>
#include  <TStyle.h>
#include  <TCanvas.h>
#include  <map>
using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2){ 
		std::cout<<"pass the text file containing root file names "<<std::endl;
		exit(0);
	}
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <map>");
	lvb_qqqq_dCSV  a(argv[1]);
	TString InputTxtFile = argv[1];
	TString OutputFileName = InputTxtFile.ReplaceAll(".txt",".root");
	a.Loop(OutputFileName.Data());

}

void lvb_qqqq_dCSV::Loop(TString OutputFileName)
{
	if (fChain == 0) return;
	int prfile = 
	//   TFile* f2 = new TFile(OutputFileName.Data(),"recreate");  
	TFile* f2 = new TFile("TbtHSignal1500GeV.root","recreate");  
	Long64_t nentries = fChain->GetEntriesFast();
	cout <<"\nTotal Events = " << nentries <<endl ;

	//===========Histogram Functions ==============================	

	/*	DefineMC_NPtEta_Histo() ;	 
		dRHisto_MCRecoObject() ;   	
		dRHisto_MCObject() ;    
		Define_Mt_Histo() ;  */

	Define_NPtEta_Histo();
	dRHisto_RecoObject();
	Define_Tag_Jet_Histo() ;  
	Define_Reco_tagjetHisto() ;   
	dR_tagjetHisto() ;     
	Category_Object_Histo() ;    
	Category_Object_dRHisto() ;     
	Category_Object_MtHisto() ;  

	//==========================Event Study ===================

	Long64_t nbytes = 0, nb = 0;
	int j, st = 0 , muon = 0, bjet = 0 , AK8 = 0, mc = 0, met = 0, bquark = 0, mu_bjet = 0, top = 0;
	int b_asso, q_forw, topbjet = -1;

	int b_match = 0;
	int qj = -1 ;
	float dR = 0.0 ;
	float dR1 = 0.0 ;
	double pz ;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		//  for (Long64_t jentry=0; jentry<5000;jentry++) { }
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		//if (ientry == 200) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		Clear_Vector() ;
		j = -1;
		b_asso = -1 ;
		b_top = -1 ;
		q_forw = -1;

		if (jentry ==  nentries-1 )  cout<<"\n TtH(MC) = "<< st << ", Muon pass = " << muon << ", bjet pass = " << bjet << ", AK8 pass "    << AK8 << ", mc = "  << mc <<", met = " << met << ", bquark = " << bquark << ", mu_bjet = " << mu_bjet << ", top = " << top << ", event top W = " << eventI << ", event top fat = " << eventII << " event Higgsjet = " << eventVI << " event WW = " << eventV <<endl  ;

		//===========================MC info========================
		T_top = -1;
		b_top = -1 ;
		T_higgs = -1;
		Higgs_W = -1 ;
		int check =0;

		for(int h =0; h< n_MC; h++){
			if (abs((*mc_PID)[h]) == 6  && abs((*mc_MomPID)[h]) == 8000001 && (*mc_Status)[h] <= 25 ) T_top = h ;

			if ( abs((*mc_PID)[h]) == 25 && abs((*mc_MomPID)[h]) == 8000001 && (*mc_Status)[h] <= 25 )  T_higgs = h ;

			if( abs((*mc_PID)[h]) == 24 && abs((*mc_MomPID)[h]) == 6 && abs((*mc_GMomPID)[h]) == 8000001) Higgs_W = h ;

			if( abs((*mc_PID)[h]) == 13  && abs((*mc_MomPID)[h]) == 24 && abs((*mc_GMomPID)[h]) == 6) top_Mu.push_back(h);

			if( abs((*mc_PID)[h]) == 13  && abs((*mc_MomPID)[h]) == 24 && abs((*mc_GMomPID)[h]) == 25 ) Higgs_Mu.push_back(h);

			if(abs((*mc_PID)[h]) <= 4 &&  (*mc_Status)[h] <= 25 && abs((*mc_MomPID)[h]) == 24 && abs((*mc_GMomPID)[h]) == 25)  Higgs_Jet.push_back(h);

			if(abs((*mc_PID)[h]) <= 4 &&  (*mc_Status)[h] <= 25 && abs((*mc_MomPID)[h]) == 24 && abs((*mc_GMomPID)[h]) == 6)    topW_q.push_back(h);

			if(abs((*mc_PID)[h]) == 5 &&  (*mc_Status)[h] <= 25 && (abs((*mc_MomPID)[h]) == 21 || abs((*mc_MomPID)[h]) <= 4))  b_asso = h ;

			if(abs((*mc_PID)[h]) == 5 &&  (*mc_Status)[h] <= 25 && abs((*mc_MomPID)[h]) == 6 )  b_top = h ;

			if(abs((*mc_PID)[h]) < 5 &&  (*mc_Status)[h] <= 25 && (abs((*mc_MomPID)[h]) == 21 || abs((*mc_MomPID)[h]) <= 4))  q_forw = h ;

		}

		/*    if (topW_q.size() == 0  ) continue ;   
		      if ( Higgs_Jet.size() == 0 ) continue ; 
		      if ( top_Mu.size() == 0) continue ;
		      if ( Higgs_Mu.size() == 0 ) continue ;
		      mc ++ ;
		      if ( T_top == -1 ) continue ;
		      if ( T_higgs == -1 ) continue ;
		      if ( b_top  == -1 ) continue ; 
		 */

		//=====================   Object selections  ================
		//======================AK8 jet selection ==================
		for(int f = 0 ; f < N_AK8Jet ; f ++)
		{
			Cut_AK8jet(f)  ;
		}
		if( n_AK8Jet.size() == 0) continue ;
		AK8 ++ ;
		//======================muonselection==========================
		for( int f = 0; f < N_Mu ; f++)
		{
			if (Cut_Muon(f) ) n_Mu.push_back(f);
		}

		if ( n_Mu.size() == 0    ) continue;		
		j = n_Mu[0] ;
		//int sizemu = (n_Mu.size() >= 2) ? 2 :  n_Mu.size() ; 
		muon ++ ;
		//========================================================
		for( int f = 0 ; f < n_Jet ; f++ )
		{
			Cut_bjet(f,1) ;   // 1 for loose, 2 for medium, 3 for tight
		}

		Fill_RecoObject() ;
		// ==================AK8 jetstudy==============================
		int size8  =  n_AK8Jet.size() ;
		int g2 = -1 ; 
		int lvl = -1 ;
		int b2 = -1;
		int t2 = -1 ;
		//if (size8 < 4) continue ;

		Wjet_selection();
		Top_selection();
		Higgs_selection() ; 
		Fatjet_selection() ;
		topbjet = -1;

		Topjet_Plots() ;		
		Higgsjet_Plots(1) ;
		Wjet_Plots(-1);
		TagJets_dRPlots() ;
		Fatjet_Plots(-1) ;

		//========= Signal Category Study starts from here========================

		if (W_boson.size() == 0)  {
			Wtag0_Category() ; // 0 W tag Category
		}

		if (W_boson.size() == 1) Wtag1_Category() ;				

		if (W_boson.size() != 2)   continue ; //
		if (b_jet.size() != 0 )  {
			WW_lvbjet_Plots() ;
			eventV ++ ;
		}		                 		     
		//==============================

	}

	f2->Write();
	f2->Close();

}
