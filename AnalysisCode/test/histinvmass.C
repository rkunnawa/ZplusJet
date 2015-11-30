#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"

using namespace std;
using namespace ROOT;

Float_t getEnergy(Float_t mass, Float_t p){
  return sqrt(pow(mass,2)+pow(p,2));
}

Float_t invmass(Float_t mass, TLorentzVector *muon1, TLorentzVector *muon2){
  
  TVector3 P1 = muon1->Vect();//3 momentum 
  TVector3 P2 = muon2->Vect();
  
  Float_t invmass = sqrt(2*pow(mass,2) + 2*muon1->Energy()*muon2->Energy() - 2*(P1.Dot(P2))); //GeV/c^2
  
  return invmass;
}

void HistInvMass(const int start = 0, const int end = 1){
  
  //Declare all histograms
  TH1F *Hinvmass0 = new TH1F("0", "Invariant Mass Distribution", 200, 0, 10);
  //Hinvmass0->Sumw2();
  TH1F *Hinvmass1 = new TH1F("1", "Invariant Mass Distribution", 200, 0, 101);
  //Hinvmass1->Sumw2();

  bool treedata = false; //true = pure tree data; false = calculated data
  
  //Muon Enhanced Proton Proton High Forest Data
  //
  //https://github.com/ilaflott/For_Ian/blob/master/4_Create_NTuples/filelists/ppMuon_data_filelist.txt
  //Copy this text file into a document called "ppMuon_data_filelist.txt" and save to src
  const string dataFileList = "ppMuon_data_filelist.txt";
  string fileList = dataFileList;
  ifstream fileStream(fileList.c_str(), ifstream::in);
  string fileName;
  fileStream >> fileName;
  
  int w = 0;
  for(; w < start; w++){fileStream>>fileName;} //skips all files before start file
  
  for(; w < end; w++){ //cycles from start file to end file
    cout<<fileName<<endl;
    if((w+1)%10==0){cout<<w+1<<endl;}
    
    TFile *file = TFile::Open(fileName.c_str());
    TTree *pf = (TTree*)file->Get("pfcandAnalyzer/pfTree");
    TTree *track = (TTree*)file->Get("ppTrack/trackTree");
    TTree *skim = (TTree*)file->Get("skimanalysis/HltTree");
    TTree *muon = (TTree*)file->Get("muonTree/HLTMuTree");
    TTree *hi =(TTree*)file->Get("hiEvtAnalyzer/HiTree");
    TTree *hlt = (TTree*)file->Get("hltanalysis/HltTree");
    
    Float_t mass = .1056583715;//mass of muon

    Float_t pfpt[100000];
    Float_t pfeta[100000];
    Float_t pfphi[100000];
    Int_t nPFpart;
    Int_t pfId[100000];
    Int_t charge[100000];
    Int_t pPA;
    Int_t noise;
    Int_t h = 0;
    Int_t i = 0;
    Int_t j = 0;
    Float_t mpt[100000];
    Float_t mrapidity[100000];
    Int_t pairs;
    Float_t mmass[100000];
    Int_t ch1[100000];
    Int_t ch2[100000];
    Float_t pz;
    Float_t vz;
    Float_t dipt1[100000];
    Float_t dipt2[100000];
    Float_t dieta1[100000];
    Float_t dieta2[100000];
    Float_t dieta[100000];
    Float_t diphi1[100000];
    Float_t diphi2[100000];
    Long64_t entries = 0;
    Float_t motherpt;
    Float_t motherrap;
    Float_t imass;
    Int_t mu12;
    Int_t mu7;
    Int_t mu3;
    Float_t Di_mass_sorted;
    
    pf->SetBranchAddress("pfPt",pfpt);
    pf->SetBranchAddress("pfEta",pfeta);
    pf->SetBranchAddress("pfPhi",pfphi);
    pf->SetBranchAddress("nPFpart",&nPFpart);
    pf->SetBranchAddress("pfId",&pfId);
    track->SetBranchAddress("trkCharge",charge);
    skim->SetBranchAddress("pPAcollisionEventSelectionPA", &pPA);
    skim->SetBranchAddress("pHBHENoiseFilter", &noise);
    muon->SetBranchAddress("Di_pt",mpt);
    muon->SetBranchAddress("Di_rapidity",mrapidity);
    muon->SetBranchAddress("Di_npair",&pairs);
    muon->SetBranchAddress("Di_mass",mmass);
    muon->SetBranchAddress("Di_charge1",ch1);
    muon->SetBranchAddress("Di_charge2",ch2);
    muon->SetBranchAddress("Di_pt1",dipt1);
    muon->SetBranchAddress("Di_pt2",dipt2);
    muon->SetBranchAddress("Di_eta",dieta);
    muon->SetBranchAddress("Di_eta1",dieta1);
    muon->SetBranchAddress("Di_eta2",dieta2);
    muon->SetBranchAddress("Di_phi1",diphi1);
    muon->SetBranchAddress("Di_phi2",diphi2);
    hi->SetBranchAddress("vz",&vz);
    hlt->SetBranchAddress("HLT_PAMu12_v1",&mu12);
    hlt->SetBranchAddress("HLT_PAMu7_v1",&mu7);
    hlt->SetBranchAddress("HLT_PAMu3_v1",&mu3);

    pf->AddFriend(muon);
    pf->AddFriend(track);
    pf->AddFriend(skim);
    pf->AddFriend(hi);
    pf->AddFriend(hlt);

    entries = pf->GetEntries();

    //Loop over each event
    for(h = 0; h < entries; h++){//Loop over every particle (first of pair)
      pf->GetEvent(h);
      
      if(pPA == 0 || noise == 0 || abs(vz) > 15)//cuts
	continue;

      if(treedata){//no (or minimal) calculation

	for(Int_t e = 0; e < pairs; e++){//loop through dimuon pairs
	  if(ch1[e] == ch2[e])
	    continue;
	  
	  motherpt = mpt[e];
	  motherrap = mrapidity[e];
	  imass = mmass[e];

	  //Fill all histograms
	  if(imass <= 10)
	    Hinvmass0->Fill(imass);
	  Hinvmass1->Fill(imass);
	}//end dimuon pair loop

      }//end if

      else{//all calculated data
	
	for(i = 0; i < nPFpart; i++){//loop through each particle
	  
	  if(pfId[i] != 3)//muon = 3
	    continue;
	  
	  for(j = i+1; j < nPFpart; j++){//loop through REST of particles; every pair is only checked once
	    
	    if(pfId[j] != 3)
	      continue;
	    if(charge[i] == charge[j])//pairs can only have opposite charge
	      continue;

	    TLorentzVector *muon1 = new TLorentzVector();
	    TLorentzVector *muon2 = new TLorentzVector();
	    muon1->SetPtEtaPhiE(pfpt[i], pfeta[i], pfphi[i], 0);//Pt Eta Phi Energy 4 vector
	    muon2->SetPtEtaPhiE(pfpt[j], pfeta[j], pfphi[j], 0);
	    muon1->SetPtEtaPhiE(pfpt[i], pfeta[i], pfphi[i], getEnergy(mass,muon1->P()));//remake vectors now that we have P
	    muon2->SetPtEtaPhiE(pfpt[j], pfeta[j], pfphi[j], getEnergy(mass,muon2->P()));
	    
	    TLorentzVector *mother = new TLorentzVector();
	    *mother = *muon1 + *muon2;
	    
	    motherpt = mother->Pt();
	    motherrap = mother->Rapidity();
	    imass = invmass(mass, muon1, muon2);
	    
	    //Fill all histograms
	    if(imass <= 10)
	      Hinvmass0->Fill(imass);
	    Hinvmass1->Fill(imass);

	  }//end second particle loop
	  
	}//end first particle loop
	
      }//end else
      
    }//end event loop
    fileStream >> fileName;
  }

  //Write all histograms to file  

  //TFile *plots = new TFile("invmass.root","RECREATE"); //Use if running as macro
  TFile f(Form("Output_%d.root",(start/(end-start))),"RECREATE"); //Use if running as batch job
  Hinvmass0->Write();
  Hinvmass1->Write();
}
