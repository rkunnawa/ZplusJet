#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooChi2Var.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"
#include "TMath.h"
#include "TTree.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"

using namespace std;
using namespace ROOT;
using namespace RooFit;

const double mass_l = 0.5;
const double mass_h = 1.5;

void fitMassBackup(){

  TFile *f = TFile::Open("invmasstree.root","UPDATE");
  TTree *theTree = (TTree*)f->Get("PfDiMuons"); // done.

  //--------- Re-sorting array variable Di_mass[Di_npair] into Di_mass_sorted.
  int nentries = theTree->GetEntries(); cout << nentries << endl;
  float Di_mass_sorted;
  float Di_mass[100000];
  int Di_npair;
  int npair;
  TBranch *b_Di_mass;
  TBranch *b_Di_npair;
  theTree->SetBranchAddress("PfDiMass", &Di_mass, &b_Di_mass);
  theTree->SetBranchAddress("PfIndex", &Di_npair, &b_Di_npair);
  theTree->Branch("Di_mass_sorted", &Di_mass_sorted, "Di_mass_sorted/F");

  for (int i=0; i<nentries; i++) 
    {
      theTree->GetEntry(i);
      npair = Di_npair;
      for(int iQQ =0 ; iQQ < npair ; iQQ++)
	{
	  Di_mass_sorted = Di_mass[iQQ];
	  theTree->Fill();
	}
    }
  //------------- Done. This part will have to be done for all the variables of interest (ID, kinematics, etc.)

  // A - using roofit language, declare the observable (Di_mass) and a set of simple variables to cut on (muon_pt, muon_eta, Di_pt, Di_charge ,... list to be extended to what suits the analysis).
  // notice that I use Di_mass_sorted instead of Di_mass.

  RooRealVar* mass = new RooRealVar("Di_mass_sorted","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}"); 
  RooRealVar* QQsign = new RooRealVar("Di_charge[1]", "Pair electric charge" ,-1,5);
  RooRealVar* mu1Pt = new RooRealVar("Di_pt1[1]","p_{T}^{#mu_{1}}",1.,1000); // arbitrary range, cuts are applied later
  RooRealVar* mu2Pt = new RooRealVar("Di_pt2[1]","p_{T}^{#mu_{2}}",1.,1000); // arbitrary range, cuts are applied later
  RooRealVar* mu1Eta = new RooRealVar("Di_eta1[1]","#eta^{#mu_{1}}",-2.4,2.4); 
  RooRealVar* mu2Eta = new RooRealVar("Di_eta2[1]","#eta^{#mu_{2}}",-2.4,2.4); 

  //making the RooDataSet, subset of the tree, keeping only the set of arguments of interest (the RooRealVars).
  RooDataSet* data0;
  RooArgSet vars(*mass, *QQsign, *mu1Pt, *mu2Pt, *mu1Eta, *mu2Eta); // careful, a RooArgSet has a limited number of arguments (max is 9, as far as I can remember).
  data0 = new RooDataSet("data0","data0",theTree,vars);
  data0->Print();

  //Cutting the RooDataSet to keep only the interesting events (example, muon eta cut between -1.2 and 1.2, and opposite sign dimuons)
  // double muEtaMinimum =-1.2;
  // double muEtaMaximum =1.2;
  // TString cut_ap(Form("(%.2f<Di_eta1[1] && Di_eta1[1] < %.2f) && (%.2f<Di_eta2[1] && Di_eta2[1] < %.2f)",muEtaMinimum, muEtaMaximum,muEtaMinimum, muEtaMaximum)); //  && (%.2f<mu2Eta && mu2Eta < %.2f) && (QQsign==0)

  // for the moment the cut does nothing because of this non-array issue
  TString cut_ap("");
  cout << cut_ap << endl; 
  RooDataSet* data;
  data =  (RooDataSet*) data0->reduce(Cut(cut_ap));
  data->Print();

  // B - Making the fit model
  //1. The parameters of interest = N_phi, N_omega, eventually add N_rho but there's nothing there if you use only one file.
  RooRealVar *n_phi = new RooRealVar("N_{ #phi(1020)}","n_phi",0,1000);
  RooRealVar *n_omega = new RooRealVar("N_{ #omega}","n_phi",0,1000);
  //2. The other fit parameters = masses, widths,
  RooRealVar *mean_phi = new RooRealVar("m_{ #phi(1020)}","#phi(1020) mean",1.02,0.95,1.1); //restrictive range, can be tuned .
  RooRealVar *mean_omega = new RooRealVar("m_{ #omega}","#omega mean",0.782,0.7,0.9); //restrictive range, can be tuned .
  // widths (for light vector mesons in cms expect about 12 MeV at high pt).
  RooRealVar *sigma_phi = new RooRealVar("#sigma_{#phi(1020)}","#sigma_{#phi(1020)}",0.005,0.02); 
  RooRealVar *sigma_omega = new RooRealVar("#sigma_{#omega}","#sigma_{#omega}",0.005,0.02); // 
  // 3. Nuisance parameters : parameters of the CB description (assume they are the same for phi and omega, so that you dont have too many fit params.)
  RooRealVar *alpha  = new RooRealVar("#alpha_{CB}","tail shift",0.3,6); // can study with MC
  RooRealVar *npow   = new RooRealVar("n_{CB}","power order",1,15);      // can study with MC
  // 4. The signal PDF : use CB (can try gaussians as well)
  RooCBShape  *sig_phi   = new RooCBShape ("cb_phi", "cb_phi",
					   *mass,*mean_phi,*sigma_phi,*alpha,*npow);

  RooCBShape  *sig_omega   = new RooCBShape ("cb_omega", "cb_omega",
					     *mass,*mean_omega,*sigma_omega,*alpha,*npow);
  // 5. Background: N_events, parameters, pdf... For simplicity, try first a Polynomial of order 2. 
  RooRealVar *nbkgd   = new RooRealVar("n_{Bkgd}","nbkgd",0,10000); // if your sample has N events, might be worth using this number instead of 10000...
  RooRealVar *bkg_a1  = new RooRealVar("a1_bkg", "bkg_{a1}", 0, -5, 5); 
  RooRealVar *bkg_a2  = new RooRealVar("a2_Bkg", "bkg_{a2}", 0, -2, 2);
  RooAbsPdf  *ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
					  *mass, RooArgList(*bkg_a1,*bkg_a2));
  // 6. Finally, build the total fit pdf = sig + bkgd
  RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total p.d.f.",
					       RooArgList(*sig_phi,*sig_omega,*ChebPdf),
					       RooArgList(*n_phi,*n_omega,*nbkgd));

  // C - Do the fit, and print the output in the terminal prompt.
  RooFitResult* fitObject = pdf->fitTo(*data,Save(),Extended(kTRUE)); // Fit options to be found online or in the user guide.
  pdf->Print();
  fitObject->Print("v");

  //Draw a plot of the data+fit (here, the data will be binned, but in the fit the data is not binned (thanks to roofit))
  int nbins = (mass_h-mass_l)/0.01; // e.g. 100 bins of size 0.01 GeV/c^2 between 0.5 and 1.5 GeV/c^2
  RooPlot* frame = mass->frame(Bins(nbins),Range(mass_l,mass_h)); 
  data->plotOn(frame);// data drawn first for pdf object to pick the proper normalisation!
  data0->plotOn(frame);// data drawn first for pdf object to pick the proper normalisation!
  pdf->plotOn(frame,Name("thePdf"));
  pdf->plotOn(frame,Components("ChebPdf"),Name("theBkg"),LineStyle(kDashed)); // look for options online or in the user guide.
  data->plotOn(frame) ;// drawing data points again over pdf line (looks best).
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  frame->GetXaxis()->CenterTitle(kTRUE);
  TCanvas* cm =  new TCanvas("cm", "cm",423,55,600,600);
  cm->cd();
  pdf->paramOn(frame,Layout(0.6,0.935,0.97),Format("NEAU",AutoPrecision(1)));
  frame->Draw();
}
