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
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TH1F.h>
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

using namespace std;

void PrintInvMass(){

  TFile *file = TFile::Open("invmass.root"); //Open the files where histograms are stored
  
  TH1::SetDefaultSumw2();
  
  TH1F *H0 = (TH1F*)file->Get("0"); //Get the histograms
  TH1F *H1 = (TH1F*)file->Get("1");
  TH1F *H2 = (TH1F*)file->Get("2");
  
  //Format Histograms
  TCanvas *IM0 = new TCanvas("c0", "Invariant Mass Distribution", 600, 600);
  
  IM0->SetLogy();
  H0->SetStats(0);
  H0->SetXTitle("Invariant Mass (GeV/c^2)");
  H0->SetYTitle("Count");
  H0->Draw();

  TCanvas *IM1 = new TCanvas("c1", "Invariant Mass Distribution", 600, 600);

  IM1->SetLogy();
  H1->SetStats(0);
  H1->SetXTitle("Invariant Mass (GeV/c^2)");
  H1->SetYTitle("Count");
  H1->Draw();

  TCanvas *IM2 = new TCanvas("c2", "Invariant Mass Distribution", 600, 600);

  IM2->SetLogy();
  H2->SetXTitle("Invariant Mass (GeV/c^2)");
  H2->SetYTitle("Count");
  H2->Draw();

  //Save canvases to .pdf
  IM0->Print("IMassPlots.pdf(","pdf");
  IM1->SaveAs("IMassPlots.pdf","pdf");
  IM2->SaveAs("IMassPlots.pdf)","pdf");
}
