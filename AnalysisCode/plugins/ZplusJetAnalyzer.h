#ifndef ZplusJetAnalyzer_H
#define ZplusJetAnalyzer_H


//
//  
// For CMSSW_7_5_3_patch1,  
// author: Raghav Kunnawalkam Elayavalli,
//         Oct 20th 2015 
//         Rutgers University, email: raghav.k.e at CERN dot CH 
//
//  
//

#include <memory>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandidateWithRef.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// include the basic jet for the PuPF jets. 
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
// include the pf candidates 
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
// include the voronoi subtraction
#include "DataFormats/HeavyIonEvent/interface/VoronoiBackground.h"
#include "RecoHI/HiJetAlgos/interface/UEParameters.h"
// include the centrality variables
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

// root headers

#include <TH1.h>
#include <TH2.h>
#include "TROOT.h"
#include "TTree.h"
#include "TLorentzVector.h"

const Float_t mu_mass = .1056583715;//mass of muon

using namespace edm;
using namespace std;

class ZplusJetAnalyzer : public edm::EDAnalyzer {

 public:
  explicit ZplusJetAnalyzer(const edm::ParameterSet&);
  ~ZplusJetAnalyzer() {}

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  //virtual void endJob() override;


  // Here are all the necessary jet collections, reco muons etc...  
  edm::InputTag mInputCollection;
  std::string JetType;
  std::string UEAlgo;
  int radius;
  edm::InputTag tagRecoMu;
  edm::InputTag tagVtx;
  //edm::EDGetTokenT<std::vector<reco::Vertex>>                vertexToken_;
  // edm::InputTag mInputPFCandCollection;
  double mRThreshold;
  double mRecoJetPtThreshold;
  //std::string JetCorrectionService;
  /* std::string mhltPath; */
  /* edm::EDGetTokenT<edm::TriggerResults> triggerBits_; */
  /* edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_; */
  /* edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_; */
  /* edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent_; */
  
  //edm::EDGetTokenT<reco::JetCorrector> jetCorrectorToken_;
  //int triggerBit;    
  
  bool isCaloJet;
  bool isJPTJet;
  bool isPFJet;
  edm::EDGetTokenT<reco::CaloJetCollection> caloJetsToken_;
  edm::EDGetTokenT<reco::PFJetCollection> pfJetsToken_;
  edm::EDGetTokenT<reco::BasicJetCollection> basicJetsToken_;
  edm::EDGetTokenT<reco::JPTJetCollection> jptJetsToken_;
  /* edm::EDGetTokenT<reco::PFCandidateCollection> pfCandToken_;  */

  /* edm::InputTag centralityTag_; */
  /* edm::EDGetTokenT<reco::Centrality> centralityToken; */
  /* edm::Handle<reco::Centrality> centrality_; */
  /* edm::InputTag centralityBinTag_; */
  /* edm::EDGetTokenT<int> centralityBinToken; */
  /* edm::Handle<int>centralityBin_; */
  
  // variables for the output root files
  edm::Service<TFileService> fout;
  // map< string, TH1D* > histos1D;  
  // map< string, TH2D* > histos2D;
  TTree * ZwithJets;
  // this is the basic analysis tree containing the following 

  int event, run, lumi, hiBin;
  double vz, vx, vy;
  // vectors of items per event - pt, eta and phi and maass for the Z
  // Muons - 
  std::vector <float> mu_p; 
  std::vector <float> mu_pt; 
  std::vector <float> mu_eta; 
  std::vector <float> mu_phi; 
  std::vector <int> mu_charge; 
  // Z -
  std::vector <float> Z_pt; 
  std::vector <float> Z_eta; 
  std::vector <float> Z_phi;
  std::vector <float> Z_mass;
  // Jets - 
  std::vector <float> jet_pt;  
  std::vector <float> jet_eta;  
  std::vector <float> jet_phi;  

  // leading Z information
  float LeadZ_pt;
  float LeadZ_eta;
  float LeadZ_phi;
  // jet in Z direction;
  float JetinZdir_pt;
  float JetinZdir_eta;
  float JetinZdir_phi;
  // recoil jet information; this is the jet that has delta phi  with Z > 2pi/3 
  float RecoilJet_pt;
  float RecoilJet_eta;
  float RecoilJet_phi;
  // difference between Recoil Jet and lead Z 
  float LeadZ_RecoilJet_deltapt;
  float LeadZ_RecoilJet_deltaeta;
  float LeadZ_RecoilJet_deltaphi;
  // difference between lead Jet and lead Z 
  float LeadZ_LeadJet_deltapt;
  float LeadZ_LeadJet_deltaeta;
  float LeadZ_LeadJet_deltaphi;

  float Z_LeadJet_deltaR; // for jets along the Z direction
  float Z_RecoilJet_deltaR;// for jets opposite the Z direction

  float Z_LeadJet_Aj;
  float Z_RecoilJet_Aj;
  
  // methods to calculate the invariant mass and energy. 

  static Float_t getEnergy(Float_t mass, Float_t p){
    return sqrt(pow(mass,2)+pow(p,2));
  }
  
  static Float_t invmass(Float_t mass, TLorentzVector *muon1, TLorentzVector *muon2){
  
    TVector3 P1 = muon1->Vect();//3 momentum 
    TVector3 P2 = muon2->Vect();
  
    Float_t invmass = sqrt(2*pow(mass,2) + 2*muon1->Energy()*muon2->Energy() - 2*(P1.Dot(P2))); //GeV/c^2
  
    return invmass;
  }

  static float deltaPhi(float phi1, float phi2){
    float dphi = fabs(phi1 - phi2);
    if(dphi > M_PI )dphi -= 2*M_PI;
    return dphi;
  }
  
};

#endif
