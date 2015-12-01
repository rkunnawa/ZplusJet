// -*- C++ -*-
//
// Package:    ZplusJetAnalyzer
// Class:      ZplusJetAnalyzer
// 
/**\class ZplusJetAnalyzer ZplusJetAnalyzer.cc ZplusJet/AnalysisCode/src/ZplusJetAnalyzer.cc

 Description: 

 Implementation:
     
*/
//
// Original Author:  Raghav Kunnawalkam Elayavalli
//                   Rutgers University, NJ 
//         Created:  Wednesday Oct 7 19:23:23 EST 2015
// $Id$
//
//

#include "ZplusJet/AnalysisCode/plugins/ZplusJetAnalyzer.h"

using namespace std;
using namespace reco;
using namespace edm;

// declare the constructors 

ZplusJetAnalyzer::ZplusJetAnalyzer(const edm::ParameterSet& iConfig) :
  mInputCollection               (iConfig.getParameter<edm::InputTag>       ("jet")),
  JetType                        (iConfig.getUntrackedParameter<std::string>("JetType")),
  UEAlgo                         (iConfig.getUntrackedParameter<std::string>("UEAlgo")),
  radius                         (iConfig.getUntrackedParameter<int>        ("radius")),
  tagRecoMu                      (iConfig.getParameter<edm::InputTag>       ("muons")),
  tagVtx                         (iConfig.getParameter<edm::InputTag>       ("vertices")),
  // mInputPFCandCollection         (iConfig.getParameter<edm::InputTag>       ("PFcands")),
  mRThreshold                    (iConfig.getParameter<double>              ("RThreshold")),
  mRecoJetPtThreshold            (iConfig.getParameter<double>              ("mRecoJetPtThreshold"))
  //JetCorrectionService           (iConfig.getParameter<std::string>         ("JetCorrections")),
  // mhltPath                       (iConfig.getParameter<std::string>         ("hltpath")),
  // triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  // triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  // triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  // triggerEvent_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("hltTrigger"))),
{
  
  std::string inputCollectionLabel(mInputCollection.label());

  isCaloJet = (std::string("calo")==JetType);
  isJPTJet  = (std::string("jpt") ==JetType);
  isPFJet   = (std::string("pf")  ==JetType);
  
  if (isCaloJet) caloJetsToken_  = consumes<reco::CaloJetCollection>(mInputCollection);
  if (isJPTJet)  jptJetsToken_   = consumes<reco::JPTJetCollection>(mInputCollection);
  if (isPFJet)   {
    if(std::string("Pu")==UEAlgo) basicJetsToken_    = consumes<reco::BasicJetCollection>(mInputCollection);
    if(std::string("Vs")==UEAlgo || std::string("pp")==UEAlgo) pfJetsToken_    = consumes<reco::PFJetCollection>(mInputCollection);
  }

  // pfCandToken_ = consumes<reco::PFCandidateCollection>(mInputPFCandCollection);
  // centralityTag_ = iConfig.getParameter<InputTag>("centralitycollection");
  // centralityToken = consumes<reco::Centrality>(centralityTag_);

  // centralityBinTag_ = (iConfig.getParameter<edm::InputTag> ("centralitybincollection"));
  // centralityBinToken = consumes<int>(centralityBinTag_);

}


void ZplusJetAnalyzer::beginJob() {

  // setup the histograms or tree for the next step. it should make a tree with only events which have Z in them and have all the jets in that tree
  TFileDirectory jet = fout->mkdir("ak" + UEAlgo + Form("%d",radius) + JetType + "JetplusZ");
  ZwithJets = jet.make<TTree>("ZplusjetTree","t");
  ZwithJets->Branch("event",&event,"event/I");
  ZwithJets->Branch("run",&run,"run/I");
  ZwithJets->Branch("lumi",&lumi,"lumi/I");
  ZwithJets->Branch("vx",&vx,"vx/F");
  ZwithJets->Branch("vy",&vy,"vy/F");
  ZwithJets->Branch("vz",&vz,"vz/F");
  ZwithJets->Branch("hiBin",&hiBin,"hiBin/I");
  ZwithJets->Branch("mu_p",&mu_p);
  ZwithJets->Branch("mu_pt",&mu_pt);
  ZwithJets->Branch("mu_eta",&mu_eta);
  ZwithJets->Branch("mu_phi",&mu_phi);
  ZwithJets->Branch("mu_charge",&mu_charge);
  ZwithJets->Branch("Z_pt",&Z_pt);
  ZwithJets->Branch("Z_eta",&Z_eta);
  ZwithJets->Branch("Z_phi",&Z_phi);
  ZwithJets->Branch("Z_mass",&Z_mass);
  ZwithJets->Branch("jet_pt",&jet_pt);
  ZwithJets->Branch("jet_eta",&jet_eta);
  ZwithJets->Branch("jet_phi",&jet_phi);
  ZwithJets->Branch("LeadZ_pt",&LeadZ_pt,"LeadZ_pt/F");
  ZwithJets->Branch("LeadZ_eta",&LeadZ_eta,"LeadZ_eta/F");
  ZwithJets->Branch("LeadZ_phi",&LeadZ_phi,"LeadZ_phi/F");
  ZwithJets->Branch("JetinZdir_pt",&JetinZdir_pt,"JetinZdir_pt/F");
  ZwithJets->Branch("JetinZdir_eta",&JetinZdir_eta,"JetinZdir_eta/F");
  ZwithJets->Branch("JetinZdir_phi",&JetinZdir_phi,"JetinZdir_phi/F");
  ZwithJets->Branch("RecoilJet_pt",&RecoilJet_pt,"RecoilJet_pt/F");
  ZwithJets->Branch("RecoilJet_eta",&RecoilJet_eta,"RecoilJet_eta/F");
  ZwithJets->Branch("RecoilJet_phi",&RecoilJet_phi,"RecoilJet_phi/F");
  ZwithJets->Branch("LeadZ_RecoilJet_deltapt",&LeadZ_RecoilJet_deltapt,"LeadZ_RecoilJet_deltapt/F");
  ZwithJets->Branch("LeadZ_RecoilJet_deltaeta",&LeadZ_RecoilJet_deltaeta,"LeadZ_RecoilJet_deltaeta/F");
  ZwithJets->Branch("LeadZ_RecoilJet_deltaphi",&LeadZ_RecoilJet_deltaphi,"LeadZ_RecoilJet_deltaphi/F");
  ZwithJets->Branch("LeadZ_LeadJet_deltapt",&LeadZ_LeadJet_deltapt,"LeadZ_LeadJet_deltapt/F");
  ZwithJets->Branch("LeadZ_LeadJet_deltaeta",&LeadZ_LeadJet_deltaeta,"LeadZ_LeadJet_deltaeta/F");
  ZwithJets->Branch("LeadZ_LeadJet_deltaphi",&LeadZ_LeadJet_deltaphi,"LeadZ_LeadJet_deltaphi/F");
  ZwithJets->Branch("Z_LeadJet_deltaR",&Z_LeadJet_deltaR,"Z_LeadJet_deltaR/F");
  ZwithJets->Branch("Z_RecoilJet_deltaR",&Z_RecoilJet_deltaR,"Z_RecoilJet_deltaR/F");
  ZwithJets->Branch("Z_LeadJet_Aj",&Z_LeadJet_Aj,"Z_LeadJet_Aj/F");
  ZwithJets->Branch("Z_RecoilJet_Aj",&Z_RecoilJet_Aj,"Z_RecoilJet_Aj/F");
  
}



void ZplusJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //Get run, event, centrality
  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  
  edm::Handle< vector<reco::Vertex> > vertex;
  iEvent.getByLabel(tagVtx,vertex);
  if(vertex->size() > 0){
    vx = vertex->begin()->x();
    vy = vertex->begin()->y();
    vz = vertex->begin()->z();
  } else {
    vx = -1;
    vy = -1;
    vz = -1;
  }

  if(fabs(vz) > 15) return;
  
  // get the centrality 
  // edm::Handle<reco::Centrality> cent;
  // iEvent.getByToken(centralityToken, cent);  //_centralitytag comes from the cfg
  // mHF->Fill(cent->EtHFtowerSum());
  // Float_t HF_energy = cent->EtHFtowerSum();  
  // edm::Handle<int> cbin;
  // iEvent.getByToken(centralityBinToken, cbin);
  
  // if(!cent.isValid()) return;
  // int hibin = -999;
  // if(cent.isValid())
  //   hibin = *cbin;
  
  edm::Handle< edm::View<reco::Muon> > muons;
  iEvent.getByLabel(tagRecoMu,muons);

  for (unsigned int i=0; i<muons->size(); i++) {
    edm::RefToBase<reco::Muon> muCand(muons,i);
    if (muCand.isNull()) continue;
    if (muCand->globalTrack().isNonnull() && muCand->innerTrack().isNonnull()) {
      if (muCand->isGlobalMuon() && muCand->isTrackerMuon() && fabs(muCand->combinedMuon()->eta()) < 2.0) {
	//edm::RefToBase<reco::Track> trk = edm::RefToBase<reco::Track>(muCand->innerTrack());
	edm::RefToBase<reco::Track> glb = edm::RefToBase<reco::Track>(muCand->combinedMuon());
	//const reco::HitPattern& p = trk->hitPattern();
	
	//GlbMu.nValMuHits[nGlb] = muCand->combinedMuon().get()->hitPattern().numberOfValidMuonHits();
	//GlbMu.nValTrkHits[nGlb] = muCand->innerTrack().get()->hitPattern().numberOfValidTrackerHits();

	//GlbMu.nTrkFound[nGlb] = trk->found();
	//GlbMu.glbChi2_ndof[nGlb] = glb->chi2()/glb->ndof();
	//GlbMu.trkChi2_ndof[nGlb] = trk->chi2()/trk->ndof();
	//GlbMu.pixLayerWMeas[nGlb] = p.pixelLayersWithMeasurement();
	//GlbMu.trkDxy[nGlb] = fabs(trk->dxy(vertex->begin()->position()));
	//GlbMu.trkDz[nGlb] = fabs(trk->dz(vertex->begin()->position()));

	//muon::SelectionType st = muon::selectionTypeFromString("TrackerMuonArbitrated");
	//GlbMu.isArbitrated[nGlb] = muon::isGoodMuon(*muCand.get(), st);

	// might have to add several quality cuts here:
        mu_charge.push_back(glb->charge());
        mu_pt.push_back(glb->pt());
        mu_p.push_back(glb->p());
        mu_eta.push_back(glb->eta());
        mu_phi.push_back(glb->phi());
	
	//GlbMu.dxy[nGlb] = glb->dxy(vertex->begin()->position());
	//GlbMu.dz[nGlb] = glb->dz(vertex->begin()->position());

	//GlbMu.trkLayerWMeas[nGlb] = muCand->globalTrack()->hitPattern().trackerLayersWithMeasurement();
	//GlbMu.nValPixHits[nGlb] = p.numberOfValidPixelHits();
	//GlbMu.nMatchedStations[nGlb] = muCand->numberOfMatchedStations();

      }

    }

  }// muon loop

  if(mu_p.size() < 2) return;

  // find the Z's from the global muons, also add method to find Z's from PF cand muons and check if they are the same.
  for(unsigned nmu_i = 0; nmu_i<mu_p.size(); ++nmu_i){

    for(unsigned nmu_j = nmu_i+1; nmu_j<mu_p.size(); ++nmu_j){

      if(mu_charge[nmu_i] == mu_charge[nmu_j]) continue; // only opposite pair muons

      TLorentzVector *muon1 = new TLorentzVector();
      TLorentzVector *muon2 = new TLorentzVector();

      muon1->SetPtEtaPhiE(mu_pt[nmu_i], mu_eta[nmu_i], mu_phi[nmu_i], getEnergy(mu_mass, mu_p[nmu_i]));
      muon2->SetPtEtaPhiE(mu_pt[nmu_j], mu_eta[nmu_j], mu_phi[nmu_j], getEnergy(mu_mass, mu_p[nmu_j]));

      TLorentzVector *mother = new TLorentzVector();
      *mother = *muon1 + *muon2;

      float mother_mass = invmass(mu_mass, muon1, muon2);

      if(mother_mass >= 70 && mother_mass <= 110) {
	// we have a Z in the event!  // since this is offshell mass it can be anywhere from 70 - 110 (original mass - 90)
	Z_pt.push_back(mother->Pt());
	Z_eta.push_back(mother->Eta());
	Z_phi.push_back(mother->Phi());
      }
      
    }// second muon loop

  }// muon loop

  if(Z_pt.size() == 0) return;
  
  LeadZ_pt = Z_pt[0];
  LeadZ_eta = Z_eta[0];
  LeadZ_phi = Z_phi[0];
  
  // pf candidates
  // iEvent.getByToken(pfCandToken_, pfCandidates);
  // const reco::PFCandidateCollection *pfCandidateColl = pfCandidates.product();

  // for(unsigned icand=0;icand<pfCandidateColl->size(); icand++){
  //   const reco::PFCandidate pfCandidate = pfCandidateColl->at(icand);
  //   reco::CandidateViewRef ref(pfcandidates_,icand);
  //   if(pfCandidate.pt() < 5) continue;
  // }// pf candidiate loop


  std::vector<Jet> recoJets;
  recoJets.clear();
  
  edm::Handle<CaloJetCollection>  caloJets;
  edm::Handle<JPTJetCollection>   jptJets;
  edm::Handle<PFJetCollection>    pfJets;
  edm::Handle<BasicJetCollection> basicJets;

  if (isCaloJet) iEvent.getByToken(caloJetsToken_, caloJets);
  if (isJPTJet)  iEvent.getByToken(jptJetsToken_, jptJets);
  if (isPFJet) {  
    if(std::string("Pu")==UEAlgo) iEvent.getByToken(basicJetsToken_, basicJets);
    if(std::string("Vs")==UEAlgo) iEvent.getByToken(pfJetsToken_, pfJets);
  }


  if (isCaloJet && !caloJets.isValid()) {
    return;
  }
  if (isJPTJet  && !jptJets.isValid()) {
    return;
  }
  if (isPFJet){
    if(std::string("Pu")==UEAlgo){if(!basicJets.isValid())   return;}
    if(std::string("Vs")==UEAlgo){if(!pfJets.isValid())   return;}
  }

  
  if (isCaloJet){
    for (unsigned ijet=0; ijet<caloJets->size(); ijet++) {
      recoJets.push_back((*caloJets)[ijet]);
    } 
  }

  if (isJPTJet){
    for (unsigned ijet=0; ijet<jptJets->size(); ijet++) 
      recoJets.push_back((*jptJets)[ijet]);
  }

  if (isPFJet) {
    if(std::string("Pu")==UEAlgo){
      for (unsigned ijet=0; ijet<basicJets->size();ijet++) {
	recoJets.push_back((*basicJets)[ijet]);
      }
    }
    if(std::string("Vs")==UEAlgo){
      for (unsigned ijet=0; ijet<pfJets->size(); ijet++){
	recoJets.push_back((*pfJets)[ijet]);
      }
    }
  }

  for (unsigned ijet=0; ijet<recoJets.size(); ijet++) {
    if (recoJets[ijet].pt() > mRecoJetPtThreshold) {

      if(fabs(recoJets[ijet].eta()) > 2.0) continue;
      jet_pt.push_back(recoJets[ijet].pt());
      jet_eta.push_back(recoJets[ijet].eta());
      jet_phi.push_back(recoJets[ijet].phi());

    }

  }// jet loop

  if(jet_pt.size() == 0) return;

  // find the recoil jet and the Delta R, Aj distributions
  float deltaRBest = 999.0;
  int matchJet = 0;
  int recoilJet = 0;
  float recoilDeltaPhi = 999.0;
  
  for(unsigned ijet = 0; ijet<recoJets.size(); ++ijet){
    float delR = deltaR(Z_eta[0], Z_phi[0], recoJets[ijet].eta(), recoJets[ijet].phi());
    if(delR < deltaRBest){
      deltaRBest = delR;
      matchJet = ijet;
    }
    float deltaphi = deltaPhi(Z_phi[0], recoJets[ijet].phi());
    if(deltaphi < recoilDeltaPhi){
      recoilDeltaPhi = deltaphi;
      recoilJet = ijet;
    }
  }

  JetinZdir_pt = recoJets[matchJet].pt();
  JetinZdir_eta = recoJets[matchJet].eta();
  JetinZdir_phi = recoJets[matchJet].phi();

  RecoilJet_pt = recoJets[recoilJet].pt();
  RecoilJet_eta = recoJets[recoilJet].eta();
  RecoilJet_phi = recoJets[recoilJet].phi();

  LeadZ_RecoilJet_deltapt = Z_pt[0] - recoJets[recoilJet].pt();
  LeadZ_RecoilJet_deltaeta = Z_eta[0] - recoJets[recoilJet].eta();
  LeadZ_RecoilJet_deltaphi = Z_phi[0] - recoJets[recoilJet].phi();

  LeadZ_LeadJet_deltapt = Z_pt[0] - recoJets[matchJet].pt();
  LeadZ_LeadJet_deltaeta = Z_eta[0] - recoJets[matchJet].eta();
  LeadZ_LeadJet_deltaphi = Z_phi[0] - recoJets[matchJet].phi();

  Z_LeadJet_deltaR = deltaR(Z_eta[0], Z_phi[0], recoJets[0].eta(), recoJets[1].phi());
  Z_LeadJet_deltaR = deltaR(Z_eta[0], Z_phi[0], recoJets[0].eta(), recoJets[1].phi());

  Z_LeadJet_Aj = (float)(Z_pt[0] - recoJets[0].pt())/(Z_pt[0] + recoJets[0].pt()); 
  Z_RecoilJet_Aj = (float)(Z_pt[0] - recoJets[recoilJet].pt())/(Z_pt[0] + recoJets[recoilJet].pt()); 
  
  ZwithJets->Fill();
}

//void ZplusJetAnalyzer::endJob() {
  // use this as your destructor
//}

//define this as a plug-in
DEFINE_FWK_MODULE(ZplusJetAnalyzer);
