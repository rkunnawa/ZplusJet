# config file to produce the zplus jet analysis. 
import FWCore.ParameterSet.Config as cms

#from FWCore.ParameterSet.VarParsing import VarParsing
###############################
####### Parameters ############
###############################
#options = VarParsing ('python')
#options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
	                    fileNames = cms.untracked.vstring('file:AOD_data.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.TFileService=cms.Service(
    "TFileService",
    fileName=cms.string('testoutput_ZplusJet.root')
)

process.ZplusAKPu4PF = cms.EDAnalyzer(
    'ZplusJetAnalyzer',
    src = cms.InputTag("akVs4PFJets"),
    JetType = cms.untracked.string('pf'),
    UEAlgo = cms.untracked.string('Pu'),
    radius = cms.int(4),
    muons = cms.InputTag("muons"),
    vertices = cms.InputTag("hiSelectedVertex"),
    #centralitycollection = cms.InputTag("hiCentrality"),
    #centralitybincollection = cms.InputTag("centralityBin","HFtowers"),
    #JetCorrections = cms.string(""),
    #PFcands = cms.InputTag("particleFlowTmp"),
    RThreshold = cms.double(0.3),  
    mRecoJetPtThreshold = cms.double(10)        
)


process.ZplusAKVs4PF = process.ZplusAKVs4PF.clone(
    UEAlgo = cms.untracked.string('Pu')
)


process.p = cms.Path(
    process.ZplusAKPu4PF *
    process.ZplusAKVs4PF
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
