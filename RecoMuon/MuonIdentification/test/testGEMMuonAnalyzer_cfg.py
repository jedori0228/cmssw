import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("GEMMuonAnalyzer",eras.Phase2C2_timing)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load( "DQMServices/Core/DQMStore_cfg" )
process.load('Validation.RecoMuon.associators_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2023_realistic_v1', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append("GEMMuonAnalyzer")
process.MessageLogger.categories.append("GEMMuonAnalyzer_Matching")
process.MessageLogger.destinations = cms.untracked.vstring("cout","junk")
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string("DEBUG"),
    default = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
    FwkReport = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    #GEMMuonAnalyzer   = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
    #GEMMuonAnalyzer_Matching = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
)

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:///somewhere/simevent.root') ##/somewhere/simevent.root" }
    fileNames = cms.untracked.vstring(
      '/store/relval/CMSSW_9_0_0_pre1/RelValZMM_14/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v0_2023D4PU200r2-v1/10000/02941085-E7BD-E611-9A42-0CC47A4D7604.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True), 
    #skipEvents=cms.untracked.uint32(1)
)

# Configurations for MuonTrackValidators


# Configurations for RecoMuonValidators

#import SimGeneral.MixingModule.mixNoPU_cfi
#process.TrackAssociatorByChi2ESProducer = Validation.RecoMuon.associators_cff.TrackAssociatorByChi2ESProducer.clone(chi2cut = 100.0,ComponentName = 'TrackAssociatorByChi2')


#import SimMuon.MCTruth.muonAssociatorByHitsHelper_cfi
#process.muonAssociatorByHits = SimMuon.MCTruth.muonAssociatorByHitsHelper_cfi.muonAssociatorByHitsHelper.clone(#ComponentName = "muonAssociatorByHits",
# #tpTag = 'mix:MergedTrackTruth',
# UseTracker = True,
# UseMuon = False,
# EfficiencyCut_track = cms.double(0.0),
# PurityCut_track = cms.double(0.0)
#)

from CommonTools.RecoAlgos.gemAssociator import *
process.gemMuonSel = gemmuon.clone()

import SimMuon.MCTruth.MuonAssociatorByHits_cfi
process.gemMuonAssociatorByHits = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("gemMuonSel"),
)

process.GEMMuonAnalyzer = cms.EDAnalyzer("GEMMuonAnalyzer",

  HistoFile = cms.string('OUTPUTTEMPLAT.root'),

  FakeRatePtCut = cms.double(5.0),
  MatchingWindowDelR = cms.double (.3),
  UseAssociators = cms.bool(True),
  UseDeltaR = cms.bool(False),
  doGeometryStudy = cms.bool(False),
  SampleProcess = cms.string('MuonGun'),

  doMatchingStudy = cms.bool(False),
  associators = cms.vstring('gemMuonAssociatorByHits'),
  label = cms.vstring('gemMuonSel')

)

process.p = cms.Path(
process.gemMuonSel*process.gemMuonAssociatorByHits
*process.GEMMuonAnalyzer
)





############## To make a output root file ###############

#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.output = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string(
#        'file:out_test.root'
#    ),
#    outputCommands = cms.untracked.vstring(
#        'keep  *_*_*_*',
#    ),
#)
#process.out_step     = cms.EndPath(process.output)
#process.p = cms.Path(process.gemMuonSel*process.muonAssociatorByHits*process.GEMMuonAnalyzer)
#process.p = cms.Path(process.gemMuonSel)
#process.out_step = cms.EndPath(process.output)
#process.schedule = cms.Schedule(
#    process.p,
#    process.out_step
#)
