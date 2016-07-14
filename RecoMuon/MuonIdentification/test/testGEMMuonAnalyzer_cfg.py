import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("GEMMuonAnalyzer", eras.Phase2GReco)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023tiltedReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023tilted_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load( "DQMServices/Core/DQMStore_cfg" )
process.load('Validation.RecoMuon.associators_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories.append("GEMMuonAnalyzer")
process.MessageLogger.categories.append("GEMMuonAnalyzer_Matching")
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.destinations = cms.untracked.vstring("cout","junk")
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
      '/store/user/jskim/TenMu_Pt_5_100_step3_PU/step3_RAW2DIGI_L1Reco_RECO_0.root'
    ),
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
  UseGEMEtaCoverageMuons = cms.bool(True),
  UseAssociators = cms.bool(True),

  associators = cms.vstring('gemMuonAssociatorByHits'),

  label = cms.vstring('gemMuonSel'),

  doMatchingStudy = cms.bool(True),
  maxPull = cms.vdouble(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 10.0, 100.0, 500.0, 1000.0, 10000.0, 100000.0),
  maxX_GE11 = cms.vdouble(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 10.0, 100.0, 500.0, 1000.0, 10000.0, 100000.0),
  maxY_GE11 = cms.vdouble(0.001, 0.01, 0.1, 1.0, 2.0, 3.0, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 100.0, 200., 300., 400., 500., 1000.),
  maxX_GE21 = cms.vdouble(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 10.0, 100., 200., 300., 400., 500., 1000.0, 10000.0),
  maxY_GE21 = cms.vdouble(0.001, 0.01, 0.1, 1.0, 2.0, 3.0, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 100.0, 200., 300., 400., 500., 1000.),
  minDotDir = cms.vdouble(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0),

  Current_trackerGEM_maxPull = cms.double (2.0),
  Current_maxDiffXGE11 = cms.double (1.5),
  Current_maxDiffYGE11 = cms.double (10.0),
  Current_maxDiffXGE21 = cms.double (2.5),
  Current_maxDiffYGE21 = cms.double (12.0),
  Current_minDotDir    = cms.double (0.9),
  
)

process.p = cms.Path(process.gemMuonSel*process.muonAssociatorByHits*process.GEMMuonAnalyzer)


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
