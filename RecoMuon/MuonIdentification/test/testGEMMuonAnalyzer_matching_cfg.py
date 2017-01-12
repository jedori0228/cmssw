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
    input = cms.untracked.int32(1000)
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

process.load('RecoMuon.MuonIdentification.trackerGEM_cfi')

process.PullXScan0 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullX = cms.double(0.),
)
process.PullXScan0p1 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullX = cms.double(0.1),
)
process.PullXScan0p5 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullX = cms.double(0.5),
)
process.PullXScan1 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullX = cms.double(1.),
)
process.PullXScan2 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullX = cms.double(2.),
)
process.PullXScan3 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullX = cms.double(3.),
)

process.DXScan0 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxDX = cms.double(0.),
)
process.DXScan0p5 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxDX = cms.double(0.5),
)
process.DXScan1 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxDX = cms.double(1.),
)
process.DXScan2 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxDX = cms.double(2.),
)
process.DXScan3 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxDX = cms.double(3.),
)
process.DXScan4 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxDX = cms.double(4.),
)

process.PullYScan0p5 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(0.5),
  MaxDY = cms.double(-1.),
)
process.PullYScan1 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(1.),
  MaxDY = cms.double(-1.),
)
process.PullYScan2 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(2.),
  MaxDY = cms.double(-1.),
)
process.PullYScan3 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(3.),
  MaxDY = cms.double(-1.),
)
process.PullYScan4 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(4.),
  MaxDY = cms.double(-1.),
)
process.PullYScan5 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(5.),
  MaxDY = cms.double(-1.),
)
process.PullYScan10 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(10.),
  MaxDY = cms.double(-1.),
)
process.PullYScan20 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(20.),
  MaxDY = cms.double(-1.),
)

process.DYScan1 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(-1.),
  MaxDY = cms.double(1.),
)
process.DYScan2 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(-1.),
  MaxDY = cms.double(2.),
)
process.DYScan3 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(-1.),
  MaxDY = cms.double(3.),
)
process.DYScan4 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(-1.),
  MaxDY = cms.double(4.),
)
process.DYScan5 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(-1.),
  MaxDY = cms.double(5.),
)
process.DYScan10 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(-1.),
  MaxDY = cms.double(10.),
)
process.DYScan20 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MaxPullY = cms.double(-1.),
  MaxDY = cms.double(20.),
)

process.DotDirScan0p99 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MinDotDir = cms.double(0.99),
)
process.DotDirScan0p95 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MinDotDir = cms.double(0.95),
)
process.DotDirScan0p90 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MinDotDir = cms.double(0.90),
)
process.DotDirScan0p85 = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  MinDotDir = cms.double(0.85),
)



import SimMuon.MCTruth.MuonAssociatorByHits_cfi
process.gemMuonAsso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("gemMuonSel"),
)
process.PullXScan0Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullXScan0"),
)
process.PullXScan0p1Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullXScan0p1"),
)
process.PullXScan0p5Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullXScan0p5"),
)
process.PullXScan1Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullXScan1"),
)
process.PullXScan2Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullXScan2"),
)
process.PullXScan3Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullXScan3"),
)

process.DXScan0Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DXScan0"),
)
process.DXScan0Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DXScan0"),
)
process.DXScan0p5Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DXScan0p5"),
)
process.DXScan1Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DXScan1"),
)
process.DXScan2Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DXScan2"),
)
process.DXScan3Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DXScan3"),
)
process.DXScan4Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DXScan4"),
)

process.PullYScan0p5Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan0p5"),
)
process.PullYScan1Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan1"),
)
process.PullYScan2Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan2"),
)
process.PullYScan3Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan3"),
)
process.PullYScan4Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan4"),
)
process.PullYScan5Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan5"),
)
process.PullYScan10Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan10"),
)
process.PullYScan20Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("PullYScan20"),
)

process.DYScan1Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DYScan1"),
)
process.DYScan2Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DYScan2"),
)
process.DYScan3Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DYScan3"),
)
process.DYScan4Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DYScan4"),
)
process.DYScan5Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DYScan5"),
)
process.DYScan10Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DYScan10"),
)
process.DYScan20Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DYScan20"),
)

process.DotDirScan0p99Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DotDirScan0p99"),
)
process.DotDirScan0p95Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DotDirScan0p95"),
)
process.DotDirScan0p90Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DotDirScan0p90"),
)
process.DotDirScan0p85Asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 useGEMs = cms.bool(True),
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
 pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
 stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
 tracksTag = cms.InputTag("DotDirScan0p85"),
)

process.GEMMuonAnalyzer = cms.EDAnalyzer("GEMMuonAnalyzer",

  HistoFile = cms.string('OUTPUTTEMPLAT.root'),

  FakeRatePtCut = cms.double(5.0),
  MatchingWindowDelR = cms.double (.3),
  UseAssociators = cms.bool(True),
  UseDeltaR = cms.bool(False),
  doGeometryStudy = cms.bool(False),
  SampleProcess = cms.string('MuonGun'),

  doMatchingStudy = cms.bool(True),
  PullXValues = cms.vdouble(0, 0.1, 0.5, 1.0, 2.0, 3.0),
  DXValues = cms.vdouble(0, 0.5, 1.0, 2.0, 3.0, 4.0),
  PullYValues = cms.vdouble(0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 10., 20.),
  DYValues = cms.vdouble(1.0, 2.0, 3.0, 4.0, 5.0, 10., 20.),
  DotDirValues = cms.vdouble(0.99, 0.95, 0.90, 0.85),
  associators = cms.vstring(
    'gemMuonAsso',
    'PullXScan0Asso', 'PullXScan0p1Asso', 'PullXScan0p5Asso', 'PullXScan1Asso' ,'PullXScan2Asso', 'PullXScan3Asso',
    'DXScan0Asso', 'DXScan0p5Asso', 'DXScan1Asso', 'DXScan2Asso', 'DXScan3Asso', 'DXScan4Asso',
    'PullYScan0p5Asso', 'PullYScan1Asso', 'PullYScan2Asso', 'PullYScan3Asso' ,'PullYScan4Asso', 'PullYScan5Asso', 'PullYScan10Asso', 'PullYScan20Asso',
    'DYScan1Asso', 'DYScan2Asso', 'DYScan3Asso' ,'DYScan4Asso', 'DYScan5Asso', 'DYScan10Asso', 'DYScan20Asso',
    'DotDirScan0p99Asso', 'DotDirScan0p95Asso', 'DotDirScan0p90Asso', 'DotDirScan0p85Asso',
  ),
  label = cms.vstring(
    'gemMuonSel',
    'PullXScan0', 'PullXScan0p1', 'PullXScan0p5', 'PullXScan1', 'PullXScan2', 'PullXScan3',
    'DXScan0', 'DXScan0p5', 'DXScan1', 'DXScan2', 'DXScan3', 'DXScan4',
    'PullYScan0p5', 'PullYScan1', 'PullYScan2', 'PullYScan3' ,'PullYScan4', 'PullYScan5', 'PullYScan10', 'PullYScan20',
    'DYScan1', 'DYScan2', 'DYScan3' ,'DYScan4', 'DYScan5', 'DYScan10', 'DYScan20',
    'DotDirScan0p99', 'DotDirScan0p95', 'DotDirScan0p90', 'DotDirScan0p85', 
  ),

)

process.p = cms.Path(
process.trackerGEM*
process.gemMuonSel*process.gemMuonAsso*
process.PullXScan0*process.PullXScan0Asso*
process.PullXScan0p1*process.PullXScan0p1Asso*
process.PullXScan0p5*process.PullXScan0p5Asso*
process.PullXScan1*process.PullXScan1Asso*
process.PullXScan2*process.PullXScan2Asso*
process.PullXScan3*process.PullXScan3Asso*
process.DXScan0*process.DXScan0Asso*
process.DXScan0p5*process.DXScan0p5Asso*
process.DXScan1*process.DXScan1Asso*
process.DXScan2*process.DXScan2Asso*
process.DXScan3*process.DXScan3Asso*
process.DXScan4*process.DXScan4Asso*
process.PullYScan0p5*process.PullYScan0p5Asso*
process.PullYScan1*process.PullYScan1Asso*
process.PullYScan2*process.PullYScan2Asso*
process.PullYScan3*process.PullYScan3Asso*
process.PullYScan4*process.PullYScan4Asso*
process.PullYScan5*process.PullYScan5Asso*
process.PullYScan10*process.PullYScan10Asso*
process.PullYScan20*process.PullYScan20Asso*
process.DYScan1*process.DYScan1Asso*
process.DYScan2*process.DYScan2Asso*
process.DYScan3*process.DYScan3Asso*
process.DYScan4*process.DYScan4Asso*
process.DYScan5*process.DYScan5Asso*
process.DYScan10*process.DYScan10Asso*
process.DYScan20*process.DYScan20Asso*
process.DotDirScan0p99*process.DotDirScan0p99Asso*
process.DotDirScan0p95*process.DotDirScan0p95Asso*
process.DotDirScan0p90*process.DotDirScan0p90Asso*
process.DotDirScan0p85*process.DotDirScan0p85Asso*
process.GEMMuonAnalyzer
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
#process.p = cms.Path(process.trackerGEM)
#process.out_step = cms.EndPath(process.output)
#process.schedule = cms.Schedule(
#    process.p,
#    process.out_step
#)
