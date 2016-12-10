import FWCore.ParameterSet.Config as cms

process = cms.Process("GEMMuonAnalyzer")

#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load( "DQMServices/Core/DQMStore_cfg" )
process.load('Validation.RecoMuon.associators_cff')
process.load('Validation.RecoMuon.selectors_cff')
process.load('Validation.RecoMuon.MuonTrackValidator_cfi')
process.load('Validation.RecoMuon.RecoMuonValidator_cfi')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/user/jskim/FromArchie/w_GE21_w_ME21/500um/out_STA_reco_GE21_ME21_pT100_1.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True), 
)

from RecoLocalMuon.GEMSegment.trackerGEM_cfi import *
process.trackerGEM = trackerGEM.clone(
  gemSegmentsToken = cms.InputTag("gemSegments", "", "STARECO"),
)

from CommonTools.RecoAlgos.gemAssociator import *
process.gemMuonSel = gemmuon.clone(
  MuonObj = cms.string("MatchingStudy"),
  trackTag = cms.InputTag("trackerGEM"),
)

###  Thsi is for association by hits
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
  HistoFile = cms.string('OUTPUT.root'),
  FakeRatePtCut = cms.double(5.0),
  MatchingWindowDelR = cms.double (.15),
  UseAssociators = cms.bool(True),
  doGeometryStudy = cms.bool(True),
  associators = cms.vstring('gemMuonAssociatorByHits'),
  label = cms.vstring('gemMuonSel'),
)


process.p = cms.Path(
process.trackerGEM
*process.gemMuonSel
*process.gemMuonAssociatorByHits
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
#process.p = cms.Path(process.trackerGEM)
#process.out_step = cms.EndPath(process.output)
#process.schedule = cms.Schedule(
#    process.p,
#    process.out_step
#)


# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCalMuon

#call to customisation function cust_2023HGCalMuon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCalMuon(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from Geometry.GEMGeometry.gemGeometryCustoms
from Geometry.GEMGeometry.gemGeometryCustoms import custom_GE11_8and8partitions_v2

#call to customisation function custom_GE11_8and8partitions_v2 imported from Geometry.GEMGeometry.gemGeometryCustoms
process = custom_GE11_8and8partitions_v2(process)

# End of customisation functions





