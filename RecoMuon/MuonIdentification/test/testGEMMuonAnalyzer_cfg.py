import FWCore.ParameterSet.Config as cms

process = cms.Process("GEMMuonAnalyzer")


process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDevReco_cff')
process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDev_cff')


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



process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:///somewhere/simevent.root') ##/somewhere/simevent.root" }
    fileNames = cms.untracked.vstring(
      '/store/user/jskim/SingleMuPt100_cfi_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO/out_recodone_009.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True), 
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

import SimMuon.MCTruth.MuonAssociatorByHits_cfi
process.muonAssociatorByHits = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
 UseTracker = True,
 UseMuon = False,
 EfficiencyCut_track = cms.double(0.0),
 PurityCut_track = cms.double(0.0),
)

from CommonTools.RecoAlgos.gemAssociator import *

process.gemMuonSel = gemmuon

process.GEMMuonAnalyzer = cms.EDAnalyzer("GEMMuonAnalyzer",
                              

                              HistoFolder = cms.string('OUTPUTTEMPLATE'),
                              HistoFile = cms.string('OUTPUTTEMPLATE.root'),

                              FakeRatePtCut = cms.double(5.0),
                              MatchingWindowDelR = cms.double (.15),
                              UseGEMEtaCoverageMuons = cms.bool(True),
                              UseAssociators = cms.bool(True),

                              associators = cms.vstring('muonAssociatorByHits'),

                              label = cms.VInputTag('gemMuonSel'),
                              
                              
)

#process.p = cms.Path(process.gemMuonSel*process.muonAssociatorByHits*process.GEMMuonAnalyzer)
#process.p = cms.Path(process.gemMuonSel*process.GEMMuonAnalyzer)






process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(
        'file:out_test.root'
    ),
    outputCommands = cms.untracked.vstring(
        'keep  *_*_*_*',
    ),
)
process.out_step     = cms.EndPath(process.output)

process.p = cms.Path(process.gemMuonSel*process.muonAssociatorByHits*process.GEMMuonAnalyzer)
#process.p = cms.Path(process.gemMuonSel*process.muonAssociatorByHits)
process.out_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.p,
    process.out_step
)
