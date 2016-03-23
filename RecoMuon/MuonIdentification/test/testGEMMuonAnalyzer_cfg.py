import FWCore.ParameterSet.Config as cms

process = cms.Process("GEMMuonAnalyzer")

process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')


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
      #'file:/cms/home/jskim/cmssw/CMSSW_6_2_0_SLHC27_trackerGEM_trackerMuon/src/work/assohits/out_sim.root'
      #'file:/cms/home/jskim/cmssw/CMSSW_6_2_0_SLHC27_trackerGEM_trackerMuon/src/work/out_reco_newGEO_pdigi_valid.root'
      #open('filelist_MuonGun_modify_TrackDetectorAssociator_newGEO_pdigi_valid.txt').readlines()
      #'/store/user/jskim/MuonGun_100_jobs_100_events_modify_TrackDetectorAssociator_geofixed_newGEO/out_reco_newGEO_pdigi_valid_000.root'
      #open('filelist_MuonGun_modify_TrackDetectorAssociator.txt').readlines()
      open('filelist_MuonGun_modify_TrackDetectorAssociator_geofixed_oldGEO.txt').readlines()
      #open('filelist_MinBias_modify_TrackDetectorAssociator_geofixed_oldGEO.txt').readlines()
    ), ##/somewhere/simevent.root" }
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True), 

)


from Validation.RecoMuon.selectors_cff import *
from Validation.RecoMuon.associators_cff import *
# Configurations for MuonTrackValidators
import Validation.RecoMuon.MuonTrackValidator_cfi


# Configurations for RecoMuonValidators
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from Validation.RecoMuon.RecoMuonValidator_cfi import *

#import SimGeneral.MixingModule.mixNoPU_cfi
from SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi import *
from SimMuon.MCTruth.MuonAssociatorByHits_cfi import muonAssociatorByHitsCommonParameters

#process.TrackAssociatorByChi2ESProducer = Validation.RecoMuon.associators_cff.TrackAssociatorByChi2ESProducer.clone(chi2cut = 6.0,ComponentName = 'TrackAssociatorByChi2')

###  Thsi is for association by hits
import SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi

process.muonAssociatorByHits = SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi.muonAssociatorByHitsESProducer.clone(ComponentName = 'muonAssociatorByHits',
 #tpTag = 'mix:MergedTrackTruth',
 UseTracker = True, 
 UseMuon = False,
 EfficiencyCut_track = cms.double(0.0),   
 PurityCut_track = cms.double(0.0)       
 )


process.GEMMuonAnalyzer = cms.EDAnalyzer("GEMMuonAnalyzer",
                              HistoFile = cms.string('GEMMuonAnalyzerOutput_MuonGun_old.root'),
                              FakeRatePtCut = cms.double(5.0),
                              MatchingWindowDelR = cms.double (.15),
                              UseAssociators = cms.bool(True),
                              associators = cms.vstring('muonAssociatorByHits'),
                              label = cms.VInputTag('gemMuonSel'),
                              #label = cms.VInputTag(cms.InputTag("standAloneMuons", "")),
                              #associatormap = cms.InputTag("tpToMuonTrackAssociation"),

                              # selection of GP for evaluation of efficiency
                              ptMinGP = cms.double(0.9),
                              minRapidityGP = cms.double(-2.5),
                              maxRapidityGP = cms.double(2.5),
                              tipGP = cms.double(3.5),
                              lipGP = cms.double(30.0),
                              chargedOnlyGP = cms.bool(True),
                              statusGP = cms.int32(1),
                              pdgIdGP = cms.vint32(13, -13),
                              #parametersDefiner = cms.string('LhcParametersDefinerForTP')

)

#process.p = cms.Path(process.GEMMuonAnalyzer)

from CommonTools.RecoAlgos.gemAssociator import *

process.gemMuonSel = gemmuon

#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.output = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string(
#        'file:out_reco.root'
#    ),
#)
#process.out_step = cms.EndPath(process.output)

process.p = cms.Path(process.gemMuonSel*process.GEMMuonAnalyzer)

#process.p = cms.Path(process.gemMuonSel)
#process.Schedule = cms.Schedule( process.p, process.out_step)
