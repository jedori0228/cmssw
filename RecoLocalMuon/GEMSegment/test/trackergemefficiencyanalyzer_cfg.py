import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackerGEMEfficiency")
process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDevReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')

# CSCGeometry depends on alignment ==> necessary to provide GlobalPositionRecord
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")

# Automatic addition of the customisation function from Geometry.GEMGeometry.gemGeometryCustoms
from Geometry.GEMGeometry.gemGeometryCustoms import custom_GE11_8and8partitions_v2
#call to customisation function custom_GE11_8and8partitions_v2 imported from Geometry.GEMGeometry.gemGeometryCustoms
process = custom_GE11_8and8partitions_v2(process)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#import sys

#infile = sys.argv[2]
#outroot = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #open('Sh_trackerGEM_path.txt').readlines()
        #open('HGC_trackerGEM_path.txt').readlines()
        #'file:/cms/home/jskim/cmssw/CMSSW_6_2_0_SLHC27_trackerGEM_trackerMuon/src/work/out_reco.root'
        open('filelist_MuonGun_modify_TrackDetectorAssociator_newGEO_pdigi_valid.txt').readlines()
        #open('filelist_MuonGun_modify_TrackDetectorAssociator.txt').readlines()
        #'file:/xrootd/store/user/jskim/condor/MuonGun_20_jobs_500_events_trackerGEM_trackerMuon/out_reco_0.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True),
)

process.trackergem = cms.EDAnalyzer('TrackerGEMEfficiencyAnalyzer',
                              # ----------------------------------------------------------------------
                              #RootFileName = cms.untracked.string("Sh_efficiency_0p1_matching.root"),
                              #RootFileName = cms.untracked.string("HGC_efficiency_0p1_matching.root"),
                              RootFileName = cms.untracked.string('temp.root'),
                              #RootFileName = cms.untracked.string('eff_oldGEO.root'),
                              # ----------------------------------------------------------------------
                              printSegmntInfo = cms.untracked.bool(False),
                              # ----------------------------------------------------------------------

)

process.p = cms.Path(process.trackergem)
