import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackerGEMEfficiency")
process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDevReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')

# CSCGeometry depends on alignment ==> necessary to provide GlobalPositionRecord
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#import sys

#infile = sys.argv[2]
#outroot = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #open('Sh_trackerGEM_path.txt').readlines()
        open('HGC_trackerGEM_path.txt').readlines()
        #'file:/xrootd/store/user/jskim/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/trackergem_FullScope_no_aging/step3_1000_1_XmT.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True),
)

process.trackergem = cms.EDAnalyzer('TrackerGEMEfficiencyAnalyzer',
                              # ----------------------------------------------------------------------
                              #RootFileName = cms.untracked.string("Sh_efficiency_0p1_matching.root"),
                              RootFileName = cms.untracked.string("HGC_efficiency_0p1_matching.root"),
                              #RootFileName = cms.untracked.string('temp.root'),
                              # ----------------------------------------------------------------------
                              printSegmntInfo = cms.untracked.bool(False),
                              # ----------------------------------------------------------------------

)

process.p = cms.Path(process.trackergem)
