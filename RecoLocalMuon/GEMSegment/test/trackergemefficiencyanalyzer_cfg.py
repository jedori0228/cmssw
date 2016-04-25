import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackerGEMEfficiency")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDevReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')

# CSCGeometry depends on alignment ==> necessary to provide GlobalPositionRecord
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")

# Automatic addition of the customisation function from Geometry.GEMGeometry.gemGeometryCustoms
#from Geometry.GEMGeometry.gemGeometryCustoms import custom_GE11_8and8partitions_v2
#call to customisation function custom_GE11_8and8partitions_v2 imported from Geometry.GEMGeometry.gemGeometryCustoms
#process = custom_GE11_8and8partitions_v2(process)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#import sys

#infile = sys.argv[2]
#outroot = sys.argv[3]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #open('filelist_MuonGun.txt').readlines()
        'file:../../../work/out_recodone.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True),
)

process.trackergem = cms.EDAnalyzer('TrackerGEMEfficiencyAnalyzer',
                              # ----------------------------------------------------------------------
                              #RootFileName = cms.untracked.string('eff.root'),
                              RootFileName = cms.untracked.string('eff_test.root'),
                              # ----------------------------------------------------------------------
                              printSegmntInfo = cms.untracked.bool(False),
                              # ----------------------------------------------------------------------

)

process.p = cms.Path(process.trackergem)

# Automatic addition of the customisation function from Geometry.GEMGeometry.gemGeometryCustoms
#from Geometry.GEMGeometry.gemGeometryCustoms import custom_GE11_8and8partitions_v2

#call to customisation function custom_GE11_8and8partitions_v2 imported from Geometry.GEMGeometry.gemGeometryCustoms
#process = custom_GE11_8and8partitions_v2(process)
