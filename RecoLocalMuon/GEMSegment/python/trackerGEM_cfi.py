import FWCore.ParameterSet.Config as cms


trackerGEM = cms.EDProducer("trackerGEM",
    gemSegmentsToken = cms.InputTag("gemSegments"),
    generalTracksToken = cms.InputTag("generalTracks"),
)
