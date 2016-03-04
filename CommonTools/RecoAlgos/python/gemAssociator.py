import FWCore.ParameterSet.Config as cms

#----------ME0Muon Collection Production for association by chi2
gemmuon = cms.EDProducer("GEMMuonTrackCollProducer",
                         #me0MuonTag = cms.InputTag("me0SegmentMatching"),
                         #selectionTags = cms.vstring('All'),
                         )
#--------------------
gemmuonColl_seq = cms.Sequence(
                             gemmuon
                             )
