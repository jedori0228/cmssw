import FWCore.ParameterSet.Config as cms

#----------GEMMuon Collection Production for association by hits
gemmuon = cms.EDProducer("GEMMuonTrackCollProducer",
                         )
#--------------------
gemmuonColl_seq = cms.Sequence(
                             gemmuon
                             )
