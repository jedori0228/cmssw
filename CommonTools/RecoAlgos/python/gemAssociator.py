import FWCore.ParameterSet.Config as cms

#----------GEMMuon Collection Production for association by hits
gemmuon = cms.EDProducer("GEMMuonTrackCollProducer",
  MuonObj = cms.string("GEMMuon"),
  MaxPullX = cms.double(3.),
  MaxX = cms.double(4.),
  MaxPullY = cms.double(9999.),
  MaxY = cms.double(9999.),
                         )
#--------------------
gemmuonColl_seq = cms.Sequence(
                             gemmuon
                             )
