import FWCore.ParameterSet.Config as cms

#----------GEMMuon Collection Production for association by hits
gemmuon = cms.EDProducer("GEMMuonTrackCollProducer",
  MuonObj = cms.string("GEMMuon"),
  MaxPullX = cms.double(4.),
  MaxDX = cms.double(3.),
  MaxPullY = cms.double(9999.),
  MaxDY = cms.double(9999.),
  MinDotDir = cms.double(-9999.),
                         )
#--------------------
gemmuonColl_seq = cms.Sequence(
                             gemmuon
                             )
