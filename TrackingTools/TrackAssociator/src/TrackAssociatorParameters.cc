// -*- C++ -*-
//
// Package:    TrackAssociator
// Class:      TrackAssociatorParameters
// 
/*

 Description: track associator parameters

*/
//
// Original Author:  Dmytro Kovalskyi
// $Id: TrackAssociatorParameters.cc,v 1.7 2009/09/06 16:32:06 dmytro Exp $
//
//

#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"


void TrackAssociatorParameters::loadParameters( const edm::ParameterSet& iConfig )
{
   dREcal = iConfig.getParameter<double>("dREcal");
   dRHcal = iConfig.getParameter<double>("dRHcal");
   dRMuon = iConfig.getParameter<double>("dRMuon");
   
   dREcalPreselection = iConfig.getParameter<double>("dREcalPreselection");
   dRHcalPreselection = iConfig.getParameter<double>("dRHcalPreselection");
   dRMuonPreselection = iConfig.getParameter<double>("dRMuonPreselection");
   dRPreshowerPreselection = iConfig.getParameter<double>("dRPreshowerPreselection");
   
   muonMaxDistanceX = iConfig.getParameter<double>("muonMaxDistanceX");
   muonMaxDistanceY = iConfig.getParameter<double>("muonMaxDistanceY");
   muonMaxDistanceSigmaX = iConfig.getParameter<double>("muonMaxDistanceSigmaX");
   muonMaxDistanceSigmaY = iConfig.getParameter<double>("muonMaxDistanceSigmaY");
   
   useEcal = iConfig.getParameter<bool>("useEcal");
   useHcal = iConfig.getParameter<bool>("useHcal");
   useHO   = iConfig.getParameter<bool>("useHO");
   useCalo = iConfig.getParameter<bool>("useCalo");
   useMuon = iConfig.getParameter<bool>("useMuon");
   usePreshower = iConfig.getParameter<bool>("usePreshower");
   
   theEBRecHitCollectionLabel       = iConfig.getParameter<edm::InputTag>("EBRecHitCollectionLabel");
   theEERecHitCollectionLabel       = iConfig.getParameter<edm::InputTag>("EERecHitCollectionLabel");
   theCaloTowerCollectionLabel      = iConfig.getParameter<edm::InputTag>("CaloTowerCollectionLabel");
   theHBHERecHitCollectionLabel     = iConfig.getParameter<edm::InputTag>("HBHERecHitCollectionLabel");
   theHORecHitCollectionLabel       = iConfig.getParameter<edm::InputTag>("HORecHitCollectionLabel");
   theDTRecSegment4DCollectionLabel = iConfig.getParameter<edm::InputTag>("DTRecSegment4DCollectionLabel");
   theCSCSegmentCollectionLabel     = iConfig.getParameter<edm::InputTag>("CSCSegmentCollectionLabel");
   theGEMSegmentCollectionLabel     = iConfig.getParameter<edm::InputTag>("GEMSegmentCollectionLabel");
   
   accountForTrajectoryChangeCalo   = iConfig.getParameter<bool>("accountForTrajectoryChangeCalo");
   // accountForTrajectoryChangeMuon   = iConfig.getParameter<bool>("accountForTrajectoryChangeMuon");
   
   truthMatch = iConfig.getParameter<bool>("truthMatch");
   muonMaxDistanceSigmaY = iConfig.getParameter<double>("trajectoryUncertaintyTolerance");


   maxPullXGE11_   = iConfig.getParameter<double>("maxPullXGE11");
   maxDiffXGE11_   = iConfig.getParameter<double>("maxDiffXGE11");
   maxPullYGE11_   = iConfig.getParameter<double>("maxPullYGE11");
   maxDiffYGE11_   = iConfig.getParameter<double>("maxDiffYGE11");
   maxPullXGE21_   = iConfig.getParameter<double>("maxPullXGE21");
   maxDiffXGE21_   = iConfig.getParameter<double>("maxDiffXGE21");
   maxPullYGE21_   = iConfig.getParameter<double>("maxPullYGE21");
   maxDiffYGE21_   = iConfig.getParameter<double>("maxDiffYGE21");
   maxDiffPhiDirection_ = iConfig.getParameter<double>("maxDiffPhiDirection");

}

TrackAssociatorParameters::TrackAssociatorParameters( const edm::ParameterSet& iConfig )
{
   loadParameters( iConfig );
}

