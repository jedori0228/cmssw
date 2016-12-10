#ifndef MuonIdentification_trackerGEM_h
#define MuonIdentification_trackerGEM_h

/** \class trackerGEM 
 *
 * \author Jason Lee, based off ME0Muon
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <DataFormats/MuonDetId/interface/GEMDetId.h>

#include "FWCore/ServiceRegistry/interface/Service.h"

#include <DataFormats/GEMRecHit/interface/GEMSegmentCollection.h>

#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"

class FreeTrajectoryState;
class MagneticField;
class SteppingHelixPropagator;
class trackerGEM : public edm::EDProducer {
 public:
  /// Constructor
  explicit trackerGEM(const edm::ParameterSet&);
  /// Destructor
  ~trackerGEM();
  /// Produce the GEMSegment collection
  virtual void produce(edm::Event&, const edm::EventSetup&);

    
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endJob();

  std::vector<reco::MuonChamberMatch> MatchGEM(const reco::Track&, const GEMSegmentCollection&, const SteppingHelixPropagator*);
    
  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
			     int , const AlgebraicSymMatrix66& ,
			     const MagneticField* );

  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
			     int , const AlgebraicSymMatrix55& ,
			     const MagneticField* );

  void getFromFTS(const FreeTrajectoryState& ,
		  GlobalVector& , GlobalVector& , 
		  int& , AlgebraicSymMatrix66& );

 private:


  edm::ESHandle<GEMGeometry> gemGeom;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegmentsToken_;
  edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;

};

#endif
