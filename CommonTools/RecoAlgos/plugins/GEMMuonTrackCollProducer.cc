#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "FWCore/Framework/interface/ESHandle.h"

#include <sstream>

#include <memory>
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/MuonReco/interface/MuonSelectors.h>

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>

class GEMMuonTrackCollProducer : public edm::stream::EDProducer<> {
public:
  explicit GEMMuonTrackCollProducer(const edm::ParameterSet&);
  //std::vector<double> findSimVtx(edm::Event& iEvent);
  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& ,
                             int , const AlgebraicSymMatrix66& ,
                             const MagneticField* );

  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& ,
                             int , const AlgebraicSymMatrix55& ,
                             const MagneticField* );
  void getFromFTS(const FreeTrajectoryState& ,
                  GlobalVector& , GlobalVector& ,
int& , AlgebraicSymMatrix66& );
  ~GEMMuonTrackCollProducer();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<reco::MuonCollection> TrackerGEM_Token;
  edm::EDGetTokenT<reco::VertexCollection> vertex_Token;
  edm::EDGetTokenT <reco::TrackCollection > generalTracksToken_;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  std::string MuonObj, trackTag;
  double MaxPullX, MaxDX, MaxPullY, MaxDY, MinDotDir;
};


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GEMMuonTrackCollProducer);


GEMMuonTrackCollProducer::GEMMuonTrackCollProducer(const edm::ParameterSet& parset){
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  vertex_Token = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  generalTracksToken_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
  GEMSegment_Token = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));
  MuonObj = parset.getParameter< std::string >("MuonObj");
  trackTag = parset.getParameter< std::string >("trackTag");
  if(MuonObj == "MatchingStudy") TrackerGEM_Token = consumes<reco::MuonCollection>(edm::InputTag(trackTag));
  MaxPullX = parset.getParameter< double >("MaxPullX");
  MaxDX = parset.getParameter< double >("MaxDX");
  MaxPullY = parset.getParameter< double >("MaxPullY");
  MaxDY = parset.getParameter< double >("MaxDY");
  MinDotDir = parset.getParameter< double >("MinDotDir");
  produces<reco::TrackCollection>();
}

GEMMuonTrackCollProducer::~GEMMuonTrackCollProducer() {
}

void GEMMuonTrackCollProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
  using namespace edm;
  
  std::auto_ptr<reco::TrackCollection> selectedTracks(new reco::TrackCollection);
 
  if(MuonObj=="MatchingStudy"){

    edm::Handle<reco::MuonCollection> TrackerGEM;
    iEvent.getByToken(TrackerGEM_Token, TrackerGEM);

    for(reco::MuonCollection::const_iterator recomuon=TrackerGEM->begin(); recomuon != TrackerGEM->end(); ++recomuon) {
      bool isGEMMuon = false;
      for(auto chmatch = recomuon->matches().begin(); chmatch != recomuon->matches().end(); ++chmatch){

        Double_t DelX = chmatch->x;
        Double_t DelX_over_sigma = chmatch->xErr;
        Double_t DelY = chmatch->y;
        Double_t DelY_over_sigma = chmatch->yErr;
        Double_t DotDir = chmatch->dXdZ;

        //std::cout << "DelX = " << DelX << ", MaxDX = " << MaxDX << std::endl;
        //std::cout << "PullX = " << DelX_over_sigma << ", MaxPullX = " << MaxPullX << std::endl;
        //std::cout << "DelY = " << DelY << ", MaxDY = " << MaxDY << std::endl;
        //std::cout << "PullY = " << DelY_over_sigma << ", MaxPullY = " << MaxPullY << std::endl;

        bool XMatched = (DelX < MaxDX) || (DelX_over_sigma < MaxPullX);
        bool YMatched = (DelY < MaxDY) || (DelY_over_sigma < MaxPullY);
        bool DirMatched = DotDir > MinDotDir;

        if(XMatched && YMatched && DirMatched){
          isGEMMuon = true;
          break;
        }
      } // END match loop

      if(isGEMMuon){
        if( !recomuon->innerTrack().isNonnull() ) continue;
        reco::TrackRef trackref = trackref = recomuon->innerTrack();
        const reco::Track* trk = &(*trackref);
        selectedTracks->push_back( *trk );        
      }

    } // TrackerGEM loop

  }
  else{

  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertex_Token, vertices);
  reco::Vertex vertex = vertices->at(0);

    for(reco::MuonCollection::const_iterator recomuon=recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon) {

      bool pass_obj(false);
      if(MuonObj=="RecoMuon") pass_obj = true;
      else if(MuonObj=="GEMMuon") pass_obj = recomuon->isGEMMuon() && recomuon->isTrackerMuon();
      else if(MuonObj=="LooseMuon") pass_obj = muon::isLooseMuon(*recomuon);
      else if(MuonObj=="MediumMuon") pass_obj = muon::isMediumMuon(*recomuon);
      else if(MuonObj=="TightMuon") pass_obj = muon::isTightMuon(*recomuon, vertex);
      else{}

      if(!pass_obj) continue;

      if( !recomuon->innerTrack().isNonnull() ) continue;
      reco::TrackRef trackref = trackref = recomuon->innerTrack();

      const reco::Track* trk = &(*trackref);
      // pointer to old track:
      //reco::Track* newTrk = new reco::Track(*trk);

      selectedTracks->push_back( *trk );
      //std::cout << "  track added" << std::endl;
      //selectedTrackExtras->push_back( *newExtra );
    }
  }
  iEvent.put(selectedTracks);

}

FreeTrajectoryState GEMMuonTrackCollProducer::getFTS(const GlobalVector& p3, const GlobalVector& r3, int charge, const AlgebraicSymMatrix55& cov, const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);

  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState GEMMuonTrackCollProducer::getFTS(const GlobalVector& p3, const GlobalVector& r3, int charge, const AlgebraicSymMatrix66& cov, const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);

  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void GEMMuonTrackCollProducer::getFromFTS(const FreeTrajectoryState& fts, GlobalVector& p3, GlobalVector& r3, int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  GlobalVector p3T(p3GV.x(), p3GV.y(), p3GV.z());
  GlobalVector r3T(r3GP.x(), r3GP.y(), r3GP.z());
  p3 = p3T;
  r3 = r3T;  //Yikes, was setting this to p3T instead of r3T!?!
  // p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  // r3.set(r3GP.x(), r3GP.y(), r3GP.z());

  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}
