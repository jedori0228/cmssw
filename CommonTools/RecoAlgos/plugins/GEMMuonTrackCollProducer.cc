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
  edm::EDGetTokenT<reco::VertexCollection> vertex_Token;
  edm::EDGetTokenT <reco::TrackCollection > generalTracksToken_;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  std::string MuonObj;
  double MaxPullX, MaxX, MaxPullY, MaxY;
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
  MaxPullX = parset.getParameter< double >("MaxPullX");
  MaxX = parset.getParameter< double >("MaxX");
  MaxPullY = parset.getParameter< double >("MaxPullY");
  MaxY = parset.getParameter< double >("MaxY");
  produces<reco::TrackCollection>();
}

GEMMuonTrackCollProducer::~GEMMuonTrackCollProducer() {
}

void GEMMuonTrackCollProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
  using namespace edm;
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertex_Token, vertices);
  reco::Vertex vertex = vertices->at(0);
  
  std::auto_ptr<reco::TrackCollection> selectedTracks(new reco::TrackCollection);
 
  reco::TrackRefProd rTracks = iEvent.getRefBeforePut<reco::TrackCollection>();
 
  if(MuonObj=="MatchingStudy"){

    edm::Handle<GEMSegmentCollection> gemSegmentCollection;
    iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);

    Handle<TrackCollection> generalTracks;
    iEvent.getByToken(generalTracksToken_, generalTracks);

    edm::ESHandle<GEMGeometry> gemGeom;
    iSetup.get<MuonGeometryRecord>().get(gemGeom);
    ESHandle<MagneticField> bField;
    iSetup.get<IdealMagneticFieldRecord>().get(bField);
    const SteppingHelixPropagator* shPropagator;
    shPropagator = new SteppingHelixPropagator(&*bField,alongMomentum);

    for(std::vector<Track>::const_iterator thisTrack = generalTracks->begin(); thisTrack != generalTracks->end(); ++thisTrack){

      if(thisTrack->pt() < 0.5) continue;
      if(thisTrack->p() < 2.5) continue;
      //if(std::abs(thisTrack->eta()) < 1.5) continue;
      //if(std::abs(thisTrack->eta()) > 3.0) continue;

      bool isGEMMuon = false;

      for(auto thisSegment = gemSegmentCollection->begin(); thisSegment != gemSegmentCollection->end(); ++thisSegment){

        GEMDetId id = thisSegment->specificRecHits()[0].gemId();
        int station = id.station();
        if (id.station() != station) continue;
        float zSign = thisTrack->pz() > 0 ? 1.0f : -1.0f;
        if ( zSign * id.region() < 0 ) continue;

        LocalPoint thisPosition(thisSegment->localPosition());
        LocalVector thisDirection(thisSegment->localDirection());

        //auto chamber = gemGeom->superChamber(id);
        auto chamber = gemGeom->chamber(id);
        GlobalPoint SegPos(chamber->toGlobal(thisPosition));
        GlobalVector SegDir(chamber->toGlobal(thisDirection));

        const float zValue = SegPos.z();

        Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());

        //Getting the initial variables for propagation

        int chargeReco = thisTrack->charge();
        GlobalVector p3reco, r3reco;

        p3reco = GlobalVector(thisTrack->outerPx(), thisTrack->outerPy(), thisTrack->outerPz());
        r3reco = GlobalVector(thisTrack->outerX(), thisTrack->outerY(), thisTrack->outerZ());

        AlgebraicSymMatrix66 covReco;
        //This is to fill the cov matrix correctly
        AlgebraicSymMatrix55 covReco_curv;
        covReco_curv = thisTrack->outerStateCovariance();
        FreeTrajectoryState initrecostate = getFTS(p3reco, r3reco, chargeReco, covReco_curv, shPropagator->magneticField());
        getFromFTS(initrecostate, p3reco, r3reco, chargeReco, covReco);

        //Now we propagate and get the propagated variables from the propagated state
        SteppingHelixStateInfo startrecostate(initrecostate);
        SteppingHelixStateInfo lastrecostate;

        //const SteppingHelixPropagator* shPropagator = 
        //dynamic_cast<const SteppingHelixPropagator*>(&*shProp);
        // for 62XSLHC
        //lastrecostate = shPropagator->propagate(startrecostate, *plane);
        //lastrecostate = shPropagator->propagateWithPath(startrecostate, *plane);
        // for 76X
        shPropagator->propagate(startrecostate, *plane,lastrecostate);

        FreeTrajectoryState finalrecostate;
        lastrecostate.getFreeState(finalrecostate);

        AlgebraicSymMatrix66 covFinalReco;
        GlobalVector p3FinalReco_glob, r3FinalReco_globv;
        getFromFTS(finalrecostate, p3FinalReco_glob, r3FinalReco_globv, chargeReco, covFinalReco);

        //To transform the global propagated track to local coordinates
        GlobalPoint r3FinalReco_glob(r3FinalReco_globv.x(),r3FinalReco_globv.y(),r3FinalReco_globv.z());

        LocalPoint r3FinalReco = chamber->toLocal(r3FinalReco_glob);
        LocalVector p3FinalReco=chamber->toLocal(p3FinalReco_glob);

        //The same goes for the error
        AlgebraicMatrix thisCov(4,4,0);
        for (int i = 1; i <=4; i++){
          for (int j = 1; j <=4; j++){
            thisCov(i,j) = thisSegment->parametersError()(i,j);
          }
        }

        LocalTrajectoryParameters ltp(r3FinalReco,p3FinalReco,chargeReco);
        JacobianCartesianToLocal jctl(chamber->surface(),ltp);
        AlgebraicMatrix56 jacobGlbToLoc = jctl.jacobian();

        AlgebraicMatrix55 Ctmp =  (jacobGlbToLoc * covFinalReco) * ROOT::Math::Transpose(jacobGlbToLoc);
        AlgebraicSymMatrix55 C;  // I couldn't find any other way, so I resort to the brute force
        for(int i=0; i<5; ++i) {
          for(int j=0; j<5; ++j) {
            C[i][j] = Ctmp[i][j];
          }
        }

        Double_t sigmax = sqrt(C[3][3]+thisSegment->localPositionError().xx() );
        Double_t sigmay = sqrt(C[4][4]+thisSegment->localPositionError().yy() );

        Double_t DelX = std::abs(thisPosition.x()-r3FinalReco.x());
        Double_t DelX_over_sigma = DelX/sigmax;
        Double_t DelY = std::abs(thisPosition.y()-r3FinalReco.y());
        Double_t DelY_over_sigma = DelY/sigmay;

        //std::cout << "DelX = " << DelX << ", MaxX = " << MaxX << std::endl;
        //std::cout << "PullX = " << DelX_over_sigma << ", MaxPullX = " << MaxPullX << std::endl;
        //std::cout << "DelY = " << DelY << ", MaxY = " << MaxY << std::endl;
        //std::cout << "PullY = " << DelY_over_sigma << ", MaxPullY = " << MaxPullY << std::endl;

        bool XMatched = (DelX < MaxX) || (DelX_over_sigma < MaxPullX);
        bool YMatched = (DelY < MaxY) || (DelY_over_sigma < MaxPullY);

        if(XMatched && YMatched){
          isGEMMuon = true;
          break;
        }


      } // GEMSegment loop

      if(isGEMMuon) selectedTracks->push_back( *thisTrack );

    } // General Track Loop

  }

  else{

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
