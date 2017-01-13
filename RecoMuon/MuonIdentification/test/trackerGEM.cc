#include <RecoMuon/MuonIdentification/test/trackerGEM.h>

#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

#include <DataFormats/MuonReco/interface/Muon.h>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"


#include "DataFormats/GeometrySurface/interface/LocalError.h"


#include "TLorentzVector.h"

#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"


#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include <DataFormats/GeometrySurface/interface/SimpleDiskBounds.h>

trackerGEM::trackerGEM(const edm::ParameterSet& iConfig) {
  gemSegmentsToken_ = consumes<GEMSegmentCollection >(iConfig.getParameter<edm::InputTag>("gemSegmentsToken"));
  generalTracksToken_ = consumes<reco::TrackCollection >(iConfig.getParameter<edm::InputTag>("generalTracksToken"));

  produces<std::vector<reco::Muon> >();

}

trackerGEM::~trackerGEM() {
}

void trackerGEM::produce(edm::Event& ev, const edm::EventSetup& setup) {
  using namespace edm;
  using namespace reco;
  
  ESHandle<MagneticField> bField;
  setup.get<IdealMagneticFieldRecord>().get(bField);
  const SteppingHelixPropagator* ThisshProp;
  ThisshProp = new SteppingHelixPropagator(&*bField,alongMomentum);

  Handle<GEMSegmentCollection> gemSegments;
  ev.getByToken(gemSegmentsToken_,gemSegments);

  Handle<TrackCollection > generalTracks;
  ev.getByToken(generalTracksToken_,generalTracks);

  std::unique_ptr<std::vector<Muon> > muons( new std::vector<Muon> ); 

  int TrackNumber = 0;
  for(std::vector<Track>::const_iterator thisTrack = generalTracks->begin(); thisTrack != generalTracks->end(); ++thisTrack, ++TrackNumber){
    //Initializing gem plane
    //Remove later
    if (thisTrack->pt() < 1.5) continue;
    if (std::fabs(thisTrack->eta()) < 1.5) continue;

    std::vector<reco::MuonChamberMatch> muonChamberMatches = MatchGEM(*thisTrack, *gemSegments, ThisshProp);

    TrackRef thisTrackRef(generalTracks,TrackNumber);
    	   
    // temp settting the muon to track p4
    Particle::Charge q = thisTrackRef->charge();
    Particle::LorentzVector p4(thisTrackRef->px(), thisTrackRef->py(), thisTrackRef->pz(), thisTrackRef->p());
    Particle::Point vtx(thisTrackRef->vx(),thisTrackRef->vy(), thisTrackRef->vz());

    reco::Muon MuonCandidate = reco::Muon(q, p4, vtx);

    MuonCandidate.setTrack(thisTrackRef);
    // need to make track from gem seg
    MuonCandidate.setOuterTrack(thisTrackRef);
    //MuonCandidate.setType(thisSegment->nRecHits());
    MuonCandidate.setMatches(muonChamberMatches);

    //MuonCandidate.setGlobalTrackPosAtSurface(r3FinalReco_glob);
    //MuonCandidate.setGlobalTrackMomAtSurface(p3FinalReco_glob);
    //MuonCandidate.setLocalTrackPosAtSurface(r3FinalReco);
    //MuonCandidate.setLocalTrackMomAtSurface(p3FinalReco);
    //MuonCandidate.setGlobalTrackCov(covFinalReco);
    //MuonCandidate.setLocalTrackCov(C);
    
    muons->push_back(MuonCandidate);
  }
  
  // put collection in event

  ev.put(std::move(muons));
  delete ThisshProp;
}

FreeTrajectoryState
trackerGEM::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
		   int charge, const AlgebraicSymMatrix55& cov,
		   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);
  
  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
trackerGEM::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
		   int charge, const AlgebraicSymMatrix66& cov,
		   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void trackerGEM::getFromFTS(const FreeTrajectoryState& fts,
			    GlobalVector& p3, GlobalVector& r3, 
			    int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  GlobalVector p3T(p3GV.x(), p3GV.y(), p3GV.z());
  GlobalVector r3T(r3GP.x(), r3GP.y(), r3GP.z());
  p3 = p3T;
  r3 = r3T;  
  // p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  // r3.set(r3GP.x(), r3GP.y(), r3GP.z());
  
  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}

void trackerGEM::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(gemGeom);

}
void trackerGEM::endJob()
{
}

std::vector<reco::MuonChamberMatch> trackerGEM::MatchGEM(const reco::Track& track, const GEMSegmentCollection& gemSegments, const SteppingHelixPropagator* shPropagator)
{
  std::vector<reco::MuonChamberMatch> vchamber;
  
  for(auto thisSegment = gemSegments.begin(); thisSegment != gemSegments.end(); ++thisSegment){

    GEMDetId id = thisSegment->gemDetId();
    float zSign = track.pz() > 0 ? 1.0f : -1.0f;
    if ( zSign * id.region() < 0 ) continue;

    LocalPoint thisPosition(thisSegment->localPosition());
    LocalVector thisDirection(thisSegment->localDirection());

    auto chamber = gemGeom->superChamber(id);
    GlobalPoint SegPos(chamber->toGlobal(thisPosition));
    GlobalVector SegDir(chamber->toGlobal(thisDirection));

    const float zValue = SegPos.z();

    Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());

    //Getting the initial variables for propagation

    int chargeReco = track.charge(); 
    GlobalVector p3reco, r3reco;

    p3reco = GlobalVector(track.outerPx(), track.outerPy(), track.outerPz());
    r3reco = GlobalVector(track.outerX(), track.outerY(), track.outerZ());

    AlgebraicSymMatrix66 covReco;
    //This is to fill the cov matrix correctly
    AlgebraicSymMatrix55 covReco_curv;
    covReco_curv = track.outerStateCovariance();
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
    /////////////////////////////////////////////////////////////////////////////////////////

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

    Double_t sigmax = C[3][3]+thisSegment->localPositionError().xx();
    Double_t sigmay = C[4][4]+thisSegment->localPositionError().yy();

    double DelX   = std::abs(thisPosition.x()-r3FinalReco.x());
    double DelY   = std::abs(thisPosition.y()-r3FinalReco.y());
    //double DelPhi = std::abs( p3FinalReco.phi() - thisDirection.phi() );
    double DotDir = p3FinalReco.unit().dot( thisDirection.unit() );

    reco::MuonChamberMatch matchedChamber;;
    matchedChamber.id = thisSegment->specificRecHits()[0].gemId();
    //==== save DX and DY
    matchedChamber.x = DelX;
    matchedChamber.y = DelY;
    //==== save PullX and PullY
    matchedChamber.xErr = DelX/sigmax;
    matchedChamber.yErr = DelY/sigmay;
    //==== save DelPhi
    matchedChamber.dXdZ = DotDir;

    //std::cout << "=============================================================================================================" << std::endl;
    //std::cout << "Track Position : eta = "<< r3FinalReco_glob.eta() << ", phi = " << r3FinalReco_glob.phi() << std::endl;
    //std::cout << "Segment Position : eta = "<< SegPos.eta() << ", phi = " << SegPos.phi() << std::endl;
    //std::cout << "Track Direction : x = "<< p3FinalReco.unit().x() << ", y = " << p3FinalReco.unit().y() << ", z = " << p3FinalReco.unit().z() << std::endl;
    //std::cout << "Segment Direction : x = "<< thisDirection.x() << ", y = " << thisDirection.y() << ", z = " << thisDirection.z() << std::endl;
    //std::cout << "=> DotDir = " << DotDir << std::endl;

    vchamber.push_back(matchedChamber);

  }

  return vchamber;

}


DEFINE_FWK_MODULE(trackerGEM);
