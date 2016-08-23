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

class GEMMuonTrackCollProducer : public edm::stream::EDProducer<> {
public:
  explicit GEMMuonTrackCollProducer(const edm::ParameterSet&);
  //std::vector<double> findSimVtx(edm::Event& iEvent);
  ~GEMMuonTrackCollProducer();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<reco::VertexCollection> vertex_Token;
  std::string MuonObj;
};


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GEMMuonTrackCollProducer);


GEMMuonTrackCollProducer::GEMMuonTrackCollProducer(const edm::ParameterSet& parset){
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  vertex_Token = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  MuonObj = parset.getParameter< std::string >("MuonObj");
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
  iEvent.put(selectedTracks);

}
