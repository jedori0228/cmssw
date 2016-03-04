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
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class GEMMuonTrackCollProducer : public edm::EDProducer {
public:
  explicit GEMMuonTrackCollProducer(const edm::ParameterSet&);
  //std::vector<double> findSimVtx(edm::Event& iEvent);
  ~GEMMuonTrackCollProducer();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
};


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GEMMuonTrackCollProducer);


GEMMuonTrackCollProducer::GEMMuonTrackCollProducer(const edm::ParameterSet& iConfig){
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
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
  
  std::auto_ptr<reco::TrackCollection> selectedTracks(new reco::TrackCollection);
 
  reco::TrackRefProd rTracks = iEvent.getRefBeforePut<reco::TrackCollection>();
  


  for(reco::MuonCollection::const_iterator recomuon=recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon) {
    if (!recomuon->isGEMMuon()) continue;
    reco::TrackRef trackref;    

    if (recomuon->innerTrack().isNonnull()) trackref = recomuon->innerTrack();

      const reco::Track* trk = &(*trackref);
      // pointer to old track:
      //reco::Track* newTrk = new reco::Track(*trk);

      selectedTracks->push_back( *trk );
      //selectedTrackExtras->push_back( *newExtra );
  }
  iEvent.put(selectedTracks);

}

