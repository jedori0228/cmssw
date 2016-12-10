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
#include <DataFormats/MuonReco/interface/MuonSelectors.h>

class GEMMuonTrackCollProducer : public edm::EDProducer {
public:
  explicit GEMMuonTrackCollProducer(const edm::ParameterSet&);
  //std::vector<double> findSimVtx(edm::Event& iEvent);
  ~GEMMuonTrackCollProducer();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<reco::MuonCollection> TrackerGEM_Token;
  edm::EDGetTokenT<reco::VertexCollection> vertex_Token;
  std::string MuonObj;
  edm::InputTag trackTag;
  double MaxPullX, MaxDX, MaxPullY, MaxDY, MinDotDir;
};


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GEMMuonTrackCollProducer);


GEMMuonTrackCollProducer::GEMMuonTrackCollProducer(const edm::ParameterSet& iConfig){
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  vertex_Token = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  MuonObj = iConfig.getParameter< std::string >("MuonObj");
  trackTag = iConfig.getParameter< edm::InputTag >("trackTag");
  if(MuonObj == "MatchingStudy") TrackerGEM_Token = consumes<reco::MuonCollection>(trackTag);
  MaxPullX = iConfig.getParameter< double >("MaxPullX");
  MaxDX = iConfig.getParameter< double >("MaxDX");
  MaxPullY = iConfig.getParameter< double >("MaxPullY");
  MaxDY = iConfig.getParameter< double >("MaxDY");
  MinDotDir = iConfig.getParameter< double >("MinDotDir");
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


