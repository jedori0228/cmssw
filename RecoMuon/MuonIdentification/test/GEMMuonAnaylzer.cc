#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 
#include "FWCore/Utilities/interface/InputTag.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include "TLorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaR.h"

//Associator for chi2: Including header files
#include "SimTracker/TrackAssociatorProducers/plugins/TrackAssociatorByChi2Impl.h"
#include "SimTracker/TrackAssociatorProducers/plugins/TrackAssociatorByHitsImpl.h"

//#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
//#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>

#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/plugins/CosmicParametersDefinerForTPESProducer.h"

#include "CommonTools/CandAlgos/interface/GenParticleCustomSelector.h"

#include "Fit/FitResult.h"
#include "TF1.h" 


#include "TMath.h"
#include "TLorentzVector.h"

#include "TH1.h" 
#include <TH2.h>
#include "TFile.h"
#include <TProfile.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <TLatex.h>
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <DataFormats/MuonDetId/interface/GEMDetId.h>


#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"
#include "TGraph.h"

#include <sstream>    
#include <string>

#include <iostream>
#include <fstream>
#include <sys/stat.h>

class GEMMuonAnalyzer : public edm::EDAnalyzer {
public:
  explicit GEMMuonAnalyzer(const edm::ParameterSet&);
  ~GEMMuonAnalyzer();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void endRun(edm::Run const&, edm::EventSetup const&);

  //protected:
  
  private:

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
  edm::EDGetTokenT <reco::TrackCollection > generalTracksToken_;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > track_Collection_Token;
  std::vector<edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator>> associators_Token; 


  bool UseAssociators;
  bool UseGEMEtaCoverageMuons;
  const TrackAssociatorByChi2Impl* associatorByChi2;



  std::vector<std::string> associators;
  std::vector<edm::InputTag> label;
  std::vector<const reco::TrackToTrackingParticleAssociator*> associator;
  //std::vector<const reco::MuonToTrackingParticleAssociator*> associator;

  //Histos for plotting
  TString histoFolder;
  TFile* histoFile; 

  double  FakeRatePtCut, MatchingWindowDelR;

  double Nevents;

  TH1F* Nevents_h;

  TH1F *GenMuon_Eta; TH1F *GenMuon_Pt; TH1F* MatchedGEMMuon_Eta; TH1F* MatchedGEMMuon_Pt;
  TH1F *VertexDiff_h;
  TH2F *PtDiff_s; TProfile *PtDiff_p; TH1F *PtDiff_h; TH1F *QOverPtDiff_h;
  TH2F *PDiff_s; TProfile *PDiff_p; TH1F *PDiff_h;

};

GEMMuonAnalyzer::GEMMuonAnalyzer(const edm::ParameterSet& iConfig) 
{
  histoFile = new TFile(iConfig.getParameter<std::string>("HistoFile").c_str(), "recreate");
  histoFolder = iConfig.getParameter<std::string>("HistoFolder").c_str();
  UseGEMEtaCoverageMuons = iConfig.getParameter< bool >("UseGEMEtaCoverageMuons");
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");

  FakeRatePtCut   = iConfig.getParameter<double>("FakeRatePtCut");
  MatchingWindowDelR   = iConfig.getParameter<double>("MatchingWindowDelR");

  //Associator for chi2: getting parameters
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  associators = iConfig.getParameter< std::vector<std::string> >("associators");

  label = iConfig.getParameter< std::vector<edm::InputTag> >("label");
  edm::InputTag genParticlesTag ("genParticles");
  genParticlesToken_ = consumes<reco::GenParticleCollection>(genParticlesTag);
  edm::InputTag trackingParticlesTag ("mix","MergedTrackTruth");
  trackingParticlesToken_ = consumes<TrackingParticleCollection>(trackingParticlesTag);
  edm::InputTag generalTracksTag ("generalTracks");
  generalTracksToken_ = consumes<reco::TrackCollection>(generalTracksTag);
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  GEMSegment_Token = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));

  //Getting tokens and doing consumers for track associators
  for (unsigned int www=0;www<label.size();www++){
    track_Collection_Token.push_back(consumes<edm::View<reco::Track> >(label[www]));
  }

  if (UseAssociators) {
    for (auto const& thisassociator :associators) {
      associators_Token.push_back(consumes<reco::TrackToTrackingParticleAssociator>(edm::InputTag(thisassociator)));
    }
  }


  std::cout<<"Contructor end"<<std::endl;
}



void GEMMuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

  //Making the directory to write plot pngs to
  mkdir(histoFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  //Histos for plotting
  Nevents_h = new TH1F("Nevents_h"      , "Nevents"   , 2, 0, 2 );

  Nevents=0;



}


GEMMuonAnalyzer::~GEMMuonAnalyzer(){}

void
GEMMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  Nevents_h->Fill(1);
  using namespace edm;


  if (UseAssociators) {
    edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
    //edm::Handle<reco::MuonToTrackingParticleAssociator> theAssociator;
    for (unsigned int w=0;w<associators.size();w++) {
      iEvent.getByToken(associators_Token[w], theAssociator);
      associator.push_back( theAssociator.product() );
    }
  }


  using namespace reco;
  Handle<GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);
  const GenParticleCollection genParticlesForChi2 = *(genParticles.product());

  unsigned int gensize=genParticles->size();

  Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByToken(trackingParticlesToken_, trackingParticles);


  if (UseGEMEtaCoverageMuons){
    //Section to turn off signal muons in the endcaps, to approximate a nu gun
    for(unsigned int i=0; i<gensize; ++i) {
      const reco::GenParticle& CurrentParticle=(*genParticles)[i];
      if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){  
	      if ( fabs( CurrentParticle.eta() ) < 1.5 || fabs( CurrentParticle.eta() ) > 2.4 ) {
	      return;
	      }
      }
    }      
  }

  Nevents++;


  Handle <TrackCollection > generalTracks;
  iEvent.getByToken (generalTracksToken_, generalTracks);

  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);

  edm::Handle<GEMSegmentCollection> gemSegmentCollection;
  iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);

  edm::ESHandle<GEMGeometry> gemGeom;
  iSetup.get<MuonGeometryRecord>().get(gemGeom);

  ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);
  ESHandle<Propagator> shProp;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", shProp);
  

  //==== Make a vector of recoMuons whose isGEMMuon() and isTrackerMuon() is fired
  std::vector<bool> IsMatched;
  reco::MuonCollection GEMMuonColl;
  for(reco::MuonCollection::const_iterator recomuon=recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon) {
    if(recomuon->isGEMMuon() && recomuon->isTrackerMuon()){
      GEMMuonColl.push_back( *recomuon );
      IsMatched.push_back(false);
    }
  }

  //==== Finding GEMMuons that match gen muons, plotting the closest of those

  for(unsigned int i=0; i<gensize; i++){

    const reco::GenParticle& CurrentParticle=(*genParticles)[i];

    if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){   

      double LowestDelR = 9999;
      double thisDelR = 9999;
      int MatchedID = -1;
      int GEMMuonID = 0;
      double VertexDiff=-1,PtDiff=-1,QOverPtDiff=-1,PDiff=-1;

      for(reco::MuonCollection::const_iterator gemmuon = GEMMuonColl.begin(); gemmuon != GEMMuonColl.end(); ++gemmuon){
        TrackRef tkRef = gemmuon->innerTrack();
        thisDelR = reco::deltaR(CurrentParticle,*tkRef);

        if( tkRef->pt() > FakeRatePtCut ){
          if( thisDelR < MatchingWindowDelR ){
            if( thisDelR < LowestDelR ){
              LowestDelR = thisDelR;
              MatchedID = GEMMuonID;
              VertexDiff = fabs(tkRef->vz()-CurrentParticle.vz());
              QOverPtDiff = ( (tkRef->charge() /tkRef->pt()) - (CurrentParticle.charge()/CurrentParticle.pt() ) )/  (CurrentParticle.charge()/CurrentParticle.pt() );
              PtDiff = (tkRef->pt() - CurrentParticle.pt())/CurrentParticle.pt();
              PDiff = (tkRef->p() - CurrentParticle.p())/CurrentParticle.p();
            }
          }
        }

        GEMMuonID++;

      } // END gemmuon loop

      // GEMMuon matched to gen muon
      if( MatchedID != -1){

        IsMatched[MatchedID] = true;

        if ((CurrentParticle.pt() >FakeRatePtCut) ){
          MatchedGEMMuon_Eta->Fill(fabs(CurrentParticle.eta()));
          if ( (TMath::Abs(CurrentParticle.eta()) > 1.5) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedGEMMuon_Pt->Fill(CurrentParticle.pt());
          }
        }

        VertexDiff_h->Fill(VertexDiff);
        PtDiff_h->Fill(PtDiff);
        QOverPtDiff_h->Fill(QOverPtDiff);
        PtDiff_s->Fill(CurrentParticle.eta(),PtDiff);
        PDiff_h->Fill(PDiff);
        PDiff_s->Fill(CurrentParticle.eta(),PDiff);
        PDiff_p->Fill(CurrentParticle.eta(),PDiff);
      
      }

      if ( (CurrentParticle.pt() >FakeRatePtCut) ){
        GenMuon_Eta->Fill(fabs(CurrentParticle.eta()));
        if ( (fabs(CurrentParticle.eta()) > 1.5) && (fabs(CurrentParticle.eta()) < 2.4) ) {
          GenMuon_Pt->Fill(CurrentParticle.pt());
        }
      }



    } // END prompt muons selection
  } // END gen particle loop
  



}


void GEMMuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) 

{
  
}

DEFINE_FWK_MODULE(GEMMuonAnalyzer);
