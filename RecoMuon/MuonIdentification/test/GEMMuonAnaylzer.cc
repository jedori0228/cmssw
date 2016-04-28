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


  bool UseAssociators;
  bool UseGEMEtaCoverageMuons;
  const TrackAssociatorByChi2Impl* associatorByChi2;



  std::vector<std::string> associators;
  std::vector<edm::InputTag> label;

  //Histos for plotting
  TString histoFolder;
  TFile* histoFile; 

  double  FakeRatePtCut, MatchingWindowDelR;

  double Nevents;

  TH1F* Nevents_h;

  TH1F *GenMuon_Eta; TH1F *GenMuon_Pt; TH1F* MatchedGEMMuon_Eta; TH1F* MatchedGEMMuon_Pt;
  TH1F *TPMuon_Eta; TH1F *TPMuon_Pt;
  TH1F *HitsMatchedGEMMuon_Eta; TH1F *HitsMatchedGEMMuon_Pt;
  TH1F *HitsUnmatchedGEMMuon_Eta; TH1F* HitsUnmatchedGEMMuon_Pt;


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
      consumes<reco::RecoToSimCollection>(edm::InputTag(thisassociator));
      consumes<reco::SimToRecoCollection>(edm::InputTag(thisassociator));
    }
  }


  std::cout<<"Contructor end"<<std::endl;
}



void GEMMuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

  //Making the directory to write plot pngs to
  //mkdir(histoFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  //Histos for plotting
  Nevents_h = new TH1F("Nevents_h", "Nevents", 2, 0, 2 );
  GenMuon_Eta = new TH1F("GenMuon_Eta", "Muon #eta", 9, 1.5, 2.4 );
  GenMuon_Pt = new TH1F("GenMuon_Pt", "Muon p_{T}", 100, 0., 100. );
  MatchedGEMMuon_Eta = new TH1F("MatchedGEMMuon_Eta", "Muon #eta", 9, 1.5, 2.4 );
  MatchedGEMMuon_Pt =  new TH1F("MatchedGEMMuon_Pt", "Muon p_{T}", 100, 0., 100. );
  TPMuon_Eta = new TH1F("TPMuon_Eta", "Muon #eta", 9, 1.5, 2.4 );
  TPMuon_Pt = new TH1F("TPMuon_Pt", "Muon p_{T}", 100, 0. , 100. );
  HitsMatchedGEMMuon_Eta = new TH1F("HitsMatchedGEMMuon_Eta", "Muon #eta", 9, 1.5, 2.4 );
  HitsMatchedGEMMuon_Pt  = new TH1F("HitsMatchedGEMMuon_Pt", "Muon p_{T}", 100,0 , 100. );
  HitsUnmatchedGEMMuon_Eta = new TH1F("HitsUnmatchedGEMMuon_Eta", "Muon #eta", 9, 1.5, 2.4 );
  HitsUnmatchedGEMMuon_Pt  = new TH1F("HitsUnmatchedGEMMuon_Pt", "Muon p_{T}", 100,0 , 100. );


  Nevents=0;



}


GEMMuonAnalyzer::~GEMMuonAnalyzer(){}

void
GEMMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  Nevents_h->Fill(1);
  using namespace edm;
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
	      if ( fabs( CurrentParticle.eta() ) < 1.6 || fabs( CurrentParticle.eta() ) > 2.4 ) {
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
      for(reco::MuonCollection::const_iterator gemmuon = GEMMuonColl.begin(); gemmuon != GEMMuonColl.end(); ++gemmuon){
        TrackRef tkRef = gemmuon->innerTrack();
        thisDelR = reco::deltaR(CurrentParticle,*tkRef);

        if( tkRef->pt() > FakeRatePtCut ){
          if( thisDelR < MatchingWindowDelR ){
            if( thisDelR < LowestDelR ){
              LowestDelR = thisDelR;
              MatchedID = GEMMuonID;
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
          if ( (TMath::Abs(CurrentParticle.eta()) > 1.6) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedGEMMuon_Pt->Fill(CurrentParticle.pt());
          }
        }
      }

      if ( (CurrentParticle.pt() >FakeRatePtCut) ){
        GenMuon_Eta->Fill(fabs(CurrentParticle.eta()));
        if ( (fabs(CurrentParticle.eta()) > 1.6) && (fabs(CurrentParticle.eta()) < 2.4) ) {
          GenMuon_Pt->Fill(CurrentParticle.pt());
        }
      }



    } // END prompt muons selection
  } // END gen particle loop
  
  
  if (UseAssociators) {

    reco::RecoToSimCollection recSimColl;
    reco::SimToRecoCollection simRecColl;
    edm::Handle<View<Track> >  trackCollection;

	  Handle<reco::SimToRecoCollection > simtorecoCollectionH;
	  iEvent.getByLabel(associators[0],simtorecoCollectionH);
	  simRecColl= *(simtorecoCollectionH.product()); 
	
	  Handle<reco::RecoToSimCollection > recotosimCollectionH;
	  iEvent.getByLabel(associators[0],recotosimCollectionH);
	  recSimColl= *(recotosimCollectionH.product());

    unsigned int trackCollectionSize = 0;
    iEvent.getByToken(track_Collection_Token[0], trackCollection);
    trackCollectionSize = trackCollection->size();

    // denominators for efficiencies
    for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++){
      TrackingParticleRef tpr(trackingParticles, i);
      TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get()); 
      TrackingParticle::Vector momentumTP; 
      TrackingParticle::Point vertexTP;
    
      if (abs(tp->pdgId()) != 13) continue;

      bool Eta_1p6_2p4 = fabs(tp->eta()) > 1.6 && fabs(tp->eta()) < 2.4,
           Pt_5 = tp->pt() > 5;
      if( Eta_1p6_2p4 && Pt_5 ){
        bool SignalMuon = false;
        if(tp->status() != -99){ // Pythia8 gen status : home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
          //int motherid=-1;
          if ((*tp->genParticle_begin())->numberOfMothers()>0)  {
            if ((*tp->genParticle_begin())->mother()->numberOfMothers()>0){
              //motherid=(*tp->genParticle_begin())->mother()->mother()->pdgId();
            }
          }  
          //std::cout<<"Mother ID = "<<motherid<<std::endl;

          if ( ( (tp->status()==1) && ( (*tp->genParticle_begin())->numberOfMothers()==0 ) )  ||
               ( (tp->status()==1) )      )    SignalMuon=true;
      
        } // END if(tp->status() != -99)        

        if(SignalMuon){
          TPMuon_Eta->Fill(fabs(tp->eta()));
          TPMuon_Pt->Fill(tp->pt());
        }

      } // END if( Eta_1p6_2p4 && Pt_5 )


    } // END for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++)


    // loop over our tracks
    //std::cout << "trackCollectionSize = " << trackCollectionSize << std::endl;
    for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
      //std::cout << i << "th trackCollection iterator" << std::endl;
      RefToBase<Track> track(trackCollection, i);

      std::vector<std::pair<TrackingParticleRef, double> > tp;
      std::vector<std::pair<TrackingParticleRef, double> > tpforfake;
      TrackingParticleRef tpr;
      TrackingParticleRef tprforfake;

      //Check if the track is associated to any gen particle
      bool TrackIsEfficient = false;

      if(recSimColl.find(track) != recSimColl.end()){
        tp = recSimColl[track];
        if (tp.size()!=0) {
          //std::cout << " recSimColl[track] size = " << tp.size() << std::endl;
          tpr = tp.begin()->first;

          //double assocChi2 = -(tp.begin()->second);
 
          //So this track is matched to a gen particle, lets get that gen particle now

          if ( (simRecColl.find(tpr) != simRecColl.end()) ){
            std::vector<std::pair<RefToBase<Track>, double> > rt;
            if(simRecColl[tpr].size() > 0){
              //std::cout << " simRecColl[tpr] size = " << simRecColl[tpr].size() << std::endl;
              rt=simRecColl[tpr];
              RefToBase<Track> bestrecotrackforeff = rt.begin()->first;
              //Only fill the efficiency histo if the track found matches up to a gen particle's best choice
              if ( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) ) {
                TrackIsEfficient=true;
                //This section fills the numerator of the efficiency calculation...

                bool Eta_1p6_2p4 = fabs(tpr->eta()) > 1.6 && fabs(tpr->eta()) < 2.4,
                     Pt_5 = tpr->pt() > 5;
                if( Eta_1p6_2p4 && Pt_5 ){
                 
                  bool SignalMuon=false;
                  
                  if (tpr->status() !=-99){
                    //int motherid=-1;
                    if ((*tpr->genParticle_begin())->numberOfMothers()>0)  {
                      if ((*tpr->genParticle_begin())->mother()->numberOfMothers()>0){
                        //motherid=(*tpr->genParticle_begin())->mother()->mother()->pdgId();
                      }
                    }
                    //std::cout<<"Mother ID = "<<motherid<<std::endl;
                    if( ( (tpr->status()==1) && ( (*tpr->genParticle_begin())->numberOfMothers()==0 ) )  ||
                        ( (tpr->status()==1) )
                    ) SignalMuon=true;
                  } // END if (tpr->status() !=-99)
                  if(SignalMuon){
                    HitsMatchedGEMMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedGEMMuon_Pt->Fill(tpr->pt());

                  } // END if(SignalMuon)

                } // END if( Eta_1p6_2p4 && Pt_5 )            

              } // END if ( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) )
            } // END if(simRecColl[tpr].size() > 0) 
          } // END  if ( (simRecColl.find(tpr) != simRecColl.end()) )

        } // END if (tp.size()!=0)
      } // END if(recSimColl.find(track) != recSimColl.end())

      // A simple way of measuring fake rate
      if (!TrackIsEfficient) {
        bool Eta_1p6_2p4 = fabs(tpr->eta()) > 1.6 && fabs(tpr->eta()) < 2.4,
             Pt_5 = tpr->pt() > 5;
        if( Eta_1p6_2p4 && Pt_5 ){
          HitsUnmatchedGEMMuon_Eta->Fill(fabs(track->eta()));
          HitsUnmatchedGEMMuon_Pt->Fill(track->pt());
        } 

      } // END if (!TrackIsEfficient)
    } // END for(View<Track>::size_type i=0; i<trackCollectionSize; ++i)

  } // END if (UseAssociators)

}


void GEMMuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) 

{

  histoFile->cd();

  Nevents_h->Write();
  GenMuon_Eta->Write();
  GenMuon_Pt->Write();
  MatchedGEMMuon_Eta->Write();
  MatchedGEMMuon_Pt->Write();
  TPMuon_Eta->Write();
  TPMuon_Pt->Write();
  HitsMatchedGEMMuon_Eta->Write();
  HitsMatchedGEMMuon_Pt->Write();
  HitsUnmatchedGEMMuon_Eta->Write();
  HitsUnmatchedGEMMuon_Pt->Write();
  
}

DEFINE_FWK_MODULE(GEMMuonAnalyzer);
