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
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
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
#include <DataFormats/GEMRecHit/interface/GEMSegmentCollection.h>

#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

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
#include "TEfficiency.h"
#include <TProfile.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <TLatex.h>
//#include "CMSStyle.C"
//#include "tdrstyle.C"
//#include "lumi.C"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"
#include "TGraph.h"

#include <sstream>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

const int n_pt_bin = 19, n_eta_bin = 9;
double pt_bin[n_pt_bin+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
double eta_bin[n_eta_bin+1] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

class GEMMuonAnalyzer : public edm::EDAnalyzer {
public:
  explicit GEMMuonAnalyzer(const edm::ParameterSet&);
  ~GEMMuonAnalyzer();
  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& ,
           int , const AlgebraicSymMatrix66& ,
           const MagneticField* );

  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& ,
           int , const AlgebraicSymMatrix55& ,
           const MagneticField* );
    void getFromFTS(const FreeTrajectoryState& ,
      GlobalVector& , GlobalVector& ,
      int& , AlgebraicSymMatrix66& );


  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  //virtual void beginJob(const edm::EventSetup&);
  void beginJob();

  //For Track Association


  //protected:

  private:

  edm::ESHandle<GEMGeometry> gemGeom;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
  std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > track_Collection_Token;
  edm::EDGetTokenT<GEMRecHitCollection> GEMRecHit_Token;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;

  bool UseAssociators;
  bool doGeometryStudy;

  std::vector<std::string> associators;
  //std::string parametersDefiner;
  
  std::vector<std::string> label;

  TH2F *GEMRecHit_GE11_GlobalPosition_scattered, *GEMSegment_GE11_GlobalPosition_scattered;
  TH2F *GEMRecHit_GE21_GlobalPosition_scattered, *GEMSegment_GE21_GlobalPosition_scattered;
  TH2F *GEMRecHit_GE11_LocalPosition_scattered, *GEMSegment_GE11_LocalPosition_scattered;
  TH2F *GEMRecHit_GE21_LocalPosition_scattered, *GEMSegment_GE21_LocalPosition_scattered;
  TH2F *GEMRecHit_GE11_odd_XZplane, *GEMRecHit_GE11_even_XZplane;
  TH2F *GEMRecHit_GE21_odd_XZplane, *GEMRecHit_GE21_even_XZplane;

  int n_GEMMuon, n_RecoMuon, n_LooseMuon, n_MediumMuon, n_TightMuon;
  int n_GEMMuon_ptcut, n_RecoMuon_ptcut, n_LooseMuon_ptcut, n_MediumMuon_ptcut, n_TightMuon_ptcut;

  TH1F *N_GEMMuon_h, *N_GEMMuon_ptcut_h;

  TH1F *TPMuon_Eta, *TPMuon_Pt, *TPMuon_Phi;
  TH1F *HitsGEMMuon_Eta, *HitsGEMMuon_Pt, *HitsGEMMuon_Phi;
  TH1F *HitsMatchedGEMMuon_Eta, *HitsMatchedGEMMuon_Pt, *HitsMatchedGEMMuon_Phi;
  TH1F *HitsUnmatchedGEMMuon_Eta,* HitsUnmatchedGEMMuon_Pt,* HitsUnmatchedGEMMuon_Phi;

  std::string rootFileName;
  std::unique_ptr<TFile> outputfile;


};

GEMMuonAnalyzer::GEMMuonAnalyzer(const edm::ParameterSet& iConfig)
{
  outputfile.reset(TFile::Open(iConfig.getParameter<std::string>("HistoFile").c_str(), "RECREATE"));
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");

  UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  doGeometryStudy = iConfig.getParameter< bool >("doGeometryStudy");
  associators = iConfig.getParameter< std::vector<std::string> >("associators");
  label = iConfig.getParameter< std::vector<std::string> >("label");

  edm::InputTag trackingParticlesTag ("mix", "MergedTrackTruth");
  trackingParticlesToken_ = consumes<TrackingParticleCollection>(trackingParticlesTag);
  GEMRecHit_Token = consumes<GEMRecHitCollection>(edm::InputTag("gemRecHits", "", "STARECO"));
  GEMSegment_Token = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments", "", "STARECO"));

  GEMRecHit_GE11_GlobalPosition_scattered = new TH2F("GEMRecHit_GE11_GlobalPosition_scattered", "GEMRecHit GE11 GlobalPosition", 800, -400., 400., 800, -400., 400.);
  GEMSegment_GE11_GlobalPosition_scattered = new TH2F("GEMSegment_GE11_GlobalPosition_scattered", "GEMSegment GE11 GlobalPosition", 800, -400., 400., 800, -400., 400.);
  GEMRecHit_GE21_GlobalPosition_scattered = new TH2F("GEMRecHit_GE21_GlobalPosition_scattered", "GEMRecHit GE21 GlobalPosition", 800, -400., 400., 800, -400., 400.);
  GEMSegment_GE21_GlobalPosition_scattered = new TH2F("GEMSegment_GE21_GlobalPosition_scattered", "GEMSegment GE21 GlobalPosition", 800, -400., 400., 800, -400., 400.);
  GEMRecHit_GE11_LocalPosition_scattered = new TH2F("GEMRecHit_GE11_LocalPosition_scattered", "GEMRecHit GE11 LocalPosition", 200, -100., 100., 20, -10., 10.);
  GEMSegment_GE11_LocalPosition_scattered = new TH2F("GEMSegment_GE11_LocalPosition_scattered", "GEMSegment GE11 LocalPosition", 200, -100., 100., 400, -200., 200.);
  GEMRecHit_GE21_LocalPosition_scattered = new TH2F("GEMRecHit_GE21_LocalPosition_scattered", "GEMRecHit GE21 LocalPosition", 200, -100., 100., 20, -10., 10.);
  GEMSegment_GE21_LocalPosition_scattered = new TH2F("GEMSegment_GE21_LocalPosition_scattered", "GEMSegment GE21 LocalPosition", 200, -100., 100., 400, -200., 200.);
  GEMRecHit_GE11_odd_XZplane = new TH2F("GEMRecHit_GE11_odd_XZplane", "GEMRecHit GE11 Odd Chamber XZ-plane",    20*10, 560., 580., 800, -400., 400.);
  GEMRecHit_GE11_even_XZplane = new TH2F("GEMRecHit_GE11_even_XZplane", "GEMRecHit GE11 Even Chamber XZ-plane", 20*10, 560., 580., 800, -400., 400.);
  GEMRecHit_GE21_odd_XZplane = new TH2F("GEMRecHit_GE21_odd_XZplane", "GEMRecHit GE21 Odd XZ-plane",    20*10, 790., 810., 800, -400., 400.);
  GEMRecHit_GE21_even_XZplane = new TH2F("GEMRecHit_GE21_even_XZplane", "GEMRecHit GE21 Even XZ-plane", 20*10, 790., 810., 800, -400., 400.);

  N_GEMMuon_h = new TH1F("N_GEMMuon_h", "Nevents", 1, 0, 1 );
  N_GEMMuon_ptcut_h = new TH1F("N_GEMMuon_ptcut_h", "Nevents", 1, 0, 1 );

  TPMuon_Eta = new TH1F("TPMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  TPMuon_Pt = new TH1F("TPMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  TPMuon_Phi = new TH1F("TPMuon_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi() );

  HitsGEMMuon_Eta = new TH1F("HitsGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsGEMMuon_Pt  = new TH1F("HitsGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsGEMMuon_Phi = new TH1F("HitsGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsMatchedGEMMuon_Eta = new TH1F("HitsMatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsMatchedGEMMuon_Pt  = new TH1F("HitsMatchedGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedGEMMuon_Phi = new TH1F("HitsMatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsUnmatchedGEMMuon_Eta = new TH1F("HitsUnmatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedGEMMuon_Pt  = new TH1F("HitsUnmatchedGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedGEMMuon_Phi = new TH1F("HitsUnmatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );

  n_GEMMuon = 0;
  n_RecoMuon = 0;
  n_LooseMuon = 0;
  n_MediumMuon = 0;
  n_TightMuon = 0;
  n_GEMMuon_ptcut = 0;
  n_RecoMuon_ptcut = 0;
  n_LooseMuon_ptcut = 0;
  n_MediumMuon_ptcut = 0;
  n_TightMuon_ptcut = 0;

  for (unsigned int www=0;www<label.size();www++){
    track_Collection_Token.push_back(consumes<edm::View<reco::Track> >( edm::InputTag( label[www] ) ));
  }

  if (UseAssociators) {
    for (auto const& thisassociator :associators) {
      //std::cout << thisassociator << std::endl;
      consumes<reco::RecoToSimCollection>(edm::InputTag(thisassociator));
      consumes<reco::SimToRecoCollection>(edm::InputTag(thisassociator));
    }
  }

}

void GEMMuonAnalyzer::beginJob()
{


}

GEMMuonAnalyzer::~GEMMuonAnalyzer(){


  outputfile->cd();

  GEMRecHit_GE11_GlobalPosition_scattered->Write();
  GEMSegment_GE11_GlobalPosition_scattered->Write();
  GEMRecHit_GE21_GlobalPosition_scattered->Write();
  GEMSegment_GE21_GlobalPosition_scattered->Write();
  GEMRecHit_GE11_LocalPosition_scattered->Write();
  GEMSegment_GE11_LocalPosition_scattered->Write();
  GEMRecHit_GE21_LocalPosition_scattered->Write();
  GEMSegment_GE21_LocalPosition_scattered->Write();
  GEMRecHit_GE11_odd_XZplane->Write();
  GEMRecHit_GE11_even_XZplane->Write();
  GEMRecHit_GE21_odd_XZplane->Write();
  GEMRecHit_GE21_even_XZplane->Write();

  N_GEMMuon_h->SetBinContent(1, n_GEMMuon);
  N_GEMMuon_ptcut_h->SetBinContent(1, n_GEMMuon_ptcut);
  N_GEMMuon_h->Write();
  N_GEMMuon_ptcut_h->Write();

  TPMuon_Eta->Write();
  TPMuon_Pt->Write();
  TPMuon_Phi->Write();

  //==== GEMMuon
  HitsGEMMuon_Eta->Write();
  HitsGEMMuon_Pt->Write();
  HitsGEMMuon_Phi->Write();
  HitsMatchedGEMMuon_Eta->Write();
  HitsMatchedGEMMuon_Pt->Write();
  HitsMatchedGEMMuon_Phi->Write();
  // Efficeicny
  TEfficiency* HitsEff_GEMMuon_Eta = new TEfficiency(*HitsMatchedGEMMuon_Eta, *TPMuon_Eta);
  TEfficiency* HitsEff_GEMMuon_Pt = new TEfficiency(*HitsMatchedGEMMuon_Pt, *TPMuon_Pt);
  TEfficiency* HitsEff_GEMMuon_Phi = new TEfficiency(*HitsMatchedGEMMuon_Phi, *TPMuon_Phi);
  HitsEff_GEMMuon_Eta->SetName("HitsEff_GEMMuon_Eta");
  HitsEff_GEMMuon_Pt->SetName("HitsEff_GEMMuon_Pt");
  HitsEff_GEMMuon_Phi->SetName("HitsEff_GEMMuon_Phi");
  HitsEff_GEMMuon_Eta->Write();
  HitsEff_GEMMuon_Pt->Write();
  HitsEff_GEMMuon_Phi->Write();
  // Fake
  HitsUnmatchedGEMMuon_Eta->Write();
  HitsUnmatchedGEMMuon_Pt->Write();
  HitsUnmatchedGEMMuon_Phi->Write();


  outputfile->Close();

}

void
GEMMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{
  using namespace edm;
  using namespace reco;

  iSetup.get<MuonGeometryRecord>().get(gemGeom);
  edm::Handle<GEMRecHitCollection> gemRecHitCollection;
  iEvent.getByToken(GEMRecHit_Token, gemRecHitCollection);
  edm::Handle<GEMSegmentCollection> gemSegmentCollection;
  iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);

  std::cout << "doGeometryStudy = " << doGeometryStudy << std::endl;
  if(doGeometryStudy){

    edm::LogVerbatim("GEMMuonAnalyzer") << "########### Geometry Study ###########";

    //==== GEMRecHit study

    for(auto thisRecHit = gemRecHitCollection->begin(); thisRecHit != gemRecHitCollection->end(); ++thisRecHit) {
      GEMDetId id = thisRecHit->gemId();
      auto roll = gemGeom->etaPartition(id);
      auto RecHitLP = thisRecHit->localPosition();
      auto RecHitGP = roll->toGlobal(RecHitLP);
      edm::LogVerbatim("GEMMuonAnalyzer") << "rechit station = " << id.station() << std::endl;
      edm::LogVerbatim("GEMMuonAnalyzer") << "Rechit id = " << id << " => (" << RecHitGP.x() << ", " << RecHitGP.y() << ", " << RecHitGP.z() << ")" << std::endl;
      if( id.station() == 1 ){
        GEMRecHit_GE11_LocalPosition_scattered->Fill(RecHitLP.x(), RecHitLP.y());
        GEMRecHit_GE11_GlobalPosition_scattered->Fill(RecHitGP.x(), RecHitGP.y());
        if(id.chamber()%2 == 0){
          GEMRecHit_GE11_even_XZplane->Fill(RecHitGP.z(), RecHitGP.x());
        }
        else{
          GEMRecHit_GE11_odd_XZplane->Fill(RecHitGP.z(), RecHitGP.x());
        }
      }
      if( id.station() == 3 ){
        GEMRecHit_GE21_LocalPosition_scattered->Fill(RecHitLP.x(), RecHitLP.y());
        GEMRecHit_GE21_GlobalPosition_scattered->Fill(RecHitGP.x(), RecHitGP.y());
        if(id.chamber()%2 == 0){
          GEMRecHit_GE21_even_XZplane->Fill(RecHitGP.z(), RecHitGP.x());
        }
        else{
          GEMRecHit_GE21_odd_XZplane->Fill(RecHitGP.z(), RecHitGP.x());
        }
      }
    }

    //==== GEMSegment study

    for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
      GEMDetId id = gems->gemDetId();
      auto chamb = gemGeom->superChamber(id);
      auto segLP = gems->localPosition();
      auto segGP = chamb->toGlobal(segLP);
      auto rechits = gems->specificRecHits();
      edm::LogVerbatim("GEMMuonAnalyzer") << "gemseg station = " << id.station() << std::endl;
      edm::LogVerbatim("GEMMuonAnalyzer") << "# of rechits used in this Segment = " << rechits.size() << std::endl;
      //for(auto rechit = rechits.begin(); rechit != rechits.end(); ++rechit){
      //  edm::LogVerbatim("GEMMuonAnalyzer") "  " << rechit->gemId() << "(" << segGP.x() << ", " << segGP.y() << ", " << segGP.z() << ")" << std::endl;
      //}
      if( id.station() == 1 ){
        GEMSegment_GE11_LocalPosition_scattered->Fill(segLP.x(), segLP.y());
        GEMSegment_GE11_GlobalPosition_scattered->Fill(segGP.x(), segGP.y());
      }
      if( id.station() == 3 ){
        GEMSegment_GE21_LocalPosition_scattered->Fill(segLP.x(), segLP.y());
        GEMSegment_GE21_GlobalPosition_scattered->Fill(segGP.x(), segGP.y());
      }
    }

  }



  //==== Track Association by hits:

  if(UseAssociators){

    Handle<TrackingParticleCollection> trackingParticles;
    iEvent.getByToken(trackingParticlesToken_, trackingParticles);

    for(unsigned int www=0;www<label.size();www++){     

      //=========================
      //==== prepare RecoTracks
      //========================= 

      Handle<reco::SimToRecoCollection > simtorecoCollectionH;
      iEvent.getByLabel(associators[www],simtorecoCollectionH);
      reco::SimToRecoCollection simRecColl= *(simtorecoCollectionH.product());

      Handle<reco::RecoToSimCollection > recotosimCollectionH;
      iEvent.getByLabel(associators[www],recotosimCollectionH);
      reco::RecoToSimCollection recSimColl= *(recotosimCollectionH.product());

      edm::Handle<View<Track> >  trackCollection;
      iEvent.getByToken(track_Collection_Token[www], trackCollection);

      //=================
      //==== Efficiency
      //=================

      for(TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++){

        TrackingParticleRef tpr(trackingParticles, i);
        TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get()); 
        TrackingParticle::Vector momentumTP; 
        TrackingParticle::Point vertexTP;
    
        if (abs(tp->pdgId()) != 13) continue;

        bool Eta_1p6_2p4 = abs(tp->eta()) > 1.6 && abs(tp->eta()) < 2.4;
        bool Pt_5 = tp->pt() > 5;

        if( Eta_1p6_2p4 ){
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
            //==== Fill the Denominator
            //==== Should be filled only once (www=0)
            if(www==0){
              if(Pt_5) TPMuon_Eta->Fill(fabs(tp->eta()));
              TPMuon_Pt->Fill(tp->pt());
              if(Pt_5) TPMuon_Phi->Fill(tp->phi());
            }

            if( (simRecColl.find(tpr) == simRecColl.end()) || (simRecColl[tpr].size() == 0) ){
              edm::LogVerbatim("GEMMuonAnalyzer") << "No SimToReco found for this TrackingParticle";
            }
            else{
              std::vector<std::pair<RefToBase<Track>, double> > rt = simRecColl[tpr];
              RefToBase<Track> rtr = rt.begin()->first;
              //std::cout << "This SimToReco :" << std::endl;
              for(std::vector<std::pair<RefToBase<Track>, double> >::const_iterator itit=rt.begin(); itit!=rt.end(); itit++){
                //std::cout << "  quality = "<<itit->second<<std::endl;
              }
              if( (recSimColl.find(rtr) == recSimColl.end()) || (recSimColl[rtr].size() ==0) ){
                edm::LogVerbatim("GEMMuonAnalyzer") << "SimToReco found, but no RecoToSim for the best SimToReco";
              }
              else{
                std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[rtr];
                TrackingParticleRef bestTPforEff = tp.begin()->first;
                if( bestTPforEff == tpr ){
                  edm::LogVerbatim("GEMMuonAnalyzer") << "Found matched RecoTrack";

                  if(label[www]=="gemMuonSel"){
                    if(Pt_5) HitsMatchedGEMMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedGEMMuon_Pt->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedGEMMuon_Phi->Fill(tpr->phi());
                  }

                }
              }
            }

          } // END Signal Muon
        } // END if( Eta_1p6_2p )


      } // END for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++)

      //=================
      //==== Fake stduy
      //=================

      //==== loop over our tracks
      for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
        RefToBase<Track> track(trackCollection, i);

        bool Eta_1p6_2p4 = fabs(track->eta()) > 1.6 && fabs(track->eta()) < 2.4;
        bool Pt_5 = track->pt() > 5;

        if( Eta_1p6_2p4 ){

          //==== count reco tracks
          if(label[www]=="gemMuonSel"){
            n_GEMMuon++;
            HitsGEMMuon_Pt->Fill(track->pt());
            if( Pt_5 ){
              n_GEMMuon_ptcut++;
              HitsGEMMuon_Eta->Fill(fabs(track->eta()));
              HitsGEMMuon_Phi->Fill(track->phi()); 
            }
          }

          //==== Check if the track is associated to any gen particle
          bool isFake = false;
          if( (recSimColl.find(track) == recSimColl.end()) || (recSimColl[track].size() == 0) ){
            edm::LogVerbatim("GEMMuonAnalyzer") << "No RecoToSim found for this RecoTrack";
            isFake = true;
          }
          //==== gen particle is found
          else{
            std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[track];
            TrackingParticleRef tpr = tp.begin()->first;
            //std::cout << " recSimColl[track] size = " << tp.size() << std::endl;
            if( (simRecColl.find(tpr) == simRecColl.end())  ||  (simRecColl[tpr].size() == 0)  ) {
              edm::LogVerbatim("GEMMuonAnalyzer") << "RecoToSim found, but no SimToReco for the best RecoToSim";
              isFake = true;
            }
            else{
              std::vector<std::pair<RefToBase<Track>, double> > rt = simRecColl[tpr];
              RefToBase<Track> bestrecotrackforeff = rt.begin()->first;
              //std::cout << " simRecColl[tpr] size = " << simRecColl[tpr].size() << std::endl;

              //std::cout << "  found matched genparticle => pdgid = " << tpr->pdgId() << std::endl;
              bool SignalMuon = false;
              if(tpr->status() !=-99){
                if ((*tpr->genParticle_begin())->numberOfMothers()>0)  {
                  if ((*tpr->genParticle_begin())->mother()->numberOfMothers()>0){
                    //motherid=(*tpr->genParticle_begin())->mother()->mother()->pdgId();
                  }
                }
                //std::cout<<"Mother ID = "<<motherid<<std::endl;
                if( ( (tpr->status()==1) && ( (*tpr->genParticle_begin())->numberOfMothers()==0 ) )  ||
                    ( (tpr->status()==1) ) ) SignalMuon=true;
              }
              if( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) && SignalMuon ) {
                edm::LogVerbatim("GEMMuonAnalyzer") << "Found matched TrackingParticle";
              }
              else{
                isFake = true;
              }

            }
          }

          if( isFake ){
            if(label[www]=="gemMuonSel"){
               if(Pt_5) HitsUnmatchedGEMMuon_Eta->Fill(fabs(track->eta()));
               HitsUnmatchedGEMMuon_Pt->Fill(track->pt());
               if(Pt_5) HitsUnmatchedGEMMuon_Phi->Fill(track->phi());
            }
          } //==== END if(isFake)

        } // ==== END if( Eta_1p6_2p4 )

      } // END track loop


    } // END for (unsigned int www=0;www<label.size();www++)

  } // END if (UseAssociators)


}

void GEMMuonAnalyzer::endJob()
{

}

FreeTrajectoryState
GEMMuonAnalyzer::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix55& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);

  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
GEMMuonAnalyzer::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix66& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);

  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void GEMMuonAnalyzer::getFromFTS(const FreeTrajectoryState& fts,
            GlobalVector& p3, GlobalVector& r3,
            int& charge, AlgebraicSymMatrix66& cov){
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

DEFINE_FWK_MODULE(GEMMuonAnalyzer);
