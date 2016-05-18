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
#include "TEfficiency.h"

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
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void endRun(edm::Run const&, edm::EventSetup const&);
  std::string DoubleToString(std::string prefix, double dd);

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
  bool doMatchingStudy;
  const TrackAssociatorByChi2Impl* associatorByChi2;



  std::vector<std::string> associators;
  std::vector<edm::InputTag> label;

  //Histos for plotting
  TFile* histoFile; 

  double  FakeRatePtCut, MatchingWindowDelR;

  TH1F* Nevents_h;

  TH1F *GenMuon_Eta; TH1F *GenMuon_Pt; TH1F* MatchedGEMMuon_Eta; TH1F* MatchedGEMMuon_Pt;
  TH1F *TPMuon_Eta; TH1F *TPMuon_Pt;
  TH1F *HitsMatchedGEMMuon_Eta; TH1F *HitsMatchedGEMMuon_Pt;
  TH1F *HitsUnmatchedGEMMuon_Eta; TH1F* HitsUnmatchedGEMMuon_Pt;

  TH1F *DelX_GE11, *DelX_over_sigma_GE11, *DelY_GE11, *DelY_over_sigma_GE11, *DotDir_GE11;
  TH1F *DelX_GE21, *DelX_over_sigma_GE21, *DelY_GE21, *DelY_over_sigma_GE21, *DotDir_GE21;

  std::vector<double> maxPull, maxX_GE11, maxY_GE11, maxX_GE21, maxY_GE21, minDotDir;
  TH1F *maxPull_values, *maxX_GE11_values, *maxY_GE11_values, *maxX_GE21_values, *maxY_GE21_values, *minDotDir_values;
  std::map< std::string, TH1F* > map_maxXPull_GE11, map_maxYPull_GE11, map_maxXPull_GE21, map_maxYPull_GE21, map_maxX_GE11, map_maxY_GE11, map_maxX_GE21, map_maxY_GE21, map_minDotDir_GE11, map_minDotDir_GE21;

};

GEMMuonAnalyzer::GEMMuonAnalyzer(const edm::ParameterSet& iConfig) 
{
  histoFile = new TFile(iConfig.getParameter<std::string>("HistoFile").c_str(), "recreate");
  UseGEMEtaCoverageMuons = iConfig.getParameter< bool >("UseGEMEtaCoverageMuons");
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  doMatchingStudy = iConfig.getParameter< bool >("doMatchingStudy");

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
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));

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

  if(doMatchingStudy){
    generalTracksToken_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
    GEMSegment_Token = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));

    maxPull = iConfig.getParameter< std::vector<double> >("maxPull");
    maxX_GE11 = iConfig.getParameter< std::vector<double> >("maxX_GE11");
    maxY_GE11 = iConfig.getParameter< std::vector<double> >("maxY_GE11");
    maxX_GE21 = iConfig.getParameter< std::vector<double> >("maxX_GE21");
    maxY_GE21 = iConfig.getParameter< std::vector<double> >("maxY_GE21");
    minDotDir = iConfig.getParameter< std::vector<double> >("minDotDir");

  }


  std::cout<<"Contructor end"<<std::endl;
}



void GEMMuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

  const int n_pt_bin = 19, n_eta_bin = 9;
  double pt_bin[n_pt_bin+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  double eta_bin[n_eta_bin+1] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

  //Histos for plotting
  Nevents_h = new TH1F("Nevents_h", "Nevents", 2, 0, 2 );
  GenMuon_Eta = new TH1F("GenMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  GenMuon_Pt = new TH1F("GenMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  MatchedGEMMuon_Eta = new TH1F("MatchedGEMMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  MatchedGEMMuon_Pt =  new TH1F("MatchedGEMMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  TPMuon_Eta = new TH1F("TPMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  TPMuon_Pt = new TH1F("TPMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedGEMMuon_Eta = new TH1F("HitsMatchedGEMMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  HitsMatchedGEMMuon_Pt  = new TH1F("HitsMatchedGEMMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedGEMMuon_Eta = new TH1F("HitsUnmatchedGEMMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedGEMMuon_Pt  = new TH1F("HitsUnmatchedGEMMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );

  if(doMatchingStudy){
    DelX_GE11 = new TH1F("DelX_GE11", "DelX_GE11", 5./0.1, 0., 5.);
    DelX_over_sigma_GE11 = new TH1F("DelX_over_sigma_GE11", "DelX_over_sigma_GE11", 5./0.1, 0., 5.);
    DelY_GE11 = new TH1F("DelY_GE11", "DelY_GE11", 20./0.1, 0., 20.);
    DelY_over_sigma_GE11 = new TH1F("DelY_over_sigma_GE11", "DelY_over_sigma_GE11", 5./0.1, 0., 5.);
    DotDir_GE11 = new TH1F("DotDir_GE11", "DotDir_GE11", 1.5/0.01, 0., 1.5);

    DelX_GE21 = new TH1F("DelX_GE21", "DelX_GE21", 5./0.1, 0., 5.);
    DelX_over_sigma_GE21 = new TH1F("DelX_over_sigma_GE21", "DelX_over_sigma_GE21", 5./0.1, 0., 5.);
    DelY_GE21 = new TH1F("DelY_GE21", "DelY_GE21", 20./0.1, 0., 20.);
    DelY_over_sigma_GE21 = new TH1F("DelY_over_sigma_GE21", "DelY_over_sigma_GE21", 5./0.1, 0., 5.);
    DotDir_GE21 = new TH1F("DotDir_GE21", "DotDir_GE21", 1.5/0.01, 0., 1.5);

    maxPull_values = new TH1F("maxPull_values", "maxPull_values", 100, 0, 100);
    maxX_GE11_values = new TH1F("maxX_GE11_values", "maxX_GE11_values", 100, 0, 100);
    maxY_GE11_values = new TH1F("maxY_GE11_values", "maxY_GE11_values", 100, 0, 100);
    maxX_GE21_values = new TH1F("maxX_GE21_values", "maxX_GE21_values", 100, 0, 100);
    maxY_GE21_values = new TH1F("maxY_GE21_values", "maxY_GE21_values", 100, 0, 100);
    minDotDir_values = new TH1F("minDotDir_values", "minDotDir_values", 100, 0, 100);

    for(unsigned int aaa=0; aaa<maxPull.size(); aaa++){
      maxPull_values->SetBinContent(aaa+1, maxPull.at(aaa));
      map_maxXPull_GE11[DoubleToString("maxXPull_GE11", maxPull.at(aaa))] = new TH1F(DoubleToString("maxXPull_GE11", maxPull.at(aaa)).data(), DoubleToString("maxXPull_GE11", maxPull.at(aaa)).data(), 2, 0, 2);
      map_maxYPull_GE11[DoubleToString("maxYPull_GE11", maxPull.at(aaa))] = new TH1F(DoubleToString("maxYPull_GE11", maxPull.at(aaa)).data(), DoubleToString("maxYPull_GE11", maxPull.at(aaa)).data(), 2, 0, 2);
      map_maxXPull_GE21[DoubleToString("maxXPull_GE21", maxPull.at(aaa))] = new TH1F(DoubleToString("maxXPull_GE21", maxPull.at(aaa)).data(), DoubleToString("maxXPull_GE21", maxPull.at(aaa)).data(), 2, 0, 2);
      map_maxYPull_GE21[DoubleToString("maxYPull_GE21", maxPull.at(aaa))] = new TH1F(DoubleToString("maxYPull_GE21", maxPull.at(aaa)).data(), DoubleToString("maxYPull_GE21", maxPull.at(aaa)).data(), 2, 0, 2);
    }
    for(unsigned int aaa=0; aaa<maxX_GE11.size(); aaa++){
      maxX_GE11_values->SetBinContent(aaa+1, maxX_GE11.at(aaa));
      std::string thistitle = DoubleToString("maxX_GE11", maxX_GE11.at(aaa));
      map_maxX_GE11[thistitle] = new TH1F(thistitle.data(), thistitle.data(), 2, 0, 2);
    }
    for(unsigned int aaa=0; aaa<maxY_GE11.size(); aaa++){
      maxY_GE11_values->SetBinContent(aaa+1, maxY_GE11.at(aaa));
      std::string thistitle = DoubleToString("maxY_GE11", maxY_GE11.at(aaa));
      map_maxY_GE11[thistitle] = new TH1F(thistitle.data(), thistitle.data(), 2, 0, 2);
    }
    for(unsigned int aaa=0; aaa<maxX_GE21.size(); aaa++){
      maxX_GE21_values->SetBinContent(aaa+1, maxX_GE21.at(aaa));
      std::string thistitle = DoubleToString("maxX_GE21", maxX_GE21.at(aaa));
      map_maxX_GE21[thistitle] = new TH1F(thistitle.data(), thistitle.data(), 2, 0, 2);
    }
    for(unsigned int aaa=0; aaa<maxY_GE21.size(); aaa++){
      maxY_GE21_values->SetBinContent(aaa+1, maxY_GE21.at(aaa));
      std::string thistitle = DoubleToString("maxY_GE21", maxY_GE21.at(aaa));
      map_maxY_GE21[thistitle] = new TH1F(thistitle.data(), thistitle.data(), 2, 0, 2);
    }
    for(unsigned int aaa=0; aaa<minDotDir.size(); aaa++){
      minDotDir_values->SetBinContent(aaa+1, minDotDir.at(aaa));
      map_minDotDir_GE11[DoubleToString("minDotDir_GE11", minDotDir.at(aaa))] = new TH1F(DoubleToString("minDotDir_GE11", minDotDir.at(aaa)).data(), DoubleToString("minDotDir_GE11", minDotDir.at(aaa)).data(), 2, 0, 2);
      map_minDotDir_GE21[DoubleToString("minDotDir_GE21", minDotDir.at(aaa))] = new TH1F(DoubleToString("minDotDir_GE21", minDotDir.at(aaa)).data(), DoubleToString("minDotDir_GE21", minDotDir.at(aaa)).data(), 2, 0, 2);
    }
    
  }



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

  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);

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
      if(recSimColl.find(track) == recSimColl.end()) continue; //FIXME
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




  //==== Matching Study
  
  if(doMatchingStudy){

    Handle <TrackCollection > generalTracks;
    iEvent.getByToken (generalTracksToken_, generalTracks);
    edm::Handle<GEMSegmentCollection> gemSegmentCollection;
    iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);

    edm::ESHandle<GEMGeometry> gemGeom;
    iSetup.get<MuonGeometryRecord>().get(gemGeom);
    ESHandle<MagneticField> bField;
    iSetup.get<IdealMagneticFieldRecord>().get(bField);
    const SteppingHelixPropagator* shPropagator;
    shPropagator = new SteppingHelixPropagator(&*bField,alongMomentum);

    for (std::vector<Track>::const_iterator thisTrack = generalTracks->begin();
         thisTrack != generalTracks->end(); ++thisTrack){

      if (thisTrack->pt() < 1.5) continue;
      if (std::abs(thisTrack->eta()) < 1.5) continue;

      for (auto thisSegment = gemSegmentCollection->begin(); thisSegment != gemSegmentCollection->end();
           ++thisSegment){

        GEMDetId id = thisSegment->specificRecHits()[0].gemId();
        int station = id.station();
        if (id.station() != station) continue;
        float zSign = thisTrack->pz() > 0 ? 1.0f : -1.0f;
        if ( zSign * id.region() < 0 ) continue;

        LocalPoint thisPosition(thisSegment->localPosition());
        LocalVector thisDirection(thisSegment->localDirection());

        auto chamber = gemGeom->superChamber(id);
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
        Double_t DotDir = p3FinalReco.unit().dot(thisDirection);

        if(station == 1){
          DelX_GE11->Fill(DelX);
          DelX_over_sigma_GE11->Fill(DelX_over_sigma);
          DelY_GE11->Fill(DelY);
          DelY_over_sigma_GE11->Fill(DelY_over_sigma);
          DotDir_GE11->Fill(DotDir);
          for(unsigned int aaa=0; aaa<maxPull.size(); aaa++){
            int matchX = DelX_over_sigma<maxPull.at(aaa)?1:0;    
            int matchY = DelY_over_sigma<maxPull.at(aaa)?1:0;
            map_maxXPull_GE11[DoubleToString("maxXPull_GE11", maxPull.at(aaa))]->Fill(matchX);
            map_maxYPull_GE11[DoubleToString("maxYPull_GE11", maxPull.at(aaa))]->Fill(matchY);
          }
          for(unsigned int aaa=0; aaa<maxX_GE11.size(); aaa++){
            int match = DelX<maxX_GE11.at(aaa)?1:0;
            std::string thistitle = DoubleToString("maxX_GE11", maxX_GE11.at(aaa));
            map_maxX_GE11[thistitle]->Fill(match);
          }
          for(unsigned int aaa=0; aaa<maxY_GE11.size(); aaa++){
            int match = DelY<maxY_GE11.at(aaa)?1:0;
            std::string thistitle = DoubleToString("maxY_GE11", maxY_GE11.at(aaa));
            map_maxY_GE11[thistitle]->Fill(match);
          }
          for(unsigned int aaa=0; aaa<minDotDir.size(); aaa++){
            int match = DotDir>minDotDir.at(aaa)?1:0;
            map_minDotDir_GE11[DoubleToString("minDotDir_GE11", minDotDir.at(aaa))]->Fill(match);
          }
        }
        if(station == 3){
          DelX_GE21->Fill(DelX);
          DelX_over_sigma_GE21->Fill(DelX_over_sigma);
          DelY_GE21->Fill(DelY);
          DelY_over_sigma_GE21->Fill(DelY_over_sigma);
          DotDir_GE21->Fill(DotDir);
          for(unsigned int aaa=0; aaa<maxPull.size(); aaa++){
            int matchX = DelX_over_sigma<maxPull.at(aaa)?1:0;
            int matchY = DelY_over_sigma<maxPull.at(aaa)?1:0;
            map_maxXPull_GE21[DoubleToString("maxXPull_GE21", maxPull.at(aaa))]->Fill(matchX);
            map_maxYPull_GE21[DoubleToString("maxYPull_GE21", maxPull.at(aaa))]->Fill(matchY);
          }
          for(unsigned int aaa=0; aaa<maxX_GE21.size(); aaa++){
            int match = DelX<maxX_GE21.at(aaa)?1:0;
            std::string thistitle = DoubleToString("maxX_GE21", maxX_GE21.at(aaa));
            map_maxX_GE21[thistitle]->Fill(match);
          }
          for(unsigned int aaa=0; aaa<maxY_GE21.size(); aaa++){
            int match = DelY<maxY_GE21.at(aaa)?1:0;
            std::string thistitle = DoubleToString("maxY_GE21", maxY_GE21.at(aaa));
            map_maxY_GE21[thistitle]->Fill(match);
          }
          for(unsigned int aaa=0; aaa<minDotDir.size(); aaa++){
            int match = DotDir>minDotDir.at(aaa)?1:0;
            map_minDotDir_GE21[DoubleToString("minDotDir_GE21", minDotDir.at(aaa))]->Fill(match);
          }
        }


      } // END gemSegment loop


    } // END general track loop


  } // END if(doMatchingStudy)

}


void GEMMuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) 

{

  histoFile->cd();

  Nevents_h->Write();

  /* gen-reco delta R matching */
  GenMuon_Eta->Write();
  GenMuon_Pt->Write();
  MatchedGEMMuon_Eta->Write();
  MatchedGEMMuon_Pt->Write();
  // Efficiecny
  TEfficiency* Eff_Eta = new TEfficiency(*MatchedGEMMuon_Eta, *GenMuon_Eta);
  TEfficiency* Eff_Pt = new TEfficiency(*MatchedGEMMuon_Pt, *GenMuon_Pt);
  Eff_Eta->SetName("Eff_Eta");
  Eff_Pt->SetName("Eff_Pt");
  Eff_Eta->Write();
  Eff_Pt->Write();
  
  /* Association by hits */
  TPMuon_Eta->Write();
  TPMuon_Pt->Write();
  HitsMatchedGEMMuon_Eta->Write();
  HitsMatchedGEMMuon_Pt->Write();
  // Efficeicny
  TEfficiency* HitsEff_Eta = new TEfficiency(*HitsMatchedGEMMuon_Eta, *TPMuon_Eta);
  TEfficiency* HitsEff_Pt = new TEfficiency(*HitsMatchedGEMMuon_Pt, *TPMuon_Pt);
  HitsEff_Eta->SetName("HitsEff_Eta");
  HitsEff_Pt->SetName("HitsEff_Pt");
  HitsEff_Eta->Write();
  HitsEff_Pt->Write();
  // Fake
  HitsUnmatchedGEMMuon_Eta->Write();
  HitsUnmatchedGEMMuon_Pt->Write();

  /* Matching Study */
  if(doMatchingStudy){

    DelX_GE11->Write();
    DelX_over_sigma_GE11->Write();
    DelY_GE11->Write();
    DelY_over_sigma_GE11->Write();
    DotDir_GE11->Write();
    DelX_GE21->Write();
    DelX_over_sigma_GE21->Write();
    DelY_GE21->Write();
    DelY_over_sigma_GE21->Write();
    DotDir_GE21->Write();

    maxPull_values->Write();
    maxX_GE11_values->Write();
    maxY_GE11_values->Write();
    maxX_GE21_values->Write();
    maxY_GE21_values->Write();
    minDotDir_values->Write();

    for(unsigned int aaa=0; aaa<maxPull.size(); aaa++){
      map_maxXPull_GE11[DoubleToString("maxXPull_GE11", maxPull.at(aaa))]->Write();
      map_maxYPull_GE11[DoubleToString("maxYPull_GE11", maxPull.at(aaa))]->Write();
      map_maxXPull_GE21[DoubleToString("maxXPull_GE21", maxPull.at(aaa))]->Write();
      map_maxYPull_GE21[DoubleToString("maxYPull_GE21", maxPull.at(aaa))]->Write();
    }
    for(unsigned int aaa=0; aaa<maxX_GE11.size(); aaa++){
      std::string thistitle = DoubleToString("maxX_GE11", maxX_GE11.at(aaa));
      map_maxX_GE11[thistitle]->Write();
    }
    for(unsigned int aaa=0; aaa<maxY_GE11.size(); aaa++){
      std::string thistitle = DoubleToString("maxY_GE11", maxY_GE11.at(aaa));
      map_maxY_GE11[thistitle]->Write();
    }
    for(unsigned int aaa=0; aaa<maxX_GE21.size(); aaa++){
      std::string thistitle = DoubleToString("maxX_GE21", maxX_GE21.at(aaa));
      map_maxX_GE21[thistitle]->Write();
    }
    for(unsigned int aaa=0; aaa<maxY_GE21.size(); aaa++){
      std::string thistitle = DoubleToString("maxY_GE21", maxY_GE21.at(aaa));
      map_maxY_GE21[thistitle]->Write();
    }
    for(unsigned int aaa=0; aaa<minDotDir.size(); aaa++){
      map_minDotDir_GE11[DoubleToString("minDotDir_GE11", minDotDir.at(aaa))]->Write();
      map_minDotDir_GE21[DoubleToString("minDotDir_GE21", minDotDir.at(aaa))]->Write();
    }

  } // END if(doMatchingStudy)




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

std::string GEMMuonAnalyzer::DoubleToString(std::string prefix, double dd){
  std::ostringstream os;
  os << dd;
  return prefix+"_"+os.str();
}

DEFINE_FWK_MODULE(GEMMuonAnalyzer);
