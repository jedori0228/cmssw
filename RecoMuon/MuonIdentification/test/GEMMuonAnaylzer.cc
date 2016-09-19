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

  //==== quick efficiency check
  Long64_t n_genmuon;
  Long64_t n_dR_matched_GEMmuon;
  Long64_t n_AssoByHits_matched_GEMmuon;

  edm::ESHandle<GEMGeometry> gemGeom;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
  edm::EDGetTokenT <reco::TrackCollection > generalTracksToken_;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<GEMRecHitCollection> GEMRecHit_Token;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > track_Collection_Token;

  bool UseAssociators;
  bool UseGEMEtaCoverageMuons;
  bool doMatchingStudy;
  const TrackAssociatorByChi2Impl* associatorByChi2;

  std::vector<std::string> associators;
  //std::vector<edm::InputTag> label;
  std::vector<std::string> label;

  //Histos for plotting
  TFile* histoFile; 

  double  FakeRatePtCut, MatchingWindowDelR;

  TH1F* Nevents_h;

  TH1F *GenMuon_Eta; TH1F *GenMuon_Pt; TH1F *GenMuon_Phi;
  TH1F *MatchedRecoMuon_Eta; TH1F *MatchedRecoMuon_Pt; TH1F *MatchedRecoMuon_Phi;
  TH1F *MatchedRecoMuon_not_GEMMuon_Eta; TH1F *MatchedRecoMuon_not_GEMMuon_Pt; TH1F *MatchedRecoMuon_not_GEMMuon_Phi;
  TH1F *MatchedRecoMuon_not_GEMMuon_no_gemseg_Eta; TH1F *MatchedRecoMuon_not_GEMMuon_no_gemseg_Pt; TH1F *MatchedRecoMuon_not_GEMMuon_no_gemseg_Phi;
  TH1F *MatchedGEMMuon_Eta; TH1F *MatchedGEMMuon_Pt; TH1F *MatchedGEMMuon_Phi;
  TH1F *MatchedGEMRecHit_Eta; TH1F *MatchedGEMRecHit_Pt; TH1F *MatchedGEMRecHit_Phi;
  TH1F *MatchedGEMRecHit_GE11_layer1_Eta; TH1F *MatchedGEMRecHit_GE11_layer1_Pt; TH1F *MatchedGEMRecHit_GE11_layer1_Phi;
  TH1F *MatchedGEMRecHit_GE11_layer2_Eta; TH1F *MatchedGEMRecHit_GE11_layer2_Pt; TH1F *MatchedGEMRecHit_GE11_layer2_Phi;
  TH1F *MatchedGEMRecHit_GE21_layer1_Eta; TH1F *MatchedGEMRecHit_GE21_layer1_Pt; TH1F *MatchedGEMRecHit_GE21_layer1_Phi;
  TH1F *MatchedGEMRecHit_GE21_layer2_Eta; TH1F *MatchedGEMRecHit_GE21_layer2_Pt; TH1F *MatchedGEMRecHit_GE21_layer2_Phi;
  TH1F *MatchedGEMRecHit_GE11_two_Eta; TH1F *MatchedGEMRecHit_GE11_two_Pt; TH1F *MatchedGEMRecHit_GE11_two_Phi;
  TH1F *MatchedGEMRecHit_GE21_two_Eta; TH1F *MatchedGEMRecHit_GE21_two_Pt; TH1F *MatchedGEMRecHit_GE21_two_Phi;
  TH1I *MatchedClusteredGEMRecHit_GE11_dBunchX; TH1I *MatchedClusteredGEMRecHit_GE21_dBunchX;
  TH1F *MatchedGEMSegment_GE11_Eta; TH1F *MatchedGEMSegment_GE11_Pt; TH1F *MatchedGEMSegment_GE11_Phi;
  TH1F *MatchedGEMSegment_GE21_Eta; TH1F *MatchedGEMSegment_GE21_Pt; TH1F *MatchedGEMSegment_GE21_Phi;
  TH1F *MatchedGEMSegment_Eta; TH1F *MatchedGEMSegment_Pt; TH1F *MatchedGEMSegment_Phi;
  TH1F *TPMuon_Eta; TH1F *TPMuon_Pt; TH1F *TPMuon_Phi;
  TH1F *HitsMatchedGEMMuon_Eta; TH1F *HitsMatchedGEMMuon_Pt; TH1F *HitsMatchedGEMMuon_Phi;
  TH1F *HitsUnmatchedGEMMuon_Eta; TH1F* HitsUnmatchedGEMMuon_Pt; TH1F* HitsUnmatchedGEMMuon_Phi;
  TH1F *HitsMatchedRecoMuon_Eta; TH1F *HitsMatchedRecoMuon_Pt; TH1F *HitsMatchedRecoMuon_Phi;
  TH1F *HitsUnmatchedRecoMuon_Eta; TH1F* HitsUnmatchedRecoMuon_Pt; TH1F* HitsUnmatchedRecoMuon_Phi;
  TH1F *HitsMatchedLooseMuon_Eta; TH1F *HitsMatchedLooseMuon_Pt; TH1F *HitsMatchedLooseMuon_Phi;
  TH1F *HitsUnmatchedLooseMuon_Eta; TH1F* HitsUnmatchedLooseMuon_Pt; TH1F* HitsUnmatchedLooseMuon_Phi;
  TH1F *HitsMatchedMediumMuon_Eta; TH1F *HitsMatchedMediumMuon_Pt; TH1F *HitsMatchedMediumMuon_Phi;
  TH1F *HitsUnmatchedMediumMuon_Eta; TH1F* HitsUnmatchedMediumMuon_Pt; TH1F* HitsUnmatchedMediumMuon_Phi;
  TH1F *HitsMatchedTightMuon_Eta; TH1F *HitsMatchedTightMuon_Pt; TH1F *HitsMatchedTightMuon_Phi;
  TH1F *HitsUnmatchedTightMuon_Eta; TH1F* HitsUnmatchedTightMuon_Pt; TH1F* HitsUnmatchedTightMuon_Phi;

  TH2F *GEMRecHit_GE11_GlobalPosition_scattered, *GEMSegment_GE11_GlobalPosition_scattered;
  TH2F *GEMRecHit_GE21_GlobalPosition_scattered, *GEMSegment_GE21_GlobalPosition_scattered;
  TH2F *GEMRecHit_GE11_LocalPosition_scattered, *GEMSegment_GE11_LocalPosition_scattered;
  TH2F *GEMRecHit_GE21_LocalPosition_scattered, *GEMSegment_GE21_LocalPosition_scattered;
  TH2F *GEMRecHit_GE11_odd_XZplane, *GEMRecHit_GE11_even_XZplane; 
  TH2F *GEMRecHit_GE21_odd_XZplane, *GEMRecHit_GE21_even_XZplane;

  TH1F *DelX_GE11, *DelX_over_sigma_GE11, *DelY_GE11, *DelY_over_sigma_GE11, *DotDir_GE11;
  TH1F *DelX_GE21, *DelX_over_sigma_GE21, *DelY_GE21, *DelY_over_sigma_GE21, *DotDir_GE21;

  std::vector<double> maxPull, maxX_GE11, maxY_GE11, maxX_GE21, maxY_GE21, minDotDir;
  TH1F *maxPull_values, *maxX_GE11_values, *maxY_GE11_values, *maxX_GE21_values, *maxY_GE21_values, *minDotDir_values;
  std::map< std::string, TH1F* > map_maxXPull_GE11, map_maxYPull_GE11, map_maxXPull_GE21, map_maxYPull_GE21, map_maxX_GE11, map_maxY_GE11, map_maxX_GE21, map_maxY_GE21, map_minDotDir_GE11, map_minDotDir_GE21;

  double Current_trackerGEM_maxPull, Current_maxDiffXGE11, Current_maxDiffYGE11, Current_maxDiffXGE21, Current_maxDiffYGE21, Current_minDotDir;

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

  //label = iConfig.getParameter< std::vector<edm::InputTag> >("label");
  label = iConfig.getParameter< std::vector<std::string> >("label");
  edm::InputTag genParticlesTag ("genParticles");
  genParticlesToken_ = consumes<reco::GenParticleCollection>(genParticlesTag);
  edm::InputTag trackingParticlesTag ("mix","MergedTrackTruth");
  trackingParticlesToken_ = consumes<TrackingParticleCollection>(trackingParticlesTag);
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  GEMRecHit_Token = consumes<GEMRecHitCollection>(edm::InputTag("gemRecHits"));
  GEMSegment_Token = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));

  //Getting tokens and doing consumers for track associators
  for (unsigned int www=0;www<label.size();www++){
    track_Collection_Token.push_back( consumes<edm::View<reco::Track> >( edm::InputTag(label[www]) ) );
  }

  if (UseAssociators) {
    for (auto const& thisassociator :associators) {
      //std::cout << thisassociator << std::endl;
      consumes<reco::RecoToSimCollection>(edm::InputTag(thisassociator));
      consumes<reco::SimToRecoCollection>(edm::InputTag(thisassociator));
    }
  }

  if(doMatchingStudy){
    generalTracksToken_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));

    maxPull = iConfig.getParameter< std::vector<double> >("maxPull");
    maxX_GE11 = iConfig.getParameter< std::vector<double> >("maxX_GE11");
    maxY_GE11 = iConfig.getParameter< std::vector<double> >("maxY_GE11");
    maxX_GE21 = iConfig.getParameter< std::vector<double> >("maxX_GE21");
    maxY_GE21 = iConfig.getParameter< std::vector<double> >("maxY_GE21");
    minDotDir = iConfig.getParameter< std::vector<double> >("minDotDir");

  }

   Current_trackerGEM_maxPull = iConfig.getParameter<double>("Current_trackerGEM_maxPull");
   Current_maxDiffXGE11   = iConfig.getParameter<double>("Current_maxDiffXGE11");
   Current_maxDiffYGE11   = iConfig.getParameter<double>("Current_maxDiffYGE11");
   Current_maxDiffXGE21   = iConfig.getParameter<double>("Current_maxDiffXGE21");
   Current_maxDiffYGE21   = iConfig.getParameter<double>("Current_maxDiffYGE21");
   Current_minDotDir      = iConfig.getParameter<double>("Current_minDotDir");

  n_genmuon = 0;
  n_dR_matched_GEMmuon = 0;
  n_AssoByHits_matched_GEMmuon = 0;

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
  GenMuon_Phi = new TH1F("GenMuon_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi());
  MatchedRecoMuon_Eta = new TH1F("MatchedRecoMuon_Eta", "RecoMuon #eta", n_eta_bin, eta_bin );
  MatchedRecoMuon_Pt =  new TH1F("MatchedRecoMuon_Pt", "RecoMuon p_{T}", n_pt_bin, pt_bin );
  MatchedRecoMuon_Phi =  new TH1F("MatchedRecoMuon_Phi", "RecoMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  MatchedRecoMuon_not_GEMMuon_Eta = new TH1F("MatchedRecoMuon_not_GEMMuon_Eta", "RecoMuon NOT GEMMuon #eta", n_eta_bin, eta_bin );
  MatchedRecoMuon_not_GEMMuon_Pt =  new TH1F("MatchedRecoMuon_not_GEMMuon_Pt", "RecoMuon NOT GEMMuon p_{T}", n_pt_bin, pt_bin );
  MatchedRecoMuon_not_GEMMuon_Phi = new TH1F("MatchedRecoMuon_not_GEMMuon_Phi", "RecoMuon NOT GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  MatchedRecoMuon_not_GEMMuon_no_gemseg_Eta = new TH1F("MatchedRecoMuon_not_GEMMuon_no_gemseg_Eta", "RecoMuon NOT GEMMuon, NO GEMSegment matched #eta", n_eta_bin, eta_bin );
  MatchedRecoMuon_not_GEMMuon_no_gemseg_Pt =  new TH1F("MatchedRecoMuon_not_GEMMuon_no_gemseg_Pt", "RecoMuon NOT GEMMuon, NO GEMSegment matched p_{T}", n_pt_bin, pt_bin );
  MatchedRecoMuon_not_GEMMuon_no_gemseg_Phi = new TH1F("MatchedRecoMuon_not_GEMMuon_no_gemseg_Phi", "RecoMuon NOT GEMMuon, NO GEMSegment matched #phi", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMMuon_Eta = new TH1F("MatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  MatchedGEMMuon_Pt =  new TH1F("MatchedGEMMuon_Pt", "GEMMuon p_{T}", n_pt_bin, pt_bin );
  MatchedGEMMuon_Phi = new TH1F("MatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMRecHit_Eta = new TH1F("MatchedGEMRecHit_Eta", "Muon #eta", n_eta_bin, eta_bin );
  MatchedGEMRecHit_Pt =  new TH1F("MatchedGEMRecHit_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  MatchedGEMRecHit_Phi = new TH1F("MatchedGEMRecHit_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMRecHit_GE11_layer1_Eta = new TH1F("MatchedGEMRecHit_GE11_layer1_Eta", "Muon #eta GE11 layer 1", n_eta_bin, eta_bin );
  MatchedGEMRecHit_GE11_layer1_Pt =  new TH1F("MatchedGEMRecHit_GE11_layer1_Pt", "Muon p_{T}  GE11 layer 1", n_pt_bin, pt_bin );
  MatchedGEMRecHit_GE11_layer1_Phi = new TH1F("MatchedGEMRecHit_GE11_layer1_Phi", "Muon #phi GE11 layer 1", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMRecHit_GE11_layer2_Eta = new TH1F("MatchedGEMRecHit_GE11_layer2_Eta", "Muon #eta GE11 layer 2", n_eta_bin, eta_bin );
  MatchedGEMRecHit_GE11_layer2_Pt =  new TH1F("MatchedGEMRecHit_GE11_layer2_Pt", "Muon p_{T}  GE11 layer 2", n_pt_bin, pt_bin );
  MatchedGEMRecHit_GE11_layer2_Phi = new TH1F("MatchedGEMRecHit_GE11_layer2_Phi", "Muon #phi GE11 layer 2", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMRecHit_GE21_layer1_Eta = new TH1F("MatchedGEMRecHit_GE21_layer1_Eta", "Muon #eta GE21 layer 1", n_eta_bin, eta_bin );
  MatchedGEMRecHit_GE21_layer1_Pt =  new TH1F("MatchedGEMRecHit_GE21_layer1_Pt", "Muon p_{T}  GE21 layer 1", n_pt_bin, pt_bin );
  MatchedGEMRecHit_GE21_layer1_Phi = new TH1F("MatchedGEMRecHit_GE21_layer1_Phi", "Muon #phi GE21 layer 1", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMRecHit_GE21_layer2_Eta = new TH1F("MatchedGEMRecHit_GE21_layer2_Eta", "Muon #eta GE21 layer 2", n_eta_bin, eta_bin );
  MatchedGEMRecHit_GE21_layer2_Pt =  new TH1F("MatchedGEMRecHit_GE21_layer2_Pt", "Muon p_{T}  GE21 layer 2", n_pt_bin, pt_bin );
  MatchedGEMRecHit_GE21_layer2_Phi = new TH1F("MatchedGEMRecHit_GE21_layer2_Phi", "Muon #phi GE21 layer 2", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMRecHit_GE11_two_Eta = new TH1F("MatchedGEMRecHit_GE11_two_Eta", "Muon #eta GE11", n_eta_bin, eta_bin );
  MatchedGEMRecHit_GE11_two_Pt =  new TH1F("MatchedGEMRecHit_GE11_two_Pt", "Muon p_{T} GE11", n_pt_bin, pt_bin );
  MatchedGEMRecHit_GE11_two_Phi = new TH1F("MatchedGEMRecHit_GE11_two_Phi", "Muon #phi GE11", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMRecHit_GE21_two_Eta = new TH1F("MatchedGEMRecHit_GE21_two_Eta", "Muon #eta GE21", n_eta_bin, eta_bin );
  MatchedGEMRecHit_GE21_two_Pt =  new TH1F("MatchedGEMRecHit_GE21_two_Pt", "Muon p_{T} GE21", n_pt_bin, pt_bin );
  MatchedGEMRecHit_GE21_two_Phi = new TH1F("MatchedGEMRecHit_GE21_two_Phi", "Muon #phi GE21", 36, -TMath::Pi(), TMath::Pi() );
  MatchedClusteredGEMRecHit_GE11_dBunchX = new TH1I("MatchedClusteredGEMRecHit_GE11_dBunchX", "#DeltaBunchX GE11", 10, 0, 10);
  MatchedClusteredGEMRecHit_GE21_dBunchX = new TH1I("MatchedClusteredGEMRecHit_GE21_dBunchX", "#DeltaBunchX GE21", 10, 0, 10);
  MatchedGEMSegment_GE11_Eta = new TH1F("MatchedGEMSegment_GE11_Eta", "Muon #eta", n_eta_bin, eta_bin );
  MatchedGEMSegment_GE11_Pt =  new TH1F("MatchedGEMSegment_GE11_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  MatchedGEMSegment_GE11_Phi = new TH1F("MatchedGEMSegment_GE11_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMSegment_GE21_Eta = new TH1F("MatchedGEMSegment_GE21_Eta", "Muon #eta", n_eta_bin, eta_bin );
  MatchedGEMSegment_GE21_Pt =  new TH1F("MatchedGEMSegment_GE21_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  MatchedGEMSegment_GE21_Phi = new TH1F("MatchedGEMSegment_GE21_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi() );
  MatchedGEMSegment_Eta = new TH1F("MatchedGEMSegment_Eta", "Muon #eta", n_eta_bin, eta_bin );
  MatchedGEMSegment_Pt =  new TH1F("MatchedGEMSegment_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  MatchedGEMSegment_Phi = new TH1F("MatchedGEMSegment_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi() );
  TPMuon_Eta = new TH1F("TPMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  TPMuon_Pt = new TH1F("TPMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  TPMuon_Phi = new TH1F("TPMuon_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsMatchedGEMMuon_Eta = new TH1F("HitsMatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsMatchedGEMMuon_Pt  = new TH1F("HitsMatchedGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedGEMMuon_Phi = new TH1F("HitsMatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsUnmatchedGEMMuon_Eta = new TH1F("HitsUnmatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedGEMMuon_Pt  = new TH1F("HitsUnmatchedGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedGEMMuon_Phi = new TH1F("HitsUnmatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsMatchedRecoMuon_Eta = new TH1F("HitsMatchedRecoMuon_Eta", "RecoMuon #eta", n_eta_bin, eta_bin );
  HitsMatchedRecoMuon_Pt  = new TH1F("HitsMatchedRecoMuon_Pt", "RecoMuon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedRecoMuon_Phi = new TH1F("HitsMatchedRecoMuon_Phi", "RecoMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsUnmatchedRecoMuon_Eta = new TH1F("HitsUnmatchedRecoMuon_Eta", "RecoMuon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedRecoMuon_Pt  = new TH1F("HitsUnmatchedRecoMuon_Pt", "RecoMuon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedRecoMuon_Phi = new TH1F("HitsUnmatchedRecoMuon_Phi", "RecoMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsMatchedLooseMuon_Eta = new TH1F("HitsMatchedLooseMuon_Eta", "LooseMuon #eta", n_eta_bin, eta_bin );
  HitsMatchedLooseMuon_Pt  = new TH1F("HitsMatchedLooseMuon_Pt", "LooseMuon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedLooseMuon_Phi = new TH1F("HitsMatchedLooseMuon_Phi", "LooseMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsUnmatchedLooseMuon_Eta = new TH1F("HitsUnmatchedLooseMuon_Eta", "LooseMuon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedLooseMuon_Pt  = new TH1F("HitsUnmatchedLooseMuon_Pt", "LooseMuon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedLooseMuon_Phi = new TH1F("HitsUnmatchedLooseMuon_Phi", "LooseMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsMatchedMediumMuon_Eta = new TH1F("HitsMatchedMediumMuon_Eta", "MediumMuon #eta", n_eta_bin, eta_bin );
  HitsMatchedMediumMuon_Pt  = new TH1F("HitsMatchedMediumMuon_Pt", "MediumMuon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedMediumMuon_Phi = new TH1F("HitsMatchedMediumMuon_Phi", "MediumMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsUnmatchedMediumMuon_Eta = new TH1F("HitsUnmatchedMediumMuon_Eta", "MediumMuon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedMediumMuon_Pt  = new TH1F("HitsUnmatchedMediumMuon_Pt", "MediumMuon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedMediumMuon_Phi = new TH1F("HitsUnmatchedMediumMuon_Phi", "MediumMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsMatchedTightMuon_Eta = new TH1F("HitsMatchedTightMuon_Eta", "TightMuon #eta", n_eta_bin, eta_bin );
  HitsMatchedTightMuon_Pt  = new TH1F("HitsMatchedTightMuon_Pt", "TightMuon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedTightMuon_Phi = new TH1F("HitsMatchedTightMuon_Phi", "TightMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsUnmatchedTightMuon_Eta = new TH1F("HitsUnmatchedTightMuon_Eta", "TightMuon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedTightMuon_Pt  = new TH1F("HitsUnmatchedTightMuon_Pt", "TightMuon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedTightMuon_Phi = new TH1F("HitsUnmatchedTightMuon_Phi", "TightMuon #phi", 36, -TMath::Pi(), TMath::Pi() );

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
  std::vector<bool> GENMuon_matched;
  for(unsigned int i=0; i<gensize; i++) GENMuon_matched.push_back(false);

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

  //==== GEMRecHit study
  iSetup.get<MuonGeometryRecord>().get(gemGeom);
  edm::Handle<GEMRecHitCollection> gemRecHitCollection;
  iEvent.getByToken(GEMRecHit_Token, gemRecHitCollection);

  for(auto thisRecHit = gemRecHitCollection->begin(); thisRecHit != gemRecHitCollection->end(); ++thisRecHit) {
    GEMDetId id = thisRecHit->gemId();
    auto roll = gemGeom->etaPartition(id);
    auto RecHitLP = thisRecHit->localPosition();
    auto RecHitGP = roll->toGlobal(RecHitLP);
    //std::cout << "rechit station = " << id.station() << std::endl;
    //std::cout << "Rechit id = " << id << " => (" << RecHitGP.x() << ", " << RecHitGP.y() << ", " << RecHitGP.z() << ")" << std::endl;
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
  edm::Handle<GEMSegmentCollection> gemSegmentCollection;
  iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);

  for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
    GEMDetId id = gems->gemDetId();
    auto chamb = gemGeom->superChamber(id);
    auto segLP = gems->localPosition();
    auto segGP = chamb->toGlobal(segLP);
    auto rechits = gems->specificRecHits();
    //std::cout << "gemseg station = " << id.station() << std::endl;
    //std::cout << "# of rechits used in this Segment = " << rechits.size() << std::endl;
    //for(auto rechit = rechits.begin(); rechit != rechits.end(); ++rechit){
    //  std::cout << "  " << rechit->gemId() << "(" << segGP.x() << ", " << segGP.y() << ", " << segGP.z() << ")" << std::endl;
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



  //==== DeltaR matching 

  edm::LogVerbatim("GEMMuonAnalyzer") << "########### DeltaR Matching ###########";

  std::vector<int> matched_RecoMuonID;

  for(unsigned int i=0; i<gensize; i++){

    const reco::GenParticle& CurrentParticle=(*genParticles)[i];

    if ( (CurrentParticle.status()==1) && ( (CurrentParticle.pdgId()==13)  || (CurrentParticle.pdgId()==-13) ) ){   

      n_genmuon++;

      edm::LogVerbatim("GEMMuonAnalyzer") << "[this GenMuon]";
      TLorentzVector igenP4;
      igenP4.SetPtEtaPhiM(CurrentParticle.pt(), CurrentParticle.eta(), CurrentParticle.phi(), CurrentParticle.mass());

      //==== Fill Numerator

      if ( (CurrentParticle.pt() >FakeRatePtCut) ){
        GenMuon_Eta->Fill(fabs(CurrentParticle.eta()));
        if ( (fabs(CurrentParticle.eta()) > 1.6) && (fabs(CurrentParticle.eta()) < 2.4) ) {
          GenMuon_Pt->Fill(CurrentParticle.pt());
          GenMuon_Phi->Fill(CurrentParticle.phi());
        }
      }

      //==== RecoMuon

      double RecoMuon_LowestDelR = 9999, RecoMuon_NotGEMMuon_LowestDelR = 9999, RecoMuon_GEMMuon_LowestDelR = 9999;
      double RecoMuon_thisDelR = 9999, RecoMuon_NotGEMMuon_thisDelR = 9999, RecoMuon_GEMMuon_thisDelR = 9999;
      bool RecoMuon_isMatched = false, RecoMuon_NotGEMMuon_isMatched = false, RecoMuon_GEMMuon_isMatched = false;
      int RecoMuonID = 0;
      for(reco::MuonCollection::const_iterator recomuon = recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon){
        TLorentzVector LV_recomuon;
        LV_recomuon.SetPxPyPzE(recomuon->px(), recomuon->py(), recomuon->pz(), recomuon->energy());

        RecoMuon_thisDelR = igenP4.DeltaR(LV_recomuon);
        if( recomuon->pt() > FakeRatePtCut){
          if( RecoMuon_thisDelR < MatchingWindowDelR ){

            //==== matched RecoMuon
            matched_RecoMuonID.push_back(RecoMuonID);
            RecoMuon_isMatched = true;
            if( RecoMuon_thisDelR < RecoMuon_LowestDelR ){
              RecoMuon_LowestDelR = RecoMuon_thisDelR;
            }

            //==== matched RecoMuon, but NOT GEMMuon
            if( !recomuon->isGEMMuon() ){
              RecoMuon_NotGEMMuon_isMatched = true;
              if( RecoMuon_thisDelR < RecoMuon_NotGEMMuon_LowestDelR ){
                RecoMuon_NotGEMMuon_LowestDelR = RecoMuon_thisDelR;
              }
            }

            //==== matched RecoMuon and (GEMMuon && TrackerMuon)
            if( recomuon->isGEMMuon() && recomuon->isTrackerMuon() ){
              RecoMuon_GEMMuon_isMatched = true;
              if( RecoMuon_thisDelR < RecoMuon_GEMMuon_LowestDelR ){
                RecoMuon_GEMMuon_LowestDelR = RecoMuon_thisDelR;
              }
            }

          }
        }


        RecoMuonID++;

      } // END RecoMuon loop
      //==== RecoMuon matched to gen muon
      if( RecoMuon_isMatched ){
        if ((CurrentParticle.pt() >FakeRatePtCut) ){
          MatchedRecoMuon_Eta->Fill(fabs(CurrentParticle.eta()));
          if ( (TMath::Abs(CurrentParticle.eta()) > 1.6) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedRecoMuon_Pt->Fill(CurrentParticle.pt());
            MatchedRecoMuon_Phi->Fill(CurrentParticle.phi());
          }
        }
      }
      if( RecoMuon_NotGEMMuon_isMatched ){
        if ((CurrentParticle.pt() >FakeRatePtCut) ){
          MatchedRecoMuon_not_GEMMuon_Eta->Fill(fabs(CurrentParticle.eta()));
          if ( (TMath::Abs(CurrentParticle.eta()) > 1.6) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedRecoMuon_not_GEMMuon_Pt->Fill(CurrentParticle.pt());
            MatchedRecoMuon_not_GEMMuon_Phi->Fill(CurrentParticle.phi());
          }
        }
      }
      if( RecoMuon_GEMMuon_isMatched ){
        GENMuon_matched.at(i) = true;
        if ((CurrentParticle.pt() >FakeRatePtCut) ){
          MatchedGEMMuon_Eta->Fill(fabs(CurrentParticle.eta()));
          if ( (TMath::Abs(CurrentParticle.eta()) > 1.6) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedGEMMuon_Pt->Fill(CurrentParticle.pt());
            MatchedGEMMuon_Phi->Fill(CurrentParticle.phi());
          }
        }
      }

      //==== GEMRecHit

      edm::LogVerbatim("GEMMuonAnalyzer") << "==GEMRecHit Loop==";
      double GEMRecHit_LowestDelR = 9999;
      double GEMRecHit_thisDelR = 9999;
      bool GEMRecHit_isMatched = false, isGE11layer1Matched = false, isGE11layer2Matched = false, isGE21layer1Matched = false, isGE21layer2Matched = false;
      int GEMRecHitID = 0;
      bool GE11doublematched = false, GE21doublematched = false;
      std::vector<GEMRecHit> matched_rh;
      for (auto thisRecHit = gemRecHitCollection->begin(); thisRecHit != gemRecHitCollection->end(); ++thisRecHit) {
        GEMDetId id = thisRecHit->gemId();
        auto roll = gemGeom->etaPartition(id);
        auto RecHitLP = thisRecHit->localPosition();
        auto RecHitGP = roll->toGlobal(RecHitLP);
        TLorentzVector RecHit;
        //RecHit.SetPtEtaPhiM(1, RecHitGP.eta(), RecHitGP.phi(), 0);
        RecHit.SetPxPyPzE(RecHitGP.x(), RecHitGP.y(), RecHitGP.z(), 1);
        GEMRecHit_thisDelR = igenP4.DeltaR(RecHit);
        edm::LogVerbatim("GEMMuonAnalyzer") << "this RecHit id = " << id << " => deltaR = " << GEMRecHit_thisDelR;
        if( GEMRecHit_thisDelR < MatchingWindowDelR ){
          matched_rh.push_back(*thisRecHit);
          GEMRecHit_isMatched = true;
          if(id.station() == 1 && id.layer() == 1) isGE11layer1Matched = true;
          if(id.station() == 1 && id.layer() == 2) isGE11layer2Matched = true;
          if(id.station() == 3 && id.layer() == 1) isGE21layer1Matched = true;
          if(id.station() == 3 && id.layer() == 2) isGE21layer2Matched = true;
          edm::LogVerbatim("GEMMuonAnalyzer") << "-> matched";
          if( GEMRecHit_thisDelR < GEMRecHit_LowestDelR ){
            GEMRecHit_LowestDelR = GEMRecHit_thisDelR;
          }
        }

        GEMRecHitID++;

      } // END GEMRecHit loop

      //==== Cluster rechits who has same SuperChamber ID
      std::map< GEMDetId, std::vector<GEMRecHit> > map_scid_to_rh;
      for(auto idit = matched_rh.begin(); idit != matched_rh.end(); ++idit){
        GEMDetId scid = idit->gemId().superChamberId();
        map_scid_to_rh[scid].push_back(*idit);
      }

      edm::LogVerbatim("GEMMuonAnalyzer") << std::endl << "RecHits are now clustered";
      //==== check both layer1 and layer2 are used
      for(auto mapit = map_scid_to_rh.begin(); mapit != map_scid_to_rh.end(); ++mapit){

        edm::LogVerbatim("GEMMuonAnalyzer") << "[this SC]";
        std::vector<GEMRecHit> thisRHs_from_SC = mapit->second;
        std::map< int, std::vector<GEMRecHit> > map_roll_to_rh;
        //==== Among the rechits who has same SuperChamber ID,
        //==== now cluster them who has same roll ID
        for(auto itit = thisRHs_from_SC.begin(); itit != thisRHs_from_SC.end(); ++itit){
          map_roll_to_rh[itit->gemId().roll()].push_back(*itit);
        }

        for(auto itit = map_roll_to_rh.begin(); itit != map_roll_to_rh.end(); ++itit){
          edm::LogVerbatim("GEMMuonAnalyzer") << "  [this roll]";
          std::vector<GEMRecHit> thisRHs_from_roll = itit->second;
          bool GE11layer1_fired(false), GE11layer2_fired(false), GE21layer1_fired(false), GE21layer2_fired(false);
          for(auto ititit = thisRHs_from_roll.begin(); ititit != thisRHs_from_roll.end(); ++ititit){
            if(ititit->gemId().station() == 1 && ititit->gemId().layer() == 1) GE11layer1_fired = true;
            if(ititit->gemId().station() == 1 && ititit->gemId().layer() == 2) GE11layer2_fired = true;
            if(ititit->gemId().station() == 3 && ititit->gemId().layer() == 1) GE21layer1_fired = true;
            if(ititit->gemId().station() == 3 && ititit->gemId().layer() == 2) GE21layer2_fired = true;
            edm::LogVerbatim("GEMMuonAnalyzer") << "  " << ititit->gemId();
          }
          bool GE11_segment_exist = false, GE21_segment_exist = false;
          if(GE11layer1_fired && GE11layer2_fired){
            GE11doublematched = true;
            edm::LogVerbatim("GEMMuonAnalyzer") << "  ==> GE11 double matched";
            MatchedClusteredGEMRecHit_GE11_dBunchX->Fill( abs(thisRHs_from_roll.at(0).BunchX() - thisRHs_from_roll.at(1).BunchX())  );
            //==== find if GEMSegment is made
            for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
              auto rechits = gems->specificRecHits();
              for(auto itrh = rechits.begin(); itrh != rechits.end(); ++itrh){
                GEMDetId rhid = itrh->gemId();
                if(rhid == thisRHs_from_roll.at(0).gemId()){
                  edm::LogVerbatim("GEMMuonAnalyzer") << "  ==> Segment exists ("<<thisRHs_from_roll.size()<<"): " << rhid;
                  GE11_segment_exist = true;
                  break;
                }
              }
            }
            if(!GE11_segment_exist) edm::LogVerbatim("GEMMuonAnalyzer") << "  ==> Segment DO NOT exist ("<<thisRHs_from_roll.size()<<")";
          }
          if(GE21layer1_fired && GE21layer2_fired){
            GE21doublematched = true;
            edm::LogVerbatim("GEMMuonAnalyzer") << "  ==> GE21 double matched";
            MatchedClusteredGEMRecHit_GE21_dBunchX->Fill( abs(thisRHs_from_roll.at(0).BunchX() - thisRHs_from_roll.at(1).BunchX())  );
            //==== find if GEMSegment is made
            for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
              auto rechits = gems->specificRecHits();
              for(auto itrh = rechits.begin(); itrh != rechits.end(); ++itrh){
                GEMDetId rhid = itrh->gemId();
                if(rhid == thisRHs_from_roll.at(0).gemId()){
                  edm::LogVerbatim("GEMMuonAnalyzer") << "  ==> Segment exists ("<<thisRHs_from_roll.size()<<"): " << rhid;
                  GE21_segment_exist = true;
                  break;
                }
              }
            }
            if(!GE21_segment_exist) edm::LogVerbatim("GEMMuonAnalyzer") << "  ==> Segment DO NOT exist ("<<thisRHs_from_roll.size()<<")";
          }

        } // END Roll Loop


      } // END SuperChamber Id loop

      //==== GEMRecHit matched to gen muon
      if( GEMRecHit_isMatched ){
        if ((CurrentParticle.pt() >FakeRatePtCut) ){

          //==== Eta
          MatchedGEMRecHit_Eta->Fill(fabs(CurrentParticle.eta()));
          if(isGE11layer1Matched) MatchedGEMRecHit_GE11_layer1_Eta->Fill(fabs(CurrentParticle.eta()));
          if(isGE11layer2Matched) MatchedGEMRecHit_GE11_layer2_Eta->Fill(fabs(CurrentParticle.eta()));
          if(isGE21layer1Matched) MatchedGEMRecHit_GE21_layer1_Eta->Fill(fabs(CurrentParticle.eta()));
          if(isGE21layer2Matched) MatchedGEMRecHit_GE21_layer2_Eta->Fill(fabs(CurrentParticle.eta()));
          if(GE11doublematched) MatchedGEMRecHit_GE11_two_Eta->Fill(fabs(CurrentParticle.eta()));
          if(GE21doublematched) MatchedGEMRecHit_GE21_two_Eta->Fill(fabs(CurrentParticle.eta()));

          //==== pt, phi
          if ( (TMath::Abs(CurrentParticle.eta()) > 1.6) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedGEMRecHit_Pt->Fill(CurrentParticle.pt());
            if(isGE11layer1Matched) MatchedGEMRecHit_GE11_layer1_Pt->Fill(CurrentParticle.pt());
            if(isGE11layer2Matched) MatchedGEMRecHit_GE11_layer2_Pt->Fill(CurrentParticle.pt());
            if(isGE21layer1Matched) MatchedGEMRecHit_GE21_layer1_Pt->Fill(CurrentParticle.pt());
            if(isGE21layer2Matched) MatchedGEMRecHit_GE21_layer2_Pt->Fill(CurrentParticle.pt());
            if(GE11doublematched) MatchedGEMRecHit_GE11_two_Pt->Fill(CurrentParticle.pt());
            if(GE21doublematched) MatchedGEMRecHit_GE21_two_Pt->Fill(CurrentParticle.pt());
            MatchedGEMRecHit_Phi->Fill(CurrentParticle.phi());
            if(isGE11layer1Matched) MatchedGEMRecHit_GE11_layer1_Phi->Fill(CurrentParticle.phi());
            if(isGE11layer2Matched) MatchedGEMRecHit_GE11_layer2_Phi->Fill(CurrentParticle.phi());
            if(isGE21layer1Matched) MatchedGEMRecHit_GE21_layer1_Phi->Fill(CurrentParticle.phi());
            if(isGE21layer2Matched) MatchedGEMRecHit_GE21_layer2_Phi->Fill(CurrentParticle.phi());
            if(GE11doublematched) MatchedGEMRecHit_GE11_two_Phi->Fill(CurrentParticle.phi());
            if(GE21doublematched) MatchedGEMRecHit_GE21_two_Phi->Fill(CurrentParticle.phi());
          }

        }
      }


      //==== GEMSegment

      edm::LogVerbatim("GEMMuonAnalyzer") << std::endl << "==GEMSegment Loop==";
      double GEMSegment_LowestDelR = 9999;
      double GEMSegment_thisDelR = 9999;
      bool GEMSegment_isMatched = false, isGE11Matched = false, isGE21Matched = false;
      int GEMSegID = 0;
      for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
        GEMDetId id = gems->gemDetId();
        auto chamb = gemGeom->superChamber(id);
        auto segLP = gems->localPosition();
        auto segGP = chamb->toGlobal(segLP);
        edm::LogVerbatim("GEMMuonAnalyzer") << "[This Segment is made of..]";
        auto rechits = gems->specificRecHits();
        for(auto itrh = rechits.begin(); itrh != rechits.end(); ++itrh){
          GEMDetId rhid = itrh->gemId();
          TLorentzVector tmp_segment;
          tmp_segment.SetPxPyPzE(segGP.x(), segGP.y(), segGP.z(), 1);
          edm::LogVerbatim("GEMMuonAnalyzer") << rhid << "=> deltaR = " << igenP4.DeltaR(tmp_segment);
        }

        TLorentzVector segment;
        //segment.SetPtEtaPhiM(1, segGP.eta(), segGP.phi(), 0);
        segment.SetPxPyPzE(segGP.x(), segGP.y(), segGP.z(), 1);
        GEMSegment_thisDelR = igenP4.DeltaR(segment);
        if( GEMSegment_thisDelR < MatchingWindowDelR ){
          GEMSegment_isMatched = true;
          edm::LogVerbatim("GEMMuonAnalyzer") << "-> matched";
          if(id.station() == 1) isGE11Matched = true;
          if(id.station() == 3) isGE21Matched = true;
          if( GEMSegment_thisDelR < GEMSegment_LowestDelR ){
            GEMSegment_LowestDelR = GEMSegment_thisDelR;
          }
        }

        GEMSegID++;

      }

      //==== GEMSegment matched to gen muon
      if( GEMSegment_isMatched ){

        if ((CurrentParticle.pt() >FakeRatePtCut) ){

          MatchedGEMSegment_Eta->Fill(fabs(CurrentParticle.eta()));
          if(isGE11Matched) MatchedGEMSegment_GE11_Eta->Fill(fabs(CurrentParticle.eta()));
          if(isGE21Matched) MatchedGEMSegment_GE21_Eta->Fill(fabs(CurrentParticle.eta()));

          if ( (TMath::Abs(CurrentParticle.eta()) > 1.6) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedGEMSegment_Pt->Fill(CurrentParticle.pt());
            if(isGE11Matched) MatchedGEMSegment_GE11_Pt->Fill(CurrentParticle.pt());
            if(isGE21Matched) MatchedGEMSegment_GE21_Pt->Fill(CurrentParticle.pt());
            MatchedGEMSegment_Phi->Fill(CurrentParticle.phi());
            if(isGE11Matched) MatchedGEMSegment_GE11_Phi->Fill(CurrentParticle.phi());
            if(isGE21Matched) MatchedGEMSegment_GE21_Phi->Fill(CurrentParticle.phi());
          }

        }
      }

      //==== RecoMuon matched, but NOT (GEMMuon && both stations segment)
      if( RecoMuon_NotGEMMuon_isMatched && !isGE11Matched && !isGE21Matched ){
        if ((CurrentParticle.pt() >FakeRatePtCut) ){
          MatchedRecoMuon_not_GEMMuon_no_gemseg_Eta->Fill(fabs(CurrentParticle.eta()));
          if ( (TMath::Abs(CurrentParticle.eta()) > 1.6) && (TMath::Abs(CurrentParticle.eta()) < 2.4) )  {
            MatchedRecoMuon_not_GEMMuon_no_gemseg_Pt->Fill(CurrentParticle.pt());
            MatchedRecoMuon_not_GEMMuon_no_gemseg_Phi->Fill(CurrentParticle.phi());
          }
        }
      }

    } // END prompt muons selection
  } // END gen particle loop
  
  
  if (UseAssociators) {

    //std::cout << "TrackingParticle size = " << trackingParticles->size() << std::endl;
    for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++){
      TrackingParticleRef tpr(trackingParticles, i);
      TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
      TrackingParticle::Vector momentumTP;
      TrackingParticle::Point vertexTP;

      if (abs(tp->pdgId()) != 13) continue;

      bool Eta_1p6_2p4 = fabs(tp->eta()) > 1.6 && fabs(tp->eta()) < 2.4;
      bool Pt_5 = tp->pt() > 5;
      if( Eta_1p6_2p4 && Pt_5 ){
        bool SignalMuon = false;
        if(tp->status() != -99){
          //==== Pythia8 gen status : home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
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
          TPMuon_Eta->Fill(fabs(tp->eta()));
          TPMuon_Pt->Fill(tp->pt());
          TPMuon_Phi->Fill(tp->phi());
          //==== Looking at SimToReco
          for (unsigned int www=0;www<label.size();www++){

            Handle<reco::SimToRecoCollection > simtorecoCollectionH;
            iEvent.getByLabel(associators[www],simtorecoCollectionH);
            reco::SimToRecoCollection simRecColl= *(simtorecoCollectionH.product());

            Handle<reco::RecoToSimCollection > recotosimCollectionH;
            iEvent.getByLabel(associators[www],recotosimCollectionH);
            reco::RecoToSimCollection recSimColl= *(recotosimCollectionH.product());

            edm::Handle<View<Track> >  trackCollection;
            iEvent.getByToken(track_Collection_Token[www], trackCollection);

            if( (simRecColl.find(tpr) == simRecColl.end()) || (simRecColl[tpr].size() == 0) ){
              edm::LogVerbatim("GEMMuonAnalyzer") << "No SimToReco found for this TrackingParticle";
            }
            else{
              std::vector<std::pair<RefToBase<Track>, double> > rt = simRecColl[tpr];
              RefToBase<Track> rtr = rt.begin()->first;
              std::cout << "This SimToReco :" << std::endl;
              for(std::vector<std::pair<RefToBase<Track>, double> >::const_iterator itit=rt.begin(); itit!=rt.end(); itit++){
                std::cout << "  quality = "<<itit->second<<std::endl;
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
                    n_AssoByHits_matched_GEMmuon++;
                    HitsMatchedGEMMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedGEMMuon_Pt->Fill(tpr->pt());
                    HitsMatchedGEMMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="recoMuonSel"){
                    HitsMatchedRecoMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedRecoMuon_Pt->Fill(tpr->pt());
                    HitsMatchedRecoMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="looseMuonSel"){
                    HitsMatchedLooseMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedLooseMuon_Pt->Fill(tpr->pt());
                    HitsMatchedLooseMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="mediumMuonSel"){
                    HitsMatchedMediumMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedMediumMuon_Pt->Fill(tpr->pt());
                    HitsMatchedMediumMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="tightMuonSel"){
                    HitsMatchedTightMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedTightMuon_Pt->Fill(tpr->pt());
                    HitsMatchedTightMuon_Phi->Fill(tpr->phi());
                  }

                }
              }
            }

          } //==== END of label loop

        } //==== END if(SignalMuon)

      } // END if( Eta_1p6_2p4 && Pt_5 )


    } // END for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++)

    for (unsigned int www=0;www<label.size();www++){

      reco::RecoToSimCollection recSimColl;
      reco::SimToRecoCollection simRecColl;
      edm::Handle<View<Track> >  trackCollection;

      Handle<reco::SimToRecoCollection > simtorecoCollectionH;
      iEvent.getByLabel(associators[www],simtorecoCollectionH);
      simRecColl= *(simtorecoCollectionH.product());

      Handle<reco::RecoToSimCollection > recotosimCollectionH;
      iEvent.getByLabel(associators[www],recotosimCollectionH);
      recSimColl= *(recotosimCollectionH.product());

      unsigned int trackCollectionSize = 0;
      iEvent.getByToken(track_Collection_Token[www], trackCollection);
      trackCollectionSize = trackCollection->size();

      //==== loop over (GEMMuon/RecoMuon/...) tracks
      //std::cout << "Obj = " << label[www] << ", trackCollectionSize = " << trackCollectionSize << std::endl;
      for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
        //std::cout << i << "th trackCollection iterator" << std::endl;
        RefToBase<Track> track(trackCollection, i);

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
            bool Eta_1p6_2p4 = fabs(tpr->eta()) > 1.6 && fabs(tpr->eta()) < 2.4;
            bool Pt_5 = tpr->pt() > 5;
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
            //if( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) && Eta_1p6_2p4 && Pt_5 && SignalMuon ) {
            if( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) && SignalMuon ) {
              edm::LogVerbatim("GEMMuonAnalyzer") << "Found matched TrackingParticle";
            }
            else{
              isFake = true;
            }
          }
        }

        if(isFake) {
          if(label[www]=="gemMuonSel"){
             HitsUnmatchedGEMMuon_Eta->Fill(fabs(track->eta()));
             HitsUnmatchedGEMMuon_Pt->Fill(track->pt());
             HitsUnmatchedGEMMuon_Phi->Fill(track->phi());
          }
          if(label[www]=="recoMuonSel"){
            HitsUnmatchedRecoMuon_Eta->Fill(fabs(track->eta()));
            HitsUnmatchedRecoMuon_Pt->Fill(track->pt());
            HitsUnmatchedRecoMuon_Phi->Fill(track->phi());
          }
          if(label[www]=="looseMuonSel"){
            HitsUnmatchedLooseMuon_Eta->Fill(fabs(track->eta()));
            HitsUnmatchedLooseMuon_Pt->Fill(track->pt());
            HitsUnmatchedLooseMuon_Phi->Fill(track->phi());
          }
          if(label[www]=="mediumMuonSel"){
            HitsUnmatchedMediumMuon_Eta->Fill(fabs(track->eta()));
            HitsUnmatchedMediumMuon_Pt->Fill(track->pt());
            HitsUnmatchedMediumMuon_Phi->Fill(track->phi());
          }
          if(label[www]=="tightMuonSel"){
            HitsUnmatchedTightMuon_Eta->Fill(fabs(track->eta()));
            HitsUnmatchedTightMuon_Pt->Fill(track->pt());
            HitsUnmatchedTightMuon_Phi->Fill(track->phi());
          }
        } //==== END if(isFake)

      }//==== END for(View<Track>::size_type i=0; i<trackCollectionSize; ++i)



    } // END muon obj loop (GEMMuon, RecoMuon..)
  } // END if (UseAssociators)






  //==== Matching Study

  if(doMatchingStudy){

    edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "########### GEMMuon Matching Study ###########";

    for(unsigned int i=0; i<gensize; i++){
      if( !GENMuon_matched.at(i) ){
        const reco::GenParticle& CurrentParticle=(*genParticles)[i];
        edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "Matching Failed GENMuon : (pt, eta, phi) = ("<<CurrentParticle.pt()<<","<<CurrentParticle.eta()<<","<<CurrentParticle.phi()<<")";
      }
    }
    edm::LogVerbatim("GEMMuonAnalyzer_Matching");

    Long64_t aaa=0;
    for(reco::MuonCollection::const_iterator recomuon=recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon, aaa++) {

      for(auto chmatch = recomuon->matches().begin(); chmatch != recomuon->matches().end(); ++chmatch){
        if( chmatch->id.subdetId() != 4 ) continue;
        if( chmatch->segmentMatches.size() == 0 ) continue;

        for(auto segmatch = chmatch->segmentMatches.begin(); segmatch != chmatch->segmentMatches.end(); ++segmatch){

          int station = segmatch->gemSegmentRef->specificRecHits()[0].gemId().station();

          Double_t sigmax = sqrt( chmatch->xErr*chmatch->xErr + segmatch->xErr*segmatch->xErr );
          Double_t sigmay = sqrt( chmatch->yErr*chmatch->yErr + segmatch->yErr*segmatch->yErr );

          GlobalVector Dir_ch(chmatch->dXdZ, chmatch->dYdZ, 1);
          GlobalVector Dir_seg(segmatch->dXdZ, segmatch->dYdZ, 1);

          Double_t DelX = std::abs(chmatch->x - segmatch->x);
          Double_t DelX_over_sigma = DelX/sigmax;
          Double_t DelY = std::abs(chmatch->y - segmatch->y);
          Double_t DelY_over_sigma = DelY/sigmay;
          Double_t DotDir = Dir_ch.unit().dot( Dir_seg.unit() );

          if( std::find( matched_RecoMuonID.begin(), matched_RecoMuonID.end(), aaa ) != matched_RecoMuonID.end() ) {
            GEMSegmentRef thisGEMSegRef = segmatch->gemSegmentRef;
            bool XMatched(false), YMatched(false), DirMatched(false);
            if(station == 1){
              XMatched = DelX < Current_maxDiffXGE11 || DelX_over_sigma < Current_trackerGEM_maxPull;
              YMatched = DelY < Current_maxDiffYGE11 || DelY_over_sigma < Current_trackerGEM_maxPull;
            }
            if(station == 3){
              XMatched = DelX < Current_maxDiffXGE21 || DelX_over_sigma < Current_trackerGEM_maxPull;
              YMatched = DelY < Current_maxDiffYGE21 || DelY_over_sigma < Current_trackerGEM_maxPull;
            }
            DirMatched = DotDir > Current_minDotDir;

            edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "this GEMSegment id = " << thisGEMSegRef->gemDetId();
            edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  DelX = " << DelX;
            edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  DelX/sigma = " << DelX_over_sigma;
            if(XMatched) edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  -> XMatch";
            edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  DelY = " << DelY;
            edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  DelY/sigma = " << DelY_over_sigma;
            if(YMatched) edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  -> YMatch";
            edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  DotDir = " << DotDir;
            if(DirMatched) edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "  -> DirMatched";
            if(XMatched && YMatched && DirMatched) edm::LogVerbatim("GEMMuonAnalyzer_Matching") << "==> GEMMuon Matched";

          }

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


        } // END segment match loop

      } // END MatchedChamber loop

    } // END RecoMuon loop  


  } // END domMatchingStudy

 

  if(doMatchingStudy){

    Handle <TrackCollection > generalTracks;
    iEvent.getByToken (generalTracksToken_, generalTracks);

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

  //==== quick eff print
  edm::LogVerbatim("GEMMuonAnalyzer");
  edm::LogVerbatim("GEMMuonAnalyzer") << "=================================";
  edm::LogVerbatim("GEMMuonAnalyzer") << "# of gen muon = " << n_genmuon;
  edm::LogVerbatim("GEMMuonAnalyzer") << "# of dR matched GEM muon = " << n_dR_matched_GEMmuon;
  edm::LogVerbatim("GEMMuonAnalyzer") << "# of AssoByHits matched GEM muon = " << n_AssoByHits_matched_GEMmuon;

  histoFile->cd();

  Nevents_h->Write();

  /* gen-reco delta R matching */
  GenMuon_Eta->Write();
  GenMuon_Pt->Write();
  GenMuon_Phi->Write();
  MatchedRecoMuon_Eta->Write();
  MatchedRecoMuon_Pt->Write();
  MatchedRecoMuon_Phi->Write();
  MatchedRecoMuon_not_GEMMuon_Eta->Write();
  MatchedRecoMuon_not_GEMMuon_Pt->Write();
  MatchedRecoMuon_not_GEMMuon_Phi->Write();
  MatchedRecoMuon_not_GEMMuon_no_gemseg_Eta->Write();
  MatchedRecoMuon_not_GEMMuon_no_gemseg_Pt->Write();
  MatchedRecoMuon_not_GEMMuon_no_gemseg_Phi->Write();
  MatchedGEMMuon_Eta->Write();
  MatchedGEMMuon_Pt->Write();
  MatchedGEMMuon_Phi->Write();
  MatchedGEMRecHit_Eta->Write();
  MatchedGEMRecHit_Pt->Write();
  MatchedGEMRecHit_Phi->Write();
  MatchedGEMRecHit_GE11_layer1_Eta->Write();
  MatchedGEMRecHit_GE11_layer1_Pt->Write();
  MatchedGEMRecHit_GE11_layer1_Phi->Write();
  MatchedGEMRecHit_GE11_layer2_Eta->Write();
  MatchedGEMRecHit_GE11_layer2_Pt->Write();
  MatchedGEMRecHit_GE11_layer2_Phi->Write();
  MatchedGEMRecHit_GE21_layer1_Eta->Write();
  MatchedGEMRecHit_GE21_layer1_Pt->Write();
  MatchedGEMRecHit_GE21_layer1_Phi->Write();
  MatchedGEMRecHit_GE21_layer2_Eta->Write();
  MatchedGEMRecHit_GE21_layer2_Pt->Write();
  MatchedGEMRecHit_GE21_layer2_Phi->Write();
  MatchedGEMRecHit_GE11_two_Eta->Write();
  MatchedGEMRecHit_GE11_two_Pt->Write();
  MatchedGEMRecHit_GE11_two_Phi->Write();
  MatchedGEMRecHit_GE21_two_Eta->Write();
  MatchedGEMRecHit_GE21_two_Pt->Write();
  MatchedGEMRecHit_GE21_two_Phi->Write();
  MatchedClusteredGEMRecHit_GE11_dBunchX->Write();
  MatchedClusteredGEMRecHit_GE21_dBunchX->Write();
  MatchedGEMSegment_GE11_Eta->Write();
  MatchedGEMSegment_GE11_Pt->Write();
  MatchedGEMSegment_GE11_Phi->Write();
  MatchedGEMSegment_GE21_Eta->Write();
  MatchedGEMSegment_GE21_Pt->Write();
  MatchedGEMSegment_GE21_Phi->Write();
  MatchedGEMSegment_Eta->Write();
  MatchedGEMSegment_Pt->Write();
  MatchedGEMSegment_Phi->Write();
  //==== RecoMuon Efficiency
  TEfficiency* Eff_RecoMuon_Eta = new TEfficiency(*MatchedRecoMuon_Eta, *GenMuon_Eta);
  TEfficiency* Eff_RecoMuon_Pt = new TEfficiency(*MatchedRecoMuon_Pt, *GenMuon_Pt);
  TEfficiency* Eff_RecoMuon_Phi = new TEfficiency(*MatchedRecoMuon_Phi, *GenMuon_Phi);
  Eff_RecoMuon_Eta->SetName("Eff_RecoMuon_Eta");
  Eff_RecoMuon_Pt->SetName("Eff_RecoMuon_Pt");
  Eff_RecoMuon_Phi->SetName("Eff_RecoMuon_Phi");
  Eff_RecoMuon_Eta->Write();
  Eff_RecoMuon_Pt->Write();
  Eff_RecoMuon_Phi->Write();
  //==== NotGEMMuon Efficiency
  TEfficiency* Eff_NotGEMMuon_Eta = new TEfficiency(*MatchedRecoMuon_not_GEMMuon_Eta, *GenMuon_Eta);
  TEfficiency* Eff_NotGEMMuon_Pt = new TEfficiency(*MatchedRecoMuon_not_GEMMuon_Pt, *GenMuon_Pt);
  TEfficiency* Eff_NotGEMMuon_Phi = new TEfficiency(*MatchedRecoMuon_not_GEMMuon_Phi, *GenMuon_Phi);
  Eff_NotGEMMuon_Eta->SetName("Eff_NotGEMMuon_Eta");
  Eff_NotGEMMuon_Pt->SetName("Eff_NotGEMMuon_Pt");
  Eff_NotGEMMuon_Phi->SetName("Eff_NotGEMMuon_Phi");
  Eff_NotGEMMuon_Eta->Write();
  Eff_NotGEMMuon_Pt->Write();
  Eff_NotGEMMuon_Phi->Write();
  //==== NOTGEMMuon, no GEMSegment matched Efficiency
  TEfficiency* Eff_NotGEMMuon_no_gemseg_Eta = new TEfficiency(*MatchedRecoMuon_not_GEMMuon_no_gemseg_Eta, *GenMuon_Eta);
  TEfficiency* Eff_NotGEMMuon_no_gemseg_Pt = new TEfficiency(*MatchedRecoMuon_not_GEMMuon_no_gemseg_Pt, *GenMuon_Pt);
  TEfficiency* Eff_NotGEMMuon_no_gemseg_Phi = new TEfficiency(*MatchedRecoMuon_not_GEMMuon_no_gemseg_Phi, *GenMuon_Phi);
  Eff_NotGEMMuon_no_gemseg_Eta->SetName("Eff_NotGEMMuon_no_gemseg_Eta");
  Eff_NotGEMMuon_no_gemseg_Pt->SetName("Eff_NotGEMMuon_no_gemseg_Pt");
  Eff_NotGEMMuon_no_gemseg_Phi->SetName("Eff_NotGEMMuon_no_gemseg_Phi");
  Eff_NotGEMMuon_no_gemseg_Eta->Write();
  Eff_NotGEMMuon_no_gemseg_Pt->Write();
  Eff_NotGEMMuon_no_gemseg_Phi->Write();
  //==== GEMMuon Efficiency
  TEfficiency* Eff_GEMMuon_Eta = new TEfficiency(*MatchedGEMMuon_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMMuon_Pt = new TEfficiency(*MatchedGEMMuon_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMMuon_Phi = new TEfficiency(*MatchedGEMMuon_Phi, *GenMuon_Phi);
  Eff_GEMMuon_Eta->SetName("Eff_GEMMuon_Eta");
  Eff_GEMMuon_Pt->SetName("Eff_GEMMuon_Pt");
  Eff_GEMMuon_Phi->SetName("Eff_GEMMuon_Phi");
  Eff_GEMMuon_Eta->Write();
  Eff_GEMMuon_Pt->Write();
  Eff_GEMMuon_Phi->Write();
  //==== GEMMuon/RecoMuon Efficiency
  TEfficiency* Eff_GEMMuon_over_RecoMuon_Eta = new TEfficiency(*MatchedGEMMuon_Eta, *MatchedRecoMuon_Eta);
  TEfficiency* Eff_GEMMuon_over_RecoMuon_Pt = new TEfficiency(*MatchedGEMMuon_Pt, *MatchedRecoMuon_Pt);
  TEfficiency* Eff_GEMMuon_over_RecoMuon_Phi = new TEfficiency(*MatchedGEMMuon_Phi, *MatchedRecoMuon_Phi);
  Eff_GEMMuon_over_RecoMuon_Eta->SetName("Eff_GEMMuon_over_RecoMuon_Eta");
  Eff_GEMMuon_over_RecoMuon_Pt->SetName("Eff_GEMMuon_over_RecoMuon_Pt");
  Eff_GEMMuon_over_RecoMuon_Phi->SetName("Eff_GEMMuon_over_RecoMuon_Phi");
  Eff_GEMMuon_over_RecoMuon_Eta->Write();
  Eff_GEMMuon_over_RecoMuon_Pt->Write();
  Eff_GEMMuon_over_RecoMuon_Phi->Write();
  //==== Rechit Efficiency
  TEfficiency* Eff_GEMRecHit_Eta = new TEfficiency(*MatchedGEMRecHit_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMRecHit_Pt = new TEfficiency(*MatchedGEMRecHit_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMRecHit_Phi = new TEfficiency(*MatchedGEMRecHit_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMRecHit_GE11_layer1_Eta = new TEfficiency(*MatchedGEMRecHit_GE11_layer1_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMRecHit_GE11_layer1_Pt = new TEfficiency(*MatchedGEMRecHit_GE11_layer1_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMRecHit_GE11_layer1_Phi = new TEfficiency(*MatchedGEMRecHit_GE11_layer1_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMRecHit_GE11_layer2_Eta = new TEfficiency(*MatchedGEMRecHit_GE11_layer2_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMRecHit_GE11_layer2_Pt = new TEfficiency(*MatchedGEMRecHit_GE11_layer2_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMRecHit_GE11_layer2_Phi = new TEfficiency(*MatchedGEMRecHit_GE11_layer2_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMRecHit_GE21_layer1_Eta = new TEfficiency(*MatchedGEMRecHit_GE21_layer1_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMRecHit_GE21_layer1_Pt = new TEfficiency(*MatchedGEMRecHit_GE21_layer1_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMRecHit_GE21_layer1_Phi = new TEfficiency(*MatchedGEMRecHit_GE21_layer1_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMRecHit_GE21_layer2_Eta = new TEfficiency(*MatchedGEMRecHit_GE21_layer2_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMRecHit_GE21_layer2_Pt = new TEfficiency(*MatchedGEMRecHit_GE21_layer2_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMRecHit_GE21_layer2_Phi = new TEfficiency(*MatchedGEMRecHit_GE21_layer2_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMRecHit_GE11_two_Eta = new TEfficiency(*MatchedGEMRecHit_GE11_two_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMRecHit_GE11_two_Pt = new TEfficiency(*MatchedGEMRecHit_GE11_two_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMRecHit_GE11_two_Phi = new TEfficiency(*MatchedGEMRecHit_GE11_two_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMRecHit_GE21_two_Eta = new TEfficiency(*MatchedGEMRecHit_GE21_two_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMRecHit_GE21_two_Pt = new TEfficiency(*MatchedGEMRecHit_GE21_two_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMRecHit_GE21_two_Phi = new TEfficiency(*MatchedGEMRecHit_GE21_two_Phi, *GenMuon_Phi);
  Eff_GEMRecHit_Eta->SetName("Eff_GEMRecHit_Eta");
  Eff_GEMRecHit_Pt->SetName("Eff_GEMRecHit_Pt");
  Eff_GEMRecHit_Phi->SetName("Eff_GEMRecHit_Phi");
  Eff_GEMRecHit_GE11_layer1_Eta->SetName("Eff_GEMRecHit_GE11_layer1_Eta");
  Eff_GEMRecHit_GE11_layer1_Pt->SetName("Eff_GEMRecHit_GE11_layer1_Pt");
  Eff_GEMRecHit_GE11_layer1_Phi->SetName("Eff_GEMRecHit_GE11_layer1_Phi");
  Eff_GEMRecHit_GE11_layer2_Eta->SetName("Eff_GEMRecHit_GE11_layer2_Eta");
  Eff_GEMRecHit_GE11_layer2_Pt->SetName("Eff_GEMRecHit_GE11_layer2_Pt");
  Eff_GEMRecHit_GE11_layer2_Phi->SetName("Eff_GEMRecHit_GE11_layer2_Phi");
  Eff_GEMRecHit_GE21_layer1_Eta->SetName("Eff_GEMRecHit_GE21_layer1_Eta");
  Eff_GEMRecHit_GE21_layer1_Pt->SetName("Eff_GEMRecHit_GE21_layer1_Pt");
  Eff_GEMRecHit_GE21_layer1_Phi->SetName("Eff_GEMRecHit_GE21_layer1_Phi");
  Eff_GEMRecHit_GE21_layer2_Eta->SetName("Eff_GEMRecHit_GE21_layer2_Eta");
  Eff_GEMRecHit_GE21_layer2_Pt->SetName("Eff_GEMRecHit_GE21_layer2_Pt");
  Eff_GEMRecHit_GE21_layer2_Phi->SetName("Eff_GEMRecHit_GE21_layer2_Phi");
  Eff_GEMRecHit_GE11_two_Eta->SetName("Eff_GEMRecHit_GE11_two_Eta");
  Eff_GEMRecHit_GE11_two_Pt->SetName("Eff_GEMRecHit_GE11_two_Pt");
  Eff_GEMRecHit_GE11_two_Phi->SetName("Eff_GEMRecHit_GE11_two_Phi");
  Eff_GEMRecHit_GE21_two_Eta->SetName("Eff_GEMRecHit_GE21_two_Eta");
  Eff_GEMRecHit_GE21_two_Pt->SetName("Eff_GEMRecHit_GE21_two_Pt");
  Eff_GEMRecHit_GE21_two_Phi->SetName("Eff_GEMRecHit_GE21_two_Phi");
  Eff_GEMRecHit_Eta->Write();
  Eff_GEMRecHit_Pt->Write();
  Eff_GEMRecHit_Phi->Write();
  Eff_GEMRecHit_GE11_layer1_Eta->Write();
  Eff_GEMRecHit_GE11_layer1_Pt->Write();
  Eff_GEMRecHit_GE11_layer1_Phi->Write();
  Eff_GEMRecHit_GE11_layer2_Eta->Write();
  Eff_GEMRecHit_GE11_layer2_Pt->Write();
  Eff_GEMRecHit_GE11_layer2_Phi->Write();
  Eff_GEMRecHit_GE21_layer1_Eta->Write();
  Eff_GEMRecHit_GE21_layer1_Pt->Write();
  Eff_GEMRecHit_GE21_layer1_Phi->Write();
  Eff_GEMRecHit_GE21_layer2_Eta->Write();
  Eff_GEMRecHit_GE21_layer2_Pt->Write();
  Eff_GEMRecHit_GE21_layer2_Phi->Write();
  Eff_GEMRecHit_GE11_two_Eta->Write();
  Eff_GEMRecHit_GE11_two_Pt->Write();
  Eff_GEMRecHit_GE11_two_Phi->Write();
  Eff_GEMRecHit_GE21_two_Eta->Write();
  Eff_GEMRecHit_GE21_two_Pt->Write();
  Eff_GEMRecHit_GE21_two_Phi->Write();
  // gemsegment
  TEfficiency* Eff_GEMSegment_GE11_Eta = new TEfficiency(*MatchedGEMSegment_GE11_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMSegment_GE11_Pt = new TEfficiency(*MatchedGEMSegment_GE11_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMSegment_GE11_Phi = new TEfficiency(*MatchedGEMSegment_GE11_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMSegment_GE21_Eta = new TEfficiency(*MatchedGEMSegment_GE21_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMSegment_GE21_Pt = new TEfficiency(*MatchedGEMSegment_GE21_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMSegment_GE21_Phi = new TEfficiency(*MatchedGEMSegment_GE21_Phi, *GenMuon_Phi);
  TEfficiency* Eff_GEMSegment_Eta = new TEfficiency(*MatchedGEMSegment_Eta, *GenMuon_Eta);
  TEfficiency* Eff_GEMSegment_Pt = new TEfficiency(*MatchedGEMSegment_Pt, *GenMuon_Pt);
  TEfficiency* Eff_GEMSegment_Phi = new TEfficiency(*MatchedGEMSegment_Phi, *GenMuon_Phi);
  Eff_GEMSegment_GE11_Eta->SetName("Eff_GEMSegment_GE11_Eta");
  Eff_GEMSegment_GE11_Pt->SetName("Eff_GEMSegment_GE11_Pt");
  Eff_GEMSegment_GE11_Phi->SetName("Eff_GEMSegment_GE11_Phi");
  Eff_GEMSegment_GE21_Eta->SetName("Eff_GEMSegment_GE21_Eta");
  Eff_GEMSegment_GE21_Pt->SetName("Eff_GEMSegment_GE21_Pt");
  Eff_GEMSegment_GE21_Phi->SetName("Eff_GEMSegment_GE21_Phi");
  Eff_GEMSegment_Eta->SetName("Eff_GEMSegment_Eta");
  Eff_GEMSegment_Pt->SetName("Eff_GEMSegment_Pt");
  Eff_GEMSegment_Phi->SetName("Eff_GEMSegment_Phi");
  Eff_GEMSegment_GE11_Eta->Write();
  Eff_GEMSegment_GE11_Pt->Write();
  Eff_GEMSegment_GE11_Phi->Write();
  Eff_GEMSegment_GE21_Eta->Write();
  Eff_GEMSegment_GE21_Pt->Write();
  Eff_GEMSegment_GE21_Phi->Write();
  Eff_GEMSegment_Eta->Write();
  Eff_GEMSegment_Pt->Write();
  Eff_GEMSegment_Phi->Write();
 
  /* Association by hits */

  TPMuon_Eta->Write();
  TPMuon_Pt->Write();
  TPMuon_Phi->Write();

  //==== GEMMuon
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

  //==== RecoMuon
  HitsMatchedRecoMuon_Eta->Write();
  HitsMatchedRecoMuon_Pt->Write();
  HitsMatchedRecoMuon_Phi->Write();
  // Efficeicny
  TEfficiency* HitsEff_RecoMuon_Eta = new TEfficiency(*HitsMatchedRecoMuon_Eta, *TPMuon_Eta);
  TEfficiency* HitsEff_RecoMuon_Pt = new TEfficiency(*HitsMatchedRecoMuon_Pt, *TPMuon_Pt);
  TEfficiency* HitsEff_RecoMuon_Phi = new TEfficiency(*HitsMatchedRecoMuon_Phi, *TPMuon_Phi);
  HitsEff_RecoMuon_Eta->SetName("HitsEff_RecoMuon_Eta");
  HitsEff_RecoMuon_Pt->SetName("HitsEff_RecoMuon_Pt");
  HitsEff_RecoMuon_Phi->SetName("HitsEff_RecoMuon_Phi");
  HitsEff_RecoMuon_Eta->Write();
  HitsEff_RecoMuon_Pt->Write();
  HitsEff_RecoMuon_Phi->Write();
  // Fake
  HitsUnmatchedRecoMuon_Eta->Write();
  HitsUnmatchedRecoMuon_Pt->Write();
  HitsUnmatchedRecoMuon_Phi->Write();

  //==== LooseMuon
  HitsMatchedLooseMuon_Eta->Write();
  HitsMatchedLooseMuon_Pt->Write();
  HitsMatchedLooseMuon_Phi->Write();
  // Efficeicny
  TEfficiency* HitsEff_LooseMuon_Eta = new TEfficiency(*HitsMatchedLooseMuon_Eta, *TPMuon_Eta);
  TEfficiency* HitsEff_LooseMuon_Pt = new TEfficiency(*HitsMatchedLooseMuon_Pt, *TPMuon_Pt);
  TEfficiency* HitsEff_LooseMuon_Phi = new TEfficiency(*HitsMatchedLooseMuon_Phi, *TPMuon_Phi);
  HitsEff_LooseMuon_Eta->SetName("HitsEff_LooseMuon_Eta");
  HitsEff_LooseMuon_Pt->SetName("HitsEff_LooseMuon_Pt");
  HitsEff_LooseMuon_Phi->SetName("HitsEff_LooseMuon_Phi");
  HitsEff_LooseMuon_Eta->Write();
  HitsEff_LooseMuon_Pt->Write();
  HitsEff_LooseMuon_Phi->Write();
  // Fake
  HitsUnmatchedLooseMuon_Eta->Write();
  HitsUnmatchedLooseMuon_Pt->Write();
  HitsUnmatchedLooseMuon_Phi->Write();

  //==== MediumMuon
  HitsMatchedMediumMuon_Eta->Write();
  HitsMatchedMediumMuon_Pt->Write();
  HitsMatchedMediumMuon_Phi->Write();
  // Efficeicny
  TEfficiency* HitsEff_MediumMuon_Eta = new TEfficiency(*HitsMatchedMediumMuon_Eta, *TPMuon_Eta);
  TEfficiency* HitsEff_MediumMuon_Pt = new TEfficiency(*HitsMatchedMediumMuon_Pt, *TPMuon_Pt);
  TEfficiency* HitsEff_MediumMuon_Phi = new TEfficiency(*HitsMatchedMediumMuon_Phi, *TPMuon_Phi);
  HitsEff_MediumMuon_Eta->SetName("HitsEff_MediumMuon_Eta");
  HitsEff_MediumMuon_Pt->SetName("HitsEff_MediumMuon_Pt");
  HitsEff_MediumMuon_Phi->SetName("HitsEff_MediumMuon_Phi");
  HitsEff_MediumMuon_Eta->Write();
  HitsEff_MediumMuon_Pt->Write();
  HitsEff_MediumMuon_Phi->Write();
  // Fake
  HitsUnmatchedMediumMuon_Eta->Write();
  HitsUnmatchedMediumMuon_Pt->Write();
  HitsUnmatchedMediumMuon_Phi->Write();

  //==== TightMuon
  HitsMatchedTightMuon_Eta->Write();
  HitsMatchedTightMuon_Pt->Write();
  HitsMatchedTightMuon_Phi->Write();
  // Efficeicny
  TEfficiency* HitsEff_TightMuon_Eta = new TEfficiency(*HitsMatchedTightMuon_Eta, *TPMuon_Eta);
  TEfficiency* HitsEff_TightMuon_Pt = new TEfficiency(*HitsMatchedTightMuon_Pt, *TPMuon_Pt);
  TEfficiency* HitsEff_TightMuon_Phi = new TEfficiency(*HitsMatchedTightMuon_Phi, *TPMuon_Phi);
  HitsEff_TightMuon_Eta->SetName("HitsEff_TightMuon_Eta");
  HitsEff_TightMuon_Pt->SetName("HitsEff_TightMuon_Pt");
  HitsEff_TightMuon_Phi->SetName("HitsEff_TightMuon_Phi");
  HitsEff_TightMuon_Eta->Write();
  HitsEff_TightMuon_Pt->Write();
  HitsEff_TightMuon_Phi->Write();
  // Fake
  HitsUnmatchedTightMuon_Eta->Write();
  HitsUnmatchedTightMuon_Pt->Write();
  HitsUnmatchedTightMuon_Phi->Write();

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
