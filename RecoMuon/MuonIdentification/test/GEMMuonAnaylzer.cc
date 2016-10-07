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
  std::string DoubleToString(double dd);

  //protected:
  
  private:

  //==== quick efficiency check
  Long64_t n_genmuon;
  Long64_t n_dR_matched_GEMmuon;
  Long64_t n_AssoByHits_matched_GEMmuon;

  edm::ESHandle<GEMGeometry> gemGeom;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<GEMRecHitCollection> GEMRecHit_Token;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > track_Collection_Token;

  bool UseAssociators;
  bool UseDeltaR;
  bool doGeometryStudy;
  const TrackAssociatorByChi2Impl* associatorByChi2;

  std::vector<std::string> associators;
  //std::vector<edm::InputTag> label;
  std::vector<std::string> label;
  bool doMatchingStudy;
  std::vector<double> PullXValues, DXValues, PullYValues, DYValues, DotDirValues;

  //Histos for plotting
  TFile* histoFile; 

  double  FakeRatePtCut, MatchingWindowDelR;

  TH1F *Nevents_h, *N_GEMMuon_h, *N_RecoMuon_h, *N_LooseMuon_h, *N_MediumMuon_h, *N_TightMuon_h,
                   *N_GEMMuon_ptcut_h, *N_RecoMuon_ptcut_h, *N_LooseMuon_ptcut_h, *N_MediumMuon_ptcut_h, *N_TightMuon_ptcut_h;
  int n_GEMMuon, n_RecoMuon, n_LooseMuon, n_MediumMuon, n_TightMuon;
  int n_GEMMuon_ptcut, n_RecoMuon_ptcut, n_LooseMuon_ptcut, n_MediumMuon_ptcut, n_TightMuon_ptcut;

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

  std::vector<int> n_GEMMuon_PullX, n_GEMMuon_ptcut_PullX;
  std::map< double, TH1F* > N_GEMMuon_PullX_h, N_GEMMuon_ptcut_PullX_h;
  std::map< double, TH1F* > HitsMatchedPullX_Eta, HitsMatchedPullX_Pt, HitsMatchedPullX_Phi,
                            HitsUnmatchedPullX_Eta, HitsUnmatchedPullX_Pt, HitsUnmatchedPullX_Phi;
  std::vector<int> n_GEMMuon_DX, n_GEMMuon_ptcut_DX;
  std::map< double, TH1F* > N_GEMMuon_DX_h, N_GEMMuon_ptcut_DX_h;
  std::map< double, TH1F* > HitsMatchedDX_Eta, HitsMatchedDX_Pt, HitsMatchedDX_Phi,
                            HitsUnmatchedDX_Eta, HitsUnmatchedDX_Pt, HitsUnmatchedDX_Phi;
  std::vector<int> n_GEMMuon_PullY, n_GEMMuon_ptcut_PullY;
  std::map< double, TH1F* > N_GEMMuon_PullY_h, N_GEMMuon_ptcut_PullY_h;
  std::map< double, TH1F* > HitsMatchedPullY_Eta, HitsMatchedPullY_Pt, HitsMatchedPullY_Phi,
                            HitsUnmatchedPullY_Eta, HitsUnmatchedPullY_Pt, HitsUnmatchedPullY_Phi;
  std::vector<int> n_GEMMuon_DY, n_GEMMuon_ptcut_DY;
  std::map< double, TH1F* > N_GEMMuon_DY_h, N_GEMMuon_ptcut_DY_h;
  std::map< double, TH1F* > HitsMatchedDY_Eta, HitsMatchedDY_Pt, HitsMatchedDY_Phi,
                            HitsUnmatchedDY_Eta, HitsUnmatchedDY_Pt, HitsUnmatchedDY_Phi;
  std::vector<int> n_GEMMuon_DotDir, n_GEMMuon_ptcut_DotDir;
  std::map< double, TH1F* > N_GEMMuon_DotDir_h, N_GEMMuon_ptcut_DotDir_h;
  std::map< double, TH1F* > HitsMatchedDotDir_Eta, HitsMatchedDotDir_Pt, HitsMatchedDotDir_Phi,
                            HitsUnmatchedDotDir_Eta, HitsUnmatchedDotDir_Pt, HitsUnmatchedDotDir_Phi;

};

GEMMuonAnalyzer::GEMMuonAnalyzer(const edm::ParameterSet& iConfig) 
{
  histoFile = new TFile(iConfig.getParameter<std::string>("HistoFile").c_str(), "recreate");
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  UseDeltaR = iConfig.getParameter< bool >("UseDeltaR");
  doGeometryStudy = iConfig.getParameter< bool >("doGeometryStudy");

  FakeRatePtCut   = iConfig.getParameter<double>("FakeRatePtCut");
  MatchingWindowDelR   = iConfig.getParameter<double>("MatchingWindowDelR");

  //Associator for chi2: getting parameters
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  associators = iConfig.getParameter< std::vector<std::string> >("associators");

  //label = iConfig.getParameter< std::vector<edm::InputTag> >("label");
  label = iConfig.getParameter< std::vector<std::string> >("label");
  doMatchingStudy = iConfig.getParameter< bool >("doMatchingStudy");
  if(doMatchingStudy){
    PullXValues = iConfig.getParameter< std::vector<double> >("PullXValues");
    DXValues = iConfig.getParameter< std::vector<double> >("DXValues");
    PullYValues = iConfig.getParameter< std::vector<double> >("PullYValues");
    DYValues = iConfig.getParameter< std::vector<double> >("DYValues");
    DotDirValues = iConfig.getParameter< std::vector<double> >("DotDirValues");
  }
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

  n_genmuon = 0;
  n_dR_matched_GEMmuon = 0;
  n_AssoByHits_matched_GEMmuon = 0;

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
  n_GEMMuon_PullX.clear();
  n_GEMMuon_ptcut_PullX.clear();
  n_GEMMuon_PullY.clear();
  n_GEMMuon_ptcut_PullY.clear();
  n_GEMMuon_PullY.clear();
  n_GEMMuon_ptcut_PullY.clear();
  n_GEMMuon_DY.clear();
  n_GEMMuon_ptcut_DY.clear();
  n_GEMMuon_DotDir.clear();
  n_GEMMuon_ptcut_DotDir.clear();

  for(unsigned int i=0; i<PullXValues.size(); i++){
    n_GEMMuon_PullX.push_back(0);
    n_GEMMuon_ptcut_PullX.push_back(0);
  }
  for(unsigned int i=0; i<DXValues.size(); i++){
    n_GEMMuon_DX.push_back(0);
    n_GEMMuon_ptcut_DX.push_back(0);
  }
  for(unsigned int i=0; i<PullYValues.size(); i++){
    n_GEMMuon_PullY.push_back(0);
    n_GEMMuon_ptcut_PullY.push_back(0);
  }
  for(unsigned int i=0; i<DYValues.size(); i++){
    n_GEMMuon_DY.push_back(0);
    n_GEMMuon_ptcut_DY.push_back(0);
  }
  for(unsigned int i=0; i<DotDirValues.size(); i++){
    n_GEMMuon_DotDir.push_back(0);
    n_GEMMuon_ptcut_DotDir.push_back(0);
  }
  

  std::cout<<"Contructor end"<<std::endl;
}



void GEMMuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

  const int n_pt_bin = 19, n_eta_bin = 9;
  double pt_bin[n_pt_bin+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  double eta_bin[n_eta_bin+1] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

  //Histos for plotting
  Nevents_h = new TH1F("Nevents_h", "Nevents", 2, 0, 2 );
  N_GEMMuon_h = new TH1F("N_GEMMuon_h", "Nevents", 1, 0, 1 );
  N_RecoMuon_h = new TH1F("N_RecoMuon_h", "Nevents", 1, 0, 1 );
  N_LooseMuon_h = new TH1F("N_LooseMuon_h", "Nevents", 1, 0, 1 );
  N_MediumMuon_h = new TH1F("N_MediumMuon_h", "Nevents", 1, 0, 1 );
  N_TightMuon_h = new TH1F("N_TightMuon_h", "Nevents", 1, 0, 1 );
  N_GEMMuon_ptcut_h = new TH1F("N_GEMMuon_ptcut_h", "Nevents", 1, 0, 1 );
  N_RecoMuon_ptcut_h = new TH1F("N_RecoMuon_ptcut_h", "Nevents", 1, 0, 1 );
  N_LooseMuon_ptcut_h = new TH1F("N_LooseMuon_ptcut_h", "Nevents", 1, 0, 1 );
  N_MediumMuon_ptcut_h = new TH1F("N_MediumMuon_ptcut_h", "Nevents", 1, 0, 1 );
  N_TightMuon_ptcut_h = new TH1F("N_TightMuon_ptcut_h", "Nevents", 1, 0, 1 );
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

  for(unsigned int i=0; i<PullXValues.size(); i++){
    double aaa = PullXValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsMatchedPullX_Eta[aaa] = new TH1F("HitsMatchedPullX_Eta"+saaa, "PullX #eta", n_eta_bin, eta_bin );
    HitsMatchedPullX_Pt[aaa]  = new TH1F("HitsMatchedPullX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedPullX_Phi[aaa] = new TH1F("HitsMatchedPullX_Phi"+saaa, "PullX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedPullX_Eta[aaa] = new TH1F("HitsUnmatchedPullX_Eta"+saaa, "PullX #eta", n_eta_bin, eta_bin );
    HitsUnmatchedPullX_Pt[aaa]  = new TH1F("HitsUnmatchedPullX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedPullX_Phi[aaa] = new TH1F("HitsUnmatchedPullX_Phi"+saaa, "PullX #phi", 36, -TMath::Pi(), TMath::Pi() );
    N_GEMMuon_PullX_h[aaa] = new TH1F("N_GEMMuon_PullX_h"+saaa, "Nevents", 1, 0, 1 );
    N_GEMMuon_ptcut_PullX_h[aaa] = new TH1F("N_GEMMuon_ptcut_PullX_h"+saaa, "Nevents", 1, 0, 1 );
  }
  for(unsigned int i=0; i<DXValues.size(); i++){
    double aaa = DXValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsMatchedDX_Eta[aaa] = new TH1F("HitsMatchedDX_Eta"+saaa, "DX #eta", n_eta_bin, eta_bin );
    HitsMatchedDX_Pt[aaa]  = new TH1F("HitsMatchedDX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedDX_Phi[aaa] = new TH1F("HitsMatchedDX_Phi"+saaa, "DX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedDX_Eta[aaa] = new TH1F("HitsUnmatchedDX_Eta"+saaa, "DX #eta", n_eta_bin, eta_bin );
    HitsUnmatchedDX_Pt[aaa]  = new TH1F("HitsUnmatchedDX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedDX_Phi[aaa] = new TH1F("HitsUnmatchedDX_Phi"+saaa, "DX #phi", 36, -TMath::Pi(), TMath::Pi() );
    N_GEMMuon_DX_h[aaa] = new TH1F("N_GEMMuon_DX_h"+saaa, "Nevents", 1, 0, 1 );
    N_GEMMuon_ptcut_DX_h[aaa] = new TH1F("N_GEMMuon_ptcut_DX_h"+saaa, "Nevents", 1, 0, 1 );
  }
  for(unsigned int i=0; i<PullYValues.size(); i++){
    double aaa = PullYValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsMatchedPullY_Eta[aaa] = new TH1F("HitsMatchedPullY_Eta"+saaa, "PullY #eta", n_eta_bin, eta_bin );
    HitsMatchedPullY_Pt[aaa]  = new TH1F("HitsMatchedPullY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedPullY_Phi[aaa] = new TH1F("HitsMatchedPullY_Phi"+saaa, "PullY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedPullY_Eta[aaa] = new TH1F("HitsUnmatchedPullY_Eta"+saaa, "PullY #eta", n_eta_bin, eta_bin );
    HitsUnmatchedPullY_Pt[aaa]  = new TH1F("HitsUnmatchedPullY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedPullY_Phi[aaa] = new TH1F("HitsUnmatchedPullY_Phi"+saaa, "PullY #phi", 36, -TMath::Pi(), TMath::Pi() );
    N_GEMMuon_PullY_h[aaa] = new TH1F("N_GEMMuon_PullY_h"+saaa, "Nevents", 1, 0, 1 );
    N_GEMMuon_ptcut_PullY_h[aaa] = new TH1F("N_GEMMuon_ptcut_PullY_h"+saaa, "Nevents", 1, 0, 1 );
  }
  for(unsigned int i=0; i<DYValues.size(); i++){
    double aaa = DYValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsMatchedDY_Eta[aaa] = new TH1F("HitsMatchedDY_Eta"+saaa, "DY #eta", n_eta_bin, eta_bin );
    HitsMatchedDY_Pt[aaa]  = new TH1F("HitsMatchedDY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedDY_Phi[aaa] = new TH1F("HitsMatchedDY_Phi"+saaa, "DY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedDY_Eta[aaa] = new TH1F("HitsUnmatchedDY_Eta"+saaa, "DY #eta", n_eta_bin, eta_bin );
    HitsUnmatchedDY_Pt[aaa]  = new TH1F("HitsUnmatchedDY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedDY_Phi[aaa] = new TH1F("HitsUnmatchedDY_Phi"+saaa, "DY #phi", 36, -TMath::Pi(), TMath::Pi() );
    N_GEMMuon_DY_h[aaa] = new TH1F("N_GEMMuon_DY_h"+saaa, "Nevents", 1, 0, 1 );
    N_GEMMuon_ptcut_DY_h[aaa] = new TH1F("N_GEMMuon_ptcut_DY_h"+saaa, "Nevents", 1, 0, 1 );
  }
  for(unsigned int i=0; i<DotDirValues.size(); i++){
    double aaa = DotDirValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsMatchedDotDir_Eta[aaa] = new TH1F("HitsMatchedDotDir_Eta"+saaa, "DotDir #eta", n_eta_bin, eta_bin );
    HitsMatchedDotDir_Pt[aaa]  = new TH1F("HitsMatchedDotDir_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedDotDir_Phi[aaa] = new TH1F("HitsMatchedDotDir_Phi"+saaa, "DotDir #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedDotDir_Eta[aaa] = new TH1F("HitsUnmatchedDotDir_Eta"+saaa, "DotDir #eta", n_eta_bin, eta_bin );
    HitsUnmatchedDotDir_Pt[aaa]  = new TH1F("HitsUnmatchedDotDir_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedDotDir_Phi[aaa] = new TH1F("HitsUnmatchedDotDir_Phi"+saaa, "DotDir #phi", 36, -TMath::Pi(), TMath::Pi() );
    N_GEMMuon_DotDir_h[aaa] = new TH1F("N_GEMMuon_DotDir_h"+saaa, "Nevents", 1, 0, 1 );
    N_GEMMuon_ptcut_DotDir_h[aaa] = new TH1F("N_GEMMuon_ptcut_DotDir_h"+saaa, "Nevents", 1, 0, 1 );
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

  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);

  iSetup.get<MuonGeometryRecord>().get(gemGeom);
  edm::Handle<GEMRecHitCollection> gemRecHitCollection;
  iEvent.getByToken(GEMRecHit_Token, gemRecHitCollection);
  edm::Handle<GEMSegmentCollection> gemSegmentCollection;
  iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);

  if(doGeometryStudy){

    //==== GEMRecHit study

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

  }

  //==== DeltaR matching 

  if(UseDeltaR){

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

        double RecoMuon_LowestDelR = 9999, RecoMuon_NotGEMMuon_LowestDelR = 9999, RecoMuon_GEMMuon_h_LowestDelR = 9999;
        double RecoMuon_thisDelR = 9999;
        bool RecoMuon_isMatched = false, RecoMuon_NotGEMMuon_isMatched = false, RecoMuoN_GEMMuon_h_isMatched = false;
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
                RecoMuoN_GEMMuon_h_isMatched = true;
                if( RecoMuon_thisDelR < RecoMuon_GEMMuon_h_LowestDelR ){
                  RecoMuon_GEMMuon_h_LowestDelR = RecoMuon_thisDelR;
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
        if( RecoMuoN_GEMMuon_h_isMatched ){
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

  } 
  
  if(UseAssociators) {

    unsigned int www_PullX = 0, www_DX = 0, www_PullY = 0, www_DY = 0, www_DotDir = 0;
    for (unsigned int www=0;www<label.size();www++){

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

      //=======================
      //==== Efficiency study
      //=======================

      //std::cout << "TrackingParticle size = " << trackingParticles->size() << std::endl;
      for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++){
        TrackingParticleRef tpr(trackingParticles, i);
        TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
        TrackingParticle::Vector momentumTP;
        TrackingParticle::Point vertexTP;

        if (abs(tp->pdgId()) != 13) continue;

        bool Eta_1p6_2p4 = fabs(tp->eta()) > 1.6 && fabs(tp->eta()) < 2.4;
        bool Pt_5 = tp->pt() > 5;

        if( Eta_1p6_2p4 ){
          bool SignalMuon = false;

          int motherId(0), grandmaId(0), greatgrandmaId(0);
          if( tp->genParticles().size()>0 && (*tp->genParticle_begin())->numberOfMothers()>0 ) {
            motherId = abs( (*tp->genParticle_begin())->mother()->pdgId() );
            if( (*tp->genParticle_begin())->mother()->numberOfMothers()>0 ) {
              grandmaId = abs( (*tp->genParticle_begin())->mother()->mother()->pdgId() );
              if( (*tp->genParticle_begin())->mother()->mother()->numberOfMothers()>0 ) {
                greatgrandmaId = abs( (*tp->genParticle_begin())->mother()->mother()->mother()->pdgId() );
              }
            }
          }
          // Is it a signal muon? 
          SignalMuon = ( grandmaId==23 || greatgrandmaId==23 ||
                       (motherId==13 && grandmaId==13 && greatgrandmaId==13)   );

          if(SignalMuon){
            //std::cout
            //<< i << '\t'
            //<< tp->status() << '\t'
            //<< motherId << '\t'
            //<< grandmaId << '\t'
            //<< greatgrandmaId << '\t'
            //<< std::endl;

            //==== Fill the Denominator
            //==== Should be filled only once (www=0)
            if( www == 0 ){
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
                    n_AssoByHits_matched_GEMmuon++;
                    if(Pt_5) HitsMatchedGEMMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedGEMMuon_Pt->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedGEMMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="recoMuonSel"){
                    if(Pt_5) HitsMatchedRecoMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedRecoMuon_Pt->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedRecoMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="looseMuonSel"){
                    if(Pt_5) HitsMatchedLooseMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedLooseMuon_Pt->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedLooseMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="mediumMuonSel"){
                    if(Pt_5) HitsMatchedMediumMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedMediumMuon_Pt->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedMediumMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www]=="tightMuonSel"){
                    if(Pt_5) HitsMatchedTightMuon_Eta->Fill(fabs(tpr->eta()));
                    HitsMatchedTightMuon_Pt->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedTightMuon_Phi->Fill(tpr->phi());
                  }
                  if(label[www].find("PullXScan") != std::string::npos){
                    double this_cut = PullXValues.at(www_PullX);
                    if(Pt_5) HitsMatchedPullX_Eta[this_cut]->Fill(fabs(tpr->eta()));
                    HitsMatchedPullX_Pt[this_cut]->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedPullX_Phi[this_cut]->Fill(tpr->phi());
                  }
                  if(label[www].find("DXScan") != std::string::npos){
                    double this_cut = DXValues.at(www_DX);
                    if(Pt_5) HitsMatchedDX_Eta[this_cut]->Fill(fabs(tpr->eta()));
                    HitsMatchedDX_Pt[this_cut]->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedDX_Phi[this_cut]->Fill(tpr->phi());
                  }
                  if(label[www].find("PullYScan") != std::string::npos){
                    double this_cut = PullYValues.at(www_PullY);
                    if(Pt_5) HitsMatchedPullY_Eta[this_cut]->Fill(fabs(tpr->eta()));
                    HitsMatchedPullY_Pt[this_cut]->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedPullY_Phi[this_cut]->Fill(tpr->phi());
                  }
                  if(label[www].find("DYScan") != std::string::npos){
                    double this_cut = DYValues.at(www_DY);
                    if(Pt_5) HitsMatchedDY_Eta[this_cut]->Fill(fabs(tpr->eta()));
                    HitsMatchedDY_Pt[this_cut]->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedDY_Phi[this_cut]->Fill(tpr->phi());
                  }
                  if(label[www].find("DotDirScan") != std::string::npos){
                    double this_cut = DotDirValues.at(www_DotDir);
                    if(Pt_5) HitsMatchedDotDir_Eta[this_cut]->Fill(fabs(tpr->eta()));
                    HitsMatchedDotDir_Pt[this_cut]->Fill(tpr->pt());
                    if(Pt_5) HitsMatchedDotDir_Phi[this_cut]->Fill(tpr->phi());
                  }

                }
              }
            }

          } //==== END if(SignalMuon)
        } // END if( Eta_1p6_2p4 )


      } // END TrackingParticle Loop

      //=================
      //==== Fake study
      //=================

      //==== loop over (GEMMuon/RecoMuon/...) tracks
      for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
        //std::cout << i << "th trackCollection iterator" << std::endl;
        RefToBase<Track> track(trackCollection, i);

        bool Eta_1p6_2p4 = fabs(track->eta()) > 1.6 && fabs(track->eta()) < 2.4;
        bool Pt_5 = track->pt() > 5;

        if( Eta_1p6_2p4 ){

          //==== count reco tracks
          if(label[www]=="gemMuonSel"){
            n_GEMMuon++;
            if( Pt_5 ) n_GEMMuon_ptcut++;
          }
          if(label[www]=="recoMuonSel"){
            n_RecoMuon++;
            if( Pt_5 ) n_RecoMuon_ptcut++;
          }
          if(label[www]=="looseMuonSel"){
            n_LooseMuon++;
            if( Pt_5 ) n_LooseMuon_ptcut++;
          }
          if(label[www]=="mediumMuonSel"){
            n_MediumMuon++;
            if( Pt_5 ) n_MediumMuon_ptcut++;
          }
          if(label[www]=="tightMuonSel"){
            n_TightMuon++;
            if( Pt_5 ) n_TightMuon_ptcut++;
          }
          if(label[www].find("PullXScan") != std::string::npos){
            n_GEMMuon_PullX[www_PullX]++;
            if( Pt_5 ) n_GEMMuon_ptcut_PullX[www_PullX]++;
          }
          if(label[www].find("DXScan") != std::string::npos){
            n_GEMMuon_DX[www_DX]++;
            if( Pt_5 ) n_GEMMuon_ptcut_DX[www_DX]++;
          }
          if(label[www].find("PullYScan") != std::string::npos){
            n_GEMMuon_PullY[www_PullY]++;
            if( Pt_5 ) n_GEMMuon_ptcut_PullY[www_PullY]++;
          }
          if(label[www].find("DYScan") != std::string::npos){
            n_GEMMuon_DY[www_DY]++;
            if( Pt_5 ) n_GEMMuon_ptcut_DY[www_DY]++;
          }
          if(label[www].find("DotDirScan") != std::string::npos){
            n_GEMMuon_DotDir[www_DotDir]++;
            if( Pt_5 ) n_GEMMuon_ptcut_DotDir[www_DotDir]++;
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

              int motherId(0), grandmaId(0), greatgrandmaId(0);
              if( tpr->genParticles().size()>0 && (*tpr->genParticle_begin())->numberOfMothers()>0 ) {
                motherId = abs( (*tpr->genParticle_begin())->mother()->pdgId() );
                if( (*tpr->genParticle_begin())->mother()->numberOfMothers()>0 ) {
                  grandmaId = abs( (*tpr->genParticle_begin())->mother()->mother()->pdgId() );
                  if( (*tpr->genParticle_begin())->mother()->mother()->numberOfMothers()>0 ) {
                    greatgrandmaId = abs( (*tpr->genParticle_begin())->mother()->mother()->mother()->pdgId() );
                  }
                }
              }
              // Is it a signal muon? 
              SignalMuon = ( grandmaId==23 || greatgrandmaId==23 ||
                           (motherId==13 && grandmaId==13 && greatgrandmaId==13)   );

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
            if(label[www]=="recoMuonSel"){
              if(Pt_5) HitsUnmatchedRecoMuon_Eta->Fill(fabs(track->eta()));
              HitsUnmatchedRecoMuon_Pt->Fill(track->pt());
              if(Pt_5) HitsUnmatchedRecoMuon_Phi->Fill(track->phi());
            }
            if(label[www]=="looseMuonSel"){
              if(Pt_5) HitsUnmatchedLooseMuon_Eta->Fill(fabs(track->eta()));
              HitsUnmatchedLooseMuon_Pt->Fill(track->pt());
              if(Pt_5) HitsUnmatchedLooseMuon_Phi->Fill(track->phi());
            }
            if(label[www]=="mediumMuonSel"){
              if(Pt_5) HitsUnmatchedMediumMuon_Eta->Fill(fabs(track->eta()));
              HitsUnmatchedMediumMuon_Pt->Fill(track->pt());
              if(Pt_5) HitsUnmatchedMediumMuon_Phi->Fill(track->phi());
            }
            if(label[www]=="tightMuonSel"){
              if(Pt_5) HitsUnmatchedTightMuon_Eta->Fill(fabs(track->eta()));
              HitsUnmatchedTightMuon_Pt->Fill(track->pt());
              if(Pt_5) HitsUnmatchedTightMuon_Phi->Fill(track->phi());
            }
            if(label[www].find("PullXScan") != std::string::npos){
              double this_cut = PullXValues.at(www_PullX);
              if(Pt_5) HitsUnmatchedPullX_Eta[this_cut]->Fill(fabs(track->eta()));
              HitsUnmatchedPullX_Pt[this_cut]->Fill(track->pt());
              if(Pt_5) HitsUnmatchedPullX_Phi[this_cut]->Fill(track->phi());
            }
            if(label[www].find("DXScan") != std::string::npos){
              double this_cut = DXValues.at(www_DX);
              if(Pt_5) HitsUnmatchedDX_Eta[this_cut]->Fill(fabs(track->eta()));
              HitsUnmatchedDX_Pt[this_cut]->Fill(track->pt());
              if(Pt_5) HitsUnmatchedDX_Phi[this_cut]->Fill(track->phi());
            }
            if(label[www].find("PullYScan") != std::string::npos){
              double this_cut = PullYValues.at(www_PullY);
              if(Pt_5) HitsUnmatchedPullY_Eta[this_cut]->Fill(fabs(track->eta()));
              HitsUnmatchedPullY_Pt[this_cut]->Fill(track->pt());
              if(Pt_5) HitsUnmatchedPullY_Phi[this_cut]->Fill(track->phi());
            }
            if(label[www].find("DYScan") != std::string::npos){
              double this_cut = DYValues.at(www_DY);
              if(Pt_5) HitsUnmatchedDY_Eta[this_cut]->Fill(fabs(track->eta()));
              HitsUnmatchedDY_Pt[this_cut]->Fill(track->pt());
              if(Pt_5) HitsUnmatchedDY_Phi[this_cut]->Fill(track->phi());
            }
            if(label[www].find("DotDirScan") != std::string::npos){
              double this_cut = DotDirValues.at(www_DotDir);
              if(Pt_5) HitsUnmatchedDotDir_Eta[this_cut]->Fill(fabs(track->eta()));
              HitsUnmatchedDotDir_Pt[this_cut]->Fill(track->pt());
              if(Pt_5) HitsUnmatchedDotDir_Phi[this_cut]->Fill(track->phi());
            }

          } //==== END if(isFake)


        } // ==== END if( Eta_1p6_2p4 )


      } // END track loop

      if(label[www].find("PullXScan") != std::string::npos) www_PullX++;
      if(label[www].find("DXScan") != std::string::npos) www_DX++;
      if(label[www].find("PullYScan") != std::string::npos) www_PullY++;
      if(label[www].find("DYScan") != std::string::npos) www_DY++;
      if(label[www].find("DotDirScan") != std::string::npos) www_DotDir++;

    } // END label Loop


  } // END Use Associator






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

  N_GEMMuon_h->SetBinContent(1, n_GEMMuon);
  N_RecoMuon_h->SetBinContent(1,n_RecoMuon);
  N_LooseMuon_h->SetBinContent(1, n_LooseMuon);
  N_MediumMuon_h->SetBinContent(1, n_MediumMuon);
  N_TightMuon_h->SetBinContent(1, n_TightMuon);
  N_GEMMuon_ptcut_h->SetBinContent(1, n_GEMMuon_ptcut);
  N_RecoMuon_ptcut_h->SetBinContent(1,n_RecoMuon_ptcut);
  N_LooseMuon_ptcut_h->SetBinContent(1, n_LooseMuon_ptcut);
  N_MediumMuon_ptcut_h->SetBinContent(1, n_MediumMuon_ptcut);
  N_TightMuon_ptcut_h->SetBinContent(1, n_TightMuon_ptcut);

  N_GEMMuon_h->Write();
  N_RecoMuon_h->Write();
  N_LooseMuon_h->Write();
  N_MediumMuon_h->Write();
  N_TightMuon_h->Write();
  N_GEMMuon_ptcut_h->Write();
  N_RecoMuon_ptcut_h->Write();
  N_LooseMuon_ptcut_h->Write();
  N_MediumMuon_ptcut_h->Write();
  N_TightMuon_ptcut_h->Write();

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

  TH1F *hist_PullXValues = new TH1F("PullXValues", "", PullXValues.size(), 0, PullXValues.size());
  for(unsigned int i=0; i<PullXValues.size(); i++){
    double aaa = PullXValues.at(i);
    hist_PullXValues->SetBinContent(i+1, aaa);   
    TString saaa = DoubleToString(aaa);
    HitsMatchedPullX_Eta[aaa]->Write();
    HitsMatchedPullX_Pt[aaa]->Write();
    HitsMatchedPullX_Phi[aaa]->Write();
    //==== Efficiency
    TEfficiency* HitsEff_Eta = new TEfficiency(*HitsMatchedPullX_Eta[aaa], *TPMuon_Eta);
    TEfficiency* HitsEff_Pt = new TEfficiency(*HitsMatchedPullX_Pt[aaa], *TPMuon_Pt);
    TEfficiency* HitsEff_Phi = new TEfficiency(*HitsMatchedPullX_Phi[aaa], *TPMuon_Phi);
    HitsEff_Eta->SetName("HitsEff_PullX_Eta_"+saaa);
    HitsEff_Pt->SetName("HitsEff_PullX_Pt_"+saaa);
    HitsEff_Phi->SetName("HitsEff_PullX_Phi_"+saaa);
    HitsEff_Eta->Write();
    HitsEff_Pt->Write();
    HitsEff_Phi->Write();

    HitsUnmatchedPullX_Eta[aaa]->Write();
    HitsUnmatchedPullX_Pt[aaa]->Write();
    HitsUnmatchedPullX_Phi[aaa]->Write();

    N_GEMMuon_PullX_h[aaa]->SetBinContent(1, n_GEMMuon_PullX.at(i));
    N_GEMMuon_PullX_h[aaa]->Write();
    N_GEMMuon_ptcut_PullX_h[aaa]->SetBinContent(1, n_GEMMuon_ptcut_PullX.at(i));
    N_GEMMuon_ptcut_PullX_h[aaa]->Write();
  }
  hist_PullXValues->Write();

  TH1F *hist_DXValues = new TH1F("DXValues", "", DXValues.size(), 0, DXValues.size());
  for(unsigned int i=0; i<DXValues.size(); i++){
    double aaa = DXValues.at(i);
    hist_DXValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsMatchedDX_Eta[aaa]->Write();
    HitsMatchedDX_Pt[aaa]->Write();
    HitsMatchedDX_Phi[aaa]->Write();
    //==== Efficiency
    TEfficiency* HitsEff_Eta = new TEfficiency(*HitsMatchedDX_Eta[aaa], *TPMuon_Eta);
    TEfficiency* HitsEff_Pt = new TEfficiency(*HitsMatchedDX_Pt[aaa], *TPMuon_Pt);
    TEfficiency* HitsEff_Phi = new TEfficiency(*HitsMatchedDX_Phi[aaa], *TPMuon_Phi);
    HitsEff_Eta->SetName("HitsEff_DX_Eta_"+saaa);
    HitsEff_Pt->SetName("HitsEff_DX_Pt_"+saaa);
    HitsEff_Phi->SetName("HitsEff_DX_Phi_"+saaa);
    HitsEff_Eta->Write();
    HitsEff_Pt->Write();
    HitsEff_Phi->Write();

    HitsUnmatchedDX_Eta[aaa]->Write();
    HitsUnmatchedDX_Pt[aaa]->Write();
    HitsUnmatchedDX_Phi[aaa]->Write();

    N_GEMMuon_DX_h[aaa]->SetBinContent(1, n_GEMMuon_DX.at(i));
    N_GEMMuon_DX_h[aaa]->Write();
    N_GEMMuon_ptcut_DX_h[aaa]->SetBinContent(1, n_GEMMuon_ptcut_DX.at(i));
    N_GEMMuon_ptcut_DX_h[aaa]->Write();
  }
  hist_DXValues->Write();

  TH1F *hist_PullYValues = new TH1F("PullYValues", "", PullYValues.size(), 0, PullYValues.size());
  for(unsigned int i=0; i<PullYValues.size(); i++){
    double aaa = PullYValues.at(i);
    hist_PullYValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsMatchedPullY_Eta[aaa]->Write();
    HitsMatchedPullY_Pt[aaa]->Write();
    HitsMatchedPullY_Phi[aaa]->Write();
    //==== Efficiency
    TEfficiency* HitsEff_Eta = new TEfficiency(*HitsMatchedPullY_Eta[aaa], *TPMuon_Eta);
    TEfficiency* HitsEff_Pt = new TEfficiency(*HitsMatchedPullY_Pt[aaa], *TPMuon_Pt);
    TEfficiency* HitsEff_Phi = new TEfficiency(*HitsMatchedPullY_Phi[aaa], *TPMuon_Phi);
    HitsEff_Eta->SetName("HitsEff_PullY_Eta_"+saaa);
    HitsEff_Pt->SetName("HitsEff_PullY_Pt_"+saaa);
    HitsEff_Phi->SetName("HitsEff_PullY_Phi_"+saaa);
    HitsEff_Eta->Write();
    HitsEff_Pt->Write();
    HitsEff_Phi->Write();

    HitsUnmatchedPullY_Eta[aaa]->Write();
    HitsUnmatchedPullY_Pt[aaa]->Write();
    HitsUnmatchedPullY_Phi[aaa]->Write();

    N_GEMMuon_PullY_h[aaa]->SetBinContent(1, n_GEMMuon_PullY.at(i));
    N_GEMMuon_PullY_h[aaa]->Write();
    N_GEMMuon_ptcut_PullY_h[aaa]->SetBinContent(1, n_GEMMuon_ptcut_PullY.at(i));
    N_GEMMuon_ptcut_PullY_h[aaa]->Write();
  }
  hist_PullYValues->Write();

  TH1F *hist_DYValues = new TH1F("DYValues", "", DYValues.size(), 0, DYValues.size());
  for(unsigned int i=0; i<DYValues.size(); i++){
    double aaa = DYValues.at(i);
    hist_DYValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsMatchedDY_Eta[aaa]->Write();
    HitsMatchedDY_Pt[aaa]->Write();
    HitsMatchedDY_Phi[aaa]->Write();
    //==== Efficiency
    TEfficiency* HitsEff_Eta = new TEfficiency(*HitsMatchedDY_Eta[aaa], *TPMuon_Eta);
    TEfficiency* HitsEff_Pt = new TEfficiency(*HitsMatchedDY_Pt[aaa], *TPMuon_Pt);
    TEfficiency* HitsEff_Phi = new TEfficiency(*HitsMatchedDY_Phi[aaa], *TPMuon_Phi);
    HitsEff_Eta->SetName("HitsEff_DY_Eta_"+saaa);
    HitsEff_Pt->SetName("HitsEff_DY_Pt_"+saaa);
    HitsEff_Phi->SetName("HitsEff_DY_Phi_"+saaa);
    HitsEff_Eta->Write();
    HitsEff_Pt->Write();
    HitsEff_Phi->Write();

    HitsUnmatchedDY_Eta[aaa]->Write();
    HitsUnmatchedDY_Pt[aaa]->Write();
    HitsUnmatchedDY_Phi[aaa]->Write();

    N_GEMMuon_DY_h[aaa]->SetBinContent(1, n_GEMMuon_DY.at(i));
    N_GEMMuon_DY_h[aaa]->Write();
    N_GEMMuon_ptcut_DY_h[aaa]->SetBinContent(1, n_GEMMuon_ptcut_DY.at(i));
    N_GEMMuon_ptcut_DY_h[aaa]->Write();
  }
  hist_DYValues->Write();

  TH1F *hist_DotDirValues = new TH1F("DotDirValues", "", DotDirValues.size(), 0, DotDirValues.size());
  for(unsigned int i=0; i<DotDirValues.size(); i++){
    double aaa = DotDirValues.at(i);
    hist_DotDirValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsMatchedDotDir_Eta[aaa]->Write();
    HitsMatchedDotDir_Pt[aaa]->Write();
    HitsMatchedDotDir_Phi[aaa]->Write();
    //==== Efficiency
    TEfficiency* HitsEff_Eta = new TEfficiency(*HitsMatchedDotDir_Eta[aaa], *TPMuon_Eta);
    TEfficiency* HitsEff_Pt = new TEfficiency(*HitsMatchedDotDir_Pt[aaa], *TPMuon_Pt);
    TEfficiency* HitsEff_Phi = new TEfficiency(*HitsMatchedDotDir_Phi[aaa], *TPMuon_Phi);
    HitsEff_Eta->SetName("HitsEff_DotDir_Eta_"+saaa);
    HitsEff_Pt->SetName("HitsEff_DotDir_Pt_"+saaa);
    HitsEff_Phi->SetName("HitsEff_DotDir_Phi_"+saaa);
    HitsEff_Eta->Write();
    HitsEff_Pt->Write();
    HitsEff_Phi->Write();

    HitsUnmatchedDotDir_Eta[aaa]->Write();
    HitsUnmatchedDotDir_Pt[aaa]->Write();
    HitsUnmatchedDotDir_Phi[aaa]->Write();

    N_GEMMuon_DotDir_h[aaa]->SetBinContent(1, n_GEMMuon_DotDir.at(i));
    N_GEMMuon_DotDir_h[aaa]->Write();
    N_GEMMuon_ptcut_DotDir_h[aaa]->SetBinContent(1, n_GEMMuon_ptcut_DotDir.at(i));
    N_GEMMuon_ptcut_DotDir_h[aaa]->Write();
  }
  hist_DotDirValues->Write();

  histoFile->cd();

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

std::string GEMMuonAnalyzer::DoubleToString(double dd){
  std::ostringstream os;
  os << dd;
  return os.str();
}

DEFINE_FWK_MODULE(GEMMuonAnalyzer);
