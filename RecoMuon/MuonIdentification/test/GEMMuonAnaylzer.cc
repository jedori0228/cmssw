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

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
  edm::EDGetTokenT<GEMRecHitCollection> GEMRecHit_Token;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > track_Collection_Token;

  bool UseAssociators;
  bool doGeometryStudy;
  std::string SampleProcess;

  std::vector<std::string> associators;
  std::vector<std::string> label;
  bool doMatchingStudy;
  std::vector<double> PullXValues, DXValues, PullYValues, DYValues, DotDirValues;

  //====Histos for plotting
  TFile* histoFile; 

  //===============
  //==== Geometry
  //===============

  TH2F *GEMRecHit_GE11_GlobalPosition_scattered, *GEMSegment_GE11_GlobalPosition_scattered;
  TH2F *GEMRecHit_GE21_GlobalPosition_scattered, *GEMSegment_GE21_GlobalPosition_scattered;
  TH2F *GEMRecHit_GE11_LocalPosition_scattered, *GEMSegment_GE11_LocalPosition_scattered;
  TH2F *GEMRecHit_GE21_LocalPosition_scattered, *GEMSegment_GE21_LocalPosition_scattered;
  TH2F *GEMRecHit_GE11_odd_XZplane, *GEMRecHit_GE11_even_XZplane;
  TH2F *GEMRecHit_GE21_odd_XZplane, *GEMRecHit_GE21_even_XZplane;

  TH1F *Nevents_h, *N_GEMMuon_dist_h;

  //==========================
  //==== Association by hits
  //==========================

  //==== Tracking Paraticle (Eff denominator)
  TH1F *TPMuon_Eta, *TPMuon_Pt, *TPMuon_Phi;
  //==== Reco Tracks (Fake denominator)
  TH1F *HitsGEMMuon_Eta, *HitsGEMMuon_Pt, *HitsGEMMuon_Phi;
  //==== Eff/Fake numerator
  TH1F *HitsMatchedGEMMuon_Eta, *HitsMatchedGEMMuon_Pt, *HitsMatchedGEMMuon_Phi;
  TH1F *HitsBkgMatchedGEMMuon_Eta, *HitsBkgMatchedGEMMuon_Pt, *HitsBkgMatchedGEMMuon_Phi;
  TH1F *HitsUnmatchedGEMMuon_Eta, *HitsUnmatchedGEMMuon_Pt, *HitsUnmatchedGEMMuon_Phi;

  //=====================
  //==== Matching study
  //=====================
  
  std::map< double, TH1F* > HitsPullX_Eta, HitsPullX_Pt, HitsPullX_Phi,
                            HitsMatchedPullX_Eta, HitsMatchedPullX_Pt, HitsMatchedPullX_Phi,
                            HitsBkgMatchedPullX_Eta, HitsBkgMatchedPullX_Pt, HitsBkgMatchedPullX_Phi,
                            HitsUnmatchedPullX_Eta, HitsUnmatchedPullX_Pt, HitsUnmatchedPullX_Phi;
  std::map< double, TH1F* > HitsDX_Eta, HitsDX_Pt, HitsDX_Phi,
                            HitsMatchedDX_Eta, HitsMatchedDX_Pt, HitsMatchedDX_Phi,
                            HitsBkgMatchedDX_Eta, HitsBkgMatchedDX_Pt, HitsBkgMatchedDX_Phi,
                            HitsUnmatchedDX_Eta, HitsUnmatchedDX_Pt, HitsUnmatchedDX_Phi;
  std::map< double, TH1F* > HitsPullY_Eta, HitsPullY_Pt, HitsPullY_Phi,
                            HitsMatchedPullY_Eta, HitsMatchedPullY_Pt, HitsMatchedPullY_Phi,
                            HitsBkgMatchedPullY_Eta, HitsBkgMatchedPullY_Pt, HitsBkgMatchedPullY_Phi,
                            HitsUnmatchedPullY_Eta, HitsUnmatchedPullY_Pt, HitsUnmatchedPullY_Phi;
  std::map< double, TH1F* > HitsDY_Eta, HitsDY_Pt, HitsDY_Phi,
                            HitsMatchedDY_Eta, HitsMatchedDY_Pt, HitsMatchedDY_Phi,
                            HitsBkgMatchedDY_Eta, HitsBkgMatchedDY_Pt, HitsBkgMatchedDY_Phi,
                            HitsUnmatchedDY_Eta, HitsUnmatchedDY_Pt, HitsUnmatchedDY_Phi;
  std::map< double, TH1F* > HitsDotDir_Eta, HitsDotDir_Pt, HitsDotDir_Phi,
                            HitsMatchedDotDir_Eta, HitsMatchedDotDir_Pt, HitsMatchedDotDir_Phi,
                            HitsBkgMatchedDotDir_Eta, HitsBkgMatchedDotDir_Pt, HitsBkgMatchedDotDir_Phi,
                            HitsUnmatchedDotDir_Eta, HitsUnmatchedDotDir_Pt, HitsUnmatchedDotDir_Phi;

};

GEMMuonAnalyzer::GEMMuonAnalyzer(const edm::ParameterSet& iConfig) 
{
  histoFile = new TFile(iConfig.getParameter<std::string>("HistoFile").c_str(), "recreate");
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  doGeometryStudy = iConfig.getParameter< bool >("doGeometryStudy");
  SampleProcess = iConfig.getParameter< std::string >("SampleProcess");

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
  edm::InputTag trackingParticlesTag ("mix","MergedTrackTruth");
  trackingParticlesToken_ = consumes<TrackingParticleCollection>(trackingParticlesTag);
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

  std::cout<<"Contructor end"<<std::endl;
}



void GEMMuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

  const int n_pt_bin = 19, n_eta_bin = 9;
  double pt_bin[n_pt_bin+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  double eta_bin[n_eta_bin+1] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

  Nevents_h = new TH1F("Nevents_h", "Nevents", 2, 0, 2 );
  N_GEMMuon_dist_h = new TH1F("N_GEMMuon_dist_h", "Nevents", 2000, 0, 2000 );

  //================
  //===== Geometry
  //================

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

  //==========================
  //==== Association by hits
  //==========================

  TPMuon_Eta = new TH1F("TPMuon_Eta", "Muon #eta", n_eta_bin, eta_bin );
  TPMuon_Pt = new TH1F("TPMuon_Pt", "Muon p_{T}", n_pt_bin, pt_bin );
  TPMuon_Phi = new TH1F("TPMuon_Phi", "Muon #phi", 36, -TMath::Pi(), TMath::Pi() );

  HitsGEMMuon_Eta = new TH1F("HitsGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsGEMMuon_Pt  = new TH1F("HitsGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsGEMMuon_Phi = new TH1F("HitsGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsMatchedGEMMuon_Eta = new TH1F("HitsMatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsMatchedGEMMuon_Pt  = new TH1F("HitsMatchedGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsMatchedGEMMuon_Phi = new TH1F("HitsMatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsBkgMatchedGEMMuon_Eta = new TH1F("HitsBkgMatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsBkgMatchedGEMMuon_Pt  = new TH1F("HitsBkgMatchedGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsBkgMatchedGEMMuon_Phi = new TH1F("HitsBkgMatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );
  HitsUnmatchedGEMMuon_Eta = new TH1F("HitsUnmatchedGEMMuon_Eta", "GEMMuon #eta", n_eta_bin, eta_bin );
  HitsUnmatchedGEMMuon_Pt  = new TH1F("HitsUnmatchedGEMMuon_Pt", "GENMuon p_{T}", n_pt_bin, pt_bin );
  HitsUnmatchedGEMMuon_Phi = new TH1F("HitsUnmatchedGEMMuon_Phi", "GEMMuon #phi", 36, -TMath::Pi(), TMath::Pi() );

  //=====================
  //==== Matching study
  //=====================

  for(unsigned int i=0; i<PullXValues.size(); i++){
    double aaa = PullXValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsPullX_Eta[aaa] = new TH1F("HitsPullX_Eta"+saaa, "PullX #eta", n_eta_bin, eta_bin );
    HitsPullX_Pt[aaa]  = new TH1F("HitsPullX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsPullX_Phi[aaa] = new TH1F("HitsPullX_Phi"+saaa, "PullX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsMatchedPullX_Eta[aaa] = new TH1F("HitsMatchedPullX_Eta"+saaa, "PullX #eta", n_eta_bin, eta_bin );
    HitsMatchedPullX_Pt[aaa]  = new TH1F("HitsMatchedPullX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedPullX_Phi[aaa] = new TH1F("HitsMatchedPullX_Phi"+saaa, "PullX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsBkgMatchedPullX_Eta[aaa] = new TH1F("HitsBkgMatchedPullX_Eta"+saaa, "PullX #eta", n_eta_bin, eta_bin );
    HitsBkgMatchedPullX_Pt[aaa]  = new TH1F("HitsBkgMatchedPullX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsBkgMatchedPullX_Phi[aaa] = new TH1F("HitsBkgMatchedPullX_Phi"+saaa, "PullX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedPullX_Eta[aaa] = new TH1F("HitsUnmatchedPullX_Eta"+saaa, "PullX #eta", n_eta_bin, eta_bin );
    HitsUnmatchedPullX_Pt[aaa]  = new TH1F("HitsUnmatchedPullX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedPullX_Phi[aaa] = new TH1F("HitsUnmatchedPullX_Phi"+saaa, "PullX #phi", 36, -TMath::Pi(), TMath::Pi() );
  }
  for(unsigned int i=0; i<DXValues.size(); i++){
    double aaa = DXValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsDX_Eta[aaa] = new TH1F("HitsDX_Eta"+saaa, "DX #eta", n_eta_bin, eta_bin );
    HitsDX_Pt[aaa]  = new TH1F("HitsDX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsDX_Phi[aaa] = new TH1F("HitsDX_Phi"+saaa, "DX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsMatchedDX_Eta[aaa] = new TH1F("HitsMatchedDX_Eta"+saaa, "DX #eta", n_eta_bin, eta_bin );
    HitsMatchedDX_Pt[aaa]  = new TH1F("HitsMatchedDX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedDX_Phi[aaa] = new TH1F("HitsMatchedDX_Phi"+saaa, "DX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsBkgMatchedDX_Eta[aaa] = new TH1F("HitsBkgMatchedDX_Eta"+saaa, "DX #eta", n_eta_bin, eta_bin );
    HitsBkgMatchedDX_Pt[aaa]  = new TH1F("HitsBkgMatchedDX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsBkgMatchedDX_Phi[aaa] = new TH1F("HitsBkgMatchedDX_Phi"+saaa, "DX #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedDX_Eta[aaa] = new TH1F("HitsUnmatchedDX_Eta"+saaa, "DX #eta", n_eta_bin, eta_bin );
    HitsUnmatchedDX_Pt[aaa]  = new TH1F("HitsUnmatchedDX_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedDX_Phi[aaa] = new TH1F("HitsUnmatchedDX_Phi"+saaa, "DX #phi", 36, -TMath::Pi(), TMath::Pi() );
  }
  for(unsigned int i=0; i<PullYValues.size(); i++){
    double aaa = PullYValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsPullY_Eta[aaa] = new TH1F("HitsPullY_Eta"+saaa, "PullY #eta", n_eta_bin, eta_bin );
    HitsPullY_Pt[aaa]  = new TH1F("HitsPullY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsPullY_Phi[aaa] = new TH1F("HitsPullY_Phi"+saaa, "PullY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsMatchedPullY_Eta[aaa] = new TH1F("HitsMatchedPullY_Eta"+saaa, "PullY #eta", n_eta_bin, eta_bin );
    HitsMatchedPullY_Pt[aaa]  = new TH1F("HitsMatchedPullY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedPullY_Phi[aaa] = new TH1F("HitsMatchedPullY_Phi"+saaa, "PullY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsBkgMatchedPullY_Eta[aaa] = new TH1F("HitsBkgMatchedPullY_Eta"+saaa, "PullY #eta", n_eta_bin, eta_bin );
    HitsBkgMatchedPullY_Pt[aaa]  = new TH1F("HitsBkgMatchedPullY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsBkgMatchedPullY_Phi[aaa] = new TH1F("HitsBkgMatchedPullY_Phi"+saaa, "PullY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedPullY_Eta[aaa] = new TH1F("HitsUnmatchedPullY_Eta"+saaa, "PullY #eta", n_eta_bin, eta_bin );
    HitsUnmatchedPullY_Pt[aaa]  = new TH1F("HitsUnmatchedPullY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedPullY_Phi[aaa] = new TH1F("HitsUnmatchedPullY_Phi"+saaa, "PullY #phi", 36, -TMath::Pi(), TMath::Pi() );
  }
  for(unsigned int i=0; i<DYValues.size(); i++){
    double aaa = DYValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsDY_Eta[aaa] = new TH1F("HitsDY_Eta"+saaa, "DY #eta", n_eta_bin, eta_bin );
    HitsDY_Pt[aaa]  = new TH1F("HitsDY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsDY_Phi[aaa] = new TH1F("HitsDY_Phi"+saaa, "DY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsMatchedDY_Eta[aaa] = new TH1F("HitsMatchedDY_Eta"+saaa, "DY #eta", n_eta_bin, eta_bin );
    HitsMatchedDY_Pt[aaa]  = new TH1F("HitsMatchedDY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedDY_Phi[aaa] = new TH1F("HitsMatchedDY_Phi"+saaa, "DY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsBkgMatchedDY_Eta[aaa] = new TH1F("HitsBkgMatchedDY_Eta"+saaa, "DY #eta", n_eta_bin, eta_bin );
    HitsBkgMatchedDY_Pt[aaa]  = new TH1F("HitsBkgMatchedDY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsBkgMatchedDY_Phi[aaa] = new TH1F("HitsBkgMatchedDY_Phi"+saaa, "DY #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedDY_Eta[aaa] = new TH1F("HitsUnmatchedDY_Eta"+saaa, "DY #eta", n_eta_bin, eta_bin );
    HitsUnmatchedDY_Pt[aaa]  = new TH1F("HitsUnmatchedDY_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedDY_Phi[aaa] = new TH1F("HitsUnmatchedDY_Phi"+saaa, "DY #phi", 36, -TMath::Pi(), TMath::Pi() );
  }
  for(unsigned int i=0; i<DotDirValues.size(); i++){
    double aaa = DotDirValues.at(i);
    TString saaa = "_"+DoubleToString(aaa);
    HitsDotDir_Eta[aaa] = new TH1F("HitsDotDir_Eta"+saaa, "DotDir #eta", n_eta_bin, eta_bin );
    HitsDotDir_Pt[aaa]  = new TH1F("HitsDotDir_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsDotDir_Phi[aaa] = new TH1F("HitsDotDir_Phi"+saaa, "DotDir #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsMatchedDotDir_Eta[aaa] = new TH1F("HitsMatchedDotDir_Eta"+saaa, "DotDir #eta", n_eta_bin, eta_bin );
    HitsMatchedDotDir_Pt[aaa]  = new TH1F("HitsMatchedDotDir_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsMatchedDotDir_Phi[aaa] = new TH1F("HitsMatchedDotDir_Phi"+saaa, "DotDir #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsBkgMatchedDotDir_Eta[aaa] = new TH1F("HitsBkgMatchedDotDir_Eta"+saaa, "DotDir #eta", n_eta_bin, eta_bin );
    HitsBkgMatchedDotDir_Pt[aaa]  = new TH1F("HitsBkgMatchedDotDir_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsBkgMatchedDotDir_Phi[aaa] = new TH1F("HitsBkgMatchedDotDir_Phi"+saaa, "DotDir #phi", 36, -TMath::Pi(), TMath::Pi() );
    HitsUnmatchedDotDir_Eta[aaa] = new TH1F("HitsUnmatchedDotDir_Eta"+saaa, "DotDir #eta", n_eta_bin, eta_bin );
    HitsUnmatchedDotDir_Pt[aaa]  = new TH1F("HitsUnmatchedDotDir_Pt"+saaa, "GENMuon p_{T}", n_pt_bin, pt_bin );
    HitsUnmatchedDotDir_Phi[aaa] = new TH1F("HitsUnmatchedDotDir_Phi"+saaa, "DotDir #phi", 36, -TMath::Pi(), TMath::Pi() );
  }


}


GEMMuonAnalyzer::~GEMMuonAnalyzer(){}

void
GEMMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  Nevents_h->Fill(1);
  using namespace edm;
  using namespace reco;

  Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByToken(trackingParticlesToken_, trackingParticles);

  iSetup.get<MuonGeometryRecord>().get(gemGeom);
  edm::Handle<GEMRecHitCollection> gemRecHitCollection;
  iEvent.getByToken(GEMRecHit_Token, gemRecHitCollection);
  edm::Handle<GEMSegmentCollection> gemSegmentCollection;
  iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);

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

  if(UseAssociators) {

    edm::LogVerbatim("GEMMuonAnalyzer") << std::endl << "==== Associator ====";
 
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

        bool SignalMuon = false;
        int motherId(0), grandmaId(0), greatgrandmaId(0);

        //====  Z->MuMu sample
        if(SampleProcess == "ZMM"){
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
        }
        //==== MuonGun sample
        else{
          if(tp->status() != -99){ // Pythia8 gen status : home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
            if( tp->genParticles().size()>0 && (*tp->genParticle_begin())->numberOfMothers()>0 ) {
              motherId = abs( (*tp->genParticle_begin())->mother()->pdgId() );
              if( (*tp->genParticle_begin())->mother()->numberOfMothers()>0 ) {
                grandmaId = abs( (*tp->genParticle_begin())->mother()->mother()->pdgId() );
                if( (*tp->genParticle_begin())->mother()->mother()->numberOfMothers()>0 ) {
                  greatgrandmaId = abs( (*tp->genParticle_begin())->mother()->mother()->mother()->pdgId() );
                }
              }
            }
            //std::cout<<"Mother ID = "<<motherid<<std::endl;

            if ( ( (tp->status()==1) && ( (*tp->genParticle_begin())->numberOfMothers()==0 ) ) ) SignalMuon=true;
          }
        }

        //std::cout
        //<< i << '\t'
        //<< tp->status() << '\t'
        //<< motherId << '\t'
        //<< grandmaId << '\t'
        //<< greatgrandmaId << '\t'
        //<< std::endl;

        if(SignalMuon){

          //==== Fill the Denominator
          //==== Should be filled only once (www=0)
          if( www == 0 ){
            if(Pt_5) TPMuon_Eta->Fill(fabs(tp->eta()));
            if(Eta_1p6_2p4) TPMuon_Pt->Fill(tp->pt());
            if(Pt_5 && Eta_1p6_2p4) TPMuon_Phi->Fill(tp->phi());
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
                  if(Eta_1p6_2p4) HitsMatchedGEMMuon_Pt->Fill(tpr->pt());
                  if(Pt_5 && Eta_1p6_2p4) HitsMatchedGEMMuon_Phi->Fill(tpr->phi());
                }
                if(label[www].find("PullXScan") != std::string::npos){
                  double this_cut = PullXValues.at(www_PullX);
                  if(Pt_5) HitsMatchedPullX_Eta[this_cut]->Fill(fabs(tpr->eta()));
                  if(Eta_1p6_2p4) HitsMatchedPullX_Pt[this_cut]->Fill(tpr->pt());
                  if(Pt_5 && Eta_1p6_2p4) HitsMatchedPullX_Phi[this_cut]->Fill(tpr->phi());
                }
                if(label[www].find("DXScan") != std::string::npos){
                  double this_cut = DXValues.at(www_DX);
                  if(Pt_5) HitsMatchedDX_Eta[this_cut]->Fill(fabs(tpr->eta()));
                  if(Eta_1p6_2p4) HitsMatchedDX_Pt[this_cut]->Fill(tpr->pt());
                  if(Pt_5 && Eta_1p6_2p4) HitsMatchedDX_Phi[this_cut]->Fill(tpr->phi());
                }
                if(label[www].find("PullYScan") != std::string::npos){
                  double this_cut = PullYValues.at(www_PullY);
                  if(Pt_5) HitsMatchedPullY_Eta[this_cut]->Fill(fabs(tpr->eta()));
                  if(Eta_1p6_2p4) HitsMatchedPullY_Pt[this_cut]->Fill(tpr->pt());
                  if(Pt_5 && Eta_1p6_2p4) HitsMatchedPullY_Phi[this_cut]->Fill(tpr->phi());
                }
                if(label[www].find("DYScan") != std::string::npos){
                  double this_cut = DYValues.at(www_DY);
                  if(Pt_5) HitsMatchedDY_Eta[this_cut]->Fill(fabs(tpr->eta()));
                  if(Eta_1p6_2p4) HitsMatchedDY_Pt[this_cut]->Fill(tpr->pt());
                  if(Pt_5 && Eta_1p6_2p4) HitsMatchedDY_Phi[this_cut]->Fill(tpr->phi());
                }
                if(label[www].find("DotDirScan") != std::string::npos){
                  double this_cut = DotDirValues.at(www_DotDir);
                  if(Pt_5) HitsMatchedDotDir_Eta[this_cut]->Fill(fabs(tpr->eta()));
                  if(Eta_1p6_2p4) HitsMatchedDotDir_Pt[this_cut]->Fill(tpr->pt());
                  if(Pt_5 && Eta_1p6_2p4) HitsMatchedDotDir_Phi[this_cut]->Fill(tpr->phi());
                }


              } //==== Matched
            } //==== Matched
          } //==== Matched



        } //==== This TP is Signal Muon





      } // END TrackingParticle Loop

      //=================
      //==== Fake study
      //=================

      int n_this_tracks = 0;

      //==== loop over (GEMMuon/RecoMuon/...) tracks
      for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
        //std::cout << i << "th trackCollection iterator" << std::endl;
        RefToBase<Track> track(trackCollection, i);

        bool Eta_1p6_2p4 = fabs(track->eta()) > 1.6 && fabs(track->eta()) < 2.4;
        bool Pt_5 = track->pt() > 5;

        n_this_tracks++;

        //==== denominator distributions
        //==== count total reco tracks
        if(label[www]=="gemMuonSel"){
          if( Eta_1p6_2p4 ) HitsGEMMuon_Pt->Fill(track->pt());
          if( Pt_5 ) HitsGEMMuon_Eta->Fill(fabs(track->eta()));
          if( Pt_5 && Eta_1p6_2p4 ) HitsGEMMuon_Phi->Fill(track->phi());
        }
        if(label[www].find("PullXScan") != std::string::npos){
          double this_cut = PullXValues.at(www_PullX);
          if( Eta_1p6_2p4 ) HitsPullX_Pt[this_cut]->Fill(track->pt());
          if( Pt_5 ) HitsPullX_Eta[this_cut]->Fill(fabs(track->eta()));
          if( Pt_5 && Eta_1p6_2p4 ) HitsPullX_Phi[this_cut]->Fill(track->phi());
        }
        if(label[www].find("DXScan") != std::string::npos){
          double this_cut = DXValues.at(www_DX);
          if( Eta_1p6_2p4 ) HitsDX_Pt[this_cut]->Fill(track->pt());
          if( Pt_5 ) HitsDX_Eta[this_cut]->Fill(fabs(track->eta()));
          if( Pt_5 && Eta_1p6_2p4 ) HitsDX_Phi[this_cut]->Fill(track->phi());
        }
        if(label[www].find("PullYScan") != std::string::npos){
          double this_cut = PullYValues.at(www_PullY);
          if( Eta_1p6_2p4 ) HitsPullY_Pt[this_cut]->Fill(track->pt());
          if( Pt_5 ) HitsPullY_Eta[this_cut]->Fill(fabs(track->eta()));
          if( Pt_5 && Eta_1p6_2p4 ) HitsPullY_Phi[this_cut]->Fill(track->phi());
        }
        if(label[www].find("DYScan") != std::string::npos){
          double this_cut = DYValues.at(www_DY);
          if( Eta_1p6_2p4 ) HitsDY_Pt[this_cut]->Fill(track->pt());
          if( Pt_5 ) HitsDY_Eta[this_cut]->Fill(fabs(track->eta()));
          if( Pt_5 && Eta_1p6_2p4 ) HitsDY_Phi[this_cut]->Fill(track->phi());
        }
        if(label[www].find("DotDirScan") != std::string::npos){
          double this_cut = DotDirValues.at(www_DotDir);
          if( Eta_1p6_2p4 ) HitsDotDir_Pt[this_cut]->Fill(track->pt());
          if( Pt_5 ) HitsDotDir_Eta[this_cut]->Fill(fabs(track->eta()));
          if( Pt_5 && Eta_1p6_2p4 ) HitsDotDir_Phi[this_cut]->Fill(track->phi());
        }

        //==== Check if the track is associated to any gen particle
        bool isFake = false;
        bool isBkg = false;

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

            //====  Z->MuMu sample
            if(SampleProcess == "ZMM"){
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
            }
            //==== MuonGun sample
            else{
              if(tpr->status() != -99){ // Pythia8 gen status : home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
                //int motherid=-1;
                if ((*tpr->genParticle_begin())->numberOfMothers()>0)  {
                  if ((*tpr->genParticle_begin())->mother()->numberOfMothers()>0){
                    //motherid=(*tpr->genParticle_begin())->mother()->mother()->pdgId();
                  }
                }
                //std::cout<<"Mother ID = "<<motherid<<std::endl;

                if ( ( (tpr->status()==1) && ( (*tpr->genParticle_begin())->numberOfMothers()==0 ) ) )  SignalMuon=true;
              }
            }

            if( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) && SignalMuon ) {
              edm::LogVerbatim("GEMMuonAnalyzer") << "This RecoTrack is matcehd to SIGNAL Muon";
            }
            else{
              edm::LogVerbatim("GEMMuonAnalyzer") << "This RecoTrack is NOT matcehd to SIGNAL Muon => Fake";
              isFake = true;
            }

            if( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) && !SignalMuon ){
              edm::LogVerbatim("GEMMuonAnalyzer") << "This RecoTrack is matcehd to BKG Muon";
              isBkg = true;
            }



          }
        }

        if( isFake ){
          if(label[www]=="gemMuonSel"){
             if(Pt_5) HitsUnmatchedGEMMuon_Eta->Fill(fabs(track->eta()));
             if(Eta_1p6_2p4) HitsUnmatchedGEMMuon_Pt->Fill(track->pt());
             if(Pt_5 && Eta_1p6_2p4) HitsUnmatchedGEMMuon_Phi->Fill(track->phi());
          }
          if(label[www].find("PullXScan") != std::string::npos){
            double this_cut = PullXValues.at(www_PullX);
            if(Pt_5) HitsUnmatchedPullX_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsUnmatchedPullX_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsUnmatchedPullX_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("DXScan") != std::string::npos){
            double this_cut = DXValues.at(www_DX);
            if(Pt_5) HitsUnmatchedDX_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsUnmatchedDX_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsUnmatchedDX_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("PullYScan") != std::string::npos){
            double this_cut = PullYValues.at(www_PullY);
            if(Pt_5) HitsUnmatchedPullY_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsUnmatchedPullY_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsUnmatchedPullY_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("DYScan") != std::string::npos){
            double this_cut = DYValues.at(www_DY);
            if(Pt_5) HitsUnmatchedDY_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsUnmatchedDY_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsUnmatchedDY_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("DotDirScan") != std::string::npos){
            double this_cut = DotDirValues.at(www_DotDir);
            if(Pt_5) HitsUnmatchedDotDir_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsUnmatchedDotDir_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsUnmatchedDotDir_Phi[this_cut]->Fill(track->phi());
          }

        } //==== END if(isFake)

        if( isBkg ){
          if(label[www]=="gemMuonSel"){
             if(Pt_5) HitsBkgMatchedGEMMuon_Eta->Fill(fabs(track->eta()));
             if(Eta_1p6_2p4) HitsBkgMatchedGEMMuon_Pt->Fill(track->pt());
             if(Pt_5 && Eta_1p6_2p4) HitsBkgMatchedGEMMuon_Phi->Fill(track->phi());
          }
          if(label[www].find("PullXScan") != std::string::npos){
            double this_cut = PullXValues.at(www_PullX);
            if(Pt_5) HitsBkgMatchedPullX_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsBkgMatchedPullX_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsBkgMatchedPullX_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("DXScan") != std::string::npos){
            double this_cut = DXValues.at(www_DX);
            if(Pt_5) HitsBkgMatchedDX_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsBkgMatchedDX_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsBkgMatchedDX_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("PullYScan") != std::string::npos){
            double this_cut = PullYValues.at(www_PullY);
            if(Pt_5) HitsBkgMatchedPullY_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsBkgMatchedPullY_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsBkgMatchedPullY_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("DYScan") != std::string::npos){
            double this_cut = DYValues.at(www_DY);
            if(Pt_5) HitsBkgMatchedDY_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsBkgMatchedDY_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsBkgMatchedDY_Phi[this_cut]->Fill(track->phi());
          }
          if(label[www].find("DotDirScan") != std::string::npos){
            double this_cut = DotDirValues.at(www_DotDir);
            if(Pt_5) HitsBkgMatchedDotDir_Eta[this_cut]->Fill(fabs(track->eta()));
            if(Eta_1p6_2p4) HitsBkgMatchedDotDir_Pt[this_cut]->Fill(track->pt());
            if(Pt_5 && Eta_1p6_2p4) HitsBkgMatchedDotDir_Phi[this_cut]->Fill(track->phi());
          }

        } //==== END if(isBkg)


      } // END track loop

      //==== # of track distribution

      if(label[www]=="gemMuonSel"){
        N_GEMMuon_dist_h->Fill(n_this_tracks);
      }
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

  N_GEMMuon_dist_h->Write();
 
  //==========================
  //==== Association by hits
  //==========================

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
  HitsBkgMatchedGEMMuon_Eta->Write();
  HitsBkgMatchedGEMMuon_Pt->Write();
  HitsBkgMatchedGEMMuon_Phi->Write();
  //==== Singal Efficeicny
  TEfficiency* HitsEff_GEMMuon_Eta = new TEfficiency(*HitsMatchedGEMMuon_Eta, *TPMuon_Eta);
  TEfficiency* HitsEff_GEMMuon_Pt = new TEfficiency(*HitsMatchedGEMMuon_Pt, *TPMuon_Pt);
  TEfficiency* HitsEff_GEMMuon_Phi = new TEfficiency(*HitsMatchedGEMMuon_Phi, *TPMuon_Phi);
  HitsEff_GEMMuon_Eta->SetName("HitsEff_GEMMuon_Eta");
  HitsEff_GEMMuon_Pt->SetName("HitsEff_GEMMuon_Pt");
  HitsEff_GEMMuon_Phi->SetName("HitsEff_GEMMuon_Phi");
  HitsEff_GEMMuon_Eta->Write();
  HitsEff_GEMMuon_Pt->Write();
  HitsEff_GEMMuon_Phi->Write();
  //==== Fake
  HitsUnmatchedGEMMuon_Eta->Write();
  HitsUnmatchedGEMMuon_Pt->Write();
  HitsUnmatchedGEMMuon_Phi->Write();

  //===============
  //==== Geometry
  //===============

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
    HitsPullX_Eta[aaa]->Write();
    HitsPullX_Pt[aaa]->Write();
    HitsPullX_Phi[aaa]->Write();
    HitsMatchedPullX_Eta[aaa]->Write();
    HitsMatchedPullX_Pt[aaa]->Write();
    HitsMatchedPullX_Phi[aaa]->Write();
    HitsBkgMatchedPullX_Eta[aaa]->Write();
    HitsBkgMatchedPullX_Pt[aaa]->Write();
    HitsBkgMatchedPullX_Phi[aaa]->Write();
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

  }
  hist_PullXValues->Write();

  TH1F *hist_DXValues = new TH1F("DXValues", "", DXValues.size(), 0, DXValues.size());
  for(unsigned int i=0; i<DXValues.size(); i++){
    double aaa = DXValues.at(i);
    hist_DXValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsDX_Eta[aaa]->Write();
    HitsDX_Pt[aaa]->Write();
    HitsDX_Phi[aaa]->Write();
    HitsMatchedDX_Eta[aaa]->Write();
    HitsMatchedDX_Pt[aaa]->Write();
    HitsMatchedDX_Phi[aaa]->Write();
    HitsBkgMatchedDX_Eta[aaa]->Write();
    HitsBkgMatchedDX_Pt[aaa]->Write();
    HitsBkgMatchedDX_Phi[aaa]->Write();
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

  }
  hist_DXValues->Write();

  TH1F *hist_PullYValues = new TH1F("PullYValues", "", PullYValues.size(), 0, PullYValues.size());
  for(unsigned int i=0; i<PullYValues.size(); i++){
    double aaa = PullYValues.at(i);
    hist_PullYValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsPullY_Eta[aaa]->Write();
    HitsPullY_Pt[aaa]->Write();
    HitsPullY_Phi[aaa]->Write();
    HitsMatchedPullY_Eta[aaa]->Write();
    HitsMatchedPullY_Pt[aaa]->Write();
    HitsMatchedPullY_Phi[aaa]->Write();
    HitsBkgMatchedPullY_Eta[aaa]->Write();
    HitsBkgMatchedPullY_Pt[aaa]->Write();
    HitsBkgMatchedPullY_Phi[aaa]->Write();
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

  }
  hist_PullYValues->Write();

  TH1F *hist_DYValues = new TH1F("DYValues", "", DYValues.size(), 0, DYValues.size());
  for(unsigned int i=0; i<DYValues.size(); i++){
    double aaa = DYValues.at(i);
    hist_DYValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsDY_Eta[aaa]->Write();
    HitsDY_Pt[aaa]->Write();
    HitsDY_Phi[aaa]->Write();
    HitsMatchedDY_Eta[aaa]->Write();
    HitsMatchedDY_Pt[aaa]->Write();
    HitsMatchedDY_Phi[aaa]->Write();
    HitsBkgMatchedDY_Eta[aaa]->Write();
    HitsBkgMatchedDY_Pt[aaa]->Write();
    HitsBkgMatchedDY_Phi[aaa]->Write();
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

  }
  hist_DYValues->Write();

  TH1F *hist_DotDirValues = new TH1F("DotDirValues", "", DotDirValues.size(), 0, DotDirValues.size());
  for(unsigned int i=0; i<DotDirValues.size(); i++){
    double aaa = DotDirValues.at(i);
    hist_DotDirValues->SetBinContent(i+1, aaa);
    TString saaa = DoubleToString(aaa);
    HitsDotDir_Eta[aaa]->Write();
    HitsDotDir_Pt[aaa]->Write();
    HitsDotDir_Phi[aaa]->Write();
    HitsMatchedDotDir_Eta[aaa]->Write();
    HitsMatchedDotDir_Pt[aaa]->Write();
    HitsMatchedDotDir_Phi[aaa]->Write();
    HitsBkgMatchedDotDir_Eta[aaa]->Write();
    HitsBkgMatchedDotDir_Pt[aaa]->Write();
    HitsBkgMatchedDotDir_Phi[aaa]->Write();
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

  }
  hist_DotDirValues->Write();

  histoFile->cd();

}

std::string GEMMuonAnalyzer::DoubleToString(double dd){
  std::ostringstream os;
  os << dd;
  return os.str();
}

DEFINE_FWK_MODULE(GEMMuonAnalyzer);
