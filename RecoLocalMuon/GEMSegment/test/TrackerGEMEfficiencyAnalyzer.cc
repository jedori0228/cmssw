// -*- C++ -*-
//
// Package:    TrackerGEMEfficiencyAnalyzer
// Class:      TrackerGEMEfficiencyAnalyzer
// 
/**\class TrackerGEMEfficiencyAnalyzer TrackerGEMEfficiencyAnalyzer.cc MyAnalyzers/TrackerGEMEfficiencyAnalyzer/src/TrackerGEMEfficiencyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>


// root include files
#include "TFile.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/GEMRecHit/interface/GEMSegmentCollection.h>
#include <DataFormats/MuonReco/interface/MuonSelectors.h>

#include <DataFormats/MuonReco/interface/MuonChamberMatch.h>
#include <DataFormats/MuonReco/interface/MuonSegmentMatch.h>

#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

//
// class declaration
//

const int n_pt_bin = 19, n_eta_bin = 9;
double pt_bin[n_pt_bin+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
double eta_bin[n_eta_bin+1] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

class TrackerGEMEfficiencyAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackerGEMEfficiencyAnalyzer(const edm::ParameterSet&);
      ~TrackerGEMEfficiencyAnalyzer();



   private:
      virtual int FindWhichBin(double value, double *bin, int n_bin);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void MakeEfficiency(TString hname, int nbins, double xbins[]);
      virtual void MakeEfficiency(TString hname, int nbins, double xmin, double xmax);
      virtual TEfficiency* GetEfficiency(TString hname);
      virtual void FillEfficiency(TString histname, bool pass, double value, int nbins, double xbins[]);
      virtual void FillEfficiency(TString histname, bool pass, double value, int nbins, double xmin, double xmax);
      virtual void WriteEfficiencies();
      virtual void MakeHistograms(TString hname, int nbins, double xbins[]);
      virtual void MakeHistograms(TString hname, int nbins, double xmin, double xmax);
      virtual TH1* GetHist(TString hname);
      virtual void FillHist(TString histname, double value, int nbins, double xbins[]);
      virtual void FillHist(TString histname, double value, int nbins, double xmin, double xmax);
      virtual void WriteHists();

      // ----------member data ---------------------------

  edm::ESHandle<GEMGeometry> gemGeom;

  edm::EDGetTokenT<reco::GenParticleCollection> GENParticle_Token;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  edm::EDGetTokenT<reco::VertexCollection> vertex_Token;

  std::string rootFileName;
  std::unique_ptr<TFile> outputfile;

  map<TString, TEfficiency*> mapefficiency;
  map<TString, TH1*> maphist;

/*
  //2d plot//
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHitsToken_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegmentsToken_;
  int n_rechit_st[3], n_segment_st[3];
  TH2F *RecHit_global_station_1, *RecHit_global_station_2, *RecHit_global_station_3;
  TH2F *Segment_RecHit_global_station_1, *Segment_RecHit_global_station_2, *Segment_RecHit_global_station_3;
  TH2F *Segment_FirstHit_global_station_1, *Segment_FirstHit_global_station_2, *Segment_FirstHit_global_station_3;
*/
};

//
// constants, enums and typedefs
//
// constructors and destructor
//
TrackerGEMEfficiencyAnalyzer::TrackerGEMEfficiencyAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  GENParticle_Token = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
  GEMSegment_Token = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));
  vertex_Token = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));

  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
  outputfile.reset(TFile::Open(rootFileName.c_str(), "RECREATE"));

/*
  //2d plot//
  gemRecHitsToken_ = consumes<GEMRecHitCollection >(edm::InputTag("gemRecHits"));
  gemSegmentsToken_ = consumes<GEMSegmentCollection >(edm::InputTag("gemSegments"));
  n_rechit_st[0] = 0; n_rechit_st[1] = 0; n_rechit_st[2] = 0;
  n_segment_st[0] = 0; n_segment_st[1] = 0; n_segment_st[2] = 0;
  RecHit_global_station_1 = new TH2F("RecHit_global_station_1", "", 1000, -500, 500, 1000, -500, 500);
  RecHit_global_station_2 = new TH2F("RecHit_global_station_2", "", 1000, -500, 500, 1000, -500, 500);
  RecHit_global_station_3 = new TH2F("RecHit_global_station_3", "", 1000, -500, 500, 1000, -500, 500);
  Segment_RecHit_global_station_1 = new TH2F("Segment_RecHit_global_station_1", "", 1000, -500, 500, 1000, -500, 500);
  Segment_RecHit_global_station_2 = new TH2F("Segment_RecHit_global_station_2", "", 1000, -500, 500, 1000, -500, 500);
  Segment_RecHit_global_station_3 = new TH2F("Segment_RecHit_global_station_3", "", 1000, -500, 500, 1000, -500, 500);
  Segment_FirstHit_global_station_1 = new TH2F("Segment_FirstHit_global_station_1", "", 1000, -500, 500, 1000, -500, 500);
  Segment_FirstHit_global_station_2 = new TH2F("Segment_FirstHit_global_station_2", "", 1000, -500, 500, 1000, -500, 500);
  Segment_FirstHit_global_station_3 = new TH2F("Segment_FirstHit_global_station_3", "", 1000, -500, 500, 1000, -500, 500);
*/
}


TrackerGEMEfficiencyAnalyzer::~TrackerGEMEfficiencyAnalyzer()
{

  outputfile->cd();

  WriteEfficiencies();
  WriteHists();

/*
  //2d plot//
  RecHit_global_station_1->Write();
  RecHit_global_station_2->Write();
  RecHit_global_station_3->Write();
  Segment_RecHit_global_station_1->Write();
  Segment_RecHit_global_station_2->Write();
  Segment_RecHit_global_station_3->Write();
  Segment_FirstHit_global_station_1->Write();
  Segment_FirstHit_global_station_2->Write();
  Segment_FirstHit_global_station_3->Write();
*/
  outputfile->Close();
/*
  //2d plot//
  std::cout << "[GEMRecHit]" << std::endl;
  for(int i=0; i<3; i++) std::cout << " station " << i+1 << " : " << n_rechit_st[i] << std::endl;
  std::cout << "[GEMSegment]" << std::endl;
  for(int i=0; i<3; i++) std::cout << " station " << i+1 << " : " << n_segment_st[i] << std::endl;
*/
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackerGEMEfficiencyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
/*
  //2d plot//
  using namespace edm;
  using namespace reco;
  using namespace std;
  iSetup.get<MuonGeometryRecord>().get(gemGeom);

  Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHitsToken_,gemRecHits);

  Handle<GEMSegmentCollection> gemSegments;
  iEvent.getByToken(gemSegmentsToken_,gemSegments);

  for (auto thisRecHit = gemRecHits->begin(); thisRecHit != gemRecHits->end();
       ++thisRecHit){
    //std::cout << "RecHit station = " << thisRecHit->gemId().station() << std::endl;
    n_rechit_st[thisRecHit->gemId().station()-1]++;

    GEMDetId id = thisRecHit->gemId();
    //std::cout << "id = " << id << std::endl;

    //GEMRecHit  local position is w.r.t. roll
    //GEMSegment local position is w.r.t. chamber
    auto roll = gemGeom->etaPartition(id);
    if(!roll){
      std::cout << "no roll for id = " << id << std::endl;
      continue; // FIXME no GEMEtaPartitions for some rolls..
    }
    LocalPoint thisLocalPos = thisRecHit->localPosition();

    //std::cout
    //<< "station = " << thisRecHit->gemId().station() << " : " << endl
    //<< "(" << thisLocalPos.x() << ", " << thisLocalPos.y() << ", " << thisLocalPos.z() << ")" << std::endl;

    GlobalPoint thisGlobalPos(roll->toGlobal(thisLocalPos));
    //std::cout << " global position obtained" << std::endl;


    if(thisRecHit->gemId().station() == 1){
      RecHit_global_station_1->Fill(thisGlobalPos.x(), thisGlobalPos.y());
      //std::cout << "(" << thisLocalPos.x() << ", " << thisLocalPos.y() << ", " << thisLocalPos.z() << ")" << std::endl;
    }
    else if(thisRecHit->gemId().station() == 2){
      RecHit_global_station_2->Fill(thisGlobalPos.x(), thisGlobalPos.y());
    }
    else if(thisRecHit->gemId().station() == 3){
      RecHit_global_station_3->Fill(thisGlobalPos.x(), thisGlobalPos.y());
    }
    else{}
  }

  for (auto thisSegment = gemSegments->begin(); thisSegment != gemSegments->end();
       ++thisSegment){
    //std::cout << "A Segment with :" << std::endl;
    //LocalPoint thisSegmentLocalPos = thisSegment->localPosition();
    //std::cout << "(" << thisSegmentLocalPos.x() << ", " << thisSegmentLocalPos.y() << ", " << thisSegmentLocalPos.z() << ")" << std::endl;
    int i=0;
    for(auto thisRecHit =  thisSegment->specificRecHits().begin();
             thisRecHit != thisSegment->specificRecHits().end();
        ++thisRecHit, i++){
      //std::cout << " station = " << thisSegment->specificRecHits()[i].gemId().station() << std::endl;
      n_segment_st[thisSegment->specificRecHits()[i].gemId().station()-1]++;

      GEMDetId id = thisRecHit->gemId();
      auto roll = gemGeom->etaPartition(id);
      if(!roll){
        std::cout << "no roll for id = " << id << std::endl;
        continue; // FIXME no GEMEtaPartitions for some rolls..
      }
      LocalPoint thisLocalPos = thisRecHit->localPosition();
      GlobalPoint thisGlobalPos(roll->toGlobal(thisLocalPos));

      if(thisRecHit->gemId().station() == 1){
        Segment_RecHit_global_station_1->Fill(thisGlobalPos.x(), thisGlobalPos.y());
      }
      else if(thisRecHit->gemId().station() == 2){
        Segment_RecHit_global_station_2->Fill(thisGlobalPos.x(), thisGlobalPos.y());
      }
      else if(thisRecHit->gemId().station() == 3){
        Segment_RecHit_global_station_3->Fill(thisGlobalPos.x(), thisGlobalPos.y());
      }
      else{}

      if( i == 0 ){
        if(thisRecHit->gemId().station() == 1){
          Segment_FirstHit_global_station_1->Fill(thisGlobalPos.x(), thisGlobalPos.y());
        }
        else if(thisRecHit->gemId().station() == 2){
          Segment_FirstHit_global_station_2->Fill(thisGlobalPos.x(), thisGlobalPos.y());
        }
        else if(thisRecHit->gemId().station() == 3){
          Segment_FirstHit_global_station_3->Fill(thisGlobalPos.x(), thisGlobalPos.y());
        }
        else{}
      }


    }
  }
*/



  iSetup.get<MuonGeometryRecord>().get(gemGeom);

  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(GENParticle_Token, genParticles);
  edm::Handle<GEMSegmentCollection> gemSegmentCollection;
  iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertex_Token, vertices);
  reco::Vertex vertex = vertices->at(0);

  /////////////////////////////////////////////////////
  /// USE HAVE TO USE MUON GUN SAMPLE FOR THIS CODE ///
  /////////////////////////////////////////////////////

  double matching_deltaR = 0.1;

  /// loop over gen particles = muons ///
  for(reco::GenParticleCollection::const_iterator genpart=genParticles->begin(); genpart != genParticles->end(); ++genpart) {
    if(fabs(genpart->pdgId()) != 13) continue;

    TLorentzVector igenP4;
    igenP4.SetPtEtaPhiM(genpart->pt(), genpart->eta(), genpart->phi(), genpart->mass());
    //cout << endl << "[gen]" << endl; igenP4.Print();

    /// loop over RecoMuon ///
    bool RecoMuon_isMatched = false, RecoMuon_isSAMuon = false, RecoMuon_isLooseMuon = false, RecoMuon_isTightMuon = false, 
         RecoMuon_isTMOneStationLoose = false, RecoMuon_isTMOneStationTight = false,
         RecoMuon_isTMOneStationAngLoose = false, RecoMuon_isTMOneStationAngTight = false,
         RecoMuon_isTMLastStationLoose = false, RecoMuon_isTMLastStationTight = false,
         RecoMuon_isTMLastStationAngLoose = false, RecoMuon_isTMLastStationAngTight = false,
         RecoMuon_isGEMMuon = false, RecoMuon_isTrackerMuon = false,
         RecoMuon_isGEMMuon_or_isTrackerMuon = false, RecoMuon_isGEMMuon_and_isTrackerMuon = false; 
    double deltaR_reco_temp = matching_deltaR;
    //const reco::Muon* matched_recoMuon = NULL;
    for(reco::MuonCollection::const_iterator recomuon=recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon) { 
      TLorentzVector irecoP4;
      irecoP4.SetPtEtaPhiM(recomuon->pt(), recomuon->eta(), recomuon->phi(), recomuon->mass());
      //cout << "[reco]" << endl; irecoP4.Print();
      if( igenP4.DeltaR(irecoP4) < deltaR_reco_temp ){
        RecoMuon_isMatched = true;
        //matched_recoMuon = &(*recomuon);
        if( recomuon->globalTrack().isNonnull() ){
          RecoMuon_isSAMuon = (recomuon->isStandAloneMuon())
                              //& (recomuon->numberOfMatchedStations() >= 1) 
                              & (recomuon->outerTrack()->hitPattern().muonStationsWithValidHits() > 1)
                              & ( std::abs(recomuon->time().timeAtIpInOut) < ( 12.5 + std::abs( recomuon->time().timeAtIpInOutErr ) ) );
        }
        RecoMuon_isLooseMuon = muon::isLooseMuon(*recomuon);
        RecoMuon_isTightMuon = muon::isTightMuon(*recomuon, vertex);
        //bool ip = fabs(recomuon->muonBestTrack()->dxy(vertex.position())) < 0.2 && fabs(recomuon->muonBestTrack()->dz(vertex.position())) < 0.5;
        
        //cout
        //<< "fabs(recomuon->muonBestTrack()->dxy(vertex.position())) = " << fabs(recomuon->muonBestTrack()->dxy(vertex.position())) << endl
        //<< "fabs(recomuon->muonBestTrack()->dz(vertex.position())) = " << fabs(recomuon->muonBestTrack()->dz(vertex.position())) << endl << endl;
        
        RecoMuon_isTMOneStationLoose = muon::isGoodMuon(*recomuon, muon::TMOneStationLoose);
        RecoMuon_isTMOneStationTight = muon::isGoodMuon(*recomuon, muon::TMOneStationTight);
        RecoMuon_isTMOneStationAngLoose = muon::isGoodMuon(*recomuon, muon::TMOneStationAngLoose);
        RecoMuon_isTMOneStationAngTight = muon::isGoodMuon(*recomuon, muon::TMOneStationAngTight);
        RecoMuon_isTMLastStationLoose = muon::isGoodMuon(*recomuon, muon::TMLastStationLoose);
        RecoMuon_isTMLastStationTight = muon::isGoodMuon(*recomuon, muon::TMLastStationTight);
        RecoMuon_isTMLastStationAngLoose = muon::isGoodMuon(*recomuon, muon::TMLastStationAngLoose);
        RecoMuon_isTMLastStationAngTight = muon::isGoodMuon(*recomuon, muon::TMLastStationAngTight);
        RecoMuon_isTrackerMuon = recomuon->isTrackerMuon();

        RecoMuon_isGEMMuon = recomuon->isGEMMuon();        
        RecoMuon_isGEMMuon_or_isTrackerMuon = recomuon->isGEMMuon() || recomuon->isTrackerMuon();
        RecoMuon_isGEMMuon_and_isTrackerMuon = recomuon->isGEMMuon() && recomuon->isTrackerMuon();

        deltaR_reco_temp = igenP4.DeltaR(irecoP4);
      }
    }
  
    bool Eta_1p6_2p4 = abs(igenP4.Eta()) > 1.6 && abs(igenP4.Eta()) < 2.4,
         Pt_5 = igenP4.Pt() > 5; 

    /// loop over GEMSegment ///
    bool GEMSegment_isMatched = false;
    double deltaR_gemseg_temp = matching_deltaR;
    for (auto gems = gemSegmentCollection->begin(); gems != gemSegmentCollection->end(); ++gems) {
      GEMDetId id = gems->gemDetId();
      auto chamb = gemGeom->superChamber(id);
      auto segLP = gems->localPosition();
      auto segGP = chamb->toGlobal(segLP);
      TLorentzVector segment;
      segment.SetPtEtaPhiM(1, segGP.eta(), segGP.phi(), 0);
      if(igenP4.DeltaR(segment) < deltaR_gemseg_temp ){
        GEMSegment_isMatched = true;
        deltaR_gemseg_temp = igenP4.DeltaR(segment);
      }
    } // end of GEMSegment loop

    //onebin//
    if(Eta_1p6_2p4 && Pt_5){
      FillEfficiency("SAMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isSAMuon, 0, 1, 0, 1);
      FillEfficiency("LooseMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isLooseMuon, 0, 1, 0, 1);
      FillEfficiency("TightMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isTightMuon, 0, 1, 0, 1);
      FillEfficiency("TMOneStationLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationLoose, 0, 1, 0, 1);
      FillEfficiency("TMOneStationTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationTight, 0, 1, 0, 1);
      FillEfficiency("TMOneStationAngLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose, 0, 1, 0, 1);
      FillEfficiency("TMOneStationAngTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight, 0, 1, 0, 1);
      FillEfficiency("TMLastStationLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationLoose, 0, 1, 0, 1);
      FillEfficiency("TMLastStationTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationTight, 0, 1, 0, 1);
      FillEfficiency("TMLastStationAngLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose, 0, 1, 0, 1);
      FillEfficiency("TMLastStationAngTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight, 0, 1, 0, 1);
      FillEfficiency("GEMMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isGEMMuon, 0, 1, 0, 1);
      FillEfficiency("TrackerMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isTrackerMuon, 0, 1, 0, 1);
      FillEfficiency("GEMMuon_or_TrackerMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isGEMMuon_or_isTrackerMuon, 0, 1, 0, 1);
      FillEfficiency("GEMMuon_and_TrackerMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isGEMMuon_and_isTrackerMuon, 0, 1, 0, 1);
      FillEfficiency("GEMSegment_eff_onebin", GEMSegment_isMatched, 0, 1, 0, 1);
      FillEfficiency("SAMuon_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isSAMuon) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("LooseMuon_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isLooseMuon) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TightMuon_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTightMuon) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMOneStationLoose_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMOneStationLoose) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMOneStationTight_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMOneStationTight) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMOneStationAngLoose_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMOneStationAngTight_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMLastStationLoose_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMLastStationLoose) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMLastStationTight_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMLastStationTight) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMLastStationAngLoose_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose) | RecoMuon_isGEMMuon , 0, 1, 0, 1);
      FillEfficiency("TMLastStationAngTight_or_GEMMuon_eff_onebin", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight) | RecoMuon_isGEMMuon , 0, 1, 0, 1);

    }

    //pt//
    if(Eta_1p6_2p4){
      FillEfficiency("SAMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isSAMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("LooseMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isLooseMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TightMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isTightMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationAngLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationAngTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationAngLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationAngTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("GEMMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isGEMMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TrackerMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isTrackerMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("GEMMuon_or_TrackerMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isGEMMuon_or_isTrackerMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("GEMMuon_and_TrackerMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isGEMMuon_and_isTrackerMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("GEMSegment_eff_pt", GEMSegment_isMatched, igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("SAMuon_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isSAMuon) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("LooseMuon_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isLooseMuon) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TightMuon_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTightMuon) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationLoose_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationLoose) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationTight_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationTight) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationAngLoose_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMOneStationAngTight_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationLoose_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationLoose) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationTight_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationTight) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationAngLoose_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
      FillEfficiency("TMLastStationAngTight_or_GEMMuon_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight) || RecoMuon_isGEMMuon , igenP4.Pt(), n_pt_bin, pt_bin);
    }

    //eta//
    if(Pt_5){
      FillEfficiency("SAMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isSAMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("LooseMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isLooseMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TightMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isTightMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationAngLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationAngTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationAngLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationAngTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("GEMMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isGEMMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TrackerMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isTrackerMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("GEMMuon_or_TrackerMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isGEMMuon_or_isTrackerMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("GEMMuon_and_TrackerMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isGEMMuon_and_isTrackerMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("GEMSegment_eff_eta", GEMSegment_isMatched, igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("SAMuon_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isSAMuon) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("LooseMuon_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isLooseMuon) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TightMuon_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTightMuon) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationLoose_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationLoose) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationTight_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationTight) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationAngLoose_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMOneStationAngTight_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationLoose_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationLoose) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationTight_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationTight) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationAngLoose_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
      FillEfficiency("TMLastStationAngTight_or_GEMMuon_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight) || RecoMuon_isGEMMuon , igenP4.Eta(), n_eta_bin, eta_bin);
    }

  }
}

int TrackerGEMEfficiencyAnalyzer::FindWhichBin(double value, double *bin, int n_bin){
  // bin[0], bin[1], ... , bin[n_bin] : n_bin+1 elements !! //
  if(value < bin[0] || value >= bin[n_bin]){
    cout << "value = " << value << " : not in range [" << bin[0] << ", " << bin[n_bin] << ")" << endl;
    return n_bin+1;
  }
  else{
    for(int i=0; i<n_bin; i++){
      if( bin[i] <= value && value < bin[i+1] ) return i;
    }
    return n_bin+2;
  }
}

void TrackerGEMEfficiencyAnalyzer::MakeEfficiency(TString hname, int nbins, double xbins[]){
  mapefficiency[hname] =  new TEfficiency(hname.Data(),hname.Data(),nbins,xbins);
}

void TrackerGEMEfficiencyAnalyzer::MakeEfficiency(TString hname, int nbins, double xmin, double xmax){
  mapefficiency[hname] =  new TEfficiency(hname.Data(),hname.Data(),nbins,xmin,xmax);
}

TEfficiency* TrackerGEMEfficiencyAnalyzer::GetEfficiency(TString hname){
  TEfficiency* h = NULL;
  std::map<TString, TEfficiency*>::iterator mapit = mapefficiency.find(hname);
  if(mapit != mapefficiency.end()) return mapit->second;

  return h;
}

void TrackerGEMEfficiencyAnalyzer::FillEfficiency(TString histname, bool pass, double value, int nbins, double xbins[]){
  if(GetEfficiency(histname)) GetEfficiency(histname)->Fill(pass, value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeEfficiency(histname, nbins, xbins);
    if(GetEfficiency(histname)) GetEfficiency(histname)->Fill(pass, value);
  }
}

void TrackerGEMEfficiencyAnalyzer::FillEfficiency(TString histname, bool pass, double value, int nbins, double xmin, double xmax){
  if(GetEfficiency(histname)) GetEfficiency(histname)->Fill(pass, value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeEfficiency(histname, nbins, xmin, xmax);
    if(GetEfficiency(histname)) GetEfficiency(histname)->Fill(pass, value);
  }
}

void TrackerGEMEfficiencyAnalyzer::WriteEfficiencies(){
  for(map<TString, TEfficiency*>::iterator mapit = mapefficiency.begin(); mapit != mapefficiency.end(); mapit++){
    mapit->second->Write();
  }
}

void TrackerGEMEfficiencyAnalyzer::MakeHistograms(TString hname, int nbins, double xbins[]){
  maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xbins);
}

void TrackerGEMEfficiencyAnalyzer::MakeHistograms(TString hname, int nbins, double xmin, double xmax){
  maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xmin,xmax);
}

TH1* TrackerGEMEfficiencyAnalyzer::GetHist(TString hname){
  TH1* h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit->second;

  return h;
}

void TrackerGEMEfficiencyAnalyzer::FillHist(TString histname, double value, int nbins, double xbins[]){
  if(GetHist(histname)) GetHist(histname)->Fill(value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xbins);
    if(GetHist(histname)) GetHist(histname)->Fill(value);
  }
}

void TrackerGEMEfficiencyAnalyzer::FillHist(TString histname, double value, int nbins, double xmin, double xmax){
  if(GetHist(histname)) GetHist(histname)->Fill(value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xmin, xmax);
    if(GetHist(histname)) GetHist(histname)->Fill(value);
  }
}

void TrackerGEMEfficiencyAnalyzer::WriteHists(){
  for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit++){
    mapit->second->Write();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackerGEMEfficiencyAnalyzer);

