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
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/GEMRecHit/interface/GEMSegmentCollection.h>
#include <DataFormats/MuonReco/interface/MuonSelectors.h>

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
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void MakeHistograms(TString hname, int nbins, double xbins[]);
      virtual void MakeHistograms(TString hname, int nbins, double xmin, double xmax);
      virtual TEfficiency* GetHist(TString hname);
      virtual void FillHist(TString histname, bool pass, double value, int nbins, double xbins[]);
      virtual void FillHist(TString histname, bool pass, double value, int nbins, double xmin, double xmax);
      virtual void WriteHists();

      // ----------member data ---------------------------

  edm::ESHandle<GEMGeometry> gemGeom;

  edm::EDGetTokenT<reco::GenParticleCollection> GENParticle_Token;
  edm::EDGetTokenT<reco::MuonCollection> TrackerGEM_Token;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  edm::EDGetTokenT<GEMSegmentCollection> GEMSegment_Token;
  edm::EDGetTokenT<reco::VertexCollection> vertex_Token;

  std::string rootFileName;
  std::unique_ptr<TFile> outputfile;

  map<TString, TEfficiency*> maphist;

};

//
// constants, enums and typedefs
//
// constructors and destructor
//
TrackerGEMEfficiencyAnalyzer::TrackerGEMEfficiencyAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  TrackerGEM_Token = consumes<reco::MuonCollection>(edm::InputTag("trackerGEM"));
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  GENParticle_Token = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
  GEMSegment_Token = consumes<GEMSegmentCollection>(edm::InputTag("gemSegments"));
  vertex_Token = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));

  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
  outputfile.reset(TFile::Open(rootFileName.c_str(), "RECREATE"));

}


TrackerGEMEfficiencyAnalyzer::~TrackerGEMEfficiencyAnalyzer()
{

  outputfile->cd();

  WriteHists();

  outputfile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackerGEMEfficiencyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(gemGeom);
  edm::Handle<reco::MuonCollection> trkMuons;
  iEvent.getByToken(TrackerGEM_Token, trkMuons);
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(GENParticle_Token, genParticles);
  edm::Handle<GEMSegmentCollection> gemSegmentCollection;
  iEvent.getByToken(GEMSegment_Token, gemSegmentCollection);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertex_Token, vertices);

  /////////////////////////////////////////////////////
  /// USE HAVE TO USE MUON GUN SAMPLE FOR THIS CODE ///
  /////////////////////////////////////////////////////

  double matching_deltaR = 0.1;

  /// loop over gen particles = muons ///
  for(reco::GenParticleCollection::const_iterator genpart=genParticles->begin(); genpart != genParticles->end(); ++genpart) {
    if(fabs(genpart->pdgId()) != 13) continue;

    TLorentzVector igenP4;
    igenP4.SetPtEtaPhiM(genpart->pt(), genpart->eta(), genpart->phi(), genpart->mass());

    /// loop over trackerGEM ///
    bool TrackerGEM_isMatched = false, TrackerGEM_isDoubleMatched = false;
    double deltaR_trk_temp = matching_deltaR;
    for(reco::MuonCollection::const_iterator trkmuon=trkMuons->begin(); trkmuon != trkMuons->end(); ++trkmuon) {
      TLorentzVector itrkP4;
      itrkP4.SetPtEtaPhiM(trkmuon->pt(), trkmuon->eta(), trkmuon->phi(), trkmuon->mass());
      if( igenP4.DeltaR(itrkP4) < deltaR_trk_temp ){
        TrackerGEM_isMatched = true;
        TrackerGEM_isDoubleMatched = trkmuon->matches().size() == 2 ? true : false;
        deltaR_trk_temp = igenP4.DeltaR(itrkP4);
      }
    } // end of trackerGEM loop

    /// loop over RecoMuon ///
    bool RecoMuon_isMatched = false, RecoMuon_isSAMuon = false, RecoMuon_isLooseMuon = false, RecoMuon_isTightMuon = false, 
         RecoMuon_isTMOneStationLoose = false, RecoMuon_isTMOneStationTight = false,
         RecoMuon_isTMOneStationAngLoose = false, RecoMuon_isTMOneStationAngTight = false,
         RecoMuon_isTMLastStationLoose = false, RecoMuon_isTMLastStationTight = false,
         RecoMuon_isTMLastStationAngLoose = false, RecoMuon_isTMLastStationAngTight = false; 
    double deltaR_reco_temp = matching_deltaR;
    for(reco::MuonCollection::const_iterator recomuon=recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon) { 
      TLorentzVector irecoP4;
      irecoP4.SetPtEtaPhiM(recomuon->pt(), recomuon->eta(), recomuon->phi(), recomuon->mass());    
      if( igenP4.DeltaR(irecoP4) < deltaR_reco_temp ){
        RecoMuon_isMatched = true;
        if( recomuon->globalTrack().isNonnull() ){
          RecoMuon_isSAMuon = (recomuon->isStandAloneMuon())
                              //& (recomuon->numberOfMatchedStations() >= 1) 
                              & (recomuon->outerTrack()->hitPattern().muonStationsWithValidHits() > 1)
                              & ( std::abs(recomuon->time().timeAtIpInOut) < ( 12.5 + std::abs( recomuon->time().timeAtIpInOutErr ) ) );
        }
        RecoMuon_isLooseMuon = muon::isLooseMuon(*recomuon);
        RecoMuon_isTightMuon = muon::isTightMuon(*recomuon, *(vertices->begin()));
        RecoMuon_isTMOneStationLoose = muon::isGoodMuon(*recomuon, muon::TMOneStationLoose);
        RecoMuon_isTMOneStationTight = muon::isGoodMuon(*recomuon, muon::TMOneStationTight);
        RecoMuon_isTMOneStationAngLoose = muon::isGoodMuon(*recomuon, muon::TMOneStationAngLoose);
        RecoMuon_isTMOneStationAngTight = muon::isGoodMuon(*recomuon, muon::TMOneStationAngTight);
        RecoMuon_isTMLastStationLoose = muon::isGoodMuon(*recomuon, muon::TMLastStationLoose);
        RecoMuon_isTMLastStationTight = muon::isGoodMuon(*recomuon, muon::TMLastStationTight);
        RecoMuon_isTMLastStationAngLoose = muon::isGoodMuon(*recomuon, muon::TMLastStationAngLoose);
        RecoMuon_isTMLastStationAngTight = muon::isGoodMuon(*recomuon, muon::TMLastStationAngTight);
        deltaR_reco_temp = igenP4.DeltaR(irecoP4);
      }
    }
   
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

    bool Eta_1p6_2p4 = abs(igenP4.Eta()) > 1.6 && abs(igenP4.Eta()) < 2.4,
                Pt_5 = igenP4.Pt() > 5;
   
    //onebin//
    if(Eta_1p6_2p4 && Pt_5){
      FillHist("TrackerGEM_eff_onebin", TrackerGEM_isMatched, 0, 1, 0, 1);
      FillHist("TrackerGEM_double_eff_onebin", TrackerGEM_isDoubleMatched, 0, 1, 0, 1);
      FillHist("SAMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isSAMuon, 0, 1, 0, 1);
      FillHist("LooseMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isLooseMuon, 0, 1, 0, 1);
      FillHist("TightMuon_eff_onebin", RecoMuon_isMatched && RecoMuon_isTightMuon, 0, 1, 0, 1);
      FillHist("TMOneStationLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationLoose, 0, 1, 0, 1);
      FillHist("TMOneStationTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationTight, 0, 1, 0, 1);
      FillHist("TMOneStationAngLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose, 0, 1, 0, 1);
      FillHist("TMOneStationAngTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight, 0, 1, 0, 1);
      FillHist("TMLastStationLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationLoose, 0, 1, 0, 1);
      FillHist("TMLastStationTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationTight, 0, 1, 0, 1);
      FillHist("TMLastStationAngLoose_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose, 0, 1, 0, 1);
      FillHist("TMLastStationAngTight_eff_onebin", RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight, 0, 1, 0, 1);
      FillHist("GEMSegment_eff_onebin", GEMSegment_isMatched, 0, 1, 0, 1);
      FillHist("SAMuon_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isSAMuon) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("LooseMuon_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isLooseMuon) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TightMuon_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTightMuon) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMOneStationLoose_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationLoose) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMOneStationTight_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationTight) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMOneStationAngLoose_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMOneStationAngTight_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMLastStationLoose_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationLoose) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMLastStationTight_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationTight) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMLastStationAngLoose_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose) | TrackerGEM_isMatched , 0, 1, 0, 1);
      FillHist("TMLastStationAngTight_or_TrackerGEM_onebin_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight) | TrackerGEM_isMatched , 0, 1, 0, 1);

    }

    //pt//
    if(Eta_1p6_2p4){
      FillHist("TrackerGEM_eff_pt", TrackerGEM_isMatched, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TrackerGEM_double_eff_pt", TrackerGEM_isDoubleMatched, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("SAMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isSAMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("LooseMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isLooseMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TightMuon_eff_pt", RecoMuon_isMatched && RecoMuon_isTightMuon, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationAngLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationAngTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationAngLoose_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationAngTight_eff_pt", RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("GEMSegment_eff_pt", GEMSegment_isMatched, igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("SAMuon_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isSAMuon) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("LooseMuon_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isLooseMuon) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TightMuon_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTightMuon) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationLoose_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationLoose) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationTight_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationTight) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationAngLoose_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMOneStationAngTight_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationLoose_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationLoose) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationTight_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationTight) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationAngLoose_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
      FillHist("TMLastStationAngTight_or_TrackerGEM_eff_pt", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight) || TrackerGEM_isMatched , igenP4.Pt(), n_pt_bin, pt_bin);
    }

    //eta//
    if(Pt_5){
      FillHist("TrackerGEM_eff_eta", TrackerGEM_isMatched, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TrackerGEM_double_eff_eta", TrackerGEM_isDoubleMatched, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("SAMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isSAMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("LooseMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isLooseMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TightMuon_eff_eta", RecoMuon_isMatched && RecoMuon_isTightMuon, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationAngLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationAngTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationAngLoose_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationAngTight_eff_eta", RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("GEMSegment_eff_eta", GEMSegment_isMatched, igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("SAMuon_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isSAMuon) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("LooseMuon_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isLooseMuon) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TightMuon_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTightMuon) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationLoose_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationLoose) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationTight_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationTight) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationAngLoose_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngLoose) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMOneStationAngTight_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMOneStationAngTight) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationLoose_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationLoose) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationTight_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationTight) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationAngLoose_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngLoose) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
      FillHist("TMLastStationAngTight_or_TrackerGEM_eff_eta", (RecoMuon_isMatched && RecoMuon_isTMLastStationAngTight) || TrackerGEM_isMatched , igenP4.Eta(), n_eta_bin, eta_bin);
    }

  }

}

void TrackerGEMEfficiencyAnalyzer::MakeHistograms(TString hname, int nbins, double xbins[]){
  maphist[hname] =  new TEfficiency(hname.Data(),hname.Data(),nbins,xbins);
}

void TrackerGEMEfficiencyAnalyzer::MakeHistograms(TString hname, int nbins, double xmin, double xmax){
  maphist[hname] =  new TEfficiency(hname.Data(),hname.Data(),nbins,xmin,xmax);
}

TEfficiency* TrackerGEMEfficiencyAnalyzer::GetHist(TString hname){
  TEfficiency* h = NULL;
  std::map<TString, TEfficiency*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit->second;

  return h;
}

void TrackerGEMEfficiencyAnalyzer::FillHist(TString histname, bool pass, double value, int nbins, double xbins[]){
  if(GetHist(histname)) GetHist(histname)->Fill(pass, value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xbins);
    if(GetHist(histname)) GetHist(histname)->Fill(pass, value);
  }
}

void TrackerGEMEfficiencyAnalyzer::FillHist(TString histname, bool pass, double value, int nbins, double xmin, double xmax){
  if(GetHist(histname)) GetHist(histname)->Fill(pass, value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xmin, xmax);
    if(GetHist(histname)) GetHist(histname)->Fill(pass, value);
  }
}

void TrackerGEMEfficiencyAnalyzer::WriteHists(){
  for(map<TString, TEfficiency*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit++){
    mapit->second->Write();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackerGEMEfficiencyAnalyzer);

