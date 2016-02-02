// -*- C++ -*-
//
// Package:    TrackerGEMFakeRateAnalyzer
// Class:      TrackerGEMFakeRateAnalyzer
// 
/**\class TrackerGEMFakeRateAnalyzer TrackerGEMFakeRateAnalyzer.cc MyAnalyzers/TrackerGEMFakeRateAnalyzer/src/TrackerGEMFakeRateAnalyzer.cc

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


// root include files
#include "TFile.h"
#include "TH1I.h"
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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

using namespace std;

const int n_pt_bin = 19, n_eta_bin = 9;
double pt_bin[n_pt_bin+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100},
       eta_bin[n_eta_bin+1] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

//
// class declaration
//

class TrackerGEMFakeRateAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrackerGEMFakeRateAnalyzer(const edm::ParameterSet&);
      ~TrackerGEMFakeRateAnalyzer();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual int FindWhichBin(double value, double *bin, int n_bin);
      virtual void MakeHistograms(TString hname, int nbins, double xbins[]);
      virtual void MakeHistograms(TString hname, int nbins, double xmin, double xmax);
      virtual TH1* GetHist(TString hname);
      virtual void FillHist(TString histname, double value, int nbins, double xbins[]);
      virtual void FillHist(TString histname, double value, int nbins, double xmin, double xmax);
      virtual void WriteHists();
      virtual TString append_number(TString string, int number);
      // ----------member data ---------------------------

  edm::EDGetTokenT<reco::GenParticleCollection> GENParticle_Token;
  edm::EDGetTokenT<reco::MuonCollection> TrackerGEM_Token;
  edm::EDGetTokenT<reco::MuonCollection> recoMuon_Token;
  edm::EDGetTokenT<reco::VertexCollection> vertex_Token;

  std::string rootFileName;
  std::unique_ptr<TFile> outputfile;

  map<TString, TH1*> maphist;

};

//
// constants, enums and typedefs
//
// constructors and destructor
//
TrackerGEMFakeRateAnalyzer::TrackerGEMFakeRateAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  TrackerGEM_Token = consumes<reco::MuonCollection>(edm::InputTag("trackerGEM"));
  GENParticle_Token = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
  recoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));
  vertex_Token = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));

  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
  outputfile.reset(TFile::Open(rootFileName.c_str(), "RECREATE"));

}


TrackerGEMFakeRateAnalyzer::~TrackerGEMFakeRateAnalyzer()
{

  outputfile->cd();

  WriteHists();

  /// variables and bin information ///

  TH1D *hist_onebin_bin = new TH1D("hist_onebin_bin", "", 1+1, 0, 1+1);
  TH1D *hist_pt_bin = new TH1D("hist_pt_bin", "", n_pt_bin+1, 0, n_pt_bin+1);
  TH1D *hist_eta_bin = new TH1D("hist_eta_bin", "", n_eta_bin+1, 0, n_eta_bin+1);
  for(int i=0; i<1+1;i++){
    hist_onebin_bin->Fill(i, i);
  }
  for(int i=0;i<n_pt_bin+1;i++){
    hist_pt_bin->Fill(i, pt_bin[i]);
  }
  for(int i=0;i<n_eta_bin+1;i++){
    hist_eta_bin->Fill(i, eta_bin[i]);
  }
  hist_onebin_bin->Write();
  hist_pt_bin->Write();
  hist_eta_bin->Write();  

  outputfile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackerGEMFakeRateAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::MuonCollection> trkMuons;
  iEvent.getByToken(TrackerGEM_Token, trkMuons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(GENParticle_Token, genParticles);
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(recoMuon_Token, recoMuons);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertex_Token, vertices);

  map<TString, int> N_fake_onebin, N_fake_pt, N_fake_eta;

	double matching_deltaR = 0.3;

  N_fake_onebin["TrackerGEM"] = 0;
  for(int i_bin=0;i_bin<n_pt_bin;i_bin++) N_fake_pt[append_number("TrackerGEM", i_bin)] = 0;
  for(int i_bin=0;i_bin<n_eta_bin;i_bin++) N_fake_eta[append_number("TrackerGEM", i_bin)] = 0;
  N_fake_onebin["TrackerGEM_double"] = 0;
  for(int i_bin=0;i_bin<n_pt_bin;i_bin++) N_fake_pt[append_number("TrackerGEM_double", i_bin)] = 0;
  for(int i_bin=0;i_bin<n_eta_bin;i_bin++) N_fake_eta[append_number("TrackerGEM_double", i_bin)] = 0;

  /// loop over trackerGEM ///
  for(reco::MuonCollection::const_iterator trkmuon=trkMuons->begin(); trkmuon != trkMuons->end(); ++trkmuon) {
    if( trkmuon->pt() < pt_bin[0]        || trkmuon->pt() >= pt_bin[n_pt_bin] ||
        abs(trkmuon->eta()) < eta_bin[0] || abs(trkmuon->eta()) >= eta_bin[n_eta_bin] ) continue; // region of interest
    bool isDoubleMatched = trkmuon->matches().size() == 2? true : false;
    TLorentzVector itrkP4;
    itrkP4.SetPtEtaPhiM(trkmuon->pt(), trkmuon->eta(), trkmuon->phi(), trkmuon->mass());
    bool isFake = true;
    /// loop over gen particles ///
    for(reco::GenParticleCollection::const_iterator genpart=genParticles->begin(); genpart != genParticles->end(); ++genpart) {
      TLorentzVector igenP4;
      igenP4.SetPtEtaPhiM(genpart->pt(), genpart->eta(), genpart->phi(), genpart->mass());
      if( TMath::Abs( genpart->pdgId() ) == 13 && itrkP4.DeltaR(igenP4) < matching_deltaR ){
        isFake = false;
        break;
      }
    } // genpart loop end
    if(isFake){
      int this_pt_bin = FindWhichBin(itrkP4.Pt(), pt_bin, n_pt_bin), this_eta_bin = FindWhichBin(abs( itrkP4.Eta() ), eta_bin, n_eta_bin);
      bool Eta_1p6_2p4 = abs(itrkP4.Eta()) > 1.6 && abs(itrkP4.Eta()) < 2.4,
                  Pt_5 = itrkP4.Pt() > 5;

      if(Eta_1p6_2p4 && Pt_5) N_fake_onebin["TrackerGEM"]++;
      if(Eta_1p6_2p4) N_fake_pt[append_number("TrackerGEM", this_pt_bin)]++;
      if(Pt_5) N_fake_eta[append_number("TrackerGEM", this_eta_bin)]++;
      if(isDoubleMatched){
        if(Eta_1p6_2p4 && Pt_5) N_fake_onebin["TrackerGEM_double"]++;
        if(Eta_1p6_2p4) N_fake_pt[append_number("TrackerGEM_double", this_pt_bin)]++;
        if(Pt_5) N_fake_eta[append_number("TrackerGEM_double", this_eta_bin)]++;
      }
    }
  } // trackGEM loop end
  FillHist("n_fake_TrackerGEM_onebin", N_fake_onebin["TrackerGEM"], 100, 0, 100);
  FillHist("n_fake_TrackerGEM_double_onebin", N_fake_onebin["TrackerGEM_double"], 100, 0, 100);
  for(int i=0;i<n_pt_bin; i++){
    FillHist(append_number("n_fake_TrackerGEM_pt", i), N_fake_pt[append_number("TrackerGEM", i)], 100, 0, 100);
    FillHist(append_number("n_fake_TrackerGEM_double_pt", i), N_fake_pt[append_number("TrackerGEM_double", i)], 100, 0, 100);
  }
  for(int i=0;i<n_eta_bin; i++){
    FillHist(append_number("n_fake_TrackerGEM_eta", i), N_fake_eta[append_number("TrackerGEM", i)], 100, 0, 100);
    FillHist(append_number("n_fake_TrackerGEM_double_eta", i), N_fake_eta[append_number("TrackerGEM_double", i)], 100, 0, 100);
  }

  const int n_obj = 11;
  TString obj[n_obj] = {"SAMuon", "LooseMuon", "TightMuon",
                        "TMOneStationLoose", "TMOneStationTight",
                        "TMOneStationAngLoose", "TMOneStationAngTight",
                        "TMLastStationLoose", "TMLastStationTight",
                        "TMLastStationAngLoose", "TMLastStationAngTight"};
  TString ortrk = "_or_TrackerGEM";
  // initializing N_fake's
  for(int i_obj=0;i_obj<n_obj;i_obj++){
    N_fake_onebin[obj[i_obj]] = 0;
    for(int i_bin=0;i_bin<n_pt_bin;i_bin++){
      N_fake_pt[append_number(obj[i_obj], i_bin)] = 0;
    }
    for(int i_bin=0;i_bin<n_eta_bin;i_bin++){
      N_fake_eta[append_number(obj[i_obj], i_bin)] = 0;
    }
  }

  /// loop over reco muon ///
  for(reco::MuonCollection::const_iterator recomuon = recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon){
    if( recomuon->pt() < pt_bin[0]        || recomuon->pt() >= pt_bin[n_pt_bin] ||
        abs(recomuon->eta()) < eta_bin[0] || abs(recomuon->eta()) >= eta_bin[n_eta_bin] ) continue; // region of interest
    TLorentzVector irecoP4;
    irecoP4.SetPtEtaPhiM(recomuon->pt(), recomuon->eta(), recomuon->phi(), recomuon->mass());
    map<TString, bool> recoMuon_isObj;
    recoMuon_isObj["SAMuon"] = false;
    recoMuon_isObj["LooseMuon"] = muon::isLooseMuon(*recomuon);
    recoMuon_isObj["TightMuon"] = muon::isTightMuon(*recomuon, *(vertices->begin()));
    recoMuon_isObj["TMOneStationLoose"] = muon::isGoodMuon(*recomuon, muon::TMOneStationLoose);
    recoMuon_isObj["TMOneStationTight"] = muon::isGoodMuon(*recomuon, muon::TMOneStationTight);
    recoMuon_isObj["TMOneStationAngLoose"] = muon::isGoodMuon(*recomuon, muon::TMOneStationAngLoose);
    recoMuon_isObj["TMOneStationAngTight"] = muon::isGoodMuon(*recomuon, muon::TMOneStationAngTight);
    recoMuon_isObj["TMLastStationLoose"] = muon::isGoodMuon(*recomuon, muon::TMLastStationLoose);
    recoMuon_isObj["TMLastStationTight"] = muon::isGoodMuon(*recomuon, muon::TMLastStationTight);
    recoMuon_isObj["TMLastStationAngLoose"] = muon::isGoodMuon(*recomuon, muon::TMLastStationAngLoose);
    recoMuon_isObj["TMLastStationAngTight"] = muon::isGoodMuon(*recomuon, muon::TMLastStationAngTight);
    if( recomuon->globalTrack().isNonnull() ){
      recoMuon_isObj["SAMuon"] = (recomuon->isStandAloneMuon())
                                 //& (recomuon->numberOfMatchedStations() >= 1) 
                                 & (recomuon->outerTrack()->hitPattern().muonStationsWithValidHits() > 1)
                                 & ( std::abs(recomuon->time().timeAtIpInOut) < ( 12.5 + std::abs( recomuon->time().timeAtIpInOutErr ) ) );
    }

    bool isFake = true;
    /// loop over gen particles ///
    for(reco::GenParticleCollection::const_iterator genpart=genParticles->begin(); genpart != genParticles->end(); ++genpart) {
      TLorentzVector igenP4;
      igenP4.SetPtEtaPhiM(genpart->pt(), genpart->eta(), genpart->phi(), genpart->mass());
      if( TMath::Abs( genpart->pdgId() ) == 13 && irecoP4.DeltaR(igenP4) < matching_deltaR ){
        isFake = false;
        break;
      }
    } // genpart loop end
    if(isFake){
      int this_pt_bin = FindWhichBin(irecoP4.Pt(), pt_bin, n_pt_bin), this_eta_bin = FindWhichBin(abs( irecoP4.Eta() ), eta_bin, n_eta_bin);
      bool Eta_1p6_2p4 = abs(irecoP4.Eta()) > 1.6 && abs(irecoP4.Eta()) < 2.4,
                  Pt_5 = irecoP4.Pt() > 5;
      for(int i_obj=0;i_obj<n_obj;i_obj++){
        if( recoMuon_isObj[obj[i_obj]] ){
          if(Eta_1p6_2p4 && Pt_5) N_fake_onebin[obj[i_obj]]++;
          if(Eta_1p6_2p4) N_fake_pt[append_number(obj[i_obj], this_pt_bin)]++;
          if(Pt_5) N_fake_eta[append_number(obj[i_obj], this_eta_bin)]++;
        }
      } 
    }
  }
  for(int i_obj=0;i_obj<n_obj;i_obj++){
    N_fake_onebin[obj[i_obj]+ortrk] = N_fake_onebin["TrackerGEM"]+N_fake_onebin[obj[i_obj]];
    FillHist("n_fake_"+obj[i_obj]+"_onebin", N_fake_onebin[obj[i_obj]], 100, 0, 100);
    FillHist("n_fake_"+obj[i_obj]+ortrk+"_onebin", N_fake_onebin[obj[i_obj]+ortrk], 100, 0, 100);
    for(int i_bin=0;i_bin<n_pt_bin;i_bin++){
      N_fake_pt[append_number(obj[i_obj]+ortrk, i_bin)] = N_fake_pt[append_number("TrackerGEM", i_bin)]+N_fake_pt[append_number(obj[i_obj], i_bin)];
      FillHist(append_number("n_fake_"+obj[i_obj]+"_pt", i_bin), N_fake_pt[append_number(obj[i_obj], i_bin)], 100, 0, 100);
      FillHist(append_number("n_fake_"+obj[i_obj]+ortrk+"_pt", i_bin), N_fake_pt[append_number(obj[i_obj]+ortrk, i_bin)], 100, 0, 100);
    }
    for(int i_bin=0;i_bin<n_eta_bin;i_bin++){
      N_fake_eta[append_number(obj[i_obj]+ortrk, i_bin)] = N_fake_eta[append_number("TrackerGEM", i_bin)]+N_fake_eta[append_number(obj[i_obj], i_bin)];
      FillHist(append_number("n_fake_"+obj[i_obj]+"_eta", i_bin), N_fake_eta[append_number(obj[i_obj], i_bin)], 100, 0, 100);
      FillHist(append_number("n_fake_"+obj[i_obj]+ortrk+"_eta", i_bin), N_fake_eta[append_number(obj[i_obj]+ortrk, i_bin)], 100, 0, 100);
    }
  }
  

}

int TrackerGEMFakeRateAnalyzer::FindWhichBin(double value, double *bin, int n_bin){
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

void TrackerGEMFakeRateAnalyzer::MakeHistograms(TString hname, int nbins, double xbins[]){
  maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xbins);
}

void TrackerGEMFakeRateAnalyzer::MakeHistograms(TString hname, int nbins, double xmin, double xmax){
  maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xmin,xmax);
}

TH1* TrackerGEMFakeRateAnalyzer::GetHist(TString hname){
  TH1* h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit->second;

  return h;
}

void TrackerGEMFakeRateAnalyzer::FillHist(TString histname, double value, int nbins, double xbins[]){
  if(GetHist(histname)) GetHist(histname)->Fill(value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xbins);
    if(GetHist(histname)) GetHist(histname)->Fill(value);
  }
}

void TrackerGEMFakeRateAnalyzer::FillHist(TString histname, double value, int nbins, double xmin, double xmax){
  if(GetHist(histname)) GetHist(histname)->Fill(value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xmin, xmax);
    if(GetHist(histname)) GetHist(histname)->Fill(value);
  }
}

void TrackerGEMFakeRateAnalyzer::WriteHists(){
  for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit++){
    mapit->second->Write();
  }
}

TString TrackerGEMFakeRateAnalyzer::append_number(TString string, int number){
  return string+"_"+TString::Itoa(number, 10);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackerGEMFakeRateAnalyzer);
