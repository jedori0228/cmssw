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
  //Associator for chi2: objects
  //edm::InputTag associatormap;
  //

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
  edm::EDGetTokenT<reco::TrackCollection > generalTracksToken_;

  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > track_Collection_Token;

  bool UseAssociators;

  std::vector<std::string> associators;
  std::vector<const TrackAssociatorBase*> associator;
  GenParticleCustomSelector gpSelector;
  //std::string parametersDefiner;
  
  std::vector<edm::InputTag> label;

  virtual void MakeHistograms(TString hname, int nbins, double xbins[]);
  virtual void MakeHistograms(TString hname, int nbins, double xmin, double xmax);
  virtual TH1* GetHist(TString hname);
  virtual void FillHist(TString histname, double value, int nbins, double xbins[]);
  virtual void FillHist(TString histname, double value, int nbins, double xmin, double xmax);
  virtual void WriteHists();

  std::string rootFileName;
  std::unique_ptr<TFile> outputfile;
  map<TString, TH1*> maphist;

  double  FakeRatePtCut, MatchingWindowDelR;

};

GEMMuonAnalyzer::GEMMuonAnalyzer(const edm::ParameterSet& iConfig)
{
  outputfile.reset(TFile::Open(iConfig.getParameter<std::string>("HistoFile").c_str(), "RECREATE"));
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");

  FakeRatePtCut   = iConfig.getParameter<double>("FakeRatePtCut");
  MatchingWindowDelR   = iConfig.getParameter<double>("MatchingWindowDelR");

  //Associator for chi2: getting parametters
  //associatormap = iConfig.getParameter< edm::InputTag >("associatormap");
  UseAssociators = iConfig.getParameter< bool >("UseAssociators");
  associators = iConfig.getParameter< std::vector<std::string> >("associators");

  gpSelector = GenParticleCustomSelector(iConfig.getParameter<double>("ptMinGP"),
           iConfig.getParameter<double>("minRapidityGP"),
           iConfig.getParameter<double>("maxRapidityGP"),
           iConfig.getParameter<double>("tipGP"),
           iConfig.getParameter<double>("lipGP"),
           iConfig.getParameter<bool>("chargedOnlyGP"),
           iConfig.getParameter<int>("statusGP"),
           iConfig.getParameter<std::vector<int> >("pdgIdGP"));
  //parametersDefiner =iConfig.getParameter<std::string>("parametersDefiner");

  label = iConfig.getParameter< std::vector<edm::InputTag> >("label");

  genParticlesToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
  edm::InputTag trackingParticlesTag ("mix", "MergedTrackTruth");
  trackingParticlesToken_ = consumes<TrackingParticleCollection>(trackingParticlesTag);
  generalTracksToken_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));

  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));

  for (unsigned int www=0;www<label.size();www++){
    track_Collection_Token.push_back(consumes<edm::View<reco::Track> >(label[www]));
  }

}

void GEMMuonAnalyzer::beginJob()
{


}

GEMMuonAnalyzer::~GEMMuonAnalyzer(){


  outputfile->cd();

  WriteHists();

  outputfile->Close();

}

void
GEMMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{
  using namespace edm;
  using namespace reco;

  if (UseAssociators) {
    edm::ESHandle<TrackAssociatorBase> theAssociator;
    for (unsigned int w=0;w<associators.size();w++) {
      iSetup.get<TrackAssociatorRecord>().get(associators[w],theAssociator);
      associator.push_back( theAssociator.product() );
    }
  }

  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);
  Handle<TrackingParticleCollection> trackingParticles;
  iEvent.getByToken(trackingParticlesToken_, trackingParticles);
  Handle <reco::TrackCollection> generalTracks;
  iEvent.getByToken (generalTracksToken_, generalTracks);
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);


  //vector<reco::Track> GEMMuonBestTracks;
  //for(reco::MuonCollection::const_iterator recomuon=recoMuons->begin(); recomuon != recoMuons->end(); ++recomuon) {
  //  if(recomuon->isGEMMuon()) GEMMuonBestTracks->push_back(recoMuons->);
  //}

  ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);
  ESHandle<Propagator> shProp;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", shProp);

  //Track Association by hits:

  if (UseAssociators) {
    for (unsigned int ww=0;ww<associators.size();ww++){
for (unsigned int www=0;www<label.size();www++){     

      reco::RecoToSimCollection recSimColl;
	    reco::SimToRecoCollection simRecColl;
	    edm::Handle<View<Track> >  trackCollection;

      unsigned int trackCollectionSize = 0;

	    if( !iEvent.getByToken(track_Collection_Token[www], trackCollection) ){ //FIXME we should replace this to GEMMuons
	      recSimColl.post_insert();
	      simRecColl.post_insert();
        std::cout << "failed to get trackCollection" << std::endl;
	  
	    }
	    else {
	      trackCollectionSize = trackCollection->size();
        std::cout << "trackCollectionSize = " << trackCollectionSize << std::endl;

	      recSimColl = associator[ww]->associateRecoToSim(trackCollection, trackingParticles, &iEvent, &iSetup);
	      simRecColl = associator[ww]->associateSimToReco(trackCollection, trackingParticles, &iEvent, &iSetup);

	    }



      // denominators for efficiencies
      for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++){

	      TrackingParticleRef tpr(trackingParticles, i);
	      TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get()); 
	      TrackingParticle::Vector momentumTP; 
	      TrackingParticle::Point vertexTP;
	  
	      if (abs(tp->pdgId()) != 13) continue;

        bool Eta_1p6_2p4 = abs(tp->eta()) > 1.6 && abs(tp->eta()) < 2.4,
             Pt_5 = tp->pt() > 5;

        if( Eta_1p6_2p4 && Pt_5 ){
          bool SignalMuon = false;

          if(tp->status() != -99){ // Pythia8 gen status : home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
	          int motherid=-1;
	          if ((*tp->genParticle_begin())->numberOfMothers()>0)  {
		          if ((*tp->genParticle_begin())->mother()->numberOfMothers()>0){
		            motherid=(*tp->genParticle_begin())->mother()->mother()->pdgId();
		          }
            }  
            std::cout<<"Mother ID = "<<motherid<<std::endl;

  	        if ( ( (tp->status()==1) && ( (*tp->genParticle_begin())->numberOfMothers()==0 ) )  ||
  		           ( (tp->status()==1) )      )    SignalMuon=true;
	    
          } // END if(tp->status() != -99)	      

          if(SignalMuon){
            FillHist("TPMuon_Eta", fabs(tp->eta()), 9, 1.5, 2.4);
          }

        } // END if( Eta_1p6_2p4 && Pt_5 )


      } // END for (TrackingParticleCollection::size_type i=0; i<trackingParticles->size(); i++)


      // loop over our tracks
      for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
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
  	        tpr = tp.begin()->first;

  	        //double assocChi2 = -(tp.begin()->second);
 
            //So this track is matched to a gen particle, lets get that gen particle now

            if ( (simRecColl.find(tpr) != simRecColl.end()) ){
              std::vector<std::pair<RefToBase<Track>, double> > rt;
              if(simRecColl[tpr].size() > 0){
                rt=simRecColl[tpr];
  		          RefToBase<Track> bestrecotrackforeff = rt.begin()->first;
                //Only fill the efficiency histo if the track found matches up to a gen particle's best choice
                if ( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) ) {
                  TrackIsEfficient=true;
                  //This section fills the numerator of the efficiency calculation...
                  
                  bool Eta_1p6_2p4 = abs(tpr->eta()) > 1.6 && abs(tpr->eta()) < 2.4,
                       Pt_5 = tpr->pt() > 5;
                  if( Eta_1p6_2p4 && Pt_5 ){
                   
                    bool SignalMuon=false;
                    
  		              if (tpr->status() !=-99){
  			              int motherid=-1;
  			              if ((*tpr->genParticle_begin())->numberOfMothers()>0)  {
  			                if ((*tpr->genParticle_begin())->mother()->numberOfMothers()>0){
  			                  motherid=(*tpr->genParticle_begin())->mother()->mother()->pdgId();
  			                }
  			              }
                      std::cout<<"Mother ID = "<<motherid<<std::endl;
  			              if ( 
  			                  ( (tpr->status()==1) && ( (*tpr->genParticle_begin())->numberOfMothers()==0 ) )  ||
  			                  ( (tpr->status()==1)  ) ) SignalMuon=true;

  		              } // END if (tpr->status() !=-99)
                    if(SignalMuon){
                      FillHist("Chi2MatchedME0Muon_Eta", fabs(tpr->eta()), 9, 1.5, 2.4 );

                    } // END if(SignalMuon)

                  } // END if( Eta_1p6_2p4 && Pt_5 )            

                } // END if ( (bestrecotrackforeff == track ) && (abs(tpr->pdgId()) == 13) )
              } // END if(simRecColl[tpr].size() > 0) 
            } // END  if ( (simRecColl.find(tpr) != simRecColl.end()) )

          } // END if (tp.size()!=0)
        } // END if(recSimColl.find(track) != recSimColl.end())

        // A simple way of measuring fake rate
        if (!TrackIsEfficient) {

          bool Eta_1p6_2p4 = abs(tpr->eta()) > 1.6 && abs(tpr->eta()) < 2.4,
               Pt_5 = tpr->pt() > 5;
          if( Eta_1p6_2p4 && Pt_5 ){
            FillHist("Chi2UnmatchedME0Muon_Eta", fabs(track->eta()), 9, 1.5,2.4 );
          } 

        } // END if (!TrackIsEfficient)

      } // END for(View<Track>::size_type i=0; i<trackCollectionSize; ++i)


}// END for (unsigned int www=0;www<label.size();www++)

    } // END for (unsigned int ww=0;ww<associators.size();ww++)
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

void GEMMuonAnalyzer::MakeHistograms(TString hname, int nbins, double xbins[]){
  maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xbins);
}

void GEMMuonAnalyzer::MakeHistograms(TString hname, int nbins, double xmin, double xmax){
  maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xmin,xmax);
}

TH1* GEMMuonAnalyzer::GetHist(TString hname){
  TH1* h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit->second;

  return h;
}

void GEMMuonAnalyzer::FillHist(TString histname, double value, int nbins, double xbins[]){
  if(GetHist(histname)) GetHist(histname)->Fill(value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xbins);
    if(GetHist(histname)) GetHist(histname)->Fill(value);
  }
}

void GEMMuonAnalyzer::FillHist(TString histname, double value, int nbins, double xmin, double xmax){
  if(GetHist(histname)) GetHist(histname)->Fill(value);
  else{
    if (nbins < 0) {
      exit(0);
    }
    MakeHistograms(histname, nbins, xmin, xmax);
    if(GetHist(histname)) GetHist(histname)->Fill(value);
  }
}

void GEMMuonAnalyzer::WriteHists(){
  for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit++){
    mapit->second->Write();
  }
}

DEFINE_FWK_MODULE(GEMMuonAnalyzer);
