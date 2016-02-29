#include "FWCore/Framework/interface/Event.h"

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

#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <DataFormats/GEMRecHit/interface/GEMSegmentCollection.h>

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

//Associator for chi2: Including header files
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
  bool UseAssociators;
  const TrackAssociatorByChi2* associatorByChi2;
  std::vector<std::string> associators;
  std::vector<const TrackAssociatorBase*> associator;
  edm::EDGetTokenT<reco::MuonCollection> RecoMuon_Token;
  GenParticleCustomSelector gpSelector;
  //std::string parametersDefiner;

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

  double Nevents;

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
  RecoMuon_Token = consumes<reco::MuonCollection>(edm::InputTag("muons"));

}

void GEMMuonAnalyzer::beginJob()
{

  Nevents=0;


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

  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(RecoMuon_Token, recoMuons);
  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel<GenParticleCollection>("genParticles", genParticles);
  const GenParticleCollection genParticlesForChi2 = *(genParticles.product());

  Nevents++;

  Handle <TrackCollection > generalTracks;
  iEvent.getByLabel <TrackCollection> ("generalTracks", generalTracks);

  Handle<GEMSegmentCollection> OurSegments;
  iEvent.getByLabel("gemSegments","",OurSegments);

  ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);
  ESHandle<Propagator> shProp;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", shProp);

  //For Track Association:

  //std::cout<<"ON BEGIN JOB:"<<std::endl;
  if (UseAssociators) {
    //std::cout<<"Inside if:"<<std::endl;
    edm::ESHandle<TrackAssociatorBase> theAssociator;
    //std::cout<<"associators size = "<<associators.size()<<std::endl;
    for (unsigned int w=0;w<associators.size();w++) {
      //std::cout<<"On step "<<w<<std::endl;
      iSetup.get<TrackAssociatorRecord>().get(associators[w],theAssociator);
      //std::cout<<"On step "<<w<<std::endl;
      associator.push_back( theAssociator.product() );
      //std::cout<<"On step "<<w<<std::endl;
    }
    //std::cout<<"Got this many associators: "<<associator.size()<<std::endl;
  }

  for (unsigned int ww=0;ww<associators.size();ww++){
  
    associatorByChi2 = dynamic_cast<const TrackAssociatorByChi2*>(associator[ww]);

    if (associatorByChi2==0) continue;

    reco::RecoToGenCollection recSimColl;
    reco::GenToRecoCollection simRecColl;

    edm::Handle<View<Track> > trackCollection;
    iEvent.getByLabel( , trackCollection) // FIXME

	  recSimColl=associatorByChi2->associateRecoToGen(trackCollection,
		  		                                          genParticles,
			 		                                          &iEvent,
					                                          &iSetup);
	  simRecColl=associatorByChi2->associateGenToReco(trackCollection,
					                                          genParticles,
					                                          &iEvent,
					                                          &iSetup);
    for (GenParticleCollection::size_type i=0; i<genParticlesForChi2.size(); i++){

      double quality = 0.;

	    GenParticleRef tpr(genParticles, i);
	    GenParticle* tp=const_cast<GenParticle*>(tpr.get());
	    TrackingParticle::Vector momentumTP; 
	    TrackingParticle::Point vertexTP;

	    //Collision like particle
	    if(! gpSelector(*tp)) continue;
	    momentumTP = tp->momentum();
	    vertexTP = tp->vertex();

      std::vector<std::pair<RefToBase<Track>, double> > rt;

	    //Check if the gen particle has been associated to any reco track
	    if(simRecColl.find(tpr) != simRecColl.end()){
        rt = (std::vector<std::pair<RefToBase<Track>, double> >) simRecColl[tpr];
        //It has, so we check that the pair TrackRef/double pair collection (vector of pairs) is not empty
        if (rt.size()!=0) {
	        //It is not empty, so there is at least one real track that the gen particle is matched to
	    
	        //We take the first element of the vector, .begin(), and the trackRef from it, ->first, this is our best possible track match
	        RefToBase<Track> assoc_recoTrack = rt.begin()->first;
	        std::cout<<"-----------------------------associated Track #"<<assoc_recoTrack.key()<<std::endl;

	        quality = rt.begin()->second;
	        std::cout
          << "TrackingParticle #" <<tpr.key() << " with pt=" << sqrt(momentumTP.perp2()) << " associated with quality:" << quality <<std::endl;

    	    //Also, seeing as we have found a gen particle that is matched to a track, it is efficient, so we put it in the numerator of the efficiency plot
	        //if (( sqrt(momentumTP.perp2()) > FakeRatePtCut) && (TMath::Abs(tp->eta()) < 2.8) )Chi2MatchedME0Muon_Eta->Fill(tp->eta());
	        //if ( ( assoc_recoTrack->pt() > FakeRatePtCut) && (TMath::Abs(tp->eta()) < 2.8) )Chi2MatchedME0Muon_Eta->Fill(tp->eta());
        }
	    }//END if(simRecColl.find(tpr) != simRecColl.end())
    }//END for (GenParticleCollection::size_type i=0; i<genParticlesForChi2.size(); i++)

    for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
	    //bool Track_is_matched = false; 
	    RefToBase<Track> track(trackCollection, i);

      bool Eta_1p6_2p4 = abs(track->eta()) > 1.6 && abs(track->eta()) < 2.4,
           Pt_5 = track->pt() > 5; 

	    //std::vector<std::pair<TrackingParticleRef, double> > tp;
	    std::vector<std::pair<GenParticleRef, double> > tp;
	    std::vector<std::pair<GenParticleRef, double> > tpforfake;
	    //TrackingParticleRef tpr;
	    GenParticleRef tpr;
	    GenParticleRef tprforfake;

	    //Check if the track is associated to any gen particle
	    if(recSimColl.find(track) != recSimColl.end()){
	  
	      tp = recSimColl[track];
	      if (tp.size()!=0) {
	        //Track_is_matched = true;
	        tpr = tp.begin()->first;

	        double assocChi2 = -(tp.begin()->second);
	   
	        //So this track is matched to a gen particle, lets get that gen particle now
	        if (  (simRecColl.find(tpr) != simRecColl.end())    ){
	          std::vector<std::pair<RefToBase<Track>, double> > rt;
	          std::cout<<"Comparing gen and reco tracks"<<std::endl;
	          if  (simRecColl[tpr].size() > 0){
		          rt=simRecColl[tpr];
		          RefToBase<Track> bestrecotrackforeff = rt.begin()->first;
		          //Only fill the efficiency histo if the track found matches up to a gen particle's best choice
		          if (bestrecotrackforeff == track) {
		            if ( Eta_1p6_2p4 && Pt_5 )FillHist("Chi2MatchedME0Muon_Eta", fabs(tpr->eta()), 9, 1.5, 2.4)  ;
		            if ( Eta_1p6_2p4 && Pt_5 )FillHist("AssociatedChi2_h", assocChi2,50,0,50);
		            if ( Eta_1p6_2p4 && Pt_5 )FillHist("AssociatedChi2_Prob_h", TMath::Prob((assocChi2)*5,5), 50, 0, 1);
		            std::cout<<"assocChi2 = "<<assocChi2<<std::endl;
		          }
	          }
	        }
	      }
	    }
    }
      



 

  } // associators loop

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
