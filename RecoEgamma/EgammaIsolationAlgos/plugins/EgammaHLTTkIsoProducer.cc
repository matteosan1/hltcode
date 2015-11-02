#include "RecoEgamma/EgammaIsolationAlgos/plugins/EgammaHLTTkIsoProducer.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"

template<typename T1>
EgammaHLTTkIsoProducer<T1>::EgammaHLTTkIsoProducer(const edm::ParameterSet& config): 
  emObjectProducer_           (consumes<T1Collection>(config.getParameter<edm::InputTag>("candidateProducer"))),
  trackProducer_              (consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("trackProducer"))),
  beamSpotProducer_           (consumes<reco::BeamSpot>(config.getParameter<edm::InputTag>("beamSpotProducer"))),
  egTrkIsoPtMin_              (config.getParameter<double>("egTrkIsoPtMin")),
  egTrkIsoConeSize_           (config.getParameter<double>("egTrkIsoConeSize")),
  egTrkIsoZSpan_              (config.getParameter<double>("egTrkIsoZSpan")),
  egTrkIsoRSpan_              (config.getParameter<double>("egTrkIsoRSpan")),
  egTrkIsoVetoConeSizeBarrel_ (config.getParameter<double>("egTrkIsoVetoConeSizeBarrel")),
  egTrkIsoVetoConeSizeEndcap_ (config.getParameter<double>("egTrkIsoVetoConeSizeEndcap")),
  egTrkIsoStripBarrel_        (config.getParameter<double>("egTrkIsoStripBarrel")),
  egTrkIsoStripEndcap_        (config.getParameter<double>("egTrkIsoStripEndcap")) {

  produces <edm::ValueMap<float>>();
}

template<typename T1>
EgammaHLTTkIsoProducer<T1>::~EgammaHLTTkIsoProducer()
{}

template<typename T1>
void EgammaHLTTkIsoProducer<T1>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  //edm::ParameterSetDescription desc;
  //desc.add<edm::InputTag>("candidateProducer", edm::InputTag("gedGsfElectrons"));
  //desc.add<edm::InputTag>("pfClusterProducer", edm::InputTag("particleFlowClusterECAL")); 
  //desc.add<double>("drMax", 0.3);
  //desc.add<double>("drVetoBarrel", 0.0);
  //desc.add<double>("drVetoEndcap", 0.0);
  //desc.add<double>("etaStripBarrel", 0.0);
  //desc.add<double>("etaStripEndcap", 0.0);
  //desc.add<double>("energyBarrel", 0.0);
  //desc.add<double>("energyEndcap", 0.0);
  //descriptions.add(defaultModuleLabel<EgammaHLTTkIsoProducer<T1>>(), desc);
}

template<typename T1>
void EgammaHLTTkIsoProducer<T1>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<T1Collection> emObjectHandle;
  iEvent.getByToken(emObjectProducer_, emObjectHandle);

  std::auto_ptr<edm::ValueMap<float> > isoMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*isoMap);
  std::vector<float> retV(emObjectHandle->size(),0);

  edm::Handle<reco::TrackCollection> tkHandle;
  iEvent.getByToken(trackProducer_, tkHandle);
  const reco::TrackCollection* trackCollection = tkHandle.product();

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotProducer_,recoBeamSpotHandle);
  const reco::BeamSpot::Point& beamSpotPosition = recoBeamSpotHandle->position(); 
  
  ElectronTkIsolation isoAlgo(egTrkIsoConeSize_, egTrkIsoVetoConeSizeBarrel_, egTrkIsoVetoConeSizeEndcap_, egTrkIsoStripBarrel_, egTrkIsoStripEndcap_, 
			      egTrkIsoPtMin_, egTrkIsoZSpan_, egTrkIsoRSpan_, trackCollection, beamSpotPosition);
  
  for(unsigned int iReco = 0; iReco < emObjectHandle->size(); iReco++) {
    T1Ref candRef(emObjectHandle, iReco);
    //const reco::Track* eleTrk = &*candRef->gsfTrack();
    retV[iReco] = isoAlgo.getIso(&*candRef).second;
  }
  
  filler.insert(emObjectHandle,retV.begin(),retV.end());
  filler.fill();
  
  iEvent.put(isoMap);
}

typedef EgammaHLTTkIsoProducer<reco::GsfElectron> ElectronHLTTkIsoProducer;
typedef EgammaHLTTkIsoProducer<pat::Electron> PATElectronHLTTkIsoProducer;

DEFINE_FWK_MODULE(ElectronHLTTkIsoProducer);
DEFINE_FWK_MODULE(PATElectronHLTTkIsoProducer);
