#ifndef EgammaIsolationProducers_EgammaHLTTkIsoProducer_h
#define EgammaIsolationProducers_EgammaHLTTkIsoProducer_h

//*****************************************************************************
// File:      EgammaHLTTkIsoProducer.h
// ----------------------------------------------------------------------------
// OrigAuth:  Matteo Sani
// Institute: UCSD
//*****************************************************************************

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

template<typename T1>
class EgammaHLTTkIsoProducer : public edm::stream::EDProducer<> {
 public:

  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;

  explicit EgammaHLTTkIsoProducer(const edm::ParameterSet&);
  ~EgammaHLTTkIsoProducer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
  virtual void produce(edm::Event&, const edm::EventSetup&);
 private:

  const edm::EDGetTokenT<T1Collection> emObjectProducer_;
  const edm::EDGetTokenT<reco::TrackCollection> trackProducer_; 
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotProducer_;

  const double egTrkIsoPtMin_              ;
  const double egTrkIsoConeSize_           ;
  const double egTrkIsoZSpan_              ;
  const double egTrkIsoRSpan_              ;
  const double egTrkIsoVetoConeSizeBarrel_ ;
  const double egTrkIsoVetoConeSizeEndcap_ ;
  const double egTrkIsoStripBarrel_        ;
  const double egTrkIsoStripEndcap_        ;
};

#endif
