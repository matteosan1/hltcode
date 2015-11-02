//*****************************************************************************
// File:      EgammaRecHitIsolation.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer, hacked by Sam Harper (ie the ugly stuff is mine)
// Institute: IIHE-VUB, RAL
//=============================================================================
//*****************************************************************************

#include "RecoEgamma/EgammaIsolationAlgos/interface/EcalPFClusterIsolation.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include <DataFormats/Math/interface/deltaR.h>

template<typename T1>
EcalPFClusterIsolation<T1>::EcalPFClusterIsolation(double drMax,
						   double drVetoBarrel,
						   double drVetoEndcap,
						   double etaStripBarrel,
						   double etaStripEndcap,
						   double energyBarrel,
						   double energyEndcap):
  drMax_(drMax),
  drVetoBarrel_(drVetoBarrel),
  drVetoEndcap_(drVetoEndcap),
  etaStripBarrel_(etaStripBarrel),
  etaStripEndcap_(etaStripEndcap),
  energyBarrel_(energyBarrel),
  energyEndcap_(energyEndcap)
{}

template<typename T1>
EcalPFClusterIsolation<T1>::~EcalPFClusterIsolation() 
{}

template<typename T1>
double EcalPFClusterIsolation<T1>::getSum(const T1Ref candRef, edm::Handle<reco::PFClusterCollection> clusterHandle) {
  
  drVeto2_ = -1.;
  float etaStrip = -1;
 
  if (fabs(candRef->eta()) < 1.479) {
    drVeto2_ = drVetoBarrel_*drVetoBarrel_;
    etaStrip = etaStripBarrel_;
  } else {
    drVeto2_ = drVetoEndcap_*drVetoEndcap_;
    etaStrip = etaStripEndcap_;
  }
   
  float etSum = 0;
  for (size_t i=0; i<clusterHandle->size(); i++) {
    reco::PFClusterRef pfclu(clusterHandle, i);

    if (fabs(candRef->eta()) < 1.479) {
      if (fabs(pfclu->pt()) < energyBarrel_)
	continue;
    } else {
      if (fabs(pfclu->energy()) < energyEndcap_)
	continue;
    }
    
    float dEta = fabs(candRef->eta() - pfclu->eta());
    if(dEta < etaStrip) continue;
    float dR2 = 0;
    if (not computedRVeto(candRef, pfclu, dR2)) continue;

    etSum += pfclu->pt();
  }

  return etSum;
}

template<typename T1>
bool EcalPFClusterIsolation<T1>::computedRVeto(T1Ref candRef, reco::PFClusterRef pfclu, float& dR2) {

  //float dR2 = deltaR2(candRef->eta(), candRef->phi(), pfclu->eta(), pfclu->phi());
  dR2 = deltaR2(candRef->eta(), candRef->phi(), pfclu->eta(), pfclu->phi());
  if(dR2 > (drMax_*drMax_))
    return false;

  if (candRef->superCluster().isNonnull()) {
    // Exclude clusters that are part of the candidate
    for (reco::CaloCluster_iterator it = candRef->superCluster()->clustersBegin(); it != candRef->superCluster()->clustersEnd(); ++it) {
      if ((*it)->seed() == pfclu->seed()) {
	return false;
      }
    }
  }

  return true;
}


template<>
bool EcalPFClusterIsolation<reco::GsfElectron>::computedRVeto(reco::GsfElectronRef candRef, reco::PFClusterRef pfclu, float& dR2) {

  dR2 = deltaR2(candRef->eta(), candRef->phi(), pfclu->eta(), pfclu->phi());
  if(dR2 > (drMax_*drMax_))
    return false;

  if (candRef->parentSuperCluster().isNonnull()) {
    // Exclude clusters that are part of the candidate
    for (reco::CaloCluster_iterator it = candRef->parentSuperCluster()->clustersBegin(); it != candRef->parentSuperCluster()->clustersEnd(); ++it) {
      if ((*it)->seed() == pfclu->seed()) {
	return false;
      }
    }
  }

  return true;
}

template<>
bool EcalPFClusterIsolation<reco::RecoChargedCandidate>::computedRVeto(T1Ref candRef, reco::PFClusterRef pfclu, float& dR2) {

  dR2 = deltaR2(candRef->eta(), candRef->phi(), pfclu->eta(), pfclu->phi());
  if(dR2 > (drMax_*drMax_) || dR2 < drVeto2_)
    return false;
  else
    return true;
}

template<>
double EcalPFClusterIsolation<reco::GsfElectron>::getSum(const reco::GsfElectronRef candRef, edm::Handle<reco::PFClusterCollection> clusterHandle) {
  

  if (fabs(candRef->superCluster()->position().eta())>2.5)
    return 0;
  //std::cout << "candRef: " << candRef->superCluster()->rawEnergy() << " " << candRef->superCluster()->position().eta() << " " << candRef->superCluster()->position().phi() << std::endl;
  // std::cout << "CLUSTER:" << std::endl;
  //  for (unsigned int i=0; i<candRef->superCluster()->size(); i++) 
  //    std::cout << candRef->superCluster()->printHitAndFraction(i) << std::endl;


  drVeto2_ = -1.;
  float etaStrip = -1;
 
  if (fabs(candRef->eta()) < 1.479) {
    drVeto2_ = drVetoBarrel_*drVetoBarrel_;
    etaStrip = etaStripBarrel_;
  } else {
    drVeto2_ = drVetoEndcap_*drVetoEndcap_;
    etaStrip = etaStripEndcap_;
  }

  //std::cout << "Cones: " << drVeto2_ << " " << etaStrip << " " << energyBarrel_ << " " << energyEndcap_ << " " << drMax_ << std::endl;

  float etSum = 0;
  for (size_t i=0; i<clusterHandle->size(); i++) {
    reco::PFClusterRef pfclu(clusterHandle, i);


    if (fabs(candRef->eta()) < 1.479) {
      if (fabs(pfclu->pt()) < energyBarrel_)
	continue;
    } else {
      if (fabs(pfclu->energy()) < energyEndcap_)
	continue;
    }
    
    float dEta = fabs(candRef->eta() - pfclu->eta());
    if(dEta < etaStrip) continue;

    float dR2 = 0;
    //auto dp = pfclu->phi()-candRef->superCluster()->position().phi(); 
    //if (dp>Float(M_PI)) 
    //  dp -= Float(2*M_PI);  
    //else if (dp < -Float(M_PI))
    //  dp += Float(2*M_PI);
    //
    //auto de = pfclu->eta()-candRef->superCluster()->position().eta();
    //float dr2 = de*de+dp+dp;

    //if (dr2 < 0.3*0.3)
    //  std::cout << "dRPFSI:" << pfclu->eta()-candRef->superCluster()->position().eta() << " " << dp << " " << pfclu->energy() << " " << std::endl; 
    if (not computedRVeto(candRef, pfclu, dR2)) continue;
//std::cout << "PFClusters: " << pfclu->energy() << " " << pfclu->eta() << " " << pfclu->phi() << " " << dR2 << std::endl;
//std::cout << "HITS:" << std::endl;
//for (unsigned int i=0; i<pfclu->size(); i++) 
//	 std::cout << pfclu->printHitAndFraction(i) << std::endl;

    //std::cout << "dRPFNO:" << pfclu->eta()-candRef->superCluster()->position().eta() << " " << dp << " " << pfclu->energy() << " " << std::endl; 
    etSum += pfclu->pt();
  }

  //std::cout << "iso: " << etSum << std::endl;
  return etSum;
}

template<>
double EcalPFClusterIsolation<reco::RecoEcalCandidate>::getSum(const reco::RecoEcalCandidateRef candRef, edm::Handle<reco::PFClusterCollection> clusterHandle) {
  

  if (fabs(candRef->superCluster()->position().eta())>2.5)
    return 0;
  //std::cout << "candRef: " << candRef->superCluster()->rawEnergy() << " " << candRef->superCluster()->position().eta() << " " << candRef->superCluster()->position().phi() << std::endl;
  //  std::cout << "CLUSTER:" << std::endl;
  //  for (unsigned int i=0; i<candRef->superCluster()->size(); i++) 
  //    std::cout << candRef->superCluster()->printHitAndFraction(i) << std::endl;
    

  drVeto2_ = -1.;
  float etaStrip = -1;
 
  if (fabs(candRef->eta()) < 1.479) {
    drVeto2_ = drVetoBarrel_*drVetoBarrel_;
    etaStrip = etaStripBarrel_;
  } else {
    drVeto2_ = drVetoEndcap_*drVetoEndcap_;
    etaStrip = etaStripEndcap_;
  }

  //std::cout << "Cones: " << drVeto2_ << " " << etaStrip << " " << energyBarrel_ << " " << energyEndcap_ << " " << drMax_ << std::endl;   

  float etSum = 0;
  for (size_t i=0; i<clusterHandle->size(); i++) {
    reco::PFClusterRef pfclu(clusterHandle, i);
    
    if (fabs(candRef->eta()) < 1.479) {
      if (fabs(pfclu->pt()) < energyBarrel_)
	continue;
    } else {
      if (fabs(pfclu->energy()) < energyEndcap_)
	continue;
    }
    
    float dEta = fabs(candRef->eta() - pfclu->eta());
    if(dEta < etaStrip) continue;
    
    float dR2 = 0;
    if (not computedRVeto(candRef, pfclu, dR2)) continue;

    //std::cout << "PFClusters: " << pfclu->energy() << " " << pfclu->eta() << " " << pfclu->phi() << " " << dR2 << std::endl;
    //std::cout << "HITS:" << std::endl;
    //for (unsigned int i=0; i<pfclu->size(); i++) 
    //  std::cout << pfclu->printHitAndFraction(i) << std::endl;

    etSum += pfclu->pt();
  }
  
  //std::cout << "iso: " << etSum << std::endl;
  return etSum;
}

template class EcalPFClusterIsolation<reco::RecoEcalCandidate>;
template class EcalPFClusterIsolation<reco::RecoChargedCandidate>;
template class EcalPFClusterIsolation<reco::Photon>;
template class EcalPFClusterIsolation<reco::GsfElectron>;
