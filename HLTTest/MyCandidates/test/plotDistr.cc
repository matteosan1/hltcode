// -*- C++ -*-
//
// Package:    MyCandidates
// Class:      MyCandidates
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Thu Feb 14 14:06:52 CET 2013
// $Id$
//
//


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <string>
#include <iostream>

class plotDistr : public edm::EDAnalyzer {
public:
  explicit plotDistr(const edm::ParameterSet&);
  ~plotDistr();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  int findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::ElectronCollection> hltEleH);
  void mcTruth(edm::Handle<reco::GenParticleCollection> genParticleH);

private:
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& setup);
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag hitEBLabel_, hitEELabel_;
  edm::InputTag rhitEBLabel_, rhitEELabel_;
  edm::InputTag trigResultsTag_;
  HLTConfigProvider hltConfig_; 
  std::vector<std::string> pathNames_;

  std::string outputFileName;
  TFile* f;
  TTree* t;
  
  Float_t rho;
  Int_t passSel[3];
  Int_t n, npf, gp_n, reco_n, nuns;
  Float_t etawidth[10], phiwidth[10];
  Float_t etawidthpf[10], phiwidthpf[10];
  Float_t eraw[10];
  Float_t e[10];
  Float_t e5x5pf[10];
  Float_t et[10];
  Float_t se[10];
  Float_t erawpf[10];
  Float_t epf[10];  
  Float_t eoppf[10];  
  Float_t etpf[10];  
  Float_t sepf[10];  
  Float_t eta[10];
  Float_t etapf[10];  
  Float_t euns[20];
  Float_t etuns[20];
  Float_t etauns[20];  
  Float_t seta[10];
  Float_t setapf[10];  
  Float_t phi[10];
  Float_t phipf[10];  
  Float_t phiuns[20];  
  Float_t sphi[10];
  Float_t sphipf[10];  
  Float_t ecal[10];
  Float_t ecalpf[10];
  Float_t sieie[10];
  Float_t sieiepf[10];
  //Float_t eop[10];
  //Float_t eoppf[10];
  Float_t deta[10];
  Float_t detapf[10];
  Float_t dphi[10];
  Float_t dphipf[10];
  Float_t tkptpf[10];
  Float_t tketapf[10];
  Float_t tkphipf[10];
  Float_t tkpt[10];
  Float_t tketa[10];
  Float_t tkphi[10];
  Float_t hoe[10];
  Float_t hcal[10];
  Float_t tkiso[10];
  Float_t hoepf[10];
  Float_t hcalpf[10];
  Float_t tkisopf[10];
  Float_t chiso[10];
  Float_t phiso[10];
  Float_t neiso[10];
  Int_t nchi2pf[10];
  Float_t chi2pf[10][100];
  Int_t mishitspf[10];
  Int_t hitspf[10];

  Float_t gp_pt[10];
  Float_t gp_eta[10];
  Float_t gp_phi[10];

  Float_t reco_e5x5[10];
  Float_t reco_et[10];
  Float_t reco_pt[10];
  Float_t reco_eta[10];
  Float_t reco_phi[10];
  Float_t reco_ecaliso[10];
  Float_t reco_tkiso[10];
  Float_t reco_sieie[10];
  Float_t reco_deta[10];
  Float_t reco_dphi[10];
  Float_t reco_eop[10];
  Float_t reco_hoe[10];
  Int_t reco_nchi2[10];
  Float_t reco_chi2[10][100];
  
  float truePU;
  int bxPU[16];
  Int_t nvtx;

  bool isData;
  bool saveReco;
  bool saveUnseeded;
  bool newClustering;
  bool oldClustering;

  TH1F* timeEB, *timeEE, *chi2EB, *chi2EE;
  TH2F* timeEB2D, *timeEE2D;
};

void plotDistr::beginRun(const edm::Run& run,const edm::EventSetup& setup) {
  //std::cout <<"begining run "<<std::endl;
  //bool changed = false;
  //hltConfig_.init(run,setup,trigResultsTag_.process(),changed); //as we need the orginal HLT config...
  //std::cout <<"table name "<<hltConfig_.tableName()<<std::endl;
}


plotDistr::plotDistr(const edm::ParameterSet& iConfig) {
  outputFileName  = iConfig.getParameter<std::string>("OutputFileName");
  isData          = iConfig.getParameter<bool>("isData");
  newClustering   = iConfig.getParameter<bool>("activateNewClustering");
  oldClustering   = iConfig.getParameter<bool>("activateOldClustering");
  saveReco        = iConfig.getParameter<bool>("saveReco");
  saveUnseeded    = iConfig.getParameter<bool>("saveUnseeded");
  pathNames_      = iConfig.getParameter<std::vector<std::string>>("trgSelection");
  trigResultsTag_ = iConfig.getParameter<edm::InputTag>("trgResults");
  hitEBLabel_     = iConfig.getParameter<edm::InputTag>("hitEBLabel");
  hitEELabel_     = iConfig.getParameter<edm::InputTag>("hitEELabel");
  rhitEBLabel_     = iConfig.getParameter<edm::InputTag>("rhitEBLabel");
  rhitEELabel_     = iConfig.getParameter<edm::InputTag>("rhitEELabel");

  timeEB   = new TH1F("timeEB",   "", 400, -100, 100);
  timeEE   = new TH1F("timeEE",   "", 400, -100, 100);
  timeEB2D = new TH2F("timeEB2D", "", 400, -100, 100, 100, 0, 100);
  timeEE2D = new TH2F("timeEE2D", "", 400, -100, 100, 100, 0, 100);

  chi2EB   = new TH1F("chi2EB", "", 1000, 0, 1000);
  chi2EE   = new TH1F("chi2EE", "", 1000, 0, 1000);
}

plotDistr::~plotDistr() 
{}

int plotDistr::findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::ElectronCollection> hltEleH) {

  int index = -1;
  for (unsigned int i=0; i<hltEleH->size(); i++) {
    reco::ElectronRef cand(hltEleH, i);
    if (cand->superCluster() == ref->superCluster()) 
      return i;
  }

  return index;
}

void plotDistr::mcTruth(edm::Handle<reco::GenParticleCollection> gpH) {
  
  gp_n = 0;
  
  for(size_t i = 0; i < gpH->size(); ++i) {
    
    const reco::GenParticleRef gp(gpH, i);
    
    if ((gp->status() >=20 and gp->status() <=29) and abs(gp->pdgId()) == 11) {
      if (gp->pt() > 0.) {
	gp_pt[gp_n]  = gp->pt();
	gp_eta[gp_n] = gp->eta();
	gp_phi[gp_n] = gp->phi();
	gp_n++;
      }
    }
  }
}

void plotDistr::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {
  
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  //const CaloGeometry *geometry = pGeometry.product();
  
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  const CaloTopology *topology = pTopology.product();
  
  //edm::Handle<edm::TriggerResults> trigResultsHandle;
  //event.getByLabel(trigResultsTag_,trigResultsHandle);

  //const edm::TriggerResults& trigResults = *trigResultsHandle;
  //const edm::TriggerNames& trigNames = event.triggerNames(trigResults);  

  //for(size_t pathNr=0;pathNr<pathNames_.size();pathNr++){
  //  passSel[pathNr] = 0;
  //  size_t pathIndex = trigNames.triggerIndex(pathNames_[pathNr]);
  //  if(pathIndex<trigResults.size() &&  trigResults.accept(pathIndex)) 
  //    passSel[pathNr] = 1;
  //}
  
  truePU = 9999.;
  for (int i=0; i<16; i++)
    bxPU[i] = 9999;

  edm::Handle<reco::GenParticleCollection> gpH;
  if (!isData) {
    event.getByLabel("genParticles", gpH);
    mcTruth(gpH);
    
    edm::Handle<std::vector<PileupSummaryInfo>> puH;
    event.getByLabel("addPileupInfo", puH);
    truePU = (*puH)[0].getTrueNumInteractions();
    for (unsigned int j=0; j<puH->size(); j++)
      bxPU[j] = (*puH)[j].getPU_NumInteractions();
  }

  edm::Handle<EcalRecHitCollection> rhitEBH, rhitEEH;
  event.getByLabel(rhitEBLabel_, rhitEBH);
  event.getByLabel(rhitEELabel_, rhitEEH);

  edm::Handle<EcalRecHitCollection> hitEBH, hitEEH;
  event.getByLabel(hitEBLabel_, hitEBH);
  event.getByLabel(hitEELabel_, hitEEH);

  EcalRecHitCollection ec;
  EcalRecHitCollection ec2;


  if (!hitEBH.failedToGet()) {
    ec = *(hitEBH.product());
    //for(EcalRecHitCollection::const_iterator recHit = (*ec).begin(); recHit != (*ec).end(); ++recHit) {
    //  timeEB->Fill(recHit->time());
    //  timeEB2D->Fill(recHit->time(), recHit->energy());
    //  chi2EB->Fill(recHit->chi2());
    //}
  }
  if (!hitEEH.failedToGet()) {
    ec2 = *(hitEEH.product());
    //for(EcalRecHitCollection::const_iterator recHit = (*ec2).begin(); recHit != (*ec2).end(); ++recHit) {
    ////if (fabs(recHit->time()) < 10.) {
    //timeEE->Fill(recHit->time());
    //timeEE2D->Fill(recHit->time(), recHit->energy());
    //chi2EE->Fill(recHit->chi2());
    //}
  }

  edm::Handle<std::vector<reco::Electron> > eH;
  edm::Handle<std::vector<reco::RecoEcalCandidate> > cH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> sieieMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> ecalMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> hcalMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> hoeMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> tkisoMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> detaMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> dphiMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> eopMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> mishitsMapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> chi2MapH;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> hitsMapH;
  edm::Handle<reco::VertexCollection> vtxH;

  if (saveUnseeded) {
    event.getByLabel(edm::InputTag("hltEgammaCandidatesUnseeded","", "TEST"), cH);

    nuns = 0;
    if (!cH.failedToGet()) {

      for (unsigned int i=0; i<cH->size(); i++) {
	if (nuns == 20)
	  continue;
	euns[nuns] = 9999.;
	etuns[nuns] = 9999.;   
	etauns[nuns] = 9999.;
	phiuns[nuns] = 9999.;
	
	reco::RecoEcalCandidateRef c(cH, i);
	euns[nuns] = c->superCluster()->rawEnergy();
	etuns[nuns] = c->superCluster()->rawEnergy()*sin(c->theta());
	etauns[nuns] = c->eta();
	phiuns[nuns] = c->phi();
	nuns++;
      }      
    }
  }

  if (saveReco) {
    edm::Handle<std::vector<reco::GsfElectron> > elH;
    event.getByLabel("gedGsfElectrons", elH);

    edm::Handle<edm::ValueMap<float> > vmH;
    event.getByLabel("electronEcalPFClusterIsolationProducer", vmH);
    edm::ValueMap<float> vm = *(vmH.product());

    edm::Handle<edm::ValueMap<float> > vmTkH;
    event.getByLabel("hltTkIsoProducer", vmTkH);
    edm::ValueMap<float> vmTk = *(vmTkH.product());

    reco_n = 0;
    for (unsigned int i=0; i<elH->size(); i++) {
      if (reco_n == 9)
	continue;
      reco_et[reco_n] = 9999.;
      reco_e5x5[reco_n] = 9999.;
      reco_eta[reco_n] = 9999.;
      reco_phi[reco_n] = 9999.;
      reco_pt[reco_n] = 9999.;
      reco_ecaliso[reco_n] = 9999.;
      reco_tkiso[reco_n] = 9999.;
      reco_sieie[reco_n] = 9999.;
      reco_deta [reco_n] = 9999.;
      reco_dphi [reco_n] = 9999.;
      reco_eop  [reco_n] = 9999.;
      reco_hoe  [reco_n] = 9999.;

      reco::GsfElectronRef el(elH, i);
      reco_et[reco_n] = el->superCluster()->rawEnergy()*sin(el->superCluster()->position().theta());
      reco_e5x5[reco_n] = el->full5x5_e5x5();
      reco_pt[reco_n] = el->pt();
      reco_eta[reco_n] = el->superCluster()->eta();
      reco_phi[reco_n] = el->superCluster()->phi();
      reco_ecaliso[reco_n] = vm[el]; 
      reco_tkiso[reco_n] = vmTk[el]; 
      reco_sieie[reco_n] = el->full5x5_sigmaIetaIeta();
      reco_deta [reco_n] = el->deltaEtaSuperClusterTrackAtVtx();
      reco_dphi [reco_n] = el->deltaPhiSuperClusterTrackAtVtx();
      reco_eop  [reco_n] = el->eSuperClusterOverP();
      reco_hoe  [reco_n] = el->full5x5_hcalOverEcal();
      
      const std::vector< std::pair<DetId, float> > hits = el->superCluster()->hitsAndFractions();
      int subdet = el->superCluster()->seed()->hitsAndFractions()[0].first.subdetId();

      reco_nchi2[reco_n] = 0;
      for (unsigned int h=0; h<el->superCluster()->size(); h++) {
	if (reco_nchi2[reco_n] == 100) {
	  std::cout << "Too many hits !" << std::endl;
	  break; 
	}
	EcalRecHitCollection::const_iterator hit;
	if (subdet == EcalBarrel) 
	  hit = rhitEBH->find(hits[h].first);
	else
	  hit =  rhitEEH->find(hits[h].first);
	
	if ((hit->checkFlag(EcalRecHit::kGood) || hit->checkFlag(EcalRecHit::kOutOfTime)) && (hit->energy() > 1)) {
	  reco_chi2[reco_n][reco_nchi2[reco_n]] = hit->chi2();
	  reco_nchi2[reco_n]++;
	}
      }

      reco_n++;
    }
  }

  //if (!isData) {
  event.getByLabel(edm::InputTag("hltPixelVertices"), vtxH);
  if (!vtxH.failedToGet())
    nvtx = vtxH->size();
  else
    nvtx = 0;

  edm::Handle<double> rhoH;
  event.getByLabel(edm::InputTag("hltFixedGridRhoFastjetAllCaloForMuons"), rhoH);
  if (!rhoH.failedToGet())
    rho = *(rhoH.product());
  else
    rho = 9999.;

  if (newClustering) {
    event.getByLabel( edm::InputTag("hltEgammaCandidates", "", "TEST"), cH);
    
    npf = 0;
    if (!cH.failedToGet()) {
    
      const reco::RecoEcalCandidateIsolationMap* sieieMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaClusterShape:sigmaIEtaIEta5x5"), sieieMapH);
      if (!sieieMapH.failedToGet())  
	sieieMapPF = sieieMapH.product();
      
      const reco::RecoEcalCandidateIsolationMap* ecalMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaEcalPFClusterIso"), ecalMapH);
      if (!ecalMapH.failedToGet()) 
	ecalMapPF = ecalMapH.product(); 
      
      const reco::RecoEcalCandidateIsolationMap* hcalMapPF  = 0;
      event.getByLabel(edm::InputTag("hltEgammaHcalPFClusterIso"), hcalMapH);
      if (!hcalMapH.failedToGet()) 
	hcalMapPF = hcalMapH.product(); 
      
      const reco::RecoEcalCandidateIsolationMap* hoeMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaHoverE"), hoeMapH);
      if (!hoeMapH.failedToGet()) 
	hoeMapPF = hoeMapH.product(); 
      
      const reco::RecoEcalCandidateIsolationMap* tkisoMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaEleGsfTrackIso"), tkisoMapH);
      if (!tkisoMapH.failedToGet())
	tkisoMapPF = tkisoMapH.product();          
      
      const reco::RecoEcalCandidateIsolationMap* detaMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaGsfTrackVars:Deta"), detaMapH);
      if (!detaMapH.failedToGet())
	detaMapPF = detaMapH.product();
      
      const reco::RecoEcalCandidateIsolationMap* dphiMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaGsfTrackVars:Dphi"), dphiMapH);
      if (!dphiMapH.failedToGet()) 
	dphiMapPF = dphiMapH.product();
      
      const reco::RecoEcalCandidateIsolationMap* eopMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaGsfTrackVars:OneOESuperMinusOneOP"), eopMapH);
      if (!eopMapH.failedToGet()) 
	eopMapPF = eopMapH.product();

//const reco::RecoEcalCandidateIsolationMap* chi2MapPF = 0;
//event.getByLabel(edm::InputTag("hltEgammaGsfTrackVars:Chi2"), chi2MapH);
//if (!chi2MapH.failedToGet()) 
//	chi2MapPF = chi2MapH.product();

      const reco::RecoEcalCandidateIsolationMap* mishitsMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaGsfTrackVars:MissingHits"), mishitsMapH);
      if (!mishitsMapH.failedToGet()) 
	mishitsMapPF = mishitsMapH.product();

      const reco::RecoEcalCandidateIsolationMap* hitsMapPF = 0;
      event.getByLabel(edm::InputTag("hltEgammaGsfTrackVars:ValidHits"), hitsMapH);
      if (!hitsMapH.failedToGet()) 
	hitsMapPF = hitsMapH.product();
      
      for (unsigned int i=0; i<cH->size(); i++) {
	//std::cout << "Canddiates" << i << std::endl;
	if (npf == 9)
	  continue;
	e5x5pf[npf] = 9999.;
	etawidthpf[npf] = 9999.;
	phiwidthpf[npf] = 9999.;
	epf[npf] = 9999.;
	eoppf[npf] = 9999.;
	etpf[npf] = 9999.;
	sepf[npf] = 9999.;
	erawpf[npf] = 9999.;
	etapf[npf] = 9999.;
	phipf[npf] = 9999.;
	setapf[npf] = 9999.;
	sphipf[npf] = 9999.;
	sieiepf[npf] = 9999.;
	ecalpf[npf] = 9999;
	detapf[npf] = 9999.;
	dphipf[npf] = 9999.;
	tkptpf[npf] = 9999.;
	tketapf[npf] = 9999.;
	tkphipf[npf] = 9999.;
	hcalpf[npf] = 9999.;
	hoepf[npf] = 9999.;
	tkisopf[npf] = 9999.;
	mishitspf[npf] = 9999.;
	//chi2pf[npf] = 9999.;
	hitspf[npf] = 9999.;
	
	reco::RecoEcalCandidateRef c(cH, i);
	//std::cout << c->eta() << " " <<  c->phi() << " " << c->et() << std::endl;
	//if (c->energy()*sin(c->theta()) < 20.)
	//  continue;
	
	std::cout << "ECCOMI" << std::endl;
	const std::vector< std::pair<DetId, float> > hits = c->superCluster()->hitsAndFractions();
	int subdet = c->superCluster()->seed()->hitsAndFractions()[0].first.subdetId();
	std::cout << "ECCOMI2" << std::endl;
	nchi2pf[npf] = 0;
	for (unsigned int h=0; h<hits.size(); h++) {
	    if (nchi2pf[npf] == 100) {
	      std::cout << "Too many hits HLT !" << std::endl;
	      break; 
	    }
	    
	    EcalRecHitCollection::const_iterator hit;
	    if (subdet == EcalBarrel)
	      hit = hitEBH->find(hits[h].first);
	    else
	      hit = hitEEH->find(hits[h].first);
	    
	    if ((hit->checkFlag(EcalRecHit::kGood) || hit->checkFlag(EcalRecHit::kOutOfTime)) and (hit->energy() > 1)) {
	      chi2pf[npf][nchi2pf[npf]] = hit->chi2();
	      nchi2pf[npf]++;
	    }
	  }
	       

	//erawpf[npf] = c->superCluster()->rawEnergy();
	//etawidthpf[npf] = c->superCluster()->etaWidth();
	//phiwidthpf[npf] = c->superCluster()->phiWidth();
	std::cout << "ECCOMI3" << std::endl;
	epf[npf] = c->superCluster()->rawEnergy();//*sin(c->theta());
	if (subdet == EcalBarrel)
	  e5x5pf[npf] = noZS::EcalClusterTools::e5x5(*(c->superCluster()->seed()), &ec, topology); 
	else
	  e5x5pf[npf] = noZS::EcalClusterTools::e5x5(*(c->superCluster()->seed()), &ec2, topology); 
	etpf[npf] = c->superCluster()->rawEnergy()*sin(c->superCluster()->position().theta());
	//sepf[npf] = c->superCluster()->seed()->energy();//*sin(c->theta());
	etapf[npf] = c->superCluster()->position().eta();
	phipf[npf] = c->superCluster()->position().phi();
	//setapf[npf] = c->superCluster()->seed()->eta();
	//sphipf[npf] = c->superCluster()->seed()->phi();

	if (sieieMapPF != 0)
	  sieiepf[npf] = (*sieieMapPF)[c];
	if (ecalMapPF != 0)
	  ecalpf[npf] = (*ecalMapPF)[c]; 
	
	if (hoeMapPF != 0) 
	  hoepf[npf] = (*hoeMapPF)[c]; 
	
	if (hcalMapPF != 0) 
	  hcalpf[npf] = (*hcalMapPF)[c]; 
	
	if (detaMapPF != 0)
	  detapf[npf] = ((*detaMapPF)[c]);
	
	if (dphiMapPF != 0)
	  dphipf[npf] = fabs((*dphiMapPF)[c]);
	
	if (eopMapPF != 0)
	  eoppf[npf] = fabs((*eopMapPF)[c]);
	
	if (mishitsMapPF != 0)
	  mishitspf[npf] = (*mishitsMapPF)[c];

	if (hitsMapPF != 0)
	  hitspf[npf] = (*hitsMapPF)[c];

	//if (chi2MapPF != 0)
	//  chi2pf[npf] = (*chi2MapPF)[c];
	
	if (tkisoMapPF != 0)
	  tkisopf[npf] = fabs((*tkisoMapPF)[c]);
	
	npf++;
      }
    }
  } 

  t->Fill();
}
 
void plotDistr::beginJob() {
  f = new TFile(outputFileName.c_str(), "recreate");
  t = new TTree("tree", "tree");
  
  if (oldClustering) {
    t->Branch("n",  &n, "n/I");
    t->Branch("ewidth", &etawidth, "ewidth[n]/F");
    t->Branch("pwidth", &phiwidth, "pwidth[n]/F");
    t->Branch("e", &e, "e[n]/F");
    t->Branch("et", &et, "et[n]/F");
    t->Branch("se", &se, "se[n]/F");
    t->Branch("eraw", &eraw, "eraw[n]/F");
    t->Branch("eta", &eta, "eta[n]/F");
    t->Branch("phi", &phi, "phi[n]/F");
    t->Branch("seta", &seta, "seta[n]/F");
    t->Branch("sphi", &sphi, "sphi[n]/F");
    t->Branch("sieie", &sieie, "sieie[n]/F");
    t->Branch("ecal", &ecal, "ecal[n]/F");
    //t->Branch("eop", &eop, "eop[n]/F");
    t->Branch("dphi", &dphi, "dphi[n]/F");
    t->Branch("deta", &deta, "deta[n]/F");
    t->Branch("tkpt",    &tkpt,  "tkpt[n]/F");
    t->Branch("tketa",   &tketa, "tketa[n]/F");
    t->Branch("tkphi",   &tkphi, "tkphi[n]/F");
    t->Branch("hoe",   &hoe, "hoe[n]/F");
    t->Branch("hcal",   &hcal, "hcal[n]/F");
    t->Branch("tkiso",   &tkiso, "tkiso[n]/F");
  }
  
  t->Branch("passHLT", &passSel, "passHLT[3]/I");
  t->Branch("nvtx", &nvtx, "nvtx/I");
  t->Branch("rho",  &rho, "rho/F");
  
  if (newClustering) {
    t->Branch("npf",  &npf, "npf/I");
    //t->Branch("ewidthpf", &etawidthpf, "ewidthpf[npf]/F");
    //t->Branch("pwidthpf", &phiwidthpf, "pwidthpf[npf]/F");
    t->Branch("epf", &epf, "epf[npf]/F");
    t->Branch("e5x5pf", &e5x5pf, "e5x5pf[npf]/F");
    t->Branch("etpf", &etpf, "etpf[npf]/F");
    t->Branch("sepf", &sepf, "sepf[npf]/F");
    //t->Branch("erawpf", &erawpf, "erawpf[npf]/F");
    t->Branch("etapf", &etapf, "etapf[npf]/F");
    t->Branch("phipf", &phipf, "phipf[npf]/F");
    //t->Branch("setapf", &setapf, "setapf[npf]/F");
    //t->Branch("sphipf", &sphipf, "sphipf[npf]/F");
    t->Branch("sieiepf", &sieiepf, "sieiepf[npf]/F");
    t->Branch("ecalpf", &ecalpf, "ecalpf[npf]/F");
    //t->Branch("eop", &eoppf, "eoppf[npf]/F");
    t->Branch("dphipf", &dphipf, "dphipf[npf]/F");
    t->Branch("detapf", &detapf, "detapf[npf]/F");
    t->Branch("tkptpf",  &tkptpf, "tkptpf[npf]/F");
    t->Branch("tketapf", &tketapf, "tketapf[npf]/F");
    t->Branch("tkphipf", &tkphipf, "tkphipf[npf]/F");
    t->Branch("hoepf",   &hoepf, "hoepf[npf]/F");
    t->Branch("hcalpf",   &hcalpf, "hcalpf[npf]/F");
    t->Branch("tkisopf",   &tkisopf, "tkisopf[npf]/F");
    t->Branch("chiso",   &chiso, "chiso[npf]/F");
    t->Branch("phiso",   &phiso, "phiso[npf]/F");
    t->Branch("neiso",   &neiso, "neiso[npf]/F");
    t->Branch("eoppf", &eoppf, "eoppf[npf]/F");
    t->Branch("nchi2pf", &nchi2pf, "nchi2pf[npf]/I");
    t->Branch("chi2pf", &chi2pf, "chi2pf[npf][100]/F");
    t->Branch("mishitspf", &mishitspf, "mishitspf[npf]/I");
    t->Branch("hitspf", &hitspf, "hitspf[npf]/I");
  }
  
  if (saveUnseeded) {   
    t->Branch("nuns",  &nuns, "nuns/I");
    t->Branch("euns", &euns, "euns[nuns]/F");
    t->Branch("etuns", &etuns, "etuns[nuns]/F");
    t->Branch("etauns", &etauns, "etauns[nuns]/F");
    t->Branch("phiuns", &phiuns, "phiuns[nuns]/F");
  }

  if (saveReco) {
    t->Branch("recon",   &reco_n,   "recon/I");
    t->Branch("recoet",  &reco_et,  "recoet[recon]/F");
    t->Branch("recoe5x5",  &reco_e5x5,  "recoe5x5[recon]/F");
    t->Branch("recopt",  &reco_pt,  "recopt[recon]/F");
    t->Branch("recoeta", &reco_eta, "recoeta[recon]/F");
    t->Branch("recophi", &reco_phi, "recophi[recon]/F");
    t->Branch("recoecaliso", &reco_ecaliso, "recoecaliso[recon]/F"); 
    t->Branch("recotkiso", &reco_tkiso, "recotkiso[recon]/F"); 
    t->Branch("recosieie", &reco_sieie, "recosieie[recon]/F");
    t->Branch("recodeta", &reco_deta, "recodeta[recon]/F");
    t->Branch("recodphi", &reco_dphi, "recodphi[recon]/F");
    t->Branch("recoeop", &reco_eop, "recoeop[recon]/F");
    t->Branch("recohoe", &reco_hoe, "recohoe[recon]/F");
    t->Branch("reconchi2", &reco_nchi2, "reconchi2[recon]/I");
    t->Branch("recochi2", &reco_chi2, "recochi2[recon][100]/F");
  }

  if (!isData) {
    t->Branch("gpn",   &gp_n,   "gpn/I");
    t->Branch("gppt",  &gp_pt,  "gppt[gpn]/F");
    t->Branch("gpeta", &gp_eta, "gpeta[gpn]/F");
    t->Branch("gpphi", &gp_phi, "gpphi[gpn]/F");
    
    t->Branch("truePU", &truePU, "truePU/F");
    t->Branch("bxPU", &bxPU, "bxPU[16]/I");
  }
}

void plotDistr::endJob() {
  f->cd();
  t->Write();
  timeEB->Write();
  timeEE->Write();
  timeEB2D->Write();
  timeEE2D->Write();
  chi2EB->Write();
  chi2EE->Write();
  f->Close();
}

void plotDistr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(plotDistr);
