#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "TROOT.h"
#include "TFile.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Event.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <iostream>
#include <vector>
#endif

void recoTree() {

  TFile* file = TFile::Open("eleIso.root");

  fwlite::Event ev(file);
  
  for(ev.toBegin(); ! ev.atEnd(); ++ev) {
    edm::Handle<std::vector<reco::GsfElectron> > objs;
    ev.getByLabel("gedGsfElectrons", objs);
    //std::vector<reco::GsfElectron> eles = *(objs.product());

    fwlite::Handle<edm::ValueMap<float> > vmH;
    vmH.getByLabel(ev, "electronEcalPFClusterIsolationProducer");
    edm::ValueMap<float> vm = *(vmH.product());

    for (unsigned int i=0; i<objs->size(); i++) {
	reco::GsfElectronRef ref(objs, i);
      std::cout << ref->pt() << " " << vm[ref] << std::endl;
    }
  }
}
