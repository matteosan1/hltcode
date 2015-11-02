from DataFormats.FWLite import Handle, Events
import ROOT

photonH = Handle('std::vector<reco::GsfElectron>')
photonLabel = ('gedGsfElectrons')

vmapH = Handle('edm::ValueMap<float>')
vmapLabel = ('electronEcalPFClusterIsolationProducer')

events = Events("eleIso.root")

for ev in events:
    ev.getByLabel(photonLabel, photonH)
    eles = photonH.product()

    ev.getByLabel(vmapLabel, vmapH)
    vmap = vmapH.product()

    for i in xrange(len(eles)):
        print eles[i].pt(), vmap(eles[i])
        
