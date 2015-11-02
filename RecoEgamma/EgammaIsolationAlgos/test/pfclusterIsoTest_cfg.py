import FWCore.ParameterSet.Config as cms

process = cms.Process("eleIso")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_74_V11'

process.load("Configuration.EventContent.EventContent_cff")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_7_4_8_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V11_mulTrh-v1/00000/001C7DBD-583C-E511-9107-0025905A60B2.root',
        '/store/relval/CMSSW_7_4_8_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V11_mulTrh-v1/00000/124BE41A-493C-E511-B909-0026189438DF.root',
        '/store/relval/CMSSW_7_4_8_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V11_mulTrh-v1/00000/AC6B87BB-953C-E511-8A84-002618943918.root',
        '/store/relval/CMSSW_7_4_8_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V11_mulTrh-v1/00000/C438BB08-5D3C-E511-9726-0025905A6084.root',
        '/store/relval/CMSSW_7_4_8_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V11_mulTrh-v1/00000/C89CDBD4-443C-E511-B817-0026189438A2.root',
        '/store/relval/CMSSW_7_4_8_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V11_mulTrh-v1/00000/F838C1BD-953C-E511-A59C-0025905A7786.root',
    ),
                            )
process.out = cms.OutputModule("PoolOutputModule",
                               process.FEVTSIMEventContent,
                               fileName = cms.untracked.string('file:eleIso.root')
)

process.out.outputCommands.append('drop *')
process.out.outputCommands.append('keep *_gedGsfElectrons_*_*')
process.out.outputCommands.append('keep *_gedGsfElectronCores_*_*')
process.out.outputCommands.append('keep particleFlowSuperClusterECAL_particleFlowSuperClusterECALBarrel_*_*')
process.out.outputCommands.append('keep particleFlowSuperClusterECAL_particleFlowSuperClusterECALEndcapWithPreshower_*_*')
#process.out.outputCommands.append('keep *_slimmedElectrons_*_*')
#process.out.outputCommands.append('keep *_photons_*_*')
process.out.outputCommands.append('keep *_*_*_eleIso')
process.out.outputCommands.append('keep *_genParticles_*_*')
process.out.outputCommands.append('keep *_addPileupInfo_*_*')
#process.out.outputCommands.append('keep EcalRecHitsSorted_*_*_*')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff")
process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
process.particleFlowRecHitECAL.producers[0].src=cms.InputTag("reducedEcalRecHitsEB")
process.particleFlowRecHitECAL.producers[1].src=cms.InputTag("reducedEcalRecHitsEE")
process.particleFlowRecHitPS.producers[0].src=cms.InputTag("reducedEcalRecHitsES")
process.particleFlowRecHitHBHE.producers[0].src=cms.InputTag("reducedHcalRecHits","hbhereco") 

process.load("RecoEgamma.EgammaIsolationAlgos.pfClusterIsolation_cfi")

process.plotDistr = cms.EDAnalyzer("plotDistr",
                                   OutputFileName = cms.string("ele_reco.root"),
                                   hitEBLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
                                   hitEELabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
                                   isData = cms.bool(False),
                                   activateNewClustering = cms.bool(False),
                                   activateOldClustering = cms.bool(False),
                                   saveReco = cms.bool(True),
                                   saveUnseeded = cms.bool(True),
                                   trgResults = cms.InputTag("TriggerResults","","TEST"),
                                   trgSelection = cms.vstring("HLT_Ele27WP85_Gsf_v1", "HLT_Ele20WP60_Ele8_Mass55_v1","HLT_Ele25WP60_SC4_Mass5_v1",)
                                   )

process.p1 = cms.Path(
#    process.particleFlowClusterWithoutHO +
    process.pfClusterIsolationSequence +
    process.plotDistr
)

#process.outpath = cms.EndPath(process.out)
