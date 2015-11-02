import FWCore.ParameterSet.Config as cms

process = cms.Process("eleIso")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_74_V7'

process.load("Configuration.EventContent.EventContent_cff")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_7_4_12/RelValZEE_13/GEN-SIM-RECO/74X_mcRun2_asymptotic_v2-v1/00000/28632B14-515C-E511-A44C-003048FFD736.root'
        #'/store/relval/CMSSW_7_4_12/RelValZEE_13/MINIAODSIM/PU25ns_74X_mcRun2_asymptotic_v2_v2-v1/00000/26CE697F-735E-E511-8256-0025905A6088.root'

    )
)
process.out = cms.OutputModule("PoolOutputModule",
                               process.FEVTSIMEventContent,
                               fileName = cms.untracked.string('file:eleIso.root')
)

process.out.outputCommands.append('keep *')
#process.out.outputCommands.append('keep *_gsfElectrons_*_*')
#process.out.outputCommands.append('keep *_photons_*_*')
#process.out.outputCommands.append('keep *_*_*_eleIso')
#process.out.outputCommands.append('keep EcalRecHitsSorted_*_*_*')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

#   m CommonTools.RecoAlgos.recoTrackSelectorPSet_cfi import recoTrackSelectorPSet
 
process.recoTrackSelector = cms.EDProducer("RecoTrackSelector",
                                           src = cms.InputTag("generalTracks"),
                                           maxChi2 = cms.double(10000.0),
                                           tip = cms.double(120.0),
                                           minRapidity = cms.double(-5.0),
                                           lip = cms.double(300.0),
                                           ptMin = cms.double(0.1),
                                           maxRapidity = cms.double(5.0),
                                           quality = cms.vstring('loose'),
                                           algorithm = cms.vstring('initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep'
                                                                   'mixedTripletStep'),
                                           minHit = cms.int32(3),
                                           min3DHit = cms.int32(0),
                                           beamSpot = cms.InputTag("offlineBeamSpot"),
                                           copyExtras = cms.untracked.bool(True), ## copies also extras and rechits on RECO
                                           copyTrajectories = cms.untracked.bool(False) # don't set this to true on AOD!
#                                           src = cms.InputTag("generalTracks"),
#                                           maxChi2 = cms.double(10000.0),
#                                           tip = cms.double(120.0),
#                                           minRapidity = cms.double(-2.5),
#                                           lip = cms.double(300.0),
#                                           ptMin = cms.double(1.0),
#                                           maxRapidity = cms.double(2.5),
#                                           quality = cms.vstring('loose'),
#                                           algorithm = cms.vstring('initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep'
#                                                                   'mixedTripletStep'),
#                                           originalAlgorithm = cms.vstring(),
#                                           algorithmMaskContains = cms.vstring(),
#                                           minLayer = cms.int32(3),
#                                           min3DLayer = cms.int32(0),
#                                           minHit = cms.int32(0),
#                                           minPixelHit = cms.int32(0),
#                                           beamSpot = cms.InputTag("offlineBeamSpot"),
#                                           usePV = cms.bool(False),
#                                           vertexTag = cms.InputTag('offlinePrimaryVertices')                                       
                                           )

process.hltTkIsoProducer = cms.EDProducer('ElectronHLTTkIsoProducer',
                                          candidateProducer = cms.InputTag("gedGsfElectrons"),
                                          trackProducer = cms.InputTag("recoTrackSelector"),
                                          beamSpotProducer = cms.InputTag("offlineBeamSpot"),
                                          egTrkIsoPtMin = cms.double(1.),
                                          egTrkIsoConeSize = cms.double(0.3),
                                          egTrkIsoZSpan = cms.double(0.15),
                                          egTrkIsoRSpan = cms.double(99999.),
                                          egTrkIsoVetoConeSizeBarrel = cms.double(0.03),
                                          egTrkIsoVetoConeSizeEndcap = cms.double(0.03),
                                          egTrkIsoStripBarrel = cms.double(0.03),
                                          egTrkIsoStripEndcap = cms.double(0.03)
                                          )



process.p1 = cms.Path(
    process.recoTrackSelector + 
    process.hltTkIsoProducer
)

process.outpath = cms.EndPath(process.out)
