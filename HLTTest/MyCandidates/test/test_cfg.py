import FWCore.ParameterSet.Config as cms

process = cms.Process("ExREG")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V10::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/afiqaize/public/files/rootfiles/output_data.root')
                            fileNames = cms.untracked.vstring('file:../../../output_sig_ElEl_final.root')
                            #fileNames = cms.untracked.vstring('file:../../../outputA.root')
                            )

process.plotDistr = cms.EDAnalyzer("plotDistr",
                                   OutputFileName = cms.string("test_sig.root"),
                                   isData = cms.bool(False),
                                   activateNewClustering = cms.bool(True),
                                   activateOldClustering = cms.bool(False),
                                   saveReco = cms.bool(False),
                                   saveUnseeded = cms.bool(True),
                                   trgResults = cms.InputTag("TriggerResults","","TEST"),
                                   trgSelection = cms.vstring("HLT_Ele20WP60_Ele8_Mass55_v1",)
                                   )
process.p = cms.Path(process.plotDistr)

                                        
