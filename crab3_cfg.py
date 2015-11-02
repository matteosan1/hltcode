from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName = 'crab_data'
config.General.workArea = 'mc'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hlt25ns_mc_multifit.py'
#config.JobType.psetName = 'hlt25ns_mc_singlefit.py'
#config.JobType.psetName = 'hlt25ns_mc_std.py'
config.JobType.outputFiles = ['multiFit.root']

#config.Data.inputDataset = '/RelValZEE_13/CMSSW_7_4_12-PU25ns_74X_mcRun2_asymptotic_v2_v2-v1/GEN-SIM-RECO'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/GEN-SIM-RAW'
#/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15Digi74-HFscaleFlat10to30Asympt25nstsg_MCRUN2_74_V9-v1/GEN-SIM-RAW'
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15DR75-Asympt25nsFlat10to30RawReco_HCALDebug_75X_mcRun2_asymptotic_v8-v1/GEN-SIM-RECO'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '193093-193999' # '193093-194075'

#config.Data.outLFNDirBase = '/store/user/sani/hlt/DYNoPrefit_stdSC'
#config.Data.outLFNDirBase = '/store/user/sani/hlt/DYSingleFit_stdSC'
#config.Data.outLFNDirBase = '/store/user/sani/hlt/hcalScaleDY_std'
config.Data.outLFNDirBase = '/store/user/sani/hlt/hcalScaleDY_mutifit'

config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_May2015_Data_analysis'
config.Data.useParent = True

config.Site.storageSite = 'T2_US_UCSD'
