from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'data'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hlt_data.py'
config.JobType.outputFiles = ['multiFit.root'] # to be changed, as to be the same name given in the ntuple producer (myntuple.root)

#config.Data.inputDataset = '/DoubleEG/Run2015D-ZElectron-PromptReco-v3/RAW-RECO'
config.Data.inputDataset = '/DoubleEG/Run2015D-ZElectron-PromptReco-v4/RAW-RECO' # for the moment you can run just on this one since it has most of the RunD statistics

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'

config.Data.outLFNDirBase = '/store/user/sani/hlt/....' # to be changed, must be a valid path to your storage area

config.Data.publication = False
config.Data.useParent = False

config.Site.storageSite = 'T2_US_UCSD' # same as before, should point to your favorite T2
