dataset_1st = '/GluGluHToBB_M125_13TeV_powheg_pythia8/RunIISummer17MiniAOD-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/MINIAODSIM'
dataset_2nd = '/GluGluHToBB_M125_13TeV_powheg_pythia8/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/GEN-SIM-RAW'
name = 'HLT_EDMTuple_DoubleBTag_Hbb_Signal_v0p5'


from CRABClient.UserUtilities import config
config = config()

config.General.requestName = name+"_"+dataset_1st.split('/')[1]
config.General.workArea = 'crab_'+name
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.numCores = 4
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HLT_tunedcontent_dump.py'

config.Data.inputDataset = dataset_1st
config.Data.secondaryInputDataset = dataset_2nd
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 2000
config.Data.outLFNDirBase = '/store/user/htholen/' + name
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName
# config.Data.allowNonValidInputDataset = True
# config.Data.ignoreLocality = True

config.Site.blacklist = ['T0_*']
# config.Site.whitelist = ['T2_DE_DESY']
config.Site.storageSite = "T2_DE_DESY"
