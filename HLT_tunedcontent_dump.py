import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST20180124162053")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:../../../GluGluHToBB_M125_13TeV_powheg_pythia8_miniaodsim_file1.root'),
    inputCommands = cms.untracked.vstring('keep *'),
    secondaryFileNames = cms.untracked.vstring(
        'file:../../../GluGluHToBB_M125_13TeV_powheg_pythia8_gen_sim_raw_file1.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/3E3CD4C8-1DAA-E711-820D-0025905B8568.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/825512FF-55AA-E711-A815-0CC47A7C340C.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/86C8D820-56AA-E711-A77C-0CC47A78A340.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150001/4427DE0B-5DAA-E711-95DB-0CC47A78A3F4.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/389530B6-62AA-E711-ADC0-0CC47A4D7698.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/50458A94-49AA-E711-8E53-0CC47A7C360E.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/507F7CAF-48AA-E711-8EDB-0025905B8610.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/DE3735B0-48AA-E711-A3DD-0025905A60F4.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/52643A20-38AA-E711-BCB5-0025905B8596.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/6A77374F-31AA-E711-BE95-0CC47A4D7692.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A4D42A6B-17AA-E711-9C02-0025905B85E8.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150001/88026120-56AA-E711-BC66-0025905A6066.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/00A54084-5CAA-E711-9BD4-003048FFD734.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/0E3B8255-14AA-E711-A144-0025905B85CA.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/38FB9483-5CAA-E711-A06E-0CC47A4D7636.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/DE096B85-5CAA-E711-865D-0025905A6104.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/2EE87096-30AA-E711-A47B-0CC47A745284.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/6CBF4329-3DAA-E711-A6FD-0CC47A4D764C.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/78F378A8-16AA-E711-90CC-0025905B85EC.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A4C20F3B-38AA-E711-B9DB-0025905B85CA.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/006B81F8-1DAA-E711-B9A7-0CC47A4C8ECE.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/441A1A58-37AA-E711-AFDC-0025905A6082.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/841461EA-0AAA-E711-806F-0025905B85FC.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/D40EC8A2-37AA-E711-9EEC-0025905A60AA.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/0CAF2DF7-68AA-E711-B628-0025905A48BC.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/5464C86A-2BAA-E711-9453-0CC47A4C8E82.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A25B1374-37AA-E711-A4C0-0025905A610C.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150001/1AC7C892-5CAA-E711-B48F-0CC47A7C34C8.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/98089770-37AA-E711-97DD-0CC47A7C34E6.root', 
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/GluGluHToBB_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150001/1831CC66-63AA-E711-871A-0CC47A4C8ECA.root'
    )
)
process.HLTConfigVersion = cms.PSet(
    tableName = cms.string('/users/xcoubez2018/doubleb/HLT_DoubleB_FirstFullImplementation/V4')
)

process.HLTIter0GroupedCkfTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0PSetTrajectoryFilterIT')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useHitsSplitting = cms.bool(False),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter0HighPtTkMuPSetTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(4),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0HighPtTkMuPSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTIter0HighPtTkMuPSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter0IterL3FromL1MuonGroupedCkfTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(10.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter0IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    foundHitBonus = cms.double(1000.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0IterL3FromL1MuonGroupedCkfTrajectoryFilterIT')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(True),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(1.0),
    maxCand = cms.int32(5),
    minNrOfHitsForRebuild = cms.int32(2),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0IterL3FromL1MuonGroupedCkfTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter0IterL3MuonGroupedCkfTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(10.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter0IterL3MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    foundHitBonus = cms.double(1000.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0IterL3MuonGroupedCkfTrajectoryFilterIT')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(True),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(1.0),
    maxCand = cms.int32(5),
    minNrOfHitsForRebuild = cms.int32(2),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0IterL3MuonGroupedCkfTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter0PSetTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator9'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter0PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTIter0PSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(4),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter1GroupedCkfTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter1ESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator16'),
    foundHitBonus = cms.double(5.0),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter1PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter1PSetTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter1ESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator16'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter1PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter1PSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter2GroupedCkfTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator16'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter2PSetTrajectoryFilterIT')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter2PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useHitsSplitting = cms.bool(False),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter2HighPtTkMuPSetTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter2HighPtTkMuESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter2HighPtTkMuPSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTIter2HighPtTkMuPSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(3),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter2IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter2HighPtTkMuESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    foundHitBonus = cms.double(1000.0),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter2IterL3FromL1MuonPSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter2IterL3FromL1MuonPSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(3),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter2IterL3MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter2HighPtTkMuESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    foundHitBonus = cms.double(1000.0),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter2IterL3MuonPSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTIter2IterL3MuonPSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(3),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter2PSetTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter2ESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator16'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter2PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTIter2PSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(1),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter3PSetTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter3ESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator16'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter3PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTIter3PSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTIter4PSetTrajectoryBuilderIT = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltIter4ESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator16'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    minNrOfHitsForRebuild = cms.untracked.int32(4),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTIter4PSetTrajectoryFilterIT')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTIter4PSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(6),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetCkf3HitTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(-1),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetCkfTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(-1),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetCkfTrajectoryFilterIterL3OI = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(10.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(-1),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(3.0),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetDetachedCkfTrajectoryBuilderForHI = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2MeasurementEstimator9'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedCkfTrajectoryFilterForHI')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(0.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialForHI'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOppositeForHI'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedCkfTrajectoryFilterForHI')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetDetachedCkfTrajectoryBuilderForHIGlobalPt8 = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2MeasurementEstimator9'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedCkfTrajectoryFilterForHIGlobalPt8')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(0.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialForHI'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOppositeForHI'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedCkfTrajectoryFilterForHIGlobalPt8')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetDetachedCkfTrajectoryFilterForHI = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(0.701),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(6),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetDetachedCkfTrajectoryFilterForHIGlobalPt8 = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(0.701),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(8.0),
    minimumNumberOfHits = cms.int32(6),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetDetachedQuadStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPDetachedQuadStepChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedQuadStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedQuadStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetDetachedQuadStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetDetachedStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetDetachedStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedStepTrajectoryFilterBase')
    ))
)

process.HLTPSetDetachedStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(2),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetDetachedTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPDetachedTripletStepChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedTripletStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetDetachedTripletStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetDetachedTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetGroupedCkfTrajectoryBuilderIterL3ForOI = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltSiStripClusters'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    deltaEta = cms.double(-1.0),
    deltaPhi = cms.double(-1.0),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    foundHitBonus = cms.double(1000.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetCkfTrajectoryFilterIterL3OI')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    propagatorProximity = cms.string('SteppingHelixPropagatorAny'),
    requireSeedHitsInRebuild = cms.bool(False),
    rescaleErrorIfFail = cms.double(1.0),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetCkfTrajectoryFilterIterL3OI')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True),
    useSeedLayer = cms.bool(False)
)

process.HLTPSetHighPtTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetHighPtTripletStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetHighPtTripletStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetHighPtTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(5),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetInitialCkfTrajectoryBuilderForHI = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    intermediateCleaning = cms.bool(False),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialForHI'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOppositeForHI'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetInitialCkfTrajectoryFilterForHI')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTPSetInitialCkfTrajectoryFilterForHI = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(6),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetInitialStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPInitialStepChi2ChargeMeasurementEstimator30'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetInitialStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(True),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(1),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetInitialStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetInitialStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetInitialStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(2),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetJetCoreStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetJetCoreStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(50),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetJetCoreStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetJetCoreStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetLowPtQuadStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetLowPtQuadStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(4),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetLowPtQuadStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetLowPtQuadStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetLowPtStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetLowPtStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(4),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetLowPtStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetLowPtStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(1),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetLowPtTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPLowPtTripletStepChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetLowPtTripletStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(4),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetLowPtTripletStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetLowPtTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetMixedStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeTightMeasurementEstimator16'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMixedStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialForMixedStep'),
    propagatorOpposite = cms.string('PropagatorWithMaterialForMixedStepOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMixedStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetMixedStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.4),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetMixedStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.05),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetMixedTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPMixedTripletStepChi2ChargeMeasurementEstimator16'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMixedTripletStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialForMixedStep'),
    propagatorOpposite = cms.string('PropagatorWithMaterialForMixedStepOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMixedTripletStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetMixedTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.4),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetMuTrackJpsiEffTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMuTrackJpsiEffTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTPSetMuTrackJpsiEffTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(9),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(1.0),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetMuTrackJpsiTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMuTrackJpsiTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTPSetMuTrackJpsiTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(8),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(10.0),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetMuonCkfTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('MuonCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    deltaEta = cms.double(-1.0),
    deltaPhi = cms.double(-1.0),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    intermediateCleaning = cms.bool(False),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    propagatorProximity = cms.string('SteppingHelixPropagatorAny'),
    rescaleErrorIfFail = cms.double(1.0),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMuonCkfTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSeedLayer = cms.bool(False)
)

process.HLTPSetMuonCkfTrajectoryBuilderSeedHit = cms.PSet(
    ComponentType = cms.string('MuonCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    deltaEta = cms.double(-1.0),
    deltaPhi = cms.double(-1.0),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    intermediateCleaning = cms.bool(False),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    propagatorProximity = cms.string('SteppingHelixPropagatorAny'),
    rescaleErrorIfFail = cms.double(1.0),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMuonCkfTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSeedLayer = cms.bool(True)
)

process.HLTPSetMuonCkfTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(-1),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetMuonTrackingRegionBuilder8356 = cms.PSet(
    DeltaEta = cms.double(0.2),
    DeltaPhi = cms.double(0.2),
    DeltaR = cms.double(0.2),
    DeltaZ = cms.double(15.9),
    EtaR_UpperLimit_Par1 = cms.double(0.25),
    EtaR_UpperLimit_Par2 = cms.double(0.15),
    Eta_fixed = cms.bool(False),
    Eta_min = cms.double(0.1),
    MeasurementTrackerName = cms.InputTag("hltESPMeasurementTracker"),
    OnDemand = cms.int32(-1),
    PhiR_UpperLimit_Par1 = cms.double(0.6),
    PhiR_UpperLimit_Par2 = cms.double(0.2),
    Phi_fixed = cms.bool(False),
    Phi_min = cms.double(0.1),
    Pt_fixed = cms.bool(False),
    Pt_min = cms.double(1.5),
    Rescale_Dz = cms.double(3.0),
    Rescale_eta = cms.double(3.0),
    Rescale_phi = cms.double(3.0),
    UseVertex = cms.bool(False),
    Z_fixed = cms.bool(True),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    input = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    maxRegions = cms.int32(2),
    precise = cms.bool(True),
    vertexCollection = cms.InputTag("pixelVertices")
)

process.HLTPSetPixelLessStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPPixelLessStepChi2ChargeMeasurementEstimator16'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelLessStepTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(4),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelLessStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetPixelLessStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetPixelLessStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.05),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetPixelPairCkfTrajectoryBuilderForHI = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator9ForHI'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelPairCkfTrajectoryFilterForHI')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialForHI'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOppositeForHI'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelPairCkfTrajectoryFilterForHI')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetPixelPairCkfTrajectoryBuilderForHIGlobalPt8 = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator9ForHI'),
    foundHitBonus = cms.double(5.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelPairCkfTrajectoryFilterForHIGlobalPt8')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialForHI'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOppositeForHI'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelPairCkfTrajectoryFilterForHIGlobalPt8')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HLTPSetPixelPairCkfTrajectoryFilterForHI = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(100),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(1.0),
    minimumNumberOfHits = cms.int32(6),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetPixelPairCkfTrajectoryFilterForHIGlobalPt8 = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(100),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(8.0),
    minimumNumberOfHits = cms.int32(6),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetPixelPairStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPPixelPairStepChi2ChargeMeasurementEstimator9'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelPairStepTrajectoryFilterInOut')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPixelPairStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(False)
)

process.HLTPSetPixelPairStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetPixelPairStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(2),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetPixelPairStepTrajectoryFilterInOut = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(1),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetPvClusterComparer = cms.PSet(
    track_chi2_max = cms.double(9999999.0),
    track_prob_min = cms.double(-1.0),
    track_pt_max = cms.double(10.0),
    track_pt_min = cms.double(2.5)
)

process.HLTPSetPvClusterComparerForBTag = cms.PSet(
    track_chi2_max = cms.double(20.0),
    track_prob_min = cms.double(-1.0),
    track_pt_max = cms.double(20.0),
    track_pt_min = cms.double(0.1)
)

process.HLTPSetPvClusterComparerForIT = cms.PSet(
    track_chi2_max = cms.double(20.0),
    track_prob_min = cms.double(-1.0),
    track_pt_max = cms.double(20.0),
    track_pt_min = cms.double(1.0)
)

process.HLTPSetTobTecStepInOutTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetTobTecStepInOutTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetTobTecStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('hltESPTobTecStepChi2ChargeMeasurementEstimator16'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetTobTecStepInOutTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7),
    minNrOfHitsForRebuild = cms.int32(4),
    propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
    propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetTobTecStepTrajectoryFilter')
    ),
    updator = cms.string('hltESPKFUpdator'),
    useSameTrajFilter = cms.bool(False)
)

process.HLTPSetTobTecStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetTobTecStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetTrajectoryBuilderForElectrons = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    intermediateCleaning = cms.bool(False),
    lostHitPenalty = cms.double(90.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('hltESPFwdElectronPropagator'),
    propagatorOpposite = cms.string('hltESPBwdElectronPropagator'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetTrajectoryFilterForElectrons')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTPSetTrajectoryBuilderForGsfElectrons = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    estimator = cms.string('hltESPChi2ChargeMeasurementEstimator2000'),
    intermediateCleaning = cms.bool(False),
    lostHitPenalty = cms.double(90.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('hltESPFwdElectronPropagator'),
    propagatorOpposite = cms.string('hltESPBwdElectronPropagator'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('HLTPSetTrajectoryFilterForElectrons')
    ),
    updator = cms.string('hltESPKFUpdator')
)

process.HLTPSetTrajectoryFilterForElectrons = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(-1),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(-1),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(2.0),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetTrajectoryFilterIT = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.3),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetTrajectoryFilterL3 = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(1000000000),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.5),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTPSetbJetRegionalTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(999.0),
    maxNumberOfHits = cms.int32(8),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(1.0),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.HLTSeedFromConsecutiveHitsCreator = cms.PSet(
    ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string(''),
    propagator = cms.string('PropagatorWithMaterial')
)

process.HLTSeedFromConsecutiveHitsCreatorIT = cms.PSet(
    ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf')
)

process.HLTSeedFromConsecutiveHitsTripletOnlyCreator = cms.PSet(
    ComponentName = cms.string('SeedFromConsecutiveHitsTripletOnlyCreator'),
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf')
)

process.HLTSeedFromProtoTracks = cms.PSet(
    ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf')
)

process.HLTSiStripClusterChargeCutForHI = cms.PSet(
    value = cms.double(2069.0)
)

process.HLTSiStripClusterChargeCutLoose = cms.PSet(
    value = cms.double(1620.0)
)

process.HLTSiStripClusterChargeCutNone = cms.PSet(
    value = cms.double(-1.0)
)

process.HLTSiStripClusterChargeCutTight = cms.PSet(
    value = cms.double(1945.0)
)

process.HLTSiStripClusterChargeCutTiny = cms.PSet(
    value = cms.double(800.0)
)

process.datasets = cms.PSet(
    AlCaLumiPixels = cms.vstring(
        'AlCa_LumiPixels_Random_v4', 
        'AlCa_LumiPixels_ZeroBias_v8'
    ),
    AlCaP0 = cms.vstring(
        'AlCa_EcalEtaEBonly_v12', 
        'AlCa_EcalEtaEEonly_v12', 
        'AlCa_EcalPi0EBonly_v12', 
        'AlCa_EcalPi0EEonly_v12'
    ),
    AlCaPhiSym = cms.vstring('AlCa_EcalPhiSym_v8'),
    BTagCSV = cms.vstring(
        'HLT_AK8PFJet330_PFAK8BTagCSV_p17_v3', 
        'HLT_AK8PFJet330_PFAK8BTagCSV_p1_v3', 
        'HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets100_CaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets200_CaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets350_CaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets40_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0_v7', 
        'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_v7', 
        'HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_SingleJet30_Mu12_SinglePFJet40_v9'
    ),
    BTagMu = cms.vstring(
        'HLT_BTagMu_AK4DiJet110_Mu5_v11', 
        'HLT_BTagMu_AK4DiJet170_Mu5_v10', 
        'HLT_BTagMu_AK4DiJet20_Mu5_v11', 
        'HLT_BTagMu_AK4DiJet40_Mu5_v11', 
        'HLT_BTagMu_AK4DiJet70_Mu5_v11', 
        'HLT_BTagMu_AK4Jet300_Mu5_v11', 
        'HLT_BTagMu_AK8DiJet170_Mu5_v7', 
        'HLT_BTagMu_AK8Jet300_Mu5_v11'
    ),
    Charmonium = cms.vstring(
        'HLT_Dimuon0_Jpsi3p5_Muon2_v5', 
        'HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v6', 
        'HLT_Dimuon0_Jpsi_L1_NoOS_v6', 
        'HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v6', 
        'HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v6', 
        'HLT_Dimuon0_Jpsi_NoVertexing_v7', 
        'HLT_Dimuon0_Jpsi_v7', 
        'HLT_Dimuon0_LowMass_L1_0er1p5R_v6', 
        'HLT_Dimuon0_LowMass_L1_0er1p5_v7', 
        'HLT_Dimuon0_LowMass_L1_4R_v6', 
        'HLT_Dimuon0_LowMass_L1_4_v7', 
        'HLT_Dimuon0_LowMass_v7', 
        'HLT_Dimuon10_PsiPrime_Barrel_Seagulls_v6', 
        'HLT_Dimuon18_PsiPrime_noCorrL1_v4', 
        'HLT_Dimuon18_PsiPrime_v13', 
        'HLT_Dimuon20_Jpsi_Barrel_Seagulls_v6', 
        'HLT_Dimuon25_Jpsi_noCorrL1_v4', 
        'HLT_Dimuon25_Jpsi_v13', 
        'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v3', 
        'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi_v4', 
        'HLT_DoubleMu4_3_Bs_v13', 
        'HLT_DoubleMu4_3_Jpsi_Displaced_v14', 
        'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v6', 
        'HLT_DoubleMu4_JpsiTrk_Displaced_v14', 
        'HLT_DoubleMu4_Jpsi_Displaced_v6', 
        'HLT_DoubleMu4_Jpsi_NoVertexing_v6', 
        'HLT_DoubleMu4_PsiPrimeTrk_Displaced_v14', 
        'HLT_Mu7p5_L2Mu2_Jpsi_v9', 
        'HLT_Mu7p5_Track2_Jpsi_v10', 
        'HLT_Mu7p5_Track3p5_Jpsi_v10', 
        'HLT_Mu7p5_Track7_Jpsi_v10'
    ),
    Commissioning = cms.vstring(
        'HLT_IsoTrackHB_v3', 
        'HLT_IsoTrackHE_v3', 
        'HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_v2'
    ),
    DQMOnlineBeamspot = cms.vstring(
        'HLT_HT300_Beamspot_v9', 
        'HLT_HT450_Beamspot_v9', 
        'HLT_ZeroBias_Beamspot_v2'
    ),
    DisplacedJet = cms.vstring(
        'HLT_HT400_DisplacedDijet40_DisplacedTrack_v11', 
        'HLT_HT425_v8', 
        'HLT_HT430_DisplacedDijet40_DisplacedTrack_v11', 
        'HLT_HT430_DisplacedDijet60_DisplacedTrack_v11', 
        'HLT_HT430_DisplacedDijet80_DisplacedTrack_v11', 
        'HLT_HT550_DisplacedDijet60_Inclusive_v11', 
        'HLT_HT550_DisplacedDijet80_Inclusive_v9', 
        'HLT_HT650_DisplacedDijet60_Inclusive_v11', 
        'HLT_HT650_DisplacedDijet80_Inclusive_v12', 
        'HLT_HT750_DisplacedDijet80_Inclusive_v12'
    ),
    DoubleEG = cms.vstring(
        'HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v3', 
        'HLT_DiSC30_18_EIso_AND_HE_Mass70_v12', 
        'HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12', 
        'HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13', 
        'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12', 
        'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13', 
        'HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v8', 
        'HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v8', 
        'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v12', 
        'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v12', 
        'HLT_DoubleEle25_CaloIdL_MW_v2', 
        'HLT_DoubleEle27_CaloIdL_MW_v2', 
        'HLT_DoubleEle33_CaloIdL_MW_v15', 
        'HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v18', 
        'HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v18', 
        'HLT_DoublePhoton33_CaloIdL_v5', 
        'HLT_DoublePhoton70_v5', 
        'HLT_DoublePhoton85_v13', 
        'HLT_ECALHT800_v9', 
        'HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v1', 
        'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v7', 
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17', 
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17', 
        'HLT_Ele27_Ele37_CaloIdL_MW_v2', 
        'HLT_Ele28_HighEta_SC20_Mass55_v11', 
        'HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v2', 
        'HLT_TriplePhoton_20_20_20_CaloIdLV2_v2', 
        'HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_v3', 
        'HLT_TriplePhoton_30_30_10_CaloIdLV2_v3', 
        'HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v3'
    ),
    DoubleMuon = cms.vstring(
        'HLT_DoubleL2Mu50_v2', 
        'HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v8', 
        'HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v8', 
        'HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v8', 
        'HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v8', 
        'HLT_DoubleMu43NoFiltersNoVtx_v3', 
        'HLT_DoubleMu48NoFiltersNoVtx_v3', 
        'HLT_DoubleMu4_Mass8_DZ_PFHT350_v7', 
        'HLT_DoubleMu8_Mass8_PFHT350_v7', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v13', 
        'HLT_Mu17_TrkIsoVVL_v11', 
        'HLT_Mu17_v11', 
        'HLT_Mu18_Mu9_DZ_v3', 
        'HLT_Mu18_Mu9_SameSign_DZ_v3', 
        'HLT_Mu18_Mu9_SameSign_v3', 
        'HLT_Mu18_Mu9_v3', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v2', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v2', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v2', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v2', 
        'HLT_Mu19_TrkIsoVVL_v2', 
        'HLT_Mu19_v2', 
        'HLT_Mu20_Mu10_DZ_v3', 
        'HLT_Mu20_Mu10_SameSign_DZ_v3', 
        'HLT_Mu20_Mu10_SameSign_v3', 
        'HLT_Mu20_Mu10_v3', 
        'HLT_Mu23_Mu12_DZ_v3', 
        'HLT_Mu23_Mu12_SameSign_DZ_v3', 
        'HLT_Mu23_Mu12_SameSign_v3', 
        'HLT_Mu23_Mu12_v3', 
        'HLT_Mu37_TkMu27_v3', 
        'HLT_Mu8_TrkIsoVVL_v11', 
        'HLT_Mu8_v11', 
        'HLT_TripleMu_10_5_5_DZ_v9', 
        'HLT_TripleMu_12_10_5_v9', 
        'HLT_TripleMu_5_3_3_Mass3p8to60_DCA_v2', 
        'HLT_TripleMu_5_3_3_Mass3p8to60_DZ_v7', 
        'HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_v5', 
        'HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_v11', 
        'HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_v12'
    ),
    DoubleMuonLowMass = cms.vstring(
        'HLT_Dimuon0_LowMass_L1_TM530_v5', 
        'HLT_DoubleMu3_TkMu_DsTau3Mu_v2', 
        'HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v5', 
        'HLT_DoubleMu3_Trk_Tau3mu_v11', 
        'HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v14', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v2', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_v2', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_v2', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_v2'
    ),
    EcalLaser = cms.vstring('HLT_EcalCalibration_v4'),
    EmptyBX = cms.vstring(
        'HLT_L1NotBptxOR_v3', 
        'HLT_L1UnpairedBunchBptxMinus_v2', 
        'HLT_L1UnpairedBunchBptxPlus_v2'
    ),
    EphemeralHLTPhysics1 = cms.vstring('HLT_Physics_part0_v7'),
    EphemeralHLTPhysics2 = cms.vstring('HLT_Physics_part1_v7'),
    EphemeralHLTPhysics3 = cms.vstring('HLT_Physics_part2_v7'),
    EphemeralHLTPhysics4 = cms.vstring('HLT_Physics_part3_v7'),
    EphemeralHLTPhysics5 = cms.vstring('HLT_Physics_part4_v7'),
    EphemeralHLTPhysics6 = cms.vstring('HLT_Physics_part5_v7'),
    EphemeralHLTPhysics7 = cms.vstring('HLT_Physics_part6_v7'),
    EphemeralHLTPhysics8 = cms.vstring('HLT_Physics_part7_v7'),
    EphemeralZeroBias1 = cms.vstring('HLT_ZeroBias_part0_v6'),
    EphemeralZeroBias2 = cms.vstring('HLT_ZeroBias_part1_v6'),
    EphemeralZeroBias3 = cms.vstring('HLT_ZeroBias_part2_v6'),
    EphemeralZeroBias4 = cms.vstring('HLT_ZeroBias_part3_v6'),
    EphemeralZeroBias5 = cms.vstring('HLT_ZeroBias_part4_v6'),
    EphemeralZeroBias6 = cms.vstring('HLT_ZeroBias_part5_v6'),
    EphemeralZeroBias7 = cms.vstring('HLT_ZeroBias_part6_v6'),
    EphemeralZeroBias8 = cms.vstring('HLT_ZeroBias_part7_v6'),
    EventDisplay = cms.vstring(
        'HLT_AK4PFJet100_v17', 
        'HLT_DoublePhoton85_v13', 
        'HLT_PFJet500_v19'
    ),
    ExpressAlignment = cms.vstring(
        'HLT_HT300_Beamspot_v9', 
        'HLT_HT450_Beamspot_v9', 
        'HLT_ZeroBias_Beamspot_v2'
    ),
    ExpressPhysics = cms.vstring(
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17', 
        'HLT_IsoMu20_v13', 
        'HLT_IsoMu24_v11', 
        'HLT_IsoMu27_v14', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14', 
        'HLT_Physics_v7', 
        'HLT_Random_v3', 
        'HLT_ZeroBias_FirstCollisionAfterAbortGap_v5', 
        'HLT_ZeroBias_IsolatedBunches_v5', 
        'HLT_ZeroBias_v6'
    ),
    FSQJet1 = cms.vstring(
        'HLT_DiPFJet15_NoCaloMatched_v14', 
        'HLT_DiPFJet25_NoCaloMatched_v14'
    ),
    FSQJet2 = cms.vstring(
        'HLT_DiPFJet15_FBEta3_NoCaloMatched_v15', 
        'HLT_DiPFJet25_FBEta3_NoCaloMatched_v15', 
        'HLT_DiPFJetAve15_HFJEC_v15', 
        'HLT_DiPFJetAve25_HFJEC_v15', 
        'HLT_DiPFJetAve35_HFJEC_v15'
    ),
    HINCaloJets = cms.vstring(
        'HLT_AK4CaloJet100_v10', 
        'HLT_AK4CaloJet120_v9', 
        'HLT_AK4CaloJet30_v11', 
        'HLT_AK4CaloJet40_v10', 
        'HLT_AK4CaloJet50_v10', 
        'HLT_AK4CaloJet80_v10'
    ),
    HINPFJets = cms.vstring(
        'HLT_AK4PFJet100_v17', 
        'HLT_AK4PFJet120_v16', 
        'HLT_AK4PFJet30_v17', 
        'HLT_AK4PFJet50_v17', 
        'HLT_AK4PFJet80_v17'
    ),
    HLTMonitor = cms.vstring(
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17', 
        'HLT_Ele32_WPTight_Gsf_v13', 
        'HLT_HT400_DisplacedDijet40_DisplacedTrack_v11', 
        'HLT_HT550_DisplacedDijet60_Inclusive_v11', 
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13', 
        'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0_v7', 
        'HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v7', 
        'HLT_PFHT510_v15', 
        'HLT_PFJet260_v18', 
        'HLT_PFJet320_v18', 
        'HLT_PFMET130_PFMHT130_IDTight_v18', 
        'HLT_PFMET140_PFMHT140_IDTight_v18'
    ),
    HLTPhysics = cms.vstring('HLT_Physics_v7'),
    HTMHT = cms.vstring(
        'HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10', 
        'HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10', 
        'HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10', 
        'HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10', 
        'HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10', 
        'HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10', 
        'HLT_Rsq0p35_v13', 
        'HLT_Rsq0p40_v13', 
        'HLT_RsqMR300_Rsq0p09_MR200_4jet_v13', 
        'HLT_RsqMR300_Rsq0p09_MR200_v13', 
        'HLT_RsqMR320_Rsq0p09_MR200_4jet_v13', 
        'HLT_RsqMR320_Rsq0p09_MR200_v13'
    ),
    HcalNZS = cms.vstring(
        'HLT_HcalNZS_v12', 
        'HLT_HcalPhiSym_v14'
    ),
    HighMultiplicityEOF = cms.vstring(
        'HLT_FullTrack_Multiplicity100_v2', 
        'HLT_FullTrack_Multiplicity130_v2', 
        'HLT_FullTrack_Multiplicity155_v3', 
        'HLT_FullTrack_Multiplicity85_v3'
    ),
    HighPtLowerPhotons = cms.vstring(
        'HLT_SinglePhoton10_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton20_Eta3p1ForPPRef_v8'
    ),
    HighPtPhoton30AndZ = cms.vstring(
        'HLT_SinglePhoton30_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton40_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton50_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton60_Eta3p1ForPPRef_v8'
    ),
    IsolatedBunch = cms.vstring('HLT_HcalIsolatedbunch_v4'),
    JetHT = cms.vstring(
        'HLT_AK8PFHT750_TrimMass50_v10', 
        'HLT_AK8PFHT800_TrimMass50_v10', 
        'HLT_AK8PFHT850_TrimMass50_v9', 
        'HLT_AK8PFHT900_TrimMass50_v9', 
        'HLT_AK8PFJet140_v13', 
        'HLT_AK8PFJet200_v13', 
        'HLT_AK8PFJet260_v14', 
        'HLT_AK8PFJet320_v14', 
        'HLT_AK8PFJet360_TrimMass30_v16', 
        'HLT_AK8PFJet380_TrimMass30_v9', 
        'HLT_AK8PFJet400_TrimMass30_v10', 
        'HLT_AK8PFJet400_v14', 
        'HLT_AK8PFJet40_v14', 
        'HLT_AK8PFJet420_TrimMass30_v9', 
        'HLT_AK8PFJet450_v14', 
        'HLT_AK8PFJet500_v14', 
        'HLT_AK8PFJet550_v9', 
        'HLT_AK8PFJet60_v13', 
        'HLT_AK8PFJet80_v13', 
        'HLT_AK8PFJetFwd140_v12', 
        'HLT_AK8PFJetFwd200_v12', 
        'HLT_AK8PFJetFwd260_v13', 
        'HLT_AK8PFJetFwd320_v13', 
        'HLT_AK8PFJetFwd400_v13', 
        'HLT_AK8PFJetFwd40_v13', 
        'HLT_AK8PFJetFwd450_v13', 
        'HLT_AK8PFJetFwd500_v13', 
        'HLT_AK8PFJetFwd60_v12', 
        'HLT_AK8PFJetFwd80_v12', 
        'HLT_CaloJet500_NoJetID_v11', 
        'HLT_CaloJet550_NoJetID_v6', 
        'HLT_DiPFJetAve100_HFJEC_v13', 
        'HLT_DiPFJetAve140_v11', 
        'HLT_DiPFJetAve160_HFJEC_v13', 
        'HLT_DiPFJetAve200_v11', 
        'HLT_DiPFJetAve220_HFJEC_v14', 
        'HLT_DiPFJetAve260_v12', 
        'HLT_DiPFJetAve300_HFJEC_v14', 
        'HLT_DiPFJetAve320_v12', 
        'HLT_DiPFJetAve400_v12', 
        'HLT_DiPFJetAve40_v12', 
        'HLT_DiPFJetAve500_v12', 
        'HLT_DiPFJetAve60_HFJEC_v13', 
        'HLT_DiPFJetAve60_v12', 
        'HLT_DiPFJetAve80_HFJEC_v13', 
        'HLT_DiPFJetAve80_v11', 
        'HLT_PFHT1050_v16', 
        'HLT_PFHT180_v15', 
        'HLT_PFHT250_v15', 
        'HLT_PFHT350MinPFJet15_v7', 
        'HLT_PFHT350_v17', 
        'HLT_PFHT370_v15', 
        'HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v7', 
        'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v6', 
        'HLT_PFHT380_SixPFJet32_v7', 
        'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v7', 
        'HLT_PFHT430_SixPFJet40_v9', 
        'HLT_PFHT430_v15', 
        'HLT_PFHT510_v15', 
        'HLT_PFHT590_v15', 
        'HLT_PFHT680_v15', 
        'HLT_PFHT780_v15', 
        'HLT_PFHT890_v15', 
        'HLT_PFJet140_v17', 
        'HLT_PFJet200_v17', 
        'HLT_PFJet260_v18', 
        'HLT_PFJet320_v18', 
        'HLT_PFJet400_v18', 
        'HLT_PFJet40_v19', 
        'HLT_PFJet450_v19', 
        'HLT_PFJet500_v19', 
        'HLT_PFJet550_v9', 
        'HLT_PFJet60_v19', 
        'HLT_PFJet80_v18', 
        'HLT_PFJetFwd140_v16', 
        'HLT_PFJetFwd200_v16', 
        'HLT_PFJetFwd260_v17', 
        'HLT_PFJetFwd320_v17', 
        'HLT_PFJetFwd400_v17', 
        'HLT_PFJetFwd40_v17', 
        'HLT_PFJetFwd450_v17', 
        'HLT_PFJetFwd500_v17', 
        'HLT_PFJetFwd60_v17', 
        'HLT_PFJetFwd80_v16', 
        'HLT_QuadPFJet103_88_75_15_v3', 
        'HLT_QuadPFJet105_88_76_15_v3', 
        'HLT_QuadPFJet111_90_80_15_v3', 
        'HLT_QuadPFJet98_83_71_15_v3'
    ),
    L1Accept = cms.vstring(
        'DST_Physics_v7', 
        'DST_ZeroBias_v2'
    ),
    L1MinimumBias = cms.vstring('HLT_L1MinimumBiasHF_OR_v2'),
    MET = cms.vstring(
        'HLT_CaloMET100_HBHECleaned_v3', 
        'HLT_CaloMET100_NotCleaned_v3', 
        'HLT_CaloMET110_NotCleaned_v3', 
        'HLT_CaloMET250_HBHECleaned_v3', 
        'HLT_CaloMET250_NotCleaned_v3', 
        'HLT_CaloMET300_HBHECleaned_v3', 
        'HLT_CaloMET350_HBHECleaned_v3', 
        'HLT_CaloMET70_HBHECleaned_v3', 
        'HLT_CaloMET80_HBHECleaned_v3', 
        'HLT_CaloMET80_NotCleaned_v3', 
        'HLT_CaloMET90_HBHECleaned_v3', 
        'HLT_CaloMET90_NotCleaned_v3', 
        'HLT_CaloMHT90_v3', 
        'HLT_DiJet110_35_Mjj650_PFMET110_v7', 
        'HLT_DiJet110_35_Mjj650_PFMET120_v7', 
        'HLT_DiJet110_35_Mjj650_PFMET130_v7', 
        'HLT_L1ETMHadSeeds_v1', 
        'HLT_MET105_IsoTrk50_v8', 
        'HLT_MET120_IsoTrk50_v8', 
        'HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v18', 
        'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v18', 
        'HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_v17', 
        'HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_v17', 
        'HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET100_PFMHT100_IDTight_PFHT60_v7', 
        'HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET110_PFMHT110_IDTight_v18', 
        'HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7', 
        'HLT_PFMET120_PFMHT120_IDTight_v18', 
        'HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET130_PFMHT130_IDTight_v18', 
        'HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET140_PFMHT140_IDTight_v18', 
        'HLT_PFMET200_HBHECleaned_v7', 
        'HLT_PFMET200_HBHE_BeamHaloCleaned_v7', 
        'HLT_PFMET200_NotCleaned_v7', 
        'HLT_PFMET250_HBHECleaned_v7', 
        'HLT_PFMET300_HBHECleaned_v7', 
        'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v7', 
        'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v18', 
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7', 
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v18', 
        'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v17', 
        'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17', 
        'HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v7', 
        'HLT_PFMETTypeOne110_PFMHT110_IDTight_v10', 
        'HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v7', 
        'HLT_PFMETTypeOne120_PFMHT120_IDTight_v10', 
        'HLT_PFMETTypeOne130_PFMHT130_IDTight_v10', 
        'HLT_PFMETTypeOne140_PFMHT140_IDTight_v9', 
        'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v7', 
        'HLT_TripleJet110_35_35_Mjj650_PFMET110_v7', 
        'HLT_TripleJet110_35_35_Mjj650_PFMET120_v7', 
        'HLT_TripleJet110_35_35_Mjj650_PFMET130_v7'
    ),
    MonteCarlo = cms.vstring(
        'MC_AK4CaloJetsFromPV_v7', 
        'MC_AK4CaloJets_v9', 
        'MC_AK4PFJets_v15', 
        'MC_AK8CaloHT_v8', 
        'MC_AK8PFHT_v14', 
        'MC_AK8PFJets_v15', 
        'MC_AK8TrimPFJets_v15', 
        'MC_CaloBTagCSV_v7', 
        'MC_CaloHT_v8', 
        'MC_CaloMET_JetIdCleaned_v9', 
        'MC_CaloMET_v8', 
        'MC_CaloMHT_v8', 
        'MC_Diphoton10_10_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass10_v12', 
        'MC_DoubleEle5_CaloIdL_MW_v13', 
        'MC_DoubleMuNoFiltersNoVtx_v7', 
        'MC_DoubleMu_TrkIsoVVL_DZ_v10', 
        'MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_v13', 
        'MC_Ele5_WPTight_Gsf_v6', 
        'MC_IsoMu_v13', 
        'MC_PFBTagCSV_v9', 
        'MC_PFHT_v14', 
        'MC_PFMET_v15', 
        'MC_PFMHT_v14', 
        'MC_ReducedIterativeTracking_v11'
    ),
    MuOnia = cms.vstring(
        'HLT_Dimuon0_Upsilon_L1_4p5NoOS_v6', 
        'HLT_Dimuon0_Upsilon_L1_4p5_v7', 
        'HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v6', 
        'HLT_Dimuon0_Upsilon_L1_4p5er2p0_v7', 
        'HLT_Dimuon0_Upsilon_L1_5M_v6', 
        'HLT_Dimuon0_Upsilon_L1_5_v7', 
        'HLT_Dimuon0_Upsilon_Muon_L1_TM0_v5', 
        'HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v5', 
        'HLT_Dimuon0_Upsilon_NoVertexing_v6', 
        'HLT_Dimuon10_Upsilon_Barrel_Seagulls_v6', 
        'HLT_Dimuon12_Upsilon_eta1p5_v13', 
        'HLT_Dimuon14_Phi_Barrel_Seagulls_v6', 
        'HLT_Dimuon24_Phi_noCorrL1_v4', 
        'HLT_Dimuon24_Upsilon_noCorrL1_v4', 
        'HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v2', 
        'HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v2', 
        'HLT_Mu20_TkMu0_Phi_v6', 
        'HLT_Mu25_TkMu0_Onia_v6', 
        'HLT_Mu25_TkMu0_Phi_v6', 
        'HLT_Mu30_TkMu0_Onia_v6', 
        'HLT_Mu7p5_L2Mu2_Upsilon_v9', 
        'HLT_Mu7p5_Track2_Upsilon_v10', 
        'HLT_Mu7p5_Track3p5_Upsilon_v10', 
        'HLT_Mu7p5_Track7_Upsilon_v10', 
        'HLT_Trimuon5_3p5_2_Upsilon_Muon_v4', 
        'HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v2'
    ),
    MuonEG = cms.vstring(
        'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v15', 
        'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v15', 
        'HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v6', 
        'HLT_DoubleMu20_7_Mass0to30_L1_DM4_v6', 
        'HLT_DoubleMu20_7_Mass0to30_Photon23_v6', 
        'HLT_Mu12_DoublePhoton20_v3', 
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13', 
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5', 
        'HLT_Mu17_Photon30_IsoCaloId_v4', 
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13', 
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5', 
        'HLT_Mu27_Ele37_CaloIdL_MW_v3', 
        'HLT_Mu37_Ele27_CaloIdL_MW_v3', 
        'HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v4', 
        'HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v4', 
        'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v16', 
        'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v16', 
        'HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v17', 
        'HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v17', 
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11', 
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9'
    ),
    NoBPTX = cms.vstring(
        'HLT_L2Mu10_NoVertex_NoBPTX3BX_v5', 
        'HLT_L2Mu10_NoVertex_NoBPTX_v6', 
        'HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_v5', 
        'HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_v4', 
        'HLT_UncorrectedJetE30_NoBPTX3BX_v5', 
        'HLT_UncorrectedJetE30_NoBPTX_v5', 
        'HLT_UncorrectedJetE60_NoBPTX3BX_v5', 
        'HLT_UncorrectedJetE70_NoBPTX3BX_v5'
    ),
    OnlineMonitor = cms.vstring( (
        'HLT_AK4CaloJet100_v10', 
        'HLT_AK4CaloJet120_v9', 
        'HLT_AK4CaloJet30_v11', 
        'HLT_AK4CaloJet40_v10', 
        'HLT_AK4CaloJet50_v10', 
        'HLT_AK4CaloJet80_v10', 
        'HLT_AK4PFJet100_v17', 
        'HLT_AK4PFJet120_v16', 
        'HLT_AK4PFJet30_v17', 
        'HLT_AK4PFJet50_v17', 
        'HLT_AK4PFJet80_v17', 
        'HLT_AK8PFHT750_TrimMass50_v10', 
        'HLT_AK8PFHT800_TrimMass50_v10', 
        'HLT_AK8PFHT850_TrimMass50_v9', 
        'HLT_AK8PFHT900_TrimMass50_v9', 
        'HLT_AK8PFJet140_v13', 
        'HLT_AK8PFJet200_v13', 
        'HLT_AK8PFJet260_v14', 
        'HLT_AK8PFJet320_v14', 
        'HLT_AK8PFJet330_PFAK8BTagCSV_p17_v3', 
        'HLT_AK8PFJet330_PFAK8BTagCSV_p1_v3', 
        'HLT_AK8PFJet360_TrimMass30_v16', 
        'HLT_AK8PFJet380_TrimMass30_v9', 
        'HLT_AK8PFJet400_TrimMass30_v10', 
        'HLT_AK8PFJet400_v14', 
        'HLT_AK8PFJet40_v14', 
        'HLT_AK8PFJet420_TrimMass30_v9', 
        'HLT_AK8PFJet450_v14', 
        'HLT_AK8PFJet500_v14', 
        'HLT_AK8PFJet550_v9', 
        'HLT_AK8PFJet60_v13', 
        'HLT_AK8PFJet80_v13', 
        'HLT_AK8PFJetFwd140_v12', 
        'HLT_AK8PFJetFwd200_v12', 
        'HLT_AK8PFJetFwd260_v13', 
        'HLT_AK8PFJetFwd320_v13', 
        'HLT_AK8PFJetFwd400_v13', 
        'HLT_AK8PFJetFwd40_v13', 
        'HLT_AK8PFJetFwd450_v13', 
        'HLT_AK8PFJetFwd500_v13', 
        'HLT_AK8PFJetFwd60_v12', 
        'HLT_AK8PFJetFwd80_v12', 
        'HLT_BTagMu_AK4DiJet110_Mu5_v11', 
        'HLT_BTagMu_AK4DiJet170_Mu5_v10', 
        'HLT_BTagMu_AK4DiJet20_Mu5_v11', 
        'HLT_BTagMu_AK4DiJet40_Mu5_v11', 
        'HLT_BTagMu_AK4DiJet70_Mu5_v11', 
        'HLT_BTagMu_AK4Jet300_Mu5_v11', 
        'HLT_BTagMu_AK8DiJet170_Mu5_v7', 
        'HLT_BTagMu_AK8Jet300_Mu5_v11', 
        'HLT_CaloJet500_NoJetID_v11', 
        'HLT_CaloJet550_NoJetID_v6', 
        'HLT_CaloMET100_HBHECleaned_v3', 
        'HLT_CaloMET100_NotCleaned_v3', 
        'HLT_CaloMET110_NotCleaned_v3', 
        'HLT_CaloMET250_HBHECleaned_v3', 
        'HLT_CaloMET250_NotCleaned_v3', 
        'HLT_CaloMET300_HBHECleaned_v3', 
        'HLT_CaloMET350_HBHECleaned_v3', 
        'HLT_CaloMET70_HBHECleaned_v3', 
        'HLT_CaloMET80_HBHECleaned_v3', 
        'HLT_CaloMET80_NotCleaned_v3', 
        'HLT_CaloMET90_HBHECleaned_v3', 
        'HLT_CaloMET90_NotCleaned_v3', 
        'HLT_CaloMHT90_v3', 
        'HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v3', 
        'HLT_DiJet110_35_Mjj650_PFMET110_v7', 
        'HLT_DiJet110_35_Mjj650_PFMET120_v7', 
        'HLT_DiJet110_35_Mjj650_PFMET130_v7', 
        'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v15', 
        'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v15', 
        'HLT_DiPFJet15_FBEta3_NoCaloMatched_v15', 
        'HLT_DiPFJet15_NoCaloMatched_v14', 
        'HLT_DiPFJet25_FBEta3_NoCaloMatched_v15', 
        'HLT_DiPFJet25_NoCaloMatched_v14', 
        'HLT_DiPFJetAve100_HFJEC_v13', 
        'HLT_DiPFJetAve140_v11', 
        'HLT_DiPFJetAve15_HFJEC_v15', 
        'HLT_DiPFJetAve160_HFJEC_v13', 
        'HLT_DiPFJetAve200_v11', 
        'HLT_DiPFJetAve220_HFJEC_v14', 
        'HLT_DiPFJetAve25_HFJEC_v15', 
        'HLT_DiPFJetAve260_v12', 
        'HLT_DiPFJetAve300_HFJEC_v14', 
        'HLT_DiPFJetAve320_v12', 
        'HLT_DiPFJetAve35_HFJEC_v15', 
        'HLT_DiPFJetAve400_v12', 
        'HLT_DiPFJetAve40_v12', 
        'HLT_DiPFJetAve500_v12', 
        'HLT_DiPFJetAve60_HFJEC_v13', 
        'HLT_DiPFJetAve60_v12', 
        'HLT_DiPFJetAve80_HFJEC_v13', 
        'HLT_DiPFJetAve80_v11', 
        'HLT_DiSC30_18_EIso_AND_HE_Mass70_v12', 
        'HLT_Dimuon0_Jpsi3p5_Muon2_v5', 
        'HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v6', 
        'HLT_Dimuon0_Jpsi_L1_NoOS_v6', 
        'HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v6', 
        'HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v6', 
        'HLT_Dimuon0_Jpsi_NoVertexing_v7', 
        'HLT_Dimuon0_Jpsi_v7', 
        'HLT_Dimuon0_LowMass_L1_0er1p5R_v6', 
        'HLT_Dimuon0_LowMass_L1_0er1p5_v7', 
        'HLT_Dimuon0_LowMass_L1_4R_v6', 
        'HLT_Dimuon0_LowMass_L1_4_v7', 
        'HLT_Dimuon0_LowMass_L1_TM530_v5', 
        'HLT_Dimuon0_LowMass_v7', 
        'HLT_Dimuon0_Upsilon_L1_4p5NoOS_v6', 
        'HLT_Dimuon0_Upsilon_L1_4p5_v7', 
        'HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v6', 
        'HLT_Dimuon0_Upsilon_L1_4p5er2p0_v7', 
        'HLT_Dimuon0_Upsilon_L1_5M_v6', 
        'HLT_Dimuon0_Upsilon_L1_5_v7', 
        'HLT_Dimuon0_Upsilon_Muon_L1_TM0_v5', 
        'HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v5', 
        'HLT_Dimuon0_Upsilon_NoVertexing_v6', 
        'HLT_Dimuon10_PsiPrime_Barrel_Seagulls_v6', 
        'HLT_Dimuon10_Upsilon_Barrel_Seagulls_v6', 
        'HLT_Dimuon12_Upsilon_eta1p5_v13', 
        'HLT_Dimuon14_Phi_Barrel_Seagulls_v6', 
        'HLT_Dimuon18_PsiPrime_noCorrL1_v4', 
        'HLT_Dimuon18_PsiPrime_v13', 
        'HLT_Dimuon20_Jpsi_Barrel_Seagulls_v6', 
        'HLT_Dimuon24_Phi_noCorrL1_v4', 
        'HLT_Dimuon24_Upsilon_noCorrL1_v4', 
        'HLT_Dimuon25_Jpsi_noCorrL1_v4', 
        'HLT_Dimuon25_Jpsi_v13', 
        'HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12', 
        'HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13', 
        'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12', 
        'HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13', 
        'HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v8', 
        'HLT_Diphoton30_18_PVrealAND_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v8', 
        'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v12', 
        'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v12', 
        'HLT_DoubleEle25_CaloIdL_MW_v2', 
        'HLT_DoubleEle27_CaloIdL_MW_v2', 
        'HLT_DoubleEle33_CaloIdL_MW_v15', 
        'HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v18', 
        'HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v18', 
        'HLT_DoubleL2Mu50_v2', 
        'HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v6', 
        'HLT_DoubleMu20_7_Mass0to30_L1_DM4_v6', 
        'HLT_DoubleMu20_7_Mass0to30_Photon23_v6', 
        'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v3', 
        'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi_v4', 
        'HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v8', 
        'HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v8', 
        'HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v8', 
        'HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v8', 
        'HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v2', 
        'HLT_DoubleMu3_TkMu_DsTau3Mu_v2', 
        'HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v5', 
        'HLT_DoubleMu3_Trk_Tau3mu_v11', 
        'HLT_DoubleMu43NoFiltersNoVtx_v3', 
        'HLT_DoubleMu48NoFiltersNoVtx_v3', 
        'HLT_DoubleMu4_3_Bs_v13', 
        'HLT_DoubleMu4_3_Jpsi_Displaced_v14', 
        'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v6', 
        'HLT_DoubleMu4_JpsiTrk_Displaced_v14', 
        'HLT_DoubleMu4_Jpsi_Displaced_v6', 
        'HLT_DoubleMu4_Jpsi_NoVertexing_v6', 
        'HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v14', 
        'HLT_DoubleMu4_Mass8_DZ_PFHT350_v7', 
        'HLT_DoubleMu4_PsiPrimeTrk_Displaced_v14', 
        'HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v2', 
        'HLT_DoubleMu8_Mass8_PFHT350_v7', 
        'HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets100_CaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets200_CaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets350_CaloBTagCSV_p33_v7', 
        'HLT_DoublePFJets40_CaloBTagCSV_p33_v7', 
        'HLT_DoublePhoton33_CaloIdL_v5', 
        'HLT_DoublePhoton70_v5', 
        'HLT_DoublePhoton85_v13', 
        'HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v10', 
        'HLT_ECALHT800_v9', 
        'HLT_Ele115_CaloIdVT_GsfTrkIdT_v12', 
        'HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v16', 
        'HLT_Ele135_CaloIdVT_GsfTrkIdT_v5', 
        'HLT_Ele145_CaloIdVT_GsfTrkIdT_v6', 
        'HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v1', 
        'HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v1', 
        'HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v7', 
        'HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v14', 
        'HLT_Ele15_IsoVVVL_PFHT450_v14', 
        'HLT_Ele15_IsoVVVL_PFHT600_v18', 
        'HLT_Ele15_WPLoose_Gsf_v1', 
        'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v7', 
        'HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v14', 
        'HLT_Ele17_WPLoose_Gsf_v1', 
        'HLT_Ele200_CaloIdVT_GsfTrkIdT_v6', 
        'HLT_Ele20_WPLoose_Gsf_v4', 
        'HLT_Ele20_WPTight_Gsf_v4', 
        'HLT_Ele20_eta2p1_WPLoose_Gsf_v4', 
        'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v16', 
        'HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v16', 
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17', 
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v11', 
        'HLT_Ele250_CaloIdVT_GsfTrkIdT_v11', 
        'HLT_Ele27_Ele37_CaloIdL_MW_v2', 
        'HLT_Ele27_WPTight_Gsf_v14', 
        'HLT_Ele28_HighEta_SC20_Mass55_v11', 
        'HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v11', 
        'HLT_Ele300_CaloIdVT_GsfTrkIdT_v11', 
        'HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v11', 
        'HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7', 
        'HLT_Ele32_WPTight_Gsf_v13', 
        'HLT_Ele35_WPTight_Gsf_L1EGMT_v3', 
        'HLT_Ele35_WPTight_Gsf_v7', 
        'HLT_Ele38_WPTight_Gsf_v7', 
        'HLT_Ele40_WPTight_Gsf_v7', 
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v16', 
        'HLT_Ele50_IsoVVVL_PFHT450_v14', 
        'HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v14', 
        'HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v16', 
        'HLT_FullTrack_Multiplicity100_v2', 
        'HLT_FullTrack_Multiplicity130_v2', 
        'HLT_FullTrack_Multiplicity155_v3', 
        'HLT_FullTrack_Multiplicity85_v3', 
        'HLT_HT400_DisplacedDijet40_DisplacedTrack_v11', 
        'HLT_HT425_v8', 
        'HLT_HT430_DisplacedDijet40_DisplacedTrack_v11', 
        'HLT_HT430_DisplacedDijet60_DisplacedTrack_v11', 
        'HLT_HT430_DisplacedDijet80_DisplacedTrack_v11', 
        'HLT_HT550_DisplacedDijet60_Inclusive_v11', 
        'HLT_HT550_DisplacedDijet80_Inclusive_v9', 
        'HLT_HT650_DisplacedDijet60_Inclusive_v11', 
        'HLT_HT650_DisplacedDijet80_Inclusive_v12', 
        'HLT_HT750_DisplacedDijet80_Inclusive_v12', 
        'HLT_HcalNZS_v12', 
        'HLT_HcalPhiSym_v14', 
        'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v10', 
        'HLT_IsoMu20_v13', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_v13', 
        'HLT_IsoMu24_v11', 
        'HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1_v3', 
        'HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1_v3', 
        'HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1_v3', 
        'HLT_IsoMu27_v14', 
        'HLT_IsoMu30_v2', 
        'HLT_IsoTrackHB_v3', 
        'HLT_IsoTrackHE_v3', 
        'HLT_L1ETMHadSeeds_v1', 
        'HLT_L1MinimumBiasHF0OR_v3', 
        'HLT_L1MinimumBiasHF_OR_v2', 
        'HLT_L1NotBptxOR_v3', 
        'HLT_L1SingleMu18_v3', 
        'HLT_L1SingleMu25_v2', 
        'HLT_L1UnpairedBunchBptxMinus_v2', 
        'HLT_L1UnpairedBunchBptxPlus_v2', 
        'HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_v2', 
        'HLT_L1_DoubleJet30_Mass_Min400_Mu10_v1', 
        'HLT_L2Mu10_NoVertex_NoBPTX3BX_v5', 
        'HLT_L2Mu10_NoVertex_NoBPTX_v6', 
        'HLT_L2Mu10_v7', 
        'HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_v5', 
        'HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_v4', 
        'HLT_L2Mu50_v2', 
        'HLT_MET105_IsoTrk50_v8', 
        'HLT_MET120_IsoTrk50_v8', 
        'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v9', 
        'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v10', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v10', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v6', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_v6', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v6', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v10', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v10', 
        'HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v18', 
        'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v18', 
        'HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_v17', 
        'HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_v17', 
        'HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v12', 
        'HLT_Mu12_DoublePFJets100_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets200_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets350_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets40_CaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagCSV_p33_v7', 
        'HLT_Mu12_DoublePhoton20_v3', 
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13', 
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5', 
        'HLT_Mu12_v1', 
        'HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v7', 
        'HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v13', 
        'HLT_Mu15_IsoVVVL_PFHT450_v13', 
        'HLT_Mu15_IsoVVVL_PFHT600_v17', 
        'HLT_Mu15_v1', 
        'HLT_Mu17_Photon30_IsoCaloId_v4', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v13', 
        'HLT_Mu17_TrkIsoVVL_v11', 
        'HLT_Mu17_v11', 
        'HLT_Mu18_Mu9_DZ_v3', 
        'HLT_Mu18_Mu9_SameSign_DZ_v3', 
        'HLT_Mu18_Mu9_SameSign_v3', 
        'HLT_Mu18_Mu9_v3', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v2', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v2', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v2', 
        'HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v2', 
        'HLT_Mu19_TrkIsoVVL_v2', 
        'HLT_Mu19_v2', 
        'HLT_Mu20_Mu10_DZ_v3', 
        'HLT_Mu20_Mu10_SameSign_DZ_v3', 
        'HLT_Mu20_Mu10_SameSign_v3', 
        'HLT_Mu20_Mu10_v3', 
        'HLT_Mu20_TkMu0_Phi_v6', 
        'HLT_Mu20_v11', 
        'HLT_Mu23_Mu12_DZ_v3', 
        'HLT_Mu23_Mu12_SameSign_DZ_v3', 
        'HLT_Mu23_Mu12_SameSign_v3', 
        'HLT_Mu23_Mu12_v3', 
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13', 
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5', 
        'HLT_Mu25_TkMu0_Onia_v6', 
        'HLT_Mu25_TkMu0_Phi_v6', 
        'HLT_Mu27_Ele37_CaloIdL_MW_v3', 
        'HLT_Mu27_v12', 
        'HLT_Mu30_TkMu0_Onia_v6', 
        'HLT_Mu37_Ele27_CaloIdL_MW_v3', 
        'HLT_Mu37_TkMu27_v3', 
        'HLT_Mu3_PFJet40_v14', 
        'HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v4', 
        'HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v4', 
        'HLT_Mu50_IsoVVVL_PFHT450_v13', 
        'HLT_Mu50_v12', 
        'HLT_Mu55_v2', 
        'HLT_Mu7p5_L2Mu2_Jpsi_v9', 
        'HLT_Mu7p5_L2Mu2_Upsilon_v9', 
        'HLT_Mu7p5_Track2_Jpsi_v10', 
        'HLT_Mu7p5_Track2_Upsilon_v10', 
        'HLT_Mu7p5_Track3p5_Jpsi_v10', 
        'HLT_Mu7p5_Track3p5_Upsilon_v10', 
        'HLT_Mu7p5_Track7_Jpsi_v10', 
        'HLT_Mu7p5_Track7_Upsilon_v10', 
        'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v16', 
        'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v16', 
        'HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v17', 
        'HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v17', 
        'HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v13', 
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11', 
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9', 
        'HLT_Mu8_TrkIsoVVL_v11', 
        'HLT_Mu8_v11', 
        'HLT_OldMu100_v3', 
        'HLT_PFHT1050_v16', 
        'HLT_PFHT180_v15', 
        'HLT_PFHT250_v15', 
        'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0_v7', 
        'HLT_PFHT300PT30_QuadPFJet_75_60_45_40_v7', 
        'HLT_PFHT350MinPFJet15_v7', 
        'HLT_PFHT350_v17', 
        'HLT_PFHT370_v15', 
        'HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v7', 
        'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v6', 
        'HLT_PFHT380_SixPFJet32_v7', 
        'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v7', 
        'HLT_PFHT430_SixPFJet40_v9', 
        'HLT_PFHT430_v15', 
        'HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10', 
        'HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10', 
        'HLT_PFHT510_v15', 
        'HLT_PFHT590_v15', 
        'HLT_PFHT680_v15', 
        'HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10', 
        'HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10', 
        'HLT_PFHT780_v15', 
        'HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10', 
        'HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10', 
        'HLT_PFHT890_v15', 
        'HLT_PFJet140_v17', 
        'HLT_PFJet200_v17', 
        'HLT_PFJet260_v18', 
        'HLT_PFJet320_v18', 
        'HLT_PFJet400_v18', 
        'HLT_PFJet40_v19', 
        'HLT_PFJet450_v19', 
        'HLT_PFJet500_v19', 
        'HLT_PFJet550_v9', 
        'HLT_PFJet60_v19', 
        'HLT_PFJet80_v18', 
        'HLT_PFJetFwd140_v16', 
        'HLT_PFJetFwd200_v16', 
        'HLT_PFJetFwd260_v17', 
        'HLT_PFJetFwd320_v17', 
        'HLT_PFJetFwd400_v17', 
        'HLT_PFJetFwd40_v17', 
        'HLT_PFJetFwd450_v17', 
        'HLT_PFJetFwd500_v17', 
        'HLT_PFJetFwd60_v17', 
        'HLT_PFJetFwd80_v16', 
        'HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET100_PFMHT100_IDTight_PFHT60_v7', 
        'HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET110_PFMHT110_IDTight_v18', 
        'HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7', 
        'HLT_PFMET120_PFMHT120_IDTight_v18', 
        'HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET130_PFMHT130_IDTight_v18', 
        'HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v7', 
        'HLT_PFMET140_PFMHT140_IDTight_v18', 
        'HLT_PFMET200_HBHECleaned_v7', 
        'HLT_PFMET200_HBHE_BeamHaloCleaned_v7', 
        'HLT_PFMET200_NotCleaned_v7', 
        'HLT_PFMET250_HBHECleaned_v7', 
        'HLT_PFMET300_HBHECleaned_v7', 
        'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v7', 
        'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v18', 
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7', 
        'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v18', 
        'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v17', 
        'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17', 
        'HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v7', 
        'HLT_PFMETTypeOne110_PFMHT110_IDTight_v10', 
        'HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v7', 
        'HLT_PFMETTypeOne120_PFMHT120_IDTight_v10', 
        'HLT_PFMETTypeOne130_PFMHT130_IDTight_v10', 
        'HLT_PFMETTypeOne140_PFMHT140_IDTight_v9', 
        'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v7', 
        'HLT_Photon120_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon120_v12', 
        'HLT_Photon150_v5', 
        'HLT_Photon165_R9Id90_HE10_IsoM_v14', 
        'HLT_Photon175_v13', 
        'HLT_Photon200_v12', 
        'HLT_Photon20_HoverELoose_v9', 
        'HLT_Photon20_v1', 
        'HLT_Photon25_v3', 
        'HLT_Photon300_NoHE_v11', 
        'HLT_Photon30_HoverELoose_v9', 
        'HLT_Photon33_v4', 
        'HLT_Photon40_HoverELoose_v9', 
        'HLT_Photon50_HoverELoose_v9', 
        'HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_v3', 
        'HLT_Photon50_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon50_v12', 
        'HLT_Photon60_HoverELoose_v9', 
        'HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_v9', 
        'HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_v4', 
        'HLT_Photon60_R9Id90_CaloIdL_IsoL_v4', 
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v3', 
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_v3', 
        'HLT_Photon75_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon75_v12', 
        'HLT_Photon90_CaloIdL_PFHT700_v14', 
        'HLT_Photon90_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon90_v12', 
        'HLT_Physics_v7', 
        'HLT_QuadPFJet103_88_75_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet103_88_75_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_QuadPFJet103_88_75_15_v3', 
        'HLT_QuadPFJet105_88_76_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet105_88_76_15_v3', 
        'HLT_QuadPFJet105_90_76_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_QuadPFJet111_90_80_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet111_90_80_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_QuadPFJet111_90_80_15_v3', 
        'HLT_QuadPFJet98_83_71_15_BTagCSV_p013_VBF2_v7', 
        'HLT_QuadPFJet98_83_71_15_DoubleBTagCSV_p013_p08_VBF1_v7', 
        'HLT_QuadPFJet98_83_71_15_v3', 
        'HLT_Random_v3', 
        'HLT_Rsq0p35_v13', 
        'HLT_Rsq0p40_v13', 
        'HLT_RsqMR300_Rsq0p09_MR200_4jet_v13', 
        'HLT_RsqMR300_Rsq0p09_MR200_v13', 
        'HLT_RsqMR320_Rsq0p09_MR200_4jet_v13', 
        'HLT_RsqMR320_Rsq0p09_MR200_v13', 
        'HLT_SingleJet30_Mu12_SinglePFJet40_v9', 
        'HLT_SinglePhoton10_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton20_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton30_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton40_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton50_Eta3p1ForPPRef_v8', 
        'HLT_SinglePhoton60_Eta3p1ForPPRef_v8', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v2', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_v2', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1_v2', 
        'HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_v2', 
        'HLT_TkMu100_v2', 
        'HLT_Trimuon5_3p5_2_Upsilon_Muon_v4', 
        'HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v2', 
        'HLT_TripleJet110_35_35_Mjj650_PFMET110_v7', 
        'HLT_TripleJet110_35_35_Mjj650_PFMET120_v7', 
        'HLT_TripleJet110_35_35_Mjj650_PFMET130_v7', 
        'HLT_TripleMu_10_5_5_DZ_v9', 
        'HLT_TripleMu_12_10_5_v9', 
        'HLT_TripleMu_5_3_3_Mass3p8to60_DCA_v2', 
        'HLT_TripleMu_5_3_3_Mass3p8to60_DZ_v7', 
        'HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v2', 
        'HLT_TriplePhoton_20_20_20_CaloIdLV2_v2', 
        'HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_v3', 
        'HLT_TriplePhoton_30_30_10_CaloIdLV2_v3', 
        'HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v3', 
        'HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_v5', 
        'HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_v11', 
        'HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_v12', 
        'HLT_UncorrectedJetE30_NoBPTX3BX_v5', 
        'HLT_UncorrectedJetE30_NoBPTX_v5', 
        'HLT_UncorrectedJetE60_NoBPTX3BX_v5', 
        'HLT_UncorrectedJetE70_NoBPTX3BX_v5', 
        'HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v7', 
        'HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v7', 
        'HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v7', 
        'HLT_ZeroBias_FirstBXAfterTrain_v3', 
        'HLT_ZeroBias_FirstCollisionAfterAbortGap_v5', 
        'HLT_ZeroBias_FirstCollisionInTrain_v4', 
        'HLT_ZeroBias_IsolatedBunches_v5', 
        'HLT_ZeroBias_LastCollisionInTrain_v3', 
        'HLT_ZeroBias_v6'
     ) ),
    ParkingHT = cms.vstring(
        'DST_CaloJet40_BTagScouting_v13', 
        'DST_CaloJet40_CaloBTagScouting_v12', 
        'DST_CaloJet40_CaloScouting_PFScouting_v13', 
        'DST_HT250_CaloBTagScouting_v8', 
        'DST_HT250_CaloScouting_v9', 
        'DST_HT410_BTagScouting_v14', 
        'DST_HT410_PFScouting_v14', 
        'DST_L1HTT_BTagScouting_v13', 
        'DST_L1HTT_CaloBTagScouting_v12', 
        'DST_L1HTT_CaloScouting_PFScouting_v13', 
        'DST_ZeroBias_BTagScouting_v13', 
        'DST_ZeroBias_CaloScouting_PFScouting_v12'
    ),
    ParkingMuon = cms.vstring(
        'DST_DoubleMu3_noVtx_CaloScouting_v5', 
        'DST_L1DoubleMu_BTagScouting_v14', 
        'DST_L1DoubleMu_CaloScouting_PFScouting_v13'
    ),
    ParkingScoutingMonitor = cms.vstring(
        'DST_CaloJet40_BTagScouting_v13', 
        'DST_CaloJet40_CaloBTagScouting_v12', 
        'DST_CaloJet40_CaloScouting_PFScouting_v13', 
        'DST_DoubleMu3_noVtx_CaloScouting_Monitoring_v5', 
        'DST_DoubleMu3_noVtx_CaloScouting_v5', 
        'DST_HT250_CaloBTagScouting_v8', 
        'DST_HT250_CaloScouting_v9', 
        'DST_HT410_BTagScouting_v14', 
        'DST_HT410_PFScouting_v14', 
        'DST_L1DoubleMu_BTagScouting_v14', 
        'DST_L1DoubleMu_CaloScouting_PFScouting_v13', 
        'DST_L1HTT_BTagScouting_v13', 
        'DST_L1HTT_CaloBTagScouting_v12', 
        'DST_L1HTT_CaloScouting_PFScouting_v13', 
        'DST_ZeroBias_BTagScouting_v13', 
        'DST_ZeroBias_CaloScouting_PFScouting_v12', 
        'HLT_Ele115_CaloIdVT_GsfTrkIdT_v12', 
        'HLT_Ele35_WPTight_Gsf_v7', 
        'HLT_IsoMu27_v14', 
        'HLT_Mu50_v12', 
        'HLT_PFHT1050_v16', 
        'HLT_Photon200_v12'
    ),
    RPCMonitor = cms.vstring('AlCa_RPCMuonNormalisation_v13'),
    ScoutingCaloCommissioning = cms.vstring(
        'DST_CaloJet40_CaloBTagScouting_v12', 
        'DST_L1HTT_CaloBTagScouting_v12', 
        'DST_ZeroBias_CaloScouting_PFScouting_v12'
    ),
    ScoutingCaloHT = cms.vstring(
        'DST_HT250_CaloBTagScouting_v8', 
        'DST_HT250_CaloScouting_v9'
    ),
    ScoutingCaloMuon = cms.vstring(
        'DST_DoubleMu3_noVtx_CaloScouting_Monitoring_v5', 
        'DST_DoubleMu3_noVtx_CaloScouting_v5'
    ),
    ScoutingPFCommissioning = cms.vstring(
        'DST_CaloJet40_BTagScouting_v13', 
        'DST_CaloJet40_CaloScouting_PFScouting_v13', 
        'DST_L1DoubleMu_BTagScouting_v14', 
        'DST_L1DoubleMu_CaloScouting_PFScouting_v13', 
        'DST_L1HTT_BTagScouting_v13', 
        'DST_L1HTT_CaloScouting_PFScouting_v13', 
        'DST_ZeroBias_BTagScouting_v13', 
        'DST_ZeroBias_CaloScouting_PFScouting_v12'
    ),
    ScoutingPFHT = cms.vstring(
        'DST_HT410_BTagScouting_v14', 
        'DST_HT410_PFScouting_v14'
    ),
    SingleElectron = cms.vstring(
        'HLT_Ele115_CaloIdVT_GsfTrkIdT_v12', 
        'HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v16', 
        'HLT_Ele135_CaloIdVT_GsfTrkIdT_v5', 
        'HLT_Ele145_CaloIdVT_GsfTrkIdT_v6', 
        'HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v1', 
        'HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v7', 
        'HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v14', 
        'HLT_Ele15_IsoVVVL_PFHT450_v14', 
        'HLT_Ele15_IsoVVVL_PFHT600_v18', 
        'HLT_Ele15_WPLoose_Gsf_v1', 
        'HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v14', 
        'HLT_Ele17_WPLoose_Gsf_v1', 
        'HLT_Ele200_CaloIdVT_GsfTrkIdT_v6', 
        'HLT_Ele20_WPLoose_Gsf_v4', 
        'HLT_Ele20_WPTight_Gsf_v4', 
        'HLT_Ele20_eta2p1_WPLoose_Gsf_v4', 
        'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v16', 
        'HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v16', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v11', 
        'HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v11', 
        'HLT_Ele250_CaloIdVT_GsfTrkIdT_v11', 
        'HLT_Ele27_WPTight_Gsf_v14', 
        'HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v11', 
        'HLT_Ele300_CaloIdVT_GsfTrkIdT_v11', 
        'HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v11', 
        'HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7', 
        'HLT_Ele32_WPTight_Gsf_v13', 
        'HLT_Ele35_WPTight_Gsf_L1EGMT_v3', 
        'HLT_Ele35_WPTight_Gsf_v7', 
        'HLT_Ele38_WPTight_Gsf_v7', 
        'HLT_Ele40_WPTight_Gsf_v7', 
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v16', 
        'HLT_Ele50_IsoVVVL_PFHT450_v14', 
        'HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v14', 
        'HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v16'
    ),
    SingleMuon = cms.vstring(
        'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v10', 
        'HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v10', 
        'HLT_IsoMu20_v13', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v10', 
        'HLT_IsoMu24_eta2p1_v13', 
        'HLT_IsoMu24_v11', 
        'HLT_IsoMu27_v14', 
        'HLT_IsoMu30_v2', 
        'HLT_L1SingleMu18_v3', 
        'HLT_L1SingleMu25_v2', 
        'HLT_L1_DoubleJet30_Mass_Min400_Mu10_v1', 
        'HLT_L2Mu10_v7', 
        'HLT_L2Mu50_v2', 
        'HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v12', 
        'HLT_Mu12_v1', 
        'HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v7', 
        'HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v13', 
        'HLT_Mu15_IsoVVVL_PFHT450_v13', 
        'HLT_Mu15_IsoVVVL_PFHT600_v17', 
        'HLT_Mu15_v1', 
        'HLT_Mu20_v11', 
        'HLT_Mu27_v12', 
        'HLT_Mu3_PFJet40_v14', 
        'HLT_Mu50_IsoVVVL_PFHT450_v13', 
        'HLT_Mu50_v12', 
        'HLT_Mu55_v2', 
        'HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v13', 
        'HLT_OldMu100_v3', 
        'HLT_TkMu100_v2'
    ),
    SinglePhoton = cms.vstring(
        'HLT_Photon120_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon120_v12', 
        'HLT_Photon150_v5', 
        'HLT_Photon165_R9Id90_HE10_IsoM_v14', 
        'HLT_Photon175_v13', 
        'HLT_Photon200_v12', 
        'HLT_Photon20_HoverELoose_v9', 
        'HLT_Photon20_v1', 
        'HLT_Photon25_v3', 
        'HLT_Photon300_NoHE_v11', 
        'HLT_Photon30_HoverELoose_v9', 
        'HLT_Photon33_v4', 
        'HLT_Photon40_HoverELoose_v9', 
        'HLT_Photon50_HoverELoose_v9', 
        'HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_v3', 
        'HLT_Photon50_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon50_v12', 
        'HLT_Photon60_HoverELoose_v9', 
        'HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_v9', 
        'HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_v4', 
        'HLT_Photon60_R9Id90_CaloIdL_IsoL_v4', 
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v3', 
        'HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_v3', 
        'HLT_Photon75_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon75_v12', 
        'HLT_Photon90_CaloIdL_PFHT700_v14', 
        'HLT_Photon90_R9Id90_HE10_IsoM_v13', 
        'HLT_Photon90_v12'
    ),
    Tau = cms.vstring(
        'HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v10', 
        'HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v10', 
        'HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v10', 
        'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v10', 
        'HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1_v3', 
        'HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1_v3', 
        'HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1_v3', 
        'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v9', 
        'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v10', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v10', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v6', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_v6', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v6', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v10', 
        'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v10', 
        'HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v7', 
        'HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v7', 
        'HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v7'
    ),
    TestEnablesEcalHcal = cms.vstring(
        'HLT_EcalCalibration_v4', 
        'HLT_HcalCalibration_v5'
    ),
    TestEnablesEcalHcalDQM = cms.vstring(
        'HLT_EcalCalibration_v4', 
        'HLT_HcalCalibration_v5'
    ),
    ZeroBias = cms.vstring(
        'HLT_Random_v3', 
        'HLT_ZeroBias_FirstBXAfterTrain_v3', 
        'HLT_ZeroBias_FirstCollisionAfterAbortGap_v5', 
        'HLT_ZeroBias_FirstCollisionInTrain_v4', 
        'HLT_ZeroBias_IsolatedBunches_v5', 
        'HLT_ZeroBias_LastCollisionInTrain_v3', 
        'HLT_ZeroBias_v6'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1800)
)

process.options = cms.untracked.PSet(
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(4),
    sizeOfStackForThreadsInKB = cms.untracked.uint32(10240),
    wantSummary = cms.untracked.bool(True)
)

process.streams = cms.PSet(
    ALCALUMIPIXELS = cms.vstring('AlCaLumiPixels'),
    ALCAP0 = cms.vstring('AlCaP0'),
    ALCAPHISYM = cms.vstring('AlCaPhiSym'),
    Calibration = cms.vstring('TestEnablesEcalHcal'),
    DQM = cms.vstring('OnlineMonitor'),
    DQMCalibration = cms.vstring('TestEnablesEcalHcalDQM'),
    DQMEventDisplay = cms.vstring('EventDisplay'),
    DQMOnlineBeamspot = cms.vstring('DQMOnlineBeamspot'),
    EcalCalibration = cms.vstring('EcalLaser'),
    Express = cms.vstring('ExpressPhysics'),
    ExpressAlignment = cms.vstring('ExpressAlignment'),
    HLTMonitor = cms.vstring('HLTMonitor'),
    NanoDST = cms.vstring('L1Accept'),
    Parking = cms.vstring(
        'ParkingHT', 
        'ParkingMuon'
    ),
    PhysicsCommissioning = cms.vstring(
        'Commissioning', 
        'HLTPhysics', 
        'HcalNZS', 
        'HighPtLowerPhotons', 
        'HighPtPhoton30AndZ', 
        'IsolatedBunch', 
        'MonteCarlo', 
        'NoBPTX', 
        'ZeroBias'
    ),
    PhysicsEGamma = cms.vstring(
        'DoubleEG', 
        'SingleElectron', 
        'SinglePhoton'
    ),
    PhysicsEndOfFill = cms.vstring(
        'EmptyBX', 
        'FSQJet1', 
        'FSQJet2', 
        'HINCaloJets', 
        'HINPFJets', 
        'HighMultiplicityEOF', 
        'L1MinimumBias'
    ),
    PhysicsHLTPhysics1 = cms.vstring(
        'EphemeralHLTPhysics1', 
        'EphemeralHLTPhysics2'
    ),
    PhysicsHLTPhysics2 = cms.vstring(
        'EphemeralHLTPhysics3', 
        'EphemeralHLTPhysics4'
    ),
    PhysicsHLTPhysics3 = cms.vstring(
        'EphemeralHLTPhysics5', 
        'EphemeralHLTPhysics6'
    ),
    PhysicsHLTPhysics4 = cms.vstring(
        'EphemeralHLTPhysics7', 
        'EphemeralHLTPhysics8'
    ),
    PhysicsHadronsTaus = cms.vstring(
        'BTagCSV', 
        'BTagMu', 
        'DisplacedJet', 
        'HTMHT', 
        'JetHT', 
        'MET', 
        'Tau'
    ),
    PhysicsMuons = cms.vstring(
        'Charmonium', 
        'DoubleMuon', 
        'DoubleMuonLowMass', 
        'MuOnia', 
        'MuonEG', 
        'SingleMuon'
    ),
    PhysicsParkingScoutingMonitor = cms.vstring('ParkingScoutingMonitor'),
    PhysicsZeroBias1 = cms.vstring(
        'EphemeralZeroBias1', 
        'EphemeralZeroBias2'
    ),
    PhysicsZeroBias2 = cms.vstring(
        'EphemeralZeroBias3', 
        'EphemeralZeroBias4'
    ),
    PhysicsZeroBias3 = cms.vstring(
        'EphemeralZeroBias5', 
        'EphemeralZeroBias6'
    ),
    PhysicsZeroBias4 = cms.vstring(
        'EphemeralZeroBias7', 
        'EphemeralZeroBias8'
    ),
    RPCMON = cms.vstring('RPCMonitor'),
    ScoutingCaloMuon = cms.vstring(
        'ScoutingCaloCommissioning', 
        'ScoutingCaloHT', 
        'ScoutingCaloMuon'
    ),
    ScoutingPF = cms.vstring(
        'ScoutingPFCommissioning', 
        'ScoutingPFHT'
    )
)

process.transferSystem = cms.PSet(
    default = cms.PSet(
        default = cms.vstring('Tier0'),
        emulator = cms.vstring('Lustre'),
        streamLookArea = cms.PSet(

        ),
        test = cms.vstring('Lustre')
    ),
    destinations = cms.vstring(
        'Tier0', 
        'DQM', 
        'ECAL', 
        'EventDisplay', 
        'Lustre', 
        'None'
    ),
    streamA = cms.PSet(
        default = cms.vstring('Tier0'),
        emulator = cms.vstring('Lustre'),
        test = cms.vstring('Lustre')
    ),
    streamCalibration = cms.PSet(
        default = cms.vstring('Tier0'),
        emulator = cms.vstring('None'),
        test = cms.vstring('Lustre')
    ),
    streamDQM = cms.PSet(
        default = cms.vstring('DQM'),
        emulator = cms.vstring('None'),
        test = cms.vstring(
            'DQM', 
            'Lustre'
        )
    ),
    streamDQMCalibration = cms.PSet(
        default = cms.vstring('DQM'),
        emulator = cms.vstring('None'),
        test = cms.vstring(
            'DQM', 
            'Lustre'
        )
    ),
    streamEcalCalibration = cms.PSet(
        default = cms.vstring('ECAL'),
        emulator = cms.vstring('None'),
        test = cms.vstring('ECAL')
    ),
    streamEventDisplay = cms.PSet(
        default = cms.vstring(
            'EventDisplay', 
            'Tier0'
        ),
        emulator = cms.vstring('None'),
        test = cms.vstring(
            'EventDisplay', 
            'Lustre'
        )
    ),
    streamExpressCosmics = cms.PSet(
        default = cms.vstring('Tier0'),
        emulator = cms.vstring('Lustre'),
        test = cms.vstring('Lustre')
    ),
    streamLookArea = cms.PSet(
        default = cms.vstring('DQM'),
        emulator = cms.vstring('None'),
        test = cms.vstring(
            'DQM', 
            'Lustre'
        )
    ),
    streamNanoDST = cms.PSet(
        default = cms.vstring('Tier0'),
        emulator = cms.vstring('None'),
        test = cms.vstring('Lustre')
    ),
    streamRPCMON = cms.PSet(
        default = cms.vstring('Tier0'),
        emulator = cms.vstring('None'),
        test = cms.vstring('Lustre')
    ),
    streamTrackerCalibration = cms.PSet(
        default = cms.vstring('Tier0'),
        emulator = cms.vstring('None'),
        test = cms.vstring('Lustre')
    ),
    transferModes = cms.vstring(
        'default', 
        'test', 
        'emulator'
    )
)

process.hltAK4CaloAbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4CaloHLT'),
    level = cms.string('L3Absolute')
)


process.hltAK4CaloCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("hltAK4CaloFastJetCorrector", "hltAK4CaloRelativeCorrector", "hltAK4CaloAbsoluteCorrector", "hltAK4CaloResidualCorrector")
)


process.hltAK4CaloFastJetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4CaloHLT'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("hltFixedGridRhoFastjetAllCalo")
)


process.hltAK4CaloJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    DxyTrVtxMax = cms.double(0.0),
    DzTrVtxMax = cms.double(0.0),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    MaxVtxZ = cms.double(15.0),
    MinVtxNdof = cms.int32(5),
    R0 = cms.double(-1.0),
    Rho_EtaMax = cms.double(4.4),
    UseOnlyOnePV = cms.bool(False),
    UseOnlyVertexTracks = cms.bool(False),
    beta = cms.double(-1.0),
    correctShape = cms.bool(False),
    csRParam = cms.double(-1.0),
    csRho_EtaMax = cms.double(-1.0),
    dRMax = cms.double(-1.0),
    dRMin = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(True),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    gridMaxRapidity = cms.double(-1.0),
    gridSpacing = cms.double(-1.0),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.3),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string(''),
    jetPtMin = cms.double(1.0),
    jetType = cms.string('CaloJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxDepth = cms.int32(-1),
    maxInputs = cms.uint32(1),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    muCut = cms.double(-1.0),
    muMax = cms.double(-1.0),
    muMin = cms.double(-1.0),
    nExclude = cms.uint32(0),
    nFilt = cms.int32(-1),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(),
    puPtMin = cms.double(10.0),
    puWidth = cms.double(0.0),
    rFilt = cms.double(-1.0),
    rFiltFactor = cms.double(-1.0),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.4),
    rcut_factor = cms.double(-1.0),
    restrictInputs = cms.bool(False),
    src = cms.InputTag("hltTowerMakerForAll"),
    srcPVs = cms.InputTag("NotUsed"),
    subjetPtMin = cms.double(-1.0),
    subtractorName = cms.string(''),
    sumRecHits = cms.bool(False),
    trimPtFracMin = cms.double(-1.0),
    useCMSBoostedTauSeedingAlgorithm = cms.bool(False),
    useConstituentSubtraction = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    useDynamicFiltering = cms.bool(False),
    useExplicitGhosts = cms.bool(False),
    useFiltering = cms.bool(False),
    useKtPruning = cms.bool(False),
    useMassDropTagger = cms.bool(False),
    usePruning = cms.bool(False),
    useSoftDrop = cms.bool(False),
    useTrimming = cms.bool(False),
    verbosity = cms.int32(0),
    voronoiRfact = cms.double(0.9),
    writeCompound = cms.bool(False),
    writeJetsWithConst = cms.bool(False),
    yCut = cms.double(-1.0),
    yMax = cms.double(-1.0),
    yMin = cms.double(-1.0),
    zcut = cms.double(-1.0)
)


process.hltAK4CaloJetsCorrected = cms.EDProducer("CorrectedCaloJetProducer",
    correctors = cms.VInputTag("hltAK4CaloCorrector"),
    src = cms.InputTag("hltAK4CaloJets")
)


process.hltAK4CaloJetsCorrectedIDPassed = cms.EDProducer("CorrectedCaloJetProducer",
    correctors = cms.VInputTag("hltAK4CaloCorrector"),
    src = cms.InputTag("hltAK4CaloJetsIDPassed")
)


process.hltAK4CaloJetsIDPassed = cms.EDProducer("HLTCaloJetIDProducer",
    JetIDParams = cms.PSet(
        ebRecHitsColl = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
        eeRecHitsColl = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
        hbheRecHitsColl = cms.InputTag("hltHbhereco"),
        hfRecHitsColl = cms.InputTag("hltHfreco"),
        hoRecHitsColl = cms.InputTag("hltHoreco"),
        useRecHits = cms.bool(True)
    ),
    jetsInput = cms.InputTag("hltAK4CaloJets"),
    max_EMF = cms.double(999.0),
    min_EMF = cms.double(1e-06),
    min_N90 = cms.int32(-2),
    min_N90hits = cms.int32(2)
)


process.hltAK4CaloJetsPF = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    DxyTrVtxMax = cms.double(0.0),
    DzTrVtxMax = cms.double(0.0),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    MaxVtxZ = cms.double(15.0),
    MinVtxNdof = cms.int32(5),
    R0 = cms.double(-1.0),
    Rho_EtaMax = cms.double(4.4),
    UseOnlyOnePV = cms.bool(False),
    UseOnlyVertexTracks = cms.bool(False),
    beta = cms.double(-1.0),
    correctShape = cms.bool(False),
    csRParam = cms.double(-1.0),
    csRho_EtaMax = cms.double(-1.0),
    dRMax = cms.double(-1.0),
    dRMin = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    gridMaxRapidity = cms.double(-1.0),
    gridSpacing = cms.double(-1.0),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.3),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string(''),
    jetPtMin = cms.double(1.0),
    jetType = cms.string('CaloJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxDepth = cms.int32(-1),
    maxInputs = cms.uint32(1),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(0),
    muCut = cms.double(-1.0),
    muMax = cms.double(-1.0),
    muMin = cms.double(-1.0),
    nExclude = cms.uint32(0),
    nFilt = cms.int32(-1),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(),
    puPtMin = cms.double(10.0),
    puWidth = cms.double(0.0),
    rFilt = cms.double(-1.0),
    rFiltFactor = cms.double(-1.0),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.4),
    rcut_factor = cms.double(-1.0),
    restrictInputs = cms.bool(False),
    src = cms.InputTag("hltTowerMakerForAll"),
    srcPVs = cms.InputTag("NotUsed"),
    subjetPtMin = cms.double(-1.0),
    subtractorName = cms.string(''),
    sumRecHits = cms.bool(False),
    trimPtFracMin = cms.double(-1.0),
    useCMSBoostedTauSeedingAlgorithm = cms.bool(False),
    useConstituentSubtraction = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    useDynamicFiltering = cms.bool(False),
    useExplicitGhosts = cms.bool(False),
    useFiltering = cms.bool(False),
    useKtPruning = cms.bool(False),
    useMassDropTagger = cms.bool(False),
    usePruning = cms.bool(False),
    useSoftDrop = cms.bool(False),
    useTrimming = cms.bool(False),
    verbosity = cms.int32(0),
    voronoiRfact = cms.double(-9.0),
    writeCompound = cms.bool(False),
    writeJetsWithConst = cms.bool(False),
    yCut = cms.double(-1.0),
    yMax = cms.double(-1.0),
    yMin = cms.double(-1.0),
    zcut = cms.double(-1.0)
)


process.hltAK4CaloRelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4CaloHLT'),
    level = cms.string('L2Relative')
)


process.hltAK4CaloResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4CaloHLT'),
    level = cms.string('L2L3Residual')
)


process.hltAK4Iter1TrackJets4Iter2 = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    DxyTrVtxMax = cms.double(0.2),
    DzTrVtxMax = cms.double(0.5),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    MaxVtxZ = cms.double(30.0),
    MinVtxNdof = cms.int32(0),
    R0 = cms.double(-1.0),
    Rho_EtaMax = cms.double(4.4),
    UseOnlyOnePV = cms.bool(True),
    UseOnlyVertexTracks = cms.bool(False),
    beta = cms.double(-1.0),
    correctShape = cms.bool(False),
    csRParam = cms.double(-1.0),
    csRho_EtaMax = cms.double(-1.0),
    dRMax = cms.double(-1.0),
    dRMin = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    gridMaxRapidity = cms.double(-1.0),
    gridSpacing = cms.double(-1.0),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.1),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string(''),
    jetPtMin = cms.double(7.5),
    jetType = cms.string('TrackJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxDepth = cms.int32(-1),
    maxInputs = cms.uint32(1),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    muCut = cms.double(-1.0),
    muMax = cms.double(-1.0),
    muMin = cms.double(-1.0),
    nExclude = cms.uint32(0),
    nFilt = cms.int32(-1),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(),
    puPtMin = cms.double(0.0),
    puWidth = cms.double(0.0),
    rFilt = cms.double(-1.0),
    rFiltFactor = cms.double(-1.0),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.4),
    rcut_factor = cms.double(-1.0),
    restrictInputs = cms.bool(False),
    src = cms.InputTag("hltIter1TrackRefsForJets4Iter2"),
    srcPVs = cms.InputTag("hltTrimmedPixelVertices"),
    subjetPtMin = cms.double(-1.0),
    subtractorName = cms.string(''),
    sumRecHits = cms.bool(False),
    trimPtFracMin = cms.double(-1.0),
    useCMSBoostedTauSeedingAlgorithm = cms.bool(False),
    useConstituentSubtraction = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    useDynamicFiltering = cms.bool(False),
    useExplicitGhosts = cms.bool(False),
    useFiltering = cms.bool(False),
    useKtPruning = cms.bool(False),
    useMassDropTagger = cms.bool(False),
    usePruning = cms.bool(False),
    useSoftDrop = cms.bool(False),
    useTrimming = cms.bool(False),
    verbosity = cms.int32(0),
    voronoiRfact = cms.double(0.9),
    writeCompound = cms.bool(False),
    writeJetsWithConst = cms.bool(False),
    yCut = cms.double(-1.0),
    yMax = cms.double(-1.0),
    yMin = cms.double(-1.0),
    zcut = cms.double(-1.0)
)


process.hltAK4PFAbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFHLT'),
    level = cms.string('L3Absolute')
)


process.hltAK4PFCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("hltAK4PFFastJetCorrector", "hltAK4PFRelativeCorrector", "hltAK4PFAbsoluteCorrector", "hltAK4PFResidualCorrector")
)


process.hltAK4PFFastJetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PFHLT'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("hltFixedGridRhoFastjetAll")
)


process.hltAK4PFJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    DxyTrVtxMax = cms.double(0.0),
    DzTrVtxMax = cms.double(0.0),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    MaxVtxZ = cms.double(15.0),
    MinVtxNdof = cms.int32(0),
    R0 = cms.double(-1.0),
    Rho_EtaMax = cms.double(4.4),
    UseOnlyOnePV = cms.bool(False),
    UseOnlyVertexTracks = cms.bool(False),
    beta = cms.double(-1.0),
    correctShape = cms.bool(False),
    csRParam = cms.double(-1.0),
    csRho_EtaMax = cms.double(-1.0),
    dRMax = cms.double(-1.0),
    dRMin = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(True),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    gridMaxRapidity = cms.double(-1.0),
    gridSpacing = cms.double(-1.0),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.0),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string(''),
    jetPtMin = cms.double(0.0),
    jetType = cms.string('PFJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxDepth = cms.int32(-1),
    maxInputs = cms.uint32(1),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(0),
    muCut = cms.double(-1.0),
    muMax = cms.double(-1.0),
    muMin = cms.double(-1.0),
    nExclude = cms.uint32(0),
    nFilt = cms.int32(-1),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(),
    puPtMin = cms.double(10.0),
    puWidth = cms.double(0.0),
    rFilt = cms.double(-1.0),
    rFiltFactor = cms.double(-1.0),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.4),
    rcut_factor = cms.double(-1.0),
    restrictInputs = cms.bool(False),
    src = cms.InputTag("hltParticleFlow"),
    srcPVs = cms.InputTag("hltPixelVertices"),
    subjetPtMin = cms.double(-1.0),
    subtractorName = cms.string(''),
    sumRecHits = cms.bool(False),
    trimPtFracMin = cms.double(-1.0),
    useCMSBoostedTauSeedingAlgorithm = cms.bool(False),
    useConstituentSubtraction = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    useDynamicFiltering = cms.bool(False),
    useExplicitGhosts = cms.bool(False),
    useFiltering = cms.bool(False),
    useKtPruning = cms.bool(False),
    useMassDropTagger = cms.bool(False),
    usePruning = cms.bool(False),
    useSoftDrop = cms.bool(False),
    useTrimming = cms.bool(False),
    verbosity = cms.int32(0),
    voronoiRfact = cms.double(-9.0),
    writeCompound = cms.bool(False),
    writeJetsWithConst = cms.bool(False),
    yCut = cms.double(-1.0),
    yMax = cms.double(-1.0),
    yMin = cms.double(-1.0),
    zcut = cms.double(-1.0)
)


process.hltAK4PFJetsCorrected = cms.EDProducer("CorrectedPFJetProducer",
    correctors = cms.VInputTag("hltAK4PFCorrector"),
    src = cms.InputTag("hltAK4PFJets")
)


process.hltAK4PFJetsLooseID = cms.EDProducer("HLTPFJetIDProducer",
    CEF = cms.double(0.99),
    CHF = cms.double(0.0),
    NCH = cms.int32(0),
    NEF = cms.double(0.99),
    NHF = cms.double(0.99),
    NTOT = cms.int32(1),
    jetsInput = cms.InputTag("hltAK4PFJets"),
    maxCF = cms.double(99.0),
    maxEta = cms.double(1e+99),
    minPt = cms.double(20.0)
)


process.hltAK4PFJetsLooseIDCorrected = cms.EDProducer("CorrectedPFJetProducer",
    correctors = cms.VInputTag("hltAK4PFCorrector"),
    src = cms.InputTag("hltAK4PFJetsLooseID")
)


process.hltAK4PFJetsTightID = cms.EDProducer("HLTPFJetIDProducer",
    CEF = cms.double(0.99),
    CHF = cms.double(0.0),
    NCH = cms.int32(0),
    NEF = cms.double(0.99),
    NHF = cms.double(0.9),
    NTOT = cms.int32(1),
    jetsInput = cms.InputTag("hltAK4PFJets"),
    maxCF = cms.double(99.0),
    maxEta = cms.double(1e+99),
    minPt = cms.double(20.0)
)


process.hltAK4PFJetsTightIDCorrected = cms.EDProducer("CorrectedPFJetProducer",
    correctors = cms.VInputTag("hltAK4PFCorrector"),
    src = cms.InputTag("hltAK4PFJetsTightID")
)


process.hltAK4PFRelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFHLT'),
    level = cms.string('L2Relative')
)


process.hltAK4PFResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFHLT'),
    level = cms.string('L2L3Residual')
)


process.hltAK8CaloAbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK8CaloHLT'),
    level = cms.string('L3Absolute')
)


process.hltAK8CaloCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("hltAK8CaloFastJetCorrector", "hltAK8CaloRelativeCorrector", "hltAK8CaloAbsoluteCorrector", "hltAK8CaloResidualCorrector")
)


process.hltAK8CaloFastJetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK8CaloHLT'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("hltFixedGridRhoFastjetAllCalo")
)


process.hltAK8CaloJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    DxyTrVtxMax = cms.double(0.0),
    DzTrVtxMax = cms.double(0.0),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    MaxVtxZ = cms.double(15.0),
    MinVtxNdof = cms.int32(5),
    R0 = cms.double(-1.0),
    Rho_EtaMax = cms.double(4.4),
    UseOnlyOnePV = cms.bool(False),
    UseOnlyVertexTracks = cms.bool(False),
    beta = cms.double(-1.0),
    correctShape = cms.bool(False),
    csRParam = cms.double(-1.0),
    csRho_EtaMax = cms.double(-1.0),
    dRMax = cms.double(-1.0),
    dRMin = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(True),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    gridMaxRapidity = cms.double(-1.0),
    gridSpacing = cms.double(-1.0),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.3),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string(''),
    jetPtMin = cms.double(1.0),
    jetType = cms.string('CaloJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxDepth = cms.int32(-1),
    maxInputs = cms.uint32(1),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    muCut = cms.double(-1.0),
    muMax = cms.double(-1.0),
    muMin = cms.double(-1.0),
    nExclude = cms.uint32(0),
    nFilt = cms.int32(-1),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(),
    puPtMin = cms.double(10.0),
    puWidth = cms.double(0.0),
    rFilt = cms.double(-1.0),
    rFiltFactor = cms.double(-1.0),
    rParam = cms.double(0.8),
    radiusPU = cms.double(0.8),
    rcut_factor = cms.double(-1.0),
    restrictInputs = cms.bool(False),
    src = cms.InputTag("hltTowerMakerForAll"),
    srcPVs = cms.InputTag("NotUsed"),
    subjetPtMin = cms.double(-1.0),
    subtractorName = cms.string(''),
    sumRecHits = cms.bool(False),
    trimPtFracMin = cms.double(-1.0),
    useCMSBoostedTauSeedingAlgorithm = cms.bool(False),
    useConstituentSubtraction = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    useDynamicFiltering = cms.bool(False),
    useExplicitGhosts = cms.bool(False),
    useFiltering = cms.bool(False),
    useKtPruning = cms.bool(False),
    useMassDropTagger = cms.bool(False),
    usePruning = cms.bool(False),
    useSoftDrop = cms.bool(False),
    useTrimming = cms.bool(False),
    verbosity = cms.int32(0),
    voronoiRfact = cms.double(0.9),
    writeCompound = cms.bool(False),
    writeJetsWithConst = cms.bool(False),
    yCut = cms.double(-1.0),
    yMax = cms.double(-1.0),
    yMin = cms.double(-1.0),
    zcut = cms.double(-1.0)
)


process.hltAK8CaloJetsCorrected = cms.EDProducer("CorrectedCaloJetProducer",
    correctors = cms.VInputTag("hltAK8CaloCorrector"),
    src = cms.InputTag("hltAK8CaloJets")
)


process.hltAK8CaloJetsCorrectedIDPassed = cms.EDProducer("CorrectedCaloJetProducer",
    correctors = cms.VInputTag("hltAK8CaloCorrector"),
    src = cms.InputTag("hltAK8CaloJetsIDPassed")
)


process.hltAK8CaloJetsIDPassed = cms.EDProducer("HLTCaloJetIDProducer",
    JetIDParams = cms.PSet(
        ebRecHitsColl = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
        eeRecHitsColl = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
        hbheRecHitsColl = cms.InputTag("hltHbhereco"),
        hfRecHitsColl = cms.InputTag("hltHfreco"),
        hoRecHitsColl = cms.InputTag("hltHoreco"),
        useRecHits = cms.bool(True)
    ),
    jetsInput = cms.InputTag("hltAK8CaloJets"),
    max_EMF = cms.double(999.0),
    min_EMF = cms.double(1e-06),
    min_N90 = cms.int32(-2),
    min_N90hits = cms.int32(2)
)


process.hltAK8CaloJetsPF = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    DxyTrVtxMax = cms.double(0.0),
    DzTrVtxMax = cms.double(0.0),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    MaxVtxZ = cms.double(15.0),
    MinVtxNdof = cms.int32(5),
    R0 = cms.double(-1.0),
    Rho_EtaMax = cms.double(4.4),
    UseOnlyOnePV = cms.bool(False),
    UseOnlyVertexTracks = cms.bool(False),
    beta = cms.double(-1.0),
    correctShape = cms.bool(False),
    csRParam = cms.double(-1.0),
    csRho_EtaMax = cms.double(-1.0),
    dRMax = cms.double(-1.0),
    dRMin = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    gridMaxRapidity = cms.double(-1.0),
    gridSpacing = cms.double(-1.0),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.3),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string(''),
    jetPtMin = cms.double(1.0),
    jetType = cms.string('CaloJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxDepth = cms.int32(-1),
    maxInputs = cms.uint32(1),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(0),
    muCut = cms.double(-1.0),
    muMax = cms.double(-1.0),
    muMin = cms.double(-1.0),
    nExclude = cms.uint32(0),
    nFilt = cms.int32(-1),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(),
    puPtMin = cms.double(10.0),
    puWidth = cms.double(0.0),
    rFilt = cms.double(-1.0),
    rFiltFactor = cms.double(-1.0),
    rParam = cms.double(0.8),
    radiusPU = cms.double(0.8),
    rcut_factor = cms.double(-1.0),
    restrictInputs = cms.bool(False),
    src = cms.InputTag("hltTowerMakerForAll"),
    srcPVs = cms.InputTag("NotUsed"),
    subjetPtMin = cms.double(-1.0),
    subtractorName = cms.string(''),
    sumRecHits = cms.bool(False),
    trimPtFracMin = cms.double(-1.0),
    useCMSBoostedTauSeedingAlgorithm = cms.bool(False),
    useConstituentSubtraction = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    useDynamicFiltering = cms.bool(False),
    useExplicitGhosts = cms.bool(False),
    useFiltering = cms.bool(False),
    useKtPruning = cms.bool(False),
    useMassDropTagger = cms.bool(False),
    usePruning = cms.bool(False),
    useSoftDrop = cms.bool(False),
    useTrimming = cms.bool(False),
    verbosity = cms.int32(0),
    voronoiRfact = cms.double(-9.0),
    writeCompound = cms.bool(False),
    writeJetsWithConst = cms.bool(False),
    yCut = cms.double(-1.0),
    yMax = cms.double(-1.0),
    yMin = cms.double(-1.0),
    zcut = cms.double(-1.0)
)


process.hltAK8CaloRelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK8CaloHLT'),
    level = cms.string('L2Relative')
)


process.hltAK8CaloResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK8CaloHLT'),
    level = cms.string('L2L3Residual')
)


process.hltAK8PFAbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK8PFHLT'),
    level = cms.string('L3Absolute')
)


process.hltAK8PFCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("hltAK8PFFastJetCorrector", "hltAK8PFRelativeCorrector", "hltAK8PFAbsoluteCorrector", "hltAK8PFResidualCorrector")
)


process.hltAK8PFFastJetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK8PFHLT'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("hltFixedGridRhoFastjetAll")
)


process.hltAK8PFJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(5),
    DxyTrVtxMax = cms.double(0.0),
    DzTrVtxMax = cms.double(0.0),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    MaxVtxZ = cms.double(15.0),
    MinVtxNdof = cms.int32(0),
    R0 = cms.double(-1.0),
    Rho_EtaMax = cms.double(4.4),
    UseOnlyOnePV = cms.bool(False),
    UseOnlyVertexTracks = cms.bool(False),
    beta = cms.double(-1.0),
    correctShape = cms.bool(False),
    csRParam = cms.double(-1.0),
    csRho_EtaMax = cms.double(-1.0),
    dRMax = cms.double(-1.0),
    dRMin = cms.double(-1.0),
    doAreaDiskApprox = cms.bool(True),
    doAreaFastjet = cms.bool(False),
    doFastJetNonUniform = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(False),
    doRhoFastjet = cms.bool(False),
    gridMaxRapidity = cms.double(-1.0),
    gridSpacing = cms.double(-1.0),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.0),
    jetAlgorithm = cms.string('AntiKt'),
    jetCollInstanceName = cms.string(''),
    jetPtMin = cms.double(0.0),
    jetType = cms.string('PFJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxDepth = cms.int32(-1),
    maxInputs = cms.uint32(1),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(0),
    muCut = cms.double(-1.0),
    muMax = cms.double(-1.0),
    muMin = cms.double(-1.0),
    nExclude = cms.uint32(0),
    nFilt = cms.int32(-1),
    nSigmaPU = cms.double(1.0),
    puCenters = cms.vdouble(),
    puPtMin = cms.double(10.0),
    puWidth = cms.double(0.0),
    rFilt = cms.double(-1.0),
    rFiltFactor = cms.double(-1.0),
    rParam = cms.double(0.8),
    radiusPU = cms.double(0.8),
    rcut_factor = cms.double(-1.0),
    restrictInputs = cms.bool(False),
    src = cms.InputTag("hltParticleFlow"),
    srcPVs = cms.InputTag("hltPixelVertices"),
    subjetPtMin = cms.double(-1.0),
    subtractorName = cms.string(''),
    sumRecHits = cms.bool(False),
    trimPtFracMin = cms.double(-1.0),
    useCMSBoostedTauSeedingAlgorithm = cms.bool(False),
    useConstituentSubtraction = cms.bool(False),
    useDeterministicSeed = cms.bool(True),
    useDynamicFiltering = cms.bool(False),
    useExplicitGhosts = cms.bool(False),
    useFiltering = cms.bool(False),
    useKtPruning = cms.bool(False),
    useMassDropTagger = cms.bool(False),
    usePruning = cms.bool(False),
    useSoftDrop = cms.bool(False),
    useTrimming = cms.bool(False),
    verbosity = cms.int32(0),
    voronoiRfact = cms.double(-9.0),
    writeCompound = cms.bool(False),
    writeJetsWithConst = cms.bool(False),
    yCut = cms.double(-1.0),
    yMax = cms.double(-1.0),
    yMin = cms.double(-1.0),
    zcut = cms.double(-1.0)
)


process.hltAK8PFJetsCorrected = cms.EDProducer("CorrectedPFJetProducer",
    correctors = cms.VInputTag("hltAK8PFCorrector"),
    src = cms.InputTag("hltAK8PFJets")
)


process.hltAK8PFJetsCorrectedMatchedToCaloJets300 = cms.EDProducer("PFJetsMatchedToFilteredCaloJetsProducer",
    CaloJetFilter = cms.InputTag("hltAK8SingleCaloJet300"),
    DeltaR = cms.double(0.5),
    PFJetSrc = cms.InputTag("hltAK8PFJetsCorrected"),
    TriggerType = cms.int32(85)
)


process.hltAK8PFJetsLooseID = cms.EDProducer("HLTPFJetIDProducer",
    CEF = cms.double(0.99),
    CHF = cms.double(0.0),
    NCH = cms.int32(0),
    NEF = cms.double(0.99),
    NHF = cms.double(0.99),
    NTOT = cms.int32(1),
    jetsInput = cms.InputTag("hltAK8PFJets"),
    maxCF = cms.double(99.0),
    maxEta = cms.double(1e+99),
    minPt = cms.double(20.0)
)


process.hltAK8PFJetsLooseIDCorrected = cms.EDProducer("CorrectedPFJetProducer",
    correctors = cms.VInputTag("hltAK8PFCorrector"),
    src = cms.InputTag("hltAK8PFJetsLooseID")
)


process.hltAK8PFJetsTightID = cms.EDProducer("HLTPFJetIDProducer",
    CEF = cms.double(0.99),
    CHF = cms.double(0.0),
    NCH = cms.int32(0),
    NEF = cms.double(0.99),
    NHF = cms.double(0.9),
    NTOT = cms.int32(1),
    jetsInput = cms.InputTag("hltAK8PFJets"),
    maxCF = cms.double(99.0),
    maxEta = cms.double(1e+99),
    minPt = cms.double(20.0)
)


process.hltAK8PFJetsTightIDCorrected = cms.EDProducer("CorrectedPFJetProducer",
    correctors = cms.VInputTag("hltAK8PFCorrector"),
    src = cms.InputTag("hltAK8PFJetsTightID")
)


process.hltAK8PFRelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK8PFHLT'),
    level = cms.string('L2Relative')
)


process.hltAK8PFResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK8PFHLT'),
    level = cms.string('L2L3Residual')
)


process.hltBLifetimeAssociatorPF = cms.EDProducer("JetTracksAssociatorAtVertex",
    coneSize = cms.double(0.4),
    jets = cms.InputTag("hltPFJetForBtag"),
    pvSrc = cms.InputTag(""),
    tracks = cms.InputTag("hltPFMuonMerging"),
    useAssigned = cms.bool(False)
)


process.hltBLifetimeAssociatorPFAK8 = cms.EDProducer("JetTracksAssociatorAtVertex",
    coneSize = cms.double(0.4),
    jets = cms.InputTag("hltPFJetForBtagAK8"),
    pvSrc = cms.InputTag(""),
    tracks = cms.InputTag("hltPFMuonMerging"),
    useAssigned = cms.bool(False)
)


process.hltBLifetimeTagInfosPF = cms.EDProducer("TrackIPProducer",
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jetTracks = cms.InputTag("hltBLifetimeAssociatorPF"),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    minimumNumberOfPixelHits = cms.int32(2),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("hltVerticesPFFilter"),
    useTrackQuality = cms.bool(False)
)


process.hltBLifetimeTagInfosPFAK8 = cms.EDProducer("TrackIPProducer",
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jetTracks = cms.InputTag("hltBLifetimeAssociatorPFAK8"),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    minimumNumberOfPixelHits = cms.int32(2),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("hltVerticesPFFilter"),
    useTrackQuality = cms.bool(False)
)


process.hltBoostedDBInclusiveMergedVerticesPFAK8 = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.0),
    secondaryVertices = cms.InputTag("hltBoostedDBTrackVertexArbitratorPFAK8")
)


process.hltBoostedDBInclusiveSecondaryVerticesPFAK8 = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.7),
    minSignificance = cms.double(2.0),
    secondaryVertices = cms.InputTag("hltBoostedDBInclusiveVertexFinderPFAK8")
)


process.hltBoostedDBInclusiveVertexFinderPFAK8 = cms.EDProducer("InclusiveCandidateVertexFinder",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterizer = cms.PSet(
        clusterMaxDistance = cms.double(0.05),
        clusterMaxSignificance = cms.double(4.5),
        clusterMinAngleCosine = cms.double(0.5),
        distanceRatio = cms.double(20.0),
        seedMax3DIPSignificance = cms.double(9999.0),
        seedMax3DIPValue = cms.double(9999.0),
        seedMin3DIPSignificance = cms.double(1.2),
        seedMin3DIPValue = cms.double(0.005)
    ),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxNTracks = cms.uint32(30),
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3.0),
    minHits = cms.uint32(8),
    minPt = cms.double(0.8),
    primaryVertices = cms.InputTag("hltVerticesPFFilter"),
    tracks = cms.InputTag("hltParticleFlow"),
    useDirectVertexFitter = cms.bool(True),
    useVertexReco = cms.bool(True),
    vertexMinAngleCosine = cms.double(0.95),
    vertexMinDLen2DSig = cms.double(2.5),
    vertexMinDLenSig = cms.double(0.5),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        primcut = cms.double(1.0),
        seccut = cms.double(3.0),
        smoothing = cms.bool(True)
    )
)


process.hltBoostedDBJetTagsPFAK8 = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('hltBoostedDoubleSecondaryVertexAK8Computer'),
    tagInfos = cms.VInputTag("hltBoostedDBSVAK8TagInfosPF")
)


process.hltBoostedDBLifetimeTagInfosPFAK8 = cms.EDProducer("CandIPProducer",
    candidates = cms.InputTag("hltParticleFlow"),
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jets = cms.InputTag("hltPFJetForBtagAK8"),
    maxDeltaR = cms.double(0.4),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    minimumNumberOfPixelHits = cms.int32(2),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("hltVerticesPFFilter"),
    useTrackQuality = cms.bool(False)
)


process.hltBoostedDBSVAK8TagInfosPF = cms.EDProducer("BoostedDoubleSVProducer",
    R0 = cms.double(0.8),
    beta = cms.double(1.0),
    maxSVDeltaRToJet = cms.double(0.7),
    svTagInfos = cms.InputTag("hltBoostedDBSecondaryVertexTagInfosPFAK8"),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackSelection = cms.PSet(
        jetDeltaRMax = cms.double(0.8),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False),
        variableJTAParsi = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5)
        )
    )
)


process.hltBoostedDBSecondaryVertexTagInfosPFAK8 = cms.EDProducer("CandSecondaryVertexProducer",
    beamSpotTag = cms.InputTag("hltOnlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("hltBoostedDBInclusiveMergedVerticesPFAK8"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("hltBoostedDBLifetimeTagInfosPFAK8"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500.0),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3.0),
        min_pT = cms.double(120.0),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(2),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(3),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(True),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(3.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.65),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.4),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.hltBoostedDBTrackVertexArbitratorPFAK8 = cms.EDProducer("CandidateVertexArbitrator",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxTimeSignificance = cms.double(3.5),
    primaryVertices = cms.InputTag("hltVerticesPFFilter"),
    secondaryVertices = cms.InputTag("hltBoostedDBInclusiveSecondaryVerticesPFAK8"),
    sigCut = cms.double(5.0),
    trackMinLayers = cms.int32(4),
    trackMinPixels = cms.int32(1),
    trackMinPt = cms.double(0.4),
    tracks = cms.InputTag("hltParticleFlow")
)


process.hltCaloJetsSix25ForHt = cms.EDProducer("HLTCaloJetCollectionProducer",
    HLTObject = cms.InputTag("hltCaloJetFilterSixC25"),
    TriggerTypes = cms.vint32(86)
)


process.hltCombinedSecondaryVertexBJetTagsCalo = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('hltCombinedSecondaryVertexV2'),
    tagInfos = cms.VInputTag("hltImpactParameterTagInfos", "hltInclusiveSecondaryVertexFinderTagInfos")
)


process.hltCombinedSecondaryVertexBJetTagsPF = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('hltCombinedSecondaryVertexV2'),
    tagInfos = cms.VInputTag("hltBLifetimeTagInfosPF", "hltSecondaryVertexTagInfosPF")
)


process.hltCombinedSecondaryVertexBJetTagsPFAK8 = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('hltCombinedSecondaryVertexV2'),
    tagInfos = cms.VInputTag("hltBLifetimeTagInfosPFAK8", "hltSecondaryVertexTagInfosPFAK8")
)


process.hltCsc2DRecHits = cms.EDProducer("CSCRecHitDProducer",
    CSCDebug = cms.untracked.bool(False),
    CSCNoOfTimeBinsForDynamicPedestal = cms.int32(2),
    CSCStripClusterChargeCut = cms.double(25.0),
    CSCStripClusterSize = cms.untracked.int32(3),
    CSCStripPeakThreshold = cms.double(10.0),
    CSCStripxtalksOffset = cms.double(0.03),
    CSCUseCalibrations = cms.bool(True),
    CSCUseGasGainCorrections = cms.bool(False),
    CSCUseReducedWireTimeWindow = cms.bool(False),
    CSCUseStaticPedestals = cms.bool(False),
    CSCUseTimingCorrections = cms.bool(True),
    CSCWireClusterDeltaT = cms.int32(1),
    CSCWireTimeWindowHigh = cms.int32(15),
    CSCWireTimeWindowLow = cms.int32(0),
    CSCstripWireDeltaTime = cms.int32(8),
    ConstSyst_ME12 = cms.double(0.0),
    ConstSyst_ME13 = cms.double(0.0),
    ConstSyst_ME1a = cms.double(0.022),
    ConstSyst_ME1b = cms.double(0.007),
    ConstSyst_ME21 = cms.double(0.0),
    ConstSyst_ME22 = cms.double(0.0),
    ConstSyst_ME31 = cms.double(0.0),
    ConstSyst_ME32 = cms.double(0.0),
    ConstSyst_ME41 = cms.double(0.0),
    NoiseLevel_ME12 = cms.double(9.0),
    NoiseLevel_ME13 = cms.double(8.0),
    NoiseLevel_ME1a = cms.double(7.0),
    NoiseLevel_ME1b = cms.double(8.0),
    NoiseLevel_ME21 = cms.double(9.0),
    NoiseLevel_ME22 = cms.double(9.0),
    NoiseLevel_ME31 = cms.double(9.0),
    NoiseLevel_ME32 = cms.double(9.0),
    NoiseLevel_ME41 = cms.double(9.0),
    UseAverageTime = cms.bool(False),
    UseFivePoleFit = cms.bool(True),
    UseParabolaFit = cms.bool(False),
    XTasymmetry_ME12 = cms.double(0.0),
    XTasymmetry_ME13 = cms.double(0.0),
    XTasymmetry_ME1a = cms.double(0.0),
    XTasymmetry_ME1b = cms.double(0.0),
    XTasymmetry_ME21 = cms.double(0.0),
    XTasymmetry_ME22 = cms.double(0.0),
    XTasymmetry_ME31 = cms.double(0.0),
    XTasymmetry_ME32 = cms.double(0.0),
    XTasymmetry_ME41 = cms.double(0.0),
    readBadChambers = cms.bool(True),
    readBadChannels = cms.bool(False),
    stripDigiTag = cms.InputTag("hltMuonCSCDigis","MuonCSCStripDigi"),
    wireDigiTag = cms.InputTag("hltMuonCSCDigis","MuonCSCWireDigi")
)


process.hltCscSegments = cms.EDProducer("CSCSegmentProducer",
    algo_psets = cms.VPSet(cms.PSet(
        algo_name = cms.string('CSCSegAlgoST'),
        algo_psets = cms.VPSet(
            cms.PSet(
                BPMinImprovement = cms.double(10000.0),
                BrutePruning = cms.bool(True),
                CSCDebug = cms.untracked.bool(False),
                CorrectTheErrors = cms.bool(True),
                Covariance = cms.double(0.0),
                ForceCovariance = cms.bool(False),
                ForceCovarianceAll = cms.bool(False),
                NormChi2Cut2D = cms.double(20.0),
                NormChi2Cut3D = cms.double(10.0),
                Pruning = cms.bool(True),
                SeedBig = cms.double(0.0015),
                SeedSmall = cms.double(0.0002),
                curvePenalty = cms.double(2.0),
                curvePenaltyThreshold = cms.double(0.85),
                dPhiFineMax = cms.double(0.025),
                dRPhiFineMax = cms.double(8.0),
                dXclusBoxMax = cms.double(4.0),
                dYclusBoxMax = cms.double(8.0),
                hitDropLimit4Hits = cms.double(0.6),
                hitDropLimit5Hits = cms.double(0.8),
                hitDropLimit6Hits = cms.double(0.3333),
                maxDPhi = cms.double(999.0),
                maxDTheta = cms.double(999.0),
                maxRatioResidualPrune = cms.double(3.0),
                maxRecHitsInCluster = cms.int32(20),
                minHitsPerSegment = cms.int32(3),
                onlyBestSegment = cms.bool(False),
                preClustering = cms.bool(True),
                preClusteringUseChaining = cms.bool(True),
                prePrun = cms.bool(True),
                prePrunLimit = cms.double(3.17),
                tanPhiMax = cms.double(0.5),
                tanThetaMax = cms.double(1.2),
                useShowering = cms.bool(False),
                yweightPenalty = cms.double(1.5),
                yweightPenaltyThreshold = cms.double(1.0)
            ), 
            cms.PSet(
                BPMinImprovement = cms.double(10000.0),
                BrutePruning = cms.bool(True),
                CSCDebug = cms.untracked.bool(False),
                CorrectTheErrors = cms.bool(True),
                Covariance = cms.double(0.0),
                ForceCovariance = cms.bool(False),
                ForceCovarianceAll = cms.bool(False),
                NormChi2Cut2D = cms.double(20.0),
                NormChi2Cut3D = cms.double(10.0),
                Pruning = cms.bool(True),
                SeedBig = cms.double(0.0015),
                SeedSmall = cms.double(0.0002),
                curvePenalty = cms.double(2.0),
                curvePenaltyThreshold = cms.double(0.85),
                dPhiFineMax = cms.double(0.025),
                dRPhiFineMax = cms.double(8.0),
                dXclusBoxMax = cms.double(4.0),
                dYclusBoxMax = cms.double(8.0),
                hitDropLimit4Hits = cms.double(0.6),
                hitDropLimit5Hits = cms.double(0.8),
                hitDropLimit6Hits = cms.double(0.3333),
                maxDPhi = cms.double(999.0),
                maxDTheta = cms.double(999.0),
                maxRatioResidualPrune = cms.double(3.0),
                maxRecHitsInCluster = cms.int32(24),
                minHitsPerSegment = cms.int32(3),
                onlyBestSegment = cms.bool(False),
                preClustering = cms.bool(True),
                preClusteringUseChaining = cms.bool(True),
                prePrun = cms.bool(True),
                prePrunLimit = cms.double(3.17),
                tanPhiMax = cms.double(0.5),
                tanThetaMax = cms.double(1.2),
                useShowering = cms.bool(False),
                yweightPenalty = cms.double(1.5),
                yweightPenaltyThreshold = cms.double(1.0)
            )
        ),
        chamber_types = cms.vstring(
            'ME1/a', 
            'ME1/b', 
            'ME1/2', 
            'ME1/3', 
            'ME2/1', 
            'ME2/2', 
            'ME3/1', 
            'ME3/2', 
            'ME4/1', 
            'ME4/2'
        ),
        parameters_per_chamber_type = cms.vint32(
            2, 1, 1, 1, 1, 
            1, 1, 1, 1, 1
        )
    )),
    algo_type = cms.int32(1),
    inputObjects = cms.InputTag("hltCsc2DRecHits")
)


process.hltDeepBLifetimeTagInfosPF = cms.EDProducer("CandIPProducer",
    candidates = cms.InputTag("hltParticleFlow"),
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jets = cms.InputTag("hltPFJetForBtag"),
    maxDeltaR = cms.double(0.4),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    minimumNumberOfPixelHits = cms.int32(2),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("hltVerticesPFFilter"),
    useTrackQuality = cms.bool(False)
)


process.hltDeepCombinedSecondaryVertexBJetTagsCalo = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepFlavourNoSL.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsInfosCalo"),
    toAdd = cms.PSet(
        probbb = cms.string('probb'),
        probcc = cms.string('probc')
    )
)


process.hltDeepCombinedSecondaryVertexBJetTagsInfos = cms.EDProducer("DeepNNTagInfoProducer",
    computer = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5.0),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500.0),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3.0),
            min_pT = cms.double(120.0),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(3),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5.0),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500.0),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3.0),
            min_pT = cms.double(120.0),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(2),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(3),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    svTagInfos = cms.InputTag("hltDeepSecondaryVertexTagInfosPF")
)


process.hltDeepCombinedSecondaryVertexBJetTagsInfosCalo = cms.EDProducer("TrackDeepNNTagInfoProducer",
    computer = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5.0),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500.0),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3.0),
            min_pT = cms.double(120.0),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5.0),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500.0),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3.0),
            min_pT = cms.double(120.0),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    svTagInfos = cms.InputTag("hltInclusiveSecondaryVertexFinderTagInfos")
)


process.hltDeepCombinedSecondaryVertexBJetTagsPF = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepFlavourNoSL.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsInfos"),
    toAdd = cms.PSet(
        probbb = cms.string('probb'),
        probcc = cms.string('probc')
    )
)


process.hltDeepInclusiveMergedVerticesPF = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.0),
    secondaryVertices = cms.InputTag("hltDeepTrackVertexArbitratorPF")
)


process.hltDeepInclusiveSecondaryVerticesPF = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.7),
    minSignificance = cms.double(2.0),
    secondaryVertices = cms.InputTag("hltDeepInclusiveVertexFinderPF")
)


process.hltDeepInclusiveVertexFinderPF = cms.EDProducer("InclusiveCandidateVertexFinder",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterizer = cms.PSet(
        clusterMaxDistance = cms.double(0.05),
        clusterMaxSignificance = cms.double(4.5),
        clusterMinAngleCosine = cms.double(0.5),
        distanceRatio = cms.double(20.0),
        seedMax3DIPSignificance = cms.double(9999.0),
        seedMax3DIPValue = cms.double(9999.0),
        seedMin3DIPSignificance = cms.double(1.2),
        seedMin3DIPValue = cms.double(0.005)
    ),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxNTracks = cms.uint32(30),
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3.0),
    minHits = cms.uint32(8),
    minPt = cms.double(0.8),
    primaryVertices = cms.InputTag("hltVerticesPFFilter"),
    tracks = cms.InputTag("hltParticleFlow"),
    useDirectVertexFitter = cms.bool(True),
    useVertexReco = cms.bool(True),
    vertexMinAngleCosine = cms.double(0.95),
    vertexMinDLen2DSig = cms.double(2.5),
    vertexMinDLenSig = cms.double(0.5),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        primcut = cms.double(1.0),
        seccut = cms.double(3.0),
        smoothing = cms.bool(True)
    )
)


process.hltDeepSecondaryVertexTagInfosPF = cms.EDProducer("CandSecondaryVertexProducer",
    beamSpotTag = cms.InputTag("hltOnlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("hltDeepInclusiveMergedVerticesPF"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("hltDeepBLifetimeTagInfosPF"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500.0),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3.0),
        min_pT = cms.double(120.0),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(2),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(3),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(True),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(3.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.65),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.4),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.hltDeepTrackVertexArbitratorPF = cms.EDProducer("CandidateVertexArbitrator",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxTimeSignificance = cms.double(3.5),
    primaryVertices = cms.InputTag("hltVerticesPFFilter"),
    secondaryVertices = cms.InputTag("hltDeepInclusiveSecondaryVerticesPF"),
    sigCut = cms.double(5.0),
    trackMinLayers = cms.int32(4),
    trackMinPixels = cms.int32(1),
    trackMinPt = cms.double(0.4),
    tracks = cms.InputTag("hltParticleFlow")
)


process.hltDoubletRecoveryClustersRefRemoval = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(16.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("hltTripletRecoveryClustersRefRemoval"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClusters"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltTripletRecoveryPFlowTrackSelectionHighPurity")
)


process.hltDoubletRecoveryMaskedMeasurementTrackerEvent = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltDoubletRecoveryClustersRefRemoval"),
    src = cms.InputTag("hltSiStripClusters")
)


process.hltDoubletRecoveryPFlowCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltDoubletRecoveryMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter2GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltDoubletRecoveryPFlowPixelSeeds"),
    useHitsSplitting = cms.bool(False)
)


process.hltDoubletRecoveryPFlowCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltDoubletRecovery'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltDoubletRecoveryMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltDoubletRecoveryPFlowCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltDoubletRecoveryPFlowPixelClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(50000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClusters"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltDoubletRecoveryPFlowPixelHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltDoubletRecoveryPFlowPixelClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag("hltDoubletRecoveryPixelLayerPairs"),
    trackingRegions = cms.InputTag("hltDoubletRecoveryPFlowPixelTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltDoubletRecoveryPFlowPixelSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("hltDoubletRecoveryPFlowPixelHitDoublets")
)


process.hltDoubletRecoveryPFlowPixelTrackingRegions = cms.EDProducer("PointSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(0.9),
        deltaPhi = cms.double(0.3),
        maxNRegions = cms.int32(100),
        maxNVertices = cms.int32(10),
        measurementTrackerName = cms.InputTag("hltDoubletRecoveryMaskedMeasurementTrackerEvent"),
        mode = cms.string('VerticesFixed'),
        nSigmaZBeamSpot = cms.double(3.0),
        nSigmaZVertex = cms.double(3.0),
        originRadius = cms.double(0.02),
        points = cms.PSet(
            eta = cms.vdouble(0.8),
            phi = cms.vdouble(3.0)
        ),
        precise = cms.bool(True),
        ptMin = cms.double(0.8),
        searchOpt = cms.bool(False),
        vertexCollection = cms.InputTag("hltTrimmedPixelVertices"),
        whereToUseMeasurementTracker = cms.string('ForSiStrips'),
        zErrorBeamSpot = cms.double(15.0),
        zErrorVetex = cms.double(0.1)
    )
)


process.hltDoubletRecoveryPFlowTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dr_par2 = cms.vdouble(3.40282346639e+38, 0.3, 0.3)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dz_par2 = cms.vdouble(3.40282346639e+38, 0.35, 0.35)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltDoubletRecoveryPFlowCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltTrimmedPixelVertices")
)


process.hltDoubletRecoveryPFlowTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltDoubletRecoveryPFlowTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltDoubletRecoveryPFlowTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltDoubletRecoveryPFlowCtfWithMaterialTracks")
)


process.hltDoubletRecoveryPixelLayerPairs = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltDoubletRecoveryClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltDoubletRecoveryClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix4', 
        'BPix1+BPix2', 
        'BPix1+BPix3', 
        'BPix2+BPix3', 
        'BPix1+FPix1_pos', 
        'BPix1+FPix1_pos', 
        'BPix1+FPix1_pos', 
        'BPix2+FPix1_neg', 
        'BPix3+BPix4', 
        'BPix3+FPix1_pos', 
        'FPix1_pos+FPix2_pos'
    )
)


process.hltDt1DRecHits = cms.EDProducer("DTRecHitProducer",
    debug = cms.untracked.bool(False),
    dtDigiLabel = cms.InputTag("hltMuonDTDigis"),
    recAlgo = cms.string('DTLinearDriftFromDBAlgo'),
    recAlgoConfig = cms.PSet(
        debug = cms.untracked.bool(False),
        doVdriftCorr = cms.bool(True),
        maxTime = cms.double(420.0),
        minTime = cms.double(-3.0),
        stepTwoFromDigi = cms.bool(False),
        tTrigMode = cms.string('DTTTrigSyncFromDB'),
        tTrigModeConfig = cms.PSet(
            debug = cms.untracked.bool(False),
            doT0Correction = cms.bool(True),
            doTOFCorrection = cms.bool(True),
            doWirePropCorrection = cms.bool(True),
            tTrigLabel = cms.string(''),
            tofCorrType = cms.int32(0),
            vPropWire = cms.double(24.4),
            wirePropCorrType = cms.int32(0)
        ),
        useUncertDB = cms.bool(True)
    )
)


process.hltDt4DSegments = cms.EDProducer("DTRecSegment4DProducer",
    Reco4DAlgoConfig = cms.PSet(
        AllDTRecHits = cms.bool(True),
        Reco2DAlgoConfig = cms.PSet(
            AlphaMaxPhi = cms.double(1.0),
            AlphaMaxTheta = cms.double(0.9),
            MaxAllowedHits = cms.uint32(50),
            debug = cms.untracked.bool(False),
            hit_afterT0_resolution = cms.double(0.03),
            nSharedHitsMax = cms.int32(2),
            nUnSharedHitsMin = cms.int32(2),
            performT0SegCorrection = cms.bool(False),
            performT0_vdriftSegCorrection = cms.bool(False),
            perform_delta_rejecting = cms.bool(False),
            recAlgo = cms.string('DTLinearDriftFromDBAlgo'),
            recAlgoConfig = cms.PSet(
                debug = cms.untracked.bool(False),
                doVdriftCorr = cms.bool(True),
                maxTime = cms.double(420.0),
                minTime = cms.double(-3.0),
                stepTwoFromDigi = cms.bool(False),
                tTrigMode = cms.string('DTTTrigSyncFromDB'),
                tTrigModeConfig = cms.PSet(
                    debug = cms.untracked.bool(False),
                    doT0Correction = cms.bool(True),
                    doTOFCorrection = cms.bool(True),
                    doWirePropCorrection = cms.bool(True),
                    tTrigLabel = cms.string(''),
                    tofCorrType = cms.int32(0),
                    vPropWire = cms.double(24.4),
                    wirePropCorrType = cms.int32(0)
                ),
                useUncertDB = cms.bool(True)
            ),
            segmCleanerMode = cms.int32(2)
        ),
        Reco2DAlgoName = cms.string('DTCombinatorialPatternReco'),
        debug = cms.untracked.bool(False),
        hit_afterT0_resolution = cms.double(0.03),
        nSharedHitsMax = cms.int32(2),
        nUnSharedHitsMin = cms.int32(2),
        performT0SegCorrection = cms.bool(False),
        performT0_vdriftSegCorrection = cms.bool(False),
        perform_delta_rejecting = cms.bool(False),
        recAlgo = cms.string('DTLinearDriftFromDBAlgo'),
        recAlgoConfig = cms.PSet(
            debug = cms.untracked.bool(False),
            doVdriftCorr = cms.bool(True),
            maxTime = cms.double(420.0),
            minTime = cms.double(-3.0),
            stepTwoFromDigi = cms.bool(False),
            tTrigMode = cms.string('DTTTrigSyncFromDB'),
            tTrigModeConfig = cms.PSet(
                debug = cms.untracked.bool(False),
                doT0Correction = cms.bool(True),
                doTOFCorrection = cms.bool(True),
                doWirePropCorrection = cms.bool(True),
                tTrigLabel = cms.string(''),
                tofCorrType = cms.int32(0),
                vPropWire = cms.double(24.4),
                wirePropCorrType = cms.int32(0)
            ),
            useUncertDB = cms.bool(True)
        ),
        segmCleanerMode = cms.int32(2)
    ),
    Reco4DAlgoName = cms.string('DTCombinatorialPatternReco4D'),
    debug = cms.untracked.bool(False),
    recHits1DLabel = cms.InputTag("hltDt1DRecHits"),
    recHits2DLabel = cms.InputTag("dt2DSegments")
)


process.hltEcalDetIdToBeRecovered = cms.EDProducer("EcalDetIdToBeRecoveredProducer",
    ebDetIdToBeRecovered = cms.string('ebDetId'),
    ebFEToBeRecovered = cms.string('ebFE'),
    ebIntegrityChIdErrors = cms.InputTag("hltEcalDigis","EcalIntegrityChIdErrors"),
    ebIntegrityGainErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainErrors"),
    ebIntegrityGainSwitchErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainSwitchErrors"),
    ebSrFlagCollection = cms.InputTag("hltEcalDigis"),
    eeDetIdToBeRecovered = cms.string('eeDetId'),
    eeFEToBeRecovered = cms.string('eeFE'),
    eeIntegrityChIdErrors = cms.InputTag("hltEcalDigis","EcalIntegrityChIdErrors"),
    eeIntegrityGainErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainErrors"),
    eeIntegrityGainSwitchErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainSwitchErrors"),
    eeSrFlagCollection = cms.InputTag("hltEcalDigis"),
    integrityBlockSizeErrors = cms.InputTag("hltEcalDigis","EcalIntegrityBlockSizeErrors"),
    integrityTTIdErrors = cms.InputTag("hltEcalDigis","EcalIntegrityTTIdErrors")
)


process.hltEcalDigis = cms.EDProducer("EcalRawToDigi",
    DoRegional = cms.bool(False),
    FEDs = cms.vint32(
        601, 602, 603, 604, 605, 
        606, 607, 608, 609, 610, 
        611, 612, 613, 614, 615, 
        616, 617, 618, 619, 620, 
        621, 622, 623, 624, 625, 
        626, 627, 628, 629, 630, 
        631, 632, 633, 634, 635, 
        636, 637, 638, 639, 640, 
        641, 642, 643, 644, 645, 
        646, 647, 648, 649, 650, 
        651, 652, 653, 654
    ),
    FedLabel = cms.InputTag("listfeds"),
    InputLabel = cms.InputTag("rawDataCollector"),
    eventPut = cms.bool(True),
    feIdCheck = cms.bool(True),
    feUnpacking = cms.bool(True),
    forceToKeepFRData = cms.bool(False),
    headerUnpacking = cms.bool(True),
    memUnpacking = cms.bool(True),
    numbTriggerTSamples = cms.int32(1),
    numbXtalTSamples = cms.int32(10),
    orderedDCCIdList = cms.vint32(
        1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        11, 12, 13, 14, 15, 
        16, 17, 18, 19, 20, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        31, 32, 33, 34, 35, 
        36, 37, 38, 39, 40, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        51, 52, 53, 54
    ),
    orderedFedList = cms.vint32(
        601, 602, 603, 604, 605, 
        606, 607, 608, 609, 610, 
        611, 612, 613, 614, 615, 
        616, 617, 618, 619, 620, 
        621, 622, 623, 624, 625, 
        626, 627, 628, 629, 630, 
        631, 632, 633, 634, 635, 
        636, 637, 638, 639, 640, 
        641, 642, 643, 644, 645, 
        646, 647, 648, 649, 650, 
        651, 652, 653, 654
    ),
    silentMode = cms.untracked.bool(True),
    srpUnpacking = cms.bool(True),
    syncCheck = cms.bool(True),
    tccUnpacking = cms.bool(True)
)


process.hltEcalPreshowerDigis = cms.EDProducer("ESRawToDigi",
    ESdigiCollection = cms.string(''),
    InstanceES = cms.string(''),
    LookupTable = cms.FileInPath('EventFilter/ESDigiToRaw/data/ES_lookup_table.dat'),
    debugMode = cms.untracked.bool(False),
    sourceTag = cms.InputTag("rawDataCollector")
)


process.hltEcalPreshowerRecHit = cms.EDProducer("ESRecHitProducer",
    ESRecoAlgo = cms.int32(0),
    ESdigiCollection = cms.InputTag("hltEcalPreshowerDigis"),
    ESrechitCollection = cms.string('EcalRecHitsES'),
    algo = cms.string('ESRecHitWorker')
)


process.hltEcalRecHit = cms.EDProducer("EcalRecHitProducer",
    ChannelStatusToBeExcluded = cms.vstring(),
    EBLaserMAX = cms.double(3.0),
    EBLaserMIN = cms.double(0.5),
    EBrechitCollection = cms.string('EcalRecHitsEB'),
    EBuncalibRecHitCollection = cms.InputTag("hltEcalUncalibRecHit","EcalUncalibRecHitsEB"),
    EELaserMAX = cms.double(8.0),
    EELaserMIN = cms.double(0.5),
    EErechitCollection = cms.string('EcalRecHitsEE'),
    EEuncalibRecHitCollection = cms.InputTag("hltEcalUncalibRecHit","EcalUncalibRecHitsEE"),
    algo = cms.string('EcalRecHitWorkerSimple'),
    algoRecover = cms.string('EcalRecHitWorkerRecover'),
    cleaningConfig = cms.PSet(
        cThreshold_barrel = cms.double(4.0),
        cThreshold_double = cms.double(10.0),
        cThreshold_endcap = cms.double(15.0),
        e4e1Threshold_barrel = cms.double(0.08),
        e4e1Threshold_endcap = cms.double(0.3),
        e4e1_a_barrel = cms.double(0.04),
        e4e1_a_endcap = cms.double(0.02),
        e4e1_b_barrel = cms.double(-0.024),
        e4e1_b_endcap = cms.double(-0.0125),
        e6e2thresh = cms.double(0.04),
        ignoreOutOfTimeThresh = cms.double(1000000000.0),
        tightenCrack_e1_double = cms.double(2.0),
        tightenCrack_e1_single = cms.double(2.0),
        tightenCrack_e4e1_single = cms.double(3.0),
        tightenCrack_e6e2_double = cms.double(3.0)
    ),
    dbStatusToBeExcludedEB = cms.vint32(14, 78, 142),
    dbStatusToBeExcludedEE = cms.vint32(14, 78, 142),
    ebDetIdToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","ebDetId"),
    ebFEToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","ebFE"),
    eeDetIdToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","eeDetId"),
    eeFEToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","eeFE"),
    flagsMapDBReco = cms.PSet(
        kDead = cms.vstring('kNoDataNoTP'),
        kGood = cms.vstring(
            'kOk', 
            'kDAC', 
            'kNoLaser', 
            'kNoisy'
        ),
        kNeighboursRecovered = cms.vstring(
            'kFixedG0', 
            'kNonRespondingIsolated', 
            'kDeadVFE'
        ),
        kNoisy = cms.vstring(
            'kNNoisy', 
            'kFixedG6', 
            'kFixedG1'
        ),
        kTowerRecovered = cms.vstring('kDeadFE')
    ),
    killDeadChannels = cms.bool(True),
    laserCorrection = cms.bool(True),
    logWarningEtThreshold_EB_FE = cms.double(50.0),
    logWarningEtThreshold_EE_FE = cms.double(50.0),
    recoverEBFE = cms.bool(True),
    recoverEBIsolatedChannels = cms.bool(False),
    recoverEBVFE = cms.bool(False),
    recoverEEFE = cms.bool(True),
    recoverEEIsolatedChannels = cms.bool(False),
    recoverEEVFE = cms.bool(False),
    singleChannelRecoveryMethod = cms.string('NeuralNetworks'),
    singleChannelRecoveryThreshold = cms.double(8.0),
    skipTimeCalib = cms.bool(True),
    triggerPrimitiveDigiCollection = cms.InputTag("hltEcalDigis","EcalTriggerPrimitives")
)


process.hltEcalUncalibRecHit = cms.EDProducer("EcalUncalibRecHitProducer",
    EBdigiCollection = cms.InputTag("hltEcalDigis","ebDigis"),
    EBhitCollection = cms.string('EcalUncalibRecHitsEB'),
    EEdigiCollection = cms.InputTag("hltEcalDigis","eeDigis"),
    EEhitCollection = cms.string('EcalUncalibRecHitsEE'),
    algo = cms.string('EcalUncalibRecHitWorkerMultiFit'),
    algoPSet = cms.PSet(
        EBamplitudeFitParameters = cms.vdouble(1.138, 1.652),
        EBtimeConstantTerm = cms.double(0.6),
        EBtimeFitLimits_Lower = cms.double(0.2),
        EBtimeFitLimits_Upper = cms.double(1.4),
        EBtimeFitParameters = cms.vdouble(
            -2.015452, 3.130702, -12.3473, 41.88921, -82.83944, 
            91.01147, -50.35761, 11.05621
        ),
        EBtimeNconst = cms.double(28.5),
        EEamplitudeFitParameters = cms.vdouble(1.89, 1.4),
        EEtimeConstantTerm = cms.double(1.0),
        EEtimeFitLimits_Lower = cms.double(0.2),
        EEtimeFitLimits_Upper = cms.double(1.4),
        EEtimeFitParameters = cms.vdouble(
            -2.390548, 3.553628, -17.62341, 67.67538, -133.213, 
            140.7432, -75.41106, 16.20277
        ),
        EEtimeNconst = cms.double(31.8),
        EcalPulseShapeParameters = cms.PSet(
            EBCorrNoiseMatrixG01 = cms.vdouble(
                1.0, 0.73354, 0.64442, 0.58851, 0.55425, 
                0.53082, 0.51916, 0.51097, 0.50732, 0.50409
            ),
            EBCorrNoiseMatrixG06 = cms.vdouble(
                1.0, 0.70946, 0.58021, 0.49846, 0.45006, 
                0.41366, 0.39699, 0.38478, 0.37847, 0.37055
            ),
            EBCorrNoiseMatrixG12 = cms.vdouble(
                1.0, 0.71073, 0.55721, 0.46089, 0.40449, 
                0.35931, 0.33924, 0.32439, 0.31581, 0.30481
            ),
            EBPulseShapeCovariance = cms.vdouble(
                3.001e-06, 1.233e-05, 0.0, -4.416e-06, -4.571e-06, 
                -3.614e-06, -2.636e-06, -1.286e-06, -8.41e-07, -5.296e-07, 
                0.0, 0.0, 1.233e-05, 6.154e-05, 0.0, 
                -2.2e-05, -2.309e-05, -1.838e-05, -1.373e-05, -7.334e-06, 
                -5.088e-06, -3.745e-06, -2.428e-06, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, -4.416e-06, -2.2e-05, 0.0, 8.319e-06, 
                8.545e-06, 6.792e-06, 5.059e-06, 2.678e-06, 1.816e-06, 
                1.223e-06, 8.245e-07, 5.589e-07, -4.571e-06, -2.309e-05, 
                0.0, 8.545e-06, 9.182e-06, 7.219e-06, 5.388e-06, 
                2.853e-06, 1.944e-06, 1.324e-06, 9.083e-07, 6.335e-07, 
                -3.614e-06, -1.838e-05, 0.0, 6.792e-06, 7.219e-06, 
                6.016e-06, 4.437e-06, 2.385e-06, 1.636e-06, 1.118e-06, 
                7.754e-07, 5.556e-07, -2.636e-06, -1.373e-05, 0.0, 
                5.059e-06, 5.388e-06, 4.437e-06, 3.602e-06, 1.917e-06, 
                1.322e-06, 9.079e-07, 6.529e-07, 4.752e-07, -1.286e-06, 
                -7.334e-06, 0.0, 2.678e-06, 2.853e-06, 2.385e-06, 
                1.917e-06, 1.375e-06, 9.1e-07, 6.455e-07, 4.693e-07, 
                3.657e-07, -8.41e-07, -5.088e-06, 0.0, 1.816e-06, 
                1.944e-06, 1.636e-06, 1.322e-06, 9.1e-07, 9.115e-07, 
                6.062e-07, 4.436e-07, 3.422e-07, -5.296e-07, -3.745e-06, 
                0.0, 1.223e-06, 1.324e-06, 1.118e-06, 9.079e-07, 
                6.455e-07, 6.062e-07, 7.217e-07, 4.862e-07, 3.768e-07, 
                0.0, -2.428e-06, 0.0, 8.245e-07, 9.083e-07, 
                7.754e-07, 6.529e-07, 4.693e-07, 4.436e-07, 4.862e-07, 
                6.509e-07, 4.418e-07, 0.0, 0.0, 0.0, 
                5.589e-07, 6.335e-07, 5.556e-07, 4.752e-07, 3.657e-07, 
                3.422e-07, 3.768e-07, 4.418e-07, 6.142e-07
            ),
            EBPulseShapeTemplate = cms.vdouble(
                0.0113979, 0.758151, 1.0, 0.887744, 0.673548, 
                0.474332, 0.319561, 0.215144, 0.147464, 0.101087, 
                0.0693181, 0.0475044
            ),
            EBdigiCollection = cms.string(''),
            EECorrNoiseMatrixG01 = cms.vdouble(
                1.0, 0.72698, 0.62048, 0.55691, 0.51848, 
                0.49147, 0.47813, 0.47007, 0.46621, 0.46265
            ),
            EECorrNoiseMatrixG06 = cms.vdouble(
                1.0, 0.71217, 0.47464, 0.34056, 0.26282, 
                0.20287, 0.17734, 0.16256, 0.15618, 0.14443
            ),
            EECorrNoiseMatrixG12 = cms.vdouble(
                1.0, 0.71373, 0.44825, 0.30152, 0.21609, 
                0.14786, 0.11772, 0.10165, 0.09465, 0.08098
            ),
            EEPulseShapeCovariance = cms.vdouble(
                3.941e-05, 3.333e-05, 0.0, -1.449e-05, -1.661e-05, 
                -1.424e-05, -1.183e-05, -6.842e-06, -4.915e-06, -3.411e-06, 
                0.0, 0.0, 3.333e-05, 2.862e-05, 0.0, 
                -1.244e-05, -1.431e-05, -1.233e-05, -1.032e-05, -5.883e-06, 
                -4.154e-06, -2.902e-06, -2.128e-06, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, -1.449e-05, -1.244e-05, 0.0, 5.84e-06, 
                6.649e-06, 5.72e-06, 4.812e-06, 2.708e-06, 1.869e-06, 
                1.33e-06, 9.186e-07, 6.446e-07, -1.661e-05, -1.431e-05, 
                0.0, 6.649e-06, 7.966e-06, 6.898e-06, 5.794e-06, 
                3.157e-06, 2.184e-06, 1.567e-06, 1.084e-06, 7.575e-07, 
                -1.424e-05, -1.233e-05, 0.0, 5.72e-06, 6.898e-06, 
                6.341e-06, 5.347e-06, 2.859e-06, 1.991e-06, 1.431e-06, 
                9.839e-07, 6.886e-07, -1.183e-05, -1.032e-05, 0.0, 
                4.812e-06, 5.794e-06, 5.347e-06, 4.854e-06, 2.628e-06, 
                1.809e-06, 1.289e-06, 9.02e-07, 6.146e-07, -6.842e-06, 
                -5.883e-06, 0.0, 2.708e-06, 3.157e-06, 2.859e-06, 
                2.628e-06, 1.863e-06, 1.296e-06, 8.882e-07, 6.108e-07, 
                4.283e-07, -4.915e-06, -4.154e-06, 0.0, 1.869e-06, 
                2.184e-06, 1.991e-06, 1.809e-06, 1.296e-06, 1.217e-06, 
                8.669e-07, 5.751e-07, 3.882e-07, -3.411e-06, -2.902e-06, 
                0.0, 1.33e-06, 1.567e-06, 1.431e-06, 1.289e-06, 
                8.882e-07, 8.669e-07, 9.522e-07, 6.717e-07, 4.293e-07, 
                0.0, -2.128e-06, 0.0, 9.186e-07, 1.084e-06, 
                9.839e-07, 9.02e-07, 6.108e-07, 5.751e-07, 6.717e-07, 
                7.911e-07, 5.493e-07, 0.0, 0.0, 0.0, 
                6.446e-07, 7.575e-07, 6.886e-07, 6.146e-07, 4.283e-07, 
                3.882e-07, 4.293e-07, 5.493e-07, 7.027e-07
            ),
            EEPulseShapeTemplate = cms.vdouble(
                0.116442, 0.756246, 1.0, 0.897182, 0.686831, 
                0.491506, 0.344111, 0.245731, 0.174115, 0.123361, 
                0.0874288, 0.061957
            ),
            EEdigiCollection = cms.string(''),
            ESdigiCollection = cms.string(''),
            EcalPreMixStage1 = cms.bool(False),
            EcalPreMixStage2 = cms.bool(False),
            UseLCcorrection = cms.untracked.bool(True)
        ),
        activeBXs = cms.vint32(
            -5, -4, -3, -2, -1, 
            0, 1, 2, 3, 4
        ),
        addPedestalUncertaintyEB = cms.double(0.0),
        addPedestalUncertaintyEE = cms.double(0.0),
        ampErrorCalculation = cms.bool(False),
        amplitudeThresholdEB = cms.double(10.0),
        amplitudeThresholdEE = cms.double(10.0),
        chi2ThreshEB_ = cms.double(65.0),
        chi2ThreshEE_ = cms.double(50.0),
        doPrefitEB = cms.bool(False),
        doPrefitEE = cms.bool(False),
        dynamicPedestalsEB = cms.bool(False),
        dynamicPedestalsEE = cms.bool(False),
        ebPulseShape = cms.vdouble(
            5.2e-05, -5.26e-05, 6.66e-05, 0.1168, 0.7575, 
            1.0, 0.8876, 0.6732, 0.4741, 0.3194
        ),
        ebSpikeThreshold = cms.double(1.042),
        eePulseShape = cms.vdouble(
            5.2e-05, -5.26e-05, 6.66e-05, 0.1168, 0.7575, 
            1.0, 0.8876, 0.6732, 0.4741, 0.3194
        ),
        gainSwitchUseMaxSampleEB = cms.bool(True),
        gainSwitchUseMaxSampleEE = cms.bool(False),
        kPoorRecoFlagEB = cms.bool(True),
        kPoorRecoFlagEE = cms.bool(False),
        mitigateBadSamplesEB = cms.bool(False),
        mitigateBadSamplesEE = cms.bool(False),
        outOfTimeThresholdGain12mEB = cms.double(5.0),
        outOfTimeThresholdGain12mEE = cms.double(1000.0),
        outOfTimeThresholdGain12pEB = cms.double(5.0),
        outOfTimeThresholdGain12pEE = cms.double(1000.0),
        outOfTimeThresholdGain61mEB = cms.double(5.0),
        outOfTimeThresholdGain61mEE = cms.double(1000.0),
        outOfTimeThresholdGain61pEB = cms.double(5.0),
        outOfTimeThresholdGain61pEE = cms.double(1000.0),
        prefitMaxChiSqEB = cms.double(25.0),
        prefitMaxChiSqEE = cms.double(10.0),
        selectiveBadSampleCriteriaEB = cms.bool(False),
        selectiveBadSampleCriteriaEE = cms.bool(False),
        simplifiedNoiseModelForGainSwitch = cms.bool(True),
        timealgo = cms.string('None'),
        useLumiInfoRunHeader = cms.bool(False)
    )
)


process.hltFEDSelector = cms.EDProducer("EvFFEDSelector",
    fedList = cms.vuint32(1023, 1024),
    inputTag = cms.InputTag("rawDataCollector")
)


process.hltFastPVJetTracksAssociator = cms.EDProducer("JetTracksAssociatorAtVertex",
    coneSize = cms.double(0.4),
    jets = cms.InputTag("hltSelector4CentralJetsL1FastJet"),
    pvSrc = cms.InputTag(""),
    tracks = cms.InputTag("hltFastPVPixelTracks"),
    useAssigned = cms.bool(False)
)


process.hltFastPVPixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltFastPVPixelTracksFilter"),
    Fitter = cms.InputTag("hltFastPVPixelTracksFitter"),
    SeedingHitSets = cms.InputTag("hltFastPVPixelTracksHitQuadruplets"),
    passLabel = cms.string('')
)


process.hltFastPVPixelTracksFilter = cms.EDProducer("PixelTrackFilterByKinematicsProducer",
    chi2 = cms.double(1000.0),
    nSigmaInvPtTolerance = cms.double(0.0),
    nSigmaTipMaxTolerance = cms.double(0.0),
    ptMin = cms.double(0.0),
    tipMax = cms.double(999.0)
)


process.hltFastPVPixelTracksFitter = cms.EDProducer("PixelFitterByHelixProjectionsProducer",
    scaleErrorsForBPix1 = cms.bool(False),
    scaleFactor = cms.double(0.65)
)


process.hltFastPVPixelTracksHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag(""),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltPixelLayerQuadrupletsRegForBTag"),
    trackingRegions = cms.InputTag("hltFastPVPixelTracksTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltFastPVPixelTracksHitDoubletsRecover = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag(""),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltPixelLayerQuadrupletsRegForBTag"),
    trackingRegions = cms.InputTag("hltFastPVPixelTracksTrackingRegionsRecover"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltFastPVPixelTracksHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.2),
    CAThetaCut = cms.double(0.002),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersRegForBTagCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltFastPVPixelTracksHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2.0),
        value1 = cms.double(200.0),
        value2 = cms.double(50.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltFastPVPixelTracksHitQuadrupletsRecover = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.2),
    CAThetaCut = cms.double(0.002),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersRegForBTagCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltFastPVPixelTracksHitDoubletsRecover"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(False),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2.0),
        value1 = cms.double(200.0),
        value2 = cms.double(50.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltFastPVPixelTracksMerger = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltFastPVPixelTracks", "hltFastPVPixelTracksRecover"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(False),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltFastPVPixelTracks", "hltFastPVPixelTracksRecover"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltFastPVPixelTracksRecover = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltFastPVPixelTracksRecoverFilter"),
    Fitter = cms.InputTag("hltFastPVPixelTracksRecoverFitter"),
    SeedingHitSets = cms.InputTag("hltFastPVPixelTracksHitQuadrupletsRecover"),
    passLabel = cms.string('')
)


process.hltFastPVPixelTracksRecoverFilter = cms.EDProducer("PixelTrackFilterByKinematicsProducer",
    chi2 = cms.double(1000.0),
    nSigmaInvPtTolerance = cms.double(0.0),
    nSigmaTipMaxTolerance = cms.double(0.0),
    ptMin = cms.double(0.0),
    tipMax = cms.double(999.0)
)


process.hltFastPVPixelTracksRecoverFitter = cms.EDProducer("PixelFitterByHelixProjectionsProducer",
    scaleErrorsForBPix1 = cms.bool(False),
    scaleFactor = cms.double(0.65)
)


process.hltFastPVPixelTracksTrackingRegions = cms.EDProducer("TauRegionalPixelSeedTrackingRegionEDProducer",
    RegionPSet = cms.PSet(
        JetSrc = cms.InputTag("hltSelector4CentralJetsL1FastJet"),
        deltaEtaRegion = cms.double(0.5),
        deltaPhiRegion = cms.double(0.3),
        howToUseMeasurementTracker = cms.string('Never'),
        measurementTrackerName = cms.InputTag("MeasurementTrackerEvent"),
        originHalfLength = cms.double(1.5),
        originRadius = cms.double(0.2),
        ptMin = cms.double(0.9),
        searchOpt = cms.bool(False),
        vertexSrc = cms.InputTag("hltFastPrimaryVertex")
    )
)


process.hltFastPVPixelTracksTrackingRegionsRecover = cms.EDProducer("TauRegionalPixelSeedTrackingRegionEDProducer",
    RegionPSet = cms.PSet(
        JetSrc = cms.InputTag("hltFastPVJetVertexChecker"),
        deltaEtaRegion = cms.double(0.5),
        deltaPhiRegion = cms.double(0.5),
        howToUseMeasurementTracker = cms.string('Never'),
        measurementTrackerName = cms.InputTag("MeasurementTrackerEvent"),
        originHalfLength = cms.double(20.0),
        originRadius = cms.double(0.2),
        ptMin = cms.double(0.9),
        searchOpt = cms.bool(False),
        vertexSrc = cms.InputTag("hltFastPVJetVertexChecker")
    )
)


process.hltFastPVPixelVertices = cms.EDProducer("PixelVertexProducer",
    Finder = cms.string('DivisiveVertexFinder'),
    Method2 = cms.bool(True),
    NTrkMin = cms.int32(2),
    PVcomparer = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPvClusterComparerForBTag')
    ),
    PtMin = cms.double(1.0),
    TrackCollection = cms.InputTag("hltFastPVPixelTracksMerger"),
    UseError = cms.bool(True),
    Verbosity = cms.int32(0),
    WtAverage = cms.bool(True),
    ZOffset = cms.double(10.0),
    ZSeparation = cms.double(0.07),
    beamSpot = cms.InputTag("hltOnlineBeamSpot")
)


process.hltFastPixelBLifetimeL3Associator = cms.EDProducer("JetTracksAssociatorAtVertex",
    coneSize = cms.double(0.4),
    jets = cms.InputTag("hltSelector8CentralJetsL1FastJet"),
    pvSrc = cms.InputTag(""),
    tracks = cms.InputTag("hltIter2MergedForBTag"),
    useAssigned = cms.bool(False)
)


process.hltFastPrimaryVertex = cms.EDProducer("FastPrimaryVertexWithWeightsProducer",
    EC_weight = cms.double(0.008),
    PixelCellHeightOverWidth = cms.double(1.8),
    barrel = cms.bool(True),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusters = cms.InputTag("hltSiPixelClustersRegForBTag"),
    endCap = cms.bool(True),
    jets = cms.InputTag("hltSelector4CentralJetsL1FastJet"),
    maxDeltaPhi = cms.double(0.21),
    maxDeltaPhi_EC = cms.double(0.14),
    maxJetEta = cms.double(2.6),
    maxJetEta_EC = cms.double(2.6),
    maxSizeX = cms.double(2.1),
    maxSizeY_q = cms.double(2.0),
    maxZ = cms.double(19.0),
    minJetEta_EC = cms.double(1.3),
    minJetPt = cms.double(0.0),
    minSizeY_q = cms.double(-0.6),
    njets = cms.int32(999),
    peakSizeY_q = cms.double(1.0),
    pixelCPE = cms.string('hltESPPixelCPEGeneric'),
    ptWeighting = cms.bool(True),
    ptWeighting_offset = cms.double(-1.0),
    ptWeighting_slope = cms.double(0.05),
    weightCut_step2 = cms.double(0.05),
    weightCut_step3 = cms.double(0.1),
    weight_SizeX1 = cms.double(0.88),
    weight_charge_down = cms.double(11000.0),
    weight_charge_peak = cms.double(22000.0),
    weight_charge_up = cms.double(190000.0),
    weight_dPhi = cms.double(0.13888888),
    weight_dPhi_EC = cms.double(0.064516129),
    weight_rho_up = cms.double(22.0),
    zClusterSearchArea_step2 = cms.double(3.0),
    zClusterSearchArea_step3 = cms.double(0.55),
    zClusterWidth_step1 = cms.double(2.0),
    zClusterWidth_step2 = cms.double(0.65),
    zClusterWidth_step3 = cms.double(0.3)
)


process.hltFixedGridRhoFastjetAll = cms.EDProducer("FixedGridRhoProducerFastjet",
    gridSpacing = cms.double(0.55),
    maxRapidity = cms.double(5.0),
    pfCandidatesTag = cms.InputTag("hltParticleFlow")
)


process.hltFixedGridRhoFastjetAllCalo = cms.EDProducer("FixedGridRhoProducerFastjet",
    gridSpacing = cms.double(0.55),
    maxRapidity = cms.double(5.0),
    pfCandidatesTag = cms.InputTag("hltTowerMakerForAll")
)


process.hltGtStage2Digis = cms.EDProducer("L1TRawToDigi",
    CTP7 = cms.untracked.bool(False),
    FWId = cms.uint32(0),
    FWOverride = cms.bool(False),
    FedIds = cms.vint32(1404),
    InputLabel = cms.InputTag("rawDataCollector"),
    MTF7 = cms.untracked.bool(False),
    MinFeds = cms.uint32(0),
    Setup = cms.string('stage2::GTSetup'),
    TMTCheck = cms.bool(True),
    debug = cms.untracked.bool(False),
    lenAMC13Header = cms.untracked.int32(8),
    lenAMC13Trailer = cms.untracked.int32(8),
    lenAMCHeader = cms.untracked.int32(8),
    lenAMCTrailer = cms.untracked.int32(0),
    lenSlinkHeader = cms.untracked.int32(8),
    lenSlinkTrailer = cms.untracked.int32(8)
)


process.hltGtStage2ObjectMap = cms.EDProducer("L1TGlobalProducer",
    AlgorithmTriggersUnmasked = cms.bool(True),
    AlgorithmTriggersUnprescaled = cms.bool(True),
    AlternativeNrBxBoardDaq = cms.uint32(0),
    BstLengthBytes = cms.int32(-1),
    EGammaInputTag = cms.InputTag("hltGtStage2Digis","EGamma"),
    EmulateBxInEvent = cms.int32(1),
    EtSumInputTag = cms.InputTag("hltGtStage2Digis","EtSum"),
    ExtInputTag = cms.InputTag("hltGtStage2Digis"),
    JetInputTag = cms.InputTag("hltGtStage2Digis","Jet"),
    L1DataBxInEvent = cms.int32(5),
    MuonInputTag = cms.InputTag("hltGtStage2Digis","Muon"),
    PrescaleCSVFile = cms.string('prescale_L1TGlobal.csv'),
    PrescaleSet = cms.uint32(1),
    PrintL1Menu = cms.untracked.bool(False),
    ProduceL1GtDaqRecord = cms.bool(True),
    ProduceL1GtObjectMapRecord = cms.bool(True),
    TauInputTag = cms.InputTag("hltGtStage2Digis","Tau"),
    TriggerMenuLuminosity = cms.string('startup'),
    Verbosity = cms.untracked.int32(0)
)


process.hltHbhePhase1Reco = cms.EDProducer("HBHEPhase1Reconstructor",
    algoConfigClass = cms.string(''),
    algorithm = cms.PSet(
        Class = cms.string('SimpleHBHEPhase1Algo'),
        activeBXs = cms.vint32(-1, 0, 1),
        applyPedConstraint = cms.bool(True),
        applyPulseJitter = cms.bool(False),
        applyTimeConstraint = cms.bool(True),
        applyTimeSlew = cms.bool(True),
        applyTimeSlewM3 = cms.bool(True),
        chiSqSwitch = cms.double(-10.0),
        correctForPhaseContainment = cms.bool(True),
        correctionPhaseNS = cms.double(6.0),
        deltaChiSqThresh = cms.double(0.001),
        dynamicPed = cms.bool(True),
        firstSampleShift = cms.int32(0),
        fitTimes = cms.int32(1),
        meanPed = cms.double(0.0),
        meanTime = cms.double(0.0),
        nMaxItersMin = cms.int32(500),
        nMaxItersNNLS = cms.int32(500),
        nnlsThresh = cms.double(1e-11),
        noiseHPD = cms.double(1.0),
        noiseSiPM = cms.double(1.0),
        pedSigmaHPD = cms.double(0.5),
        pedSigmaSiPM = cms.double(0.00065),
        pedestalUpperLimit = cms.double(2.7),
        pulseJitter = cms.double(1.0),
        respCorrM3 = cms.double(1.0),
        samplesToAdd = cms.int32(2),
        tdcTimeShift = cms.double(0.0),
        timeMax = cms.double(12.5),
        timeMin = cms.double(-12.5),
        timeSigmaHPD = cms.double(5.0),
        timeSigmaSiPM = cms.double(2.5),
        timeSlewPars = cms.vdouble(
            12.2999, -2.19142, 0.0, 12.2999, -2.19142, 
            0.0, 12.2999, -2.19142, 0.0
        ),
        timeSlewParsType = cms.int32(3),
        ts4Max = cms.vdouble(100.0, 45000.0),
        ts4Min = cms.double(0.0),
        ts4Thresh = cms.double(0.0),
        ts4chi2 = cms.vdouble(15.0, 15.0),
        useM2 = cms.bool(False),
        useM3 = cms.bool(True),
        useMahi = cms.bool(False)
    ),
    digiLabelQIE11 = cms.InputTag("hltHcalDigis"),
    digiLabelQIE8 = cms.InputTag("hltHcalDigis"),
    dropZSmarkedPassed = cms.bool(True),
    flagParametersQIE11 = cms.PSet(

    ),
    flagParametersQIE8 = cms.PSet(
        hitEnergyMinimum = cms.double(1.0),
        hitMultiplicityThreshold = cms.int32(17),
        nominalPedestal = cms.double(3.0),
        pulseShapeParameterSets = cms.VPSet(
            cms.PSet(
                pulseShapeParameters = cms.vdouble(
                    0.0, 100.0, -50.0, 0.0, -15.0, 
                    0.15
                )
            ), 
            cms.PSet(
                pulseShapeParameters = cms.vdouble(
                    100.0, 2000.0, -50.0, 0.0, -5.0, 
                    0.05
                )
            ), 
            cms.PSet(
                pulseShapeParameters = cms.vdouble(
                    2000.0, 1000000.0, -50.0, 0.0, 95.0, 
                    0.0
                )
            ), 
            cms.PSet(
                pulseShapeParameters = cms.vdouble(
                    -1000000.0, 1000000.0, 45.0, 0.1, 1000000.0, 
                    0.0
                )
            )
        )
    ),
    makeRecHits = cms.bool(True),
    processQIE11 = cms.bool(True),
    processQIE8 = cms.bool(True),
    pulseShapeParametersQIE11 = cms.PSet(

    ),
    pulseShapeParametersQIE8 = cms.PSet(
        LeftSlopeCut = cms.vdouble(5.0, 2.55, 2.55),
        LeftSlopeThreshold = cms.vdouble(250.0, 500.0, 100000.0),
        LinearCut = cms.vdouble(-3.0, -0.054, -0.054),
        LinearThreshold = cms.vdouble(20.0, 100.0, 100000.0),
        MinimumChargeThreshold = cms.double(20.0),
        MinimumTS4TS5Threshold = cms.double(100.0),
        R45MinusOneRange = cms.double(0.2),
        R45PlusOneRange = cms.double(0.2),
        RMS8MaxCut = cms.vdouble(-13.5, -11.5, -11.5),
        RMS8MaxThreshold = cms.vdouble(20.0, 100.0, 100000.0),
        RightSlopeCut = cms.vdouble(5.0, 4.15, 4.15),
        RightSlopeSmallCut = cms.vdouble(1.08, 1.16, 1.16),
        RightSlopeSmallThreshold = cms.vdouble(150.0, 200.0, 100000.0),
        RightSlopeThreshold = cms.vdouble(250.0, 400.0, 100000.0),
        TS3TS4ChargeThreshold = cms.double(70.0),
        TS3TS4UpperChargeThreshold = cms.double(20.0),
        TS4TS5ChargeThreshold = cms.double(70.0),
        TS4TS5LowerCut = cms.vdouble(
            -1.0, -0.7, -0.5, -0.4, -0.3, 
            0.1
        ),
        TS4TS5LowerThreshold = cms.vdouble(
            100.0, 120.0, 160.0, 200.0, 300.0, 
            500.0
        ),
        TS4TS5UpperCut = cms.vdouble(1.0, 0.8, 0.75, 0.72),
        TS4TS5UpperThreshold = cms.vdouble(70.0, 90.0, 100.0, 400.0),
        TS5TS6ChargeThreshold = cms.double(70.0),
        TS5TS6UpperChargeThreshold = cms.double(20.0),
        TriangleIgnoreSlow = cms.bool(False),
        TrianglePeakTS = cms.uint32(10000),
        UseDualFit = cms.bool(True)
    ),
    recoParamsFromDB = cms.bool(True),
    saveDroppedInfos = cms.bool(False),
    saveEffectivePedestal = cms.bool(False),
    saveInfos = cms.bool(False),
    setLegacyFlagsQIE11 = cms.bool(False),
    setLegacyFlagsQIE8 = cms.bool(False),
    setNegativeFlagsQIE11 = cms.bool(False),
    setNegativeFlagsQIE8 = cms.bool(False),
    setNoiseFlagsQIE11 = cms.bool(False),
    setNoiseFlagsQIE8 = cms.bool(True),
    setPulseShapeFlagsQIE11 = cms.bool(False),
    setPulseShapeFlagsQIE8 = cms.bool(True),
    sipmQNTStoSum = cms.int32(3),
    sipmQTSShift = cms.int32(0),
    tsFromDB = cms.bool(False)
)


process.hltHbhereco = cms.EDProducer("HBHEPlan1Combiner",
    algorithm = cms.PSet(
        Class = cms.string('SimplePlan1RechitCombiner')
    ),
    hbheInput = cms.InputTag("hltHbhePhase1Reco"),
    ignorePlan1Topology = cms.bool(False),
    usePlan1Mode = cms.bool(True)
)


process.hltHcalDigis = cms.EDProducer("HcalRawToDigi",
    ComplainEmptyData = cms.untracked.bool(False),
    ElectronicsMap = cms.string(''),
    ExpectedOrbitMessageTime = cms.untracked.int32(-1),
    FEDs = cms.untracked.vint32(),
    FilterDataQuality = cms.bool(True),
    HcalFirstFED = cms.untracked.int32(700),
    InputLabel = cms.InputTag("rawDataCollector"),
    UnpackCalib = cms.untracked.bool(True),
    UnpackTTP = cms.untracked.bool(False),
    UnpackUMNio = cms.untracked.bool(True),
    UnpackZDC = cms.untracked.bool(True),
    UnpackerMode = cms.untracked.int32(0),
    firstSample = cms.int32(0),
    lastSample = cms.int32(9),
    silent = cms.untracked.bool(True)
)


process.hltHfprereco = cms.EDProducer("HFPreReconstructor",
    digiLabel = cms.InputTag("hltHcalDigis"),
    dropZSmarkedPassed = cms.bool(True),
    forceSOI = cms.int32(-1),
    soiShift = cms.int32(0),
    sumAllTimeSlices = cms.bool(False),
    tsFromDB = cms.bool(False)
)


process.hltHfreco = cms.EDProducer("HFPhase1Reconstructor",
    PETstat = cms.PSet(
        HcalAcceptSeverityLevel = cms.int32(9),
        longETParams = cms.vdouble(
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0
        ),
        longEnergyParams = cms.vdouble(
            43.5, 45.7, 48.32, 51.36, 54.82, 
            58.7, 63.0, 67.72, 72.86, 78.42, 
            84.4, 90.8, 97.62
        ),
        long_R = cms.vdouble(0.98),
        long_R_29 = cms.vdouble(0.8),
        shortETParams = cms.vdouble(
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0
        ),
        shortEnergyParams = cms.vdouble(
            35.1773, 35.37, 35.7933, 36.4472, 37.3317, 
            38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 
            47.4813, 49.98, 52.7093
        ),
        short_R = cms.vdouble(0.8),
        short_R_29 = cms.vdouble(0.8)
    ),
    S8S1stat = cms.PSet(
        HcalAcceptSeverityLevel = cms.int32(9),
        isS8S1 = cms.bool(True),
        longETParams = cms.vdouble(
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0
        ),
        longEnergyParams = cms.vdouble(
            40.0, 100.0, 100.0, 100.0, 100.0, 
            100.0, 100.0, 100.0, 100.0, 100.0, 
            100.0, 100.0, 100.0
        ),
        long_optimumSlope = cms.vdouble(
            0.3, 0.1, 0.1, 0.1, 0.1, 
            0.1, 0.1, 0.1, 0.1, 0.1, 
            0.1, 0.1, 0.1
        ),
        shortETParams = cms.vdouble(
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0
        ),
        shortEnergyParams = cms.vdouble(
            40.0, 100.0, 100.0, 100.0, 100.0, 
            100.0, 100.0, 100.0, 100.0, 100.0, 
            100.0, 100.0, 100.0
        ),
        short_optimumSlope = cms.vdouble(
            0.3, 0.1, 0.1, 0.1, 0.1, 
            0.1, 0.1, 0.1, 0.1, 0.1, 
            0.1, 0.1, 0.1
        )
    ),
    S9S1stat = cms.PSet(
        HcalAcceptSeverityLevel = cms.int32(9),
        isS8S1 = cms.bool(False),
        longETParams = cms.vdouble(
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0
        ),
        longEnergyParams = cms.vdouble(
            43.5, 45.7, 48.32, 51.36, 54.82, 
            58.7, 63.0, 67.72, 72.86, 78.42, 
            84.4, 90.8, 97.62
        ),
        long_optimumSlope = cms.vdouble(
            -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 
            0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 
            0.135313, 0.136289, 0.0589927
        ),
        shortETParams = cms.vdouble(
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0
        ),
        shortEnergyParams = cms.vdouble(
            35.1773, 35.37, 35.7933, 36.4472, 37.3317, 
            38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 
            47.4813, 49.98, 52.7093
        ),
        short_optimumSlope = cms.vdouble(
            -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 
            0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 
            0.135313, 0.136289, 0.0589927
        )
    ),
    algoConfigClass = cms.string('HFPhase1PMTParams'),
    algorithm = cms.PSet(
        Class = cms.string('HFFlexibleTimeCheck'),
        energyWeights = cms.vdouble(
            1.0, 1.0, 1.0, 0.0, 1.0, 
            0.0, 2.0, 0.0, 2.0, 0.0, 
            2.0, 0.0, 1.0, 0.0, 0.0, 
            1.0, 0.0, 1.0, 0.0, 2.0, 
            0.0, 2.0, 0.0, 2.0, 0.0, 
            1.0
        ),
        rejectAllFailures = cms.bool(True),
        soiPhase = cms.uint32(1),
        tfallIfNoTDC = cms.double(-101.0),
        timeShift = cms.double(0.0),
        tlimits = cms.vdouble(-1000.0, 1000.0, -1000.0, 1000.0),
        triseIfNoTDC = cms.double(-100.0)
    ),
    checkChannelQualityForDepth3and4 = cms.bool(False),
    inputLabel = cms.InputTag("hltHfprereco"),
    setNoiseFlags = cms.bool(True),
    useChannelQualityFromDB = cms.bool(False)
)


process.hltHoreco = cms.EDProducer("HcalHitReconstructor",
    HFInWindowStat = cms.PSet(

    ),
    PETstat = cms.PSet(

    ),
    S8S1stat = cms.PSet(

    ),
    S9S1stat = cms.PSet(

    ),
    Subdetector = cms.string('HO'),
    correctForPhaseContainment = cms.bool(True),
    correctForTimeslew = cms.bool(True),
    correctTiming = cms.bool(False),
    correctionPhaseNS = cms.double(13.0),
    dataOOTCorrectionCategory = cms.string('Data'),
    dataOOTCorrectionName = cms.string(''),
    digiLabel = cms.InputTag("hltHcalDigis"),
    digiTimeFromDB = cms.bool(True),
    digistat = cms.PSet(

    ),
    dropZSmarkedPassed = cms.bool(True),
    firstAuxTS = cms.int32(4),
    firstSample = cms.int32(4),
    hfTimingTrustParameters = cms.PSet(

    ),
    mcOOTCorrectionCategory = cms.string('MC'),
    mcOOTCorrectionName = cms.string(''),
    recoParamsFromDB = cms.bool(True),
    samplesToAdd = cms.int32(4),
    saturationParameters = cms.PSet(
        maxADCvalue = cms.int32(127)
    ),
    setHSCPFlags = cms.bool(False),
    setNegativeFlags = cms.bool(False),
    setNoiseFlags = cms.bool(False),
    setPulseShapeFlags = cms.bool(False),
    setSaturationFlags = cms.bool(False),
    setTimingTrustFlags = cms.bool(False),
    tsFromDB = cms.bool(True),
    useLeakCorrection = cms.bool(False)
)


process.hltHtMhtCaloJetsSixC25 = cms.EDProducer("HLTHtMhtProducer",
    excludePFMuons = cms.bool(False),
    jetsLabel = cms.InputTag("hltCaloJetsSix25ForHt"),
    maxEtaJetHt = cms.double(2.6),
    maxEtaJetMht = cms.double(5.0),
    minNJetHt = cms.int32(6),
    minNJetMht = cms.int32(0),
    minPtJetHt = cms.double(25.0),
    minPtJetMht = cms.double(0.0),
    pfCandidatesLabel = cms.InputTag(""),
    usePt = cms.bool(True)
)


process.hltHtMhtPFJetsSixC32 = cms.EDProducer("HLTHtMhtProducer",
    excludePFMuons = cms.bool(False),
    jetsLabel = cms.InputTag("hltPFJetsSix32ForHt"),
    maxEtaJetHt = cms.double(2.6),
    maxEtaJetMht = cms.double(999.0),
    minNJetHt = cms.int32(6),
    minNJetMht = cms.int32(0),
    minPtJetHt = cms.double(32.0),
    minPtJetMht = cms.double(0.0),
    pfCandidatesLabel = cms.InputTag("hltParticleFlow"),
    usePt = cms.bool(True)
)


process.hltImpactParameterTagInfos = cms.EDProducer("TrackIPProducer",
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jetTracks = cms.InputTag("hltFastPixelBLifetimeL3Associator"),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    minimumNumberOfPixelHits = cms.int32(2),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("hltVerticesL3","WithBS"),
    useTrackQuality = cms.bool(False)
)


process.hltInclusiveMergedVertices = cms.EDProducer("VertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.0),
    secondaryVertices = cms.InputTag("hltTrackVertexArbitrator")
)


process.hltInclusiveMergedVerticesPF = cms.EDProducer("VertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.0),
    secondaryVertices = cms.InputTag("hltTrackVertexArbitratorPF")
)


process.hltInclusiveSecondaryVertexFinderTagInfos = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("hltOnlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("hltInclusiveMergedVertices"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("hltImpactParameterTagInfos"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500.0),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3.0),
        min_pT = cms.double(120.0),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(2),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(2),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(True),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(2.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.79),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.5),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.hltInclusiveSecondaryVertices = cms.EDProducer("VertexMerger",
    maxFraction = cms.double(0.7),
    minSignificance = cms.double(2.0),
    secondaryVertices = cms.InputTag("hltInclusiveVertexFinder")
)


process.hltInclusiveSecondaryVerticesPF = cms.EDProducer("VertexMerger",
    maxFraction = cms.double(0.7),
    minSignificance = cms.double(2.0),
    secondaryVertices = cms.InputTag("hltInclusiveVertexFinderPF")
)


process.hltInclusiveVertexFinder = cms.EDProducer("InclusiveVertexFinder",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterizer = cms.PSet(
        clusterMaxDistance = cms.double(0.05),
        clusterMaxSignificance = cms.double(4.5),
        clusterMinAngleCosine = cms.double(0.5),
        distanceRatio = cms.double(20.0),
        seedMax3DIPSignificance = cms.double(9999.0),
        seedMax3DIPValue = cms.double(9999.0),
        seedMin3DIPSignificance = cms.double(1.2),
        seedMin3DIPValue = cms.double(0.005)
    ),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxNTracks = cms.uint32(30),
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3.0),
    minHits = cms.uint32(8),
    minPt = cms.double(0.8),
    primaryVertices = cms.InputTag("hltVerticesL3"),
    tracks = cms.InputTag("hltIter2MergedForBTag"),
    useDirectVertexFitter = cms.bool(True),
    useVertexReco = cms.bool(True),
    vertexMinAngleCosine = cms.double(0.95),
    vertexMinDLen2DSig = cms.double(2.5),
    vertexMinDLenSig = cms.double(0.5),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        primcut = cms.double(1.0),
        seccut = cms.double(3.0),
        smoothing = cms.bool(True)
    )
)


process.hltInclusiveVertexFinderPF = cms.EDProducer("InclusiveVertexFinder",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterizer = cms.PSet(
        clusterMaxDistance = cms.double(0.05),
        clusterMaxSignificance = cms.double(4.5),
        clusterMinAngleCosine = cms.double(0.5),
        distanceRatio = cms.double(20.0),
        seedMax3DIPSignificance = cms.double(9999.0),
        seedMax3DIPValue = cms.double(9999.0),
        seedMin3DIPSignificance = cms.double(1.2),
        seedMin3DIPValue = cms.double(0.005)
    ),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxNTracks = cms.uint32(30),
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3.0),
    minHits = cms.uint32(8),
    minPt = cms.double(0.8),
    primaryVertices = cms.InputTag("hltVerticesPFFilter"),
    tracks = cms.InputTag("hltPFMuonMerging"),
    useDirectVertexFitter = cms.bool(True),
    useVertexReco = cms.bool(True),
    vertexMinAngleCosine = cms.double(0.95),
    vertexMinDLen2DSig = cms.double(2.5),
    vertexMinDLenSig = cms.double(0.5),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        primcut = cms.double(1.0),
        seccut = cms.double(3.0),
        smoothing = cms.bool(True)
    )
)


process.hltIter0IterL3FromL1MuonCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('none'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string('GroupedCkfTrajectoryBuilder'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter0IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks"),
    useHitsSplitting = cms.bool(True)
)


process.hltIter0IterL3FromL1MuonCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter0'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter0IterL3FromL1MuonCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltIterL3FromL1MuonPixelTracks"),
    InputVertexCollection = cms.InputTag("hltIterL3FromL1MuonTrimmedPixelVertices"),
    SeedCreatorPSet = cms.PSet(
        refToPSet_ = cms.string('HLTSeedFromProtoTracks')
    ),
    TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
    originHalfLength = cms.double(0.3),
    originRadius = cms.double(0.1),
    useEventsWithNoVertex = cms.bool(True),
    usePV = cms.bool(False),
    useProtoTrackKinematics = cms.bool(False)
)


process.hltIter0IterL3FromL1MuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
            d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
            dr_exp = cms.vint32(4, 4, 2147483647),
            dr_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dr_par2 = cms.vdouble(0.3, 0.3, 3.40282346639e+38)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 2147483647),
            dz_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dz_par2 = cms.vdouble(0.35, 0.35, 3.40282346639e+38)
        ),
        maxChi2 = cms.vdouble(3.40282346639e+38, 3.40282346639e+38, 3.40282346639e+38),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 3, 4),
        minLayers = cms.vint32(3, 3, 4),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 3, 4)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter0IterL3FromL1MuonCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltIterL3FromL1MuonTrimmedPixelVertices")
)


process.hltIter0IterL3FromL1MuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter0IterL3FromL1MuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter0IterL3FromL1MuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter0IterL3FromL1MuonCtfWithMaterialTracks")
)


process.hltIter0IterL3MuonCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('none'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string('GroupedCkfTrajectoryBuilder'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter0IterL3MuonPSetGroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks"),
    useHitsSplitting = cms.bool(True)
)


process.hltIter0IterL3MuonCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter0'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter0IterL3MuonCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter0IterL3MuonPixelSeedsFromPixelTracks = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltIterL3MuonPixelTracks"),
    InputVertexCollection = cms.InputTag("hltIterL3MuonTrimmedPixelVertices"),
    SeedCreatorPSet = cms.PSet(
        refToPSet_ = cms.string('HLTSeedFromProtoTracks')
    ),
    TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
    originHalfLength = cms.double(0.3),
    originRadius = cms.double(0.1),
    useEventsWithNoVertex = cms.bool(True),
    usePV = cms.bool(False),
    useProtoTrackKinematics = cms.bool(False)
)


process.hltIter0IterL3MuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
            d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
            dr_exp = cms.vint32(4, 4, 2147483647),
            dr_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dr_par2 = cms.vdouble(0.3, 0.3, 3.40282346639e+38)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 2147483647),
            dz_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dz_par2 = cms.vdouble(0.35, 0.35, 3.40282346639e+38)
        ),
        maxChi2 = cms.vdouble(3.40282346639e+38, 3.40282346639e+38, 3.40282346639e+38),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 3, 4),
        minLayers = cms.vint32(3, 3, 4),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 3, 4)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter0IterL3MuonCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltIterL3MuonTrimmedPixelVertices")
)


process.hltIter0IterL3MuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter0IterL3MuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter0IterL3MuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter0IterL3MuonCtfWithMaterialTracks")
)


process.hltIter0PFLowPixelSeedsFromPixelTracks = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltPixelTracks"),
    InputVertexCollection = cms.InputTag("hltTrimmedPixelVertices"),
    SeedCreatorPSet = cms.PSet(
        refToPSet_ = cms.string('HLTSeedFromProtoTracks')
    ),
    TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
    originHalfLength = cms.double(0.3),
    originRadius = cms.double(0.1),
    useEventsWithNoVertex = cms.bool(True),
    usePV = cms.bool(False),
    useProtoTrackKinematics = cms.bool(False)
)


process.hltIter0PFLowPixelSeedsFromPixelTracksForBTag = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltFastPVPixelTracksMerger"),
    InputVertexCollection = cms.InputTag("hltFastPVPixelVertices"),
    SeedCreatorPSet = cms.PSet(
        refToPSet_ = cms.string('HLTSeedFromProtoTracks')
    ),
    TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
    originHalfLength = cms.double(0.3),
    originRadius = cms.double(0.1),
    useEventsWithNoVertex = cms.bool(True),
    usePV = cms.bool(False),
    useProtoTrackKinematics = cms.bool(False)
)


process.hltIter0PFlowCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter0GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter0PFLowPixelSeedsFromPixelTracks"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter0PFlowCkfTrackCandidatesForBTag = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClustersRegForBTag"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter0GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter0PFLowPixelSeedsFromPixelTracksForBTag"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter0PFlowCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter0'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter0PFlowCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter0PFlowCtfWithMaterialTracksForBTag = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter0'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClustersRegForBTag"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter0PFlowCkfTrackCandidatesForBTag"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter0PFlowTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(0.4, 0.4, 0.4),
            dr_par2 = cms.vdouble(0.3, 0.3, 0.3)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(0.4, 0.4, 0.4),
            dz_par2 = cms.vdouble(0.35, 0.35, 0.35)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 3, 4),
        minLayers = cms.vint32(3, 3, 4),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 3, 4)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter0PFlowCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltTrimmedPixelVertices")
)


process.hltIter0PFlowTrackCutClassifierForBTag = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(0.4, 0.4, 0.4),
            dr_par2 = cms.vdouble(0.3, 0.3, 0.3)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(0.4, 0.4, 0.4),
            dz_par2 = cms.vdouble(0.35, 0.35, 0.35)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 3, 4),
        minLayers = cms.vint32(3, 3, 4),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 3, 4)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter0PFlowCtfWithMaterialTracksForBTag"),
    vertices = cms.InputTag("hltFastPVPixelVertices")
)


process.hltIter0PFlowTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter0PFlowTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter0PFlowTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter0PFlowCtfWithMaterialTracks")
)


process.hltIter0PFlowTrackSelectionHighPurityForBTag = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter0PFlowTrackCutClassifierForBTag","MVAValues"),
    originalQualVals = cms.InputTag("hltIter0PFlowTrackCutClassifierForBTag","QualityMasks"),
    originalSource = cms.InputTag("hltIter0PFlowCtfWithMaterialTracksForBTag")
)


process.hltIter1ClustersRefRemoval = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClusters"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter0PFlowTrackSelectionHighPurity")
)


process.hltIter1ClustersRefRemovalForBTag = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClustersRegForBTag"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter0PFlowTrackSelectionHighPurityForBTag")
)


process.hltIter1MaskedMeasurementTrackerEvent = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltIter1ClustersRefRemoval"),
    src = cms.InputTag("hltSiStripClusters")
)


process.hltIter1MaskedMeasurementTrackerEventForBTag = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltIter1ClustersRefRemovalForBTag"),
    src = cms.InputTag("hltSiStripClustersRegForBTag")
)


process.hltIter1Merged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIter0PFlowTrackSelectionHighPurity", "hltIter1PFlowTrackSelectionHighPurity"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIter0PFlowTrackSelectionHighPurity", "hltIter1PFlowTrackSelectionHighPurity"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIter1MergedForBTag = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIter0PFlowTrackSelectionHighPurityForBTag", "hltIter1PFlowTrackSelectionHighPurityForBTag"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIter0PFlowTrackSelectionHighPurityForBTag", "hltIter1PFlowTrackSelectionHighPurityForBTag"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIter1PFLowPixelSeedsFromPixelTracks = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltIter1PixelTracks"),
    InputVertexCollection = cms.InputTag("hltTrimmedPixelVertices"),
    SeedCreatorPSet = cms.PSet(
        refToPSet_ = cms.string('HLTSeedFromProtoTracks')
    ),
    TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
    originHalfLength = cms.double(0.3),
    originRadius = cms.double(0.1),
    useEventsWithNoVertex = cms.bool(True),
    usePV = cms.bool(False),
    useProtoTrackKinematics = cms.bool(False)
)


process.hltIter1PFlowCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltIter1MaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter1GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter1PFLowPixelSeedsFromPixelTracks"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter1PFlowCkfTrackCandidatesForBTag = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltIter1MaskedMeasurementTrackerEventForBTag"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter1GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter1PFlowPixelSeedsForBTag"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter1PFlowCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter1'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltIter1MaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter1PFlowCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter1PFlowCtfWithMaterialTracksForBTag = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter1'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltIter1MaskedMeasurementTrackerEventForBTag"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter1PFlowCkfTrackCandidatesForBTag"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter1PFlowPixelClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(50000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClusters"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltIter1PFlowPixelClusterCheckForBTag = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClustersRegForBTag"),
    MaxNumberOfCosmicClusters = cms.uint32(800000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClustersRegForBTag"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltIter1PFlowPixelHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltIter1PFlowPixelClusterCheck"),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIter1PixelLayerQuadruplets"),
    trackingRegions = cms.InputTag("hltIter1PFlowPixelTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIter1PFlowPixelHitDoubletsForBTag = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltIter1PFlowPixelClusterCheckForBTag"),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIter1PixelLayerQuadrupletsForBTag"),
    trackingRegions = cms.InputTag("hltIter1PFlowPixelTrackingRegionsForBTag"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIter1PFlowPixelHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.3),
    CAThetaCut = cms.double(0.004),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltIter1PFlowPixelHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.4),
        pt2 = cms.double(2.0),
        value1 = cms.double(1000.0),
        value2 = cms.double(100.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIter1PFlowPixelHitQuadrupletsForBTag = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.3),
    CAThetaCut = cms.double(0.004),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersRegForBTagCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltIter1PFlowPixelHitDoubletsForBTag"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2.0),
        value1 = cms.double(2000.0),
        value2 = cms.double(150.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIter1PFlowPixelSeedsForBTag = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("hltIter1PFlowPixelHitQuadrupletsForBTag")
)


process.hltIter1PFlowPixelTrackingRegions = cms.EDProducer("GlobalTrackingRegionWithVerticesEDProducer",
    RegionPSet = cms.PSet(
        VertexCollection = cms.InputTag("hltTrimmedPixelVertices"),
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        fixedError = cms.double(0.2),
        nSigmaZ = cms.double(4.0),
        originRadius = cms.double(0.05),
        precise = cms.bool(True),
        ptMin = cms.double(0.4),
        sigmaZVertex = cms.double(3.0),
        useFakeVertices = cms.bool(False),
        useFixedError = cms.bool(True),
        useFoundVertices = cms.bool(True),
        useMultipleScattering = cms.bool(False)
    )
)


process.hltIter1PFlowPixelTrackingRegionsForBTag = cms.EDProducer("CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(1.0),
        deltaPhi = cms.double(1.0),
        input = cms.InputTag("hltSelector8CentralJetsL1FastJet"),
        maxNRegions = cms.int32(100),
        maxNVertices = cms.int32(10),
        measurementTrackerName = cms.InputTag("hltIter1MaskedMeasurementTrackerEventForBTag"),
        mode = cms.string('VerticesFixed'),
        nSigmaZBeamSpot = cms.double(4.0),
        nSigmaZVertex = cms.double(4.0),
        originRadius = cms.double(0.05),
        precise = cms.bool(True),
        ptMin = cms.double(0.3),
        searchOpt = cms.bool(True),
        vertexCollection = cms.InputTag("hltFastPVPixelVertices"),
        whereToUseMeasurementTracker = cms.string('ForSiStrips'),
        zErrorBeamSpot = cms.double(15.0),
        zErrorVetex = cms.double(0.1)
    )
)


process.hltIter1PFlowTrackCutClassifierDetached = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(1.0, 1.0, 1.0),
            dr_par2 = cms.vdouble(1.0, 1.0, 1.0)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(1.0, 1.0, 1.0),
            dz_par2 = cms.vdouble(1.0, 1.0, 1.0)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.0, 0.7, 0.4),
        maxDr = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDz = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(99, 3, 3),
        min3DLayers = cms.vint32(1, 2, 3),
        minLayers = cms.vint32(5, 5, 5),
        minNVtxTrk = cms.int32(2),
        minNdof = cms.vdouble(-1.0, -1.0, -1.0),
        minPixelHits = cms.vint32(0, 0, 1)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter1PFlowCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltTrimmedPixelVertices")
)


process.hltIter1PFlowTrackCutClassifierDetachedForBTag = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(1.0, 1.0, 1.0),
            dr_par2 = cms.vdouble(1.0, 1.0, 1.0)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(1.0, 1.0, 1.0),
            dz_par2 = cms.vdouble(1.0, 1.0, 1.0)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.0, 0.7, 0.4),
        maxDr = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDz = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(99, 3, 3),
        min3DLayers = cms.vint32(1, 2, 3),
        minLayers = cms.vint32(5, 5, 5),
        minNVtxTrk = cms.int32(2),
        minNdof = cms.vdouble(-1.0, -1.0, -1.0),
        minPixelHits = cms.vint32(0, 0, 1)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter1PFlowCtfWithMaterialTracksForBTag"),
    vertices = cms.InputTag("hltFastPVPixelVertices")
)


process.hltIter1PFlowTrackCutClassifierMerged = cms.EDProducer("ClassifierMerger",
    inputClassifiers = cms.vstring(
        'hltIter1PFlowTrackCutClassifierPrompt', 
        'hltIter1PFlowTrackCutClassifierDetached'
    )
)


process.hltIter1PFlowTrackCutClassifierMergedForBTag = cms.EDProducer("ClassifierMerger",
    inputClassifiers = cms.vstring(
        'hltIter1PFlowTrackCutClassifierPromptForBTag', 
        'hltIter1PFlowTrackCutClassifierDetachedForBTag'
    )
)


process.hltIter1PFlowTrackCutClassifierPrompt = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(3, 3, 3),
            dr_par1 = cms.vdouble(3.40282346639e+38, 1.0, 0.9),
            dr_par2 = cms.vdouble(3.40282346639e+38, 1.0, 0.85)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(3, 3, 3),
            dz_par1 = cms.vdouble(3.40282346639e+38, 1.0, 0.9),
            dz_par2 = cms.vdouble(3.40282346639e+38, 1.0, 0.8)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDz = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 2)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter1PFlowCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltTrimmedPixelVertices")
)


process.hltIter1PFlowTrackCutClassifierPromptForBTag = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(3, 3, 3),
            dr_par1 = cms.vdouble(3.40282346639e+38, 1.0, 0.9),
            dr_par2 = cms.vdouble(3.40282346639e+38, 1.0, 0.85)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(3, 3, 3),
            dz_par1 = cms.vdouble(3.40282346639e+38, 1.0, 0.9),
            dz_par2 = cms.vdouble(3.40282346639e+38, 1.0, 0.8)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDz = cms.vdouble(3.40282346639e+38, 1.0, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 2)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter1PFlowCtfWithMaterialTracksForBTag"),
    vertices = cms.InputTag("hltFastPVPixelVertices")
)


process.hltIter1PFlowTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter1PFlowTrackCutClassifierMerged","MVAValues"),
    originalQualVals = cms.InputTag("hltIter1PFlowTrackCutClassifierMerged","QualityMasks"),
    originalSource = cms.InputTag("hltIter1PFlowCtfWithMaterialTracks")
)


process.hltIter1PFlowTrackSelectionHighPurityForBTag = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter1PFlowTrackCutClassifierMergedForBTag","MVAValues"),
    originalQualVals = cms.InputTag("hltIter1PFlowTrackCutClassifierMergedForBTag","QualityMasks"),
    originalSource = cms.InputTag("hltIter1PFlowCtfWithMaterialTracksForBTag")
)


process.hltIter1PixelLayerQuadruplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltIter1ClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltIter1ClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4', 
        'BPix1+BPix2+BPix3+FPix1_pos', 
        'BPix1+BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos+FPix2_pos', 
        'BPix1+BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltIter1PixelLayerQuadrupletsForBTag = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHitsRegForBTag'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltIter1ClustersRefRemovalForBTag"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHitsRegForBTag'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltIter1ClustersRefRemovalForBTag"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4', 
        'BPix1+BPix2+BPix3+FPix1_pos', 
        'BPix1+BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos+FPix2_pos', 
        'BPix1+BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltIter1PixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltPixelTracksFilter"),
    Fitter = cms.InputTag("hltPixelTracksFitter"),
    SeedingHitSets = cms.InputTag("hltIter1PFlowPixelHitQuadruplets"),
    passLabel = cms.string('')
)


process.hltIter1TrackAndTauJets4Iter2 = cms.EDProducer("TauJetSelectorForHLTTrackSeeding",
    etaMaxCaloJet = cms.double(2.7),
    etaMinCaloJet = cms.double(-2.7),
    fractionMaxChargedPUInCaloCone = cms.double(0.3),
    fractionMinCaloInTauCone = cms.double(0.7),
    inputCaloJetTag = cms.InputTag("hltAK4CaloJetsPFEt5"),
    inputTrackJetTag = cms.InputTag("hltAK4Iter1TrackJets4Iter2"),
    inputTrackTag = cms.InputTag("hltIter1Merged"),
    isolationConeSize = cms.double(0.5),
    nTrkMaxInCaloCone = cms.int32(0),
    ptMinCaloJet = cms.double(10.0),
    ptTrkMaxInCaloCone = cms.double(1.4),
    tauConeSize = cms.double(0.2)
)


process.hltIter1TrackRefsForJets4Iter2 = cms.EDProducer("ChargedRefCandidateProducer",
    particleType = cms.string('pi+'),
    src = cms.InputTag("hltIter1Merged")
)


process.hltIter2ClustersRefRemoval = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(16.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("hltIter1ClustersRefRemoval"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClusters"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter1PFlowTrackSelectionHighPurity")
)


process.hltIter2ClustersRefRemovalForBTag = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(16.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("hltIter1ClustersRefRemovalForBTag"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClustersRegForBTag"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter1PFlowTrackSelectionHighPurityForBTag")
)


process.hltIter2IterL3FromL1MuonCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltIter2IterL3FromL1MuonMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter2IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter2IterL3FromL1MuonPixelSeeds"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter2IterL3FromL1MuonClustersRefRemoval = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(16.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClusters"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity")
)


process.hltIter2IterL3FromL1MuonCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter2'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltIter2IterL3FromL1MuonMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter2IterL3FromL1MuonCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter2IterL3FromL1MuonMaskedMeasurementTrackerEvent = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltIter2IterL3FromL1MuonClustersRefRemoval"),
    src = cms.InputTag("hltSiStripClusters")
)


process.hltIter2IterL3FromL1MuonMerged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "hltIter2IterL3FromL1MuonTrackSelectionHighPurity"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity", "hltIter2IterL3FromL1MuonTrackSelectionHighPurity"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIter2IterL3FromL1MuonPixelClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(50000),
    MaxNumberOfPixelClusters = cms.uint32(10000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClusters"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltIter2IterL3FromL1MuonPixelHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltIter2IterL3FromL1MuonPixelClusterCheck"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIter2IterL3FromL1MuonPixelLayerTriplets"),
    trackingRegions = cms.InputTag("hltIterL3FromL1MuonPixelTracksTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIter2IterL3FromL1MuonPixelHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.3),
    CAPhiCut = cms.double(0.1),
    CAThetaCut = cms.double(0.015),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("hltIter2IterL3FromL1MuonPixelHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(8.0),
        value1 = cms.double(100.0),
        value2 = cms.double(6.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIter2IterL3FromL1MuonPixelLayerTriplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltIter2IterL3FromL1MuonClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltIter2IterL3FromL1MuonClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3', 
        'BPix2+BPix3+BPix4', 
        'BPix1+BPix3+BPix4', 
        'BPix1+BPix2+BPix4', 
        'BPix2+BPix3+FPix1_pos', 
        'BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos', 
        'BPix1+BPix2+FPix1_neg', 
        'BPix2+FPix1_pos+FPix2_pos', 
        'BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos', 
        'BPix1+FPix1_neg+FPix2_neg', 
        'FPix1_pos+FPix2_pos+FPix3_pos', 
        'FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltIter2IterL3FromL1MuonPixelSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("hltIter2IterL3FromL1MuonPixelHitTriplets")
)


process.hltIter2IterL3FromL1MuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
            d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
            dr_exp = cms.vint32(4, 4, 2147483647),
            dr_par1 = cms.vdouble(3.40282346639e+38, 0.4, 3.40282346639e+38),
            dr_par2 = cms.vdouble(3.40282346639e+38, 0.3, 3.40282346639e+38)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 2147483647),
            dz_par1 = cms.vdouble(3.40282346639e+38, 0.4, 3.40282346639e+38),
            dz_par2 = cms.vdouble(3.40282346639e+38, 0.35, 3.40282346639e+38)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 3.40282346639e+38),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter2IterL3FromL1MuonCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltIterL3FromL1MuonTrimmedPixelVertices")
)


process.hltIter2IterL3FromL1MuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter2IterL3FromL1MuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter2IterL3FromL1MuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter2IterL3FromL1MuonCtfWithMaterialTracks")
)


process.hltIter2IterL3MuonCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltIter2IterL3MuonMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter2IterL3MuonPSetGroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter2IterL3MuonPixelSeeds"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter2IterL3MuonClustersRefRemoval = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(16.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClusters"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity")
)


process.hltIter2IterL3MuonCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter2'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltIter2IterL3MuonMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter2IterL3MuonCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter2IterL3MuonMaskedMeasurementTrackerEvent = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltIter2IterL3MuonClustersRefRemoval"),
    src = cms.InputTag("hltSiStripClusters")
)


process.hltIter2IterL3MuonMerged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIter0IterL3MuonTrackSelectionHighPurity", "hltIter2IterL3MuonTrackSelectionHighPurity"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIter0IterL3MuonTrackSelectionHighPurity", "hltIter2IterL3MuonTrackSelectionHighPurity"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIter2IterL3MuonPixelClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(50000),
    MaxNumberOfPixelClusters = cms.uint32(10000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClusters"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltIter2IterL3MuonPixelHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltIter2IterL3MuonPixelClusterCheck"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIter2IterL3MuonPixelLayerTriplets"),
    trackingRegions = cms.InputTag("hltIterL3MuonPixelTracksTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIter2IterL3MuonPixelHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.3),
    CAPhiCut = cms.double(0.1),
    CAThetaCut = cms.double(0.015),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("hltIter2IterL3MuonPixelHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(8.0),
        value1 = cms.double(100.0),
        value2 = cms.double(6.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIter2IterL3MuonPixelLayerTriplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltIter2IterL3MuonClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltIter2IterL3MuonClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3', 
        'BPix2+BPix3+BPix4', 
        'BPix1+BPix3+BPix4', 
        'BPix1+BPix2+BPix4', 
        'BPix2+BPix3+FPix1_pos', 
        'BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos', 
        'BPix1+BPix2+FPix1_neg', 
        'BPix2+FPix1_pos+FPix2_pos', 
        'BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos', 
        'BPix1+FPix1_neg+FPix2_neg', 
        'FPix1_pos+FPix2_pos+FPix3_pos', 
        'FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltIter2IterL3MuonPixelSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("hltIter2IterL3MuonPixelHitTriplets")
)


process.hltIter2IterL3MuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
            d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
            dr_exp = cms.vint32(4, 4, 2147483647),
            dr_par1 = cms.vdouble(3.40282346639e+38, 0.4, 3.40282346639e+38),
            dr_par2 = cms.vdouble(3.40282346639e+38, 0.3, 3.40282346639e+38)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 2147483647),
            dz_par1 = cms.vdouble(3.40282346639e+38, 0.4, 3.40282346639e+38),
            dz_par2 = cms.vdouble(3.40282346639e+38, 0.35, 3.40282346639e+38)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 3.40282346639e+38),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter2IterL3MuonCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltIterL3MuonTrimmedPixelVertices")
)


process.hltIter2IterL3MuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter2IterL3MuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter2IterL3MuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter2IterL3MuonCtfWithMaterialTracks")
)


process.hltIter2MaskedMeasurementTrackerEvent = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltIter2ClustersRefRemoval"),
    src = cms.InputTag("hltSiStripClusters")
)


process.hltIter2MaskedMeasurementTrackerEventForBTag = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltIter2ClustersRefRemovalForBTag"),
    src = cms.InputTag("hltSiStripClustersRegForBTag")
)


process.hltIter2Merged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIter1Merged", "hltIter2PFlowTrackSelectionHighPurity"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIter1Merged", "hltIter2PFlowTrackSelectionHighPurity"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIter2MergedForBTag = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIter1MergedForBTag", "hltIter2PFlowTrackSelectionHighPurityForBTag"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIter1MergedForBTag", "hltIter2PFlowTrackSelectionHighPurityForBTag"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIter2PFlowCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltIter2MaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter2GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter2PFlowPixelSeeds"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter2PFlowCkfTrackCandidatesForBTag = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltIter2MaskedMeasurementTrackerEventForBTag"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter2GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltIter2PFlowPixelSeedsForBTag"),
    useHitsSplitting = cms.bool(False)
)


process.hltIter2PFlowCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter2'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltIter2MaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter2PFlowCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter2PFlowCtfWithMaterialTracksForBTag = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter2'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltIter2MaskedMeasurementTrackerEventForBTag"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter2PFlowCkfTrackCandidatesForBTag"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltIter2PFlowPixelClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(50000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClusters"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltIter2PFlowPixelClusterCheckForBTag = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClustersRegForBTag"),
    MaxNumberOfCosmicClusters = cms.uint32(800000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClustersRegForBTag"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltIter2PFlowPixelHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltIter2PFlowPixelClusterCheck"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIter2PixelLayerTriplets"),
    trackingRegions = cms.InputTag("hltIter2PFlowPixelTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIter2PFlowPixelHitDoubletsForBTag = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltIter2PFlowPixelClusterCheckForBTag"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIter2PixelLayerTripletsForBTag"),
    trackingRegions = cms.InputTag("hltIter2PFlowPixelTrackingRegionsForBTag"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIter2PFlowPixelHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.3),
    CAPhiCut = cms.double(0.1),
    CAThetaCut = cms.double(0.004),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("hltIter2PFlowPixelHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.4),
        pt2 = cms.double(8.0),
        value1 = cms.double(100.0),
        value2 = cms.double(6.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIter2PFlowPixelHitTripletsForBTag = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.3),
    CAPhiCut = cms.double(0.1),
    CAThetaCut = cms.double(0.004),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("hltIter2PFlowPixelHitDoubletsForBTag"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(8.0),
        value1 = cms.double(100.0),
        value2 = cms.double(6.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIter2PFlowPixelSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("hltIter2PFlowPixelHitTriplets")
)


process.hltIter2PFlowPixelSeedsForBTag = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("hltIter2PFlowPixelHitTripletsForBTag")
)


process.hltIter2PFlowPixelTrackingRegions = cms.EDProducer("CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(0.8),
        deltaPhi = cms.double(0.8),
        input = cms.InputTag("hltIter1TrackAndTauJets4Iter2"),
        maxNRegions = cms.int32(100),
        maxNVertices = cms.int32(10),
        measurementTrackerName = cms.InputTag("hltIter2MaskedMeasurementTrackerEvent"),
        mode = cms.string('VerticesFixed'),
        nSigmaZBeamSpot = cms.double(3.0),
        nSigmaZVertex = cms.double(4.0),
        originRadius = cms.double(0.025),
        precise = cms.bool(True),
        ptMin = cms.double(0.4),
        searchOpt = cms.bool(True),
        vertexCollection = cms.InputTag("hltTrimmedPixelVertices"),
        whereToUseMeasurementTracker = cms.string('ForSiStrips'),
        zErrorBeamSpot = cms.double(15.0),
        zErrorVetex = cms.double(0.05)
    )
)


process.hltIter2PFlowPixelTrackingRegionsForBTag = cms.EDProducer("CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(0.8),
        deltaPhi = cms.double(0.8),
        input = cms.InputTag("hltSelector8CentralJetsL1FastJet"),
        maxNRegions = cms.int32(100),
        maxNVertices = cms.int32(10),
        measurementTrackerName = cms.InputTag("hltIter2MaskedMeasurementTrackerEventForBTag"),
        mode = cms.string('VerticesFixed'),
        nSigmaZBeamSpot = cms.double(3.0),
        nSigmaZVertex = cms.double(4.0),
        originRadius = cms.double(0.025),
        precise = cms.bool(True),
        ptMin = cms.double(0.8),
        searchOpt = cms.bool(True),
        vertexCollection = cms.InputTag("hltFastPVPixelVertices"),
        whereToUseMeasurementTracker = cms.string('ForSiStrips'),
        zErrorBeamSpot = cms.double(15.0),
        zErrorVetex = cms.double(0.05)
    )
)


process.hltIter2PFlowTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dr_par2 = cms.vdouble(3.40282346639e+38, 0.3, 0.3)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dz_par2 = cms.vdouble(3.40282346639e+38, 0.35, 0.35)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter2PFlowCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltTrimmedPixelVertices")
)


process.hltIter2PFlowTrackCutClassifierForBTag = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dr_par2 = cms.vdouble(3.40282346639e+38, 0.3, 0.3)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dz_par2 = cms.vdouble(3.40282346639e+38, 0.35, 0.35)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter2PFlowCtfWithMaterialTracksForBTag"),
    vertices = cms.InputTag("hltFastPVPixelVertices")
)


process.hltIter2PFlowTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter2PFlowTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter2PFlowTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter2PFlowCtfWithMaterialTracks")
)


process.hltIter2PFlowTrackSelectionHighPurityForBTag = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter2PFlowTrackCutClassifierForBTag","MVAValues"),
    originalQualVals = cms.InputTag("hltIter2PFlowTrackCutClassifierForBTag","QualityMasks"),
    originalSource = cms.InputTag("hltIter2PFlowCtfWithMaterialTracksForBTag")
)


process.hltIter2PixelLayerTriplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltIter2ClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltIter2ClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3', 
        'BPix2+BPix3+BPix4', 
        'BPix1+BPix3+BPix4', 
        'BPix1+BPix2+BPix4', 
        'BPix2+BPix3+FPix1_pos', 
        'BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos', 
        'BPix1+BPix2+FPix1_neg', 
        'BPix2+FPix1_pos+FPix2_pos', 
        'BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos', 
        'BPix1+FPix1_neg+FPix2_neg', 
        'FPix1_pos+FPix2_pos+FPix3_pos', 
        'FPix1_neg+FPix2_neg+FPix3_neg', 
        'BPix1+BPix3+FPix1_pos', 
        'BPix1+BPix2+FPix2_pos', 
        'BPix1+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix2_neg', 
        'BPix1+FPix2_neg+FPix3_neg', 
        'BPix1+FPix1_neg+FPix3_neg', 
        'BPix1+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_pos+FPix3_pos'
    )
)


process.hltIter2PixelLayerTripletsForBTag = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHitsRegForBTag'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltIter2ClustersRefRemovalForBTag"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHitsRegForBTag'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltIter2ClustersRefRemovalForBTag"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3', 
        'BPix2+BPix3+BPix4', 
        'BPix1+BPix3+BPix4', 
        'BPix1+BPix2+BPix4', 
        'BPix2+BPix3+FPix1_pos', 
        'BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos', 
        'BPix1+BPix2+FPix1_neg', 
        'BPix2+FPix1_pos+FPix2_pos', 
        'BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos', 
        'BPix1+FPix1_neg+FPix2_neg', 
        'FPix1_pos+FPix2_pos+FPix3_pos', 
        'FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltIterL3FromL1MuonPixelLayerQuadruplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4', 
        'BPix1+BPix2+BPix3+FPix1_pos', 
        'BPix1+BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos+FPix2_pos', 
        'BPix1+BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltIterL3FromL1MuonPixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltIterL3MuonPixelTracksFilter"),
    Fitter = cms.InputTag("hltIterL3MuonPixelTracksFitter"),
    SeedingHitSets = cms.InputTag("hltIterL3FromL1MuonPixelTracksHitQuadruplets"),
    passLabel = cms.string('')
)


process.hltIterL3FromL1MuonPixelTracksHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag(""),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIterL3FromL1MuonPixelLayerQuadruplets"),
    trackingRegions = cms.InputTag("hltIterL3FromL1MuonPixelTracksTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIterL3FromL1MuonPixelTracksHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.2),
    CAThetaCut = cms.double(0.005),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltIterL3FromL1MuonPixelTracksHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2.0),
        value1 = cms.double(200.0),
        value2 = cms.double(50.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIterL3FromL1MuonPixelTracksTrackingRegions = cms.EDProducer("CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(0.35),
        deltaPhi = cms.double(0.2),
        input = cms.InputTag("hltIterL3MuonL1MuonNoL2Selector"),
        maxNRegions = cms.int32(2),
        maxNVertices = cms.int32(1),
        measurementTrackerName = cms.InputTag(""),
        mode = cms.string('BeamSpotSigma'),
        nSigmaZBeamSpot = cms.double(4.0),
        nSigmaZVertex = cms.double(3.0),
        originRadius = cms.double(0.2),
        precise = cms.bool(True),
        ptMin = cms.double(10.0),
        searchOpt = cms.bool(False),
        vertexCollection = cms.InputTag("notUsed"),
        whereToUseMeasurementTracker = cms.string('Never'),
        zErrorBeamSpot = cms.double(24.2),
        zErrorVetex = cms.double(0.2)
    )
)


process.hltIterL3FromL1MuonPixelVertices = cms.EDProducer("PixelVertexProducer",
    Finder = cms.string('DivisiveVertexFinder'),
    Method2 = cms.bool(True),
    NTrkMin = cms.int32(2),
    PVcomparer = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPvClusterComparerForIT')
    ),
    PtMin = cms.double(1.0),
    TrackCollection = cms.InputTag("hltIterL3MuonPixelTracks"),
    UseError = cms.bool(True),
    Verbosity = cms.int32(0),
    WtAverage = cms.bool(True),
    ZOffset = cms.double(5.0),
    ZSeparation = cms.double(0.05),
    beamSpot = cms.InputTag("hltOnlineBeamSpot")
)


process.hltIterL3FromL1MuonTrimmedPixelVertices = cms.EDProducer("PixelVertexCollectionTrimmer",
    PVcomparer = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPvClusterComparerForIT')
    ),
    fractionSumPt2 = cms.double(0.3),
    maxVtx = cms.uint32(100),
    minSumPt2 = cms.double(0.0),
    src = cms.InputTag("hltIterL3FromL1MuonPixelVertices")
)


process.hltIterL3MuonAndMuonFromL1Merged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIterL3MuonMerged", "hltIter2IterL3FromL1MuonMerged"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIterL3MuonMerged", "hltIter2IterL3FromL1MuonMerged"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIterL3MuonCandidates = cms.EDProducer("L3MuonCandidateProducerFromMuons",
    InputObjects = cms.InputTag("hltIterL3Muons")
)


process.hltIterL3MuonL1MuonNoL2Selector = cms.EDProducer("HLTL1MuonNoL2Selector",
    CentralBxOnly = cms.bool(True),
    InputObjects = cms.InputTag("hltGtStage2Digis","Muon"),
    L1MaxEta = cms.double(5.0),
    L1MinPt = cms.double(-1.0),
    L1MinQuality = cms.uint32(7),
    L2CandTag = cms.InputTag("hltL2MuonCandidates"),
    SeedMapTag = cms.InputTag("hltL2Muons")
)


process.hltIterL3MuonMerged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIterL3OIMuonTrackSelectionHighPurity", "hltIter2IterL3MuonMerged"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIterL3OIMuonTrackSelectionHighPurity", "hltIter2IterL3MuonMerged"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltIterL3MuonPixelLayerQuadruplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4', 
        'BPix1+BPix2+BPix3+FPix1_pos', 
        'BPix1+BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos+FPix2_pos', 
        'BPix1+BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltIterL3MuonPixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltIterL3MuonPixelTracksFilter"),
    Fitter = cms.InputTag("hltIterL3MuonPixelTracksFitter"),
    SeedingHitSets = cms.InputTag("hltIterL3MuonPixelTracksHitQuadruplets"),
    passLabel = cms.string('')
)


process.hltIterL3MuonPixelTracksFilter = cms.EDProducer("PixelTrackFilterByKinematicsProducer",
    chi2 = cms.double(1000.0),
    nSigmaInvPtTolerance = cms.double(0.0),
    nSigmaTipMaxTolerance = cms.double(0.0),
    ptMin = cms.double(0.1),
    tipMax = cms.double(1.0)
)


process.hltIterL3MuonPixelTracksFitter = cms.EDProducer("PixelFitterByHelixProjectionsProducer",
    scaleErrorsForBPix1 = cms.bool(False),
    scaleFactor = cms.double(0.65)
)


process.hltIterL3MuonPixelTracksHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag(""),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltIterL3MuonPixelLayerQuadruplets"),
    trackingRegions = cms.InputTag("hltIterL3MuonPixelTracksTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltIterL3MuonPixelTracksHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.2),
    CAThetaCut = cms.double(0.005),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltIterL3MuonPixelTracksHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2.0),
        value1 = cms.double(200.0),
        value2 = cms.double(50.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltIterL3MuonPixelTracksTrackingRegions = cms.EDProducer("MuonTrackingRegionEDProducer",
    DeltaEta = cms.double(0.2),
    DeltaPhi = cms.double(0.15),
    DeltaR = cms.double(0.025),
    DeltaZ = cms.double(24.2),
    EtaR_UpperLimit_Par1 = cms.double(0.25),
    EtaR_UpperLimit_Par2 = cms.double(0.15),
    Eta_fixed = cms.bool(True),
    Eta_min = cms.double(0.0),
    MeasurementTrackerName = cms.InputTag(""),
    OnDemand = cms.int32(-1),
    PhiR_UpperLimit_Par1 = cms.double(0.6),
    PhiR_UpperLimit_Par2 = cms.double(0.2),
    Phi_fixed = cms.bool(True),
    Phi_min = cms.double(0.0),
    Pt_fixed = cms.bool(True),
    Pt_min = cms.double(2.0),
    Rescale_Dz = cms.double(4.0),
    Rescale_eta = cms.double(3.0),
    Rescale_phi = cms.double(3.0),
    UseVertex = cms.bool(False),
    Z_fixed = cms.bool(True),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    input = cms.InputTag("hltL2SelectorForL3IO"),
    maxRegions = cms.int32(5),
    precise = cms.bool(True),
    vertexCollection = cms.InputTag("notUsed")
)


process.hltIterL3MuonPixelVertices = cms.EDProducer("PixelVertexProducer",
    Finder = cms.string('DivisiveVertexFinder'),
    Method2 = cms.bool(True),
    NTrkMin = cms.int32(2),
    PVcomparer = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPvClusterComparerForIT')
    ),
    PtMin = cms.double(1.0),
    TrackCollection = cms.InputTag("hltIterL3MuonPixelTracks"),
    UseError = cms.bool(True),
    Verbosity = cms.int32(0),
    WtAverage = cms.bool(True),
    ZOffset = cms.double(5.0),
    ZSeparation = cms.double(0.05),
    beamSpot = cms.InputTag("hltOnlineBeamSpot")
)


process.hltIterL3MuonTrimmedPixelVertices = cms.EDProducer("PixelVertexCollectionTrimmer",
    PVcomparer = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPvClusterComparerForIT')
    ),
    fractionSumPt2 = cms.double(0.3),
    maxVtx = cms.uint32(100),
    minSumPt2 = cms.double(0.0),
    src = cms.InputTag("hltIterL3MuonPixelVertices")
)


process.hltIterL3Muons = cms.EDProducer("MuonIdProducer",
    CaloExtractorPSet = cms.PSet(
        CenterConeOnCalIntersection = cms.bool(False),
        ComponentName = cms.string('CaloExtractorByAssociator'),
        DR_Max = cms.double(1.0),
        DR_Veto_E = cms.double(0.07),
        DR_Veto_H = cms.double(0.1),
        DR_Veto_HO = cms.double(0.1),
        DepositInstanceLabels = cms.vstring(
            'ecal', 
            'hcal', 
            'ho'
        ),
        DepositLabel = cms.untracked.string('Cal'),
        NoiseTow_EB = cms.double(0.04),
        NoiseTow_EE = cms.double(0.15),
        Noise_EB = cms.double(0.025),
        Noise_EE = cms.double(0.1),
        Noise_HB = cms.double(0.2),
        Noise_HE = cms.double(0.2),
        Noise_HO = cms.double(0.2),
        PrintTimeReport = cms.untracked.bool(False),
        PropagatorName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(False),
            UseMuonNavigation = cms.untracked.bool(False)
        ),
        Threshold_E = cms.double(0.2),
        Threshold_H = cms.double(0.5),
        Threshold_HO = cms.double(0.5),
        TrackAssociatorParameters = cms.PSet(
            CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
            CaloTowerCollectionLabel = cms.InputTag("Notused"),
            DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
            EBRecHitCollectionLabel = cms.InputTag("Notused"),
            EERecHitCollectionLabel = cms.InputTag("Notused"),
            HBHERecHitCollectionLabel = cms.InputTag("Notused"),
            HORecHitCollectionLabel = cms.InputTag("Notused"),
            accountForTrajectoryChangeCalo = cms.bool(False),
            dREcal = cms.double(1.0),
            dREcalPreselection = cms.double(1.0),
            dRHcal = cms.double(1.0),
            dRHcalPreselection = cms.double(1.0),
            dRMuon = cms.double(9999.0),
            dRMuonPreselection = cms.double(0.2),
            dRPreshowerPreselection = cms.double(0.2),
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            propagateAllDirections = cms.bool(True),
            trajectoryUncertaintyTolerance = cms.double(-1.0),
            truthMatch = cms.bool(False),
            useCalo = cms.bool(True),
            useEcal = cms.bool(False),
            useHO = cms.bool(False),
            useHcal = cms.bool(False),
            useMuon = cms.bool(False),
            usePreshower = cms.bool(False)
        ),
        UseRecHitsFlag = cms.bool(False)
    ),
    JetExtractorPSet = cms.PSet(
        ComponentName = cms.string('JetExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.1),
        ExcludeMuonVeto = cms.bool(True),
        JetCollectionLabel = cms.InputTag("Notused"),
        PrintTimeReport = cms.untracked.bool(False),
        PropagatorName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(False),
            UseMuonNavigation = cms.untracked.bool(False)
        ),
        Threshold = cms.double(5.0),
        TrackAssociatorParameters = cms.PSet(
            CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
            CaloTowerCollectionLabel = cms.InputTag("Notused"),
            DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
            EBRecHitCollectionLabel = cms.InputTag("Notused"),
            EERecHitCollectionLabel = cms.InputTag("Notused"),
            HBHERecHitCollectionLabel = cms.InputTag("Notused"),
            HORecHitCollectionLabel = cms.InputTag("Notused"),
            accountForTrajectoryChangeCalo = cms.bool(False),
            dREcal = cms.double(0.5),
            dREcalPreselection = cms.double(0.5),
            dRHcal = cms.double(0.5),
            dRHcalPreselection = cms.double(0.5),
            dRMuon = cms.double(9999.0),
            dRMuonPreselection = cms.double(0.2),
            dRPreshowerPreselection = cms.double(0.2),
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            propagateAllDirections = cms.bool(True),
            trajectoryUncertaintyTolerance = cms.double(-1.0),
            truthMatch = cms.bool(False),
            useCalo = cms.bool(True),
            useEcal = cms.bool(False),
            useHO = cms.bool(False),
            useHcal = cms.bool(False),
            useMuon = cms.bool(False),
            usePreshower = cms.bool(False)
        )
    ),
    MuonCaloCompatibility = cms.PSet(
        MuonTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root'),
        PionTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root'),
        allSiPMHO = cms.bool(False),
        delta_eta = cms.double(0.02),
        delta_phi = cms.double(0.02)
    ),
    TimingFillerParameters = cms.PSet(
        CSCTimingParameters = cms.PSet(
            CSCStripError = cms.double(7.0),
            CSCStripTimeOffset = cms.double(0.0),
            CSCTimeOffset = cms.double(0.0),
            CSCWireError = cms.double(8.6),
            CSCWireTimeOffset = cms.double(0.0),
            CSCsegments = cms.InputTag("hltCscSegments"),
            MatchParameters = cms.PSet(
                CSCsegments = cms.InputTag("hltCscSegments"),
                DTradius = cms.double(0.01),
                DTsegments = cms.InputTag("hltDt4DSegments"),
                TightMatchCSC = cms.bool(True),
                TightMatchDT = cms.bool(False)
            ),
            PruneCut = cms.double(100.0),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                RPCLayers = cms.bool(True)
            ),
            UseStripTime = cms.bool(True),
            UseWireTime = cms.bool(True),
            debug = cms.bool(False)
        ),
        DTTimingParameters = cms.PSet(
            DTTimeOffset = cms.double(2.7),
            DTsegments = cms.InputTag("hltDt4DSegments"),
            DoWireCorr = cms.bool(False),
            DropTheta = cms.bool(True),
            HitError = cms.double(6.0),
            HitsMin = cms.int32(5),
            MatchParameters = cms.PSet(
                CSCsegments = cms.InputTag("hltCscSegments"),
                DTradius = cms.double(0.01),
                DTsegments = cms.InputTag("hltDt4DSegments"),
                TightMatchCSC = cms.bool(True),
                TightMatchDT = cms.bool(False)
            ),
            PruneCut = cms.double(10000.0),
            RequireBothProjections = cms.bool(False),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                RPCLayers = cms.bool(True)
            ),
            UseSegmentT0 = cms.bool(False),
            debug = cms.bool(False)
        ),
        EcalEnergyCut = cms.double(0.4),
        ErrorCSC = cms.double(7.4),
        ErrorDT = cms.double(6.0),
        ErrorEB = cms.double(2.085),
        ErrorEE = cms.double(6.95),
        UseCSC = cms.bool(True),
        UseDT = cms.bool(True),
        UseECAL = cms.bool(True)
    ),
    TrackAssociatorParameters = cms.PSet(
        CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
        CaloTowerCollectionLabel = cms.InputTag("Notused"),
        DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
        EBRecHitCollectionLabel = cms.InputTag("Notused"),
        EERecHitCollectionLabel = cms.InputTag("Notused"),
        HBHERecHitCollectionLabel = cms.InputTag("Notused"),
        HORecHitCollectionLabel = cms.InputTag("Notused"),
        accountForTrajectoryChangeCalo = cms.bool(False),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        dRHcal = cms.double(9999.0),
        dRHcalPreselection = cms.double(0.2),
        dRMuon = cms.double(9999.0),
        dRMuonPreselection = cms.double(0.2),
        dRPreshowerPreselection = cms.double(0.2),
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        propagateAllDirections = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        truthMatch = cms.bool(False),
        useCalo = cms.bool(False),
        useEcal = cms.bool(False),
        useHO = cms.bool(False),
        useHcal = cms.bool(False),
        useMuon = cms.bool(True),
        usePreshower = cms.bool(False)
    ),
    TrackExtractorPSet = cms.PSet(
        BeamSpotLabel = cms.InputTag("hltOnlineBeamSpot"),
        BeamlineOption = cms.string('BeamSpotFromEvent'),
        Chi2Ndof_Max = cms.double(1e+64),
        Chi2Prob_Min = cms.double(-1.0),
        ComponentName = cms.string('TrackExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.01),
        Diff_r = cms.double(0.1),
        Diff_z = cms.double(0.2),
        NHits_Min = cms.uint32(0),
        Pt_Min = cms.double(-1.0),
        inputTrackCollection = cms.InputTag("hltIter2IterL3FromL1MuonMerged")
    ),
    TrackerKinkFinderParameters = cms.PSet(
        diagonalOnly = cms.bool(False),
        usePosition = cms.bool(False)
    ),
    addExtraSoftMuons = cms.bool(False),
    arbitrateTrackerMuons = cms.bool(True),
    arbitrationCleanerOptions = cms.PSet(
        ClusterDPhi = cms.double(0.6),
        ClusterDTheta = cms.double(0.02),
        Clustering = cms.bool(True),
        ME1a = cms.bool(True),
        Overlap = cms.bool(True),
        OverlapDPhi = cms.double(0.0786),
        OverlapDTheta = cms.double(0.02)
    ),
    debugWithTruthMatching = cms.bool(False),
    ecalDepositName = cms.string('ecal'),
    fillCaloCompatibility = cms.bool(False),
    fillEnergy = cms.bool(False),
    fillGlobalTrackQuality = cms.bool(False),
    fillGlobalTrackRefits = cms.bool(False),
    fillIsolation = cms.bool(False),
    fillMatching = cms.bool(True),
    fillTrackerKink = cms.bool(False),
    globalTrackQualityInputTag = cms.InputTag("glbTrackQual"),
    hcalDepositName = cms.string('hcal'),
    hoDepositName = cms.string('ho'),
    inputCollectionLabels = cms.VInputTag("hltIter2IterL3FromL1MuonMerged", "hltL3MuonsIterL3Links"),
    inputCollectionTypes = cms.vstring(
        'inner tracks', 
        'links'
    ),
    jetDepositName = cms.string('jets'),
    maxAbsDx = cms.double(3.0),
    maxAbsDy = cms.double(9999.0),
    maxAbsEta = cms.double(3.0),
    maxAbsPullX = cms.double(4.0),
    maxAbsPullY = cms.double(9999.0),
    minCaloCompatibility = cms.double(0.6),
    minNumberOfMatches = cms.int32(1),
    minP = cms.double(0.0),
    minPCaloMuon = cms.double(1000000000.0),
    minPt = cms.double(8.0),
    ptThresholdToFillCandidateP4WithGlobalFit = cms.double(200.0),
    runArbitrationCleaner = cms.bool(False),
    sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double(2.0),
    trackDepositName = cms.string('tracker'),
    writeIsoDeposits = cms.bool(False)
)


process.hltIterL3MuonsFromL2 = cms.EDProducer("L3TrackCombiner",
    labels = cms.VInputTag("hltL3MuonsIterL3OI", "hltL3MuonsIterL3IO")
)


process.hltIterL3MuonsFromL2LinksCombination = cms.EDProducer("L3TrackLinksCombiner",
    labels = cms.VInputTag("hltL3MuonsIterL3OI", "hltL3MuonsIterL3IO")
)


process.hltIterL3OIL3MuonCandidates = cms.EDProducer("L3MuonCandidateProducer",
    InputLinksObjects = cms.InputTag("hltIterL3OIL3MuonsLinksCombination"),
    InputObjects = cms.InputTag("hltIterL3OIL3Muons"),
    MuonPtOption = cms.string('Tracker')
)


process.hltIterL3OIL3Muons = cms.EDProducer("L3TrackCombiner",
    labels = cms.VInputTag("hltL3MuonsIterL3OI")
)


process.hltIterL3OIL3MuonsLinksCombination = cms.EDProducer("L3TrackLinksCombiner",
    labels = cms.VInputTag("hltL3MuonsIterL3OI")
)


process.hltIterL3OIMuCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('iter10'),
    Fitter = cms.string('hltESPKFFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string('hltESPMeasurementTracker'),
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterial'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIterL3OITrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.hltIterL3OIMuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(True),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
            d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
            dr_exp = cms.vint32(4, 4, 2147483647),
            dr_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dr_par2 = cms.vdouble(0.3, 0.3, 3.40282346639e+38)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 2147483647),
            dz_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dz_par2 = cms.vdouble(0.35, 0.35, 3.40282346639e+38)
        ),
        maxChi2 = cms.vdouble(3.40282346639e+38, 3.40282346639e+38, 3.40282346639e+38),
        maxChi2n = cms.vdouble(10.0, 1.0, 0.4),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
        maxLostLayers = cms.vint32(4, 3, 2),
        min3DLayers = cms.vint32(1, 2, 1),
        minLayers = cms.vint32(3, 5, 5),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 1)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIterL3OIMuCtfWithMaterialTracks"),
    vertices = cms.InputTag("Notused")
)


process.hltIterL3OIMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIterL3OIMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIterL3OIMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIterL3OIMuCtfWithMaterialTracks")
)


process.hltIterL3OISeedsFromL2Muons = cms.EDProducer("TSGForOI",
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    SF1 = cms.double(3.0),
    SF2 = cms.double(4.0),
    SF3 = cms.double(5.0),
    SF4 = cms.double(7.0),
    SF5 = cms.double(10.0),
    UseHitLessSeeds = cms.bool(True),
    UseStereoLayersInTEC = cms.bool(False),
    adjustErrorsDynamicallyForHitless = cms.bool(True),
    adjustErrorsDynamicallyForHits = cms.bool(True),
    debug = cms.untracked.bool(False),
    estimator = cms.string('hltESPChi2MeasurementEstimator100'),
    eta1 = cms.double(1.0),
    eta2 = cms.double(1.4),
    fixedErrorRescaleFactorForHitless = cms.double(10.0),
    fixedErrorRescaleFactorForHits = cms.double(3.0),
    hitsToTry = cms.int32(3),
    layersToTry = cms.int32(2),
    maxEtaForTOB = cms.double(1.8),
    maxSeeds = cms.uint32(5),
    minEtaForTEC = cms.double(0.7),
    pT1 = cms.double(13.0),
    pT2 = cms.double(30.0),
    pT3 = cms.double(70.0),
    propagatorName = cms.string('PropagatorWithMaterial'),
    src = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    tsosDiff = cms.double(0.03)
)


process.hltIterL3OITrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltSiStripClusters"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string(''),
    TrajectoryBuilder = cms.string('CkfTrajectoryBuilder'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTPSetMuonCkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('muonSeededTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    src = cms.InputTag("hltIterL3OISeedsFromL2Muons"),
    useHitsSplitting = cms.bool(False)
)


process.hltL2MuonCandidates = cms.EDProducer("L2MuonCandidateProducer",
    InputObjects = cms.InputTag("hltL2Muons","UpdatedAtVtx")
)


process.hltL2MuonSeeds = cms.EDProducer("L2MuonSeedGeneratorFromL1T",
    CentralBxOnly = cms.bool(True),
    EtaMatchingBins = cms.vdouble(0.0, 2.5),
    GMTReadoutCollection = cms.InputTag(""),
    InputObjects = cms.InputTag("hltGtStage2Digis","Muon"),
    L1MaxEta = cms.double(2.5),
    L1MinPt = cms.double(0.0),
    L1MinQuality = cms.uint32(7),
    MatchDR = cms.vdouble(0.3),
    OfflineSeedLabel = cms.untracked.InputTag("hltL2OfflineMuonSeeds"),
    Propagator = cms.string('SteppingHelixPropagatorAny'),
    ServiceParameters = cms.PSet(
        Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny'),
        RPCLayers = cms.bool(True),
        UseMuonNavigation = cms.untracked.bool(True)
    ),
    UseOfflineSeed = cms.untracked.bool(True),
    UseUnassociatedL1 = cms.bool(False)
)


process.hltL2Muons = cms.EDProducer("L2MuonProducer",
    DoSeedRefit = cms.bool(False),
    InputObjects = cms.InputTag("hltL2MuonSeeds"),
    L2TrajBuilderParameters = cms.PSet(
        BWFilterParameters = cms.PSet(
            BWSeedType = cms.string('fromGenerator'),
            CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
            DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
            EnableCSCMeasurement = cms.bool(True),
            EnableDTMeasurement = cms.bool(True),
            EnableRPCMeasurement = cms.bool(True),
            FitDirection = cms.string('outsideIn'),
            MaxChi2 = cms.double(100.0),
            MuonTrajectoryUpdatorParameters = cms.PSet(
                ExcludeRPCFromFit = cms.bool(False),
                Granularity = cms.int32(0),
                MaxChi2 = cms.double(25.0),
                RescaleError = cms.bool(False),
                RescaleErrorFactor = cms.double(100.0),
                UseInvalidHits = cms.bool(True)
            ),
            NumberOfSigma = cms.double(3.0),
            Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            RPCRecSegmentLabel = cms.InputTag("hltRpcRecHits")
        ),
        DoBackwardFilter = cms.bool(True),
        DoRefit = cms.bool(False),
        DoSeedRefit = cms.bool(False),
        FilterParameters = cms.PSet(
            CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
            DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
            EnableCSCMeasurement = cms.bool(True),
            EnableDTMeasurement = cms.bool(True),
            EnableRPCMeasurement = cms.bool(True),
            FitDirection = cms.string('insideOut'),
            MaxChi2 = cms.double(1000.0),
            MuonTrajectoryUpdatorParameters = cms.PSet(
                ExcludeRPCFromFit = cms.bool(False),
                Granularity = cms.int32(0),
                MaxChi2 = cms.double(25.0),
                RescaleError = cms.bool(False),
                RescaleErrorFactor = cms.double(100.0),
                UseInvalidHits = cms.bool(True)
            ),
            NumberOfSigma = cms.double(3.0),
            Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            RPCRecSegmentLabel = cms.InputTag("hltRpcRecHits")
        ),
        NavigationType = cms.string('Standard'),
        SeedPosition = cms.string('in'),
        SeedPropagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        SeedTransformerParameters = cms.PSet(
            Fitter = cms.string('hltESPKFFittingSmootherForL2Muon'),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            NMinRecHits = cms.uint32(2),
            Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            RescaleError = cms.double(100.0),
            UseSubRecHits = cms.bool(False)
        )
    ),
    MuonTrajectoryBuilder = cms.string('Exhaustive'),
    SeedTransformerParameters = cms.PSet(
        Fitter = cms.string('hltESPKFFittingSmootherForL2Muon'),
        MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
        NMinRecHits = cms.uint32(2),
        Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        RescaleError = cms.double(100.0),
        UseSubRecHits = cms.bool(False)
    ),
    ServiceParameters = cms.PSet(
        Propagators = cms.untracked.vstring(
            'hltESPFastSteppingHelixPropagatorAny', 
            'hltESPFastSteppingHelixPropagatorOpposite'
        ),
        RPCLayers = cms.bool(True),
        UseMuonNavigation = cms.untracked.bool(True)
    ),
    TrackLoaderParameters = cms.PSet(
        DoSmoothing = cms.bool(False),
        MuonUpdatorAtVertexParameters = cms.PSet(
            BeamSpotPosition = cms.vdouble(0.0, 0.0, 0.0),
            BeamSpotPositionErrors = cms.vdouble(0.1, 0.1, 5.3),
            MaxChi2 = cms.double(1000000.0),
            Propagator = cms.string('hltESPFastSteppingHelixPropagatorOpposite')
        ),
        Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
        TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        VertexConstraint = cms.bool(True),
        beamSpot = cms.InputTag("hltOnlineBeamSpot")
    )
)


process.hltL2OfflineMuonSeeds = cms.EDProducer("MuonSeedGenerator",
    CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
    CSC_01 = cms.vdouble(
        0.166, 0.0, 0.0, 0.031, 0.0, 
        0.0
    ),
    CSC_01_1_scale = cms.vdouble(-1.915329, 0.0),
    CSC_02 = cms.vdouble(
        0.612, -0.207, 0.0, 0.067, -0.001, 
        0.0
    ),
    CSC_03 = cms.vdouble(
        0.787, -0.338, 0.029, 0.101, -0.008, 
        0.0
    ),
    CSC_12 = cms.vdouble(
        -0.161, 0.254, -0.047, 0.042, -0.007, 
        0.0
    ),
    CSC_12_1_scale = cms.vdouble(-6.434242, 0.0),
    CSC_12_2_scale = cms.vdouble(-1.63622, 0.0),
    CSC_12_3_scale = cms.vdouble(-1.63622, 0.0),
    CSC_13 = cms.vdouble(
        0.901, -1.302, 0.533, 0.045, 0.005, 
        0.0
    ),
    CSC_13_2_scale = cms.vdouble(-6.077936, 0.0),
    CSC_13_3_scale = cms.vdouble(-1.701268, 0.0),
    CSC_14 = cms.vdouble(
        0.606, -0.181, -0.002, 0.111, -0.003, 
        0.0
    ),
    CSC_14_3_scale = cms.vdouble(-1.969563, 0.0),
    CSC_23 = cms.vdouble(
        -0.081, 0.113, -0.029, 0.015, 0.008, 
        0.0
    ),
    CSC_23_1_scale = cms.vdouble(-19.084285, 0.0),
    CSC_23_2_scale = cms.vdouble(-6.079917, 0.0),
    CSC_24 = cms.vdouble(
        0.004, 0.021, -0.002, 0.053, 0.0, 
        0.0
    ),
    CSC_24_1_scale = cms.vdouble(-6.055701, 0.0),
    CSC_34 = cms.vdouble(
        0.062, -0.067, 0.019, 0.021, 0.003, 
        0.0
    ),
    CSC_34_1_scale = cms.vdouble(-11.520507, 0.0),
    DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
    DT_12 = cms.vdouble(
        0.183, 0.054, -0.087, 0.028, 0.002, 
        0.0
    ),
    DT_12_1_scale = cms.vdouble(-3.692398, 0.0),
    DT_12_2_scale = cms.vdouble(-3.518165, 0.0),
    DT_13 = cms.vdouble(
        0.315, 0.068, -0.127, 0.051, -0.002, 
        0.0
    ),
    DT_13_1_scale = cms.vdouble(-4.520923, 0.0),
    DT_13_2_scale = cms.vdouble(-4.257687, 0.0),
    DT_14 = cms.vdouble(
        0.359, 0.052, -0.107, 0.072, -0.004, 
        0.0
    ),
    DT_14_1_scale = cms.vdouble(-5.644816, 0.0),
    DT_14_2_scale = cms.vdouble(-4.808546, 0.0),
    DT_23 = cms.vdouble(
        0.13, 0.023, -0.057, 0.028, 0.004, 
        0.0
    ),
    DT_23_1_scale = cms.vdouble(-5.320346, 0.0),
    DT_23_2_scale = cms.vdouble(-5.117625, 0.0),
    DT_24 = cms.vdouble(
        0.176, 0.014, -0.051, 0.051, 0.003, 
        0.0
    ),
    DT_24_1_scale = cms.vdouble(-7.490909, 0.0),
    DT_24_2_scale = cms.vdouble(-6.63094, 0.0),
    DT_34 = cms.vdouble(
        0.044, 0.004, -0.013, 0.029, 0.003, 
        0.0
    ),
    DT_34_1_scale = cms.vdouble(-13.783765, 0.0),
    DT_34_2_scale = cms.vdouble(-11.901897, 0.0),
    EnableCSCMeasurement = cms.bool(True),
    EnableDTMeasurement = cms.bool(True),
    EnableME0Measurement = cms.bool(False),
    ME0RecSegmentLabel = cms.InputTag("me0Segments"),
    OL_1213 = cms.vdouble(
        0.96, -0.737, 0.0, 0.052, 0.0, 
        0.0
    ),
    OL_1213_0_scale = cms.vdouble(-4.488158, 0.0),
    OL_1222 = cms.vdouble(
        0.848, -0.591, 0.0, 0.062, 0.0, 
        0.0
    ),
    OL_1222_0_scale = cms.vdouble(-5.810449, 0.0),
    OL_1232 = cms.vdouble(
        0.184, 0.0, 0.0, 0.066, 0.0, 
        0.0
    ),
    OL_1232_0_scale = cms.vdouble(-5.964634, 0.0),
    OL_2213 = cms.vdouble(
        0.117, 0.0, 0.0, 0.044, 0.0, 
        0.0
    ),
    OL_2213_0_scale = cms.vdouble(-7.239789, 0.0),
    OL_2222 = cms.vdouble(
        0.107, 0.0, 0.0, 0.04, 0.0, 
        0.0
    ),
    OL_2222_0_scale = cms.vdouble(-7.667231, 0.0),
    SMB_10 = cms.vdouble(
        1.387, -0.038, 0.0, 0.19, 0.0, 
        0.0
    ),
    SMB_10_0_scale = cms.vdouble(2.448566, 0.0),
    SMB_11 = cms.vdouble(
        1.247, 0.72, -0.802, 0.229, -0.075, 
        0.0
    ),
    SMB_11_0_scale = cms.vdouble(2.56363, 0.0),
    SMB_12 = cms.vdouble(
        2.128, -0.956, 0.0, 0.199, 0.0, 
        0.0
    ),
    SMB_12_0_scale = cms.vdouble(2.283221, 0.0),
    SMB_20 = cms.vdouble(
        1.011, -0.052, 0.0, 0.188, 0.0, 
        0.0
    ),
    SMB_20_0_scale = cms.vdouble(1.486168, 0.0),
    SMB_21 = cms.vdouble(
        1.043, -0.124, 0.0, 0.183, 0.0, 
        0.0
    ),
    SMB_21_0_scale = cms.vdouble(1.58384, 0.0),
    SMB_22 = cms.vdouble(
        1.474, -0.758, 0.0, 0.185, 0.0, 
        0.0
    ),
    SMB_22_0_scale = cms.vdouble(1.346681, 0.0),
    SMB_30 = cms.vdouble(
        0.505, -0.022, 0.0, 0.215, 0.0, 
        0.0
    ),
    SMB_30_0_scale = cms.vdouble(-3.629838, 0.0),
    SMB_31 = cms.vdouble(
        0.549, -0.145, 0.0, 0.207, 0.0, 
        0.0
    ),
    SMB_31_0_scale = cms.vdouble(-3.323768, 0.0),
    SMB_32 = cms.vdouble(
        0.67, -0.327, 0.0, 0.22, 0.0, 
        0.0
    ),
    SMB_32_0_scale = cms.vdouble(-3.054156, 0.0),
    SME_11 = cms.vdouble(
        3.295, -1.527, 0.112, 0.378, 0.02, 
        0.0
    ),
    SME_11_0_scale = cms.vdouble(1.325085, 0.0),
    SME_12 = cms.vdouble(
        0.102, 0.599, 0.0, 0.38, 0.0, 
        0.0
    ),
    SME_12_0_scale = cms.vdouble(2.279181, 0.0),
    SME_13 = cms.vdouble(
        -1.286, 1.711, 0.0, 0.356, 0.0, 
        0.0
    ),
    SME_13_0_scale = cms.vdouble(0.104905, 0.0),
    SME_21 = cms.vdouble(
        -0.529, 1.194, -0.358, 0.472, 0.086, 
        0.0
    ),
    SME_21_0_scale = cms.vdouble(-0.040862, 0.0),
    SME_22 = cms.vdouble(
        -1.207, 1.491, -0.251, 0.189, 0.243, 
        0.0
    ),
    SME_22_0_scale = cms.vdouble(-3.457901, 0.0),
    SME_31 = cms.vdouble(
        -1.594, 1.482, -0.317, 0.487, 0.097, 
        0.0
    ),
    SME_32 = cms.vdouble(
        -0.901, 1.333, -0.47, 0.41, 0.073, 
        0.0
    ),
    SME_41 = cms.vdouble(
        -0.003, 0.005, 0.005, 0.608, 0.076, 
        0.0
    ),
    SME_42 = cms.vdouble(
        -0.003, 0.005, 0.005, 0.608, 0.076, 
        0.0
    ),
    beamSpotTag = cms.InputTag("hltOnlineBeamSpot"),
    crackEtas = cms.vdouble(0.2, 1.6, 1.7),
    crackWindow = cms.double(0.04),
    deltaEtaCrackSearchWindow = cms.double(0.25),
    deltaEtaSearchWindow = cms.double(0.2),
    deltaPhiSearchWindow = cms.double(0.25),
    scaleDT = cms.bool(True)
)


process.hltL2SelectorForL3IO = cms.EDProducer("HLTMuonL2SelectorForL3IO",
    InputLinks = cms.InputTag("hltIterL3OIL3MuonsLinksCombination"),
    MaxNormalizedChi2 = cms.double(20.0),
    MaxPtDifference = cms.double(0.3),
    MinNhits = cms.int32(1),
    MinNmuonHits = cms.int32(1),
    applyL3Filters = cms.bool(False),
    l2Src = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    l3OISrc = cms.InputTag("hltIterL3OIL3MuonCandidates")
)


process.hltL3MuonsIterL3IO = cms.EDProducer("L3MuonProducer",
    L3TrajBuilderParameters = cms.PSet(
        GlbRefitterParameters = cms.PSet(
            CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
            Chi2CutCSC = cms.double(150.0),
            Chi2CutDT = cms.double(10.0),
            Chi2CutRPC = cms.double(1.0),
            DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
            DYTthrs = cms.vint32(30, 15),
            DoPredictionsOnly = cms.bool(False),
            Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
            HitThreshold = cms.int32(1),
            MuonHitsOption = cms.int32(1),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            PropDirForCosmics = cms.bool(False),
            Propagator = cms.string('hltESPSmartPropagatorAny'),
            RefitDirection = cms.string('insideOut'),
            RefitFlag = cms.bool(True),
            RefitRPCHits = cms.bool(True),
            SkipStation = cms.int32(-1),
            TrackerRecHitBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
            TrackerSkipSection = cms.int32(-1),
            TrackerSkipSystem = cms.int32(-1)
        ),
        GlobalMuonTrackMatcher = cms.PSet(
            Chi2Cut_1 = cms.double(50.0),
            Chi2Cut_2 = cms.double(50.0),
            Chi2Cut_3 = cms.double(200.0),
            DeltaDCut_1 = cms.double(40.0),
            DeltaDCut_2 = cms.double(10.0),
            DeltaDCut_3 = cms.double(15.0),
            DeltaRCut_1 = cms.double(0.1),
            DeltaRCut_2 = cms.double(0.2),
            DeltaRCut_3 = cms.double(1.0),
            Eta_threshold = cms.double(1.2),
            LocChi2Cut = cms.double(0.001),
            MinP = cms.double(2.5),
            MinPt = cms.double(1.0),
            Propagator = cms.string('hltESPSmartPropagator'),
            Pt_threshold1 = cms.double(0.0),
            Pt_threshold2 = cms.double(999999999.0),
            Quality_1 = cms.double(20.0),
            Quality_2 = cms.double(15.0),
            Quality_3 = cms.double(7.0)
        ),
        MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
        MuonTrackingRegionBuilder = cms.PSet(
            DeltaEta = cms.double(0.04),
            DeltaPhi = cms.double(0.15),
            DeltaR = cms.double(0.025),
            DeltaZ = cms.double(24.2),
            EtaR_UpperLimit_Par1 = cms.double(0.25),
            EtaR_UpperLimit_Par2 = cms.double(0.15),
            Eta_fixed = cms.bool(True),
            Eta_min = cms.double(0.1),
            MeasurementTrackerName = cms.InputTag("hltESPMeasurementTracker"),
            OnDemand = cms.int32(-1),
            PhiR_UpperLimit_Par1 = cms.double(0.6),
            PhiR_UpperLimit_Par2 = cms.double(0.2),
            Phi_fixed = cms.bool(True),
            Phi_min = cms.double(0.1),
            Pt_fixed = cms.bool(True),
            Pt_min = cms.double(3.0),
            Rescale_Dz = cms.double(4.0),
            Rescale_eta = cms.double(3.0),
            Rescale_phi = cms.double(3.0),
            UseVertex = cms.bool(False),
            Z_fixed = cms.bool(True),
            beamSpot = cms.InputTag("hltOnlineBeamSpot"),
            input = cms.InputTag("hltL2SelectorForL3IO"),
            maxRegions = cms.int32(2),
            precise = cms.bool(True),
            vertexCollection = cms.InputTag("pixelVertices")
        ),
        PCut = cms.double(2.5),
        PtCut = cms.double(1.0),
        RefitRPCHits = cms.bool(True),
        ScaleTECxFactor = cms.double(-1.0),
        ScaleTECyFactor = cms.double(-1.0),
        TrackTransformer = cms.PSet(
            DoPredictionsOnly = cms.bool(False),
            Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            Propagator = cms.string('hltESPSmartPropagatorAny'),
            RefitDirection = cms.string('insideOut'),
            RefitRPCHits = cms.bool(True),
            Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
            TrackerRecHitBuilder = cms.string('hltESPTTRHBWithTrackAngle')
        ),
        TrackerPropagator = cms.string('SteppingHelixPropagatorAny'),
        TrackerRecHitBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        matchToSeeds = cms.bool(True),
        tkTrajBeamSpot = cms.InputTag("hltOnlineBeamSpot"),
        tkTrajLabel = cms.InputTag("hltIter2IterL3MuonMerged"),
        tkTrajMaxChi2 = cms.double(9999.0),
        tkTrajMaxDXYBeamSpot = cms.double(9999.0),
        tkTrajUseVertex = cms.bool(False),
        tkTrajVertex = cms.InputTag("hltIterL3MuonPixelVertices")
    ),
    MuonCollectionLabel = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    ServiceParameters = cms.PSet(
        Propagators = cms.untracked.vstring(
            'hltESPSmartPropagatorAny', 
            'SteppingHelixPropagatorAny', 
            'hltESPSmartPropagator', 
            'hltESPSteppingHelixPropagatorOpposite'
        ),
        RPCLayers = cms.bool(True),
        UseMuonNavigation = cms.untracked.bool(True)
    ),
    TrackLoaderParameters = cms.PSet(
        DoSmoothing = cms.bool(False),
        MuonSeededTracksInstance = cms.untracked.string('L2Seeded'),
        MuonUpdatorAtVertexParameters = cms.PSet(
            BeamSpotPositionErrors = cms.vdouble(0.1, 0.1, 5.3),
            MaxChi2 = cms.double(1000000.0),
            Propagator = cms.string('hltESPSteppingHelixPropagatorOpposite')
        ),
        PutTkTrackIntoEvent = cms.untracked.bool(False),
        SmoothTkTrack = cms.untracked.bool(False),
        Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
        VertexConstraint = cms.bool(False),
        beamSpot = cms.InputTag("hltOnlineBeamSpot")
    )
)


process.hltL3MuonsIterL3Links = cms.EDProducer("MuonLinksProducerForHLT",
    InclusiveTrackerTrackCollection = cms.InputTag("hltIter2IterL3FromL1MuonMerged"),
    LinkCollection = cms.InputTag("hltIterL3MuonsFromL2LinksCombination"),
    pMin = cms.double(2.5),
    ptMin = cms.double(2.5),
    shareHitFraction = cms.double(0.19)
)


process.hltL3MuonsIterL3OI = cms.EDProducer("L3MuonProducer",
    L3TrajBuilderParameters = cms.PSet(
        GlbRefitterParameters = cms.PSet(
            CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
            Chi2CutCSC = cms.double(150.0),
            Chi2CutDT = cms.double(10.0),
            Chi2CutRPC = cms.double(1.0),
            DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
            DYTthrs = cms.vint32(30, 15),
            DoPredictionsOnly = cms.bool(False),
            Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
            HitThreshold = cms.int32(1),
            MuonHitsOption = cms.int32(1),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            PropDirForCosmics = cms.bool(False),
            Propagator = cms.string('hltESPSmartPropagatorAny'),
            RefitDirection = cms.string('insideOut'),
            RefitFlag = cms.bool(True),
            RefitRPCHits = cms.bool(True),
            SkipStation = cms.int32(-1),
            TrackerRecHitBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
            TrackerSkipSection = cms.int32(-1),
            TrackerSkipSystem = cms.int32(-1)
        ),
        GlobalMuonTrackMatcher = cms.PSet(
            Chi2Cut_1 = cms.double(50.0),
            Chi2Cut_2 = cms.double(50.0),
            Chi2Cut_3 = cms.double(200.0),
            DeltaDCut_1 = cms.double(40.0),
            DeltaDCut_2 = cms.double(10.0),
            DeltaDCut_3 = cms.double(15.0),
            DeltaRCut_1 = cms.double(0.1),
            DeltaRCut_2 = cms.double(0.2),
            DeltaRCut_3 = cms.double(1.0),
            Eta_threshold = cms.double(1.2),
            LocChi2Cut = cms.double(0.001),
            MinP = cms.double(2.5),
            MinPt = cms.double(1.0),
            Propagator = cms.string('hltESPSmartPropagator'),
            Pt_threshold1 = cms.double(0.0),
            Pt_threshold2 = cms.double(999999999.0),
            Quality_1 = cms.double(20.0),
            Quality_2 = cms.double(15.0),
            Quality_3 = cms.double(7.0)
        ),
        MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
        MuonTrackingRegionBuilder = cms.PSet(
            DeltaEta = cms.double(0.2),
            DeltaPhi = cms.double(0.15),
            DeltaR = cms.double(0.025),
            DeltaZ = cms.double(24.2),
            EtaR_UpperLimit_Par1 = cms.double(0.25),
            EtaR_UpperLimit_Par2 = cms.double(0.15),
            Eta_fixed = cms.bool(True),
            Eta_min = cms.double(0.1),
            MeasurementTrackerName = cms.InputTag("hltESPMeasurementTracker"),
            OnDemand = cms.int32(-1),
            PhiR_UpperLimit_Par1 = cms.double(0.6),
            PhiR_UpperLimit_Par2 = cms.double(0.2),
            Phi_fixed = cms.bool(True),
            Phi_min = cms.double(0.1),
            Pt_fixed = cms.bool(False),
            Pt_min = cms.double(3.0),
            Rescale_Dz = cms.double(4.0),
            Rescale_eta = cms.double(3.0),
            Rescale_phi = cms.double(3.0),
            UseVertex = cms.bool(False),
            Z_fixed = cms.bool(False),
            beamSpot = cms.InputTag("hltOnlineBeamSpot"),
            input = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
            maxRegions = cms.int32(2),
            precise = cms.bool(True),
            vertexCollection = cms.InputTag("pixelVertices")
        ),
        PCut = cms.double(2.5),
        PtCut = cms.double(1.0),
        RefitRPCHits = cms.bool(True),
        ScaleTECxFactor = cms.double(-1.0),
        ScaleTECyFactor = cms.double(-1.0),
        TrackTransformer = cms.PSet(
            DoPredictionsOnly = cms.bool(False),
            Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            Propagator = cms.string('hltESPSmartPropagatorAny'),
            RefitDirection = cms.string('insideOut'),
            RefitRPCHits = cms.bool(True),
            Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
            TrackerRecHitBuilder = cms.string('hltESPTTRHBWithTrackAngle')
        ),
        TrackerPropagator = cms.string('SteppingHelixPropagatorAny'),
        TrackerRecHitBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        tkTrajBeamSpot = cms.InputTag("hltOnlineBeamSpot"),
        tkTrajLabel = cms.InputTag("hltIterL3OIMuonTrackSelectionHighPurity"),
        tkTrajMaxChi2 = cms.double(9999.0),
        tkTrajMaxDXYBeamSpot = cms.double(9999.0),
        tkTrajUseVertex = cms.bool(False),
        tkTrajVertex = cms.InputTag("Notused")
    ),
    MuonCollectionLabel = cms.InputTag("hltL2Muons","UpdatedAtVtx"),
    ServiceParameters = cms.PSet(
        Propagators = cms.untracked.vstring(
            'hltESPSmartPropagatorAny', 
            'SteppingHelixPropagatorAny', 
            'hltESPSmartPropagator', 
            'hltESPSteppingHelixPropagatorOpposite'
        ),
        RPCLayers = cms.bool(True),
        UseMuonNavigation = cms.untracked.bool(True)
    ),
    TrackLoaderParameters = cms.PSet(
        DoSmoothing = cms.bool(True),
        MuonSeededTracksInstance = cms.untracked.string('L2Seeded'),
        MuonUpdatorAtVertexParameters = cms.PSet(
            BeamSpotPositionErrors = cms.vdouble(0.1, 0.1, 5.3),
            MaxChi2 = cms.double(1000000.0),
            Propagator = cms.string('hltESPSteppingHelixPropagatorOpposite')
        ),
        PutTkTrackIntoEvent = cms.untracked.bool(False),
        SmoothTkTrack = cms.untracked.bool(False),
        Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
        TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        VertexConstraint = cms.bool(False),
        beamSpot = cms.InputTag("hltOnlineBeamSpot")
    )
)


process.hltLightPFTracks = cms.EDProducer("LightPFTrackProducer",
    TkColList = cms.VInputTag("hltPFMuonMerging"),
    TrackQuality = cms.string('none'),
    UseQuality = cms.bool(False)
)


process.hltMergedTracks = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltTripletRecoveryMerged", "hltDoubletRecoveryPFlowTrackSelectionHighPurity"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltTripletRecoveryMerged", "hltDoubletRecoveryPFlowTrackSelectionHighPurity"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltMuonCSCDigis = cms.EDProducer("CSCDCCUnpacker",
    Debug = cms.untracked.bool(False),
    ErrorMask = cms.uint32(0),
    ExaminerMask = cms.uint32(535557110),
    FormatedEventDump = cms.untracked.bool(False),
    InputObjects = cms.InputTag("rawDataCollector"),
    PrintEventNumber = cms.untracked.bool(False),
    SuppressZeroLCT = cms.untracked.bool(True),
    UnpackStatusDigis = cms.bool(False),
    UseExaminer = cms.bool(True),
    UseFormatStatus = cms.bool(True),
    UseSelectiveUnpacking = cms.bool(True),
    VisualFEDInspect = cms.untracked.bool(False),
    VisualFEDShort = cms.untracked.bool(False),
    runDQM = cms.untracked.bool(False)
)


process.hltMuonDTDigis = cms.EDProducer("DTUnpackingModule",
    dataType = cms.string('DDU'),
    dqmOnly = cms.bool(False),
    inputLabel = cms.InputTag("rawDataCollector"),
    maxFEDid = cms.untracked.int32(779),
    minFEDid = cms.untracked.int32(770),
    readOutParameters = cms.PSet(
        debug = cms.untracked.bool(False),
        localDAQ = cms.untracked.bool(False),
        performDataIntegrityMonitor = cms.untracked.bool(False),
        rosParameters = cms.PSet(
            debug = cms.untracked.bool(False),
            localDAQ = cms.untracked.bool(False),
            performDataIntegrityMonitor = cms.untracked.bool(False),
            readDDUIDfromDDU = cms.untracked.bool(True),
            readingDDU = cms.untracked.bool(True),
            writeSC = cms.untracked.bool(True)
        )
    ),
    useStandardFEDid = cms.bool(True)
)


process.hltMuonLinks = cms.EDProducer("MuonLinksProducerForHLT",
    InclusiveTrackerTrackCollection = cms.InputTag("hltPFMuonMerging"),
    LinkCollection = cms.InputTag("hltL3MuonsIterL3Links"),
    pMin = cms.double(2.5),
    ptMin = cms.double(2.5),
    shareHitFraction = cms.double(0.8)
)


process.hltMuonRPCDigis = cms.EDProducer("RPCUnpackingModule",
    InputLabel = cms.InputTag("rawDataCollector"),
    doSynchro = cms.bool(False)
)


process.hltMuons = cms.EDProducer("MuonIdProducer",
    CaloExtractorPSet = cms.PSet(
        CenterConeOnCalIntersection = cms.bool(False),
        ComponentName = cms.string('CaloExtractorByAssociator'),
        DR_Max = cms.double(1.0),
        DR_Veto_E = cms.double(0.07),
        DR_Veto_H = cms.double(0.1),
        DR_Veto_HO = cms.double(0.1),
        DepositInstanceLabels = cms.vstring(
            'ecal', 
            'hcal', 
            'ho'
        ),
        DepositLabel = cms.untracked.string('Cal'),
        NoiseTow_EB = cms.double(0.04),
        NoiseTow_EE = cms.double(0.15),
        Noise_EB = cms.double(0.025),
        Noise_EE = cms.double(0.1),
        Noise_HB = cms.double(0.2),
        Noise_HE = cms.double(0.2),
        Noise_HO = cms.double(0.2),
        PrintTimeReport = cms.untracked.bool(False),
        PropagatorName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(False),
            UseMuonNavigation = cms.untracked.bool(False)
        ),
        Threshold_E = cms.double(0.2),
        Threshold_H = cms.double(0.5),
        Threshold_HO = cms.double(0.5),
        TrackAssociatorParameters = cms.PSet(
            CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
            CaloTowerCollectionLabel = cms.InputTag("hltTowerMakerForAll"),
            DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
            EBRecHitCollectionLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
            EERecHitCollectionLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
            HBHERecHitCollectionLabel = cms.InputTag("hltHbhereco"),
            HORecHitCollectionLabel = cms.InputTag("hltHoreco"),
            accountForTrajectoryChangeCalo = cms.bool(False),
            dREcal = cms.double(1.0),
            dREcalPreselection = cms.double(1.0),
            dRHcal = cms.double(1.0),
            dRHcalPreselection = cms.double(1.0),
            dRMuon = cms.double(9999.0),
            dRMuonPreselection = cms.double(0.2),
            dRPreshowerPreselection = cms.double(0.2),
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            propagateAllDirections = cms.bool(True),
            trajectoryUncertaintyTolerance = cms.double(-1.0),
            truthMatch = cms.bool(False),
            useCalo = cms.bool(True),
            useEcal = cms.bool(False),
            useHO = cms.bool(False),
            useHcal = cms.bool(False),
            useMuon = cms.bool(False),
            usePreshower = cms.bool(False)
        ),
        UseRecHitsFlag = cms.bool(False)
    ),
    JetExtractorPSet = cms.PSet(
        ComponentName = cms.string('JetExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.1),
        ExcludeMuonVeto = cms.bool(True),
        JetCollectionLabel = cms.InputTag("hltAK4CaloJetsPFEt5"),
        PrintTimeReport = cms.untracked.bool(False),
        PropagatorName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(False),
            UseMuonNavigation = cms.untracked.bool(False)
        ),
        Threshold = cms.double(5.0),
        TrackAssociatorParameters = cms.PSet(
            CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
            CaloTowerCollectionLabel = cms.InputTag("hltTowerMakerForAll"),
            DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
            EBRecHitCollectionLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
            EERecHitCollectionLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
            HBHERecHitCollectionLabel = cms.InputTag("hltHbhereco"),
            HORecHitCollectionLabel = cms.InputTag("hltHoreco"),
            accountForTrajectoryChangeCalo = cms.bool(False),
            dREcal = cms.double(0.5),
            dREcalPreselection = cms.double(0.5),
            dRHcal = cms.double(0.5),
            dRHcalPreselection = cms.double(0.5),
            dRMuon = cms.double(9999.0),
            dRMuonPreselection = cms.double(0.2),
            dRPreshowerPreselection = cms.double(0.2),
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            propagateAllDirections = cms.bool(True),
            trajectoryUncertaintyTolerance = cms.double(-1.0),
            truthMatch = cms.bool(False),
            useCalo = cms.bool(True),
            useEcal = cms.bool(False),
            useHO = cms.bool(False),
            useHcal = cms.bool(False),
            useMuon = cms.bool(False),
            usePreshower = cms.bool(False)
        )
    ),
    MuonCaloCompatibility = cms.PSet(
        MuonTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root'),
        PionTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root'),
        allSiPMHO = cms.bool(False),
        delta_eta = cms.double(0.02),
        delta_phi = cms.double(0.02)
    ),
    TimingFillerParameters = cms.PSet(
        CSCTimingParameters = cms.PSet(
            CSCStripError = cms.double(7.0),
            CSCStripTimeOffset = cms.double(0.0),
            CSCTimeOffset = cms.double(0.0),
            CSCWireError = cms.double(8.6),
            CSCWireTimeOffset = cms.double(0.0),
            CSCsegments = cms.InputTag("hltCscSegments"),
            MatchParameters = cms.PSet(
                CSCsegments = cms.InputTag("hltCscSegments"),
                DTradius = cms.double(0.01),
                DTsegments = cms.InputTag("hltDt4DSegments"),
                TightMatchCSC = cms.bool(True),
                TightMatchDT = cms.bool(False)
            ),
            PruneCut = cms.double(100.0),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                RPCLayers = cms.bool(True)
            ),
            UseStripTime = cms.bool(True),
            UseWireTime = cms.bool(True),
            debug = cms.bool(False)
        ),
        DTTimingParameters = cms.PSet(
            DTTimeOffset = cms.double(2.7),
            DTsegments = cms.InputTag("hltDt4DSegments"),
            DoWireCorr = cms.bool(False),
            DropTheta = cms.bool(True),
            HitError = cms.double(6.0),
            HitsMin = cms.int32(5),
            MatchParameters = cms.PSet(
                CSCsegments = cms.InputTag("hltCscSegments"),
                DTradius = cms.double(0.01),
                DTsegments = cms.InputTag("hltDt4DSegments"),
                TightMatchCSC = cms.bool(True),
                TightMatchDT = cms.bool(False)
            ),
            PruneCut = cms.double(10000.0),
            RequireBothProjections = cms.bool(False),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                RPCLayers = cms.bool(True)
            ),
            UseSegmentT0 = cms.bool(False),
            debug = cms.bool(False)
        ),
        EcalEnergyCut = cms.double(0.4),
        ErrorCSC = cms.double(7.4),
        ErrorDT = cms.double(6.0),
        ErrorEB = cms.double(2.085),
        ErrorEE = cms.double(6.95),
        UseCSC = cms.bool(True),
        UseDT = cms.bool(True),
        UseECAL = cms.bool(True)
    ),
    TrackAssociatorParameters = cms.PSet(
        CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
        CaloTowerCollectionLabel = cms.InputTag("hltTowerMakerForAll"),
        DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
        EBRecHitCollectionLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
        EERecHitCollectionLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
        HBHERecHitCollectionLabel = cms.InputTag("hltHbhereco"),
        HORecHitCollectionLabel = cms.InputTag("hltHoreco"),
        accountForTrajectoryChangeCalo = cms.bool(False),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        dRHcal = cms.double(9999.0),
        dRHcalPreselection = cms.double(0.2),
        dRMuon = cms.double(9999.0),
        dRMuonPreselection = cms.double(0.2),
        dRPreshowerPreselection = cms.double(0.2),
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        propagateAllDirections = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        truthMatch = cms.bool(False),
        useCalo = cms.bool(False),
        useEcal = cms.bool(True),
        useHO = cms.bool(True),
        useHcal = cms.bool(True),
        useMuon = cms.bool(True),
        usePreshower = cms.bool(False)
    ),
    TrackExtractorPSet = cms.PSet(
        BeamSpotLabel = cms.InputTag("hltOnlineBeamSpot"),
        BeamlineOption = cms.string('BeamSpotFromEvent'),
        Chi2Ndof_Max = cms.double(1e+64),
        Chi2Prob_Min = cms.double(-1.0),
        ComponentName = cms.string('TrackExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.01),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(0.1),
        Diff_z = cms.double(0.2),
        NHits_Min = cms.uint32(0),
        Pt_Min = cms.double(-1.0),
        inputTrackCollection = cms.InputTag("hltPFMuonMerging")
    ),
    TrackerKinkFinderParameters = cms.PSet(
        diagonalOnly = cms.bool(False),
        usePosition = cms.bool(False)
    ),
    addExtraSoftMuons = cms.bool(False),
    arbitrateTrackerMuons = cms.bool(False),
    arbitrationCleanerOptions = cms.PSet(
        ClusterDPhi = cms.double(0.6),
        ClusterDTheta = cms.double(0.02),
        Clustering = cms.bool(True),
        ME1a = cms.bool(True),
        Overlap = cms.bool(True),
        OverlapDPhi = cms.double(0.0786),
        OverlapDTheta = cms.double(0.02)
    ),
    debugWithTruthMatching = cms.bool(False),
    ecalDepositName = cms.string('ecal'),
    fillCaloCompatibility = cms.bool(True),
    fillEnergy = cms.bool(True),
    fillGlobalTrackQuality = cms.bool(False),
    fillGlobalTrackRefits = cms.bool(False),
    fillIsolation = cms.bool(True),
    fillMatching = cms.bool(True),
    fillTrackerKink = cms.bool(False),
    globalTrackQualityInputTag = cms.InputTag("glbTrackQual"),
    hcalDepositName = cms.string('hcal'),
    hoDepositName = cms.string('ho'),
    inputCollectionLabels = cms.VInputTag("hltPFMuonMerging", "hltMuonLinks", "hltL2Muons"),
    inputCollectionTypes = cms.vstring(
        'inner tracks', 
        'links', 
        'outer tracks'
    ),
    jetDepositName = cms.string('jets'),
    maxAbsDx = cms.double(3.0),
    maxAbsDy = cms.double(9999.0),
    maxAbsEta = cms.double(3.0),
    maxAbsPullX = cms.double(4.0),
    maxAbsPullY = cms.double(9999.0),
    minCaloCompatibility = cms.double(0.6),
    minNumberOfMatches = cms.int32(1),
    minP = cms.double(10.0),
    minPCaloMuon = cms.double(1000000000.0),
    minPt = cms.double(10.0),
    ptThresholdToFillCandidateP4WithGlobalFit = cms.double(200.0),
    runArbitrationCleaner = cms.bool(False),
    sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double(2.0),
    trackDepositName = cms.string('tracker'),
    writeIsoDeposits = cms.bool(False)
)


process.hltOnlineBeamSpot = cms.EDProducer("BeamSpotOnlineProducer",
    changeToCMSCoordinates = cms.bool(False),
    gtEvmLabel = cms.InputTag(""),
    maxRadius = cms.double(2.0),
    maxZ = cms.double(40.0),
    setSigmaZ = cms.double(0.0),
    src = cms.InputTag("hltScalersRawToDigi")
)


process.hltPFJetForBtag = cms.EDProducer("HLTPFJetCollectionProducer",
    HLTObject = cms.InputTag("hltPFJetForBtagSelector"),
    TriggerTypes = cms.vint32(86)
)


process.hltPFJetForBtagAK8 = cms.EDProducer("HLTPFJetCollectionProducer",
    HLTObject = cms.InputTag("hltPFJetForBtagSelectorAK8"),
    TriggerTypes = cms.vint32(86)
)


process.hltPFJetsCorrectedMatchedToCaloJets350AK8 = cms.EDProducer("PFJetsMatchedToFilteredCaloJetsProducer",
    CaloJetFilter = cms.InputTag("hltSingleCaloJet350AK8"),
    DeltaR = cms.double(0.5),
    PFJetSrc = cms.InputTag("hltAK8PFJetsCorrected"),
    TriggerType = cms.int32(85)
)


process.hltPFJetsSix32ForHt = cms.EDProducer("HLTPFJetCollectionProducer",
    HLTObject = cms.InputTag("hltPFJetFilterSixC32"),
    TriggerTypes = cms.vint32(86)
)


process.hltPFMuonMerging = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIterL3MuonAndMuonFromL1Merged", "hltMergedTracks"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIterL3MuonAndMuonFromL1Merged", "hltMergedTracks"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltParticleFlow = cms.EDProducer("PFProducer",
    GedElectronValueMap = cms.InputTag("gedGsfElectronsTmp"),
    GedPhotonValueMap = cms.InputTag("tmpGedPhotons","valMapPFEgammaCandToPhoton"),
    PFEGammaCandidates = cms.InputTag("particleFlowEGamma"),
    X0_Map = cms.string('RecoParticleFlow/PFProducer/data/allX0histos.root'),
    algoType = cms.uint32(0),
    blocks = cms.InputTag("hltParticleFlowBlock"),
    calibHF_a_EMHAD = cms.vdouble(
        1.42215, 1.00496, 0.68961, 0.81656, 0.98504, 
        0.98504, 1.00802, 1.0593, 1.4576, 1.4576
    ),
    calibHF_a_EMonly = cms.vdouble(
        0.96945, 0.96701, 0.76309, 0.82268, 0.87583, 
        0.89718, 0.98674, 1.4681, 1.458, 1.458
    ),
    calibHF_b_EMHAD = cms.vdouble(
        1.27541, 0.85361, 0.86333, 0.89091, 0.94348, 
        0.94348, 0.9437, 1.0034, 1.0444, 1.0444
    ),
    calibHF_b_HADonly = cms.vdouble(
        1.27541, 0.85361, 0.86333, 0.89091, 0.94348, 
        0.94348, 0.9437, 1.0034, 1.0444, 1.0444
    ),
    calibHF_eta_step = cms.vdouble(
        0.0, 2.9, 3.0, 3.2, 4.2, 
        4.4, 4.6, 4.8, 5.2, 5.4
    ),
    calibHF_use = cms.bool(False),
    calibPFSCEle_Fbrem_barrel = cms.vdouble(
        0.6, 6.0, -0.0255975, 0.0576727, 0.975442, 
        -0.000546394, 1.26147, 25.0, -0.02025, 0.04537, 
        0.9728, -0.0008962, 1.172
    ),
    calibPFSCEle_Fbrem_endcap = cms.vdouble(
        0.9, 6.5, -0.0692932, 0.101776, 0.995338, 
        -0.00236548, 0.874998, 1.653, -0.0750184, 0.147, 
        0.923165, 0.000474665, 1.10782
    ),
    calibPFSCEle_barrel = cms.vdouble(
        1.004, -1.536, 22.88, -1.467, 0.3555, 
        0.6227, 14.65, 2051.0, 25.0, 0.9932, 
        -0.5444, 0.0, 0.5438, 0.7109, 7.645, 
        0.2904, 0.0
    ),
    calibPFSCEle_endcap = cms.vdouble(
        1.153, -16.5975, 5.668, -0.1772, 16.22, 
        7.326, 0.0483, -4.068, 9.406
    ),
    calibrationsLabel = cms.string('HLT'),
    cleanedHF = cms.VInputTag("hltParticleFlowRecHitHF:Cleaned", "hltParticleFlowClusterHF:Cleaned"),
    coneEcalIsoForEgammaSC = cms.double(0.3),
    coneTrackIsoForEgammaSC = cms.double(0.3),
    cosmicRejectionDistance = cms.double(1.0),
    debug = cms.untracked.bool(False),
    dptRel_DispVtx = cms.double(10.0),
    dzPV = cms.double(0.2),
    egammaElectrons = cms.InputTag(""),
    electron_iso_combIso_barrel = cms.double(10.0),
    electron_iso_combIso_endcap = cms.double(10.0),
    electron_iso_mva_barrel = cms.double(-0.1875),
    electron_iso_mva_endcap = cms.double(-0.1075),
    electron_iso_pt = cms.double(10.0),
    electron_missinghits = cms.uint32(1),
    electron_noniso_mvaCut = cms.double(-0.1),
    electron_protectionsForJetMET = cms.PSet(
        maxDPhiIN = cms.double(0.1),
        maxE = cms.double(50.0),
        maxEcalEOverPRes = cms.double(0.2),
        maxEcalEOverP_1 = cms.double(0.5),
        maxEcalEOverP_2 = cms.double(0.2),
        maxEeleOverPout = cms.double(0.2),
        maxEeleOverPoutRes = cms.double(0.5),
        maxEleHcalEOverEcalE = cms.double(0.1),
        maxHcalE = cms.double(10.0),
        maxHcalEOverEcalE = cms.double(0.1),
        maxHcalEOverP = cms.double(1.0),
        maxNtracks = cms.double(3.0),
        maxTrackPOverEele = cms.double(1.0)
    ),
    eventFactorForCosmics = cms.double(10.0),
    eventFractionForCleaning = cms.double(0.5),
    eventFractionForRejection = cms.double(0.8),
    factors_45 = cms.vdouble(10.0, 100.0),
    iCfgCandConnector = cms.PSet(
        bCalibPrimary = cms.bool(False),
        bCalibSecondary = cms.bool(False),
        bCorrect = cms.bool(False),
        nuclCalibFactors = cms.vdouble(0.8, 0.15, 0.5, 0.5, 0.05)
    ),
    isolatedElectronID_mvaWeightFile = cms.string('RecoEgamma/ElectronIdentification/data/TMVA_BDTSimpleCat_17Feb2011.weights.xml'),
    maxDPtOPt = cms.double(1.0),
    maxDeltaPhiPt = cms.double(7.0),
    maxSignificance = cms.double(2.5),
    metFactorForCleaning = cms.double(4.0),
    metFactorForFakes = cms.double(4.0),
    metFactorForHighEta = cms.double(25.0),
    metFactorForRejection = cms.double(4.0),
    metSignificanceForCleaning = cms.double(3.0),
    metSignificanceForRejection = cms.double(4.0),
    minDeltaMet = cms.double(0.4),
    minEnergyForPunchThrough = cms.double(100.0),
    minHFCleaningPt = cms.double(5.0),
    minMomentumForPunchThrough = cms.double(100.0),
    minPixelHits = cms.int32(1),
    minPtForPostCleaning = cms.double(20.0),
    minSignificance = cms.double(2.5),
    minSignificanceReduction = cms.double(1.4),
    minTrackerHits = cms.int32(8),
    muon_ECAL = cms.vdouble(0.5, 0.5),
    muon_HCAL = cms.vdouble(3.0, 3.0),
    muon_HO = cms.vdouble(0.9, 0.9),
    muons = cms.InputTag("hltMuons"),
    nTrackIsoForEgammaSC = cms.uint32(2),
    nsigma_TRACK = cms.double(1.0),
    pf_GlobC_mvaWeightFile = cms.string('RecoParticleFlow/PFProducer/data/TMVARegression_BDTG_PFGlobalCorr_14Dec2011.root'),
    pf_Res_mvaWeightFile = cms.string('RecoParticleFlow/PFProducer/data/TMVARegression_BDTG_PFRes_14Dec2011.root'),
    pf_convID_mvaWeightFile = cms.string('RecoParticleFlow/PFProducer/data/MVAnalysis_BDT.weights_pfConversionAug0411.txt'),
    pf_conv_mvaCut = cms.double(0.0),
    pf_electronID_crackCorrection = cms.bool(False),
    pf_electronID_mvaWeightFile = cms.string('RecoParticleFlow/PFProducer/data/MVAnalysis_BDT.weights_PfElectrons23Jan_IntToFloat.txt'),
    pf_electron_mvaCut = cms.double(-0.1),
    pf_electron_output_col = cms.string('electrons'),
    pf_locC_mvaWeightFile = cms.string('RecoParticleFlow/PFProducer/data/TMVARegression_BDTG_PFClusterLCorr_14Dec2011.root'),
    pf_nsigma_ECAL = cms.double(0.0),
    pf_nsigma_HCAL = cms.double(1.0),
    photon_HoE = cms.double(0.05),
    photon_MinEt = cms.double(10.0),
    photon_SigmaiEtaiEta_barrel = cms.double(0.0125),
    photon_SigmaiEtaiEta_endcap = cms.double(0.034),
    photon_combIso = cms.double(10.0),
    photon_protectionsForJetMET = cms.PSet(
        sumPtTrackIso = cms.double(2.0),
        sumPtTrackIsoSlope = cms.double(0.001)
    ),
    postHFCleaning = cms.bool(False),
    postMuonCleaning = cms.bool(True),
    ptErrorScale = cms.double(8.0),
    ptFactorForHighEta = cms.double(2.0),
    pt_Error = cms.double(1.0),
    punchThroughFactor = cms.double(3.0),
    punchThroughMETFactor = cms.double(4.0),
    rejectTracks_Bad = cms.bool(False),
    rejectTracks_Step45 = cms.bool(False),
    sumEtEcalIsoForEgammaSC_barrel = cms.double(1.0),
    sumEtEcalIsoForEgammaSC_endcap = cms.double(2.0),
    sumPtTrackIsoForEgammaSC_barrel = cms.double(4.0),
    sumPtTrackIsoForEgammaSC_endcap = cms.double(4.0),
    sumPtTrackIsoForPhoton = cms.double(-1.0),
    sumPtTrackIsoSlopeForPhoton = cms.double(-1.0),
    trackQuality = cms.string('highPurity'),
    useCalibrationsFromDB = cms.bool(True),
    useEGammaElectrons = cms.bool(False),
    useEGammaFilters = cms.bool(False),
    useEGammaSupercluster = cms.bool(False),
    useHO = cms.bool(False),
    usePFConversions = cms.bool(False),
    usePFDecays = cms.bool(False),
    usePFElectrons = cms.bool(False),
    usePFNuclearInteractions = cms.bool(False),
    usePFPhotons = cms.bool(False),
    usePFSCEleCalib = cms.bool(True),
    usePhotonReg = cms.bool(False),
    useProtectionsForJetMET = cms.bool(True),
    useRegressionFromDB = cms.bool(False),
    useVerticesForNeutral = cms.bool(True),
    verbose = cms.untracked.bool(False),
    vertexCollection = cms.InputTag("hltPixelVertices")
)


process.hltParticleFlowBlock = cms.EDProducer("PFBlockProducer",
    debug = cms.untracked.bool(False),
    elementImporters = cms.VPSet(
        cms.PSet(
            DPtOverPtCuts_byTrackAlgo = cms.vdouble(
                0.5, 0.5, 0.5, 0.5, 0.5, 
                0.5
            ),
            NHitCuts_byTrackAlgo = cms.vuint32(
                3, 3, 3, 3, 3, 
                3
            ),
            importerName = cms.string('GeneralTracksImporter'),
            muonSrc = cms.InputTag("hltMuons"),
            source = cms.InputTag("hltLightPFTracks"),
            useIterativeTracking = cms.bool(False)
        ), 
        cms.PSet(
            BCtoPFCMap = cms.InputTag(""),
            importerName = cms.string('ECALClusterImporter'),
            source = cms.InputTag("hltParticleFlowClusterECALUnseeded")
        ), 
        cms.PSet(
            importerName = cms.string('GenericClusterImporter'),
            source = cms.InputTag("hltParticleFlowClusterHCAL")
        ), 
        cms.PSet(
            importerName = cms.string('GenericClusterImporter'),
            source = cms.InputTag("hltParticleFlowClusterHF")
        ), 
        cms.PSet(
            importerName = cms.string('GenericClusterImporter'),
            source = cms.InputTag("hltParticleFlowClusterPSUnseeded")
        )
    ),
    linkDefinitions = cms.VPSet(
        cms.PSet(
            linkType = cms.string('PS1:ECAL'),
            linkerName = cms.string('PreshowerAndECALLinker'),
            useKDTree = cms.bool(True)
        ), 
        cms.PSet(
            linkType = cms.string('PS2:ECAL'),
            linkerName = cms.string('PreshowerAndECALLinker'),
            useKDTree = cms.bool(True)
        ), 
        cms.PSet(
            linkType = cms.string('TRACK:ECAL'),
            linkerName = cms.string('TrackAndECALLinker'),
            useKDTree = cms.bool(True)
        ), 
        cms.PSet(
            linkType = cms.string('TRACK:HCAL'),
            linkerName = cms.string('TrackAndHCALLinker'),
            useKDTree = cms.bool(True)
        ), 
        cms.PSet(
            linkType = cms.string('ECAL:HCAL'),
            linkerName = cms.string('ECALAndHCALLinker'),
            useKDTree = cms.bool(False)
        ), 
        cms.PSet(
            linkType = cms.string('HFEM:HFHAD'),
            linkerName = cms.string('HFEMAndHFHADLinker'),
            useKDTree = cms.bool(False)
        )
    ),
    verbose = cms.untracked.bool(False)
)


process.hltParticleFlowClusterECALUncorrectedUnseeded = cms.EDProducer("PFClusterProducer",
    energyCorrector = cms.PSet(

    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string('Basic2DGenericTopoClusterizer'),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                gatheringThreshold = cms.double(0.08),
                gatheringThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                gatheringThreshold = cms.double(0.3),
                gatheringThresholdPt = cms.double(0.0)
            )
        ),
        useCornerCells = cms.bool(True)
    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowClusterizer'),
        allCellsPositionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.08),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1),
            timeResolutionCalcBarrel = cms.PSet(
                constantTerm = cms.double(0.428192),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0510871),
                noiseTerm = cms.double(1.10889),
                noiseTermLowE = cms.double(1.31883),
                threshHighE = cms.double(5.0),
                threshLowE = cms.double(0.5)
            ),
            timeResolutionCalcEndcap = cms.PSet(
                constantTerm = cms.double(0.0),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(5.72489999999),
                noiseTermLowE = cms.double(6.92683000001),
                threshHighE = cms.double(10.0),
                threshLowE = cms.double(1.0)
            )
        ),
        excludeOtherSeeds = cms.bool(True),
        maxIterations = cms.uint32(50),
        minFracTot = cms.double(1e-20),
        minFractionToKeep = cms.double(1e-07),
        positionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.08),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(9),
            timeResolutionCalcBarrel = cms.PSet(
                constantTerm = cms.double(0.428192),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0510871),
                noiseTerm = cms.double(1.10889),
                noiseTermLowE = cms.double(1.31883),
                threshHighE = cms.double(5.0),
                threshLowE = cms.double(0.5)
            ),
            timeResolutionCalcEndcap = cms.PSet(
                constantTerm = cms.double(0.0),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(5.72489999999),
                noiseTermLowE = cms.double(6.92683000001),
                threshHighE = cms.double(10.0),
                threshLowE = cms.double(1.0)
            )
        ),
        positionCalcForConvergence = cms.PSet(
            T0_EB = cms.double(7.4),
            T0_EE = cms.double(3.1),
            T0_ES = cms.double(1.2),
            W0 = cms.double(4.2),
            X0 = cms.double(0.89),
            algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
            minAllowedNormalization = cms.double(0.0),
            minFractionInCalc = cms.double(0.0)
        ),
        recHitEnergyNorms = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                recHitEnergyNorm = cms.double(0.08)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                recHitEnergyNorm = cms.double(0.3)
            )
        ),
        showerSigma = cms.double(1.5),
        stoppingTolerance = cms.double(1e-08)
    ),
    positionReCalc = cms.PSet(
        T0_EB = cms.double(7.4),
        T0_EE = cms.double(3.1),
        T0_ES = cms.double(1.2),
        W0 = cms.double(4.2),
        X0 = cms.double(0.89),
        algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
        minAllowedNormalization = cms.double(0.0),
        minFractionInCalc = cms.double(0.0)
    ),
    recHitCleaners = cms.VPSet(),
    recHitsSource = cms.InputTag("hltParticleFlowRecHitECALUnseeded"),
    seedFinder = cms.PSet(
        algoName = cms.string('LocalMaximumSeedFinder'),
        nNeighbours = cms.int32(8),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                seedingThreshold = cms.double(0.6),
                seedingThresholdPt = cms.double(0.15)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                seedingThreshold = cms.double(0.23),
                seedingThresholdPt = cms.double(0.0)
            )
        )
    )
)


process.hltParticleFlowClusterECALUnseeded = cms.EDProducer("CorrectedECALPFClusterProducer",
    energyCorrector = cms.PSet(
        algoName = cms.string('PFClusterEMEnergyCorrector'),
        applyCrackCorrections = cms.bool(False)
    ),
    inputECAL = cms.InputTag("hltParticleFlowClusterECALUncorrectedUnseeded"),
    inputPS = cms.InputTag("hltParticleFlowClusterPSUnseeded"),
    minimumPSEnergy = cms.double(0.0)
)


process.hltParticleFlowClusterHBHE = cms.EDProducer("PFClusterProducer",
    energyCorrector = cms.PSet(

    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string('Basic2DGenericTopoClusterizer'),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('HCAL_BARREL1'),
                gatheringThreshold = cms.double(0.8),
                gatheringThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('HCAL_ENDCAP'),
                gatheringThreshold = cms.double(0.8),
                gatheringThresholdPt = cms.double(0.0)
            )
        ),
        useCornerCells = cms.bool(True)
    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowClusterizer'),
        allCellsPositionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.8),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1)
        ),
        clusterTimeResFromSeed = cms.bool(False),
        excludeOtherSeeds = cms.bool(True),
        maxIterations = cms.uint32(50),
        maxNSigmaTime = cms.double(10.0),
        minChi2Prob = cms.double(0.0),
        minFracTot = cms.double(1e-20),
        minFractionToKeep = cms.double(1e-07),
        positionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.8),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(5)
        ),
        recHitEnergyNorms = cms.VPSet(
            cms.PSet(
                detector = cms.string('HCAL_BARREL1'),
                recHitEnergyNorm = cms.double(0.8)
            ), 
            cms.PSet(
                detector = cms.string('HCAL_ENDCAP'),
                recHitEnergyNorm = cms.double(0.8)
            )
        ),
        showerSigma = cms.double(10.0),
        stoppingTolerance = cms.double(1e-08),
        timeResolutionCalcBarrel = cms.PSet(
            constantTerm = cms.double(2.82),
            constantTermLowE = cms.double(4.24),
            corrTermLowE = cms.double(0.0),
            noiseTerm = cms.double(21.86),
            noiseTermLowE = cms.double(8.0),
            threshHighE = cms.double(15.0),
            threshLowE = cms.double(6.0)
        ),
        timeResolutionCalcEndcap = cms.PSet(
            constantTerm = cms.double(2.82),
            constantTermLowE = cms.double(4.24),
            corrTermLowE = cms.double(0.0),
            noiseTerm = cms.double(21.86),
            noiseTermLowE = cms.double(8.0),
            threshHighE = cms.double(15.0),
            threshLowE = cms.double(6.0)
        ),
        timeSigmaEB = cms.double(10.0),
        timeSigmaEE = cms.double(10.0)
    ),
    positionReCalc = cms.PSet(

    ),
    recHitCleaners = cms.VPSet(),
    recHitsSource = cms.InputTag("hltParticleFlowRecHitHBHE"),
    seedFinder = cms.PSet(
        algoName = cms.string('LocalMaximumSeedFinder'),
        nNeighbours = cms.int32(4),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('HCAL_BARREL1'),
                seedingThreshold = cms.double(1.0),
                seedingThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('HCAL_ENDCAP'),
                seedingThreshold = cms.double(1.1),
                seedingThresholdPt = cms.double(0.0)
            )
        )
    )
)


process.hltParticleFlowClusterHCAL = cms.EDProducer("PFMultiDepthClusterProducer",
    clustersSource = cms.InputTag("hltParticleFlowClusterHBHE"),
    energyCorrector = cms.PSet(

    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('PFMultiDepthClusterizer'),
        allCellsPositionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.8),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1)
        ),
        minFractionToKeep = cms.double(1e-07),
        nSigmaEta = cms.double(2.0),
        nSigmaPhi = cms.double(2.0)
    ),
    positionReCalc = cms.PSet(

    )
)


process.hltParticleFlowClusterHF = cms.EDProducer("PFClusterProducer",
    energyCorrector = cms.PSet(

    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string('Basic2DGenericTopoClusterizer'),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('HF_EM'),
                gatheringThreshold = cms.double(0.8),
                gatheringThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('HF_HAD'),
                gatheringThreshold = cms.double(0.8),
                gatheringThresholdPt = cms.double(0.0)
            )
        ),
        useCornerCells = cms.bool(False)
    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowClusterizer'),
        allCellsPositionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.8),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1)
        ),
        excludeOtherSeeds = cms.bool(True),
        maxIterations = cms.uint32(50),
        minFracTot = cms.double(1e-20),
        minFractionToKeep = cms.double(1e-07),
        positionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.8),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(5)
        ),
        recHitEnergyNorms = cms.VPSet(
            cms.PSet(
                detector = cms.string('HF_EM'),
                recHitEnergyNorm = cms.double(0.8)
            ), 
            cms.PSet(
                detector = cms.string('HF_HAD'),
                recHitEnergyNorm = cms.double(0.8)
            )
        ),
        showerSigma = cms.double(0.0),
        stoppingTolerance = cms.double(1e-08)
    ),
    positionReCalc = cms.PSet(

    ),
    recHitCleaners = cms.VPSet(),
    recHitsSource = cms.InputTag("hltParticleFlowRecHitHF"),
    seedFinder = cms.PSet(
        algoName = cms.string('LocalMaximumSeedFinder'),
        nNeighbours = cms.int32(0),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('HF_EM'),
                seedingThreshold = cms.double(1.4),
                seedingThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('HF_HAD'),
                seedingThreshold = cms.double(1.4),
                seedingThresholdPt = cms.double(0.0)
            )
        )
    )
)


process.hltParticleFlowClusterPSUnseeded = cms.EDProducer("PFClusterProducer",
    energyCorrector = cms.PSet(

    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string('Basic2DGenericTopoClusterizer'),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('PS1'),
                gatheringThreshold = cms.double(6e-05),
                gatheringThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('PS2'),
                gatheringThreshold = cms.double(6e-05),
                gatheringThresholdPt = cms.double(0.0)
            )
        ),
        useCornerCells = cms.bool(False)
    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowClusterizer'),
        excludeOtherSeeds = cms.bool(True),
        maxIterations = cms.uint32(50),
        minFracTot = cms.double(1e-20),
        minFractionToKeep = cms.double(1e-07),
        positionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(6e-05),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1)
        ),
        recHitEnergyNorms = cms.VPSet(
            cms.PSet(
                detector = cms.string('PS1'),
                recHitEnergyNorm = cms.double(6e-05)
            ), 
            cms.PSet(
                detector = cms.string('PS2'),
                recHitEnergyNorm = cms.double(6e-05)
            )
        ),
        showerSigma = cms.double(0.3),
        stoppingTolerance = cms.double(1e-08)
    ),
    positionReCalc = cms.PSet(

    ),
    recHitCleaners = cms.VPSet(),
    recHitsSource = cms.InputTag("hltParticleFlowRecHitPSUnseeded"),
    seedFinder = cms.PSet(
        algoName = cms.string('LocalMaximumSeedFinder'),
        nNeighbours = cms.int32(4),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('PS1'),
                seedingThreshold = cms.double(0.00012),
                seedingThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('PS2'),
                seedingThreshold = cms.double(0.00012),
                seedingThresholdPt = cms.double(0.0)
            )
        )
    )
)


process.hltParticleFlowRecHitECALUnseeded = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
        barrel = cms.PSet(

        ),
        endcap = cms.PSet(

        ),
        name = cms.string('PFRecHitECALNavigator')
    ),
    producers = cms.VPSet(
        cms.PSet(
            name = cms.string('PFEBRecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    name = cms.string('PFRecHitQTestECALMultiThreshold'),
                    thresholds = cms.vdouble(
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        1.25023, 1.25033, 1.25047, 1.25068, 1.25097, 
                        1.25139, 1.25199, 1.25286, 1.2541, 1.25587, 
                        1.25842, 1.26207, 1.26729, 1.27479, 1.28553, 
                        1.30092, 1.32299, 1.35462, 1.39995, 1.46493, 
                        1.55807, 1.69156, 1.88291, 2.15716, 2.55027, 
                        3.11371, 3.92131, 5.07887, 6.73803, 9.11615, 
                        10.0, 10.0, 10.0, 10.0, 10.0, 
                        10.0, 10.0, 10.0, 10.0, 1.25023, 
                        1.25033, 1.25047, 1.25068, 1.25097, 1.25139, 
                        1.25199, 1.25286, 1.2541, 1.25587, 1.25842, 
                        1.26207, 1.26729, 1.27479, 1.28553, 1.30092, 
                        1.32299, 1.35462, 1.39995, 1.46493, 1.55807, 
                        1.69156, 1.88291, 2.15716, 2.55027, 3.11371, 
                        3.92131, 5.07887, 6.73803, 9.11615, 10.0, 
                        10.0, 10.0, 10.0, 10.0, 10.0, 
                        10.0, 10.0, 10.0
                    )
                ), 
                cms.PSet(
                    cleaningThreshold = cms.double(2.0),
                    name = cms.string('PFRecHitQTestECAL'),
                    skipTTRecoveredHits = cms.bool(True),
                    timingCleaning = cms.bool(True),
                    topologicalCleaning = cms.bool(True)
                )
            ),
            srFlags = cms.InputTag("hltEcalDigis"),
            src = cms.InputTag("hltEcalRecHit","EcalRecHitsEB")
        ), 
        cms.PSet(
            name = cms.string('PFEERecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    name = cms.string('PFRecHitQTestECALMultiThreshold'),
                    thresholds = cms.vdouble(
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        0.25, 0.25, 0.25, 0.25, 0.25, 
                        1.25023, 1.25033, 1.25047, 1.25068, 1.25097, 
                        1.25139, 1.25199, 1.25286, 1.2541, 1.25587, 
                        1.25842, 1.26207, 1.26729, 1.27479, 1.28553, 
                        1.30092, 1.32299, 1.35462, 1.39995, 1.46493, 
                        1.55807, 1.69156, 1.88291, 2.15716, 2.55027, 
                        3.11371, 3.92131, 5.07887, 6.73803, 9.11615, 
                        10.0, 10.0, 10.0, 10.0, 10.0, 
                        10.0, 10.0, 10.0, 10.0, 1.25023, 
                        1.25033, 1.25047, 1.25068, 1.25097, 1.25139, 
                        1.25199, 1.25286, 1.2541, 1.25587, 1.25842, 
                        1.26207, 1.26729, 1.27479, 1.28553, 1.30092, 
                        1.32299, 1.35462, 1.39995, 1.46493, 1.55807, 
                        1.69156, 1.88291, 2.15716, 2.55027, 3.11371, 
                        3.92131, 5.07887, 6.73803, 9.11615, 10.0, 
                        10.0, 10.0, 10.0, 10.0, 10.0, 
                        10.0, 10.0, 10.0
                    )
                ), 
                cms.PSet(
                    cleaningThreshold = cms.double(2.0),
                    name = cms.string('PFRecHitQTestECAL'),
                    skipTTRecoveredHits = cms.bool(True),
                    timingCleaning = cms.bool(True),
                    topologicalCleaning = cms.bool(True)
                )
            ),
            srFlags = cms.InputTag("hltEcalDigis"),
            src = cms.InputTag("hltEcalRecHit","EcalRecHitsEE")
        )
    )
)


process.hltParticleFlowRecHitHBHE = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
        name = cms.string('PFRecHitHCALNavigator'),
        sigmaCut = cms.double(4.0),
        timeResolutionCalc = cms.PSet(
            constantTerm = cms.double(1.92),
            constantTermLowE = cms.double(6.0),
            corrTermLowE = cms.double(0.0),
            noiseTerm = cms.double(8.64),
            noiseTermLowE = cms.double(0.0),
            threshHighE = cms.double(8.0),
            threshLowE = cms.double(2.0)
        )
    ),
    producers = cms.VPSet(cms.PSet(
        name = cms.string('PFHBHERecHitCreator'),
        qualityTests = cms.VPSet(
            cms.PSet(
                name = cms.string('PFRecHitQTestThreshold'),
                threshold = cms.double(0.8)
            ), 
            cms.PSet(
                cleaningThresholds = cms.vdouble(0.0),
                flags = cms.vstring('Standard'),
                maxSeverities = cms.vint32(11),
                name = cms.string('PFRecHitQTestHCALChannel')
            )
        ),
        src = cms.InputTag("hltHbhereco")
    ))
)


process.hltParticleFlowRecHitHF = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
        barrel = cms.PSet(

        ),
        endcap = cms.PSet(

        ),
        name = cms.string('PFRecHitHCALNavigator')
    ),
    producers = cms.VPSet(cms.PSet(
        EMDepthCorrection = cms.double(22.0),
        HADDepthCorrection = cms.double(25.0),
        HFCalib29 = cms.double(1.07),
        LongFibre_Cut = cms.double(120.0),
        LongFibre_Fraction = cms.double(0.1),
        ShortFibre_Cut = cms.double(60.0),
        ShortFibre_Fraction = cms.double(0.01),
        name = cms.string('PFHFRecHitCreator'),
        qualityTests = cms.VPSet(
            cms.PSet(
                cleaningThresholds = cms.vdouble(0.0, 120.0, 60.0),
                flags = cms.vstring(
                    'Standard', 
                    'HFLong', 
                    'HFShort'
                ),
                maxSeverities = cms.vint32(11, 9, 9),
                name = cms.string('PFRecHitQTestHCALChannel')
            ), 
            cms.PSet(
                cuts = cms.VPSet(
                    cms.PSet(
                        depth = cms.int32(1),
                        threshold = cms.double(1.2)
                    ), 
                    cms.PSet(
                        depth = cms.int32(2),
                        threshold = cms.double(1.8)
                    )
                ),
                name = cms.string('PFRecHitQTestHCALThresholdVsDepth')
            )
        ),
        src = cms.InputTag("hltHfreco"),
        thresh_HF = cms.double(0.4)
    ))
)


process.hltParticleFlowRecHitPSUnseeded = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
        name = cms.string('PFRecHitPreshowerNavigator')
    ),
    producers = cms.VPSet(cms.PSet(
        name = cms.string('PFPSRecHitCreator'),
        qualityTests = cms.VPSet(cms.PSet(
            name = cms.string('PFRecHitQTestThreshold'),
            threshold = cms.double(7e-06)
        )),
        src = cms.InputTag("hltEcalPreshowerRecHit","EcalRecHitsES")
    ))
)


process.hltPixelLayerQuadruplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4', 
        'BPix1+BPix2+BPix3+FPix1_pos', 
        'BPix1+BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos+FPix2_pos', 
        'BPix1+BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltPixelLayerQuadrupletsRegForBTag = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHitsRegForBTag'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHitsRegForBTag'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4', 
        'BPix1+BPix2+BPix3+FPix1_pos', 
        'BPix1+BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos+FPix2_pos', 
        'BPix1+BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    )
)


process.hltPixelLayerTripletsWithClustersRemoval = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltPixelTripletsClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltPixelTripletsClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3', 
        'BPix2+BPix3+BPix4', 
        'BPix1+BPix3+BPix4', 
        'BPix1+BPix2+BPix4', 
        'BPix2+BPix3+FPix1_pos', 
        'BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos', 
        'BPix1+BPix2+FPix1_neg', 
        'BPix2+FPix1_pos+FPix2_pos', 
        'BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos', 
        'BPix1+FPix1_neg+FPix2_neg', 
        'FPix1_pos+FPix2_pos+FPix3_pos', 
        'FPix1_neg+FPix2_neg+FPix3_neg', 
        'BPix1+BPix3+FPix1_pos', 
        'BPix1+BPix2+FPix2_pos', 
        'BPix1+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix2_neg', 
        'BPix1+FPix2_neg+FPix3_neg', 
        'BPix1+FPix1_neg+FPix3_neg', 
        'BPix1+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_pos+FPix3_pos'
    )
)


process.hltPixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltPixelTracksFilter"),
    Fitter = cms.InputTag("hltPixelTracksFitter"),
    SeedingHitSets = cms.InputTag("hltPixelTracksHitQuadruplets"),
    passLabel = cms.string('')
)


process.hltPixelTracksFilter = cms.EDProducer("PixelTrackFilterByKinematicsProducer",
    chi2 = cms.double(1000.0),
    nSigmaInvPtTolerance = cms.double(0.0),
    nSigmaTipMaxTolerance = cms.double(0.0),
    ptMin = cms.double(0.1),
    tipMax = cms.double(1.0)
)


process.hltPixelTracksFitter = cms.EDProducer("PixelFitterByHelixProjectionsProducer",
    scaleErrorsForBPix1 = cms.bool(False),
    scaleFactor = cms.double(0.65)
)


process.hltPixelTracksFromTriplets = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltPixelTracksFilter"),
    Fitter = cms.InputTag("hltPixelTracksFitter"),
    SeedingHitSets = cms.InputTag("hltPixelTracksHitTriplets"),
    passLabel = cms.string('')
)


process.hltPixelTracksHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag(""),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltPixelLayerQuadruplets"),
    trackingRegions = cms.InputTag("hltPixelTracksTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltPixelTracksHitDoubletsForTriplets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag(""),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltPixelLayerTripletsWithClustersRemoval"),
    trackingRegions = cms.InputTag("hltPixelTracksTrackingRegionsForTriplets"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltPixelTracksHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.2),
    CAThetaCut = cms.double(0.002),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltPixelTracksHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2.0),
        value1 = cms.double(200.0),
        value2 = cms.double(50.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltPixelTracksHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.0),
    CAPhiCut = cms.double(0.2),
    CAThetaCut = cms.double(0.002),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("hltSiPixelClustersCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("hltPixelTracksHitDoubletsForTriplets"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(False),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2.0),
        value1 = cms.double(200.0),
        value2 = cms.double(50.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltPixelTracksMerged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltPixelTracks", "hltPixelTracksFromTriplets"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltPixelTracks", "hltPixelTracksFromTriplets"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltPixelTracksTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        nSigmaZ = cms.double(4.0),
        originRadius = cms.double(0.02),
        precise = cms.bool(True),
        ptMin = cms.double(0.8)
    )
)


process.hltPixelTracksTrackingRegionsForTriplets = cms.EDProducer("PointSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(0.4),
        deltaPhi = cms.double(0.4),
        maxNRegions = cms.int32(100),
        maxNVertices = cms.int32(10),
        measurementTrackerName = cms.InputTag("hltDoubletRecoveryMaskedMeasurementTrackerEvent"),
        mode = cms.string('BeamSpotFixed'),
        nSigmaZBeamSpot = cms.double(3.0),
        nSigmaZVertex = cms.double(3.0),
        originRadius = cms.double(0.1),
        points = cms.PSet(
            eta = cms.vdouble(2.1, -0.8),
            phi = cms.vdouble(1.8, -3.2)
        ),
        precise = cms.bool(True),
        ptMin = cms.double(0.8),
        searchOpt = cms.bool(False),
        vertexCollection = cms.InputTag("none"),
        whereToUseMeasurementTracker = cms.string('never'),
        zErrorBeamSpot = cms.double(15.0),
        zErrorVetex = cms.double(0.1)
    )
)


process.hltPixelTripletsClustersRefRemoval = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('undefQuality'),
    maxChi2 = cms.double(3000.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClusters"),
    stripClusters = cms.InputTag(""),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltPixelTracks")
)


process.hltPixelVertices = cms.EDProducer("PixelVertexProducer",
    Finder = cms.string('DivisiveVertexFinder'),
    Method2 = cms.bool(True),
    NTrkMin = cms.int32(2),
    PVcomparer = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPvClusterComparerForIT')
    ),
    PtMin = cms.double(1.0),
    TrackCollection = cms.InputTag("hltPixelTracksMerged"),
    UseError = cms.bool(True),
    Verbosity = cms.int32(0),
    WtAverage = cms.bool(True),
    ZOffset = cms.double(5.0),
    ZSeparation = cms.double(0.05),
    beamSpot = cms.InputTag("hltOnlineBeamSpot")
)


process.hltRpcRecHits = cms.EDProducer("RPCRecHitProducer",
    deadSource = cms.string('File'),
    deadvecfile = cms.FileInPath('RecoLocalMuon/RPCRecHit/data/RPCDeadVec.dat'),
    maskSource = cms.string('File'),
    maskvecfile = cms.FileInPath('RecoLocalMuon/RPCRecHit/data/RPCMaskVec.dat'),
    recAlgo = cms.string('RPCRecHitStandardAlgo'),
    recAlgoConfig = cms.PSet(

    ),
    rpcDigiLabel = cms.InputTag("hltMuonRPCDigis")
)


process.hltScalersRawToDigi = cms.EDProducer("ScalersRawToDigi",
    scalersInputTag = cms.InputTag("rawDataCollector")
)


process.hltSecondaryVertexTagInfosPF = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("hltOnlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("hltInclusiveMergedVerticesPF"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("hltBLifetimeTagInfosPF"),
    trackSelection = cms.PSet(
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(2),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(2)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(True),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(2.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.79),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.5),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.hltSecondaryVertexTagInfosPFAK8 = cms.EDProducer("SecondaryVertexProducer",
    beamSpotTag = cms.InputTag("hltOnlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("hltInclusiveMergedVerticesPF"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("hltBLifetimeTagInfosPFAK8"),
    trackSelection = cms.PSet(
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(2),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(2)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(True),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(2.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.79),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.5),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.hltSiPixelClusters = cms.EDProducer("SiPixelClusterProducer",
    ChannelThreshold = cms.int32(1000),
    ClusterThreshold = cms.int32(4000),
    ClusterThreshold_L1 = cms.int32(2000),
    MissCalibrate = cms.untracked.bool(True),
    SeedThreshold = cms.int32(1000),
    SplitClusters = cms.bool(False),
    VCaltoElectronGain = cms.int32(47),
    VCaltoElectronGain_L1 = cms.int32(50),
    VCaltoElectronOffset = cms.int32(-60),
    VCaltoElectronOffset_L1 = cms.int32(-670),
    maxNumberOfClusters = cms.int32(40000),
    payloadType = cms.string('HLT'),
    src = cms.InputTag("hltSiPixelDigis")
)


process.hltSiPixelClustersCache = cms.EDProducer("SiPixelClusterShapeCacheProducer",
    onDemand = cms.bool(False),
    src = cms.InputTag("hltSiPixelClusters")
)


process.hltSiPixelClustersRegForBTag = cms.EDProducer("SiPixelClusterProducer",
    ChannelThreshold = cms.int32(1000),
    ClusterThreshold = cms.int32(4000),
    ClusterThreshold_L1 = cms.int32(2000),
    MissCalibrate = cms.untracked.bool(True),
    SeedThreshold = cms.int32(1000),
    SplitClusters = cms.bool(False),
    VCaltoElectronGain = cms.int32(47),
    VCaltoElectronGain_L1 = cms.int32(50),
    VCaltoElectronOffset = cms.int32(-60),
    VCaltoElectronOffset_L1 = cms.int32(-670),
    maxNumberOfClusters = cms.int32(40000),
    payloadType = cms.string('HLT'),
    src = cms.InputTag("hltSiPixelDigisRegForBTag")
)


process.hltSiPixelClustersRegForBTagCache = cms.EDProducer("SiPixelClusterShapeCacheProducer",
    onDemand = cms.bool(False),
    src = cms.InputTag("hltSiPixelClustersRegForBTag")
)


process.hltSiPixelDigis = cms.EDProducer("SiPixelRawToDigi",
    CablingMapLabel = cms.string(''),
    ErrorList = cms.vint32(),
    IncludeErrors = cms.bool(False),
    InputLabel = cms.InputTag("rawDataCollector"),
    Regions = cms.PSet(

    ),
    Timing = cms.untracked.bool(False),
    UsePhase1 = cms.bool(True),
    UsePilotBlade = cms.bool(False),
    UseQualityInfo = cms.bool(False),
    UserErrorList = cms.vint32()
)


process.hltSiPixelDigisRegForBTag = cms.EDProducer("SiPixelRawToDigi",
    CablingMapLabel = cms.string(''),
    ErrorList = cms.vint32(),
    IncludeErrors = cms.bool(False),
    InputLabel = cms.InputTag("rawDataCollector"),
    Regions = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaPhi = cms.vdouble(0.5),
        inputs = cms.VInputTag("hltSelectorCentralJets20L1FastJeta"),
        maxZ = cms.vdouble(24.0)
    ),
    Timing = cms.untracked.bool(False),
    UsePhase1 = cms.bool(True),
    UsePilotBlade = cms.bool(False),
    UseQualityInfo = cms.bool(False),
    UserErrorList = cms.vint32()
)


process.hltSiPixelRecHits = cms.EDProducer("SiPixelRecHitConverter",
    CPE = cms.string('hltESPPixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("hltSiPixelClusters")
)


process.hltSiPixelRecHitsRegForBTag = cms.EDProducer("SiPixelRecHitConverter",
    CPE = cms.string('hltESPPixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("hltSiPixelClustersRegForBTag")
)


process.hltSiStripClusters = cms.EDProducer("MeasurementTrackerEventProducer",
    Phase2TrackerCluster1DProducer = cms.string(''),
    badPixelFEDChannelCollectionLabels = cms.VInputTag(),
    inactivePixelDetectorLabels = cms.VInputTag(),
    inactiveStripDetectorLabels = cms.VInputTag("hltSiStripExcludedFEDListProducer"),
    measurementTracker = cms.string('hltESPMeasurementTracker'),
    pixelCablingMapLabel = cms.string(''),
    pixelClusterProducer = cms.string('hltSiPixelClusters'),
    skipClusters = cms.InputTag(""),
    stripClusterProducer = cms.string('hltSiStripRawToClustersFacility'),
    switchOffPixelsIfEmpty = cms.bool(True)
)


process.hltSiStripClustersRegForBTag = cms.EDProducer("MeasurementTrackerEventProducer",
    Phase2TrackerCluster1DProducer = cms.string(''),
    badPixelFEDChannelCollectionLabels = cms.VInputTag(),
    inactivePixelDetectorLabels = cms.VInputTag(),
    inactiveStripDetectorLabels = cms.VInputTag("hltSiStripExcludedFEDListProducer"),
    measurementTracker = cms.string('hltESPMeasurementTracker'),
    pixelCablingMapLabel = cms.string(''),
    pixelClusterProducer = cms.string('hltSiPixelClustersRegForBTag'),
    skipClusters = cms.InputTag(""),
    stripClusterProducer = cms.string('hltSiStripRawToClustersFacility'),
    switchOffPixelsIfEmpty = cms.bool(True)
)


process.hltSiStripExcludedFEDListProducer = cms.EDProducer("SiStripExcludedFEDListProducer",
    ProductLabel = cms.InputTag("rawDataCollector")
)


process.hltSiStripRawToClustersFacility = cms.EDProducer("SiStripClusterizerFromRaw",
    Algorithms = cms.PSet(
        CommonModeNoiseSubtractionMode = cms.string('Median'),
        PedestalSubtractionFedMode = cms.bool(True),
        SiStripFedZeroSuppressionMode = cms.uint32(4),
        TruncateInSuppressor = cms.bool(True),
        doAPVRestore = cms.bool(False),
        useCMMeanMap = cms.bool(False)
    ),
    Clusterizer = cms.PSet(
        Algorithm = cms.string('ThreeThresholdAlgorithm'),
        ChannelThreshold = cms.double(2.0),
        ClusterThreshold = cms.double(5.0),
        MaxAdjacentBad = cms.uint32(0),
        MaxSequentialBad = cms.uint32(1),
        MaxSequentialHoles = cms.uint32(0),
        QualityLabel = cms.string(''),
        RemoveApvShots = cms.bool(True),
        SeedThreshold = cms.double(3.0),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
        ),
        setDetId = cms.bool(True)
    ),
    DoAPVEmulatorCheck = cms.bool(False),
    ProductLabel = cms.InputTag("rawDataCollector"),
    onDemand = cms.bool(True)
)


process.hltTowerMakerForAll = cms.EDProducer("CaloTowersCreator",
    AllowMissingInputs = cms.bool(False),
    EBGrid = cms.vdouble(),
    EBSumThreshold = cms.double(0.2),
    EBThreshold = cms.double(0.07),
    EBWeight = cms.double(1.0),
    EBWeights = cms.vdouble(),
    EEGrid = cms.vdouble(),
    EESumThreshold = cms.double(0.45),
    EEThreshold = cms.double(0.3),
    EEWeight = cms.double(1.0),
    EEWeights = cms.vdouble(),
    EcalRecHitSeveritiesToBeExcluded = cms.vstring(
        'kTime', 
        'kWeird', 
        'kBad'
    ),
    EcalSeveritiesToBeUsedInBadTowers = cms.vstring(),
    EcutTower = cms.double(-1000.0),
    HBGrid = cms.vdouble(),
    HBThreshold = cms.double(0.7),
    HBWeight = cms.double(1.0),
    HBWeights = cms.vdouble(),
    HEDGrid = cms.vdouble(),
    HEDThreshold = cms.double(0.8),
    HEDWeight = cms.double(1.0),
    HEDWeights = cms.vdouble(),
    HESGrid = cms.vdouble(),
    HESThreshold = cms.double(0.8),
    HESWeight = cms.double(1.0),
    HESWeights = cms.vdouble(),
    HF1Grid = cms.vdouble(),
    HF1Threshold = cms.double(0.5),
    HF1Weight = cms.double(1.0),
    HF1Weights = cms.vdouble(),
    HF2Grid = cms.vdouble(),
    HF2Threshold = cms.double(0.85),
    HF2Weight = cms.double(1.0),
    HF2Weights = cms.vdouble(),
    HOGrid = cms.vdouble(),
    HOThreshold0 = cms.double(3.5),
    HOThresholdMinus1 = cms.double(3.5),
    HOThresholdMinus2 = cms.double(3.5),
    HOThresholdPlus1 = cms.double(3.5),
    HOThresholdPlus2 = cms.double(3.5),
    HOWeight = cms.double(1e-99),
    HOWeights = cms.vdouble(),
    HcalAcceptSeverityLevel = cms.uint32(9),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32(9999),
    HcalPhase = cms.int32(0),
    HcalThreshold = cms.double(-1000.0),
    MomConstrMethod = cms.int32(1),
    MomEBDepth = cms.double(0.3),
    MomEEDepth = cms.double(0.0),
    MomHBDepth = cms.double(0.2),
    MomHEDepth = cms.double(0.4),
    UseEcalRecoveredHits = cms.bool(False),
    UseEtEBTreshold = cms.bool(False),
    UseEtEETreshold = cms.bool(False),
    UseHO = cms.bool(False),
    UseHcalRecoveredHits = cms.bool(False),
    UseRejectedHitsOnly = cms.bool(False),
    UseRejectedRecoveredEcalHits = cms.bool(False),
    UseRejectedRecoveredHcalHits = cms.bool(False),
    UseSymEBTreshold = cms.bool(False),
    UseSymEETreshold = cms.bool(False),
    ecalInputs = cms.VInputTag("hltEcalRecHit:EcalRecHitsEB", "hltEcalRecHit:EcalRecHitsEE"),
    hbheInput = cms.InputTag("hltHbhereco"),
    hfInput = cms.InputTag("hltHfreco"),
    hoInput = cms.InputTag("hltHoreco")
)


process.hltTrackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxTimeSignificance = cms.double(3.5),
    primaryVertices = cms.InputTag("hltVerticesL3"),
    secondaryVertices = cms.InputTag("hltInclusiveSecondaryVertices"),
    sigCut = cms.double(5.0),
    trackMinLayers = cms.int32(4),
    trackMinPixels = cms.int32(1),
    trackMinPt = cms.double(0.4),
    tracks = cms.InputTag("hltIter2MergedForBTag")
)


process.hltTrackVertexArbitratorPF = cms.EDProducer("TrackVertexArbitrator",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3.0),
    fitterTini = cms.double(256.0),
    maxTimeSignificance = cms.double(3.5),
    primaryVertices = cms.InputTag("hltVerticesPFFilter"),
    secondaryVertices = cms.InputTag("hltInclusiveSecondaryVerticesPF"),
    sigCut = cms.double(5.0),
    trackMinLayers = cms.int32(4),
    trackMinPixels = cms.int32(1),
    trackMinPt = cms.double(0.4),
    tracks = cms.InputTag("hltPFMuonMerging")
)


process.hltTriggerSummaryAOD = cms.EDProducer("TriggerSummaryProducerAOD",
    moduleLabelPatternsToMatch = cms.vstring('hlt*'),
    moduleLabelPatternsToSkip = cms.vstring(),
    processName = cms.string('@'),
    throw = cms.bool(True)
)


process.hltTriggerSummaryRAW = cms.EDProducer("TriggerSummaryProducerRAW",
    processName = cms.string('@')
)


process.hltTrimmedPixelVertices = cms.EDProducer("PixelVertexCollectionTrimmer",
    PVcomparer = cms.PSet(
        refToPSet_ = cms.string('HLTPSetPvClusterComparerForIT')
    ),
    fractionSumPt2 = cms.double(0.3),
    maxVtx = cms.uint32(100),
    minSumPt2 = cms.double(0.0),
    src = cms.InputTag("hltPixelVertices")
)


process.hltTripletRecoveryClustersRefRemoval = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(16.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("hltIter2ClustersRefRemoval"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("hltSiPixelClusters"),
    stripClusters = cms.InputTag("hltSiStripRawToClustersFacility"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter2PFlowTrackSelectionHighPurity")
)


process.hltTripletRecoveryMaskedMeasurementTrackerEvent = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
    OnDemand = cms.bool(False),
    clustersToSkip = cms.InputTag("hltTripletRecoveryClustersRefRemoval"),
    src = cms.InputTag("hltSiStripClusters")
)


process.hltTripletRecoveryMerged = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(20.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("hltIter2Merged", "hltTripletRecoveryPFlowTrackSelectionHighPurity"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(False),
    hasSelector = cms.vint32(0, 0),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("hltIter2Merged", "hltTripletRecoveryPFlowTrackSelectionHighPurity"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(False),
        tLists = cms.vint32(0, 1)
    )),
    trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.hltTripletRecoveryPFlowCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("hltTripletRecoveryMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TrajectoryBuilder = cms.string(''),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('HLTIter2GroupedCkfTrajectoryBuilderIT')
    ),
    TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(False),
    doSeedingRegionRebuilding = cms.bool(False),
    maxNSeeds = cms.uint32(100000),
    maxSeedsBeforeCleaning = cms.uint32(1000),
    src = cms.InputTag("hltTripletRecoveryPFlowPixelSeeds"),
    useHitsSplitting = cms.bool(False)
)


process.hltTripletRecoveryPFlowCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltTripletRecovery'),
    Fitter = cms.string('hltESPFittingSmootherIT'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltTripletRecoveryMaskedMeasurementTrackerEvent"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltTripletRecoveryPFlowCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hltTripletRecoveryPFlowPixelClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("hltSiStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(50000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClusters"),
    cut = cms.string(''),
    doClusterCheck = cms.bool(False),
    silentClusterCheck = cms.untracked.bool(False)
)


process.hltTripletRecoveryPFlowPixelHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("hltTripletRecoveryPFlowPixelClusterCheck"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(0),
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("hltTripletRecoveryPixelLayerTriplets"),
    trackingRegions = cms.InputTag("hltTripletRecoveryPFlowPixelTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.hltTripletRecoveryPFlowPixelHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.3),
    CAPhiCut = cms.double(0.1),
    CAThetaCut = cms.double(0.004),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("hltTripletRecoveryPFlowPixelHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(8.0),
        value1 = cms.double(100.0),
        value2 = cms.double(50.0)
    ),
    useBendingCorrection = cms.bool(True)
)


process.hltTripletRecoveryPFlowPixelSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("hltTripletRecoveryPFlowPixelHitTriplets")
)


process.hltTripletRecoveryPFlowPixelTrackingRegions = cms.EDProducer("PointSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(2.1),
        deltaPhi = cms.double(0.55),
        maxNRegions = cms.int32(100),
        maxNVertices = cms.int32(10),
        measurementTrackerName = cms.InputTag("hltTripletRecoveryMaskedMeasurementTrackerEvent"),
        mode = cms.string('VerticesFixed'),
        nSigmaZBeamSpot = cms.double(3.0),
        nSigmaZVertex = cms.double(3.0),
        originRadius = cms.double(0.05),
        points = cms.PSet(
            eta = cms.vdouble(0.4, -0.8, 4.4),
            phi = cms.vdouble(1.7, -3.1, -2.0)
        ),
        precise = cms.bool(True),
        ptMin = cms.double(0.8),
        searchOpt = cms.bool(False),
        vertexCollection = cms.InputTag("hltTrimmedPixelVertices"),
        whereToUseMeasurementTracker = cms.string('ForSiStrips'),
        zErrorBeamSpot = cms.double(15.0),
        zErrorVetex = cms.double(0.1)
    )
)


process.hltTripletRecoveryPFlowTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            dr_exp = cms.vint32(4, 4, 4),
            dr_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dr_par2 = cms.vdouble(3.40282346639e+38, 0.3, 0.3)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 4),
            dz_par1 = cms.vdouble(3.40282346639e+38, 0.4, 0.4),
            dz_par2 = cms.vdouble(3.40282346639e+38, 0.35, 0.35)
        ),
        maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 0, 0),
        minLayers = cms.vint32(3, 3, 3),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltTripletRecoveryPFlowCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltTrimmedPixelVertices")
)


process.hltTripletRecoveryPFlowTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltTripletRecoveryPFlowTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltTripletRecoveryPFlowTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltTripletRecoveryPFlowCtfWithMaterialTracks")
)


process.hltTripletRecoveryPixelLayerTriplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        skipClusters = cms.InputTag("hltTripletRecoveryClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('hltSiPixelRecHits'),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        skipClusters = cms.InputTag("hltTripletRecoveryClustersRefRemoval"),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3', 
        'BPix2+BPix3+BPix4', 
        'BPix1+BPix3+BPix4', 
        'BPix1+BPix2+BPix4', 
        'BPix2+BPix3+FPix1_pos', 
        'BPix2+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix1_pos', 
        'BPix1+BPix2+FPix1_neg', 
        'BPix2+FPix1_pos+FPix2_pos', 
        'BPix2+FPix1_neg+FPix2_neg', 
        'BPix1+FPix1_pos+FPix2_pos', 
        'BPix1+FPix1_neg+FPix2_neg', 
        'FPix1_pos+FPix2_pos+FPix3_pos', 
        'FPix1_neg+FPix2_neg+FPix3_neg', 
        'BPix1+BPix3+FPix1_pos', 
        'BPix1+BPix2+FPix2_pos', 
        'BPix1+BPix3+FPix1_neg', 
        'BPix1+BPix2+FPix2_neg', 
        'BPix1+FPix2_neg+FPix3_neg', 
        'BPix1+FPix1_neg+FPix3_neg', 
        'BPix1+FPix2_pos+FPix3_pos', 
        'BPix1+FPix1_pos+FPix3_pos'
    )
)


process.hltVerticesL3 = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        TkDAClusParameters = cms.PSet(
            Tmin = cms.double(2.4),
            Tpurge = cms.double(2.0),
            Tstop = cms.double(0.5),
            coolingFactor = cms.double(0.6),
            d0CutOff = cms.double(999.0),
            dzCutOff = cms.double(4.0),
            uniquetrkweight = cms.double(0.9),
            use_vdt = cms.untracked.bool(True),
            vertexSize = cms.double(0.15),
            zmerge = cms.double(0.01)
        ),
        algorithm = cms.string('DA_vect')
    ),
    TkFilterParameters = cms.PSet(
        algorithm = cms.string('filter'),
        maxD0Significance = cms.double(999.0),
        maxEta = cms.double(100.0),
        maxNormalizedChi2 = cms.double(20.0),
        minPixelLayersWithHits = cms.int32(2),
        minPt = cms.double(0.0),
        minSiliconLayersWithHits = cms.int32(5),
        trackQuality = cms.string('any')
    ),
    TrackLabel = cms.InputTag("hltIter2MergedForBTag"),
    beamSpotLabel = cms.InputTag("hltOnlineBeamSpot"),
    verbose = cms.untracked.bool(False),
    vertexCollections = cms.VPSet(
        cms.PSet(
            algorithm = cms.string('AdaptiveVertexFitter'),
            chi2cutoff = cms.double(3.0),
            label = cms.string(''),
            maxDistanceToBeam = cms.double(1.0),
            minNdof = cms.double(0.0),
            useBeamConstraint = cms.bool(False)
        ), 
        cms.PSet(
            algorithm = cms.string('AdaptiveVertexFitter'),
            chi2cutoff = cms.double(3.0),
            label = cms.string('WithBS'),
            maxDistanceToBeam = cms.double(1.0),
            minNdof = cms.double(0.0),
            useBeamConstraint = cms.bool(True)
        )
    )
)


process.hltVerticesPF = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        TkDAClusParameters = cms.PSet(
            Tmin = cms.double(2.4),
            Tpurge = cms.double(2.0),
            Tstop = cms.double(0.5),
            coolingFactor = cms.double(0.6),
            d0CutOff = cms.double(999.0),
            dzCutOff = cms.double(4.0),
            uniquetrkweight = cms.double(0.9),
            use_vdt = cms.untracked.bool(True),
            vertexSize = cms.double(0.15),
            zmerge = cms.double(0.01)
        ),
        algorithm = cms.string('DA_vect')
    ),
    TkFilterParameters = cms.PSet(
        algorithm = cms.string('filter'),
        maxD0Significance = cms.double(999.0),
        maxEta = cms.double(100.0),
        maxNormalizedChi2 = cms.double(20.0),
        minPixelLayersWithHits = cms.int32(2),
        minPt = cms.double(0.0),
        minSiliconLayersWithHits = cms.int32(5),
        trackQuality = cms.string('any')
    ),
    TrackLabel = cms.InputTag("hltPFMuonMerging"),
    beamSpotLabel = cms.InputTag("hltOnlineBeamSpot"),
    verbose = cms.untracked.bool(False),
    vertexCollections = cms.VPSet(
        cms.PSet(
            algorithm = cms.string('AdaptiveVertexFitter'),
            chi2cutoff = cms.double(3.0),
            label = cms.string(''),
            maxDistanceToBeam = cms.double(1.0),
            minNdof = cms.double(0.0),
            useBeamConstraint = cms.bool(False)
        ), 
        cms.PSet(
            algorithm = cms.string('AdaptiveVertexFitter'),
            chi2cutoff = cms.double(3.0),
            label = cms.string('WithBS'),
            maxDistanceToBeam = cms.double(1.0),
            minNdof = cms.double(0.0),
            useBeamConstraint = cms.bool(True)
        )
    )
)


process.hltAK4CaloJetsPFEt5 = cms.EDFilter("EtMinCaloJetSelector",
    etMin = cms.double(5.0),
    filter = cms.bool(False),
    src = cms.InputTag("hltAK4CaloJetsPF")
)


process.hltAK8CaloJetsPFEt5 = cms.EDFilter("EtMinCaloJetSelector",
    etMin = cms.double(5.0),
    filter = cms.bool(False),
    src = cms.InputTag("hltAK8CaloJetsPF")
)


process.hltAK8SingleCaloJet300 = cms.EDFilter("HLT1CaloJet",
    MaxEta = cms.double(2.5),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(1),
    MinPt = cms.double(300.0),
    inputTag = cms.InputTag("hltAK8CaloJetsCorrectedIDPassed"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(85)
)


process.hltBTagCaloCSVp44Single = cms.EDFilter("HLTCaloJetTag",
    JetTags = cms.InputTag("hltCombinedSecondaryVertexBJetTagsCalo"),
    Jets = cms.InputTag("hltSelector8CentralJetsL1FastJet"),
    MaxTag = cms.double(999999.0),
    MinJets = cms.int32(1),
    MinTag = cms.double(0.44),
    TriggerType = cms.int32(86),
    saveTags = cms.bool(True)
)


process.hltBTagCaloDeepCSVp44Single = cms.EDFilter("HLTCaloJetTag",
    JetTags = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsCalo","probb"),
    Jets = cms.InputTag("hltSelector8CentralJetsL1FastJet"),
    MaxTag = cms.double(999999.0),
    MinJets = cms.int32(1),
    MinTag = cms.double(0.21),
    TriggerType = cms.int32(86),
    saveTags = cms.bool(True)
)


process.hltBTagPFBoostedDoubleSecondaryVertexpZSingle = cms.EDFilter("HLTPFJetTag",
    JetTags = cms.InputTag("hltBoostedDBJetTagsPFAK8"),
    Jets = cms.InputTag("hltPFJetForBtagAK8"),
    MaxTag = cms.double(999999.0),
    MinJets = cms.int32(1),
    MinTag = cms.double(0.4),
    TriggerType = cms.int32(86),
    saveTags = cms.bool(True)
)


process.hltBTagPFCSVp075Double = cms.EDFilter("HLTPFJetTag",
    JetTags = cms.InputTag("hltCombinedSecondaryVertexBJetTagsPF"),
    Jets = cms.InputTag("hltPFJetForBtag"),
    MaxTag = cms.double(999999.0),
    MinJets = cms.int32(2),
    MinTag = cms.double(0.75),
    TriggerType = cms.int32(86),
    saveTags = cms.bool(True)
)


process.hltBTagPFCSVp4Single = cms.EDFilter("HLTPFJetTag",
    JetTags = cms.InputTag("hltCombinedSecondaryVertexBJetTagsPFAK8"),
    Jets = cms.InputTag("hltPFJetForBtagAK8"),
    MaxTag = cms.double(999999.0),
    MinJets = cms.int32(1),
    MinTag = cms.double(0.4),
    TriggerType = cms.int32(86),
    saveTags = cms.bool(True)
)


process.hltBTagPFDeepCSVp075Double = cms.EDFilter("HLTPFJetTag",
    JetTags = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsPF","probb"),
    Jets = cms.InputTag("hltPFJetForBtag"),
    MaxTag = cms.double(999999.0),
    MinJets = cms.int32(2),
    MinTag = cms.double(0.5),
    TriggerType = cms.int32(86),
    saveTags = cms.bool(True)
)


process.hltBoolEnd = cms.EDFilter("HLTBool",
    result = cms.bool(True)
)


process.hltBoolFalse = cms.EDFilter("HLTBool",
    result = cms.bool(False)
)


process.hltCaloJetFilterSixC25 = cms.EDFilter("HLT1CaloJet",
    MaxEta = cms.double(2.6),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(6),
    MinPt = cms.double(25.0),
    inputTag = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(86)
)


process.hltCaloSixJet25HT300 = cms.EDFilter("HLTHtMhtFilter",
    htLabels = cms.VInputTag("hltHtMhtCaloJetsSixC25"),
    meffSlope = cms.vdouble(1.0),
    mhtLabels = cms.VInputTag("hltHtMhtCaloJetsSixC25"),
    minHt = cms.vdouble(300.0),
    minMeff = cms.vdouble(0.0),
    minMht = cms.vdouble(0.0),
    saveTags = cms.bool(True)
)


process.hltFastPVJetVertexChecker = cms.EDFilter("JetVertexChecker",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    doFilter = cms.bool(False),
    jetTracks = cms.InputTag("hltFastPVJetTracksAssociator"),
    maxChi2 = cms.double(20.0),
    maxETA = cms.double(2.4),
    maxNJetsToCheck = cms.int32(2),
    maxNjetsOutput = cms.int32(2),
    maxTrackPt = cms.double(20.0),
    minPt = cms.double(0.0),
    minPtRatio = cms.double(0.1),
    newMethod = cms.bool(True),
    pvErr_x = cms.double(0.0015),
    pvErr_y = cms.double(0.0015),
    pvErr_z = cms.double(1.5)
)


process.hltFastPVPixelVertexFilter = cms.EDFilter("VertexSelector",
    cut = cms.string('!isFake && ndof > 0 && abs(z) <= 25 && position.Rho <= 2'),
    filter = cms.bool(True),
    src = cms.InputTag("hltFastPrimaryVertex")
)


process.hltFastPVPixelVerticesFilter = cms.EDFilter("VertexSelector",
    cut = cms.string('!isFake && ndof > 0 && abs(z) <= 25 && position.Rho <= 2'),
    filter = cms.bool(True),
    src = cms.InputTag("hltFastPVPixelVertices")
)


process.hltL1sHTT280to500erIorHTT250to340erQuadJet = cms.EDFilter("HLTL1TSeed",
    L1EGammaInputTag = cms.InputTag("hltGtStage2Digis","EGamma"),
    L1EtSumInputTag = cms.InputTag("hltGtStage2Digis","EtSum"),
    L1GlobalInputTag = cms.InputTag("hltGtStage2Digis"),
    L1JetInputTag = cms.InputTag("hltGtStage2Digis","Jet"),
    L1MuonInputTag = cms.InputTag("hltGtStage2Digis","Muon"),
    L1ObjectMapInputTag = cms.InputTag("hltGtStage2ObjectMap"),
    L1SeedsLogicalExpression = cms.string('L1_HTT280er OR L1_HTT300er OR L1_HTT320er OR L1_HTT340er OR L1_HTT380er OR L1_HTT400er OR L1_HTT450er OR L1_HTT500er OR L1_HTT250er_QuadJet_70_55_40_35_er2p5 OR L1_HTT280er_QuadJet_70_55_40_35_er2p5 OR L1_HTT300er_QuadJet_70_55_40_35_er2p5 OR L1_HTT320er_QuadJet_70_55_40_40_er2p4 OR L1_HTT320er_QuadJet_70_55_40_40_er2p5 OR L1_HTT340er_QuadJet_70_55_40_40_er2p5 OR L1_HTT320er_QuadJet_70_55_45_45_er2p5 OR L1_HTT340er_QuadJet_70_55_45_45_er2p5'),
    L1TauInputTag = cms.InputTag("hltGtStage2Digis","Tau"),
    saveTags = cms.bool(True)
)


process.hltL1sSingleJet170IorSingleJet180IorSingleJet200 = cms.EDFilter("HLTL1TSeed",
    L1EGammaInputTag = cms.InputTag("hltGtStage2Digis","EGamma"),
    L1EtSumInputTag = cms.InputTag("hltGtStage2Digis","EtSum"),
    L1GlobalInputTag = cms.InputTag("hltGtStage2Digis"),
    L1JetInputTag = cms.InputTag("hltGtStage2Digis","Jet"),
    L1MuonInputTag = cms.InputTag("hltGtStage2Digis","Muon"),
    L1ObjectMapInputTag = cms.InputTag("hltGtStage2ObjectMap"),
    L1SeedsLogicalExpression = cms.string('L1_SingleJet170 OR L1_SingleJet180 OR L1_SingleJet200'),
    L1TauInputTag = cms.InputTag("hltGtStage2Digis","Tau"),
    saveTags = cms.bool(True)
)


process.hltPFJetFilterSixC32 = cms.EDFilter("HLT1PFJet",
    MaxEta = cms.double(2.6),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(6),
    MinPt = cms.double(32.0),
    inputTag = cms.InputTag("hltAK4PFJetsCorrected"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(86)
)


process.hltPFJetForBtagSelector = cms.EDFilter("HLT1PFJet",
    MaxEta = cms.double(2.6),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(1),
    MinPt = cms.double(30.0),
    inputTag = cms.InputTag("hltAK4PFJetsCorrected"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(86)
)


process.hltPFJetForBtagSelectorAK8 = cms.EDFilter("HLT1PFJet",
    MaxEta = cms.double(2.6),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(1),
    MinPt = cms.double(30.0),
    inputTag = cms.InputTag("hltAK8PFJetsCorrected"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(86)
)


process.hltPFSixJet32HT380 = cms.EDFilter("HLTHtMhtFilter",
    htLabels = cms.VInputTag("hltHtMhtPFJetsSixC32"),
    meffSlope = cms.vdouble(1.0),
    mhtLabels = cms.VInputTag("hltHtMhtPFJetsSixC32"),
    minHt = cms.vdouble(380.0),
    minMeff = cms.vdouble(0.0),
    minMht = cms.vdouble(0.0),
    saveTags = cms.bool(True)
)


process.hltPreAK8PFJet400 = cms.EDFilter("HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag("hltGtStage2Digis"),
    offset = cms.uint32(0)
)


process.hltPreAK8PFJetXPFAK8DoubleBpY = cms.EDFilter("HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag("hltGtStage2Digis"),
    offset = cms.uint32(0)
)


process.hltPrePFHT380SixPFJet32DoublePFBTagCSV2p2 = cms.EDFilter("HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag("hltGtStage2Digis"),
    offset = cms.uint32(0)
)


process.hltPrePFHT380SixPFJet32DoublePFBTagDeepCSV2p2 = cms.EDFilter("HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag("hltGtStage2Digis"),
    offset = cms.uint32(0)
)


process.hltSelector4CentralJetsL1FastJet = cms.EDFilter("LargestEtCaloJetSelector",
    filter = cms.bool(False),
    maxNumber = cms.uint32(4),
    src = cms.InputTag("hltSelectorCentralJets20L1FastJeta")
)


process.hltSelector8CentralJetsL1FastJet = cms.EDFilter("LargestEtCaloJetSelector",
    filter = cms.bool(False),
    maxNumber = cms.uint32(8),
    src = cms.InputTag("hltSelectorCentralJets30L1FastJeta")
)


process.hltSelectorCentralJets20L1FastJeta = cms.EDFilter("EtaRangeCaloJetSelector",
    etaMax = cms.double(2.4),
    etaMin = cms.double(-2.4),
    src = cms.InputTag("hltSelectorJets20L1FastJet")
)


process.hltSelectorCentralJets30L1FastJeta = cms.EDFilter("EtaRangeCaloJetSelector",
    etaMax = cms.double(2.4),
    etaMin = cms.double(-2.4),
    src = cms.InputTag("hltSelectorJets30L1FastJet")
)


process.hltSelectorJets20L1FastJet = cms.EDFilter("EtMinCaloJetSelector",
    etMin = cms.double(20.0),
    filter = cms.bool(False),
    src = cms.InputTag("hltAK4CaloJetsCorrected")
)


process.hltSelectorJets30L1FastJet = cms.EDFilter("EtMinCaloJetSelector",
    etMin = cms.double(30.0),
    filter = cms.bool(False),
    src = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed")
)


process.hltSingleCaloJet350AK8 = cms.EDFilter("HLT1CaloJet",
    MaxEta = cms.double(5.0),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(1),
    MinPt = cms.double(350.0),
    inputTag = cms.InputTag("hltAK8CaloJetsCorrectedIDPassed"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(85)
)


process.hltSinglePFJet330AK8 = cms.EDFilter("HLT1PFJet",
    MaxEta = cms.double(5.0),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(1),
    MinPt = cms.double(330.0),
    inputTag = cms.InputTag("hltAK8PFJetsCorrectedMatchedToCaloJets300"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(85)
)


process.hltSinglePFJet400AK8 = cms.EDFilter("HLT1PFJet",
    MaxEta = cms.double(5.0),
    MaxMass = cms.double(-1.0),
    MinE = cms.double(-1.0),
    MinEta = cms.double(-1.0),
    MinMass = cms.double(-1.0),
    MinN = cms.int32(1),
    MinPt = cms.double(400.0),
    inputTag = cms.InputTag("hltPFJetsCorrectedMatchedToCaloJets350AK8"),
    saveTags = cms.bool(True),
    triggerType = cms.int32(85)
)


process.hltTriggerType = cms.EDFilter("HLTTriggerTypeFilter",
    SelectedTriggerType = cms.int32(1)
)


process.hltVerticesPFFilter = cms.EDFilter("VertexSelector",
    cut = cms.string('!isFake'),
    filter = cms.bool(True),
    src = cms.InputTag("hltVerticesPFSelector")
)


process.hltVerticesPFSelector = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = cms.PSet(
        maxRho = cms.double(2.0),
        maxZ = cms.double(24.0),
        minNdof = cms.double(4.0),
        pvSrc = cms.InputTag("hltVerticesPF")
    ),
    src = cms.InputTag("hltVerticesPF")
)


process.hltGetConditions = cms.EDAnalyzer("EventSetupRecordDataGetter",
    toGet = cms.VPSet(),
    verbose = cms.untracked.bool(False)
)


process.hltGetRaw = cms.EDAnalyzer("HLTGetRaw",
    RawDataCollection = cms.InputTag("rawDataCollector")
)


process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
    fileName = cms.untracked.string('DQMIO.root')
)


process.hltOutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('HLT_*')
    ),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('HLT_tunedcontent_v2-1.root'),
    outputCommands = cms.untracked.vstring(
        'drop *', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoJetedmRefToBaseProdTofloatsAssociationVector_*_*_*', 
        'keep *RefToBaseProd*_*_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep *_ak4GenJets*_*_*', 
        'keep *_ak4GenJetsNoNu*_*_*', 
        'keep *_ak4PFJets*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed*_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet*_*_*', 
        'keep *_hltAK4PFJetsCorrected*_*_*', 
        'keep *_hltPFJetForBtag*_*_*', 
        'keep *_ak8GenJets*_*_*', 
        'keep *_ak8PFJets*_*_*', 
        'keep *_ak8GenJetsNoNu*_*_*', 
        'keep *_hltAK8CaloJetsCorrectedIDPassed*_*_*', 
        'keep *_hltAK8PFJetsCorrected*_*_*', 
        'keep *_hltPFJetForBtagAK8*_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_addPileupInfo_*_*'
    )
)


process.DQMStore = cms.Service("DQMStore",
    LSbasedMode = cms.untracked.bool(False),
    collateHistograms = cms.untracked.bool(False),
    enableMultiThread = cms.untracked.bool(True),
    forceResetOnBeginLumi = cms.untracked.bool(False),
    referenceFileName = cms.untracked.string(''),
    verbose = cms.untracked.int32(0),
    verboseQT = cms.untracked.int32(0)
)


process.FastTimerService = cms.Service("FastTimerService",
    dqmLumiSectionsRange = cms.untracked.uint32(2500),
    dqmMemoryRange = cms.untracked.double(1000000.0),
    dqmMemoryResolution = cms.untracked.double(5000.0),
    dqmModuleMemoryRange = cms.untracked.double(100000.0),
    dqmModuleMemoryResolution = cms.untracked.double(500.0),
    dqmModuleTimeRange = cms.untracked.double(40.0),
    dqmModuleTimeResolution = cms.untracked.double(0.2),
    dqmPath = cms.untracked.string('HLT/TimerService'),
    dqmPathMemoryRange = cms.untracked.double(1000000.0),
    dqmPathMemoryResolution = cms.untracked.double(5000.0),
    dqmPathTimeRange = cms.untracked.double(100.0),
    dqmPathTimeResolution = cms.untracked.double(0.5),
    dqmTimeRange = cms.untracked.double(2000.0),
    dqmTimeResolution = cms.untracked.double(5.0),
    enableDQM = cms.untracked.bool(True),
    enableDQMbyLumiSection = cms.untracked.bool(True),
    enableDQMbyModule = cms.untracked.bool(False),
    enableDQMbyPath = cms.untracked.bool(False),
    enableDQMbyProcesses = cms.untracked.bool(True),
    printEventSummary = cms.untracked.bool(False),
    printJobSummary = cms.untracked.bool(True),
    printRunSummary = cms.untracked.bool(True)
)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    ),
    categories = cms.untracked.vstring(
        'FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary', 
        'TriggerSummaryProducerAOD', 
        'L1GtTrigReport', 
        'L1TGlobalSummary', 
        'HLTrigReport', 
        'FastReport'
    ),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            reportEvery = cms.untracked.int32(1)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        suppressDebug = cms.untracked.vstring(),
        suppressError = cms.untracked.vstring(),
        suppressInfo = cms.untracked.vstring(),
        suppressWarning = cms.untracked.vstring(),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(10)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(10)
        ),
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True),
        suppressDebug = cms.untracked.vstring(),
        suppressError = cms.untracked.vstring(),
        suppressInfo = cms.untracked.vstring(),
        suppressWarning = cms.untracked.vstring(),
        threshold = cms.untracked.string('INFO')
    ),
    destinations = cms.untracked.vstring(
        'warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'
    ),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True),
        suppressDebug = cms.untracked.vstring(),
        suppressError = cms.untracked.vstring(),
        suppressInfo = cms.untracked.vstring(),
        suppressWarning = cms.untracked.vstring(),
        threshold = cms.untracked.string('INFO')
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        placeholder = cms.untracked.bool(True),
        suppressDebug = cms.untracked.vstring(),
        suppressError = cms.untracked.vstring(),
        suppressInfo = cms.untracked.vstring(),
        suppressWarning = cms.untracked.vstring(),
        threshold = cms.untracked.string('INFO')
    ),
    statistics = cms.untracked.vstring('cerr'),
    suppressDebug = cms.untracked.vstring(),
    suppressError = cms.untracked.vstring(
        'hltOnlineBeamSpot', 
        'hltL3MuonCandidates', 
        'hltL3TkTracksFromL2OIState', 
        'hltPFJetCtfWithMaterialTracks', 
        'hltL3TkTracksFromL2IOHit', 
        'hltL3TkTracksFromL2OIHit'
    ),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(
        'hltOnlineBeamSpot', 
        'hltCtf3HitL1SeededWithMaterialTracks', 
        'hltL3MuonsOIState', 
        'hltPixelTracksForHighMult', 
        'hltHITPixelTracksHE', 
        'hltHITPixelTracksHB', 
        'hltCtfL1SeededWithMaterialTracks', 
        'hltRegionalTracksForL3MuonIsolation', 
        'hltSiPixelClusters', 
        'hltActivityStartUpElectronPixelSeeds', 
        'hltLightPFTracks', 
        'hltPixelVertices3DbbPhi', 
        'hltL3MuonsIOHit', 
        'hltPixelTracks', 
        'hltSiPixelDigis', 
        'hltL3MuonsOIHit', 
        'hltL1SeededElectronGsfTracks', 
        'hltL1SeededStartUpElectronPixelSeeds', 
        'hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV', 
        'hltCtfActivityWithMaterialTracks'
    ),
    threshold = cms.untracked.string('INFO'),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True),
        suppressDebug = cms.untracked.vstring(),
        suppressError = cms.untracked.vstring(),
        suppressInfo = cms.untracked.vstring(),
        suppressWarning = cms.untracked.vstring(),
        threshold = cms.untracked.string('INFO')
    )
)


process.ThroughputService = cms.Service("ThroughputService",
    dqmPath = cms.untracked.string('HLT/Throughput'),
    timeRange = cms.untracked.double(60000.0),
    timeResolution = cms.untracked.double(5.828)
)


process.AnyDirectionAnalyticalPropagator = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('AnyDirectionAnalyticalPropagator'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('anyDirection')
)


process.CSCChannelMapperESProducer = cms.ESProducer("CSCChannelMapperESProducer",
    AlgoName = cms.string('CSCChannelMapperPostls1')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(False),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CSCIndexerESProducer = cms.ESProducer("CSCIndexerESProducer",
    AlgoName = cms.string('CSCIndexerPostls1')
)


process.CSCObjectMapESProducer = cms.ESProducer("CSCObjectMapESProducer",
    appendToDataLabel = cms.string('')
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL', 
        'ZDC', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapAuto = cms.untracked.bool(False),
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz'),
    SkipHE = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.ClusterShapeHitFilterESProducer = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('ClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    )
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.DTObjectMapESProducer = cms.ESProducer("DTObjectMapESProducer",
    appendToDataLabel = cms.string('')
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.GlobalParameters = cms.ESProducer("StableParametersTrivialProducer",
    IfCaloEtaNumberBits = cms.uint32(4),
    IfMuEtaNumberBits = cms.uint32(6),
    NumberChips = cms.uint32(1),
    NumberConditionChips = cms.uint32(1),
    NumberL1CenJet = cms.uint32(4),
    NumberL1EGamma = cms.uint32(12),
    NumberL1ForJet = cms.uint32(4),
    NumberL1IsoEG = cms.uint32(4),
    NumberL1Jet = cms.uint32(12),
    NumberL1JetCounts = cms.uint32(12),
    NumberL1Mu = cms.uint32(4),
    NumberL1Muon = cms.uint32(8),
    NumberL1NoIsoEG = cms.uint32(4),
    NumberL1Tau = cms.uint32(12),
    NumberL1TauJet = cms.uint32(4),
    NumberPhysTriggers = cms.uint32(512),
    NumberPhysTriggersExtended = cms.uint32(64),
    NumberPsbBoards = cms.int32(7),
    NumberTechnicalTriggers = cms.uint32(64),
    OrderConditionChip = cms.vint32(1),
    OrderOfChip = cms.vint32(1),
    PinsOnChip = cms.uint32(512),
    PinsOnConditionChip = cms.uint32(512),
    TotalBxInEvent = cms.int32(5),
    UnitLength = cms.int32(8),
    WordLength = cms.int32(64),
    appendToDataLabel = cms.string('')
)


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.HcalTopologyIdealEP = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(True),
    appendToDataLabel = cms.string('')
)


process.MaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterial'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.MaterialPropagatorForHI = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialForHI'),
    Mass = cms.double(0.139),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.MaterialPropagatorParabolicMF = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialParabolicMf'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.OppositeMaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.OppositeMaterialPropagatorForHI = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialOppositeForHI'),
    Mass = cms.double(0.139),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.OppositeMaterialPropagatorParabolicMF = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.OppositePropagatorWithMaterialForMixedStep = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialForMixedStepOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(0.1),
    useRungeKutta = cms.bool(False)
)


process.ParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(-1),
    version = cms.string('Parabolic')
)


process.PropagatorWithMaterialForLoopers = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialForLoopers'),
    Mass = cms.double(0.1396),
    MaxDPhi = cms.double(4.0),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.PropagatorWithMaterialForMixedStep = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialForMixedStep'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(0.1),
    useRungeKutta = cms.bool(False)
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.SiStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.SiStripRegionConnectivity = cms.ESProducer("SiStripRegionConnectivity",
    EtaDivisions = cms.untracked.uint32(20),
    EtaMax = cms.untracked.double(2.5),
    PhiDivisions = cms.untracked.uint32(20)
)


process.SimpleSecondaryVertex3TrkComputer = cms.ESProducer("SimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(3),
    minVertices = cms.uint32(1),
    unBoost = cms.bool(False),
    use3d = cms.bool(True),
    useSignificance = cms.bool(True)
)


process.SteppingHelixPropagatorAny = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorAny'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('anyDirection'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.TrackerDigiGeometryESModule = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.TrackerGeometricDetESModule = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(-1)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.caloDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('CaloDetIdAssociator'),
    etaBinSize = cms.double(0.087),
    hcalRegion = cms.int32(2),
    includeBadChambers = cms.bool(False),
    includeGEM = cms.bool(False),
    includeME0 = cms.bool(False),
    nEta = cms.int32(70),
    nPhi = cms.int32(72)
)


process.cosmicsNavigationSchoolESProducer = cms.ESProducer("NavigationSchoolESProducer",
    ComponentName = cms.string('CosmicNavigationSchool'),
    SimpleMagneticField = cms.string('')
)


process.ecalDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('EcalDetIdAssociator'),
    etaBinSize = cms.double(0.02),
    hcalRegion = cms.int32(2),
    includeBadChambers = cms.bool(False),
    includeGEM = cms.bool(False),
    includeME0 = cms.bool(False),
    nEta = cms.int32(300),
    nPhi = cms.int32(360)
)


process.ecalSeverityLevel = cms.ESProducer("EcalSeverityLevelESProducer",
    dbstatusMask = cms.PSet(
        kBad = cms.vstring(
            'kNonRespondingIsolated', 
            'kDeadVFE', 
            'kDeadFE', 
            'kNoDataNoTP'
        ),
        kGood = cms.vstring('kOk'),
        kProblematic = cms.vstring(
            'kDAC', 
            'kNoLaser', 
            'kNoisy', 
            'kNNoisy', 
            'kNNNoisy', 
            'kNNNNoisy', 
            'kNNNNNoisy', 
            'kFixedG6', 
            'kFixedG1', 
            'kFixedG0'
        ),
        kRecovered = cms.vstring(),
        kTime = cms.vstring(),
        kWeird = cms.vstring()
    ),
    flagMask = cms.PSet(
        kBad = cms.vstring(
            'kFaultyHardware', 
            'kDead', 
            'kKilled'
        ),
        kGood = cms.vstring('kGood'),
        kProblematic = cms.vstring(
            'kPoorReco', 
            'kPoorCalib', 
            'kNoisy', 
            'kSaturated'
        ),
        kRecovered = cms.vstring(
            'kLeadingEdgeRecovered', 
            'kTowerRecovered'
        ),
        kTime = cms.vstring('kOutOfTime'),
        kWeird = cms.vstring(
            'kWeird', 
            'kDiWeird'
        )
    ),
    timeThresh = cms.double(2.0)
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('HcalDetIdAssociator'),
    etaBinSize = cms.double(0.087),
    hcalRegion = cms.int32(2),
    includeBadChambers = cms.bool(False),
    includeGEM = cms.bool(False),
    includeME0 = cms.bool(False),
    nEta = cms.int32(70),
    nPhi = cms.int32(72)
)


process.hcalRecAlgos = cms.ESProducer("HcalRecAlgoESProducer",
    DropChannelStatusBits = cms.vstring(
        'HcalCellMask', 
        'HcalCellOff', 
        'HcalCellDead'
    ),
    RecoveredRecHitBits = cms.vstring(),
    SeverityLevels = cms.VPSet(
        cms.PSet(
            ChannelStatus = cms.vstring(),
            Level = cms.int32(0),
            RecHitFlags = cms.vstring('TimingFromTDC')
        ), 
        cms.PSet(
            ChannelStatus = cms.vstring('HcalCellCaloTowerProb'),
            Level = cms.int32(1),
            RecHitFlags = cms.vstring()
        ), 
        cms.PSet(
            ChannelStatus = cms.vstring('HcalCellExcludeFromHBHENoiseSummary'),
            Level = cms.int32(5),
            RecHitFlags = cms.vstring()
        ), 
        cms.PSet(
            ChannelStatus = cms.vstring(),
            Level = cms.int32(8),
            RecHitFlags = cms.vstring(
                'HBHEHpdHitMultiplicity', 
                'HBHEIsolatedNoise', 
                'HBHEFlatNoise', 
                'HBHESpikeNoise', 
                'HBHETS4TS5Noise', 
                'HBHENegativeNoise', 
                'HBHEPulseFitBit', 
                'HBHEOOTPU'
            )
        ), 
        cms.PSet(
            ChannelStatus = cms.vstring(),
            Level = cms.int32(11),
            RecHitFlags = cms.vstring(
                'HFLongShort', 
                'HFS8S1Ratio', 
                'HFPET', 
                'HFSignalAsymmetry'
            )
        ), 
        cms.PSet(
            ChannelStatus = cms.vstring('HcalCellHot'),
            Level = cms.int32(15),
            RecHitFlags = cms.vstring()
        ), 
        cms.PSet(
            ChannelStatus = cms.vstring(
                'HcalCellOff', 
                'HcalCellDead'
            ),
            Level = cms.int32(20),
            RecHitFlags = cms.vstring()
        )
    ),
    appendToDataLabel = cms.string(''),
    phase = cms.uint32(1)
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer")


process.hltBoostedDoubleSecondaryVertexAK8Computer = cms.ESProducer("CandidateBoostedDoubleSecondaryVertexESProducer",
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SecondaryVertex/data/BoostedDoubleSV_AK8_BDT_v4.weights.xml.gz')
)


process.hltCombinedSecondaryVertex = cms.ESProducer("CombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring(
        'CombinedSVRecoVertex', 
        'CombinedSVPseudoVertex', 
        'CombinedSVNoVertex'
    ),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string('HLT'),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(3),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5.0),
        maxDistToAxis = cms.double(0.07),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0)
    ),
    trackSelection = cms.PSet(
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5.0),
        maxDistToAxis = cms.double(0.07),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.hltCombinedSecondaryVertexV2 = cms.ESProducer("CombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring(
        'CombinedSVIVFV2RecoVertex', 
        'CombinedSVIVFV2PseudoVertex', 
        'CombinedSVIVFV2NoVertex'
    ),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string('HLT'),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(3),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5.0),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500.0),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3.0),
        min_pT = cms.double(120.0),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5.0),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500.0),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3.0),
        min_pT = cms.double(120.0),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.hltDisplacedDijethltESPPromptTrackCountingESProducer = cms.ESProducer("PromptTrackCountingESProducer",
    deltaR = cms.double(-1.0),
    deltaRmin = cms.double(0.0),
    impactParameterType = cms.int32(1),
    maxImpactParameter = cms.double(0.1),
    maxImpactParameterSig = cms.double(999999.0),
    maximumDecayLength = cms.double(999999.0),
    maximumDistanceToJetAxis = cms.double(999999.0),
    minimumImpactParameter = cms.double(-1.0),
    nthTrack = cms.int32(-1),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True)
)


process.hltDisplacedDijethltESPTrackCounting2D1st = cms.ESProducer("TrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(1),
    max_pT = cms.double(500.0),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3.0),
    maximumDecayLength = cms.double(999999.0),
    maximumDistanceToJetAxis = cms.double(9999999.0),
    min_pT = cms.double(120.0),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(0.05),
    nthTrack = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(False),
    useVariableJTA = cms.bool(False)
)


process.hltESPAnalyticalPropagator = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('hltESPAnalyticalPropagator'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum')
)


process.hltESPBwdAnalyticalPropagator = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('hltESPBwdAnalyticalPropagator'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum')
)


process.hltESPBwdElectronPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('hltESPBwdElectronPropagator'),
    Mass = cms.double(0.000511),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.hltESPChi2ChargeLooseMeasurementEstimator16 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2ChargeLooseMeasurementEstimator16'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPChi2ChargeMeasurementEstimator16 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2ChargeMeasurementEstimator16'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPChi2ChargeMeasurementEstimator2000 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2ChargeMeasurementEstimator2000'),
    MaxChi2 = cms.double(2000.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPChi2ChargeMeasurementEstimator30 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
    MaxChi2 = cms.double(30.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPChi2ChargeMeasurementEstimator9 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2ChargeMeasurementEstimator9'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(15.0)
)


process.hltESPChi2ChargeMeasurementEstimator9ForHI = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2ChargeMeasurementEstimator9ForHI'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutForHI')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(15.0)
)


process.hltESPChi2ChargeTightMeasurementEstimator16 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2ChargeTightMeasurementEstimator16'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPChi2MeasurementEstimator100 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2MeasurementEstimator100'),
    MaxChi2 = cms.double(40.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1e+12),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(4.0)
)


process.hltESPChi2MeasurementEstimator16 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2MeasurementEstimator16'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3.0)
)


process.hltESPChi2MeasurementEstimator30 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2MeasurementEstimator30'),
    MaxChi2 = cms.double(30.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3.0)
)


process.hltESPChi2MeasurementEstimator9 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPChi2MeasurementEstimator9'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3.0)
)


process.hltESPCloseComponentsMerger5D = cms.ESProducer("CloseComponentsMergerESProducer5D",
    ComponentName = cms.string('hltESPCloseComponentsMerger5D'),
    DistanceMeasure = cms.string('hltESPKullbackLeiblerDistance5D'),
    MaxComponents = cms.int32(12)
)


process.hltESPDetachedQuadStepChi2ChargeMeasurementEstimator9 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPDetachedQuadStepChi2ChargeMeasurementEstimator9'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPDetachedQuadStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPDetachedQuadStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.13)
)


process.hltESPDetachedStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPDetachedStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.13)
)


process.hltESPDetachedTripletStepChi2ChargeMeasurementEstimator9 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPDetachedTripletStepChi2ChargeMeasurementEstimator9'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPDetachedTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPDetachedTripletStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.13)
)


process.hltESPDisplacedDijethltPromptTrackCountingESProducer = cms.ESProducer("PromptTrackCountingESProducer",
    deltaR = cms.double(-1.0),
    deltaRmin = cms.double(0.0),
    impactParameterType = cms.int32(1),
    maxImpactParameter = cms.double(0.1),
    maxImpactParameterSig = cms.double(999999.0),
    maximumDecayLength = cms.double(999999.0),
    maximumDistanceToJetAxis = cms.double(999999.0),
    minimumImpactParameter = cms.double(-1.0),
    nthTrack = cms.int32(-1),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True)
)


process.hltESPDisplacedDijethltPromptTrackCountingESProducerLong = cms.ESProducer("PromptTrackCountingESProducer",
    deltaR = cms.double(-1.0),
    deltaRmin = cms.double(0.0),
    impactParameterType = cms.int32(1),
    maxImpactParameter = cms.double(0.2),
    maxImpactParameterSig = cms.double(999999.0),
    maximumDecayLength = cms.double(999999.0),
    maximumDistanceToJetAxis = cms.double(999999.0),
    minimumImpactParameter = cms.double(-1.0),
    nthTrack = cms.int32(-1),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True)
)


process.hltESPDisplacedDijethltTrackCounting2D1st = cms.ESProducer("TrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(1),
    max_pT = cms.double(500.0),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3.0),
    maximumDecayLength = cms.double(999999.0),
    maximumDistanceToJetAxis = cms.double(9999999.0),
    min_pT = cms.double(120.0),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(0.05),
    nthTrack = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(False),
    useVariableJTA = cms.bool(False)
)


process.hltESPDisplacedDijethltTrackCounting2D2ndLong = cms.ESProducer("TrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(1),
    max_pT = cms.double(500.0),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3.0),
    maximumDecayLength = cms.double(999999.0),
    maximumDistanceToJetAxis = cms.double(9999999.0),
    min_pT = cms.double(120.0),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(0.2),
    nthTrack = cms.int32(2),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.hltESPDummyDetLayerGeometry = cms.ESProducer("DetLayerGeometryESProducer",
    ComponentName = cms.string('hltESPDummyDetLayerGeometry')
)


process.hltESPEcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.hltESPElectronMaterialEffects = cms.ESProducer("GsfMaterialEffectsESProducer",
    BetheHeitlerCorrection = cms.int32(2),
    BetheHeitlerParametrization = cms.string('BetheHeitler_cdfmom_nC6_O5.par'),
    ComponentName = cms.string('hltESPElectronMaterialEffects'),
    EnergyLossUpdator = cms.string('GsfBetheHeitlerUpdator'),
    Mass = cms.double(0.000511),
    MultipleScatteringUpdator = cms.string('MultipleScatteringUpdator')
)


process.hltESPFastSteppingHelixPropagatorAny = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('anyDirection'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(True)
)


process.hltESPFastSteppingHelixPropagatorOpposite = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('hltESPFastSteppingHelixPropagatorOpposite'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(True)
)


process.hltESPFittingSmootherIT = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('hltESPFittingSmootherIT'),
    EstimateCut = cms.double(-1.0),
    Fitter = cms.string('hltESPTrajectoryFitterRK'),
    LogPixelProbabilityCut = cms.double(-16.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPTrajectorySmootherRK'),
    appendToDataLabel = cms.string('')
)


process.hltESPFittingSmootherRK = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(False),
    ComponentName = cms.string('hltESPFittingSmootherRK'),
    EstimateCut = cms.double(-1.0),
    Fitter = cms.string('hltESPTrajectoryFitterRK'),
    LogPixelProbabilityCut = cms.double(-16.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(False),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPTrajectorySmootherRK'),
    appendToDataLabel = cms.string('')
)


process.hltESPFlexibleKFFittingSmoother = cms.ESProducer("FlexibleKFFittingSmootherESProducer",
    ComponentName = cms.string('hltESPFlexibleKFFittingSmoother'),
    appendToDataLabel = cms.string(''),
    looperFitter = cms.string('hltESPKFFittingSmootherForLoopers'),
    standardFitter = cms.string('hltESPKFFittingSmootherWithOutliersRejectionAndRK')
)


process.hltESPFwdElectronPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('hltESPFwdElectronPropagator'),
    Mass = cms.double(0.000511),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.hltESPGlobalDetLayerGeometry = cms.ESProducer("GlobalDetLayerGeometryESProducer",
    ComponentName = cms.string('hltESPGlobalDetLayerGeometry')
)


process.hltESPGlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.hltESPGsfElectronFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('hltESPGsfElectronFittingSmoother'),
    EstimateCut = cms.double(-1.0),
    Fitter = cms.string('hltESPGsfTrajectoryFitter'),
    LogPixelProbabilityCut = cms.double(-16.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPGsfTrajectorySmoother'),
    appendToDataLabel = cms.string('')
)


process.hltESPGsfTrajectoryFitter = cms.ESProducer("GsfTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPGsfTrajectoryFitter'),
    GeometricalPropagator = cms.string('hltESPAnalyticalPropagator'),
    MaterialEffectsUpdator = cms.string('hltESPElectronMaterialEffects'),
    Merger = cms.string('hltESPCloseComponentsMerger5D'),
    RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry')
)


process.hltESPGsfTrajectorySmoother = cms.ESProducer("GsfTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPGsfTrajectorySmoother'),
    ErrorRescaling = cms.double(100.0),
    GeometricalPropagator = cms.string('hltESPBwdAnalyticalPropagator'),
    MaterialEffectsUpdator = cms.string('hltESPElectronMaterialEffects'),
    Merger = cms.string('hltESPCloseComponentsMerger5D'),
    RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry')
)


process.hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30'),
    MaxChi2 = cms.double(30.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(15.0)
)


process.hltESPInitialStepChi2ChargeMeasurementEstimator30 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPInitialStepChi2ChargeMeasurementEstimator30'),
    MaxChi2 = cms.double(30.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(15.0)
)


process.hltESPInitialStepChi2MeasurementEstimator36 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPInitialStepChi2MeasurementEstimator36'),
    MaxChi2 = cms.double(36.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3.0)
)


process.hltESPKFFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(False),
    ComponentName = cms.string('hltESPKFFittingSmoother'),
    EstimateCut = cms.double(-1.0),
    Fitter = cms.string('hltESPKFTrajectoryFitter'),
    LogPixelProbabilityCut = cms.double(-16.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(False),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPKFTrajectorySmoother'),
    appendToDataLabel = cms.string('')
)


process.hltESPKFFittingSmootherForL2Muon = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(False),
    ComponentName = cms.string('hltESPKFFittingSmootherForL2Muon'),
    EstimateCut = cms.double(-1.0),
    Fitter = cms.string('hltESPKFTrajectoryFitterForL2Muon'),
    LogPixelProbabilityCut = cms.double(-16.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(False),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPKFTrajectorySmootherForL2Muon'),
    appendToDataLabel = cms.string('')
)


process.hltESPKFFittingSmootherForLoopers = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('hltESPKFFittingSmootherForLoopers'),
    EstimateCut = cms.double(20.0),
    Fitter = cms.string('hltESPKFTrajectoryFitterForLoopers'),
    LogPixelProbabilityCut = cms.double(-14.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPKFTrajectorySmootherForLoopers'),
    appendToDataLabel = cms.string('')
)


process.hltESPKFFittingSmootherWithOutliersRejectionAndRK = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('hltESPKFFittingSmootherWithOutliersRejectionAndRK'),
    EstimateCut = cms.double(20.0),
    Fitter = cms.string('hltESPRKTrajectoryFitter'),
    LogPixelProbabilityCut = cms.double(-14.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPRKTrajectorySmoother'),
    appendToDataLabel = cms.string('')
)


process.hltESPKFTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPKFTrajectoryFitter'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.hltESPKFTrajectoryFitterForL2Muon = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPKFTrajectoryFitterForL2Muon'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.hltESPKFTrajectoryFitterForLoopers = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPKFTrajectoryFitterForLoopers'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.hltESPKFTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPKFTrajectorySmoother'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100.0),
    minHits = cms.int32(3)
)


process.hltESPKFTrajectorySmootherForL2Muon = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPKFTrajectorySmootherForL2Muon'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPFastSteppingHelixPropagatorOpposite'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100.0),
    minHits = cms.int32(3)
)


process.hltESPKFTrajectorySmootherForLoopers = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPKFTrajectorySmootherForLoopers'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(3)
)


process.hltESPKFTrajectorySmootherForMuonTrackLoader = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPSmartPropagatorAnyOpposite'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(3)
)


process.hltESPKFUpdator = cms.ESProducer("KFUpdatorESProducer",
    ComponentName = cms.string('hltESPKFUpdator')
)


process.hltESPKullbackLeiblerDistance5D = cms.ESProducer("DistanceBetweenComponentsESProducer5D",
    ComponentName = cms.string('hltESPKullbackLeiblerDistance5D'),
    DistanceMeasure = cms.string('KullbackLeibler')
)


process.hltESPL3MuKFTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPL3MuKFTrajectoryFitter'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPSmartPropagatorAny'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPLowPtQuadStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPLowPtQuadStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.16)
)


process.hltESPLowPtStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPLowPtStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.16)
)


process.hltESPLowPtTripletStepChi2ChargeMeasurementEstimator9 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPLowPtTripletStepChi2ChargeMeasurementEstimator9'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPLowPtTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPLowPtTripletStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.16)
)


process.hltESPMeasurementTracker = cms.ESProducer("MeasurementTrackerESProducer",
    ComponentName = cms.string('hltESPMeasurementTracker'),
    DebugPixelModuleQualityDB = cms.untracked.bool(False),
    DebugPixelROCQualityDB = cms.untracked.bool(False),
    DebugStripAPVFiberQualityDB = cms.untracked.bool(False),
    DebugStripModuleQualityDB = cms.untracked.bool(False),
    DebugStripStripQualityDB = cms.untracked.bool(False),
    HitMatcher = cms.string('StandardMatcher'),
    MaskBadAPVFibers = cms.bool(True),
    PixelCPE = cms.string('hltESPPixelCPEGeneric'),
    SiStripQualityLabel = cms.string(''),
    StripCPE = cms.string('hltESPStripCPEfromTrackAngle'),
    UsePixelModuleQualityDB = cms.bool(True),
    UsePixelROCQualityDB = cms.bool(True),
    UseStripAPVFiberQualityDB = cms.bool(True),
    UseStripModuleQualityDB = cms.bool(True),
    UseStripStripQualityDB = cms.bool(True),
    badStripCuts = cms.PSet(
        TEC = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TIB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TID = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TOB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        )
    )
)


process.hltESPMixedStepClusterShapeHitFilter = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('hltESPMixedStepClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    )
)


process.hltESPMixedStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPMixedStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.11)
)


process.hltESPMixedTripletStepChi2ChargeMeasurementEstimator16 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPMixedTripletStepChi2ChargeMeasurementEstimator16'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPMixedTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPMixedTripletStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.11)
)


process.hltESPMuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.hltESPMuonTransientTrackingRecHitBuilder = cms.ESProducer("MuonTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('hltESPMuonTransientTrackingRecHitBuilder')
)


process.hltESPPixelCPEGeneric = cms.ESProducer("PixelCPEGenericESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('hltESPPixelCPEGeneric'),
    DoCosmics = cms.bool(False),
    EdgeClusterErrorX = cms.double(50.0),
    EdgeClusterErrorY = cms.double(85.0),
    IrradiationBiasCorrection = cms.bool(False),
    LoadTemplatesFromDB = cms.bool(True),
    MagneticFieldRecord = cms.ESInputTag(""),
    PixelErrorParametrization = cms.string('NOTcmsim'),
    TruncatePixelCharge = cms.bool(True),
    UseErrorsFromTemplates = cms.bool(True),
    eff_charge_cut_highX = cms.double(1.0),
    eff_charge_cut_highY = cms.double(1.0),
    eff_charge_cut_lowX = cms.double(0.0),
    eff_charge_cut_lowY = cms.double(0.0),
    inflate_all_errors_no_trk_angle = cms.bool(False),
    inflate_errors = cms.bool(False),
    size_cutX = cms.double(3.0),
    size_cutY = cms.double(3.0),
    useLAAlignmentOffsets = cms.bool(False),
    useLAWidthFromDB = cms.bool(False)
)


process.hltESPPixelCPETemplateReco = cms.ESProducer("PixelCPETemplateRecoESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('hltESPPixelCPETemplateReco'),
    DoCosmics = cms.bool(False),
    DoLorentz = cms.bool(True),
    LoadTemplatesFromDB = cms.bool(True),
    UseClusterSplitter = cms.bool(False),
    speed = cms.int32(-2)
)


process.hltESPPixelLessStepChi2ChargeMeasurementEstimator16 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPPixelLessStepChi2ChargeMeasurementEstimator16'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPPixelLessStepClusterShapeHitFilter = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('hltESPPixelLessStepClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    )
)


process.hltESPPixelLessStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPPixelLessStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.11)
)


process.hltESPPixelPairStepChi2ChargeMeasurementEstimator9 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPPixelPairStepChi2ChargeMeasurementEstimator9'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1e+12),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(15.0)
)


process.hltESPPixelPairStepChi2MeasurementEstimator25 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPPixelPairStepChi2MeasurementEstimator25'),
    MaxChi2 = cms.double(25.0),
    MaxDisplacement = cms.double(100.0),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(10.0),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3.0)
)


process.hltESPPixelPairTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPPixelPairTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.19)
)


process.hltESPRKTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPRKTrajectoryFitter'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.hltESPRKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPRKTrajectorySmoother'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100.0),
    minHits = cms.int32(3)
)


process.hltESPRungeKuttaTrackerPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('hltESPRungeKuttaTrackerPropagator'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(True)
)


process.hltESPSmartPropagator = cms.ESProducer("SmartPropagatorESProducer",
    ComponentName = cms.string('hltESPSmartPropagator'),
    Epsilon = cms.double(5.0),
    MuonPropagator = cms.string('hltESPSteppingHelixPropagatorAlong'),
    PropagationDirection = cms.string('alongMomentum'),
    TrackerPropagator = cms.string('PropagatorWithMaterial')
)


process.hltESPSmartPropagatorAny = cms.ESProducer("SmartPropagatorESProducer",
    ComponentName = cms.string('hltESPSmartPropagatorAny'),
    Epsilon = cms.double(5.0),
    MuonPropagator = cms.string('SteppingHelixPropagatorAny'),
    PropagationDirection = cms.string('alongMomentum'),
    TrackerPropagator = cms.string('PropagatorWithMaterial')
)


process.hltESPSmartPropagatorAnyOpposite = cms.ESProducer("SmartPropagatorESProducer",
    ComponentName = cms.string('hltESPSmartPropagatorAnyOpposite'),
    Epsilon = cms.double(5.0),
    MuonPropagator = cms.string('SteppingHelixPropagatorAny'),
    PropagationDirection = cms.string('oppositeToMomentum'),
    TrackerPropagator = cms.string('PropagatorWithMaterialOpposite')
)


process.hltESPSoftLeptonByDistance = cms.ESProducer("LeptonTaggerByDistanceESProducer",
    distance = cms.double(0.5)
)


process.hltESPSteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('hltESPSteppingHelixPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.hltESPSteppingHelixPropagatorOpposite = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('hltESPSteppingHelixPropagatorOpposite'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.hltESPStripCPEfromTrackAngle = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('hltESPStripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.hltESPTTRHBWithTrackAngle = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('hltESPTTRHBWithTrackAngle'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    PixelCPE = cms.string('hltESPPixelCPEGeneric'),
    StripCPE = cms.string('hltESPStripCPEfromTrackAngle')
)


process.hltESPTTRHBuilderAngleAndTemplate = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('hltESPTTRHBuilderAngleAndTemplate'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    PixelCPE = cms.string('hltESPPixelCPETemplateReco'),
    StripCPE = cms.string('hltESPStripCPEfromTrackAngle')
)


process.hltESPTTRHBuilderPixelOnly = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('hltESPTTRHBuilderPixelOnly'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    PixelCPE = cms.string('hltESPPixelCPEGeneric'),
    StripCPE = cms.string('Fake')
)


process.hltESPTTRHBuilderWithoutAngle4PixelTriplets = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('hltESPTTRHBuilderWithoutAngle4PixelTriplets'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    PixelCPE = cms.string('hltESPPixelCPEGeneric'),
    StripCPE = cms.string('Fake')
)


process.hltESPTobTecStepChi2ChargeMeasurementEstimator16 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('hltESPTobTecStepChi2ChargeMeasurementEstimator16'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3.0),
    pTChargeCutThreshold = cms.double(-1.0)
)


process.hltESPTobTecStepClusterShapeHitFilter = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('hltESPTobTecStepClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutTight')
    )
)


process.hltESPTobTecStepFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(False),
    ComponentName = cms.string('hltESPTobTecStepFitterSmoother'),
    EstimateCut = cms.double(30.0),
    Fitter = cms.string('hltESPTobTecStepRKFitter'),
    LogPixelProbabilityCut = cms.double(-16.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(7),
    NoInvalidHitsBeginEnd = cms.bool(False),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPTobTecStepRKSmoother'),
    appendToDataLabel = cms.string('')
)


process.hltESPTobTecStepFittingSmootherForLoopers = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(False),
    ComponentName = cms.string('hltESPTobTecStepFitterSmootherForLoopers'),
    EstimateCut = cms.double(30.0),
    Fitter = cms.string('hltESPTobTecStepRKFitterForLoopers'),
    LogPixelProbabilityCut = cms.double(-16.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(7),
    NoInvalidHitsBeginEnd = cms.bool(False),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('hltESPTobTecStepRKSmootherForLoopers'),
    appendToDataLabel = cms.string('')
)


process.hltESPTobTecStepFlexibleKFFittingSmoother = cms.ESProducer("FlexibleKFFittingSmootherESProducer",
    ComponentName = cms.string('hltESPTobTecStepFlexibleKFFittingSmoother'),
    appendToDataLabel = cms.string(''),
    looperFitter = cms.string('hltESPTobTecStepFitterSmootherForLoopers'),
    standardFitter = cms.string('hltESPTobTecStepFitterSmoother')
)


process.hltESPTobTecStepRKTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPTobTecStepRKFitter'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(7)
)


process.hltESPTobTecStepRKTrajectoryFitterForLoopers = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPTobTecStepRKFitterForLoopers'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(7)
)


process.hltESPTobTecStepRKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPTobTecStepRKSmoother'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(7)
)


process.hltESPTobTecStepRKTrajectorySmootherForLoopers = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPTobTecStepRKSmootherForLoopers'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(7)
)


process.hltESPTobTecStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPTobTecStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.09)
)


process.hltESPTrackAlgoPriorityOrder = cms.ESProducer("TrackAlgoPriorityOrderESProducer",
    ComponentName = cms.string('hltESPTrackAlgoPriorityOrder'),
    algoOrder = cms.vstring(),
    appendToDataLabel = cms.string('')
)


process.hltESPTrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer",
    appendToDataLabel = cms.string(''),
    trackerGeometryLabel = cms.untracked.string('')
)


process.hltESPTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('hltESPTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(0.0),
    ValidHitBonus = cms.double(100.0),
    allowSharedFirstHit = cms.bool(False),
    fractionShared = cms.double(0.5)
)


process.hltESPTrajectoryFitterRK = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('hltESPTrajectoryFitterRK'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.hltESPTrajectorySmootherRK = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('hltESPTrajectorySmootherRK'),
    Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
    Updator = cms.string('hltESPKFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100.0),
    minHits = cms.int32(3)
)


process.hltPixelTracksCleanerBySharedHits = cms.ESProducer("PixelTrackCleanerBySharedHitsESProducer",
    ComponentName = cms.string('hltPixelTracksCleanerBySharedHits'),
    appendToDataLabel = cms.string(''),
    useQuadrupletAlgo = cms.bool(False)
)


process.hltTrackCleaner = cms.ESProducer("TrackCleanerESProducer",
    ComponentName = cms.string('hltTrackCleaner'),
    appendToDataLabel = cms.string('')
)


process.hoDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('HODetIdAssociator'),
    etaBinSize = cms.double(0.087),
    hcalRegion = cms.int32(2),
    includeBadChambers = cms.bool(False),
    includeGEM = cms.bool(False),
    includeME0 = cms.bool(False),
    nEta = cms.int32(30),
    nPhi = cms.int32(72)
)


process.muonDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('MuonDetIdAssociator'),
    etaBinSize = cms.double(0.125),
    hcalRegion = cms.int32(2),
    includeBadChambers = cms.bool(False),
    includeGEM = cms.bool(False),
    includeME0 = cms.bool(False),
    nEta = cms.int32(48),
    nPhi = cms.int32(48)
)


process.muonSeededTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('muonSeededTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(1.0),
    ValidHitBonus = cms.double(1000.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.1)
)


process.navigationSchoolESProducer = cms.ESProducer("NavigationSchoolESProducer",
    ComponentName = cms.string('SimpleNavigationSchool'),
    SimpleMagneticField = cms.string('ParabolicMf')
)


process.preshowerDetIdAssociator = cms.ESProducer("DetIdAssociatorESProducer",
    ComponentName = cms.string('PreshowerDetIdAssociator'),
    etaBinSize = cms.double(0.1),
    hcalRegion = cms.int32(2),
    includeBadChambers = cms.bool(False),
    includeGEM = cms.bool(False),
    includeME0 = cms.bool(False),
    nEta = cms.int32(60),
    nPhi = cms.int32(30)
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    )
)


process.siPixelTemplateDBObjectESProducer = cms.ESProducer("SiPixelTemplateDBObjectESProducer")


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.CSCChannelMapperESSource = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('CSCChannelMapperRecord')
)


process.CSCINdexerESSource = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('CSCIndexerRecord')
)


process.GlobalParametersRcdSource = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('L1TGlobalParametersRcd')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(10),
        connectionRetrialTimeOut = cms.untracked.int32(60),
        connectionTimeOut = cms.untracked.int32(0),
        enableConnectionSharing = cms.untracked.bool(True),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False),
        idleConnectionCleanupPeriod = cms.untracked.int32(10),
        messageLevel = cms.untracked.int32(0)
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('94X_mc2017_realistic_TSG_2017_12_19_13_49_40'),
    pfnPostfix = cms.untracked.string('None'),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HepPDTESSource = cms.ESSource("HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/pythiaparticle.tbl')
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    fromDDD = cms.untracked.bool(False),
    toGet = cms.untracked.vstring('GainWidths')
)


process.hltESSBTagRecord = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('JetTagComputerRecord')
)


process.hltESSEcalSeverityLevel = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalSeverityLevelAlgoRcd')
)


process.hltESSHcalSeverityLevel = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('HcalSeverityLevelComputerRcd')
)


process.HLTIterativeTrackingTripletRecovery = cms.Sequence(process.hltTripletRecoveryClustersRefRemoval+process.hltTripletRecoveryMaskedMeasurementTrackerEvent+process.hltTripletRecoveryPixelLayerTriplets+process.hltTripletRecoveryPFlowPixelTrackingRegions+process.hltTripletRecoveryPFlowPixelClusterCheck+process.hltTripletRecoveryPFlowPixelHitDoublets+process.hltTripletRecoveryPFlowPixelHitTriplets+process.hltTripletRecoveryPFlowPixelSeeds+process.hltTripletRecoveryPFlowCkfTrackCandidates+process.hltTripletRecoveryPFlowCtfWithMaterialTracks+process.hltTripletRecoveryPFlowTrackCutClassifier+process.hltTripletRecoveryPFlowTrackSelectionHighPurity)


process.HLTDoLocalPixelSequence = cms.Sequence(process.hltSiPixelDigis+process.hltSiPixelClusters+process.hltSiPixelClustersCache+process.hltSiPixelRecHits)


process.HLTFastRecopixelvertexingSequence = cms.Sequence(process.hltSelector4CentralJetsL1FastJet+process.hltFastPrimaryVertex+process.hltFastPVPixelVertexFilter+process.hltFastPVPixelTracksFilter+process.hltFastPVPixelTracksFitter+process.hltFastPVPixelTracksTrackingRegions+process.hltFastPVPixelTracksHitDoublets+process.hltFastPVPixelTracksHitQuadruplets+process.hltFastPVPixelTracks+process.hltFastPVJetTracksAssociator+process.hltFastPVJetVertexChecker+process.hltFastPVPixelTracksRecoverFilter+process.hltFastPVPixelTracksRecoverFitter+process.hltFastPVPixelTracksTrackingRegionsRecover+process.hltFastPVPixelTracksHitDoubletsRecover+process.hltFastPVPixelTracksHitQuadrupletsRecover+process.hltFastPVPixelTracksRecover+process.hltFastPVPixelTracksMerger+process.hltFastPVPixelVertices+process.hltFastPVPixelVerticesFilter)


process.HLTRecoPixelTracksSequence = cms.Sequence(process.hltPixelTracksFilter+process.hltPixelTracksFitter+process.hltPixelTracksTrackingRegions+process.hltPixelLayerQuadruplets+process.hltPixelTracksHitDoublets+process.hltPixelTracksHitQuadruplets+process.hltPixelTracks+process.hltPixelTripletsClustersRefRemoval+process.hltPixelTracksTrackingRegionsForTriplets+process.hltPixelLayerTripletsWithClustersRemoval+process.hltPixelTracksHitDoubletsForTriplets+process.hltPixelTracksHitTriplets+process.hltPixelTracksFromTriplets+process.hltPixelTracksMerged)


process.HLTIterativeTrackingIteration2ForIterL3Muon = cms.Sequence(process.hltIter2IterL3MuonClustersRefRemoval+process.hltIter2IterL3MuonMaskedMeasurementTrackerEvent+process.hltIter2IterL3MuonPixelLayerTriplets+process.hltIter2IterL3MuonPixelClusterCheck+process.hltIter2IterL3MuonPixelHitDoublets+process.hltIter2IterL3MuonPixelHitTriplets+process.hltIter2IterL3MuonPixelSeeds+process.hltIter2IterL3MuonCkfTrackCandidates+process.hltIter2IterL3MuonCtfWithMaterialTracks+process.hltIter2IterL3MuonTrackCutClassifier+process.hltIter2IterL3MuonTrackSelectionHighPurity)


process.HLTL1UnpackerSequence = cms.Sequence(process.hltGtStage2Digis+process.hltGtStage2ObjectMap)


process.HLTMuonLocalRecoSequence = cms.Sequence(process.hltMuonDTDigis+process.hltDt1DRecHits+process.hltDt4DSegments+process.hltMuonCSCDigis+process.hltCsc2DRecHits+process.hltCscSegments+process.hltMuonRPCDigis+process.hltRpcRecHits)


process.HLTDoLocalHcalSequence = cms.Sequence(process.hltHcalDigis+process.hltHbhePhase1Reco+process.hltHbhereco+process.hltHfprereco+process.hltHfreco+process.hltHoreco)


process.HLTIter1TrackAndTauJets4Iter2Sequence = cms.Sequence(process.hltIter1TrackRefsForJets4Iter2+process.hltAK4Iter1TrackJets4Iter2+process.hltIter1TrackAndTauJets4Iter2)


process.HLTBeamSpot = cms.Sequence(process.hltScalersRawToDigi+process.hltOnlineBeamSpot)


process.HLTL2muonrecoNocandSequence = cms.Sequence(process.HLTMuonLocalRecoSequence+process.hltL2OfflineMuonSeeds+process.hltL2MuonSeeds+process.hltL2Muons)


process.HLTBtagCSVSequencePFAK8 = cms.Sequence(process.hltVerticesPF+process.hltVerticesPFSelector+process.hltVerticesPFFilter+process.hltPFJetForBtagSelectorAK8+process.hltPFJetForBtagAK8+process.hltBLifetimeAssociatorPFAK8+process.hltBLifetimeTagInfosPFAK8+process.hltInclusiveVertexFinderPF+process.hltInclusiveSecondaryVerticesPF+process.hltTrackVertexArbitratorPF+process.hltInclusiveMergedVerticesPF+process.hltSecondaryVertexTagInfosPFAK8+process.hltCombinedSecondaryVertexBJetTagsPFAK8+process.hltBoostedDBJetTagsPFAK8)


process.HLTAK4CaloCorrectorProducersSequence = cms.Sequence(process.hltAK4CaloFastJetCorrector+process.hltAK4CaloRelativeCorrector+process.hltAK4CaloAbsoluteCorrector+process.hltAK4CaloResidualCorrector+process.hltAK4CaloCorrector)


process.HLTBtagBoostedDoubleSecondaryVertexSequencePFAK8 = cms.Sequence(process.hltVerticesPF+process.hltVerticesPFSelector+process.hltVerticesPFFilter+process.hltPFJetForBtagSelectorAK8+process.hltPFJetForBtagAK8+process.hltBoostedDBLifetimeTagInfosPFAK8+process.hltBoostedDBInclusiveVertexFinderPFAK8+process.hltBoostedDBInclusiveSecondaryVerticesPFAK8+process.hltBoostedDBTrackVertexArbitratorPFAK8+process.hltBoostedDBInclusiveMergedVerticesPFAK8+process.hltBoostedDBSecondaryVertexTagInfosPFAK8+process.hltBoostedDBSVAK8TagInfosPF+process.hltBoostedDBJetTagsPFAK8)


process.HLTIterativeTrackingDoubletRecovery = cms.Sequence(process.hltDoubletRecoveryClustersRefRemoval+process.hltDoubletRecoveryMaskedMeasurementTrackerEvent+process.hltDoubletRecoveryPixelLayerPairs+process.hltDoubletRecoveryPFlowPixelTrackingRegions+process.hltDoubletRecoveryPFlowPixelClusterCheck+process.hltDoubletRecoveryPFlowPixelHitDoublets+process.hltDoubletRecoveryPFlowPixelSeeds+process.hltDoubletRecoveryPFlowCkfTrackCandidates+process.hltDoubletRecoveryPFlowCtfWithMaterialTracks+process.hltDoubletRecoveryPFlowTrackCutClassifier+process.hltDoubletRecoveryPFlowTrackSelectionHighPurity)


process.HLTEndSequence = cms.Sequence(process.hltBoolEnd)


process.HLTIterativeTrackingIteration2 = cms.Sequence(process.hltIter2ClustersRefRemoval+process.hltIter2MaskedMeasurementTrackerEvent+process.hltIter2PixelLayerTriplets+process.hltIter2PFlowPixelTrackingRegions+process.hltIter2PFlowPixelClusterCheck+process.hltIter2PFlowPixelHitDoublets+process.hltIter2PFlowPixelHitTriplets+process.hltIter2PFlowPixelSeeds+process.hltIter2PFlowCkfTrackCandidates+process.hltIter2PFlowCtfWithMaterialTracks+process.hltIter2PFlowTrackCutClassifier+process.hltIter2PFlowTrackSelectionHighPurity)


process.HLTIterativeTrackingIteration1ForBTag = cms.Sequence(process.hltIter1ClustersRefRemovalForBTag+process.hltIter1MaskedMeasurementTrackerEventForBTag+process.hltIter1PixelLayerQuadrupletsForBTag+process.hltIter1PFlowPixelTrackingRegionsForBTag+process.hltIter1PFlowPixelClusterCheckForBTag+process.hltIter1PFlowPixelHitDoubletsForBTag+process.hltIter1PFlowPixelHitQuadrupletsForBTag+process.hltIter1PFlowPixelSeedsForBTag+process.hltIter1PFlowCkfTrackCandidatesForBTag+process.hltIter1PFlowCtfWithMaterialTracksForBTag+process.hltIter1PFlowTrackCutClassifierPromptForBTag+process.hltIter1PFlowTrackCutClassifierDetachedForBTag+process.hltIter1PFlowTrackCutClassifierMergedForBTag+process.hltIter1PFlowTrackSelectionHighPurityForBTag)


process.HLTDoLocalStripSequence = cms.Sequence(process.hltSiStripExcludedFEDListProducer+process.hltSiStripRawToClustersFacility+process.hltSiStripClusters)


process.HLTIterL3OImuonTkCandidateSequence = cms.Sequence(process.hltIterL3OISeedsFromL2Muons+process.hltIterL3OITrackCandidates+process.hltIterL3OIMuCtfWithMaterialTracks+process.hltIterL3OIMuonTrackCutClassifier+process.hltIterL3OIMuonTrackSelectionHighPurity+process.hltL3MuonsIterL3OI)


process.HLTRecopixelvertexingSequence = cms.Sequence(process.HLTRecoPixelTracksSequence+process.hltPixelVertices+process.hltTrimmedPixelVertices)


process.HLTBeginSequence = cms.Sequence(process.hltTriggerType+process.HLTL1UnpackerSequence+process.HLTBeamSpot)


process.HLTRecoPixelTracksSequenceForIterL3FromL1Muon = cms.Sequence(process.hltIterL3FromL1MuonPixelTracksTrackingRegions+process.hltIterL3FromL1MuonPixelLayerQuadruplets+process.hltIterL3FromL1MuonPixelTracksHitDoublets+process.hltIterL3FromL1MuonPixelTracksHitQuadruplets+process.hltIterL3FromL1MuonPixelTracks)


process.HLTPreshowerSequence = cms.Sequence(process.hltEcalPreshowerDigis+process.hltEcalPreshowerRecHit)


process.HLTDoLocalPixelSequenceRegForBTag = cms.Sequence(process.hltSelectorJets20L1FastJet+process.hltSelectorCentralJets20L1FastJeta+process.hltSiPixelDigisRegForBTag+process.hltSiPixelClustersRegForBTag+process.hltSiPixelClustersRegForBTagCache+process.hltSiPixelRecHitsRegForBTag+process.hltPixelLayerQuadrupletsRegForBTag)


process.HLTIterL3MuonRecoPixelTracksSequence = cms.Sequence(process.hltIterL3MuonPixelTracksFilter+process.hltIterL3MuonPixelTracksFitter+process.hltIterL3MuonPixelTracksTrackingRegions+process.hltIterL3MuonPixelLayerQuadruplets+process.hltIterL3MuonPixelTracksHitDoublets+process.hltIterL3MuonPixelTracksHitQuadruplets+process.hltIterL3MuonPixelTracks)


process.HLTAK8CaloCorrectorProducersSequence = cms.Sequence(process.hltAK8CaloFastJetCorrector+process.hltAK8CaloRelativeCorrector+process.hltAK8CaloAbsoluteCorrector+process.hltAK8CaloResidualCorrector+process.hltAK8CaloCorrector)


process.HLTL2muonrecoSequence = cms.Sequence(process.HLTL2muonrecoNocandSequence+process.hltL2MuonCandidates)


process.HLTBtagCSVSequencePF = cms.Sequence(process.hltVerticesPF+process.hltVerticesPFSelector+process.hltVerticesPFFilter+process.hltPFJetForBtagSelector+process.hltPFJetForBtag+process.hltBLifetimeAssociatorPF+process.hltBLifetimeTagInfosPF+process.hltInclusiveVertexFinderPF+process.hltInclusiveSecondaryVerticesPF+process.hltTrackVertexArbitratorPF+process.hltInclusiveMergedVerticesPF+process.hltSecondaryVertexTagInfosPF+process.hltCombinedSecondaryVertexBJetTagsPF)


process.HLTIterativeTrackingIteration0ForBTag = cms.Sequence(process.hltIter0PFLowPixelSeedsFromPixelTracksForBTag+process.hltIter0PFlowCkfTrackCandidatesForBTag+process.hltIter0PFlowCtfWithMaterialTracksForBTag+process.hltIter0PFlowTrackCutClassifierForBTag+process.hltIter0PFlowTrackSelectionHighPurityForBTag)


process.HLTRecopixelvertexingSequenceForIterL3FromL1Muon = cms.Sequence(process.HLTRecoPixelTracksSequenceForIterL3FromL1Muon+process.hltIterL3FromL1MuonPixelVertices+process.hltIterL3FromL1MuonTrimmedPixelVertices)


process.HLTIterativeTrackingIteration0ForIterL3FromL1Muon = cms.Sequence(process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks+process.hltIter0IterL3FromL1MuonCkfTrackCandidates+process.hltIter0IterL3FromL1MuonCtfWithMaterialTracks+process.hltIter0IterL3FromL1MuonTrackCutClassifier+process.hltIter0IterL3FromL1MuonTrackSelectionHighPurity)


process.HLTIterativeTrackingIteration0ForIterL3Muon = cms.Sequence(process.hltIter0IterL3MuonPixelSeedsFromPixelTracks+process.hltIter0IterL3MuonCkfTrackCandidates+process.hltIter0IterL3MuonCtfWithMaterialTracks+process.hltIter0IterL3MuonTrackCutClassifier+process.hltIter0IterL3MuonTrackSelectionHighPurity)


process.HLTAK8PFCorrectorProducersSequence = cms.Sequence(process.hltAK8PFFastJetCorrector+process.hltAK8PFRelativeCorrector+process.hltAK8PFAbsoluteCorrector+process.hltAK8PFResidualCorrector+process.hltAK8PFCorrector)


process.HLTIterativeTrackingIteration2ForIterL3FromL1Muon = cms.Sequence(process.hltIter2IterL3FromL1MuonClustersRefRemoval+process.hltIter2IterL3FromL1MuonMaskedMeasurementTrackerEvent+process.hltIter2IterL3FromL1MuonPixelLayerTriplets+process.hltIter2IterL3FromL1MuonPixelClusterCheck+process.hltIter2IterL3FromL1MuonPixelHitDoublets+process.hltIter2IterL3FromL1MuonPixelHitTriplets+process.hltIter2IterL3FromL1MuonPixelSeeds+process.hltIter2IterL3FromL1MuonCkfTrackCandidates+process.hltIter2IterL3FromL1MuonCtfWithMaterialTracks+process.hltIter2IterL3FromL1MuonTrackCutClassifier+process.hltIter2IterL3FromL1MuonTrackSelectionHighPurity)


process.HLTIterativeTrackingIteration2ForBTag = cms.Sequence(process.hltIter2ClustersRefRemovalForBTag+process.hltIter2MaskedMeasurementTrackerEventForBTag+process.hltIter2PixelLayerTripletsForBTag+process.hltIter2PFlowPixelTrackingRegionsForBTag+process.hltIter2PFlowPixelClusterCheckForBTag+process.hltIter2PFlowPixelHitDoubletsForBTag+process.hltIter2PFlowPixelHitTripletsForBTag+process.hltIter2PFlowPixelSeedsForBTag+process.hltIter2PFlowCkfTrackCandidatesForBTag+process.hltIter2PFlowCtfWithMaterialTracksForBTag+process.hltIter2PFlowTrackCutClassifierForBTag+process.hltIter2PFlowTrackSelectionHighPurityForBTag)


process.HLTBtagDeepCSVSequencePF = cms.Sequence(process.hltVerticesPF+process.hltVerticesPFSelector+process.hltVerticesPFFilter+process.hltPFJetForBtagSelector+process.hltPFJetForBtag+process.hltDeepBLifetimeTagInfosPF+process.hltDeepInclusiveVertexFinderPF+process.hltDeepInclusiveSecondaryVerticesPF+process.hltDeepTrackVertexArbitratorPF+process.hltDeepInclusiveMergedVerticesPF+process.hltDeepSecondaryVertexTagInfosPF+process.hltDeepCombinedSecondaryVertexBJetTagsInfos+process.hltDeepCombinedSecondaryVertexBJetTagsPF)


process.HLTAK4PFCorrectorProducersSequence = cms.Sequence(process.hltAK4PFFastJetCorrector+process.hltAK4PFRelativeCorrector+process.hltAK4PFAbsoluteCorrector+process.hltAK4PFResidualCorrector+process.hltAK4PFCorrector)


process.HLTIterativeTrackingIteration1 = cms.Sequence(process.hltIter1ClustersRefRemoval+process.hltIter1MaskedMeasurementTrackerEvent+process.hltIter1PixelLayerQuadruplets+process.hltIter1PFlowPixelTrackingRegions+process.hltIter1PFlowPixelClusterCheck+process.hltIter1PFlowPixelHitDoublets+process.hltIter1PFlowPixelHitQuadruplets+process.hltIter1PixelTracks+process.hltIter1PFLowPixelSeedsFromPixelTracks+process.hltIter1PFlowCkfTrackCandidates+process.hltIter1PFlowCtfWithMaterialTracks+process.hltIter1PFlowTrackCutClassifierPrompt+process.hltIter1PFlowTrackCutClassifierDetached+process.hltIter1PFlowTrackCutClassifierMerged+process.hltIter1PFlowTrackSelectionHighPurity)


process.HLTIterativeTrackingIteration0 = cms.Sequence(process.hltIter0PFLowPixelSeedsFromPixelTracks+process.hltIter0PFlowCkfTrackCandidates+process.hltIter0PFlowCtfWithMaterialTracks+process.hltIter0PFlowTrackCutClassifier+process.hltIter0PFlowTrackSelectionHighPurity)


process.HLTIterativeTrackingIter02 = cms.Sequence(process.HLTIterativeTrackingIteration0+process.HLTIterativeTrackingIteration1+process.hltIter1Merged+process.HLTIter1TrackAndTauJets4Iter2Sequence+process.HLTIterativeTrackingIteration2+process.hltIter2Merged+process.HLTIterativeTrackingTripletRecovery+process.hltTripletRecoveryMerged+process.HLTIterativeTrackingDoubletRecovery+process.hltMergedTracks)


process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence = cms.Sequence(process.hltEcalDigis+process.hltEcalUncalibRecHit+process.hltEcalDetIdToBeRecovered+process.hltEcalRecHit)


process.HLTIterativeTrackingIter02ForBTag = cms.Sequence(process.HLTIterativeTrackingIteration0ForBTag+process.HLTIterativeTrackingIteration1ForBTag+process.hltIter1MergedForBTag+process.HLTIterativeTrackingIteration2ForBTag+process.hltIter2MergedForBTag)


process.HLTDoLocalStripSequenceRegForBTag = cms.Sequence(process.hltSiStripExcludedFEDListProducer+process.hltSiStripRawToClustersFacility+process.hltSiStripClustersRegForBTag)


process.HLTTrackReconstructionForPF = cms.Sequence(process.HLTDoLocalPixelSequence+process.HLTRecopixelvertexingSequence+process.HLTDoLocalStripSequence+process.HLTIterativeTrackingIter02+process.hltPFMuonMerging+process.hltMuonLinks+process.hltMuons)


process.HLTParticleFlowSequence = cms.Sequence(process.HLTPreshowerSequence+process.hltParticleFlowRecHitECALUnseeded+process.hltParticleFlowRecHitHBHE+process.hltParticleFlowRecHitHF+process.hltParticleFlowRecHitPSUnseeded+process.hltParticleFlowClusterECALUncorrectedUnseeded+process.hltParticleFlowClusterPSUnseeded+process.hltParticleFlowClusterECALUnseeded+process.hltParticleFlowClusterHBHE+process.hltParticleFlowClusterHCAL+process.hltParticleFlowClusterHF+process.hltLightPFTracks+process.hltParticleFlowBlock+process.hltParticleFlow)


process.HLTIterL3MuonRecopixelvertexingSequence = cms.Sequence(process.HLTIterL3MuonRecoPixelTracksSequence+process.hltIterL3MuonPixelVertices+process.hltIterL3MuonTrimmedPixelVertices)


process.HLTIterativeTrackingIter02ForIterL3Muon = cms.Sequence(process.HLTIterativeTrackingIteration0ForIterL3Muon+process.HLTIterativeTrackingIteration2ForIterL3Muon+process.hltIter2IterL3MuonMerged)


process.HLTAK4CaloJetsCorrectionSequence = cms.Sequence(process.hltFixedGridRhoFastjetAllCalo+process.HLTAK4CaloCorrectorProducersSequence+process.hltAK4CaloJetsCorrected+process.hltAK4CaloJetsCorrectedIDPassed)


process.HLTAK8PFJetsCorrectionSequence = cms.Sequence(process.hltFixedGridRhoFastjetAll+process.HLTAK8PFCorrectorProducersSequence+process.hltAK8PFJetsCorrected+process.hltAK8PFJetsLooseIDCorrected+process.hltAK8PFJetsTightIDCorrected)


process.HLTTrackReconstructionForBTag = cms.Sequence(process.HLTDoLocalPixelSequenceRegForBTag+process.HLTFastRecopixelvertexingSequence+process.HLTDoLocalStripSequenceRegForBTag+process.HLTIterativeTrackingIter02ForBTag)


process.HLTDoCaloSequence = cms.Sequence(process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence+process.HLTDoLocalHcalSequence+process.hltTowerMakerForAll)


process.HLTAK4PFJetsCorrectionSequence = cms.Sequence(process.hltFixedGridRhoFastjetAll+process.HLTAK4PFCorrectorProducersSequence+process.hltAK4PFJetsCorrected+process.hltAK4PFJetsLooseIDCorrected+process.hltAK4PFJetsTightIDCorrected)


process.HLTIterL3IOmuonTkCandidateSequence = cms.Sequence(process.HLTIterL3MuonRecopixelvertexingSequence+process.HLTIterativeTrackingIter02ForIterL3Muon+process.hltL3MuonsIterL3IO)


process.HLTAK4CaloJetsReconstructionSequence = cms.Sequence(process.HLTDoCaloSequence+process.hltAK4CaloJets+process.hltAK4CaloJetsIDPassed)


process.HLTBtagCSVSequenceL3 = cms.Sequence(process.hltSelectorJets30L1FastJet+process.hltSelectorCentralJets30L1FastJeta+process.hltSelector8CentralJetsL1FastJet+process.HLTTrackReconstructionForBTag+process.hltVerticesL3+process.hltFastPixelBLifetimeL3Associator+process.hltImpactParameterTagInfos+process.hltInclusiveVertexFinder+process.hltInclusiveSecondaryVertices+process.hltTrackVertexArbitrator+process.hltInclusiveMergedVertices+process.hltInclusiveSecondaryVertexFinderTagInfos+process.hltCombinedSecondaryVertexBJetTagsCalo)


process.HLTAK8CaloJetsCorrectionSequence = cms.Sequence(process.hltFixedGridRhoFastjetAllCalo+process.HLTAK8CaloCorrectorProducersSequence+process.hltAK8CaloJetsCorrected+process.hltAK8CaloJetsCorrectedIDPassed)


process.HLTIterativeTrackingIter02ForIterL3FromL1Muon = cms.Sequence(process.HLTIterativeTrackingIteration0ForIterL3FromL1Muon+process.HLTIterativeTrackingIteration2ForIterL3FromL1Muon+process.hltIter2IterL3FromL1MuonMerged)


process.HLTDoCaloSequencePF = cms.Sequence(process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence+process.HLTDoLocalHcalSequence+process.hltTowerMakerForAll)


process.HLTAK8CaloJetsPrePFRecoSequence = cms.Sequence(process.HLTDoCaloSequencePF+process.hltAK8CaloJetsPF+process.hltAK4CaloJetsPF)


process.HLTBtagDeepCSVSequenceL3 = cms.Sequence(process.hltSelectorJets30L1FastJet+process.hltSelectorCentralJets30L1FastJeta+process.hltSelector8CentralJetsL1FastJet+process.HLTTrackReconstructionForBTag+process.hltVerticesL3+process.hltFastPixelBLifetimeL3Associator+process.hltImpactParameterTagInfos+process.hltInclusiveVertexFinder+process.hltInclusiveSecondaryVertices+process.hltTrackVertexArbitrator+process.hltInclusiveMergedVertices+process.hltInclusiveSecondaryVertexFinderTagInfos+process.hltDeepCombinedSecondaryVertexBJetTagsInfosCalo+process.hltDeepCombinedSecondaryVertexBJetTagsCalo)


process.HLTPreAK8PFJetsRecoSequence = cms.Sequence(process.HLTAK8CaloJetsPrePFRecoSequence+process.hltAK8CaloJetsPFEt5+process.hltAK4CaloJetsPFEt5)


process.HLTAK4CaloJetsPrePFRecoSequence = cms.Sequence(process.HLTDoCaloSequencePF+process.hltAK4CaloJetsPF)


process.HLTAK4CaloJetsSequence = cms.Sequence(process.HLTAK4CaloJetsReconstructionSequence+process.HLTAK4CaloJetsCorrectionSequence)


process.HLTIterL3IOmuonFromL1TkCandidateSequence = cms.Sequence(process.HLTRecopixelvertexingSequenceForIterL3FromL1Muon+process.HLTIterativeTrackingIter02ForIterL3FromL1Muon)


process.HLTAK8CaloJetsReconstructionSequence = cms.Sequence(process.HLTDoCaloSequence+process.hltAK8CaloJets+process.hltAK8CaloJetsIDPassed)


process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence = cms.Sequence(process.HLTIterL3OImuonTkCandidateSequence+process.hltIterL3OIL3MuonsLinksCombination+process.hltIterL3OIL3Muons+process.hltIterL3OIL3MuonCandidates+process.hltL2SelectorForL3IO+process.HLTIterL3IOmuonTkCandidateSequence+process.hltIterL3MuonsFromL2LinksCombination+process.hltIterL3MuonsFromL2)


process.HLTIterL3muonTkCandidateSequence = cms.Sequence(process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence+process.hltIterL3MuonL1MuonNoL2Selector+process.HLTIterL3IOmuonFromL1TkCandidateSequence)


process.HLTPreAK4PFJetsRecoSequence = cms.Sequence(process.HLTAK4CaloJetsPrePFRecoSequence+process.hltAK4CaloJetsPFEt5)


process.HLTAK8CaloJetsSequence = cms.Sequence(process.HLTAK8CaloJetsReconstructionSequence+process.HLTAK8CaloJetsCorrectionSequence)


process.HLTL3muonrecoNocandSequence = cms.Sequence(process.HLTIterL3muonTkCandidateSequence+process.hltIterL3MuonMerged+process.hltIterL3MuonAndMuonFromL1Merged+process.hltL3MuonsIterL3Links+process.hltIterL3Muons)


process.HLTL3muonrecoSequence = cms.Sequence(process.HLTL3muonrecoNocandSequence+process.hltIterL3MuonCandidates)


process.HLTAK8PFJetsReconstructionSequence = cms.Sequence(process.HLTL2muonrecoSequence+process.HLTL3muonrecoSequence+process.HLTTrackReconstructionForPF+process.HLTParticleFlowSequence+process.hltAK8PFJets+process.hltAK8PFJetsLooseID+process.hltAK8PFJetsTightID)


process.HLTAK8PFJetsSequence = cms.Sequence(process.HLTPreAK8PFJetsRecoSequence+process.HLTAK8PFJetsReconstructionSequence+process.HLTAK8PFJetsCorrectionSequence)


process.HLTAK4PFJetsReconstructionSequence = cms.Sequence(process.HLTL2muonrecoSequence+process.HLTL3muonrecoSequence+process.HLTTrackReconstructionForPF+process.HLTParticleFlowSequence+process.hltAK4PFJets+process.hltAK4PFJetsLooseID+process.hltAK4PFJetsTightID)


process.HLTAK4PFJetsSequence = cms.Sequence(process.HLTPreAK4PFJetsRecoSequence+process.HLTAK4PFJetsReconstructionSequence+process.HLTAK4PFJetsCorrectionSequence)


process.HLTriggerFirstPath = cms.Path(process.hltGetConditions+process.hltGetRaw+process.hltBoolFalse)


process.HLT_AK8PFJet400_v14 = cms.Path(process.HLTBeginSequence+process.hltL1sSingleJet170IorSingleJet180IorSingleJet200+process.hltPreAK8PFJet400+process.HLTAK8CaloJetsSequence+process.hltSingleCaloJet350AK8+process.HLTAK8PFJetsSequence+process.hltPFJetsCorrectedMatchedToCaloJets350AK8+process.hltSinglePFJet400AK8+process.HLTEndSequence)


process.HLT_AK8PFJet330_PFAK8BTagCSV_p1_v3 = cms.Path(process.HLTBeginSequence+process.hltL1sSingleJet170IorSingleJet180IorSingleJet200+process.hltPreAK8PFJetXPFAK8DoubleBpY+process.HLTAK8CaloJetsSequence+process.hltAK8SingleCaloJet300+process.HLTAK8PFJetsSequence+process.hltAK8PFJetsCorrectedMatchedToCaloJets300+process.hltSinglePFJet330AK8+process.HLTBtagCSVSequencePFAK8+process.hltBTagPFCSVp4Single+process.HLTEndSequence)


process.HLT_AK8PFJetX_PFAK8DoubleB_pY_v1 = cms.Path(process.HLTBeginSequence+process.hltL1sSingleJet170IorSingleJet180IorSingleJet200+process.hltPreAK8PFJetXPFAK8DoubleBpY+process.HLTAK8CaloJetsSequence+process.hltAK8SingleCaloJet300+process.HLTAK8PFJetsSequence+process.hltAK8PFJetsCorrectedMatchedToCaloJets300+process.hltSinglePFJet330AK8+process.HLTBtagBoostedDoubleSecondaryVertexSequencePFAK8+process.hltBTagPFBoostedDoubleSecondaryVertexpZSingle+process.HLTEndSequence)


process.HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v7 = cms.Path(process.HLTBeginSequence+process.hltL1sHTT280to500erIorHTT250to340erQuadJet+process.hltPrePFHT380SixPFJet32DoublePFBTagCSV2p2+process.HLTAK4CaloJetsSequence+process.hltCaloJetFilterSixC25+process.hltCaloJetsSix25ForHt+process.hltHtMhtCaloJetsSixC25+process.hltCaloSixJet25HT300+process.HLTBtagCSVSequenceL3+process.hltBTagCaloCSVp44Single+process.HLTAK4PFJetsSequence+process.hltPFJetFilterSixC32+process.hltPFJetsSix32ForHt+process.hltHtMhtPFJetsSixC32+process.hltPFSixJet32HT380+process.HLTBtagCSVSequencePF+process.hltBTagPFCSVp075Double+process.HLTEndSequence)


process.HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v6 = cms.Path(process.HLTBeginSequence+process.hltL1sHTT280to500erIorHTT250to340erQuadJet+process.hltPrePFHT380SixPFJet32DoublePFBTagDeepCSV2p2+process.HLTAK4CaloJetsSequence+process.hltCaloJetFilterSixC25+process.hltCaloJetsSix25ForHt+process.hltHtMhtCaloJetsSixC25+process.hltCaloSixJet25HT300+process.HLTBtagDeepCSVSequenceL3+process.hltBTagCaloDeepCSVp44Single+process.HLTAK4PFJetsSequence+process.hltPFJetFilterSixC32+process.hltPFJetsSix32ForHt+process.hltHtMhtPFJetsSixC32+process.hltPFSixJet32HT380+process.HLTBtagDeepCSVSequencePF+process.hltBTagPFDeepCSVp075Double+process.HLTEndSequence)


process.HLTriggerFinalPath = cms.Path(process.hltGtStage2Digis+process.hltScalersRawToDigi+process.hltFEDSelector+process.hltTriggerSummaryAOD+process.hltTriggerSummaryRAW+process.hltBoolFalse)


process.DQMOutput = cms.EndPath(process.dqmOutput)


process.Output = cms.EndPath(process.hltOutput)


process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_AK8PFJet400_v14, process.HLT_AK8PFJet330_PFAK8BTagCSV_p1_v3, process.HLT_AK8PFJetX_PFAK8DoubleB_pY_v1, process.HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v7, process.HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v6, process.HLTriggerFinalPath, process.Output ])

