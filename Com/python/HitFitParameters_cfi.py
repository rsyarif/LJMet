import FWCore.ParameterSet.Config as cms

defaultHitFitParameters = cms.PSet(
    isData                   = cms.bool(False),  # is sample data?
    hitFitLepId              = cms.int32(13),      # which lepton, 13 - muon or  11 - electron? # maybe this should be done automatically?
    hitfitDebug              = cms.untracked.bool(False),
    hitfitDefault            = cms.untracked.FileInPath("TopQuarkAnalysis/TopHitFit/data/setting/RunHitFitConfiguration.txt"),
    hitfitElectronResolution = cms.untracked.FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafElectronResolution.txt"),
    hitfitElectronObjRes     = cms.untracked.bool(False),
    hitfitMuonResolution     = cms.untracked.FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafMuonResolution.txt"),
    hitfitMuonObjRes         = cms.untracked.bool(False),
    hitfitUdscJetResolution  = cms.untracked.FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafUdscJetResolution.txt"),
    hitfitBJetResolution     = cms.untracked.FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafBJetResolution.txt"),
    hitfitJetObjRes          = cms.untracked.bool(False),
    hitfitJetCorrectionLevel = cms.untracked.string('L3Absolute'),#3Absolute'),#L7Parton'),
    #hitfitUdscJES            = cms.double(1.0),
    #hitfitBJES               = cms.double(1.0),
    hitfitMETResolution      = cms.untracked.FileInPath("TopQuarkAnalysis/TopHitFit/data/resolution/tqafKtResolution.txt"),
    hitfitMETsObjRes         = cms.untracked.bool(False),
    hitfitLepWMass           = cms.untracked.double(80.4),
    hitfitHadWMass           = cms.untracked.double(80.4),
    hitfitTopMass            = cms.untracked.double(0.0),
    hitfitNuSolution         = cms.untracked.int32(2),
    hitfitMinLeptonPt        = cms.untracked.double(15.0),
    hitfitMinJetPt           = cms.untracked.double(15.0),
    hitfitMinMET             = cms.untracked.double(0.0),
    hitfitUseNLeadJets       = cms.untracked.bool(False),
    hitfitMaxNJet            = cms.untracked.uint32(5)
    )
