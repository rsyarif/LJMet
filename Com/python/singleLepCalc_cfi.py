import FWCore.ParameterSet.Config as cms

singleLepCalc = cms.PSet(
    triggerSummary = cms.InputTag("hltTriggerSummaryAOD"),
    rhoSrc     = cms.InputTag("kt6PFJets", 'rho'),
    isWJets     = cms.bool(False),
    isTB       = cms.bool(False)
)
