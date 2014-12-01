import FWCore.ParameterSet.Config as cms

JetSubCalc = cms.PSet(
    dummy_parameter = cms.string('dummy'),
    CA8PrunedJetColl = cms.InputTag("goodPatJetsCA8PrunedPFPacked"),
    debug          = cms.bool(False)
)