import FWCore.ParameterSet.Config as cms

JetSubCalc = cms.PSet(
    dummy_parameter = cms.string('dummy'),
    CA8TopJetColl = cms.InputTag("goodPatJetsCATopTagPFPacked"),
    CA8HEPTopJetColl = cms.InputTag("goodPatJetsCAHEPTopTagPFPacked"),
    CA8PrunedJetColl = cms.InputTag("goodPatJetsCA8PrunedPFPacked"),
    CA8JetColl = cms.InputTag("goodPatJetsCA8PF"),
    bDiscriminant     = cms.string("combinedSecondaryVertexBJetTags"),
    CA8TopTagInfo	  = cms.string("CATop")

)