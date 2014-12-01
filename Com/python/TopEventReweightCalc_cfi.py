import FWCore.ParameterSet.Config as cms


TopEventReweightCalc = cms.PSet (
    isMc         = cms.bool(True),
    genParticles = cms.InputTag("genParticles"),
    ptReweight = cms.bool(False),
    reweightBSemiLeptDecyas = cms.bool(False),
    nuDecayFractionSource = cms.double(0.25),
    nuDecayFractionTarget = cms.double(0.239),
    reweightBfragmentation = cms.bool(False),
    fragSourceFile = cms.string('/src/LJMet/Dilepton/data/TopUtils_data_MC_BJES_TuneZ2star.root'),
    fragTargetFile = cms.string('/src/LJMet/Dilepton/data/TopUtils_data_MC_BJES_TuneZ2star_rbLEP.root')
)
