import FWCore.ParameterSet.Config as cms

PdfCalc = cms.PSet (
    PdfSetNames = cms.untracked.vstring("cteq66.LHgrid"),
    reducedInfo = cms.untracked.bool(True),
    defaultPdfSet = cms.untracked.string("cteq66.LHgrid")
)
