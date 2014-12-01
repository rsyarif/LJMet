import FWCore.ParameterSet.Config as cms

StopCalc = cms.PSet(
    dummy_parameter = cms.string('dummy'),
    PdfInfoTag = cms.untracked.InputTag("generator"),
    PdfSetNames = cms.untracked.vstring(
        "cteq66.LHgrid"
        , "MRST2006nnlo.LHgrid"
        , "NNPDF10_100.LHgrid"
        )
    )
