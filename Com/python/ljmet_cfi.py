import FWCore.ParameterSet.Config as cms

ljmet = cms.PSet(
    isMc      = cms.bool(True),
    verbosity = cms.int32(0),
    useHcalLaserEventFilter = cms.bool(True),
    runs                 = cms.vint32([]),
    excluded_calculators = cms.vstring()
    )
