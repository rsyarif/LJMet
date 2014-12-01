import FWCore.ParameterSet.Config as cms


electronSelector = cms.PSet(
    version = cms.string('DATA2010'),
    #quality = cms.string('OFF'),
    quality = cms.string('TOP_VERYLOOSE'),
    #quality = cms.string('EWK_Wenu'),
    electronId = cms.string('someElectronId'),
    minEt = cms.double(15.0),
    eta = cms.double(2.5),
    relIso = cms.double(0.2)
)    

