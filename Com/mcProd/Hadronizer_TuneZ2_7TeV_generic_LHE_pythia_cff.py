#!/usr/bin/env python
#########################################################################
#
# lhe_to_cmsDriver_cfg.py
#
# CMSSW config file for supplying cmsDriver with info to prepare
# simulation config files
#
# Usage:
#
# Gena Kukartsev, February 16, 2012
#

import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PythiaUEZ2Settings_cfi import *

generator = cms.EDFilter("Pythia6HadronizerFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    comEnergy = cms.double(7000.0),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring('MSEL=0          ! User defined processes',
                                        'PMAS(5,1)=4.4   ! b quark mass',
                                        'PMAS(6,1)=172.4 ! t quark mass'),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings',
                                    'processParameters'
                                    )
        )
)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string ('$Source: /local/reps/CMSSW/UserCode/LJMet/Com/python/Hadronizer_TuneZ2_7TeV_generic_LHE_pythia_cff.py,v $'),
    annotation = cms.untracked.string('runs Z2 Pythia6')
)
