import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes



#
# Parameters that can be set via command line
# when submitting Condor jobs
#
isMc_settable       = True
isSignalMc_settable = False




def FindFile(name):

    fname = 'file.txt'

    return fname



process = cms.Process("LJMetCom")



##################################################################
#
# All input files needed for the job to run
# Specify them here, and they will automatically be correctly
# transferred to Condor when needed
# NOTE: you can define as many or as few entries as you wish,
#       names are up to you
miscFiles = {}
miscFiles['jec_uncertainty']  = '../cond/Summer12_V2_DATA_AK5PF_UncertaintySources.txt'
miscFiles['btag_performance'] = '../cond/btag_performance_db062012.root'
miscFiles['json']             = '../data/json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
miscFiles['MCL1JetPar']       = '../data/START53_V7G_L1FastJet_AK5PFchs.txt'
miscFiles['MCL2JetPar']       = '../data/START53_V7G_L2Relative_AK5PFchs.txt'
miscFiles['MCL3JetPar']       = '../data/START53_V7G_L3Absolute_AK5PFchs.txt'
miscFiles['DataL1JetPar']     = '../data/FT_53_V10_AN3_L1FastJet_AK5PFchs.txt'
miscFiles['DataL2JetPar']     = '../data/FT_53_V10_AN3_L2Relative_AK5PFchs.txt'
miscFiles['DataL3JetPar']     = '../data/FT_53_V10_AN3_L3Absolute_AK5PFchs.txt'
miscFiles['DataResJetPar']    = '../data/FT_53_V10_AN3_L2L3Residual_AK5PFchs.txt'



############################################################
#
# FWLite application options
process.load('LJMet.Com.ljmet_cfi')
process.ljmet.isMc      = cms.bool(isMc_settable)
process.ljmet.verbosity = cms.int32(0)
process.ljmet.excluded_calculators = cms.vstring(
    'DileptonCalc',
    #'WprimeCalc',
    #'LjetsTopoCalcNew'
    )



############################################################
#
# common calculator options
process.load('LJMet.Com.commonCalc_cfi')
process.CommonCalc.dummy_parameter = cms.string('Dummy parameter value')



############################################################
#
# pileup calculator options
process.load('LJMet.Com.pileupCalc_cfi')
process.PileUpCalc.verbosity = process.ljmet.verbosity



############################################################
#
# PDF calculator options
process.load('LJMet.Com.pdfCalc_cfi')



############################################################
#
# Stop calculator options
process.load('LJMet.Com.stopCalc_cfi')



############################################################
#
# Wprime calculator options
process.load('LJMet.Com.wprimeCalc_cfi')
process.WprimeCalc.isWJets = cms.bool(False)



############################################################
#
# LjetsTopoCalcNew options
process.load('LJMet.Com.ljetsTopoCalcNew_cfi')
process.LjetsTopoCalcNew.useBestTop = cms.bool(True)
process.LjetsTopoCalcNew.debug      = cms.bool(True)



#############################################################
#
# Event selector options
process.event_selector = cms.PSet(

    # general settings
    selection = cms.string('StopSelector'),
    isMc = process.ljmet.isMc,
    is_signal_mc = cms.bool(isSignalMc_settable),

    # jet energy scale
    JEC_txtfile              = cms.string(miscFiles['jec_uncertainty']),
    JECup		     = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    do53xJEC                 = cms.bool(True),
    MCL1JetPar               = cms.string(miscFiles['MCL1JetPar']),
    MCL2JetPar               = cms.string(miscFiles['MCL2JetPar']),
    MCL3JetPar               = cms.string(miscFiles['MCL3JetPar']),
    DataL1JetPar             = cms.string(miscFiles['DataL1JetPar']),
    DataL2JetPar             = cms.string(miscFiles['DataL2JetPar']),
    DataL3JetPar             = cms.string(miscFiles['DataL3JetPar']),
    DataResJetPar            = cms.string(miscFiles['DataResJetPar']),

    # b tagging
    btagOP                  = cms.string("CSVM"),
    btag_cond_file    = cms.string(miscFiles['btag_performance']),
    btag_eff_label    = cms.string('TTBARDISCRIMBTAGCSV'),
    btag_sf_label     = cms.string('TTBARWPBTAGCSVM'),
    mistag_label      = cms.string('MISTAGCSVM'),
    btag_discriminant = cms.double(0.679),
    btag_min_discr    = cms.double(0.679),
    btagger           = cms.string('combinedSecondaryVertexBJetTags'),
    BTagUncertUp             = cms.bool(False),
    BTagUncertDown           = cms.bool(False),

    # cuts
    trigger_cut  = cms.bool(False),
    dump_trigger = cms.bool(True),
    trigger_path = cms.string('HLT_IsoMu24_eta2p1_v15'), #sig MC
    pv_cut         = cms.bool(True),
    hbhe_cut       = cms.bool(process.ljmet.isMc!=cms.bool(True)),
    doLaserCalFilt       = cms.bool(process.ljmet.isMc!=cms.bool(True)),
    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(20.0),
    jet_maxeta               = cms.double(5.0),
    min_jet                  = cms.int32(4),
    max_jet                  = cms.int32(4000),
    muon_cuts                = cms.bool(True),
    min_tight_muon           = cms.int32(1),
    tight_muon_minpt         = cms.double(10.0),
    tight_muon_maxeta        = cms.double(2.5),
    tight_muon_mindeltaR_jet = cms.double(0.3),
    #tight_muon_maxIpPv            = cms.double(0.2), # used to be 0.02
    tight_muon_maxMuonPvDeltaZ    = cms.double(0.5),
    max_tight_muon           = cms.int32(1),
    loose_muon_minpt         = cms.double(10.0),
    loose_muon_maxeta        = cms.double(2.5),
    loose_muon_veto          = cms.bool(True),
    electron_minEt           = cms.double(15.0),
    electron_maxeta          = cms.double(2.5),
    electron_veto            = cms.bool(True),
    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(0.0),
    btag_cuts                = cms.bool(True),
    btag_1                   = cms.bool(True),
    btag_2                   = cms.bool(True),
    btag_3                   = cms.bool(False),

    # input collections
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('goodOfflinePrimaryVertices'),
    jet_collection           = cms.InputTag('goodPatJetsPFlow'),
    muon_collection          = cms.InputTag('selectedPatMuonsPFlow'),
    electron_collection      = cms.InputTag('selectedPatElectronsPFlow'),
    met_collection           = cms.InputTag('patMETsPFlow'),
    type1corrmet_collection  = cms.InputTag('pfType1CorrectedMet'),

    )






##################################################################
#
# Input files
#
# NOTE: keep your test inputs in the python files as in
#       this example, and they will be correctly substituted with
#       specified input events when you submit to Condor
#       (
#
#       nEvents and skipEvents are for interactive use, their
#       values will be correctly reset when you submit Condor
#
#input_module = 'LJMet.Com.t1t1bar_WbG_120GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-TLBSM_cff'
#input_module = 'LJMet.Com.t1t1bar_WbG_150GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-TLBSM_cff'
#input_module = 'LJMet.Com.t1t1bar_WbG_172GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-TLBSM_cff'
#input_module = 'LJMet.Com.t1t1bar_WbG_180GeV_8TeV_Rutgers_533_PU2012Startup_Fastsim-TLBSM_cff'
input_module = 'LJMet.Com.TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_sample_cff'
#input_module = 'LJMet.Com.TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1_TLBSM_53x_v2_sample2_cff'
#input_module = 'LJMet.Com.Run2012A_Tlbsm53X_sample_cff'
#input_module = 'LJMet.Com.TT_CT10_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v2_TLBSM_53x_v2_sample_cff'
process.load(input_module)
process.inputs.nEvents    = cms.int32(10000)
process.inputs.skipEvents = cms.int32(0)



############################################################
#
# JSON
JsonFile = miscFiles['json']
myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
if not process.ljmet.isMc:
    process.inputs.lumisToProcess.extend(myList)
        
        
        
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('ljmet_tree'),
    treeName   = cms.string('ljmet'),
    #treeName_el   = cms.string('treetop_el')
)



#######################################################
#
# Object selector options
#



############################################################
#
# Primary vertex
process.load('PhysicsTools.SelectorUtils.pvSelector_cfi')
process.pvSelector.pvSrc   = cms.InputTag('goodOfflinePrimaryVertices')
process.pvSelector.minNdof = cms.double(4.0)
process.pvSelector.maxZ    = cms.double(24.0)
process.pvSelector.maxRho  = cms.double(2.0)



############################################################
#
# tight muon
process.load('PhysicsTools.SelectorUtils.pfMuonSelector_cfi') 
process.pfMuonSelector.version            = cms.string('TOPPAG12_LJETS')
process.pfMuonSelector.Chi2               = cms.double(10.0)
process.pfMuonSelector.minTrackerLayers   = cms.int32(6)
process.pfMuonSelector.minValidMuHits     = cms.int32(1)
process.pfMuonSelector.maxIp              = cms.double(0.2)
process.pfMuonSelector.minPixelHits       = cms.int32(1)
process.pfMuonSelector.minMatchedStations = cms.int32(2)
process.pfMuonSelector.maxPfRelIso        = cms.double(0.12)
process.pfMuonSelector.cutsToIgnore       = cms.vstring()



############################################################
#
# loose muon
process.looseMuonSelector = process.pfMuonSelector.clone()
process.looseMuonSelector.maxPfRelIso  = cms.double(0.2)
process.looseMuonSelector.cutsToIgnore = cms.vstring('Chi2',
                                                     'minTrackerLayers',
                                                     'minValidMuHits',
                                                     'maxIp',
                                                     'minPixelHits',
                                                     'minMatchedStations')



############################################################
#
# electron
process.load('PhysicsTools.SelectorUtils.pfElectronSelector_cfi')
process.pfElectronSelector.version = cms.string('SPRING11')
process.pfElectronSelector.electronIDused = cms.string('eidTight') # data and bg MC
#process.pfElectronSelector.electronIDused = cms.string('eidTightMC') # signal MC
process.pfElectronSelector.PFIso = cms.double(0.2)
process.pfElectronSelector.cutsToIgnore = cms.vstring('electronID',
                                                      'MVA',
                                                      'MaxMissingHits',
                                                      'D0',
                                                      'ConversionRejection')



############################################################
#
# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')
