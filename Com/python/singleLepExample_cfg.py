import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("LJMetCom")

relBase    = str('/data1/speer/tblsm/cmssw/git-transition/new-git_5_3_6/')
############################################################
#
# FWLite application options
process.load('LJMet.Com.ljmet_cfi')
process.ljmet.isMc = cms.bool(True)
process.ljmet.excluded_calculators = cms.vstring(
	'DileptonCalc',
	'StopCalc',
	'PdfCalc',
	'ChargedHiggsCalc',
	'TprimeCalc',
	'LjetsTopoCalc',
	'WprimeCalc'
	) 

# common calculator options
process.load('LJMet.Com.commonCalc_cfi')
process.CommonCalc.dummy_parameter = cms.string('Dummy parameter value')

# Stop calculator options
process.load('LJMet.Com.stopCalc_cfi')

# singleLep calculator options
process.load('LJMet.Com.singleLepCalc_cfi')
process.singleLepCalc.isWJets = cms.bool(False)

# LjetsTopoCalc options
process.load('LJMet.Com.ljetsTopoCalcNew_cfi')
process.LjetsTopoCalcNew.useBestTop = cms.bool(True)

# Stop calculator options
process.load('LJMet.Com.JetSubCalc_cfi')


############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('singleLepSelector'),

    # cuts
    debug  = cms.bool(False),

    isMc  = cms.bool(True),
    doLaserCalFilt  = cms.bool(False),
    
    trigger_cut  = cms.bool(True),
    dump_trigger = cms.bool(False),
    
    mctrigger_path_el = cms.string('HLT_Ele27_WP80_v10'), 
    mctrigger_path_mu = cms.string('HLT_IsoMu24_eta2p1_v13'), 
    trigger_path_el = cms.vstring('HLT_Ele27_WP80_v8','HLT_Ele27_WP80_v9','HLT_Ele27_WP80_v10','HLT_Ele27_WP80_v11'), 
    trigger_path_mu = cms.vstring('HLT_IsoMu24_eta2p1_v11','HLT_IsoMu24_eta2p1_v12','HLT_IsoMu24_eta2p1_v13','HLT_IsoMu24_eta2p1_v14','HLT_IsoMu24_eta2p1_v15'),

    #testing muons 
    #mctrigger_path_el = cms.string('HLT_Mu40_v12'), 
    #mctrigger_path_mu = cms.string('HLT_Mu40_v12'), 
    #trigger_path_el = cms.vstring('HLT_Mu40_v12'), 
    #trigger_path_mu = cms.vstring('HLT_Mu40_v12'),
    
    pv_cut         = cms.bool(True),
    hbhe_cut       = cms.bool(True),

    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(4.7),
    min_jet                  = cms.int32(1),
    max_jet                  = cms.int32(4000),
    leading_jet_pt           = cms.double(100.0),
    removeJetLepOverlap		 = cms.bool(True),

    muon_cuts                = cms.bool(True),
    tight_muon_minpt         = cms.double(26.0), # 26 GeV
    tight_muon_maxeta        = cms.double(2.1),
    loose_muon_minpt         = cms.double(10.0),
    loose_muon_maxeta        = cms.double(2.4),
    min_tight_muon           = cms.int32(0),

    electron_cuts                = cms.bool(True),
    tight_electron_minpt         = cms.double(30.0), # 30 GeV
    tight_electron_maxeta        = cms.double(2.5),
    loose_electron_minpt         = cms.double(20.0),
    loose_electron_maxeta        = cms.double(2.5),
    min_tight_electron           = cms.int32(0),
    
    min_tight_lepton         = cms.int32(1),
    max_tight_lepton         = cms.int32(1),
    trigger_consistent       = cms.bool(True),
    second_lepton_veto       = cms.bool(True),
    tau_veto				 = cms.bool(False),
    
    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(20.0),
    
    btag_cuts                = cms.bool(True),
    btagger                  = cms.string('combinedSecondaryVertexBJetTags'),
    btag_min_discr           = cms.double(0.679),
    btag_1                   = cms.bool(False),
    btag_2                   = cms.bool(False),
    btag_3                   = cms.bool(False),

    BTagUncertUp             = cms.bool(False),
    BTagUncertDown           = cms.bool(False),
    JECup                    = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    JEC_txtfile = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V5_DATA_UncertaintySources_AK5PF.txt'),
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('goodOfflinePrimaryVertices'),
    jet_collection           = cms.InputTag('goodPatJetsPFlow'),
    muon_collection          = cms.InputTag('selectedPatMuonsPFlow'),
    electron_collection      = cms.InputTag('selectedPatElectronsPFlow'),
    tau_collection			 = cms.InputTag('selectedPatTausPFlow'),
    met_collection           = cms.InputTag('patMETsPFlow'),
    type1corrmet_collection  = cms.InputTag('pfType1CorrectedMet'),

    do53xJEC                 = cms.bool(False),

    MCL1JetPar               = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_MC_L1FastJet_AK5PFchs.txt'),
    MCL2JetPar               = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_MC_L2Relative_AK5PFchs.txt'),
    MCL3JetPar               = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_MC_L3Absolute_AK5PFchs.txt'),

    DataL1JetPar             = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L1FastJet_AK5PFchs.txt'),
    DataL2JetPar             = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L2Relative_AK5PFchs.txt'),
    DataL3JetPar             = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L3Absolute_AK5PFchs.txt'),
    DataResJetPar            = cms.string('CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L2L3Residual_AK5PFchs.txt')
    )


#######################################################
#
# Input files
#

process.inputs = cms.PSet (
       nEvents    = cms.int32(1000),
           skipEvents = cms.int32(0),
           lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
           fileNames  = cms.vstring(
              ['file:///mnt/hadoop/store/results/B2G/TTJets_SemiLeptMGDecays_8TeV-madgraph/StoreResults-Summer12_DR53X-PU_S10_START53_V7A_ext-v1_TLBSM_53x_v3-99bd99199697666ff01397dad5652e9e/TTJets_SemiLeptMGDecays_8TeV-madgraph/USER/StoreResults-Summer12_DR53X-PU_S10_START53_V7A_ext-v1_TLBSM_53x_v3-99bd99199697666ff01397dad5652e9e/0000/00071EA3-C8E9-E211-85D1-00261894394F.root']
                   )
       )

#process.load(input_module)


# JSON
if (not process.ljmet.isMc==cms.bool(True)):
    JsonFile = 'CMSSW_BASE/src/LJMet/singletPrime/json/Jan222013ReReco_json.txt'
    myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
    process.inputs.lumisToProcess.extend(myList)
        
                
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('ljmet_tree'),
    treeName   = cms.string('ljmet'),
)


#######################################################
#
# Object selector options
#

# Primary vertex
process.load('PhysicsTools.SelectorUtils.pvSelector_cfi')
process.pvSelector.pvSrc   = cms.InputTag('goodOfflinePrimaryVertices')
process.pvSelector.minNdof = cms.double(4.0)
process.pvSelector.maxZ    = cms.double(24.0)
process.pvSelector.maxRho  = cms.double(2.0)

# Tight muon
process.load('LJMet.Com.pfMuonSelector_cfi') 
process.pfMuonSelector.version          = cms.string('TOPPAG12_LJETS')
process.pfMuonSelector.Chi2             = cms.double(10.0)
process.pfMuonSelector.NHits            = cms.int32(0)
process.pfMuonSelector.minValidMuHits       = cms.int32(1)
process.pfMuonSelector.maxIp               = cms.double(0.2)
process.pfMuonSelector.maxPfRelIso            = cms.double(0.12) # 0.12
process.pfMuonSelector.minPixelHits       = cms.int32(1)
process.pfMuonSelector.minMatchedStations = cms.int32(2)
process.pfMuonSelector.minTrackerLayers = cms.int32(6)
process.pfMuonSelector.cutsToIgnore     = cms.vstring('TrackerMuon')

# Loose muon
process.looseMuonSelector = process.pfMuonSelector.clone()
process.looseMuonSelector.PFIso        = cms.double(0.2)
process.looseMuonSelector.nLayersWithMeasurement = cms.int32(0) # not sure why it has to be like this
process.looseMuonSelector.cutsToIgnore = cms.vstring('TrackerMuon',
                                                     'Chi2',
                                                     'NHits',
                                                     'NValMuHits',
                                                     'D0',
                                                     'nPixelHits',
                                                     'nMatchedStations',
                                                     #'nLayersWithMeasurement'
                                                     )

# electron
process.load('LJMet.Com.cutbasedIDSelector_cfi')
process.cutbasedIDSelector.version = cms.string('TIGHT')
process.cutbasedIDSelector.cutsToIgnore     = cms.vstring()

# loose electron
process.looseElectronSelector = process.cutbasedIDSelector.clone()
process.looseElectronSelector.version = cms.string('VETO')
process.looseElectronSelector.cutsToIgnore     = cms.vstring()

# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')
