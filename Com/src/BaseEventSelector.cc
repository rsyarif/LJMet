/*
  Interface class for FWLite PAT analyzer-selectors
  Specific selectors must implement the () operator

  Author: Gena Kukartsev, 2010,2012
*/



#include "LJMet/Com/interface/BaseEventSelector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "LJMet/Com/interface/FileExists.h"

using namespace std;

BaseEventSelector::BaseEventSelector():
  mName(""),
  mLegend(""){
}



BaseEventSelector::~BaseEventSelector(){
}



void BaseEventSelector::BeginJob(std::map<std::string, edm::ParameterSet const > par){
  std::string _key;
  _key = "event_selector";
  bool _missing_config = false;
  if ( par.find(_key)!=par.end() ){
    if (par[_key].exists("isMc"))        mbPar["isMc"]        = par[_key].getParameter<bool>        ("isMc");
    else                                 mbPar["isMc"]        = false;

    if (par[_key].exists("btagOP"))      msPar["btagOP"]      = par[_key].getParameter<std::string> ("btagOP");
    else                                 msPar["btagOP"]      = "CSVM";

    if (par[_key].exists("JECup"))       mbPar["JECup"]       = par[_key].getParameter<bool>        ("JECup");
    else                                 mbPar["JECup"]       = false;    
    if (par[_key].exists("JECdown"))     mbPar["JECdown"]     = par[_key].getParameter<bool>        ("JECdown");
    else                                 mbPar["JECdown"]     = false;
    if (par[_key].exists("JERup"))       mbPar["JERup"]       = par[_key].getParameter<bool>        ("JERup");
    else                                 mbPar["JERup"]       = false;
    if (par[_key].exists("JERdown"))     mbPar["JERdown"]     = par[_key].getParameter<bool>        ("JERdown");
    else                                 mbPar["JERdown"]     = false;

    if (par[_key].exists("METup"))       mbPar["METup"]       = par[_key].getParameter<bool>        ("METup");
    else                                 mbPar["METup"]       = false;    
    if (par[_key].exists("METdown"))     mbPar["METdown"]     = par[_key].getParameter<bool>        ("METdown");
    else                                 mbPar["METdown"]     = false;
    if (mbPar["METup"] || mbPar["METdown"] ) {
      if (par[_key].exists("METuncert"))   mdPar["METuncert"]     = par[_key].getParameter<double> ("METuncert");
      else                                 mdPar["METuncert"]     = 0.1;
      std::cout << mLegend << "MET uncertainty: modify unclustered energy by "
                << (mbPar["METup"] ? "+" : "-")
		<< mdPar["METuncert"] << std::endl;
    }

    if (par[_key].exists("JEC_txtfile")) msPar["JEC_txtfile"] = par[_key].getParameter<std::string> ("JEC_txtfile");
    else{
      msPar["JEC_txtfile"]  = "";
      _missing_config = true;
    }
    if (par[_key].exists("JEC_source")) msPar["JEC_source"] = par[_key].getParameter<std::string> ("JEC_source");
    else                                msPar["JEC_source"]     = "Total";

    if (par[_key].exists("do53xJEC"))    mbPar["do53xJEC"]    = par[_key].getParameter<bool>        ("do53xJEC");
    else                                 mbPar["do53xJEC"]    = false;

    if (par[_key].exists("BTagUncertUp"))   mbPar["BTagUncertUp"]   = par[_key].getParameter<bool>        ("BTagUncertUp");
    else                                    mbPar["BTagUncertUp"]   = false;
    if (par[_key].exists("BTagUncertDown")) mbPar["BTagUncertDown"] = par[_key].getParameter<bool>        ("BTagUncertDown");
    else                                    mbPar["BTagUncertDown"] = false;

// This is to change the BC and light SFs individually
    if (par[_key].exists("BTagUncertBcUp"))   mbPar["BTagUncertBcUp"]   = par[_key].getParameter<bool>        ("BTagUncertBcUp");
    else                                    mbPar["BTagUncertBcUp"]   = false;
    if (par[_key].exists("BTagUncertBcDown")) mbPar["BTagUncertBcDown"] = par[_key].getParameter<bool>        ("BTagUncertBcDown");
    else                                    mbPar["BTagUncertBcDown"] = false;
    if (par[_key].exists("BTagUncertLightUp"))   mbPar["BTagUncertLightUp"]   = par[_key].getParameter<bool>        ("BTagUncertLightUp");
    else                                    mbPar["BTagUncertLightUp"]   = false;
    if (par[_key].exists("BTagUncertLightDown")) mbPar["BTagUncertLightDown"] = par[_key].getParameter<bool>        ("BTagUncertLightDown");
    else                                    mbPar["BTagUncertLightDown"] = false;

    if (par[_key].exists("MCL1JetPar")) msPar["MCL1JetPar"] = par[_key].getParameter<std::string> ("MCL1JetPar");
    else{
      msPar["MCL1JetPar"] = "../data/START53_V7G_L1FastJet_AK5PFchs.txt";
      _missing_config = true;
    }
    if (par[_key].exists("MCL2JetPar")) msPar["MCL2JetPar"] = par[_key].getParameter<std::string> ("MCL2JetPar");
    else{
      msPar["MCL2JetPar"] = "../data/START53_V7G_L2Relative_AK5PFchs.txt";
      _missing_config = true;
    }
    if (par[_key].exists("MCL3JetPar")) msPar["MCL3JetPar"] = par[_key].getParameter<std::string> ("MCL3JetPar");
    else{
      msPar["MCL3JetPar"] = "../data/START53_V7G_L3Absolute_AK5PFchs.txt";
      _missing_config = true;
    }

    if (par[_key].exists("DataL1JetPar"))  msPar["DataL1JetPar"]  = par[_key].getParameter<std::string> ("DataL1JetPar");
    else{
      msPar["DataL1JetPar"]  = "../data/FT_53_V10_AN3_L1FastJet_AK5PFchs.txt";
      _missing_config = true;
    }
    if (par[_key].exists("DataL2JetPar"))  msPar["DataL2JetPar"]  = par[_key].getParameter<std::string> ("DataL2JetPar");
    else{
      msPar["DataL2JetPar"]  = "../data/FT_53_V10_AN3_L2Relative_AK5PFchs.txt";
      _missing_config = true;
    }
    if (par[_key].exists("DataL3JetPar"))  msPar["DataL3JetPar"]  = par[_key].getParameter<std::string> ("DataL3JetPar");
    else{
      msPar["DataL3JetPar"]  = "../data/FT_53_V10_AN3_L3Absolute_AK5PFchs.txt";
      _missing_config = true;
    }
    if (par[_key].exists("DataResJetPar")) msPar["DataResJetPar"] = par[_key].getParameter<std::string> ("DataResJetPar");
    else{
      msPar["DataResJetPar"] = "../data/FT_53_V10_AN3_L2L3Residual_AK5PFchs.txt";
      _missing_config = true;
    }

    if (_missing_config){
      std::cout << mLegend
		<< "ONE OF THE FOLLOWING CONFIG OPTIONS MISSING!\n"
		<< "MCL1JetPar, MCL2JetPar, MCL3JetPar, DataL1JetPar, DataL2JetPar,\n"
		<< "DataL3JetPar, DataResJetPar" <<std::endl;
      std::cout << mLegend
		<< "USING DEFAULT VALUES" << std::endl;
    }

  }

  
  msPar["btagger"] = mBtagCond.getAlgoName(msPar["btagOP"]);
  mdPar["btag_min_discr"] = mBtagCond.getDiscriminant(msPar["btagOP"]);

  bTagCut = mdPar["btag_min_discr"];
  std::cout << "b-tag check "<<msPar["btagOP"]<<" "<< msPar["btagger"]<<" "<<mdPar["btag_min_discr"]<<std::endl;

  if ( mbPar["isMc"] && ( mbPar["JECup"] || mbPar["JECdown"])) {
    std::cout << mLegend << "Applying 53X jet energy corrections Uncertainty: "<< (mbPar["JECup"]?"Up\n":"Down\n");
    std::cout << mLegend << "   Source : "<< msPar["JEC_source"]<<std::endl;
    std::cout << mLegend << "   File   : "<< msPar["JEC_txtfile"]<<std::endl;
    fexists(msPar["JEC_txtfile"], true);
    jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters(msPar["JEC_txtfile"].c_str(), msPar["JEC_source"])));
  }

  //gSystem->Load("libFWCoreFWLite.so");
  //AutoLibraryLoader::enable();
 
  vector<JetCorrectorParameters> vPar;
  if ( mbPar["isMc"] && ( mbPar["do53xJEC"] ) ) {
    // Create the JetCorrectorParameter objects, the order does not matter.
    // START53_V7G
    fexists(msPar["MCL3JetPar"], true);
    fexists(msPar["MCL2JetPar"], true);
    fexists(msPar["MCL1JetPar"], true);
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(msPar["MCL3JetPar"]);
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(msPar["MCL2JetPar"]);
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(msPar["MCL1JetPar"]);
    // Load the JetCorrectorParameter objects into a vector,
    // IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);

    std::cout << mLegend << "Applying 53X jet energy corrections" << std::endl;

  }
  else if ( !mbPar["isMc"] && ( mbPar["do53xJEC"] ) ) {
    // Create the JetCorrectorParameter objects, the order does not matter.
    // GR_P_V43
    fexists(msPar["DataResJetPar"], true);
    fexists(msPar["DataL3JetPar"], true);
    fexists(msPar["DataL2JetPar"], true);
    fexists(msPar["DataL1JetPar"], true);
    JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(msPar["DataResJetPar"]); 
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(msPar["DataL3JetPar"]);
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(msPar["DataL2JetPar"]);
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(msPar["DataL1JetPar"]);
    // Load the JetCorrectorParameter objects into a vector,
    // IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);
    vPar.push_back(*ResJetPar);

    std::cout << mLegend << "Applying 53X jet energy corrections" << std::endl;

  }
  else{

    std::cout << mLegend 
	      << "NOT applying 53X jet energy corrections - ARE YOU SURE?" 
	      << std::endl;
    
  }
  JetCorrector = new FactorizedJetCorrector(vPar);


}



void BaseEventSelector::EndJob(){
}



void BaseEventSelector::AnalyzeEvent( edm::EventBase const & event, LjmetEventContent & ec ){
}



std::string BaseEventSelector::GetName(){
  return mName;
}



void BaseEventSelector::setName(std::string name){
  mName = name;
}



void BaseEventSelector::init(){
  //
  // private init method to be called by LjmetFactory
  // when registering the selector
  //

  mLegend = "["+mName+"]: ";
  std::cout << mLegend << "registering "
	    << mName << std::endl;

  return;
}



void BaseEventSelector::SetHistogram(std::string name,
				     int nbins, 
				     double low, 
				     double high){
  //
  // Declare a new histogram to be created for the module
  //

  mpEc->SetHistogram(mName, name, nbins, low, high);

  return;
}



void BaseEventSelector::SetHistValue(std::string name, 
				     double value){
  mpEc->SetHistValue(mName, name, value);

  return;
}



// evaluates a signed perp components of v1 relative to v2
// the sign is defined by Phi
double BaseEventSelector::GetPerp(TVector3 & v1, TVector3 & v2){
  
  double perp;
  
  double _mag = v1.Cross(v2.Unit()).Mag();
  double _phi1 = v1.Phi();
  double _phi2 = v2.Phi();
  double _dphi = _phi1 - _phi2;
  double rPI = TMath::Pi();
  if (_dphi>rPI || (_dphi>-rPI && _dphi<0)) perp = _mag;
  else perp = -_mag;
  
  return perp;
}

 

bool BaseEventSelector::AreMatched ( const reco::Candidate & c1,
                                     const reco::Candidate & c2,
                                     double DR,
                                     double DPtRel ) {
  unsigned int nPass=0;
  if (  (reco::deltaR(c1, c2) < DR)   && (fabs(c2.pt() - c1.pt())/ c2.pt()<DPtRel)) {
    nPass++;
  }
  return (nPass>0);
}



std::vector<edm::Ptr<pat::Jet> >      const & BaseEventSelector::GetAllJets()           const { return mvAllJets; }
std::vector<edm::Ptr<pat::Jet> >      const & BaseEventSelector::GetSelectedJets()      const { return mvSelJets; }
std::vector<edm::Ptr<pat::Jet> >      const & BaseEventSelector::GetLooseJets()         const { return mvSelJets; }
std::vector<edm::Ptr<pat::Jet> >      const & BaseEventSelector::GetSelectedBtagJets()  const { return mvSelBtagJets; }
std::vector<std::pair<TLorentzVector,bool>>  const & BaseEventSelector::GetCorrJetsWithBTags()  const { return mvCorrJetsWithBTags; }
std::vector<edm::Ptr<pat::Muon> >     const & BaseEventSelector::GetAllMuons()          const { return mvAllMuons; }
std::vector<edm::Ptr<pat::Muon> >     const & BaseEventSelector::GetSelectedMuons()     const { return mvSelMuons; }
std::vector<edm::Ptr<pat::Muon> >     const & BaseEventSelector::GetLooseMuons()        const { return mvLooseMuons; }
std::vector<edm::Ptr<pat::Electron> > const & BaseEventSelector::GetAllElectrons()      const { return mvAllElectrons; }
std::vector<edm::Ptr<pat::Electron> > const & BaseEventSelector::GetSelectedElectrons() const { return mvSelElectrons; }
std::vector<edm::Ptr<pat::Electron> > const & BaseEventSelector::GetLooseElectrons()    const { return mvLooseElectrons; }
edm::Ptr<pat::MET>                    const & BaseEventSelector::GetMet()               const { return mpMet; }
edm::Ptr<reco::PFMET>                    const & BaseEventSelector::GetType1CorrMet()   const { return mpType1CorrMet; }
TLorentzVector                        const & BaseEventSelector::GetCorrectedMet()      const { return correctedMET_p4; }
std::vector<unsigned >                const & BaseEventSelector::GetSelectedTriggers()  const { return mvSelTriggers; }
std::vector<edm::Ptr<reco::Vertex> >  const & BaseEventSelector::GetSelectedPVs()       const { return mvSelPVs; }

double const & 
BaseEventSelector::GetTestValue() const{
  return mTestValue; 
}



void BaseEventSelector::SetTestValue(double & test){
  mTestValue = test;
  return;
}



void BaseEventSelector::SetCorrJetsWithBTags(std::vector<std::pair<TLorentzVector,bool>> & jets){
  mvCorrJetsWithBTags = jets; 
  return;
}



void BaseEventSelector::SetCorrectedMet(TLorentzVector & met){
  correctedMET_p4 = met;
  return;
}



void BaseEventSelector::SetMc(bool isMc){
  mbIsMc = isMc;
  return;
}



bool BaseEventSelector::IsMc(){
  return mbIsMc;
}

bool BaseEventSelector::isJetTagged(const pat::Jet & jet, edm::EventBase const & event, bool applySF)
{
  bool _isTagged = false;


  if ( jet.bDiscriminator( msPar["btagger"] ) > bTagCut ) _isTagged = true;
  // DEBUG: list all available taggers with values
  //     cout << "_isTagged "<<_isTagged<<endl;

  //       const std::vector<std::pair<std::string, float> > & vpd = jet.getPairDiscri();
  //       for ( std::vector<std::pair<std::string, float> >::const_iterator i = vpd.begin();
  // 	    i != vpd.end(); ++i){
  // 	std::cout << mLegend << i->first << ", " << i->second << std::endl;
  //       }

  if (mbPar["isMc"] && applySF){
    
    TLorentzVector lvjet = correctJet(jet, event);

    double _lightSf  = mBtagCond.GetMistagScaleFactor(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
    double _lightEff = mBtagCond.GetMistagRateMC(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*_lightSf;
    if ( mbPar["BTagUncertUp"] || mbPar["BTagUncertLightUp"]) _lightSf += mBtagCond.GetMistagSFUncertUp(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
    else if ( mbPar["BTagUncertDown"] || mbPar["BTagUncertLightDown"])_lightSf -= mBtagCond.GetMistagSFUncertDown(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);

    int _jetFlavor = abs(jet.partonFlavour());
    double _btagSf  = mBtagCond.GetBtagScaleFactor(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
    double _btagEff = mBtagCond.GetBtagEfficiencyMC(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*_btagSf;
    if ( mbPar["BTagUncertUp"] || mbPar["BTagUncertBcUp"]) _btagSf += (mBtagCond.GetBtagSFUncertUp(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*(_jetFlavor==4?2:1));
    else if ( mbPar["BTagUncertDown"] || mbPar["BTagUncertBcDown"])_btagSf -= (mBtagCond.GetBtagSFUncertDown(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*(_jetFlavor==4?2:1));

    mBtagSfUtil.SetSeed(abs(static_cast<int>(sin(jet.phi())*100000)));

    // sanity check
    bool _orig_tag = _isTagged;

    mBtagSfUtil.modifyBTagsWithSF(_isTagged, _jetFlavor, _btagSf, _btagEff,
				  _lightSf, _lightEff);

    // sanity check
    if (_isTagged != _orig_tag) ++mNBtagSfCorrJets;

  } // end of btag scale factor corrections
  return _isTagged;

}

TLorentzVector BaseEventSelector::correctJet(const pat::Jet & jet, edm::EventBase const & event)
{
  // JES and JES systematics
    pat::Jet correctedJet;
    if (mbPar["do53xJEC"])
        correctedJet = jet.correctedJet(0);                 //copy original jet
   else
        correctedJet = jet;                                 //copy 52x corrected jet

  double ptscale = 1.0;
  double unc = 1.0;
  double pt = correctedJet.pt();
  double correction = 1.0;

  edm::Handle<double> rhoHandle;
  edm::InputTag rhoSrc_("kt6PFJets", "rho");
  event.getByLabel(rhoSrc_, rhoHandle);
  double rho = std::max(*(rhoHandle.product()), 0.0);

  if ( mbPar["isMc"] ){ 

    if (mbPar["do53xJEC"]) {
      // 53x Jet Energy corrections were not applied to the TLBSM 53x v2 pat-tuples
      // Therefore, we need to undo the 52x corrections and then apply the 53x ones

      double pt_raw = jet.correctedJet(0).pt();
      JetCorrector->setJetEta(jet.eta());
      JetCorrector->setJetPt(pt_raw);
      JetCorrector->setJetA(jet.jetArea());
      JetCorrector->setRho(rho); 

      try{
	correction = JetCorrector->getCorrection();
      }
      catch(...){
	std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
	std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
	std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
      }
      
      correctedJet.scaleEnergy(correction);
      pt = correctedJet.pt();

    }
    double factor = 0.0; //For Nominal Case

    if ( abs(jet.eta()) < 0.5 ) {
      factor = .052;
      if (mbPar["JERup"]) factor = 0.115;
      if (mbPar["JERdown"]) factor = -0.01;
    }
    else if ( abs(jet.eta()) < 1.1 && abs(jet.eta()) >= 0.5 ) {
      factor = 0.057;
      if (mbPar["JERup"]) factor = 0.114;
      if (mbPar["JERdown"]) factor = 0.001;
    }
    else if ( abs(jet.eta()) < 1.7 && abs(jet.eta()) >= 1.1 ) {
      factor = 0.096;
      if (mbPar["JERup"]) factor = 0.161;
      if (mbPar["JERdown"]) factor = 0.032;
    }
    else if ( abs(jet.eta()) < 2.3 && abs(jet.eta()) >= 1.7 ) {
      factor = 0.134;
      if (mbPar["JERup"]) factor = 0.228;
      if (mbPar["JERdown"]) factor = 0.042;

    }
    else if (abs(jet.eta()) < 5.0 && abs(jet.eta()) >=2.3 ) {
      factor = 0.288;
      if (mbPar["JERup"]) factor = 0.488;
      if (mbPar["JERdown"]) factor = 0.089;
    }

    const reco::GenJet * genJet = jet.genJet();
    if (genJet && genJet->pt()>15. && (abs(genJet->pt()/pt-1)<0.5)){
      double gen_pt = genJet->pt();
      double reco_pt = pt;
      double deltapt = (reco_pt - gen_pt) * factor;
      ptscale = max(0.0, (reco_pt + deltapt) / reco_pt);
    }

    if ( mbPar["JECup"] || mbPar["JECdown"]) {
      jecUnc->setJetEta(jet.eta());
      jecUnc->setJetPt(pt*ptscale);

      if (mbPar["JECup"]) { 
	try{
          unc = jecUnc->getUncertainty(true);
	}
	catch(...){ // catch all exceptions. Jet Uncertainty tool throws when binning out of range
	  std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
	  std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
	  std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
	  unc = 0.0;
	}
        unc = 1 + unc; 
      }
      else { 
	try{
          unc = jecUnc->getUncertainty(false);
	}
	catch(...){
	  std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
	  std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
	  std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
	  unc = 0.0;
	}
        unc = 1 - unc; 
      }

    }
  }
  else if (!mbPar["isMc"]) {
      
      if (mbPar["do53xJEC"]) {
	
          // 53x Jet Energy corrections were not applied to the TLBSM 53x v2 pat-tuples
          // Therefore, we need to undo the 52x corrections and then apply the 53x ones

	  double pt_raw = jet.correctedJet(0).pt();
          JetCorrector->setJetEta(jet.eta());
          JetCorrector->setJetPt(pt_raw);
          JetCorrector->setJetA(jet.jetArea());
          JetCorrector->setRho(rho); 
	
	  try{
	    correction = JetCorrector->getCorrection();
	  }
	  catch(...){
	    std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
	    std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
	    std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
	  }
	  
          correctedJet.scaleEnergy(correction);
          pt = correctedJet.pt();
      }
  }

  TLorentzVector jetP4;
  jetP4.SetPtEtaPhiM(correctedJet.pt()*unc*ptscale, correctedJet.eta(),correctedJet.phi(), correctedJet.mass() );
//   if (correctedJet.pt()> 15.)
//   {
//   cout << "correction: "<< unc<<" "<<ptscale<<endl;
//   std::cout<<"jet pt before : "<<correctedJet.pt()<<" eta: "<< correctedJet.eta()<<" phi: "<<correctedJet.phi()<<" mass: "<<correctedJet.mass()<<" energy: "<<correctedJet.energy()<<std::endl;
//   std::cout<<"jet pt after  : "<<jetP4.Pt()<<" eta: "<<jetP4.Eta()<<" phi: "<<jetP4.Phi()<<" mass: "<<jetP4.M()<<" energy: "<<jetP4.E()<<std::endl;
// }

  // sanity check - save correction of the first jet
  if (mNCorrJets==0){
    double _orig_pt = jet.correctedJet(0).pt();
    if (fabs(_orig_pt)<0.000000001){
      _orig_pt = 0.000000001;
    }
    SetHistValue("jes_correction", jetP4.Pt()/_orig_pt);
    ++mNCorrJets;
  }


  return jetP4;
}


TLorentzVector BaseEventSelector::correctMet(const pat::MET & met, edm::EventBase const & event)
{
  double correctedMET_px = met.px();
  double correctedMET_py = met.py();
// cout <<"MET correction "<<correctedMET_px <<" "<<correctedMET_py<<endl;
    for (std::vector<edm::Ptr<pat::Jet> >::const_iterator ijet = mvAllJets.begin();
         ijet != mvAllJets.end(); ++ijet){
      TLorentzVector lv = correctJet(**ijet, event);
      correctedMET_px +=  (**ijet).px() -lv.Px();
      correctedMET_py += 	(**ijet).py() -lv.Py();
// cout << endl<<correctedMET_px <<" "<<correctedMET_py<<" "<<
// (**ijet).px() <<" "<<lv.Px()<<" "<< (**ijet).py() <<" "<<lv.Py()<<" "<< (**ijet).pt() <<" "<<lv.Perp()<<endl<<endl;
    }

// cout <<"MET correction "<<correctedMET_px <<" "<<correctedMET_py<<endl;

// This is for the MET uncertainty
  if (mbPar["METup"] || mbPar["METdown"] ) {
    TVector2 uMET = unclusMET(met);
    int sign = (mbPar["METup"] ? +1 : -1);
    correctedMET_px += uMET.Px()*sign*mdPar["METuncert"];
    correctedMET_py += uMET.Py()*sign*mdPar["METuncert"];
  }
// cout <<"MET correction "<<correctedMET_px <<" "<<correctedMET_py<<endl;

  correctedMET_p4.SetPxPyPzE(correctedMET_px, correctedMET_py, 0, sqrt(correctedMET_px*correctedMET_px+correctedMET_py*correctedMET_py));
// cout << "MET before : "<<met.px()<<" "<<met.py()<<" "<<met.pt()<<endl;
// cout << "MET after  : "<<correctedMET_p4.Px()<<" "<<correctedMET_p4.Py()<<" "<<correctedMET_p4.Perp()<<endl;

  // sanity check histogram
  double _orig_met = met.pt();
  if (fabs(_orig_met)<0.000000001){
    _orig_met = 0.000000001;
  }
  SetHistValue("met_correction", correctedMET_p4.Pt()/_orig_met);


  return correctedMET_p4;
}

TVector2 BaseEventSelector::unclusMET(const pat::MET & met) const{
  double missetX = met.px();
  double missetY = met.py();
  // cout <<"Now GetUnclusScaledMET\n";
  // cout << "MET original   : "<< met.pt()<<  " "<< met.px()<<  " "<< met.py()<<  " "<<endl;
  // cout << "Jets used to calculate unclustered MET\n"<<endl;
  for (std::vector<edm::Ptr<pat::Jet> >::const_iterator ijet = mvAllJets.begin();
       ijet != mvAllJets.end(); ++ijet){
    if ((fabs((**ijet).eta())< 5.0) && ((**ijet).pt()>10.0)) {
    // cout<<" jet eta,pt "<<(**ijet).eta()<<" "<<(**ijet).pt()<<" "<<(*ijet)->px()<<" "<<(*ijet)->py()<<endl;
    missetX += (**ijet).px();
    missetY += (**ijet).py();
    }
  }
  // cout << "MET after jet  : "<<sqrt(pow(missetX,2.) + pow(missetY,2.))<<  " "<< missetX<<  " "<< missetY<<  " "<<endl;

  for (std::vector<edm::Ptr<pat::Electron> >::const_iterator iel = mvAllElectrons.begin(); iel != mvAllElectrons.end(); iel++){   
  // cout<<" e eta,pt "<<(*iel)->eta()<<" "<<(*iel)->pt()<<" "<<(*iel)->px()<<" "<<(*iel)->py()<<endl;
    missetX += (*iel)->px();
    missetY += (*iel)->py();
  }
  // cout << "MET after el   : "<<sqrt(pow(missetX,2.) + pow(missetY,2.))<<  " "<< missetX<<  " "<< missetY<<  " "<<endl;
  for (std::vector<edm::Ptr<pat::Muon> >::const_iterator imu = mvAllMuons.begin(); imu != mvAllMuons.end(); imu++){
    // cout<<" mu eta,pt "<<(*imu)->eta()<<" "<<(*imu)->pt()<<" "<<(*imu)->px()<<" "<<(*imu)->py()<<endl;
    missetX += (*imu)->px();
    missetY += (*imu)->py();
  }
  // cout << "MET after mu   : "<<sqrt(pow(missetX,2.) + pow(missetY,2.))<<  " "<< missetX<<  " "<< missetY<<  " "<<endl;
//     for (unsigned int i=0; i<taus.size(); i++){
//       cout<<" tau eta,pt "<<taus[i].p4.Eta()<<" "<<taus[i].p4.Pt()<<endl;
//       missetX += taus[i].p4.Px();
//       missetY += taus[i].p4.Py();
//     }
    return TVector2(missetX,missetY);

}




void BaseEventSelector::SetEventContent(LjmetEventContent * pEc){

  mpEc = pEc;

  return;
}



void BaseEventSelector::Init( void ){

  // init sanity check histograms
  mpEc->SetHistogram(mName, "jes_correction",100,0.8,1.2);
  mpEc->SetHistogram(mName, "met_correction",100,0.0,2.0);
  mpEc->SetHistogram(mName, "nBtagSfCorrections", 100,0,10);

  return;
}



void BaseEventSelector::BeginEvent( edm::EventBase const & event,
				    LjmetEventContent & ec ){
  //
  // Do what any event selector must do before event gets checked
  //

  mNCorrJets = 0;
  mNBtagSfCorrJets = 0;

  return;
}



void BaseEventSelector::EndEvent( edm::EventBase const & event,
				    LjmetEventContent & ec ){
  //
  // Do what any event selector must do after event processing is done
  // (but before event content gets saved to file)
  //

  SetHistValue("nBtagSfCorrections", mNBtagSfCorrJets);

  return;
}
