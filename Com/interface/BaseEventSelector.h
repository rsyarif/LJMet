#ifndef LJMet_Com_interface_BaseEventSelector_h
#define LJMet_Com_interface_BaseEventSelector_h

/*
   Interface class for FWLite PAT analyzer-selectors
   Specific selectors must implement the () operator

   Author: Gena Kukartsev, 2010,2012
*/

#include <cmath>
#include <iostream>

#include "TROOT.h"
#include "TVector3.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h"
#include "LJMet/Com/interface/BTagSFUtil.h"
#include "LJMet/Com/interface/BtagHardcodedConditions.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "FWCore/Framework/interface/Event.h"

class BaseEventSelector : public EventSelector {
  //
  // Base class for all event selector plugins
  //


  friend class LjmetFactory;


 public:
  

  BaseEventSelector();
    
    
  virtual ~BaseEventSelector();
  

  virtual void BeginJob(std::map<std::string, edm::ParameterSet const > par);
  virtual bool operator()( edm::EventBase const & event, pat::strbitset & ret) = 0;
  virtual void EndJob();



  virtual void AnalyzeEvent( edm::EventBase const & event, LjmetEventContent & ec );
  std::string  GetName();
  double       GetPerp(TVector3 & v1, TVector3 & v2);
  bool         AreMatched ( const reco::Candidate & c1,
			    const reco::Candidate & c2,
			    double DR,
			    double DPtRel );
  
  
  std::vector<edm::Ptr<pat::Jet> >            const & GetAllJets()           const;
  std::vector<edm::Ptr<pat::Jet> >            const & GetSelectedJets()      const;
  std::vector<edm::Ptr<pat::Jet> >            const & GetLooseJets()         const;
  std::vector<edm::Ptr<pat::Jet> >            const & GetSelectedBtagJets()  const;
  std::vector<std::pair<TLorentzVector,bool>> const & GetCorrJetsWithBTags()  const;
  std::vector<edm::Ptr<pat::Muon> >           const & GetAllMuons()          const;
  std::vector<edm::Ptr<pat::Muon> >           const & GetSelectedMuons()     const;
  std::vector<edm::Ptr<pat::Muon> >           const & GetLooseMuons()        const;
  std::vector<edm::Ptr<pat::Electron> >       const & GetAllElectrons()      const;
  std::vector<edm::Ptr<pat::Electron> >       const & GetSelectedElectrons() const;
  std::vector<edm::Ptr<pat::Electron> >       const & GetLooseElectrons()    const;
  edm::Ptr<pat::MET>                          const & GetMet()               const;
  edm::Ptr<reco::PFMET>                       const & GetType1CorrMet()      const;
  TLorentzVector                              const & GetCorrectedMet()      const;
  std::vector<unsigned int>                   const & GetSelectedTriggers()  const;
  std::vector<edm::Ptr<reco::Vertex> >        const & GetSelectedPVs()  const;
  double const & GetTestValue() const;


  void SetMc(bool isMc);
  bool IsMc();


  // LJMET event content setters
  void Init( void );
  void SetEventContent(LjmetEventContent * pEc);
  void SetHistogram(std::string name, int nbins, double low, double high);
  void SetHistValue(std::string name, double value);
  void SetTestValue(double & test);

  void SetCorrectedMet(TLorentzVector & met);
  void SetCorrJetsWithBTags(std::vector<std::pair<TLorentzVector,bool>> & jets);

  bool isJetTagged(const pat::Jet &jet, edm::EventBase const & event, bool applySF = true);
  TLorentzVector correctJet(const pat::Jet & jet, edm::EventBase const & event);
  TLorentzVector correctMet(const pat::MET & met, edm::EventBase const & event);

 protected:

  std::vector<edm::Ptr<pat::Jet> >      mvAllJets;
  std::vector<edm::Ptr<pat::Jet> >      mvSelJets;
  std::vector<edm::Ptr<pat::Jet> >      mvLooseJets;
  std::vector<std::pair<TLorentzVector,bool> >      mvCorrJetsWithBTags;
  std::vector<edm::Ptr<pat::Jet> >      mvSelBtagJets;
  std::vector<edm::Ptr<pat::Muon> >     mvAllMuons;
  std::vector<edm::Ptr<pat::Muon> >     mvSelMuons;
  std::vector<edm::Ptr<pat::Muon> >     mvLooseMuons;
  std::vector<edm::Ptr<pat::Electron> > mvAllElectrons;
  std::vector<edm::Ptr<pat::Electron> > mvSelElectrons;
  std::vector<edm::Ptr<pat::Electron> > mvLooseElectrons;
  edm::Ptr<pat::MET>                    mpMet;
  edm::Ptr<reco::PFMET>                 mpType1CorrMet;
  TLorentzVector                        correctedMET_p4;
  std::vector<unsigned int>             mvSelTriggers;
  std::vector<edm::Ptr<reco::Vertex> >  mvSelPVs;
  double                                mTestValue;

  // containers for config parameter values
  std::map<std::string,bool>           mbPar;
  std::map<std::string,int>            miPar;
  std::map<std::string,double>         mdPar;
  std::map<std::string,std::string>    msPar;
  std::map<std::string, edm::InputTag> mtPar;
  std::map<std::string,std::vector<std::string> > mvsPar;

  std::string mName;
  std::string mLegend;

  bool mbIsMc;
  


 private:

  void init();
  void setName(std::string name);
  void BeginEvent( edm::EventBase const & event, LjmetEventContent & ec );
  void EndEvent( edm::EventBase const & event, LjmetEventContent & ec );

  TVector2 unclusMET(const pat::MET & met) const;


  BTagSFUtil mBtagSfUtil;
  BtagHardcodedConditions mBtagCond;
  double bTagCut;
  JetCorrectionUncertainty *jecUnc;
  FactorizedJetCorrector *JetCorrector;

  LjmetEventContent * mpEc;

  int mNCorrJets;
  int mNBtagSfCorrJets;
};



#endif
