#ifndef LJMet_Com_interface_LjmetPdfWeightProducer_h
#define LJMet_Com_interface_LjmetPdfWeightProducer_h

/*
  FWLite port of ElectroWeakAnalysis/Utilities/src/PdfWeightProducer

   Author: Gena Kukartsev, 2013
*/



#include <iostream>
#include <vector>

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include <Math/VectorUtil.h>


class LjmetPdfWeightProducer {
  //
  // FWLite port of ElectroWeakAnalysis/Utilities/src/PdfWeightProducer
  //



 public:

  explicit LjmetPdfWeightProducer(const edm::ParameterSet&);
  ~LjmetPdfWeightProducer();
 
  virtual void beginJob() ;
  std::map<std::string,std::vector<double> > produce(edm::EventBase const & event);
  //virtual void endJob() ;
  
  std::string fixPOWHEG_;
  edm::InputTag genTag_;
  edm::InputTag pdfInfoTag_;
  std::vector<std::string> pdfSetNames_;
  std::vector<std::string> pdfShortNames_;


 private:

};



#endif
