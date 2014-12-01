/*
  Calculator for the pileup weights

   Author: Zaixing Mao, 2012
*/



#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"


class LjmetFactory;


class PileUpCalc : public BaseCalc{
  
public:
  
  PileUpCalc():mVerbosity(1){}
  virtual ~PileUpCalc(){}

  virtual int BeginJob(){

    // grab parameter values
    if (mPset.exists("verbosity")){
      mVerbosity = mPset.getParameter<int>("verbosity");
    } else mVerbosity = 0;

    // Distribution used for Summer2012 MC. 
    // https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios

    return 0;
  }
    virtual int ProduceEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}

  
private:

    int mVerbosity;

};



static int reg = LjmetFactory::GetInstance()->Register(new PileUpCalc(), "PileUpCalc");



int PileUpCalc::ProduceEvent(edm::EventBase const & event,
                             BaseEventSelector * selector){
    //
    // produce and store some new data for other modules
    //

    double _val = 2.34;
    selector->SetTestValue(_val);

    return 0;
}



int PileUpCalc::AnalyzeEvent(edm::EventBase const & event,
                             BaseEventSelector * selector){

    //
    //_____ Pile-up ____________________________________________
    //
    // get pile up handle
    edm::InputTag puInfoSrc("addPileupInfo");
    edm::Handle<std::vector< PileupSummaryInfo > >  hvPuInfo;

    bool isMc = selector->IsMc();
    float nTrue = -1.;
    int nInteractions = -1;
    int bunchXing = -1;

    if ( isMc ){
        event.getByLabel(puInfoSrc, hvPuInfo);
        for (std::vector<PileupSummaryInfo>::const_iterator iPu=hvPuInfo->begin();
             iPu != hvPuInfo->end();
             ++iPu) {
            if ( iPu->getBunchCrossing() == 0 ){
                bunchXing = iPu->getBunchCrossing();
                nInteractions = iPu->getPU_NumInteractions();
                nTrue = iPu->getTrueNumInteractions();
                //hists["nInteractions"] -> Fill(iPu->getPU_NumInteractions());
                //hists["nTrueInteractions"] -> Fill(iPu->getTrueNumInteractions());
                break;
            }
        }
    }
    SetValue("bunchXing", bunchXing);
    SetValue("nInteractions", nInteractions);
    SetValue("nTrueInteractions", nTrue);

  return 0;
}
