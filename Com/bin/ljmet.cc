
//
// FWLite PAT analyzer 
//
// Gena Kukartsev, March 2012  
// 

#include <cmath>
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "LJMet/Com/interface/BaseEventSelector.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "Math/GenVector/Cartesian2D.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"
#include "LJMet/Com/interface/HcalLaserEventFilter2012Standalone.h"
#include "LJMet/Com/interface/FileExists.h"



//===============================================================>
//
// forward declarations
//



bool JsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        const edm::EventBase &event);



///////////////////////////
// ///////////////////// //
// // Main Subroutine // //
// ///////////////////// //
///////////////////////////

int main (int argc, char* argv[]) 
{

  // legend for self ID in messages
  std::string legend = "[";
  legend.append(argv[0]);
  legend.append("]: ");


  // processing the config file
  std::cout << legend << "getting parameters from config file" << std::endl;
  std::cout << legend << "the following parameter sets found:" << std::endl;
  // Get the python configuration
  PythonProcessDesc builder(argv[1]);
  boost::shared_ptr<edm::ProcessDesc> b = builder.processDesc();
  boost::shared_ptr<edm::ParameterSet> parameters = b->getProcessPSet();
  parameters->registerIt();


  // get some config parameters
  // other parameters are accessed directly in the event selector
  edm::ParameterSet const& selectorParams = parameters->getParameter<edm::ParameterSet>("event_selector");
  edm::ParameterSet const& ljmetParams      = parameters->getParameter<edm::ParameterSet>("ljmet");
  edm::ParameterSet const& inputs         = parameters->getParameter<edm::ParameterSet>("inputs");
  edm::ParameterSet const& outputs        = parameters->getParameter<edm::ParameterSet>("outputs");
  int const maxEvents = inputs.getParameter<int>("nEvents");
  int const nEventsToSkip = inputs.getParameter<int>("skipEvents");


  // output file name base
  std::string _outputName = outputs.getParameter<std::string>("outputName");


  // log file
  std::string _logName = _outputName+".log";
  fstream _logfile;
  _logfile.open(_logName, fstream::out);


  // usage
  if ( argc < 2 ) {
    std::cout << legend << "usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }


  // greeting message
  std::cout << legend << "hello from " << argv[0] << "!" << std::endl;


  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();


  // harvest all parameter sets from config file
  std::map<std::string, edm::ParameterSet const> mPar;
  std::vector<std::string> vParNames = parameters->getParameterNames();
  for (std::vector<std::string>::const_iterator name = vParNames.begin();
       name != vParNames.end();
       ++name){

    if ( name->find("@")==std::string::npos ){
      // skip weird entries starting with @
      std::cout << legend << *name << std::endl;
      edm::ParameterSet const _ps = parameters->getParameter<edm::ParameterSet>(*name);
      mPar.insert( std::pair<std::string, edm::ParameterSet const>(*name, _ps) );
    }
  }
  //std::cout << std::endl;


  // TFileService for saving the output ROOT tree
  std::cout << legend << "setting up TFileService" << std::endl;
  fwlite::TFileService fs = fwlite::TFileService( _outputName+".root" );

  
  // output tree
  std::cout << legend << "Creating output tree" << std::endl;
  std::string const _treename = outputs.getParameter<std::string>("treeName");
  TTree * _tree = fs.make<TTree>(_treename.c_str(), _treename.c_str(), 64000000);
    

  // book histograms
  TFileDirectory theDir = fs.mkdir( "histos" );
  std::cout << legend << "booking histograms" << std::endl;
  std::map<std::string, TH1*> hists;
  hists["nevents"] = theDir.make<TH1I>( "nevents",
					"nevents", 
					1, 0, 2 ) ;
//  hists["nInteractions"] = theDir.make<TH1D>( "hist_nInteractions",
//					      "nInteractions", 
//					      400, 0, 400 ) ;
//  hists["nTrueInteractions"] = theDir.make<TH1D>( "hist_nTrueInteractions",
//						  "nTrueInteractions", 
//						  400, 0, 400 ) ;

  // internal LJMet event content
  LjmetEventContent ec(mPar);
  ec.SetTree(_tree);

  // maps for branches
  // each entry will create and fill a branch
  //std::map<string,int> m_int_branch;
  //std::map<string,double> m_double_branch;
  //std::map<string,std::vector<double> > m_vector_branch;


  // The factory for event selector and calculator plugins
  LjmetFactory * factory = LjmetFactory::GetInstance();


  // choose event selector
  std::cout << legend << "instantiating the event selector" << std::endl;
  std::string selection = selectorParams.getParameter<std::string>("selection");
  BaseEventSelector * theSelector = 0;
  theSelector = factory->GetEventSelector(selection);

  // sanity check histograms from the selector
  theSelector->SetEventContent(&ec);
  theSelector->Init();

  theSelector->BeginJob(mPar);


  // set excluded calculators
  if (mPar.find("ljmet")!=mPar.end()){

    std::vector<std::string> vLjmetParams = mPar["ljmet"].getParameterNames();

    if (std::find(vLjmetParams.begin(), vLjmetParams.end(), "excluded_calculators")!=vLjmetParams.end()){
      std::vector<std::string> vExcl = 
	mPar["ljmet"].getParameter<std::vector<std::string> >("excluded_calculators");

      factory->SetExcludedCalcs(vExcl);
    }

  }


  // send config parameters to calculators
  factory->SetAllCalcConfig(mPar);


  
  // Run BeginJob() for calculators
  factory->BeginJobAllCalc();



  // create histograms
  std::map<std::string,std::map<std::string,LjmetEventContent::HistMetadata> > & mh = ec.GetHistMap();
  std::map<std::string,std::map<std::string,LjmetEventContent::HistMetadata> >::iterator iMod;
  std::map<std::string,LjmetEventContent::HistMetadata>::iterator iHist;
  for (iMod=mh.begin();iMod!=mh.end();++iMod){

    TFileDirectory _dir = fs.mkdir( iMod->first.c_str() );

    for (iHist=iMod->second.begin();iHist!=iMod->second.end();++iHist){
      std::cout << legend
		<< "Creating " << iMod->first << "/"
		<< iHist->second.GetName() << std::endl;
      iHist->second.SetHist( _dir.make<TH1F>(iHist->second.GetName().c_str(),
					     iHist->second.GetName().c_str(),
					     iHist->second.GetNBins(),
					     iHist->second.GetXMin(),
					     iHist->second.GetXMax() ) );
    }
  }



  
  // data/MC flag
  // FIXME: may be superseded by auto data/MC detection
  bool isMc    = ljmetParams.getParameter<bool>("isMc");
  theSelector->SetMc(isMc);
  if (isMc) std::cout << legend << "Using MC truth info" << std::endl;

  //=============================================================>
  //
  // JSON file processing
  //
  std::vector<edm::LuminosityBlockRange> vJson;
  if ( (!isMc) && (inputs.exists("lumisToProcess")) ) {
      std::vector<edm::LuminosityBlockRange> const & lumisTemp =
	inputs.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> > ("lumisToProcess");
      vJson.resize( lumisTemp.size() );
      copy( lumisTemp.begin(), lumisTemp.end(), vJson.begin() );
  }


  // list of allowed runs: if not empty,
  // only events from listed runs are considered
  // Superseded by the JSON file machinery
  std::vector<int> const & runs = ljmetParams.getParameter<std::vector<int> >("runs");

  HcalLaserEventFilter2012 *hcal = 0;
  
  if ( (!isMc) && (inputs.exists("useHcalLaserEventFilter")) && 
	(inputs.getParameter<bool>("useHcalLaserEventFilter")==true) ) {
    std::string name;
    if (inputs.exists("hcalLaserEventFilterList")) {
      name = inputs.getParameter<std::string>("hcalLaserEventFilterList");
    } else {
      char * release = getenv ("CMSSW_RELEASE_BASE");
      name = string(release) + "/src/EventFilter/HcalRawToDigi/data/HCALLaser2012AllDatasets.txt.gz";
    }
    cout << "Will apply HcalLaserEventFilter with list from "<< name<<endl;
    fexists(name,true);
    hcal = new HcalLaserEventFilter2012(name);
  }

  // This object 'event' is used both to get all information from the
  // event as well as to store histograms, etc.
  std::cout << legend << "Setting up chain event" << std::endl;
  fwlite::ChainEvent ev ( inputs.getParameter<std::vector<std::string> > ("fileNames") );



  //=============================================================>
  //
  // event loop
  //
  std::cout << legend << "Begin loop over events" << std::endl;
  int nev = 0;
  bool firstEvent = true;
  for (ev.toBegin(); 
       !ev.atEnd() && nev!=maxEvents; 
       ++ev, ++nev, firstEvent=false) {
    
    // skip specified number of events
    if (firstEvent && nEventsToSkip != 0){
      if (nEventsToSkip < ev.size()){
	std::cout << "Skipping " << nEventsToSkip << "events..." <<  std::endl;	
	ev.to(nEventsToSkip);
	nev += nEventsToSkip;
      }
      else{
	std::cout << legend << "Cannot skip " << nEventsToSkip << "events, it is more than I have: " << ev.size() << std::endl;	
      }
    }
    
    
    // current event
    edm::EventBase const & event = ev;

    // count event before any selection
    hists["nevents"]->Fill(1);

    // progress printout
    if ( nev % 100 == 0 ) std::cout << legend << nev << " events processed. Processing run " << event.id().run() << ", event " << event.id().event() << std::endl;


    if ( (!isMc) ){

    // Check if the event is in the JSON file or in the HCAL list
    // check if the run needs to be processed
      if (! JsonContainsEvent (vJson, event) ) continue;
      else if ( runs.size() > 0 &&
		find( runs.begin(),
		      runs.end(),
		      event.id().run() ) == runs.end() ) continue;
       if (hcal!=0){
         if (!hcal->filter (event.id().run(), event.id().luminosityBlock(), event.id().event())) continue;
      }
    }


    //
    //_____ Run private begin-of-event methods ___________________
    //
    factory->RunBeginEvent(event, ec);



    // run producers
    factory->RunAllProducers(event, theSelector);



    // event selection
    pat::strbitset ret = theSelector->getBitTemplate();
    bool passed = (*theSelector)( event, ret );


    if ( passed ) {

      //
      //_____ Run all variable calculators now ___________________
      //
      factory->RunAllCalculators(event, theSelector, ec);
      

      //
      //_____ Run selector-specific code if any___________________
      //
      theSelector->AnalyzeEvent(event, ec);
      


      //
      //_____ Run private end-of-event methods ___________________
      //
      factory->RunEndEvent(event, ec);



      //
      //_____Fill output file ____________________________________
      //
      ec.Fill();

    } // end if statement for final cut requirements
    


  } // end loop over events
  
  
  std::cout << legend << "Selection" << std::endl;
  theSelector->print(std::cout);
  theSelector->print(_logfile);
  

  _logfile.close();



  // Run EndJob() for calculators
  factory->EndJobAllCalc();



  // EndJob() for the selector
  theSelector->EndJob();



  delete theSelector;
  
  return 0;
}



bool JsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        const edm::EventBase &event)
{
  //
  // check if the event is in the good lumi list (JSON)
  // return result
  //
  
  if ( jsonVec.empty() ) return true;
  

  bool (* funcPtr) (edm::LuminosityBlockRange const &,
		    edm::LuminosityBlockID const &) = &edm::contains;
  
  edm::LuminosityBlockID lumiID (event.id().run(), 
				 event.id().luminosityBlock());
  
  std::vector< edm::LuminosityBlockRange >::const_iterator iter = 
    std::find_if (jsonVec.begin(), jsonVec.end(),
		  boost::bind(funcPtr, _1, lumiID) );


  return jsonVec.end() != iter;
}
