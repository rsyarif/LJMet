#ifndef PhysicsTools_PatUtils_interface_PFMuonSelector_h
#define PhysicsTools_PatUtils_interface_PFMuonSelector_h

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>

class PFMuonSelector : public Selector<pat::Muon> {

 public: // interface

  bool verbose_;
  
  enum Version_t { SPRING11, N_VERSIONS };

  PFMuonSelector() {}

  PFMuonSelector( edm::ParameterSet const & parameters ) {

    verbose_ = false;
    
    std::string versionStr = parameters.getParameter<std::string>("version");

    Version_t version = N_VERSIONS;

    if ( versionStr == "SPRING11" ) {
      version = SPRING11;
    }
    else {
      throw cms::Exception("InvalidInput") << "Expect version to be one of SPRING11" << std::endl;
    }

    initialize( version, 
		parameters.getParameter<bool>  ("GlobalMuon"),
		parameters.getParameter<bool>  ("TrackerMuon"),
		parameters.getParameter<double>("Chi2"),
		parameters.getParameter<double>("D0")  ,
		parameters.getParameter<int>   ("NHits")   ,
		parameters.getParameter<int>   ("NValMuHits"),
		parameters.getParameter<double>("PFIso"),
		parameters.getParameter<int>   ("nPixelHits"),
        parameters.getParameter<int>   ("nMatchedStations"),
		parameters.getParameter<int>   ("nLayersWithMeasurement")
		);
    if ( parameters.exists("cutsToIgnore") )
      setIgnoredCuts( parameters.getParameter<std::vector<std::string> >("cutsToIgnore") );
	
    retInternal_ = getBitTemplate();

  }



  void initialize( Version_t version,
		   double globalMuon = true,
		   double trackerMuon = true,
		   double chi2 = 10.0,
		   double d0 = 0.02,
		   int nhits = 11,
		   int nValidMuonHits = 0,
		   double pfiso = 0.15,
		   int minPixelHits = 1,
           int minNMatches = 1,
           int minNLayers = 6
 )
  {
    version_ = version; 

    push_back("GlobalMuon",       globalMuon);
    push_back("TrackerMuon",      trackerMuon);
    push_back("Chi2",             chi2   );
    push_back("D0",               d0     );
    push_back("NHits",            nhits  );
    push_back("NValMuHits",       nValidMuonHits  );
    push_back("PFIso",            pfiso );
    push_back("nPixelHits",       minPixelHits);
    push_back("nMatchedStations", minNMatches);
    push_back("nLayersWithMeasurement", minNLayers);

    set("GlobalMuon");
    set("TrackerMuon");
    set("Chi2");
    set("D0");
    set("NHits");
    set("NValMuHits");
    set("PFIso");   
    set("nPixelHits");
    set("nMatchedStations");  
    set("nLayersWithMeasurement");  

    indexChi2_          = index_type(&bits_, "Chi2"            );
    indexD0_            = index_type(&bits_, "D0"              );
    indexNHits_         = index_type(&bits_, "NHits"           );
    indexNValMuHits_    = index_type(&bits_, "NValMuHits"      );
    indexPFIso_         = index_type(&bits_, "PFIso"           );
    indexPixHits_       = index_type(&bits_, "nPixelHits"      );
    indexStations_      = index_type(&bits_, "nMatchedStations");
    indexLayers_        = index_type(&bits_, "nLayersWithMeasurement");

  }

  // Allow for multiple definitions of the cuts. 
  bool operator()( const pat::Muon & muon, pat::strbitset & ret ) 
  { 
    if (version_ == SPRING11 ) return spring11Cuts(muon, ret);
    else {
      return false;
    }
  }

  using Selector<pat::Muon>::operator();

  bool spring11Cuts( const pat::Muon & muon, pat::strbitset & ret)
  {
    ret.set(false);
    
    bool isGlobal  = muon.isGlobalMuon();
    bool isTracker = muon.isTrackerMuon();

    double norm_chi2 = 9999999.0;
    if ( muon.globalTrack().isNonnull() && muon.globalTrack().isAvailable() )
      norm_chi2 = muon.normChi2();
    double corr_d0 = 999999.0;
    if ( muon.globalTrack().isNonnull() && muon.globalTrack().isAvailable() )    
      corr_d0 = muon.dB();

    int nhits = static_cast<int>( muon.innerTrack()->numberOfValidHits() );
    int nValidMuonHits = 0;
    if ( muon.globalTrack().isNonnull() && muon.globalTrack().isAvailable() )
      nValidMuonHits = static_cast<int> (muon.globalTrack()->hitPattern().numberOfValidMuonHits());

    double chIso = muon.userIsolation(pat::PfChargedHadronIso);
    double nhIso = muon.userIsolation(pat::PfNeutralHadronIso);
    double gIso  = muon.userIsolation(pat::PfGammaIso);
    double puIso  = muon.userIsolation(pat::PfPUChargedHadronIso);
    double pt    = muon.pt() ;

    double pfIso = (chIso + std::max(0.,nhIso + gIso - 0.5*puIso))/pt;
    //double pfIso = (chIso + nhIso + gIso) / pt;

    int nPixelHits = 0;
    if ( muon.innerTrack().isNonnull() && muon.innerTrack().isAvailable() )
      nPixelHits = muon.innerTrack()->hitPattern().pixelLayersWithMeasurement();

    int nMatchedStations = muon.numberOfMatches();

    int nLayersWithMeasurement =0;
    if ( muon.innerTrack().isNonnull() && muon.innerTrack().isAvailable() )
      nLayersWithMeasurement = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement();

    if ( isGlobal  || ignoreCut("GlobalMuon")  )  passCut(ret, "GlobalMuon" );
    if ( isTracker  || ignoreCut("TrackerMuon")  )  passCut(ret, "TrackerMuon" );
    if ( norm_chi2        <  cut(indexChi2_,   double()) || ignoreCut(indexChi2_)    ) passCut(ret, indexChi2_   );
    if ( nhits            >= cut(indexNHits_,  int()   ) || ignoreCut(indexNHits_)   ) passCut(ret, indexNHits_  );
    if ( nValidMuonHits   >= cut(indexNValMuHits_,int()) || ignoreCut(indexNValMuHits_)) passCut(ret, indexNValMuHits_  );
    if ( fabs(corr_d0)    <  cut(indexD0_,     double()) || ignoreCut(indexD0_)      ) passCut(ret, indexD0_     );
    if ( pfIso            <  cut(indexPFIso_, double())  || ignoreCut(indexPFIso_)  ) passCut(ret, indexPFIso_ );
    if ( nPixelHits       >= cut(indexPixHits_,int())    || ignoreCut(indexPixHits_))  passCut(ret, indexPixHits_);
    if ( nMatchedStations >= cut(indexStations_,int())  || ignoreCut(indexStations_))  passCut(ret, indexStations_);
    if ( nLayersWithMeasurement >= cut(indexLayers_,int())  || ignoreCut(indexLayers_))  passCut(ret, indexLayers_);

    setIgnored(ret);
    
    return (bool)ret;
  }

  
  
 private: // member variables
  
  Version_t version_;

  index_type indexChi2_;
  index_type indexD0_;
  index_type indexNHits_;
  index_type indexNValMuHits_;
  index_type indexPFIso_;
  index_type indexPixHits_;
  index_type indexStations_;
  index_type indexLayers_;


};

#endif
