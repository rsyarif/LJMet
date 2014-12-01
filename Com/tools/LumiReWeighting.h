#ifndef PhysicsTools_Utilities_interface_LumiReWeighting_cc
#define PhysicsTools_Utilities_interface_LumiReWeighting_cc

#include "TH1.h"
#include "TFile.h"
#include <string>
#include <iostream>

class LumiReWeighting {
public:
  LumiReWeighting( std::string generatedFile,
		   std::string dataFile,
		   std::string GenHistName = "pileup",
		   std::string DataHistName = "pileup" ) :
      generatedFileName_( generatedFile), 
      dataFileName_     ( dataFile ), 
      GenHistName_        ( GenHistName ), 
      DataHistName_        ( DataHistName )
      {
	generatedFile_ =  new TFile(generatedFileName_.c_str()) ; //MC distribution
	dataFile_      =  new TFile(dataFileName_.c_str()) ;	   //Data distribution

	Data_distr_ = (TH1F*) dataFile_->Get( DataHistName_.c_str() )->Clone();
	MC_distr_   = (TH1F*) generatedFile_->Get( GenHistName_.c_str() )->Clone();

	// MC * data/MC = data, so the weights are data/MC:

	// normalize both histograms first
	Data_distr_->Scale( 1.0/ Data_distr_->Integral() );
	MC_distr_->Scale( 1.0/ MC_distr_->Integral() );
	std::cout << Data_distr_->Integral()<<" "<<MC_distr_->Integral()<<std::endl;

	weights_ =  static_cast<TH1*>(Data_distr_->Clone()) ;

	weights_->SetName("lumiWeights");

	TH1* den = dynamic_cast<TH1*>(MC_distr_->Clone());

	//den->Scale(1.0/ den->Integral());
	std::cout << weights_->Integral()<<" "<<den->Integral()<<std::endl;

	weights_->Divide( den );  // so now the average weight should be 1.0
	std::cout << weights_->Integral()<<" "<<den->Integral()<<std::endl;

	std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;

	int NBins = weights_->GetNbinsX();

	for(int ibin = 1; ibin<NBins+1; ++ibin){
	  std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
	}
};

  double weight( int npv ) {
    int bin = weights_->GetXaxis()->FindBin( npv );
    return weights_->GetBinContent( bin );
  }

  double weight( float npv ) {
    int bin = weights_->GetXaxis()->FindBin( npv );
    return weights_->GetBinContent( bin );
  }

//   double weight( const edm::EventBase &e ) ;
// 
//   double weightOOT( const edm::EventBase &e ) ;
// 
//   void weightOOT_init(); 

protected:

  std::string generatedFileName_;
  std::string dataFileName_;
  std::string GenHistName_;
  std::string DataHistName_;
  TFile*     generatedFile_;
  TFile*     dataFile_;
  TH1*	  weights_;

  //keep copies of normalized distributions:

  TH1*	   MC_distr_;
  TH1*	   Data_distr_;
};

#endif
