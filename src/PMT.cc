
#include "PMT.h"

using namespace std;


ClassImp( PMT )

PMT::PMT()
{}

PMT::~PMT()
{}

PMT::PMT( Int_t _nbins, Double_t _min, Double_t _max, Pedestal _ped, SPEResponse _res )
{
  nbins = _nbins;
  min = _min;
  max = _max;

  ped = _ped;
  res = _res;
    
  spectrum = new TH1D( "hspectrum", "PMT spectrum; Charge [AU]; Entries", nbins, min, max );

  spectrum->SetMarkerStyle( 20 );
  spectrum->SetMarkerSize( 0.75 );
  spectrum->SetLineColor( kBlack );
  spectrum->SetMarkerColor( kBlack );
           
}

void PMT::GenSpectrum( Int_t ntot, Double_t mu )
{
  spectrum->Reset();
   
  gRandom->SetSeed(0);
  
  for ( Int_t i=0; i<ntot; i++ )
    {
      Double_t q = ped.GenQ();
      
      Int_t npe = gRandom->Poisson( mu );
      
      for ( Int_t j=0; j<npe; j++ )
	{
	  q += res.GenQ();
	  
	}
      
      spectrum->Fill( q );
      
    }
  
  return;
  
}
