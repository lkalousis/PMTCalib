
#ifndef PMT_H
#define PMT_H

#include "iostream"
#include "iomanip"

#include "TObject.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"

#include "Pedestal.h"
#include "SPEResponse.h"

class PMT : public TObject
{
 private:
  
  TH1D *spectrum;
  
  Int_t nbins;
  Double_t min;
  Double_t max;
  
  Pedestal ped;
  SPEResponse res;
  
  
 public:
  
  PMT();
  
  virtual ~PMT();
  
  PMT( Int_t _nbins, Double_t _min, Double_t _max, Pedestal _ped, SPEResponse _res );
      
  void GenSpectrum( Int_t ntot, Double_t mu );
  
  TH1D* GetSpectrum(){ return spectrum; };
  
  void DrawSpectrum()
  {
    spectrum->SetMaximum( 2.5*spectrum->GetBinContent( spectrum->GetMaximumBin() ) );
    spectrum->Draw( "PEZ" );
      
  };
    
  ClassDef( PMT, 1 )
    
};

#endif
