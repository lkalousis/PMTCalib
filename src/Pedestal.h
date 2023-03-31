
#ifndef PEDESTAL_H
#define PEDESTAL_H

#include "TObject.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TString.h"
#include "TH1D.h"
#include "TMinuit.h"

class Pedestal : public TObject
{
 private:
    
  TF1 *pedfunc;
  TF1 *pedfit;
      
 public:
  
  Pedestal();
  
  virtual ~Pedestal();

  Double_t wbin;
  
  Double_t Q0;
  Double_t s0;

  Double_t dQ0;
  Double_t ds0;
  
  Pedestal( Double_t _Q0, Double_t _s0 );
  
  Double_t GenQ();

  void LocatePedestal( TH1 *hspec, Double_t _Q0, Double_t _s0 );
  Int_t status;
  Double_t chi2;
  
  ClassDef( Pedestal, 1 )
        
};

#endif
