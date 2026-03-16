
#ifndef OCCUPANCY_H
#define OCCUPANCY_H

#include "iostream"
#include "iomanip"
//#include "math.h"

#include "TObject.h"
#include "TMath.h"

#include "TH1D.h"
#include "TGraph.h"

class Occupancy : public TObject
{
 private:
  
  Double_t Q0;
  Double_t s0;

 public:
  
  Occupancy();
  virtual ~Occupancy();
  Occupancy( Double_t _Q0, Double_t _s0 );

  Double_t Gauss1( Float_t x );
  Float_t FindG( TH1D* _h, Float_t f );
    
  ClassDef( Occupancy, 1 )
  
};

#endif
